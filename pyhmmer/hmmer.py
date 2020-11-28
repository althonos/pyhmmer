# coding: utf-8

import contextlib
import ctypes
import queue
import time
import threading
import typing
import os
import multiprocessing

from .easel import Alphabet, DigitalSequence, TextSequence, SequenceFile, SSIWriter
from .plan7 import Background, Pipeline, TopHits, HMM, HMMFile, Profile


class _PipelineThread(threading.Thread):
    @staticmethod
    def _none_callback(hmm: HMM, total: int) -> None:
        pass

    def __init__(
        self,
        sequences: typing.Iterable[DigitalSequence],
        hmm_queue: "queue.Queue[typing.Optional[typing.Tuple[int, HMM]]]",
        hmm_count: multiprocessing.Value,  # type: ignore
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[HMM, int], None]],
        options: typing.Dict[str, typing.Any],
    ) -> None:
        super().__init__()
        self.options = options
        self.pipeline = Pipeline(alphabet=Alphabet.amino(), **options)
        self.hits = TopHits()
        self.sequences = sequences
        self.hmm_queue = hmm_queue
        self.hmm_count = hmm_count
        self.hits_queue = hits_queue
        self.callback = callback or self._none_callback
        self.kill_switch = kill_switch
        self.error = None  # type: typing.Optional[BaseException]

    def run(self) -> None:
        while not self.kill_switch.is_set():
            args = self.hmm_queue.get()
            if args is None:
                self.hmm_queue.task_done()
                return
            else:
                index, hmm = args
            try:
                hits = self.pipeline.search(hmm, self.sequences)
                self.hits_queue.put((index, hits))
                self.hmm_queue.task_done()
                self.callback(hmm, self.hmm_count.value)  # type: ignore
                self.pipeline.clear()
            except BaseException as exc:
                self.error = exc
                self.kill()
                return

    def kill(self) -> None:
        self.kill_switch.set()


def hmmsearch(
    hmms: typing.Iterable[HMM],
    sequences: typing.Sequence[DigitalSequence],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
    **options: typing.Any,
) -> typing.Iterator[TopHits]:
    # count the number of CPUs to use
    _cpus = cpus if cpus > 0 else multiprocessing.cpu_count()

    # create the queues to pass the HMM objects around, as well as atomic
    # values that we use to synchronize the threads
    hits_queue = queue.PriorityQueue()  # type: ignore
    hmm_queue = queue.Queue()  # type: ignore
    hmm_count = multiprocessing.Value(ctypes.c_ulong)
    kill_switch = threading.Event()

    # create and launch one pipeline thread per CPU
    threads = []
    for _ in range(_cpus):
        thread = _PipelineThread(
            sequences, hmm_queue, hmm_count, hits_queue, kill_switch, callback, options
        )
        thread.start()
        threads.append(thread)

    # queue the HMMs passed as arguments
    for index, hmm in enumerate(hmms):
        hmm_count.value += 1
        hmm_queue.put((index, hmm))

    # poison-pill the queue so that threads terminate when they
    # have consumed all the HMMs
    for _ in threads:
        hmm_queue.put(None)

    # wait for all threads to be completed
    for thread in threads:
        thread.join()
        if thread.error is not None:
            raise thread.error

    # give back results
    while not hits_queue.empty():
        yield hits_queue.get_nowait()[1]


def hmmpress(hmms: typing.Iterable[HMM], output: str) -> None:

    DEFAULT_L = 400

    with contextlib.ExitStack() as ctx:
        h3p = ctx.enter_context(open("{}.h3p".format(output), "wb"))
        h3m = ctx.enter_context(open("{}.h3m".format(output), "wb"))
        h3f = ctx.enter_context(open("{}.h3f".format(output), "wb"))
        h3i = ctx.enter_context(SSIWriter("{}.h3i".format(output)))
        fh = h3i.add_file(output, format=0)

        for nmodel, hmm in enumerate(hmms):
            # create the background model on the first iteration
            if nmodel == 0:
                bg = Background(hmm.alphabet)
                bg.L = DEFAULT_L

            # build the optimized models
            gm = Profile(hmm.M, hmm.alphabet)
            gm.configure(hmm, bg, DEFAULT_L)
            om = gm.optimized()

            # update the disk offsets of the optimized model to be written
            om.offsets.model = h3m.tell()
            om.offsets.profile = h3p.tell()
            om.offsets.filter = h3f.tell()

            # add the HMM name, and optionally the HMM accession to the index
            h3i.add_key(hmm.name, fh, h3m_offset, 0, 0)
            if hmm.accession is not None:
                h3i.add_alias(hmm.accession, hmm.name)

            # write the HMM in binary format, and the optimized profile
            hmm.write(h3m, binary=True)
            om.write(h3f, h3p)


# add a very limited CLI so that this module can be invoked in a shell:
#     $ python -m pyhmmer.hmmsearch <hmmfile> <seqdb>
if __name__ == "__main__":

    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", required=False, default=0, type=int)
    parser.add_argument("hmmfile")
    parser.add_argument("seqdb")
    args = parser.parse_args()

    with SequenceFile(args.seqdb) as seqfile:
        alphabet = seqfile.guess_alphabet()
        if alphabet is None:
            print("could not guess alphabet of input, exiting")
            sys.exit(1)

        seq = TextSequence()
        sequences = []
        while seqfile.readinto(seq) is not None:
            sequences.append(seq.digitize(alphabet))
            seq.clear()

        with HMMFile(args.hmmfile) as hmms:
            hits_list = hmmsearch(hmms, sequences, cpus=args.jobs)

            for hits in hits_list:
                for hit in hits:
                    if hit.is_reported():
                        print(
                            hit.name.decode(),
                            "-",
                            hit.domains[0].alignment.hmm_name.decode(),
                            hit.evalue,
                            hit.score,
                            hit.bias,
                            sep="\t",
                        )
