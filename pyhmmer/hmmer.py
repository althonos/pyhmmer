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


def _hmmsearch_singlethreaded(
    hmms: typing.Iterable[HMM],
    sequences: typing.Sequence[DigitalSequence],
    callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
    **options  # type: typing.Any
) -> typing.Iterator:
    # create the queues to pass the HMM objects around, as well as atomic
    # values that we use to synchronize the threads
    hits_queue = queue.PriorityQueue()  # type: ignore
    hmm_queue = queue.Queue()  # type: ignore
    hmm_count = multiprocessing.Value(ctypes.c_ulong)
    kill_switch = threading.Event()

    # create the thread (to recycle code)
    thread = _PipelineThread(
        sequences,
        hmm_queue,
        hmm_count,
        hits_queue,
        kill_switch,
        callback,
        options
    )

    # queue the HMMs passed as arguments
    for index, hmm in enumerate(hmms):
        hmm_count.value += 1
        hmm_queue.put((index, hmm))

    # poison-pill the queue so that threads terminate when they
    # have consumed all the HMMs
    hmm_queue.put(None)

    # launch the thread code, but in the main thread
    thread.run()
    if thread.error is not None:
        raise thread.error

    # give back results
    while not hits_queue.empty():
        yield hits_queue.get_nowait()[1]

def _hmmsearch_multithreaded(
    hmms: typing.Iterable[HMM],
    sequences: typing.Sequence[DigitalSequence],
    cpus: int,
    callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
    **options  # type: typing.Any
) -> typing.Iterator:
    # create the queues to pass the HMM objects around, as well as atomic
    # values that we use to synchronize the threads
    hits_queue = queue.PriorityQueue()  # type: ignore
    hmm_queue = queue.Queue()  # type: ignore
    hmm_count = multiprocessing.Value(ctypes.c_ulong)
    kill_switch = threading.Event()

    # create and launch one pipeline thread per CPU
    threads = []
    for _ in range(cpus):
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


def hmmsearch(
    hmms: typing.Iterable[HMM],
    sequences: typing.Sequence[DigitalSequence],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
    **options  # type: typing.Any
) -> typing.Iterator[TopHits]:
    # count the number of CPUs to use
    _cpus = cpus if cpus > 0 else multiprocessing.cpu_count()

    if _cpus > 1:
        return _hmmsearch_multithreaded(hmms, sequences, _cpus, callback, **options)
    else:
        return _hmmsearch_singlethreaded(hmms, sequences, callback, **options)



def hmmpress(
    hmms: typing.Iterable[HMM],
    output: typing.Union[str, "os.PathLike[str]"],
) -> int:

    DEFAULT_L = 400
    path = os.fspath(output)
    nmodel = 0

    with contextlib.ExitStack() as ctx:
        h3p = ctx.enter_context(open("{}.h3p".format(path), "wb"))
        h3m = ctx.enter_context(open("{}.h3m".format(path), "wb"))
        h3f = ctx.enter_context(open("{}.h3f".format(path), "wb"))
        h3i = ctx.enter_context(SSIWriter("{}.h3i".format(path)))
        fh = h3i.add_file(path, format=0)

        for hmm in hmms:
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
            h3i.add_key(hmm.name, fh, om.offsets.model, 0, 0)
            if hmm.accession is not None:
                h3i.add_alias(hmm.accession, hmm.name)

            # write the HMM in binary format, and the optimized profile
            hmm.write(h3m, binary=True)
            om.write(h3f, h3p)
            nmodel += 1

    # return the number of written HMMs
    return nmodel


# add a very limited CLI so that this module can be invoked in a shell:
#     $ python -m pyhmmer.hmmsearch <hmmfile> <seqdb>
if __name__ == "__main__":

    import argparse
    import sys


    def _hmmsearch(args: argparse.Namespace) -> int:
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

        return 0


    def _hmmpress(args: argparse.Namespace) -> int:
        for ext in ["h3m", "h3i", "h3f", "h3p"]:
            path = "{}.{}".format(args.hmmfile, ext)
            if os.path.exists(path):
                if args.force:
                    os.remove(path)
                else:
                    raise FileExistsError(path)

        with HMMFile(args.hmmfile) as hmms:
            hmmpress(hmms, args.hmmfile)

        return 0


    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", required=False, default=0, type=int)
    subparsers = parser.add_subparsers(dest="cmd", help='HMMER command to run', required=True)

    parser_hmmsearch = subparsers.add_parser("hmmsearch")
    parser_hmmsearch.add_argument("hmmfile")
    parser_hmmsearch.add_argument("seqdb")

    parser_hmmpress = subparsers.add_parser("hmmpress")
    parser_hmmpress.add_argument("hmmfile")
    parser_hmmpress.add_argument("-f", "--force", action="store_true")

    args = parser.parse_args()
    if args.cmd == "hmmsearch":
        sys.exit(_hmmsearch(args))
    elif args.cmd == "hmmpress":
        sys.exit(_hmmpress(args))
