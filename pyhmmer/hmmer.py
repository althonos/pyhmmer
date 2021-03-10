# coding: utf-8
"""Reimplementation of HMMER binaries with the PyHMMER API.
"""

import abc
import contextlib
import collections
import ctypes
import queue
import time
import threading
import typing
import os
import multiprocessing

import psutil

from .easel import Alphabet, DigitalSequence, DigitalMSA, TextSequence, SequenceFile, SSIWriter
from .plan7 import Builder, Background, Pipeline, TopHits, HMM, HMMFile, Profile
from .utils import peekable

# the query type for the pipeline
_Q = typing.TypeVar("_Q")

# --- Pipeline threads -------------------------------------------------------

class _PipelineThread(typing.Generic[_Q], threading.Thread):
    @staticmethod
    def _none_callback(hmm: HMM, total: int) -> None:
        pass

    def __init__(
        self,
        sequences: typing.Iterable[DigitalSequence],
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, _Q]]]",
        query_count: multiprocessing.Value,  # type: ignore
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[_Q, int], None]],
        options: typing.Dict[str, typing.Any],
    ) -> None:
        super().__init__()
        self.options = options
        self.pipeline = Pipeline(alphabet=Alphabet.amino(), **options)
        self.sequences = sequences
        self.query_queue = query_queue
        self.query_count = query_count
        self.hits_queue = hits_queue
        self.callback = callback or self._none_callback
        self.kill_switch = kill_switch
        self.error: typing.Optional[BaseException] = None

    def run(self) -> None:
        while not self.kill_switch.is_set():
            args = self.query_queue.get()
            if args is None:
                self.query_queue.task_done()
                return
            else:
                index, query = args
            try:
                hits = self.search(query)
                self.hits_queue.put((index, hits))
                self.query_queue.task_done()
                self.callback(query, self.query_count.value)  # type: ignore
                self.pipeline.clear()
            except BaseException as exc:
                self.error = exc
                self.kill()
                return

    def kill(self) -> None:
        self.kill_switch.set()

    @abc.abstractmethod
    def search(self, query: _Q) -> TopHits:
        return NotImplemented


class _HMMPipelineThread(_PipelineThread[HMM]):
    def search(self, query: HMM) -> TopHits:
        return self.pipeline.search_hmm(query, self.sequences)


class _SequencePipelineThread(_PipelineThread[DigitalSequence]):
    def __init__(
        self,
        sequences: typing.Iterable[DigitalSequence],
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, DigitalSequence]]]",
        query_count: multiprocessing.Value,  # type: ignore
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[DigitalSequence, int], None]],
        options: typing.Dict[str, typing.Any],
        builder: Builder,
    ) -> None:
        super().__init__(
            sequences,
            query_queue,
            query_count,
            hits_queue,
            kill_switch,
            callback,
            options,
        )
        self.builder = builder

    def search(self, query: DigitalSequence) -> TopHits:
        return self.pipeline.search_seq(query, self.sequences, self.builder)


class _MSAPipelineThread(_PipelineThread[DigitalMSA]):
    def __init__(
        self,
        sequences: typing.Iterable[DigitalSequence],
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, DigitalMSA]]]",
        query_count: multiprocessing.Value,  # type: ignore
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[DigitalMSA, int], None]],
        options: typing.Dict[str, typing.Any],
        builder: Builder,
    ) -> None:
        super().__init__(
            sequences,
            query_queue,
            query_count,
            hits_queue,
            kill_switch,
            callback,
            options,
        )
        self.builder = builder

    def search(self, query: DigitalMSA) -> TopHits:
        return self.pipeline.search_msa(query, self.sequences, self.builder)


# --- Search runners ---------------------------------------------------------

@abc.abstractmethod
class _Search(typing.Generic[_Q], abc.ABC):

    def __init__(
        self,
        queries: typing.Iterable[_Q],
        sequences: typing.Collection[DigitalSequence],
        cpus: int = 0,
        callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
        **options,
    ) -> None:
        self.queries = queries
        self.sequences = sequences
        self.cpus = cpus
        self.callback = callback
        self.options = options

    @abc.abstractmethod
    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, _Q]]]",
        query_count: "multiprocessing.Value[int]",
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
    ) -> _PipelineThread[_Q]:
        return NotImplemented  # type: ignore

    def _single_threaded(self) -> typing.Iterator[TopHits]:
        # create the queues to pass the HMM objects around, as well as atomic
        # values that we use to synchronize the threads
        query_queue = queue.Queue()  # type: ignore
        query_count = multiprocessing.Value(ctypes.c_ulong)
        hits_queue = queue.PriorityQueue()  # type: ignore
        kill_switch = threading.Event()

        # create the thread (to recycle code)
        thread = self._new_thread(query_queue, query_count, hits_queue, kill_switch)

        # queue the HMMs passed as arguments
        for index, query in enumerate(self.queries):
            query_count.value += 1
            query_queue.put((index, query))

        # poison-pill the queue so that threads terminate when they
        # have consumed all the HMMs
        query_queue.put(None)

        # launch the thread code, but in the main thread
        thread.run()
        if thread.error is not None:
            raise thread.error

        # give back results
        while not hits_queue.empty():
            yield hits_queue.get_nowait()[1]

    def _multi_threaded(self) -> typing.Iterator[TopHits]:
        # create the queues to pass the HMM objects around, as well as atomic
        # values that we use to synchronize the threads
        hits_queue = queue.PriorityQueue()  # type: ignore
        query_queue = queue.Queue()  # type: ignore
        query_count = multiprocessing.Value(ctypes.c_ulong)
        kill_switch = threading.Event()

        # create and launch one pipeline thread per CPU
        threads = []
        for _ in range(self.cpus):
            thread = self._new_thread(query_queue, query_count, hits_queue, kill_switch)
            thread.start()
            threads.append(thread)

        # queue the HMMs passed as arguments
        for index, query in enumerate(self.queries):
            query_count.value += 1
            query_queue.put((index, query))

        # poison-pill the queue so that threads terminate when they
        # have consumed all the HMMs
        for _ in threads:
            query_queue.put(None)

        # wait for all threads to be completed
        for thread in threads:
            thread.join()
            if thread.error is not None:
                raise thread.error

        # give back results
        while not hits_queue.empty():
            yield hits_queue.get_nowait()[1]

    def run(self) -> typing.Iterator[TopHits]:
        if self.cpus == 1:
            return self._single_threaded()
        else:
            return self._multi_threaded()


class _HMMSearch(_Search[HMM]):

    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, HMM]]]",
        query_count: "multiprocessing.Value[int]",
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
    ) -> _HMMPipelineThread:
        return _HMMPipelineThread(
            self.sequences,
            query_queue,
            query_count,
            hits_queue,
            kill_switch,
            self.callback,
            self.options,
        )


class _SequenceSearch(_Search[DigitalSequence]):

    def __init__(
        self,
        builder: Builder,
        queries: typing.Iterable[DigitalSequence],
        sequences: typing.Collection[DigitalSequence],
        cpus: int = 0,
        callback: typing.Optional[typing.Callable[[DigitalSequence, int], None]] = None,
        **options,
    ) -> None:
        super().__init__(queries, sequences, cpus, callback, **options)
        self.builder = builder

    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, DigitalSequence]]]",
        query_count: "multiprocessing.Value[int]",
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
    ) -> _SequencePipelineThread:
        return _SequencePipelineThread(
            self.sequences,
            query_queue,
            query_count,
            hits_queue,
            kill_switch,
            self.callback,
            self.options,
            self.builder.copy(),
        )


class _MSASearch(_Search[DigitalMSA]):

    def __init__(
        self,
        builder: Builder,
        queries: typing.Iterable[DigitalMSA],
        sequences: typing.Collection[DigitalSequence],
        cpus: int = 0,
        callback: typing.Optional[typing.Callable[[DigitalMSA, int], None]] = None,
        **options,
    ) -> None:
        super().__init__(queries, sequences, cpus, callback, **options)
        self.builder = builder

    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, DigitalMSA]]]",
        query_count: "multiprocessing.Value[int]",
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
    ) -> _MSAPipelineThread:
        return _MSAPipelineThread(
            self.sequences,
            query_queue,
            query_count,
            hits_queue,
            kill_switch,
            self.callback,
            self.options,
            self.builder.copy(),
        )


# --- hmmsearch --------------------------------------------------------------

def hmmsearch(
    queries: typing.Iterable[HMM],
    sequences: typing.Collection[DigitalSequence],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
    **options,  # type: typing.Any
) -> typing.Iterator[TopHits]:
    """Search HMM profiles against a sequence database.

    Arguments:
        queries (iterable of `~pyhmmer.plan7.HMM`): The query HMMs to
            search in the database.
        sequences (collection of `~pyhmmer.easel.DigitalSequence`): A
            database of sequences to query.
        cpus (`int`): The number of threads to run in parallel. Pass ``1``
            to run everything in the main thread, ``0`` to automatically
            select a suitable number (using `psutil.cpu_count`), or any
            positive number otherwise.
        callback (callable): A callback that is called everytime a query is
            processed with two arguments: the query, and the total number
            of queries. This can be used to display progress in UI.

    Yields:
        `~pyhmmer.plan7.TopHits`: An object reporting *top hits* for each
        query, in the same order the queries were passed in the input.

    Raises:
        `~pyhmmer.errors.AlphabetMismatch`: When any of the query HMMs
        and the sequences do not share the same alphabet.

    Note:
        Any additional arguments passed to the `hmmsearch` function will be
        passed transparently to the `~pyhmmer.plan7.Pipeline` to be created.

    .. versionadded:: 0.1.0

    """
    # count the number of CPUs to use
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or multiprocessing.cpu_count()
    runner = _HMMSearch(queries, sequences, _cpus, callback, **options)
    return runner.run()


# --- phmmer -----------------------------------------------------------------

@typing.overload
def phmmer(
    queries: typing.Iterable[DigitalMSA],
    sequences: typing.Collection[DigitalSequence],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[DigitalMSA, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options: typing.Any,
) -> typing.Iterator[TopHits]:
    ...

def phmmer(
    queries: typing.Iterable[DigitalSequence],
    sequences: typing.Collection[DigitalSequence],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[DigitalSequence, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options: typing.Any,
) -> typing.Iterator[TopHits]:
    """Search protein sequences against a sequence database.

    Arguments:
        queries (iterable of `DigitalSequence` or `DigitalMSA`): The query
            sequences to search in the database.
        sequences (collection of `~pyhmmer.easel.DigitalSequence`): A
            database of sequences to query.
        cpus (`int`): The number of threads to run in parallel. Pass ``1`` to
            run everything in the main thread, ``0`` to automatically
            select a suitable number (using `psutil.cpu_count`), or any
            positive number otherwise.
        callback (callable): A callback that is called everytime a query is
            processed with two arguments: the query, and the total number
            of queries. This can be used to display progress in UI.
        builder (`~pyhmmer.plan7.Builder`, optional): A builder to configure
            how the queries are converted to HMMs. Passing `None` will create
            a default instance.

    Yields:
        `~pyhmmer.plan7.TopHits`: A *top hits* instance for each query,
        in the same order the queries were passed in the input.

    Note:
        Any additional keyword arguments passed to the `phmmer` function
        will be passed transparently to the `~pyhmmer.plan7.Pipeline` to
        be created in each worker thread.

    .. versionadded:: 0.2.0

    .. versionchanged:: 0.3.0
       Allow using `DigitalMSA` queries.

    """
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or multiprocessing.cpu_count()
    _builder = Builder(Alphabet.amino()) if builder is None else builder

    try:
        _queries = peekable(queries)
        _item = _queries.peek()
    except StopIteration:
        _item = None

    if isinstance(_item, DigitalSequence):
        runner = _SequenceSearch(
            _builder, _queries, sequences, _cpus, callback, **options
        )
    else:
        runner = _MSASearch(
            _builder, _queries, sequences, _cpus, callback, **options
        )
    return runner.run()


# --- nhmmer -----------------------------------------------------------------

@typing.overload
def nhmmer(
    queries: typing.Iterable[DigitalMSA],
    sequences: typing.Collection[DigitalSequence],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[DigitalMSA, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options: typing.Any,
) -> typing.Iterator[TopHits]:
    ...

def nhmmer(
    queries: typing.Iterable[DigitalSequence],
    sequences: typing.Collection[DigitalSequence],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[DigitalSequence, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options: typing.Any,
) -> typing.Iterator[TopHits]:
    """Search protein sequences against a sequence database.

    See Also:
        The equivalent function for proteins, `~pyhmmer.hmmer.phmmer`.

    .. versionadded:: 0.3.0

    """
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or multiprocessing.cpu_count()
    _builder = Builder(Alphabet.dna()) if builder is None else builder

    try:
        _queries = peekable(queries)
        _item = _queries.peek()
    except StopIteration:
        _item = None

    if isinstance(_item, DigitalSequence):
        runner = _SequenceSearch(
            _builder, _queries, sequences, _cpus, callback, **options
        )
    else:
        runner = _MSASearch(
            _builder, _queries, sequences, _cpus, callback, **options
        )
    return runner.run()


# --- hmmpress ---------------------------------------------------------------

def hmmpress(
    hmms: typing.Iterable[HMM], output: typing.Union[str, "os.PathLike[str]"],
) -> int:
    """Press several HMMs into a database.

    Calling this function will create 4 files at the given location:
    ``{output}.h3p`` (containing the optimized profiles),
    ``{output}.h3m`` (containing the binary HMMs),
    ``{output}.h3f`` (containing the MSV parameters), and
    ``{output}.h3i`` (the SSI index mapping the previous files).

    Arguments:
        hmms (iterable of `~pyhmmer.plan7.HMM`): The HMMs to be pressed
            together in the file.
        output (`str` or `os.PathLike`): The path to an output location
            where to write the different files.

    """
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
                print("could not guess alphabet of input, exiting", file=sys.stderr)
                return 1
            seqfile.set_digital(alphabet)
            sequences = list(seqfile)

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

    def _phmmer(args: argparse.Namespace) -> int:
        with SequenceFile(args.seqdb) as seqfile:
            alphabet = seqfile.guess_alphabet() or Alphabet.amino()
            seqfile.set_digital(alphabet)
            sequences = list(seqfile)

        with SequenceFile(args.seqfile) as queries:
            queries.set_digital(alphabet)
            hits_list = phmmer(queries, sequences, cpus=args.jobs)  # type: ignore

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
                    print(f"file {path} already exists")
                    return 1

        with HMMFile(args.hmmfile) as hmms:
            hmmpress(hmms, args.hmmfile)

        return 0

    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", required=False, default=0, type=int)
    subparsers = parser.add_subparsers(
        dest="cmd", help="HMMER command to run", required=True
    )

    parser_hmmsearch = subparsers.add_parser("hmmsearch")
    parser_hmmsearch.add_argument("hmmfile")
    parser_hmmsearch.add_argument("seqdb")

    parser_phmmer = subparsers.add_parser("phmmer")
    parser_phmmer.add_argument("seqfile")
    parser_phmmer.add_argument("seqdb")

    parser_hmmpress = subparsers.add_parser("hmmpress")
    parser_hmmpress.add_argument("hmmfile")
    parser_hmmpress.add_argument("-f", "--force", action="store_true")

    args = parser.parse_args()
    if args.cmd == "hmmsearch":
        sys.exit(_hmmsearch(args))
    elif args.cmd == "hmmpress":
        sys.exit(_hmmpress(args))
    elif args.cmd == "phmmer":
        sys.exit(_phmmer(args))
