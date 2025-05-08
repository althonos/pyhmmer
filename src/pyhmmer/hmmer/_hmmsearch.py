import collections
import contextlib
import operator
import ctypes
import queue
import multiprocessing
import typing
import os
import threading

import psutil

from ..easel import Alphabet, DigitalSequence, DigitalMSA, DigitalSequenceBlock, SequenceFile
from ..plan7 import TopHits, Builder, Pipeline, HMM, Profile, OptimizedProfile, HMMPressedFile, OptimizedProfileBlock
from ..utils import singledispatchmethod, peekable
from ._base import _BaseDispatcher, _BaseWorker, _BaseChore, _AnyProfile

_SEARCHQueryType = typing.Union[_AnyProfile]
_P = typing.TypeVar("_P", HMM, Profile, OptimizedProfile)

if typing.TYPE_CHECKING:
    from ._base import Unpack, PipelineOptions, BACKEND

# --- Worker -------------------------------------------------------------------

class _SEARCHWorker(
    _BaseWorker[
        _SEARCHQueryType,
        typing.Union[DigitalSequenceBlock, "SequenceFile[DigitalSequence]"],
        "TopHits[_SEARCHQueryType]",
    ],
):
    @singledispatchmethod
    def query(self, query) -> "TopHits[Any]":  # type: ignore
        raise TypeError(
            "Unsupported query type for `hmmsearch`: {}".format(type(query).__name__)
        )

    @query.register(HMM)
    @query.register(Profile)
    @query.register(OptimizedProfile)
    def _(self, query: _AnyProfile) -> "TopHits[_AnyProfile]":  # type: ignore
        assert self.pipeline is not None
        return self.pipeline.search_hmm(query, self.targets)


class _SEARCHThread(_SEARCHWorker, threading.Thread):
    pass


class _SEARCHProcess(_SEARCHWorker, multiprocessing.Process):
    pass


# --- Dispatcher ---------------------------------------------------------------

class _SEARCHDispatcher(
    _BaseDispatcher[
        _SEARCHQueryType,
        typing.Union[DigitalSequenceBlock, "SequenceFile[DigitalSequence]"],
        "TopHits[_SEARCHQueryType]",
    ]
):
    def _new_worker(
        self,
        query_queue: "queue.Queue[typing.Optional[_BaseChore[_SEARCHQueryType, TopHits[_SEARCHQueryType]]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _SEARCHWorker:
        if isinstance(self.targets, SequenceFile):
            assert self.targets.name is not None
            targets = SequenceFile(
                self.targets.name,
                format=self.targets.format,
                digital=True,
                alphabet=self.options["alphabet"],
            )
        else:
            targets = self.targets  # type: ignore
        params = [
            targets,
            query_queue,
            query_count,
            kill_switch,
            self.callback,
            self.options,
            self.builder,
        ]
        if self.backend == "threading":
            return _SEARCHThread(*params)
        elif self.backend == "multiprocessing":
            return _SEARCHProcess(*params)
        else:
            raise ValueError(f"Invalid backend for `hmmsearch`: {self.backend!r}")


class _ReverseSEARCHDispatcher(
    _BaseDispatcher[
        "_SEARCHQueryType",
        DigitalSequenceBlock,
        "TopHits[_SEARCHQueryType]",
    ]
):
    """A ``hmmsearch`` dispatcher that parallelizes on the targets.
    """

    def __init__(
        self,
        queries: typing.Iterable["_SEARCHQueryType"],
        targets: DigitalSequenceBlock,
        cpus: int = 0,
        callback: typing.Optional[typing.Callable[["_SEARCHQueryType", int], None]] = None,
        builder: typing.Optional[Builder] = None,
        timeout: int = 1,
        backend: "BACKEND" = "threading",
        **options,  # type: Unpack[PipelineOptions]
    ) -> None:
        super().__init__(
            queries,
            targets,
            cpus,
            callback,
            builder,
            timeout,
            backend,
            **options
        )
        # only use as many CPUs as there are targets (if only a few), but
        # that may be a waste for less than N sequences per CPUs?
        self.cpus = max(1, min(cpus, len(targets)))
        # attempt to balance the chunks so that every thread gets about the
        # same number of *residues* (not the same number of *sequences*!)
        self.target_chunks = self._make_chunks(targets)

    def _make_chunks(self, targets: DigitalSequenceBlock) -> typing.List[DigitalSequenceBlock]:
        # compute chunksize from total sequence lengths
        # TODO: implement this as a Cython function with quick access to the
        #       sequence data?
        total_length = sum(len(seq) for seq in targets)
        chunksize = (total_length + self.cpus - 1) // self.cpus
        # balance sequence residues across chunks
        current_size = 0
        chunk_indices = [0]
        for i, seq in enumerate(targets):
            current_size += len(seq)
            if current_size > chunksize:
                chunk_indices.append(i)
                current_size = 0
        while len(chunk_indices) <= self.cpus:
            chunk_indices.append(len(targets))
        # NB: this does not copy data, as `DigitalSequenceBlock` are implemented
        #     as views of `DigitalSequence` objects, so slicing is cheap.
        return [targets[i:j] for i,j in zip(chunk_indices, chunk_indices[1:])]

    def _new_worker(
        self,
        query_queue,
        query_count,
        kill_switch: threading.Event,
        targets: DigitalSequenceBlock = None,
    ):
        if targets is None:
            targets = self.targets
        params = [
            targets,
            query_queue,
            query_count,
            kill_switch,
            None,
            self.options,
            self.builder,
        ]
        if self.backend == "threading":
            return _SEARCHThread(*params)
        elif self.backend == "multiprocessing":
            return _SEARCHProcess(*params)
        else:
            raise ValueError(f"Invalid backend for `hmmsearch`: {self.backend!r}")

    def _single_threaded(self) -> typing.Iterator["TopHits[_SEARCHQueryType]"]:
        for i, hits in enumerate(super()._single_threaded(), start=1):
            if self.callback is not None:
                self.callback(hits.query, i)
            yield hits

    def _multi_threaded(self) -> typing.Iterator["TopHits[_SEARCHQueryType]"]:
        with contextlib.ExitStack() as ctx:
            # single shared kill switch and query count
            if self.backend == "multiprocessing":
                manager = ctx.enter_context(multiprocessing.Manager())
                query_count = manager.Value(ctypes.c_ulong, 0)
                kill_switch = manager.Event()
            else:
                query_count = multiprocessing.Value(ctypes.c_ulong)  # type: ignore
                kill_switch = threading.Event()

            # create and launch one pipeline thread per CPU, each with its own
            # queue as they all need to get the same copy of each query
            workers = []
            queues = []
            for i in range(self.cpus):
                # create the queues to pass the query objects around, only
                # one item at a time since we synchronize query-by-query
                results: typing.Deque[_BaseChore[_Q, _R]] = collections.deque()
                if self.backend == "multiprocessing":
                    query_queue = ctx.enter_context(contextlib.closing(multiprocessing.Queue()))
                elif self.backend == "threading":
                    query_queue = queue.Queue()
                # create worker
                worker = self._new_worker(query_queue, query_count, kill_switch, targets=self.target_chunks[i])
                worker.start()
                workers.append(worker)
                queues.append(query_queue)

            # catch exceptions to kill threads in the background before exiting
            try:
                hits = None
                # process queries sequentially
                for i, query in enumerate(self.queries):
                    # yield hits obtained at the previous iteration
                    if hits is not None:
                        yield hits
                    query_count.value += 1
                    # create one chore per worker
                    chores = []
                    for worker, worker_queue in zip(workers, queues):
                        if isinstance(query, OptimizedProfile):
                            chore = self._new_chore(query.copy())
                        else:
                            chore = self._new_chore(query)
                        chores.append(chore)
                        worker_queue.put(chore)
                    # collect hits
                    partial_hits = []
                    for chore in chores:
                        partial_hits.append(chore.get())
                    # merge hits
                    hits = TopHits.merge(*partial_hits)
                    # call callback here after the hits have been merged
                    if self.callback is not None:
                        self.callback(chore.query, query_count.value)
                # now that we exhausted all queries, poison pill the
                # threads so they stop on their own gracefully
                for worker in workers:
                    worker.query_queue.put(None)
                    worker.join()
                    if self.backend == "multiprocessing":
                        worker.query_queue.close()
                # yield the final hits
                if hits is not None:
                    yield hits
            except BaseException as e:
                # make sure threads are killed to avoid being stuck,
                # e.g. after a KeyboardInterrupt, then re-raise
                try:
                    kill_switch.set()
                except queue.Full:
                    pass
                for worker in workers:
                    worker.join()
                    if self.backend == "multiprocessing":
                        worker.query_queue.close()
                        worker.query_queue.join_thread()
                raise e


# --- hmmsearch --------------------------------------------------------------

def hmmsearch(
    queries: typing.Union[_P, typing.Iterable[_P]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_P, int], None]] = None,
    backend: "BACKEND" = "threading",
    parallel: "PARALLEL" = None,
    **options,  # type: Unpack[PipelineOptions]
) -> typing.Iterator["TopHits[_P]"]:
    """Search HMM profiles against a sequence database.

    In HMMER many-to-many comparisons, a *search* is the operation of
    querying with profile HMMs a database of sequences.

    The `hmmsearch` function offers two ways of managing the database that
    will be selected based on the type of the ``sequences`` argument. If
    ``sequences`` is an `SequenceFile` object, `hmmsearch` will reopen the
    file in each thread, and load targets *iteratively* to scan with the
    query. Otherwise, it will *pre-fetch* the target sequences into a
    `DigitalSequenceBlock` collection, and share them across threads
    without copy. The *pre-fetching* gives much higher performance at the
    cost of extra  startup time and much higher memory consumption. You may
    want to check how much memory is available (for instance with
    `psutil.virtual_memory`) before trying to load a whole sequence database,
    but it is really recommended to do so whenever possible.

    Arguments:
        queries (iterable of `HMM`, `Profile` or `OptimizedProfile`): The
            query HMMs or profiles to search for in the database. Note that
            passing a single object is supported.
        sequences (iterable of `~pyhmmer.easel.DigitalSequence`): A
            database of sequences to query. If you plan on using the
            same sequences several times, consider storing them into
            a `~pyhmmer.easel.DigitalSequenceBlock` directly. If a
            `~pyhmmer.easel.SequenceFile` is given, profiles will be loaded
            iteratively from disk rather than prefetched.
        cpus (`int`): The number of threads to run in parallel. Pass ``1``
            to run everything in the main thread, ``0`` to automatically
            select a suitable number (using `psutil.cpu_count`), or any
            positive number otherwise.
        callback (callable): A callback that is called everytime a query is
            processed with two arguments: the query, and the total number
            of queries. This can be used to display progress in UI.
        backend (`str`): The parallel backend to use for workers to be
            executed. Supports ``threading`` to use thread-based parallelism,
            or ``multiprocessing`` to use process-based parallelism.
        parallel (`str`): The parallel strategy to use. Supports ``queries``
            to run queries in parallel, or ``targets`` to parallelize on
            targets while running one query at a time. If `None` given,
            use ``queries`` by default unless we can detect that there is
            a single or a small number of queries. Note that parallelization
            on ``targets`` does not work with `~pyhmmer.easel.SequenceFile`
            targets.

    Yields:
        `~pyhmmer.plan7.TopHits`: An object reporting *top hits* for each
        query, in the same order the queries were passed in the input.

    Raises:
        `~pyhmmer.errors.AlphabetMismatch`: When any of the query HMMs
            and the sequences do not share the same alphabet.
        `RuntimeError`: When attempting to use ``targets`` parallel
            strategy with targets from a `~pyhmmer.easel.SequenceFile`.

    Note:
        Any additional arguments passed to the `hmmsearch` function will be
        passed transparently to the `~pyhmmer.plan7.Pipeline` to be created.
        For instance, to run a ``hmmsearch`` using a bitscore cutoffs of
        5 instead of the default E-value cutoff, use::

            >>> hits = next(hmmsearch(thioesterase, proteins, T=5))
            >>> hits[0].score
            8.601...

        Since *version 0.11.0*, ``mypy`` should be able to detection which
        keywords can be passed to `hmmsearch` using a `TypedDict` annotation.

    .. versionadded:: 0.1.0

    .. versionadded:: 0.11.1
       ``parallel`` argument to select parallelization strategy.

    .. versionchanged:: 0.4.9
       Allow using `Profile` and `OptimizedProfile` queries.

    .. versionchanged:: 0.7.0
        Queries may now be an iterable of different types, or a single object.

    """
    cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1
    alphabet = options.get("alphabet")

    if not isinstance(queries, collections.abc.Iterable):
        queries = (queries,)

    if isinstance(sequences, SequenceFile):
        if sequences.name is None:
            raise ValueError("expected named `SequenceFile` for targets")
        if not sequences.digital:
            raise ValueError("expected digital mode `SequenceFile` for targets")
        assert sequences.alphabet is not None
        alphabet = alphabet or sequences.alphabet
        targets: typing.Union["SequenceFile[DigitalSequence]", DigitalSequenceBlock] = sequences
    elif isinstance(sequences, DigitalSequenceBlock):
        alphabet = alphabet or sequences.alphabet
        targets = sequences
    else:
        queries = peekable(queries)
        try:
            alphabet = alphabet or queries.peek().alphabet or Alphabet.amino()
            targets = DigitalSequenceBlock(alphabet, sequences)
        except StopIteration:
            alphabet = alphabet or Alphabet.amino()
            targets = DigitalSequenceBlock(alphabet)

    # attempt to optimize parallelism based on the number of queries --
    # for few queries it's probably more efficient to parallelize on targets
    # instead, but we need to have a DigitalSequenceBlock for that
    _queries_hint = operator.length_hint(queries)
    _few_queries = _queries_hint != 0 and _queries_hint < cpus
    if parallel is None:
        if _few_queries and isinstance(targets, DigitalSequenceBlock):
            parallel = "targets"
        else:
            parallel = "queries"
    if parallel == "targets" and not isinstance(targets, DigitalSequenceBlock):
        raise RuntimeError("cannot use ``targets`` parallel mode with a sequence file")

    # start the dispatcher
    if "alphabet" not in options:
        options["alphabet"] = alphabet
    dclass = _SEARCHDispatcher if parallel == "queries" else _ReverseSEARCHDispatcher
    dispatcher = dclass(
        queries=queries,
        targets=targets,
        cpus=cpus,
        backend=backend,
        callback=callback,  # type: ignore
        builder=None,
        **options,
    )
    return dispatcher.run()  # type: ignore
