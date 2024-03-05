# coding: utf-8
"""Reimplementation of HMMER binaries with the PyHMMER API.

Note:
    Functions of this module handle parallelization using threads to run
    searches in parallel for the different queries. If less queries are
    given, the number of threads will be reduced to avoid spawning idle
    threads.

"""

import abc
import contextlib
import collections
import copy
import ctypes
import itertools
import io
import multiprocessing
import os
import operator
import queue
import threading
import time
import typing

import psutil

from .easel import (
    Alphabet,
    DigitalSequence,
    DigitalMSA,
    MSA,
    MSAFile,
    TextSequence,
    SequenceFile,
    SSIWriter,
    DigitalSequenceBlock,
)
from .plan7 import (
    Builder,
    Background,
    Pipeline,
    LongTargetsPipeline,
    TopHits,
    IterationResult,
    HMM,
    HMMFile,
    HMMPressedFile,
    Profile,
    TraceAligner,
    OptimizedProfile,
    OptimizedProfileBlock,
)
from .utils import peekable, singledispatchmethod

# the query type for the pipeline
_Q = typing.TypeVar("_Q")
# the target type for the pipeline
_T = typing.TypeVar("_T")
# the result type for the pipeline
_R = typing.TypeVar("_R")
_I = typing.TypeVar("_I")

# the query types for the different tasks
_PHMMERQueryType = typing.Union[DigitalSequence, DigitalMSA]
_SEARCHQueryType = typing.Union[HMM, Profile, OptimizedProfile]
_NHMMERQueryType = typing.Union[_PHMMERQueryType, _SEARCHQueryType]
_JACKHMMERQueryType = typing.Union[DigitalSequence, _SEARCHQueryType]

# `typing.Literal`` is only available in Python 3.8 and later
if typing.TYPE_CHECKING:
    try:
        from typing import Literal
    except ImportError:
        from typing_extensions import Literal  # type: ignore

# --- Result class -----------------------------------------------------------


class _Chore(typing.Generic[_Q, _R]):
    """A chore for a worker thread.

    Attributes:
        query (`object`): The query object to be processed by the worker
            thread. Exact type depends on the pipeline type.
        event (`threading.Event`): An event flag to set when the query
            is done being processed.
        hits (`pyhmmer.plan7.TopHits`): The hits obtained after processing
            the query.
        exception (`BaseException`): An exception that occured while
            processing the query.

    """

    query: _Q
    event: threading.Event
    result: typing.Optional[_R]
    exception: typing.Optional[BaseException]

    __slots__ = ("query", "event", "result", "exception")

    def __init__(self, query: _Q) -> None:
        """Create a new chore from the given query."""
        self.query = query
        self.event = threading.Event()
        self.result = None
        self.exception = None

    def available(self) -> bool:
        """Return whether the chore is done and results are available."""
        return self.event.is_set()

    def wait(self, timeout: typing.Optional[float] = None) -> bool:
        """Wait for the chore to be done."""
        return self.event.wait(timeout)

    def get(self) -> _R:
        """Get the results of the chore, blocking if the chore was not done."""
        self.event.wait()
        if self.exception is not None:
            raise self.exception
        return typing.cast(_R, self.result)

    def complete(self, result: _R) -> None:
        """Mark the chore as done and record ``result`` as the results."""
        self.result = result
        self.event.set()

    def fail(self, exception: BaseException) -> None:
        """Mark the chore as done and record ``exception`` as the error."""
        self.exception = exception
        self.event.set()


# --- Pipeline threads -------------------------------------------------------


class _BaseWorker(typing.Generic[_Q, _T, _R], threading.Thread):
    """A generic worker thread to parallelize a pipelined search.

    Attributes:
        targets (`DigitalSequenceBlock` or `OptimizedProfileBlock`): The
            target to search for hits, either a digital sequence block while
            in search mode, or an optimized profile block while in scan mode.
        query_queue (`queue.Queue`): The queue used to pass queries
            between threads. It contains the query, its index so that the
            results can be returned in the same order, and a `_ResultBuffer`
            where to store the result when the query has been processed.
        query_count (`multiprocessing.Value`): An atomic counter storing
            the total number of queries that have currently been loaded.
            Passed to the ``callback`` so that an UI can show the total
            for a progress bar.
        kill_switch (`threading.Event`): An event flag shared between
            all worker threads, used to notify emergency exit.
        callback (`callable`, optional): An optional callback to be called
            after each query has been processed. It should accept two
            arguments: the query object that was processed, and the total
            number of queries read until now.
        options (`dict`): A dictionary of options to be passed to the
            `pyhmmer.plan7.Pipeline` object wrapped by the worker thread.
        pipeline_class (`type`): The pipeline class to use to search for
            hits. Use `~plan7.LongTargetsPipeline` for `nhmmer`, and
            `~plan7.Pipeline` everywhere else.
        builder (`~pyhmmer.plan7.Builder`, *optional*): The builder to use
            for translating sequence or alignment queries into `HMM` objects.
            May be `None` if the queries are expected to be `HMM` only.

    """

    @staticmethod
    def _none_callback(hmm: _Q, total: int) -> None:
        pass

    def __init__(
        self,
        targets: _T,
        query_queue: "queue.Queue[typing.Optional[_Chore[_Q, _R]]]",
        query_count: multiprocessing.Value,  # type: ignore
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[_Q, int], None]],
        options: typing.Dict[str, typing.Any],
        pipeline_class: typing.Type[Pipeline],
        alphabet: Alphabet,
        builder: typing.Optional[Builder] = None,
    ) -> None:
        super().__init__()
        self.options = options
        self.targets: _T = targets
        self.pipeline = pipeline_class(alphabet=alphabet, **options)
        self.query_queue: "queue.Queue[typing.Optional[_Chore[_Q, _R]]]" = query_queue
        self.query_count = query_count
        self.callback: typing.Optional[typing.Callable[[_Q, int], None]] = (
            callback or self._none_callback
        )
        self.kill_switch = kill_switch
        self.builder = builder

    def run(self) -> None:
        while not self.kill_switch.is_set():
            # attempt to get the next argument, with a timeout
            # so that the thread can periodically check if it has
            # been killed, even when no queries are available
            try:
                chore = self.query_queue.get(timeout=1)
            except queue.Empty:
                continue
            # check if arguments from the queue are a poison-pill (`None`),
            # in which case the thread will stop running
            if chore is None:
                break
            # process the query, making sure to capture any exception
            # and then mark the hits as "found" using a `threading.Event`
            try:
                hits = self.process(chore.query)
                chore.complete(hits)
            except BaseException as exc:
                self.kill()
                chore.fail(exc)
        if isinstance(self.targets, (SequenceFile, HMMPressedFile)):
            self.targets.close()

    def kill(self) -> None:
        """Set the synchronized kill switch for all threads."""
        self.kill_switch.set()

    def process(self, query: _Q) -> _R:
        """Process a single query and return the resulting hits."""
        if isinstance(self.targets, (HMMPressedFile, SequenceFile)):
            self.targets.rewind()
        hits = self.query(query)
        self.callback(query, self.query_count.value)  # type: ignore
        self.pipeline.clear()
        return hits

    @abc.abstractmethod
    def query(self, query: _Q) -> _R:
        """Run a single query against the target database."""
        return NotImplemented


class _SEARCHWorker(
    _BaseWorker[
        _SEARCHQueryType,
        typing.Union[DigitalSequenceBlock, SequenceFile],
        TopHits,
    ]
):
    @singledispatchmethod
    def query(self, query) -> TopHits:  # type: ignore
        raise TypeError(
            "Unsupported query type for `hmmsearch`: {}".format(type(query).__name__)
        )

    @query.register(HMM)
    @query.register(Profile)
    @query.register(OptimizedProfile)
    def _(self, query: typing.Union[HMM, Profile, OptimizedProfile]) -> TopHits:  # type: ignore
        return self.pipeline.search_hmm(query, self.targets)


class _PHMMERWorker(
    _BaseWorker[
        _PHMMERQueryType,
        typing.Union[DigitalSequenceBlock, SequenceFile],
        TopHits,
    ]
):
    @singledispatchmethod
    def query(self, query) -> TopHits:  # type: ignore
        raise TypeError(
            "Unsupported query type for `phmmer`: {}".format(type(query).__name__)
        )

    @query.register(DigitalSequence)
    def _(self, query: DigitalSequence) -> TopHits:  # type: ignore
        return self.pipeline.search_seq(query, self.targets, self.builder)

    @query.register(DigitalMSA)
    def _(self, query: DigitalMSA) -> TopHits:  # type: ignore
        return self.pipeline.search_msa(query, self.targets, self.builder)


class _JACKHMMERWorker(
    typing.Generic[_I],
    _BaseWorker[
        _JACKHMMERQueryType,
        DigitalSequenceBlock,
        _I,
    ],
):
    def __init__(
        self,
        targets: DigitalSequenceBlock,
        query_queue: "queue.Queue[typing.Optional[_Chore[_JACKHMMERQueryType, _I]]]",
        query_count: multiprocessing.Value,  # type: ignore
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]],
        options: typing.Dict[str, typing.Any],
        pipeline_class: typing.Type[Pipeline],
        alphabet: Alphabet,
        builder: typing.Optional[Builder] = None,
        max_iterations: typing.Optional[int] = 5,
        select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
        checkpoints: bool = False,
    ) -> None:
        super().__init__(
            targets=targets,
            query_queue=query_queue,
            query_count=query_count,
            kill_switch=kill_switch,
            callback=callback,
            options=options,
            pipeline_class=pipeline_class,
            alphabet=alphabet,
            builder=builder,
        )
        self.select_hits = select_hits
        self.max_iterations = max_iterations
        self.checkpoints = checkpoints

    @singledispatchmethod
    def query(self, query) -> typing.Union[IterationResult, typing.Iterable[IterationResult]]:  # type: ignore
        raise TypeError(
            "Unsupported query type for `jackhmmer`: {}".format(type(query).__name__)
        )

    @query.register(DigitalSequence)
    def _(self, query: DigitalSequence) -> typing.Union[IterationResult, typing.Iterable[IterationResult]]:  # type: ignore
        iterator = self.pipeline.iterate_seq(
            query, self.targets, self.builder, self.select_hits
        )
        return self._iterate(iterator, self.checkpoints)

    @query.register(HMM)
    def _(self, query: HMM) -> typing.Union[IterationResult, typing.Iterable[IterationResult]]:  # type: ignore
        iterator = self.pipeline.iterate_hmm(
            query, self.targets, self.builder, self.select_hits
        )
        return self._iterate(iterator, self.checkpoints)

    @typing.overload
    def _iterate(
        self, iterator: typing.Iterable[IterationResult], checkpoints: "Literal[False]"
    ) -> IterationResult:
        ...

    @typing.overload
    def _iterate(
        self, iterator: typing.Iterable[IterationResult], checkpoints: "Literal[True]"
    ) -> typing.Iterable[IterationResult]:
        ...

    @typing.overload
    def _iterate(
        self, iterator: typing.Iterable[IterationResult], checkpoints: bool = False
    ) -> typing.Union[IterationResult, typing.Iterable[IterationResult]]:
        ...

    def _iterate(
        self, iterator: typing.Iterable[IterationResult], checkpoints: bool = False
    ) -> typing.Union[IterationResult, typing.Iterable[IterationResult]]:
        iteration_checkpoints = []
        for iteration in itertools.islice(iterator, self.max_iterations):
            if checkpoints:
                iteration_checkpoints.append(iteration)
            if iteration.converged:
                break
        return iteration_checkpoints if checkpoints else iteration


class _NHMMERWorker(
    _BaseWorker[
        _NHMMERQueryType,
        typing.Union[DigitalSequenceBlock, SequenceFile],
        TopHits,
    ]
):
    @singledispatchmethod
    def query(self, query) -> TopHits:  # type: ignore
        raise TypeError(
            "Unsupported query type for `nhmmer`: {}".format(type(query).__name__)
        )

    @query.register(DigitalSequence)
    def _(self, query: DigitalSequence) -> TopHits:  # type: ignore
        return self.pipeline.search_seq(query, self.targets, self.builder)

    @query.register(DigitalMSA)
    def _(self, query: DigitalMSA) -> TopHits:  # type: ignore
        return self.pipeline.search_msa(query, self.targets, self.builder)

    @query.register(HMM)
    @query.register(Profile)
    @query.register(OptimizedProfile)
    def _(self, query: typing.Union[HMM, Profile, OptimizedProfile]) -> TopHits:  # type: ignore
        return self.pipeline.search_hmm(query, self.targets)


class _SCANWorker(
    _BaseWorker[
        DigitalSequence,
        typing.Union[OptimizedProfileBlock, HMMPressedFile],
        TopHits,
    ]
):
    @singledispatchmethod
    def query(self, query) -> TopHits:  # type: ignore
        raise TypeError(
            "Unsupported query type for `hmmscan`: {}".format(type(query).__name__)
        )

    @query.register(DigitalSequence)
    def _(self, query: DigitalSequence) -> TopHits:  # type: ignore
        return self.pipeline.scan_seq(query, self.targets)


# --- Search runners ---------------------------------------------------------


class _BaseDispatcher(typing.Generic[_Q, _T, _R], abc.ABC):
    def __init__(
        self,
        queries: typing.Iterable[_Q],
        targets: _T,
        cpus: int = 0,
        callback: typing.Optional[typing.Callable[[_Q, int], None]] = None,
        pipeline_class: typing.Type[Pipeline] = Pipeline,
        alphabet: Alphabet = Alphabet.amino(),
        builder: typing.Optional[Builder] = None,
        timeout: int = 1,
        **options,  # type: object
    ) -> None:
        self.queries = queries
        self.targets: _T = targets
        self.callback: typing.Optional[typing.Callable[[_Q, int], None]] = callback
        self.options = options
        self.pipeline_class = pipeline_class
        self.alphabet = alphabet
        self.builder = builder
        self.timeout = timeout

        # make sure a positive number of CPUs is requested
        if cpus <= 0:
            raise ValueError("`cpus` must be strictly positive, not {!r}".format(cpus))

        # reduce the number of threads if there are less queries (at best
        # use one thread by query)
        hint = operator.length_hint(queries, -1)
        self.cpus = 1 if hint == 0 else min(cpus, hint) if hint > 0 else cpus

    @abc.abstractmethod
    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[_Chore[_Q, _R]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _BaseWorker[_Q, _T, _R]:
        return NotImplemented

    def _single_threaded(self) -> typing.Iterator[_R]:
        # create the queues to pass the HMM objects around, as well as atomic
        # values that we use to synchronize the threads
        query_queue = queue.Queue()  # type: ignore
        query_count = multiprocessing.Value(ctypes.c_ulong)
        kill_switch = threading.Event()

        # create the thread (to recycle code)
        thread = self._new_thread(query_queue, query_count, kill_switch)

        # process each HMM iteratively and yield the result
        # immediately so that the user can iterate over the
        # TopHits one at a time
        for query in self.queries:
            query_count.value += 1
            yield thread.process(query)

        # close the targets if they were coming from a file
        if isinstance(thread.targets, (SequenceFile, HMMPressedFile)):
            thread.targets.close()

    def _multi_threaded(self) -> typing.Iterator[_R]:
        # create the queues to pass the query objects around, as well as
        # atomic values that we use to synchronize the threads
        results: typing.Deque[_Chore[_Q, _R]] = collections.deque()
        query_queue = queue.Queue(maxsize=2*self.cpus)  # type: ignore
        query_count = multiprocessing.Value(ctypes.c_ulong)
        kill_switch = threading.Event()

        # create and launch one pipeline thread per CPU
        threads = []
        for _ in range(self.cpus):
            thread = self._new_thread(query_queue, query_count, kill_switch)
            thread.start()
            threads.append(thread)

        # catch exceptions to kill threads in the background before exiting
        try:
            # alternate between feeding queries to the threads and
            # yielding back results, if available. the priority is
            # given to filling the query queue, so that no worker
            # ever idles.
            for query in self.queries:
                # prepare to add the next query to the queue
                query_count.value += 1
                chore: _Chore[_Q, _R] = _Chore(query)
                # attempt to add the query to the query queue,
                # but use a timeout so that it doesn't deadlock
                # if the background thread errored (which would
                # set the kill switch), typically because the
                # user-provided callback failed
                while not kill_switch.is_set():
                    with contextlib.suppress(queue.Full):
                        query_queue.put(chore, timeout=self.timeout)  # <-- blocks if too many chores in queue
                        results.append(chore)
                        break
                # aggressively wait for the result with a very short
                # timeout, and exit the loop if the queue is not full
                if results[0].available():
                    yield results[0].get()
                    results.popleft()
            # now that we exhausted all queries, poison pill the
            # threads so they stop on their own gracefully
            for _ in threads:
                query_queue.put(None)
            # yield all remaining results, in order
            while results:
                yield results[0].get()  # <-- blocks until result is available
                results.popleft()
        except BaseException:
            # make sure threads are killed to avoid being stuck,
            # e.g. after a KeyboardInterrupt, then re-raise
            kill_switch.set()
            raise

    def run(self) -> typing.Iterator[_R]:
        if self.cpus == 1:
            return self._single_threaded()
        else:
            return self._multi_threaded()


class _SEARCHDispatcher(
    _BaseDispatcher[
        _SEARCHQueryType,
        typing.Union[DigitalSequenceBlock, SequenceFile],
        TopHits,
    ]
):
    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[_Chore[_SEARCHQueryType, TopHits]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _SEARCHWorker:
        if isinstance(self.targets, SequenceFile):
            assert self.targets.name is not None
            targets = SequenceFile(
                self.targets.name,
                format=self.targets.format,
                digital=True,
                alphabet=self.alphabet,
            )
        else:
            targets = self.targets  # type: ignore
        return _SEARCHWorker(
            targets,
            query_queue,
            query_count,
            kill_switch,
            self.callback,
            self.options,
            self.pipeline_class,
            self.alphabet,
        )


class _PHMMERDispatcher(
    _BaseDispatcher[
        _PHMMERQueryType,
        typing.Union[DigitalSequenceBlock, SequenceFile],
        TopHits,
    ]
):
    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[_Chore[_PHMMERQueryType, TopHits]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _PHMMERWorker:
        if isinstance(self.targets, SequenceFile):
            assert self.targets.name is not None
            targets = SequenceFile(
                self.targets.name,
                format=self.targets.format,
                digital=True,
                alphabet=self.alphabet,
            )
        else:
            targets = self.targets  # type: ignore
        return _PHMMERWorker(
            targets,
            query_queue,
            query_count,
            kill_switch,
            self.callback,
            self.options,
            self.pipeline_class,
            self.alphabet,
            copy.copy(self.builder),
        )


class _JACKHMMERDispatcher(
    typing.Generic[_I],
    _BaseDispatcher[
        _JACKHMMERQueryType,
        DigitalSequenceBlock,
        _I,
    ],
):
    """A dispatcher to run JackHMMER iterative searches."""

    def __init__(
        self,
        queries: typing.Iterable[_JACKHMMERQueryType],
        targets: DigitalSequenceBlock,
        cpus: int = 0,
        callback: typing.Optional[
            typing.Callable[[_JACKHMMERQueryType, int], None]
        ] = None,
        pipeline_class: typing.Type[Pipeline] = Pipeline,
        alphabet: Alphabet = Alphabet.amino(),
        builder: typing.Optional[Builder] = None,
        timeout: int = 1,
        max_iterations: typing.Optional[int] = 5,
        select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
        checkpoints: bool = False,
        **options,  # type: object
    ) -> None:
        super().__init__(
            queries=queries,
            targets=targets,
            cpus=cpus,
            callback=callback,
            pipeline_class=pipeline_class,
            alphabet=alphabet,
            builder=builder,
            timeout=timeout,
            **options,
        )
        self.max_iterations = max_iterations
        self.select_hits = select_hits
        self.checkpoints = checkpoints

    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[_Chore[_JACKHMMERQueryType, _I]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _JACKHMMERWorker[_I]:
        return _JACKHMMERWorker(
            self.targets,
            query_queue,
            query_count,
            kill_switch,
            self.callback,
            self.options,
            self.pipeline_class,
            self.alphabet,
            copy.copy(self.builder),
            self.max_iterations,
            self.select_hits,
            self.checkpoints,
        )


class _NHMMERDispatcher(
    _BaseDispatcher[
        _NHMMERQueryType,
        typing.Union[DigitalSequenceBlock, SequenceFile],
        TopHits,
    ]
):
    def __init__(
        self,
        queries: typing.Iterable[_NHMMERQueryType],
        targets: typing.Union[DigitalSequenceBlock, SequenceFile],
        cpus: int = 0,
        callback: typing.Optional[
            typing.Callable[[_NHMMERQueryType, int], None]
        ] = None,
        pipeline_class: typing.Type[Pipeline] = LongTargetsPipeline,
        alphabet: Alphabet = Alphabet.dna(),
        builder: typing.Optional[Builder] = None,
        timeout: int = 1,
        **options,  # type: typing.Dict[str, object]
    ) -> None:
        super().__init__(
            queries=queries,
            targets=targets,
            cpus=cpus,
            callback=callback,
            pipeline_class=pipeline_class,
            alphabet=alphabet,
            builder=builder,
            timeout=timeout,
            **options,
        )

    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[_Chore[_NHMMERQueryType, TopHits]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _NHMMERWorker:
        if isinstance(self.targets, SequenceFile):
            assert self.targets.name is not None
            targets = SequenceFile(
                self.targets.name,
                format=self.targets.format,
                digital=True,
                alphabet=self.alphabet,
            )
        else:
            targets = self.targets  # type: ignore
        return _NHMMERWorker(
            targets,
            query_queue,
            query_count,
            kill_switch,
            self.callback,
            self.options,
            self.pipeline_class,
            self.alphabet,
            copy.copy(self.builder),
        )


class _SCANDispatcher(
    _BaseDispatcher[
        DigitalSequence,
        typing.Union[OptimizedProfileBlock, HMMPressedFile],
        TopHits,
    ]
):
    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[_Chore[DigitalSequence, TopHits]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _SCANWorker:
        if isinstance(self.targets, HMMPressedFile):
            assert self.targets.name is not None
            targets = HMMPressedFile(self.targets.name)
        else:
            targets = self.targets  # type: ignore
        return _SCANWorker(
            targets,
            query_queue,
            query_count,
            kill_switch,
            self.callback,
            self.options,
            self.pipeline_class,
            self.alphabet,
        )


# --- hmmsearch --------------------------------------------------------------


def hmmsearch(
    queries: typing.Union[_SEARCHQueryType, typing.Iterable[_SEARCHQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_SEARCHQueryType, int], None]] = None,
    **options,  # type: typing.Dict[str, object]
) -> typing.Iterator[TopHits]:
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
            `SequenceFile` is given, profiles will be loaded iteratively
            from disk rather than prefetched.
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
        For instance, to run a ``hmmsearch`` using a bitscore cutoffs of
        5 instead of the default E-value cutoff, use::

            >>> hits = next(hmmsearch(thioesterase, proteins, T=5))
            >>> hits[0].score
            8.601...

    .. versionadded:: 0.1.0

    .. versionchanged:: 0.4.9
       Allow using `Profile` and `OptimizedProfile` queries.

    .. versionchanged:: 0.7.0
        Queries may now be an iterable of different types, or a single object.

    """
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1

    if not isinstance(queries, collections.abc.Iterable):
        queries = (queries,)

    if isinstance(sequences, SequenceFile):
        sequence_file: SequenceFile = sequences
        if sequence_file.name is None:
            raise ValueError("expected named `SequenceFile` for targets")
        if not sequence_file.digital:
            raise ValueError("expected digital mode `SequenceFile` for targets")
        assert sequence_file.alphabet is not None
        alphabet = sequence_file.alphabet
        targets: typing.Union[SequenceFile, DigitalSequenceBlock] = sequence_file
    elif isinstance(sequences, DigitalSequenceBlock):
        alphabet = sequences.alphabet
        targets = sequences
    else:
        queries = peekable(queries)
        try:
            alphabet = queries.peek().alphabet
            targets = DigitalSequenceBlock(alphabet, sequences)
        except StopIteration:
            alphabet = Alphabet.amino()
            targets = DigitalSequenceBlock(alphabet)

    dispatcher = _SEARCHDispatcher(
        queries=queries,
        targets=targets,
        cpus=_cpus,
        callback=callback,
        alphabet=alphabet,
        builder=None,
        pipeline_class=Pipeline,
        **options,
    )
    return dispatcher.run()


# --- phmmer -----------------------------------------------------------------


def phmmer(
    queries: typing.Union[_PHMMERQueryType, typing.Iterable[_PHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_PHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: typing.Dict[str, object]
) -> typing.Iterator[TopHits]:
    """Search protein sequences against a sequence database.

    Arguments:
        queries (iterable of `DigitalSequence` or `DigitalMSA`): The query
            sequences to search for in the sequence database. Passing a
            single object is supported.
        sequences (iterable of `~pyhmmer.easel.DigitalSequence`): A database
            of sequences to query. If you plan on using the same sequences
            several times, consider storing them into a
            `~pyhmmer.easel.DigitalSequenceBlock` directly. If a
            `SequenceFile` is given, profiles will be loaded iteratively
            from disk rather than prefetched.
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

    Raises:
        `~pyhmmer.errors.AlphabetMismatch`: When any of the query sequence
            the profile or the optional builder do not share the same
            alphabet.

    Note:
        Any additional keyword arguments passed to the `phmmer` function
        will be passed transparently to the `~pyhmmer.plan7.Pipeline` to
        be created in each worker thread.

    .. versionadded:: 0.2.0

    .. versionchanged:: 0.3.0
       Allow using `DigitalMSA` queries.

    .. versionchanged:: 0.7.0
        Queries may now be an iterable of different types, or a single object.

    """
    _alphabet = Alphabet.amino()
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1
    _builder = Builder(_alphabet) if builder is None else builder

    if not isinstance(queries, collections.abc.Iterable):
        queries = (queries,)

    if isinstance(sequences, SequenceFile):
        sequence_file: SequenceFile = sequences
        if sequence_file.name is None:
            raise ValueError("expected named `SequenceFile` for targets")
        if not sequence_file.digital:
            raise ValueError("expected digital mode `SequenceFile` for targets")
        assert sequence_file.alphabet is not None
        alphabet = sequence_file.alphabet
        targets: typing.Union[SequenceFile, DigitalSequenceBlock] = sequence_file
    elif isinstance(sequences, DigitalSequenceBlock):
        alphabet = sequences.alphabet
        targets = sequences
    else:
        alphabet = _alphabet
        targets = DigitalSequenceBlock(_alphabet, sequences)

    dispatcher = _PHMMERDispatcher(
        queries=queries,
        targets=targets,
        cpus=_cpus,
        callback=callback,
        pipeline_class=Pipeline,
        alphabet=alphabet,
        builder=_builder,
        **options,
    )
    return dispatcher.run()


# --- jackhmmer -----------------------------------------------------------------


@typing.overload
def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
    checkpoints: "Literal[True]",
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: typing.Dict[str, object]
) -> typing.Iterator[typing.Iterable[IterationResult]]:
    ...


@typing.overload
def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
    checkpoints: "Literal[False]",
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: typing.Dict[str, object]
) -> typing.Iterator[IterationResult]:
    ...


@typing.overload
def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
    checkpoints: bool = False,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: typing.Dict[str, object]
) -> typing.Union[
    typing.Iterator[IterationResult], typing.Iterator[typing.Iterable[IterationResult]]
]:
    ...


def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
    checkpoints: bool = False,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: object
) -> typing.Union[
    typing.Iterator[IterationResult], typing.Iterator[typing.Iterable[IterationResult]]
]:
    """Search protein sequences against a sequence database.

    Arguments:
        queries (iterable of `DigitalSequence`): The query sequences to search
            for in the sequence database. Passing a single sequence object
            is supported.
        sequences (iterable of `~pyhmmer.easel.DigitalSequence`): A database
            of sequences to query. If you plan on using the same sequences
            several times, consider storing them into a
            `~pyhmmer.easel.DigitalSequenceBlock` directly. `jackhmmer` does
            not support passing a `~pyhmmer.easel.SequenceFile` at the
            moment.
        max_iterations (`int`): The maximum number of iterations for the
            search. Hits will be returned early if the searched converged.
        select_hits (callable, optional): A function or callable object
            for manually selecting hits during each iteration. It should
            take a single `~pyhmmer.plan7.TopHits` argument and change the
            inclusion of individual hits with the `~pyhmmer.plan7.Hit.include`
            and `~pyhmmer.plan7.Hit.drop` methods of `~pyhmmer.plan7.Hit`
            objects.
        checkpoints (`bool`): A logical flag to return the results at each
            iteration 'checkpoint'. If `True`, then an iterable of up to
            ``max_iterations`` `~pyhmmer.plan7.IterationResult` will be
            returned, rather than just the final iteration. This is similar
            to ``--chkhmm`` amd ``--chkali`` flags from HMMER3's
            ``jackhmmer`` interface.
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
        `~pyhmmer.plan7.IterationResult`: An *iteration result* instance for
        each query, in the same order the queries were passed in the input.
        If ``checkpoint`` option is `True`, all iterations will be returned
        instead of the last one.

    Raises:
        `~pyhmmer.errors.AlphabetMismatch`: When any of the query sequence
            the profile or the optional builder do not share the same
            alphabet.

    Note:
        Any additional keyword arguments passed to the `jackhmmer` function
        will be passed transparently to the `~pyhmmer.plan7.Pipeline` to
        be created in each worker thread.

    Caution:
        Default values used for ``jackhmmer`` do not correspond to the
        default parameters used for creating a pipeline in the other cases.
        If no parameter value is given as a keyword argument, `jackhmmer`
        will create the pipeline with ``incE=0.001`` and ``incdomE=0.001``,
        where a default `~pyhmmer.plan7.Pipeline` would use ``incE=0.01``
        and ``incdomE=0.01``.

    .. versionadded:: 0.8.0

    """
    _alphabet = Alphabet.amino()
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1
    _builder = Builder(_alphabet, architecture="hand") if builder is None else builder

    options.setdefault("incE", 0.001)
    options.setdefault("incdomE", 0.001)

    if not isinstance(queries, collections.abc.Iterable):
        queries = (queries,)

    if isinstance(sequences, SequenceFile):
        sequence_file: SequenceFile = sequences
        if sequence_file.name is None:
            raise ValueError("expected named `SequenceFile` for targets")
        if not sequence_file.digital:
            raise ValueError("expected digital mode `SequenceFile` for targets")
        assert sequence_file.alphabet is not None
        alphabet = sequence_file.alphabet
        targets = typing.cast(DigitalSequenceBlock, sequence_file.read_block())
    elif isinstance(sequences, DigitalSequenceBlock):
        alphabet = sequences.alphabet
        targets = sequences
    else:
        alphabet = _alphabet
        targets = DigitalSequenceBlock(_alphabet, sequences)

    dispatcher = _JACKHMMERDispatcher(  # type: ignore
        queries=queries,
        targets=targets,
        cpus=_cpus,
        callback=callback,
        pipeline_class=Pipeline,
        alphabet=alphabet,
        builder=_builder,
        max_iterations=max_iterations,
        select_hits=select_hits,
        checkpoints=checkpoints,
        **options,
    )
    return dispatcher.run()


# --- nhmmer -----------------------------------------------------------------


def nhmmer(
    queries: typing.Union[_NHMMERQueryType, typing.Iterable[_NHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_NHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: typing.Dict[str, object]
) -> typing.Iterator[TopHits]:
    """Search nucleotide sequences against a sequence database.

    Arguments:
        queries (iterable of `DigitalSequence`, `DigitalMSA`, `HMM`): The
            query sequences or profiles to search for in the sequence
            database. Passing a single object is supported.
        sequences (iterable of `~pyhmmer.easel.DigitalSequence`): A
            database of sequences to query. If you plan on using the
            same sequences several times, consider storing them into
            a `~pyhmmer.easel.DigitalSequenceBlock` directly. If a
            `SequenceFile` is given, profiles will be loaded iteratively
            from disk rather than prefetched.
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
        Any additional keyword arguments passed to the `nhmmer` function
        will be passed to the `~pyhmmer.plan7.LongTargetsPipeline` created
        in each worker thread. The ``strand`` argument can be used to
        restrict the search on the direct or reverse strand.

    Caution:
        This function is not just `phmmer` for nucleotide sequences; it
        actually uses a `~pyhmmer.plan7.LongTargetsPipeline` internally
        instead of processing each target sequence in its entirety when
        searching for hits. This avoids hitting the maximum target size
        that can be used (100,000 residues), which may be a problem for
        some larger genomes.

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.4.9
       Allow using `Profile` and `OptimizedProfile` queries.

    .. versionchanged:: 0.7.0
        Queries may now be an iterable of different types, or a single object.

    """
    _alphabet = Alphabet.dna()
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1

    if builder is None:
        _builder = Builder(
            _alphabet,
            seed=options.get("seed", 42),
            window_length=options.get("window_length"),
            window_beta=options.get("window_beta"),
        )
    else:
        _builder = builder

    if not isinstance(queries, collections.abc.Iterable):
        queries = (queries,)

    if isinstance(sequences, SequenceFile):
        sequence_file: SequenceFile = sequences
        if sequence_file.name is None:
            raise ValueError("expected named `SequenceFile` for targets")
        if not sequence_file.digital:
            raise ValueError("expected digital mode `SequenceFile` for targets")
        assert sequence_file.alphabet is not None
        alphabet = sequence_file.alphabet
        targets: typing.Union[SequenceFile, DigitalSequenceBlock] = sequence_file
    elif isinstance(sequences, DigitalSequenceBlock):
        alphabet = sequences.alphabet
        targets = sequences
    else:
        alphabet = _alphabet
        targets = DigitalSequenceBlock(_alphabet, sequences)

    dispatcher = _NHMMERDispatcher(
        queries=queries,
        targets=targets,
        cpus=_cpus,
        callback=callback,
        pipeline_class=LongTargetsPipeline,
        alphabet=alphabet,
        builder=_builder,
        **options,
    )
    return dispatcher.run()


# --- hmmpress ---------------------------------------------------------------


def hmmpress(
    hmms: typing.Iterable[HMM],
    output: typing.Union[str, "os.PathLike[str]"],
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
            om = gm.to_optimized()

            # update the disk offsets of the optimized model to be written
            om.offsets.model = h3m.tell()
            om.offsets.profile = h3p.tell()
            om.offsets.filter = h3f.tell()

            # check that hmm has a name
            if hmm.name is None:
                raise ValueError("HMMs must have a name to be pressed.")
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


# --- hmmalign ---------------------------------------------------------------


def hmmalign(
    hmm: HMM,
    sequences: typing.Iterable[DigitalSequence],
    *,
    digitize: bool = False,
    trim: bool = False,
    all_consensus_cols: bool = True,
) -> MSA:
    """Align several sequences to a reference HMM, and return the MSA.

    Arguments:
        hmm (`~pyhmmer.plan7.HMM`): The reference HMM to use for the
            alignment.
        sequences (iterable of `~pyhmmer.easel.DigitalSequence`): The
            sequences to align to the HMM. If you plan on using the
            same sequences several times, consider storing them into
            a `~pyhmmer.easel.DigitalSequenceBlock` directly.
        trim (`bool`): Trim off any residues that get assigned to
            flanking :math:`N` and :math:`C` states (in profile traces)
            or :math:`I_0` and :math:`I_m` (in core traces).
        digitize (`bool`): If set to `True`, returns a `DigitalMSA`
            instead of a `TextMSA`.
        all_consensus_cols (`bool`): Force a column to be created for
            every consensus column in the model, even if it means having
            all gap character in a column.

    Returns:
        `~pyhmmer.easel.MSA`: A multiple sequence alignment containing
        the aligned sequences, either a `TextMSA` or a `DigitalMSA`
        depending on the value of the ``digitize`` argument.

    See Also:
        The `~pyhmmer.plan7.TraceAligner` class, which lets you inspect the
        intermediate tracebacks obtained for each alignment before building
        a MSA.

    .. versionadded:: 0.4.7

    """
    aligner = TraceAligner()
    if not isinstance(sequences, DigitalSequenceBlock):
        sequences = DigitalSequenceBlock(hmm.alphabet, sequences)
    traces = aligner.compute_traces(hmm, sequences)
    return aligner.align_traces(
        hmm,
        sequences,
        traces,
        trim=trim,
        digitize=digitize,
        all_consensus_cols=all_consensus_cols,
    )


# --- hmmscan ----------------------------------------------------------------


def hmmscan(
    queries: typing.Union[DigitalSequence, typing.Iterable[DigitalSequence]],
    profiles: typing.Iterable[typing.Union[HMM, Profile, OptimizedProfile]],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[DigitalSequence, int], None]] = None,
    background: typing.Optional[Background] = None,
    **options,  # type: typing.Dict[str, object]
) -> typing.Iterator[TopHits]:
    """Scan query sequences against a profile database.

    In HMMER many-to-many comparisons, a *scan* is the operation of querying
    with sequences a database of profile HMMs. It is necessary slower than
    a *search* because reconfiguring profiles between each queries has
    additional overhead, so it's recommended to use a *search* if the order
    of the comparisons is not important.

    The `hmmscan` function offers two ways of managing the database that will
    be selected based on the type of the ``profiles`` argument. If
    ``profiles`` is an `HMMPressedFile` object, `hmmscan` will reopen the
    file in each thread, and load profiles *iteratively* to scan with the
    query. Otherwise, it will *pre-fetch* the optimized profiles into an
    `OptimizedProfileBlock` collection, and share them across queries. The
    *pre-fetching* gives much higher performance at the cost of extra
    startup time and much higher memory consumption. You may want to check
    how much memory is available (for instance with `psutil.virtual_memory`)
    before trying to load a whole pHMM database.

    Arguments:
        queries (iterable of `DigitalSequence`): The query sequences to scan
            with the database. Passing a single query is supported.
        profiles (iterable of `HMM`, `Profile` or `OptimizedProfile`): A
            database of profiles to query. If you plan on using the
            same targets several times, consider converting them into
            `OptimizedProfile` and storing them into an `OptimizedProfileBlock`
            ahead of time. If a `HMMPressedFile` is given, profiles will be
            loaded iteratively from disk rather than prefetched.
        cpus (`int`): The number of threads to run in parallel. Pass ``1``
            to run everything in the main thread, ``0`` to automatically
            select a suitable number (using `psutil.cpu_count`), or any
            positive number otherwise.
        callback (callable): A callback that is called everytime a query is
            processed with two arguments: the query, and the total number
            of queries. This can be used to display progress in UI.
        background (`pyhmmer.plan7.Background`, *optional*): A background
            object to use for configuring the profiles. If `None` given,
            create a default one.

    Yields:
        `~pyhmmer.plan7.TopHits`: An object reporting *top hits* for each
        query, in the same order the queries were passed in the input.

    Raises:
        `~pyhmmer.errors.AlphabetMismatch`: When any of the query sequence
            and the profile do not share the same alphabet.

    Note:
        Any additional keyword arguments passed to the `phmmer` function
        will be passed transparently to the `~pyhmmer.plan7.Pipeline` to
        be created in each worker thread.

    Hint:
        If reading the profiles from a pressed HMM database, make sure to
        use the `HMMFile.optimized_profiles` method so that profiles are
        read iteratively from the file during the scan loop::

            >>> with HMMFile("tests/data/hmms/db/t2pks.hmm") as hmm_file:
            ...     targets = hmm_file.optimized_profiles()
            ...     all_hits = list(hmmscan(proteins, targets, E=1e-10))
            >>> sum(len(hits) for hits in all_hits)
            26

        Otherwise, passing ``hmm_file`` as the ``profiles`` argument of
        `hmmscan` would cause the entire HMM file to be loaded in memory
        into an `OptimizedProfileBlock` otherwise.

    .. versionadded:: 0.7.0

    """
    _alphabet = Alphabet.amino()
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1
    _background = Background(_alphabet) if background is None else background
    options.setdefault("background", _background)  # type: ignore

    if not isinstance(queries, collections.abc.Iterable):
        queries = (queries,)
    if isinstance(profiles, HMMPressedFile):
        alphabet = _alphabet  # FIXME: try to detect from content instead?
        targets = profiles
    elif isinstance(profiles, OptimizedProfileBlock):
        alphabet = profiles.alphabet
        targets = profiles  # type: ignore
    else:
        alphabet = _alphabet
        block = OptimizedProfileBlock(_alphabet)
        for item in profiles:
            if isinstance(item, HMM):
                profile = Profile(item.M, item.alphabet)
                profile.configure(item, _background)
                item = profile
            if isinstance(item, Profile):
                item = item.to_optimized()
            if isinstance(item, OptimizedProfile):
                block.append(item)
            else:
                ty = type(item).__name__
                raise TypeError(
                    "Expected HMM, Profile or OptimizedProfile, found {}".format(ty)
                )
        targets = block  # type: ignore

    dispatcher = _SCANDispatcher(
        queries=queries,
        targets=targets,
        cpus=_cpus,
        callback=callback,
        pipeline_class=Pipeline,
        alphabet=alphabet,
        builder=None,
        **options,
    )
    return dispatcher.run()


# add a very limited CLI so that this module can be invoked in a shell:
#     $ python -m pyhmmer.hmmsearch <hmmfile> <seqdb>
if __name__ == "__main__":
    import argparse
    import sys

    # don't load target databases in memory if they would take more than
    # 80% of the remainining available memory (because more memory will be
    # needed afterwards to allocate the `TopHits` for each query)
    MAX_MEMORY_LOAD = 0.80

    def _hmmsearch(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.seqdb).st_size

        with SequenceFile(args.seqdb, digital=True) as sequences:
            # pre-load the database if it is small enough to fit in memory
            if database_size < available_memory * MAX_MEMORY_LOAD:
                sequences = sequences.read_block()  # type: ignore
            # load the query HMMs iteratively
            with HMMFile(args.hmmfile) as hmms:
                if hmms.is_pressed():
                    hmms = hmms.optimized_profiles()  # type: ignore
                hits_list = hmmsearch(hmms, sequences, cpus=args.jobs)  # type: ignore
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                (hit.accession or b"-").decode(),
                                (hits.query_name or b"-").decode(),
                                (hits.query_accession or b"-").decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    def _phmmer(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.seqdb).st_size

        alphabet = Alphabet.amino()
        with SequenceFile(args.seqdb, digital=True, alphabet=alphabet) as sequences:
            # pre-load the database if it is small enough to fit in memory
            if database_size < available_memory * MAX_MEMORY_LOAD:
                sequences = sequences.read_block()  # type: ignore
            # load the query sequences iteratively
            with SequenceFile(args.seqfile, digital=True, alphabet=alphabet) as queries:
                hits_list = phmmer(queries, sequences, cpus=args.jobs)  # type: ignore
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                "-",
                                hit.best_domain.alignment.hmm_accession.decode(),
                                hit.best_domain.alignment.hmm_name.decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    @contextlib.contextmanager
    def open_query_file(
        queryfile: typing.Union["os.PathLike[str]", typing.BinaryIO],
        alphabet: Alphabet,
    ) -> typing.Iterator[typing.Union[SequenceFile, HMMFile]]:
        """Open either a sequence file or an HMM file."""
        try:
            yield SequenceFile(queryfile, digital=True, alphabet=alphabet)
        except ValueError:
            yield HMMFile(queryfile)

    def _jackhmmer(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.seqdb).st_size

        alphabet = Alphabet.amino()
        with SequenceFile(args.seqdb, digital=True, alphabet=alphabet) as sequences:
            # pre-load the database if it is small enough to fit in memory
            if database_size < available_memory * MAX_MEMORY_LOAD:
                sequences = sequences.read_block()  # type: ignore
            # load the query sequences or HMMs iteratively
            with open_query_file(args.queryfile, alphabet) as queries:
                result = jackhmmer(queries, sequences, checkpoint=False, cpus=args.jobs)
                for hits in result.hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                "-",
                                hit.best_domain.alignment.hmm_accession.decode(),
                                hit.best_domain.alignment.hmm_name.decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    def _nhmmer(args: argparse.Namespace) -> int:
        # at the moment `LongTargetsPipeline` only support block targets, not files
        with SequenceFile(args.seqdb, digital=True) as seqfile:
            with SequenceFile(args.seqfile, digital=True) as queryfile:
                hits_list = nhmmer(queryfile, seqfile, cpus=args.jobs)  # type: ignore
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                "-",
                                hit.best_domain.alignment.hmm_accession.decode(),
                                hit.best_domain.alignment.hmm_name.decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    def _hmmscan(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.hmmdb).st_size

        with SequenceFile(args.seqfile, digital=True) as seqfile:
            with HMMFile(args.hmmdb) as hmms:
                # pre-load profiles is they can fit into memory
                targets = hmms.optimized_profiles() if hmms.is_pressed() else hmms
                if (
                    hmms.is_pressed()
                    and database_size < available_memory * MAX_MEMORY_LOAD
                ):
                    targets = OptimizedProfileBlock(seqfile.alphabet, targets)  # type: ignore
                # load the query sequences iteratively
                hits_list = hmmscan(seqfile, targets, cpus=args.jobs)  # type: ignore
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                (hit.accession or b"-").decode(),
                                (hits.query_name or b"-").decode(),
                                (hits.query_accession or b"-").decode(),
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

    def _hmmalign(args: argparse.Namespace) -> int:
        try:
            with SequenceFile(args.seqfile, args.informat, digital=True) as seqfile:
                sequences: typing.List[DigitalSequence] = list(seqfile)  # type: ignore
        except EOFError as err:
            print(err, file=sys.stderr)
            return 1

        with HMMFile(args.hmmfile) as hmms:
            hmm = next(hmms)
            if next(hmms, None) is not None:
                print("HMM file contains more than one HMM, exiting", file=sys.stderr)
                return 1

        msa = hmmalign(hmm, sequences, trim=args.trim)
        if args.output == "-":
            with io.BytesIO() as out:
                msa.write(out, args.outformat)
                print(out.getvalue().decode("ascii"), end="")
        else:
            with open(args.output, "wb") as out:
                msa.write(out, args.outformat)

        return 0

    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", required=False, default=0, type=int)
    subparsers = parser.add_subparsers(
        dest="cmd", help="HMMER command to run", required=True
    )

    parser_hmmsearch = subparsers.add_parser("hmmsearch")
    parser_hmmsearch.set_defaults(call=_hmmsearch)
    parser_hmmsearch.add_argument("hmmfile")
    parser_hmmsearch.add_argument("seqdb")

    parser_phmmer = subparsers.add_parser("phmmer")
    parser_phmmer.set_defaults(call=_phmmer)
    parser_phmmer.add_argument("seqfile")
    parser_phmmer.add_argument("seqdb")

    parser_phmmer = subparsers.add_parser("jackhmmer")
    parser_phmmer.set_defaults(call=_jackhmmer)
    parser_phmmer.add_argument("queryfile")  # can be sequences or HMM
    parser_phmmer.add_argument("seqdb")

    parser_nhmmer = subparsers.add_parser("nhmmer")
    parser_nhmmer.set_defaults(call=_nhmmer)
    parser_nhmmer.add_argument("seqfile")
    parser_nhmmer.add_argument("seqdb")

    parser_hmmsearch = subparsers.add_parser("hmmscan")
    parser_hmmsearch.set_defaults(call=_hmmscan)
    parser_hmmsearch.add_argument("hmmdb")
    parser_hmmsearch.add_argument("seqfile")

    parser_hmmpress = subparsers.add_parser("hmmpress")
    parser_hmmpress.set_defaults(call=_hmmpress)
    parser_hmmpress.add_argument("hmmfile")
    parser_hmmpress.add_argument("-f", "--force", action="store_true")

    parser_hmmalign = subparsers.add_parser("hmmalign")
    parser_hmmalign.set_defaults(call=_hmmalign)
    parser_hmmalign.add_argument("hmmfile", metavar="<hmmfile>")
    parser_hmmalign.add_argument(
        "seqfile",
        metavar="<seqfile>",
    )
    parser_hmmalign.add_argument(
        "-o",
        "--output",
        action="store",
        default="-",
        metavar="<f>",
        help="output alignment to file <f>, not stdout",
    )
    parser_hmmalign.add_argument(
        "--trim",
        action="store_true",
        help="trim terminal tails of nonaligned residues from alignment",
    )
    parser_hmmalign.add_argument(
        "--informat",
        action="store",
        metavar="<s>",
        help="assert <seqfile> is in format <s> (no autodetection)",
        choices=SequenceFile._FORMATS.keys(),
    )
    parser_hmmalign.add_argument(
        "--outformat",
        action="store",
        metavar="<s>",
        help="output alignment in format <s>",
        default="stockholm",
        choices=MSAFile._FORMATS.keys(),
    )

    args = parser.parse_args()
    sys.exit(args.call(args))
