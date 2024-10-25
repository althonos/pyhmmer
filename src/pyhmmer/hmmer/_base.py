import abc
import ctypes
import collections
import contextlib
import itertools
import multiprocessing.synchronize
import multiprocessing.connection
import operator
import queue
import sys
import typing
import threading

from ..easel import DigitalSequenceBlock, SequenceFile, DigitalSequence, DigitalMSA
from ..plan7 import Pipeline, Builder, HMMPressedFile, HMM, Profile, OptimizedProfile, IterationResult, OptimizedProfileBlock
from ..utils import singledispatchmethod

# --- Type annotations ---------------------------------------------------------

# the query type for the pipeline
_Q = typing.TypeVar("_Q")
# the target type for the pipeline
_T = typing.TypeVar("_T")
# the result type for the pipeline
_R = typing.TypeVar("_R")

# type aliases
_AnyProfile = typing.Union[HMM, Profile, OptimizedProfile]

# `typing.Literal`` is only available in Python 3.8 and later
if typing.TYPE_CHECKING:

    if sys.version_info >= (3, 8):
        from typing import Literal, TypedDict
    else:
        from typing_extensions import Literal, TypedDict  # type: ignore

    if sys.version_info >= (3, 11):
        from typing import Unpack
    else:
        from typing_extensions import Unpack

    from .plan7 import BIT_CUTOFFS, STRAND, Background
    from .easel import Alphabet

    BACKEND = Literal["threading", "multiprocessing"]

    class PipelineOptions(TypedDict, total=False):
        alphabet: Alphabet
        background: typing.Optional[Background]
        bias_filter: bool
        null2: bool
        seed: int
        Z: float
        domZ: typing.Optional[float]
        F1: float
        F2: float
        F3: float
        E: float
        T: typing.Optional[float]
        domE: float
        domT: typing.Optional[float]
        incE: float
        incT: typing.Optional[float]
        incdomE: float
        incdomT: typing.Optional[float]
        bit_cutoffs: typing.Optional[BIT_CUTOFFS]

    class LongTargetsPipelineOptions(PipelineOptions, total=False):
        strand: typing.Optional[STRAND]
        B1: int
        B2: int
        B3: int
        block_length: int
        window_length: typing.Optional[int]
        window_beta: typing.Optional[float]


# --- Result class -------------------------------------------------------------

class _BaseChore(typing.Generic[_Q, _R], abc.ABC):
    query: _Q

    def __init__(self, query: _Q):
        self.query = query

    # --- Dispatcher methods ---

    @abc.abstractmethod
    def available(self) -> bool:
        """Return whether the chore is done and results are available."""

    @abc.abstractmethod
    def wait(self, timeout: typing.Optional[float] = None) -> bool:
        """Wait for the chore to be done."""

    @abc.abstractmethod
    def get(self) -> _R:
        """Get the results of the chore, blocking if the chore was not done."""

    # --- Worker methods ---

    @abc.abstractmethod
    def complete(self, result: _R) -> None:
        """Mark the chore as done and record ``result`` as the results."""
       
    @abc.abstractmethod
    def fail(self, exception: BaseException) -> None:
        """Mark the chore as done and record ``exception`` as the error."""


class _ThreadChore(typing.Generic[_Q, _R], _BaseChore[_Q, _R]):
    """A chore for a worker thread.

    Attributes:
        query (`object`): The query object to be processed by the worker
            thread. Exact type depends on the pipeline type.
        event (`threading.Event`): An event flag to set when the query
            is done being processed.
        result (`pyhmmer.plan7.TopHits` or `BaseException`): The results
            obtained after processing the query, or an exception on failure.

    """

    result: typing.Union[_R, BaseException, None]
    event: threading.Event

    def __init__(self, query: _Q) -> None:
        """Create a new chore from the given query."""
        super().__init__(query)
        self.event = threading.Event()
        self.result = None

    def available(self) -> bool:
        return self.event.is_set()

    def wait(self, timeout: typing.Optional[float] = None) -> bool:
        return self.event.wait(timeout)

    def get(self) -> _R:
        self.event.wait()
        assert self.result is not None
        if isinstance(self.result, BaseException):
            raise self.result
        else:
            return self.result

    def complete(self, result: _R) -> None:
        self.result = result
        self.event.set()

    def fail(self, exception: BaseException) -> None:
        self.result = exception
        self.event.set()


class _ProcessChore(typing.Generic[_Q, _R], _BaseChore[_Q, _R]):
    """A chore for a worker process.

    Attributes:
        query (`object`): The query object to be processed by the worker
            thread. Exact type depends on the pipeline type.
        conns (`multiprocessing.Connection`): The pipe end where the 
            worker should send the result.
        connr (`multiprocessing.Connection`): The pipe end where the 
            dispatcher should receive the result.

    """

    connr: multiprocessing.connection.Connection
    conns: multiprocessing.connection.Connection

    def __init__(self, query: _Q) -> None:
        """Create a new chore from the given query."""
        super().__init__(query)
        self.connr, self.conns = multiprocessing.Pipe()

    def available(self) -> bool:
        """Return whether the chore is done and results are available."""
        return self.connr.poll(0)

    def wait(self, timeout: typing.Optional[float] = None) -> bool:
        """Wait for the chore to be done."""
        return self.connr.poll(timeout)

    def complete(self, result: _R) -> None:
        """Mark the chore as done and record ``result`` as the results."""
        self.conns.send(result)

    def fail(self, exception: BaseException) -> None:
        """Mark the chore as done and record ``exception`` as the error."""
        self.conns.send(exception)

    def get(self) -> _R:
        """Get the results of the chore, blocking if the chore was not done."""
        self.wait()
        result = self.connr.recv()
        self.connr.close()
        self.conns.close()
        assert result is not None
        if isinstance(result, BaseException):
            raise result
        else:
            return result  # type: ignore


# --- Workers ------------------------------------------------------------------

class _BaseWorker(typing.Generic[_Q, _T, _R]):
    """A generic worker to parallelize a pipelined search.

    Attributes:
        targets (`DigitalSequenceBlock` or `OptimizedProfileBlock`): The
            target to search for hits, either a digital sequence block while
            in search mode, or an optimized profile block while in scan mode.
        query_queue (`queue.Queue`): The queue used to pass queries between 
            the dispatcher and the workers. It contains the query, its index 
            so that the results can be returned in the same order, and 
            a way to pass back the results corresponding to that query
            back to the dispatcher.
        query_count (`multiprocessing.Value`): An atomic counter storing
            the total number of queries that have been loaded so far.
            Passed to the ``callback`` so that an UI can show the total
            for a progress bar.
        kill_switch (`threading.Event`): An event flag shared between
            all worker threads, used to notify emergency exit.
        callback (`callable`, optional): An optional callback to be called
            after each query has been processed. It should accept two
            arguments: the query object that was processed, and the total
            number of queries read until now. For `multiprocessing` workers,
            this callback should be picklable.
        options (`dict`): A dictionary of options to be passed to the
            `pyhmmer.plan7.Pipeline` object wrapped by the worker.
        builder (`~pyhmmer.plan7.Builder`, *optional*): The builder to use
            for translating sequence or alignment queries into `HMM` objects.
            May be `None` if the queries are expected to be `HMM` only.

    """

    pipeline_class: typing.ClassVar[typing.Type[Pipeline]] = Pipeline

    @staticmethod
    def _none_callback(hmm: _Q, total: int) -> None:
        pass

    def __init__(
        self,
        targets: _T,
        query_queue: "queue.Queue[typing.Optional[_BaseChore[_Q, _R]]]",
        query_count: multiprocessing.Value,  # type: ignore
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[_Q, int], None]],
        options: "PipelineOptions",
        builder: typing.Optional[Builder] = None,
    ) -> None:
        super().__init__()
        self.options = options
        self.targets: _T = targets
        self.pipeline: typing.Optional[Pipeline] = None
        self.pipeline_options = options
        self.query_queue: "queue.Queue[typing.Optional[_BaseChore[_Q, _R]]]" = query_queue
        self.query_count = query_count
        self.callback: typing.Optional[typing.Callable[[_Q, int], None]] = (
            callback or self._none_callback
        )
        self.kill_switch = kill_switch
        self.builder = builder

    def run(self) -> None:
        while not self.is_killed():
            # attempt to get the next argument, with a timeout
            # so that the worker can periodically check if it has
            # been killed, even when no queries are available
            try:
                chore = self.query_queue.get(timeout=1)
            except queue.Empty:
                continue
            except BrokenPipeError:
                self.kill()
                break
            except ConnectionResetError:
                self.kill()
                break
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

    if typing.TYPE_CHECKING:
        def start(self) -> None: ...

    def is_killed(self) -> bool:
        try:
            return self.kill_switch.is_set()
        except BrokenPipeError:  # the connection was closed already
            return True
        except ConnectionResetError:  # the connection was closed already
            return True
        except FileNotFoundError:  # the Event manager has been closed already
            return True

    def kill(self) -> None:
        """Set the synchronized kill switch for all threads."""
        if not self.is_killed():
            self.kill_switch.set()

    def process(self, query: _Q) -> _R:
        """Process a single query and return the resulting hits."""
        if isinstance(self.targets, (HMMPressedFile, SequenceFile)):
            self.targets.rewind()
        if self.pipeline is None:
            self.pipeline = self.pipeline_class(**self.pipeline_options)
        hits = self.query(query)
        self.callback(query, self.query_count.value)  # type: ignore
        self.pipeline.clear()
        return hits

    @abc.abstractmethod
    def query(self, query: _Q) -> _R:
        """Run a single query against the target database."""
        return NotImplemented


# --- Dispatcher ---------------------------------------------------------------

class _BaseDispatcher(typing.Generic[_Q, _T, _R], abc.ABC):
    def __init__(
        self,
        queries: typing.Iterable[_Q],
        targets: _T,
        cpus: int = 0,
        callback: typing.Optional[typing.Callable[[_Q, int], None]] = None,
        builder: typing.Optional[Builder] = None,
        timeout: int = 1,
        backend: "BACKEND" = "threading",
        **options,  # type: Unpack[PipelineOptions]
    ) -> None:
        self.queries = queries
        self.targets: _T = targets
        self.callback: typing.Optional[typing.Callable[[_Q, int], None]] = callback
        self.options = options
        self.builder = builder
        self.timeout = timeout

        # make sure a positive number of CPUs is requested
        if cpus <= 0:
            raise ValueError("`cpus` must be strictly positive, not {!r}".format(cpus))

        # reduce the number of threads if there are less queries (at best
        # use one thread by query)
        hint = operator.length_hint(queries, -1)
        self.cpus = 1 if hint == 0 else min(cpus, hint) if hint > 0 else cpus

        # force "threading" backend when running everything in main thread
        self.backend = "threading" if self.cpus == 1 else backend

    @abc.abstractmethod
    def _new_worker(
        self,
        query_queue: "queue.Queue[typing.Optional[_BaseChore[_Q, _R]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _BaseWorker[_Q, _T, _R]:
        return NotImplemented

    def _new_chore(
        self,
        query: _Q,
    ) -> _BaseChore[_Q, _R]:
        if self.backend == "threading":
            return _ThreadChore(query)
        elif self.backend == "multiprocessing":
            return _ProcessChore(query)
        else:
            raise ValueError(f"Invalid parallel backend: {self.backend!r}")

    def _single_threaded(self) -> typing.Iterator[_R]:
        # create the queues to pass the HMM objects around, as well as atomic
        # values that we use to synchronize the threads
        query_queue = queue.Queue()  # type: ignore
        query_count = multiprocessing.Value(ctypes.c_ulong)
        kill_switch = threading.Event()

        # create the thread (to recycle code)
        thread = self._new_worker(query_queue, query_count, kill_switch)

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
        with contextlib.ExitStack() as ctx:

            # create the queues to pass the query objects around, as well as
            # atomic values that we use to synchronize the threads
            results: typing.Deque[_BaseChore[_Q, _R]] = collections.deque()
            if self.backend == "multiprocessing":
                manager = ctx.enter_context(multiprocessing.Manager())
                query_queue = manager.Queue(maxsize=2*self.cpus)
                query_count = manager.Value(ctypes.c_ulong, 0)
                kill_switch = manager.Event()
            elif self.backend == "threading":
                query_queue = queue.Queue(maxsize=2*self.cpus)
                query_count = multiprocessing.Value(ctypes.c_ulong)  # type: ignore
                kill_switch = threading.Event()

            # create and launch one pipeline thread per CPU
            workers = []
            for _ in range(self.cpus):
                worker = self._new_worker(query_queue, query_count, kill_switch)
                worker.start()
                workers.append(worker)

            # catch exceptions to kill threads in the background before exiting
            try:
                # alternate between feeding queries to the threads and
                # yielding back results, if available. the priority is
                # given to filling the query queue, so that no worker
                # ever idles.
                for query in self.queries:
                    # prepare to add the next query to the queue
                    query_count.value += 1
                    chore = self._new_chore(query)
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
                for _ in workers:
                    query_queue.put(None)
                # yield all remaining results, in order
                while results:
                    yield results[0].get()  # <-- blocks until result is available
                    results.popleft()
            except BaseException:
                # make sure threads are killed to avoid being stuck,
                # e.g. after a KeyboardInterrupt, then re-raise
                try:
                    kill_switch.set()
                except queue.Full:
                    pass
                raise

    def run(self) -> typing.Iterator[_R]:
        if self.cpus == 1:
            return self._single_threaded()
        else:
            return self._multi_threaded()
