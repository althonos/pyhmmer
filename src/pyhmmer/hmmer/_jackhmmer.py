import collections
import copy
import queue
import os
import itertools
import multiprocessing
import sys
import typing
import threading

import psutil

from ..easel import Alphabet, DigitalSequence, DigitalMSA, DigitalSequenceBlock, SequenceFile
from ..plan7 import TopHits, HMM, Profile, OptimizedProfile, Pipeline, Builder, IterationResult
from ..utils import singledispatchmethod, peekable
from ._base import _BaseDispatcher, _BaseWorker, _BaseChore, _AnyProfile

_JACKHMMERQueryType = typing.Union[DigitalSequence, _AnyProfile]
_I = typing.TypeVar("_I") # generic iteration result

if typing.TYPE_CHECKING:

    if sys.version_info >= (3, 8):
        from typing import Literal
    else:
        from typing_extensions import Literal  # type: ignore

    if sys.version_info >= (3, 11):
        from typing import Unpack
    else:
        from typing_extensions import Unpack

    from ._base import PipelineOptions, BACKEND

# --- Worker -------------------------------------------------------------------

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
        query_queue: "queue.Queue[typing.Optional[_BaseChore[_JACKHMMERQueryType, _I]]]",
        query_count: multiprocessing.Value,  # type: ignore
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]],
        options: "PipelineOptions",
        builder: typing.Optional[Builder] = None,
        max_iterations: typing.Optional[int] = 5,
        select_hits: typing.Optional[typing.Callable[["TopHits[_JACKHMMERQueryType]"], None]] = None,
        checkpoints: bool = False,
    ) -> None:
        super().__init__(
            targets=targets,
            query_queue=query_queue,
            query_count=query_count,
            kill_switch=kill_switch,
            callback=callback,
            options=options,
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
        assert self.pipeline is not None
        iterator = self.pipeline.iterate_seq(
            query,
            self.targets,
            self.builder,
            self.select_hits  # type: ignore
        )
        return self._iterate(iterator, self.checkpoints)

    @query.register(HMM)
    def _(self, query: HMM) -> typing.Union[IterationResult, typing.Iterable[IterationResult]]:  # type: ignore
        assert self.pipeline is not None
        iterator = self.pipeline.iterate_hmm(
            query,
            self.targets,
            self.builder,
            self.select_hits  # type: ignore
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


class _JACKHMMERThread(_JACKHMMERWorker, threading.Thread):
    pass


class _JACKHMMERProcess(_JACKHMMERWorker, multiprocessing.Process):
    pass


# --- Dispatcher ---------------------------------------------------------------

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
        builder: typing.Optional[Builder] = None,
        timeout: int = 1,
        max_iterations: typing.Optional[int] = 5,
        select_hits: typing.Optional[typing.Callable[["TopHits[_JACKHMMERQueryType]"], None]] = None,
        checkpoints: bool = False,
        **options,  # type: Unpack[PipelineOptions]
    ) -> None:
        super().__init__(
            queries=queries,
            targets=targets,
            cpus=cpus,
            callback=callback,
            builder=builder,
            timeout=timeout,
            **options,
        )
        self.max_iterations = max_iterations
        self.select_hits = select_hits
        self.checkpoints = checkpoints

    def _new_worker(
        self,
        query_queue: "queue.Queue[typing.Optional[_BaseChore[_JACKHMMERQueryType, _I]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _JACKHMMERWorker[_I]:
        params = [
            self.targets,
            query_queue,
            query_count,
            kill_switch,
            self.callback,
            self.options,
            copy.copy(self.builder),
            self.max_iterations,
            self.select_hits,
            self.checkpoints,
        ]
        if self.backend == "threading":
            return _JACKHMMERThread(*params)
        elif self.backend == "multiprocessing":
            return _JACKHMMERProcess(*params)
        else:
            raise ValueError(f"Invalid backend for `jackhmmer`: {self.backend!r}")


# --- jackhmmer -----------------------------------------------------------------

@typing.overload
def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[["TopHits[_JACKHMMERQueryType]"], None]] = None,
    checkpoints: "Literal[True]",
    cpus: int = 0,
    backend: "BACKEND" = "threading",
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: Unpack[PipelineOptions]
) -> typing.Iterator[typing.Iterable[IterationResult]]:
    ...


@typing.overload
def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[["TopHits[_JACKHMMERQueryType]"], None]] = None,
    checkpoints: "Literal[False]",
    cpus: int = 0,
    backend: "BACKEND" = "threading",
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: Unpack[PipelineOptions]
) -> typing.Iterator[IterationResult]:
    ...


@typing.overload
def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Union[typing.Iterable[DigitalSequence], "SequenceFile[DigitalSequence]"],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[["TopHits[_JACKHMMERQueryType]"], None]] = None,
    checkpoints: bool = False,
    cpus: int = 0,
    backend: "BACKEND" = "threading",
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: Unpack[PipelineOptions]
) -> typing.Union[
    typing.Iterator[IterationResult], typing.Iterator[typing.Iterable[IterationResult]]
]:
    ...


def jackhmmer(
    queries: typing.Union[_JACKHMMERQueryType, typing.Iterable[_JACKHMMERQueryType]],
    sequences: typing.Union[typing.Iterable[DigitalSequence], "SequenceFile[DigitalSequence]"],
    *,
    max_iterations: typing.Optional[int] = 5,
    select_hits: typing.Optional[typing.Callable[["TopHits[_JACKHMMERQueryType]"], None]] = None,
    checkpoints: bool = False,
    cpus: int = 0,
    backend: "BACKEND" = "threading",
    callback: typing.Optional[typing.Callable[[_JACKHMMERQueryType, int], None]] = None,
    builder: typing.Optional[Builder] = None,
    **options,  # type: Unpack[PipelineOptions]
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
        backend (`str`): The parallel backend to use for workers to be
            executed. Supports ``threading`` to use thread-based parallelism,
            or ``multiprocessing`` to use process-based parallelism.

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
    cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1
    options.setdefault("incE", 0.001)
    options.setdefault("incdomE", 0.001)
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
        targets = typing.cast(DigitalSequenceBlock, sequences.read_block())
    elif isinstance(sequences, DigitalSequenceBlock):
        alphabet = alphabet or sequences.alphabet
        targets = sequences
    else:
        sequences = peekable(sequences)
        try:
            alphabet = alphabet or sequences.peek().alphabet or Alphabet.amino()
            targets = DigitalSequenceBlock(alphabet, sequences)
        except StopIteration:
            alphabet = alphabet or Alphabet.amino()
            targets = DigitalSequenceBlock(alphabet)

    if builder is None:
        builder = Builder(
            alphabet,
            seed=options.get("seed", 42),
            architecture="hand",
        )

    if "alphabet" not in options:
        options["alphabet"] = alphabet
    dispatcher = _JACKHMMERDispatcher(  # type: ignore
        queries=queries,
        targets=targets,
        cpus=cpus,
        backend="threading",
        callback=callback,
        builder=builder,
        max_iterations=max_iterations,
        select_hits=select_hits,
        checkpoints=checkpoints,
        **options,
    )
    return dispatcher.run()

