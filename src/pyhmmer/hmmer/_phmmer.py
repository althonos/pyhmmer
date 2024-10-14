import collections
import copy
import os
import queue
import multiprocessing
import typing
import threading

import psutil

from ..easel import Alphabet, DigitalSequence, DigitalMSA, DigitalSequenceBlock, SequenceFile
from ..plan7 import TopHits, Builder, Pipeline
from ..utils import singledispatchmethod, peekable
from ._base import _BaseDispatcher, _BaseWorker, _BaseChore

_PHMMERQueryType = typing.Union[DigitalSequence, DigitalMSA]
_M = typing.TypeVar("_M", DigitalSequence, DigitalMSA)

if typing.TYPE_CHECKING:
    from ._base import Unpack, PipelineOptions, BACKEND

# --- Worker -------------------------------------------------------------------

class _PHMMERWorker(
    _BaseWorker[
        _PHMMERQueryType,
        typing.Union[DigitalSequenceBlock, "SequenceFile[DigitalSequence]"],
        "TopHits[_PHMMERQueryType]",
    ],
    threading.Thread
):
    @singledispatchmethod
    def query(self, query) -> "TopHits[Any]":  # type: ignore
        raise TypeError(
            "Unsupported query type for `phmmer`: {}".format(type(query).__name__)
        )

    @query.register(DigitalSequence)
    def _(self, query: DigitalSequence) -> "TopHits[DigitalSequence]":  # type: ignore
        assert self.pipeline is not None
        return self.pipeline.search_seq(query, self.targets, self.builder)

    @query.register(DigitalMSA)
    def _(self, query: DigitalMSA) -> "TopHits[DigitalMSA]":  # type: ignore
        assert self.pipeline is not None
        return self.pipeline.search_msa(query, self.targets, self.builder)


class _PHMMERThread(_PHMMERWorker, threading.Thread):
    pass


class _PHMMERProcess(_PHMMERWorker, multiprocessing.Process):
    pass

# --- Dispatcher ---------------------------------------------------------------

class _PHMMERDispatcher(
    _BaseDispatcher[
        _PHMMERQueryType,
        typing.Union[DigitalSequenceBlock, "SequenceFile[DigitalSequence]"],
        "TopHits[_PHMMERQueryType]",
    ]
):
    def _new_worker(
        self,
        query_queue: "queue.Queue[typing.Optional[_BaseChore[_PHMMERQueryType, TopHits[_PHMMERQueryType]]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _PHMMERWorker:
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
            copy.copy(self.builder),
        ]
        if self.backend == "threading":
            return _PHMMERThread(*params)
        elif self.backend == "multiprocessing":
            return _PHMMERProcess(*params)
        else:
            raise ValueError(f"Invalid backend for `phmmer`: {self.backend!r}")

# --- phmmer -----------------------------------------------------------------

def phmmer(
    queries: typing.Union[_M, typing.Iterable[_M]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_M, int], None]] = None,
    backend: "BACKEND" = "threading",
    builder: typing.Optional[Builder] = None,
    **options,  # type: Unpack[PipelineOptions]
) -> typing.Iterator["TopHits[_M]"]:
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
        backend (`str`): The parallel backend to use for workers to be
            executed. Supports ``threading`` to use thread-based parallelism,
            or ``multiprocessing`` to use process-based parallelism.

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
        sequences = peekable(sequences)
        try:
            alphabet = alphabet or sequences.peek().alphabet or Alphabet.amino()
            targets = DigitalSequenceBlock(alphabet, sequences)
        except StopIteration:
            alphabet = alphabet or Alphabet.amino()
            targets = DigitalSequenceBlock(alphabet)

    if "alphabet" not in options:
        options["alphabet"] = alphabet
    dispatcher = _PHMMERDispatcher(
        queries=queries,
        targets=targets,
        cpus=cpus,
        backend=backend,
        callback=callback,  # type: ignore
        builder=builder,
        **options,
    )
    return dispatcher.run()  # type: ignore

