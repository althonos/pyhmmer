import collections
import copy
import os
import queue
import multiprocessing
import typing
import threading

import psutil

from ..easel import Alphabet, DigitalSequence, DigitalMSA, DigitalSequenceBlock, SequenceFile
from ..plan7 import TopHits, HMM, Profile, OptimizedProfile, LongTargetsPipeline, Builder
from ..utils import singledispatchmethod, peekable
from ._base import _BaseDispatcher, _BaseWorker, _BaseChore, _AnyProfile

_NHMMERQueryType = typing.Union[DigitalSequence, DigitalMSA, _AnyProfile]
_N = typing.TypeVar("_N", DigitalSequence, DigitalMSA, HMM, Profile, OptimizedProfile)

if typing.TYPE_CHECKING:
    from ._base import Unpack, LongTargetsPipelineOptions, BACKEND

# --- Worker -------------------------------------------------------------------

class _NHMMERWorker(
    _BaseWorker[
        _NHMMERQueryType,
        typing.Union[DigitalSequenceBlock, "SequenceFile[DigitalSequence]"],
        "TopHits[_NHMMERQueryType]",
    ],
    threading.Thread
):

    pipeline_class = LongTargetsPipeline

    @singledispatchmethod
    def query(self, query) -> "TopHits[Any]":  # type: ignore
        raise TypeError(
            "Unsupported query type for `nhmmer`: {}".format(type(query).__name__)
        )

    @query.register(DigitalSequence)
    def _(self, query: DigitalSequence) -> "TopHits[DigitalSequence]":  # type: ignore
        assert self.pipeline is not None
        return self.pipeline.search_seq(query, self.targets, self.builder)

    @query.register(DigitalMSA)
    def _(self, query: DigitalMSA) -> "TopHits[DigitalMSA]":  # type: ignore
        assert self.pipeline is not None
        return self.pipeline.search_msa(query, self.targets, self.builder)

    @query.register(HMM)
    @query.register(Profile)
    @query.register(OptimizedProfile)
    def _(self, query: _AnyProfile) -> "TopHits[_AnyProfile]":  # type: ignore
        assert self.pipeline is not None
        return self.pipeline.search_hmm(query, self.targets)


class _NHMMERThread(_NHMMERWorker, threading.Thread):
    pass


class _NHMMERProcess(_NHMMERWorker, multiprocessing.Process):
    pass


# --- Dispatcher ---------------------------------------------------------------

class _NHMMERDispatcher(
    _BaseDispatcher[
        _NHMMERQueryType,
        typing.Union[DigitalSequenceBlock, "SequenceFile[DigitalSequence]"],
        "TopHits[_NHMMERQueryType]",
    ]
):
    def __init__(
        self,
        queries: typing.Iterable[_NHMMERQueryType],
        targets: typing.Union[DigitalSequenceBlock, "SequenceFile[DigitalSequence]"],
        cpus: int = 0,
        callback: typing.Optional[
            typing.Callable[[_NHMMERQueryType, int], None]
        ] = None,
        builder: typing.Optional[Builder] = None,
        timeout: int = 1,
        **options,  # type: Unpack[LongTargetsPipelineOptions]
    ) -> None:
        super().__init__(
            queries=queries,
            targets=targets,
            cpus=cpus,
            callback=callback,
            builder=builder,
            timeout=timeout,
            **options,  # type: ignore
        )

    def _new_worker(
        self,
        query_queue: "queue.Queue[typing.Optional[_BaseChore[_NHMMERQueryType, TopHits[_NHMMERQueryType]]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _NHMMERWorker:
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
            return _NHMMERThread(*params)
        elif self.backend == "multiprocessing":
            return _NHMMERProcess(*params)
        else:
            raise ValueError(f"Invalid backend for `nhmmer`: {self.backend!r}")


# --- nhmmer -----------------------------------------------------------------

def nhmmer(
    queries: typing.Union[_N, typing.Iterable[_N]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_N, int], None]] = None,
    backend: "BACKEND" = "threading",
    builder: typing.Optional[Builder] = None,
    **options,  # type: Unpack[LongTargetsPipelineOptions]
) -> typing.Iterator["TopHits[_N]"]:
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
        backend (`str`): The parallel backend to use for workers to be
            executed. Supports ``threading`` to use thread-based parallelism,
            or ``multiprocessing`` to use process-based parallelism.

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
        alphabet = alphabet or Alphabet.dna()
        targets = DigitalSequenceBlock(alphabet, sequences)

    if builder is None:
        builder = Builder(
            alphabet,
            seed=options.get("seed", 42),
            window_length=options.get("window_length"),
            window_beta=options.get("window_beta"),
        )

    if "alphabet" not in options:
        options["alphabet"] = alphabet
    dispatcher = _NHMMERDispatcher(
        queries=queries,
        targets=targets,
        cpus=cpus,
        backend=backend,
        callback=callback,  # type: ignore
        builder=builder,
        **options,
    )
    return dispatcher.run()  # type: ignore

