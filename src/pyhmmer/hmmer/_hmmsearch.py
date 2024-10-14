import collections
import queue
import multiprocessing
import typing
import os
import threading

import psutil

from ..easel import Alphabet, DigitalSequence, DigitalMSA, DigitalSequenceBlock, SequenceFile
from ..plan7 import TopHits, Builder, Pipeline, HMM, Profile, OptimizedProfile, HMMPressedFile, OptimizedProfileBlock
from ..utils import singledispatchmethod, peekable
from ._base import _BaseDispatcher, _BaseWorker, _BaseChore, _AnyProfile, _ProcessChore

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
        ]
        if self.backend == "threading":
            return _SEARCHThread(*params)
        elif self.backend == "multiprocessing":
            return _SEARCHProcess(*params)
        else:
            raise ValueError(f"Invalid backend for `hmmsearch`: {self.backend!r}")


# --- hmmsearch --------------------------------------------------------------

def hmmsearch(
    queries: typing.Union[_P, typing.Iterable[_P]],
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[_P, int], None]] = None,
    backend: "BACKEND" = "threading",
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
            `SequenceFile` is given, profiles will be loaded iteratively
            from disk rather than prefetched.
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

        Since *version 0.11.0*, ``mypy`` should be able to detection which
        keywords can be passed to `hmmsearch` using a `TypedDict` annotation.

    .. versionadded:: 0.1.0

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
        # sequence_file = sequences
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

    if "alphabet" not in options:
        options["alphabet"] = alphabet
    dispatcher = _SEARCHDispatcher(
        queries=queries,
        targets=targets,
        cpus=cpus,
        backend=backend,
        callback=callback,  # type: ignore
        builder=None,
        **options,
    )
    return dispatcher.run()  # type: ignore

