import collections
import queue
import multiprocessing
import os
import typing
import threading

import psutil

from ..easel import Alphabet, DigitalSequence, DigitalMSA, DigitalSequenceBlock, SequenceFile
from ..plan7 import TopHits, Builder, Pipeline, Background, HMM, Profile, OptimizedProfile, HMMFile, HMMPressedFile, OptimizedProfileBlock
from ..utils import singledispatchmethod
from ._base import _BaseDispatcher, _BaseWorker, _BaseChore, _AnyProfile

if typing.TYPE_CHECKING:
    from ._base import Unpack, PipelineOptions, BACKEND

# --- Worker -------------------------------------------------------------------

class _SCANWorker(
    _BaseWorker[
        DigitalSequence,
        typing.Union[OptimizedProfileBlock, HMMPressedFile],
        "TopHits[DigitalSequence]",
    ],
    threading.Thread
):
    @singledispatchmethod
    def query(self, query) -> "TopHits[Any]":  # type: ignore
        raise TypeError(
            "Unsupported query type for `hmmscan`: {}".format(type(query).__name__)
        )

    @query.register(DigitalSequence)
    def _(self, query: DigitalSequence) -> "TopHits[DigitalSequence]":  # type: ignore
        assert self.pipeline is not None
        return self.pipeline.scan_seq(query, self.targets)


class _SCANThread(_SCANWorker, threading.Thread):
    pass


class _SCANProcess(_SCANWorker, multiprocessing.Process):
    pass

# --- Dispatcher ---------------------------------------------------------------

class _SCANDispatcher(
    _BaseDispatcher[
        DigitalSequence,
        typing.Union[OptimizedProfileBlock, HMMPressedFile],
        "TopHits[DigitalSequence]",
    ]
):
    def _new_worker(
        self,
        query_queue: "queue.Queue[typing.Optional[_BaseChore[DigitalSequence, TopHits[DigitalSequence]]]]",
        query_count: "multiprocessing.Value[int]",  # type: ignore
        kill_switch: threading.Event,
    ) -> _SCANWorker:
        if isinstance(self.targets, HMMPressedFile):
            assert self.targets.name is not None
            targets = HMMPressedFile(self.targets.name)
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
            return _SCANThread(*params)
        elif self.backend == "multiprocessing":
            return _SCANProcess(*params)
        else:
            raise ValueError(f"Invalid backend for `hmmsearch`: {self.backend!r}")


# --- hmmscan ----------------------------------------------------------------

def hmmscan(
    queries: typing.Union[DigitalSequence, typing.Iterable[DigitalSequence]],
    profiles: typing.Iterable[typing.Union[HMM, Profile, OptimizedProfile]],
    *,
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[DigitalSequence, int], None]] = None,
    backend: "BACKEND" = "threading",
    **options,  # type: Unpack[PipelineOptions]
) -> typing.Iterator["TopHits[DigitalSequence]"]:
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
        backend (`str`): The parallel backend to use for workers to be
            executed. Supports ``threading`` to use thread-based parallelism,
            or ``multiprocessing`` to use process-based parallelism.

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

            >>> with HMMFile("tests/data/hmms/db/PF02826.hmm") as hmm_file:
            ...     targets = hmm_file.optimized_profiles()
            ...     all_hits = list(hmmscan(proteins, targets, E=1e-10))
            >>> sum(len(hits) for hits in all_hits)
            6

        Otherwise, passing ``hmm_file`` as the ``profiles`` argument of
        `hmmscan` would cause the entire HMM file to be loaded in memory
        into an `OptimizedProfileBlock` otherwise.

    .. versionadded:: 0.7.0

    """
    cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1
    alphabet = options.get("alphabet")
    background = options.get("background")

    if not isinstance(queries, collections.abc.Iterable):
        queries = (queries,)
    if isinstance(profiles, HMMPressedFile):
        opt = profiles.read()
        profiles.rewind()
        if opt is not None:
            alphabet = alphabet or opt.alphabet
        else:
            alphabet = Alphabet.amino()
        targets = profiles
    elif isinstance(profiles, OptimizedProfileBlock):
        alphabet = alphabet or profiles.alphabet
        targets = profiles  # type: ignore
    else:
        block = None
        for item in profiles:
            alphabet = alphabet or item.alphabet
            if block is None:
                block = OptimizedProfileBlock(item.alphabet)
            if isinstance(item, HMM):
                if background is None:
                    background = Background(item.alphabet)
                profile = Profile(item.M, item.alphabet)
                profile.configure(item, background)
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
        if alphabet is None:
            alphabet = Alphabet.amino()
        if block is None:
            block = OptimizedProfileBlock(alphabet)
        targets = block  # type: ignore

    if "alphabet" not in options:
        options["alphabet"] = alphabet
    if "background" not in options and background is not None:
        options["background"] = background
    dispatcher = _SCANDispatcher(
        queries=queries,
        targets=targets,
        cpus=cpus,
        backend=backend,
        callback=callback,
        builder=None,
        **options,
    )
    return dispatcher.run()

