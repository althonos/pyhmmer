# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyhmmer/compare/v0.10.9...HEAD


## [v0.10.9] - 2024-03-12
[v0.10.9]: https://github.com/althonos/pyhmmer/compare/v0.10.8...v0.10.9

### Fixed
- Reallocation issue causing segmentation faults in `nhmmer` with more than 64 sequences ([#62](https://github.com/althonos/pyhmmer/issues/62)).


## [v0.10.8] - 2024-03-06
[v0.10.8]: https://github.com/althonos/pyhmmer/compare/v0.10.7...v0.10.8

### Added
- Getter to access the strand of a `Domain` produced by a `LongTargetsPipeline`.

### Changed
- Display model and cutoff names in `MissingCutoffs` error message, if any.
- Allow `LongTargetsPipeline` to be configured with window length and beta parameters.
- Make `nhmmer` use the window length and beta from the options when creating a `Builder`.

### Fixed
- `nhmmer` not computing E-values for non-default window lengths ([moshi4/pybarrnap#2](https://github.com/moshi4/pybarrnap/issues/2)).
- `SequenceFile` and `MSAFile` crashing with a segmentation fault when given the path to a folder rather than a file.


## [v0.10.7] - 2024-03-04
[v0.10.7]: https://github.com/althonos/pyhmmer/compare/v0.10.6...v0.10.7

### Added
- Pre-compiled wheels for PyPy 3.10.

### Fixed
- Invalid pointer cast in `__getbuffer__` method of `Matrix` and `Vector` objects.
- Remaining tests failing to run on missing `importlib-resources`.
- `pyhmmer.hmmer` dispatchers possibly dead-locking on background thread errors ([#60](https://github.com/althonos/pyhmmer/issues/60)).


## [v0.10.6] - 2024-02-20
[v0.10.6]: https://github.com/althonos/pyhmmer/compare/v0.10.5...v0.10.6

### Added
-  `armv7` and `aarch64` to the `PKGBUILD` architectures.

### Changed
- `SSIReader` and `SSIWriter` constructors now accept path-like objects.
- Skip tests dependending on `importlib.resources.files` when it is not available on the host machine.

### Fixed
- Memory leak caused by alphabet allocation in `Pipeline._scan_loop_file`.


## [v0.10.5] - 2024-02-16
[v0.10.5]: https://github.com/althonos/pyhmmer/compare/v0.10.4...v0.10.5

### Added
- `Alignment` properties to get the original lengths of the sequence and HMM being stored.
- `Hit.length` property storing the length of the hit sequence (or HMM).
- `TopHits.query_length` storing the length of the hit HMM (or query).
- `Alignment.posterior_probabilities` property showing an encoded representation of posteriors ([#59](https://github.com/althonos/pyhmmer/pull/59), by [@arajkovic](https://github.com/arajkovic)).
- `Trace.score` method to compute a trace score from a given profile and sequence.
- `Alignment.__sizeof__` implementation leveraing `p7_alidisplay_SizeOf`. 

### Fixed
- `Cutoffs` proxy objects not recording their owner to prevent deallocation.
- Avoid GIL re-acquisition in `GeneticCode.translate`.
- Query metadata not being recorded in `Hits` obtained from `daemon.Client`.
- Empty `MatrixU8` creation attempting zero-allocation.
- `VectorU8.zeros` allocating 4x more memory than required.
- Memory leak caused by string duplication in `__getbuffer__` methods of `Matrix` and `Vector` types.


## [v0.10.4] - 2023-10-29
[v0.10.4]: https://github.com/althonos/pyhmmer/compare/v0.10.3...v0.10.4

### Added
- `residue_markups` argument to `TextSequence` and `DigitalSequence` constructors.
- `__reduce__` implementation to `TextSequence`, `DigitalSequence`, `TextSequenceBlock` and `DigitalSequenceBlock`.

### Changed
- Handling of `easel` I/O methods to avoid implicit GIL acquisition for error checking.

### Fixed
- Syntax errors in type annotation files.


## [v0.10.3] - 2023-10-22
[v0.10.3]: https://github.com/althonos/pyhmmer/compare/v0.10.2...v0.10.3

### Added
- Out-of-band pickle serialization of `Bitfield` objects.
- Getters for `float` attributes and forward/backward parameters of `OptimizedProfile`.
- `InvalidHMM` error raised by `HMM.validate`.

### Changed
- Mark `HMM.zero` method as `noexcept`.
- Increase size of buffer for the query queue in the `hmmer` dispatcher.

### Fixed
- Unneeded semaphore in `pyhmmer.hmmer` message passing implementation.
- Broken assertion in `Bitfield._from_raw_bytes`.
- Relax tolerance of HMM validation in `TraceAligner.align_traces`.


## [v0.10.2] - 2023-08-20
[v0.10.2]: https://github.com/althonos/pyhmmer/compare/v0.10.1...v0.10.2

### Fixed
- Invalid buffer write in `DigitalSequenceBlock.translate` ([#50](https://github.com/althonos/pyhmmer/issues/50)).


## [v0.10.1] - 2023-08-17
[v0.10.1]: https://github.com/althonos/pyhmmer/compare/v0.10.0...v0.10.1

### Added
- `HMM.set_consensus` method to set the consensus for a method or compute it from the emission probabilities.

### Fixed
- Platform detection for MacOS and Armv7 platforms in `setup.py`.
- `pyhmmer.plan7.HMM` constructor setting a consensus string forcefully.


## [v0.10.0] - 2023-08-16
[v0.10.0]: https://github.com/althonos/pyhmmer/compare/v0.9.0...v0.10.0

### Added
- Support for compiling wheels for Aarch64 and NEON-enabled Arm platforms.

### Changed
- Updated HMMER to [`v3.4`](https://github.com/EddyRivasLab/hmmer/releases/tag/hmmer-3.4).
- Updated Easel to [`v0.49`](https://github.com/EddyRivasLab/easel/releases/tag/easel-0.49).
- Use [`cibuildwheel`](https://github.com/pypa/cibuildwheel) to build wheel distributions.

### Fixed
- Patch missing `PyInterpreterState_GetID` preventing the package from working on PyPy 3.9.


## [v0.9.0] - 2023-08-03
[v0.9.0]: https://github.com/althonos/pyhmmer/compare/v0.8.2...v0.9.0

### Added
- `TopHits.mode` property showing from which pipeline mode (*search* or *scan*) the hits were obtained.

### Changed
- Updated the code for Cython `v3.0`.

### Fixed
- `TopHits.merge` not properly handling inclusion and reporting for domains ([#46](https://github.com/althonos/pyhmmer/issues/46), [#47](https://github.com/althonos/pyhmmer/pull/47), by [@zdk123](https://github.com/zdk123)).


## [v0.8.2] - 2023-06-07
[v0.8.2]: https://github.com/althonos/pyhmmer/compare/v0.8.1...v0.8.2

### Added
- Bracket-style `repr` implementation to `HMM`, `Profile` and `OptimizedProfile` showing model alphabet, length and name.
- `MissingCutoffs` and `InvalidParameter` exceptions inheriting `ValueError`.

### Changed
- Replace `pthread` locks with `PyThread` API for synchronizing models in `OptimizedProfileBlock`.

### Fixed
- Sequence length extraction in `LongTargetsPipeline.search_hmm` ([#42](https://github.com/althonos/pyhmmer/issues/42)).
- `LongTargetsPipeline.search_msa` not building a HMM with `Builder.build_msa`.


## [v0.8.1] - 2023-05-19
[v0.8.1]: https://github.com/althonos/pyhmmer/compare/v0.8.0...v0.8.1

### Added
- `HMM.validate` method to ensure a HMM holds HMMER structural constraints.
- `plan7.Transitions` enum with transition names for indexing `HMM.transition_probabilities`.


## [v0.8.0] - 2023-05-01
[v0.8.0]: https://github.com/althonos/pyhmmer/compare/v0.7.4...v0.8.0

PyHMMER has been accepted for publication in [Bioinformatics](https://academic.oup.com/bioinformatics). Paper can be reached at [doi:10.1093/bioinformatics/btad214](https://doi.org/10.1093/bioinformatics/btad214).

### Added
- `pyhmmer.hmmer.jackhmmer` function to run several JackHMMER iterative searches in parallel using multithreading ([#35](https://github.com/althonos/pyhmmer/pull/35), by [@zdk123](https://github.com/zdk123)).
- `HMM.to_profile` shortcut method to allocate and configure a new `Profile` object.

### Fixed
- Type annotations of `Pipeline.iterate_seq` and `Pipeline.iterate_hmm`.
- Potential memory leak on exceptions raised by `HMMPressedFile.read`.
- `Offsets.profile` not recording offsets properly, causing `pyhmmer.hmmer.hmmpress` to produce invalid pressed files ([#37](https://github.com/althonos/pyhmmer/issues/36)).

### Changed
- `HMM.__init__` and `HMM.sample` now take the `Alphabet` as the first argument, for consistency with the rest of the API.
- `HMM` now require a `name` argument.

### Removed
- Deprecated `ignore_gaps` argument in `SequenceFile.__init__`.
- Deprecated `Sequence.taxonomy_id` property.


## [v0.7.4] - 2023-04-14
[v0.7.4]: https://github.com/althonos/pyhmmer/compare/v0.7.3...v0.7.4

### Added
- *Recipes* page to the documentation with code example for loading multiple HMM files ([#24](https://github.com/althonos/pyhmmer/issues/24), by [@zdk123](https://github.com/zdk123)).

### Fixed
- `TraceAligner` methods causing a segfault when passed an uninitialized HMM ([#36](https://github.com/althonos/pyhmmer/issues/36)).

### Changed
- `HMM` default constructor now always creates a valid HMM (with respects to probability arrays).
- `TraceAligner` now validates the input `HMM` before calling the HMMER code.
- Use stack allocation for all error buffers instead of creating empty `bytearray` objects where applicable.


## [v0.7.3] - 2023-03-24
[v0.7.3]: https://github.com/althonos/pyhmmer/compare/v0.7.2...v0.7.3

### Fixed
- Wrong argument type in `IterativeSearch.iterate_hmm` method ([#34](https://github.com/althonos/pyhmmer/pull/34), by [@zdk123](https://github.com/zdk123)).


## [v0.7.2] - 2023-02-17
[v0.7.2]: https://github.com/althonos/pyhmmer/compare/v0.7.1...v0.7.2

### Added
- `easel.GeneticCode` class wrapping an `ESL_GENCODE` struct for configuring translation.
- `DigitalSequence.translate` method to translate a nucleotide sequence to a protein sequence. Metadata is copied from the source sequence to its translation ([#31](https://github.com/althonos/pyhmmer/pull/31), by [@valentynbez](https://github.com/valentynbez)).

### Deprecated
- `Sequence.taxonomy_id` property, as it is not used by Easel and implementation is not consistent (see [EddyRivasLab/easel#68](https://github.com/EddyRivasLab/easel/issues/68)).


## [v0.7.1] - 2022-12-15
[v0.7.1]: https://github.com/althonos/pyhmmer/compare/v0.7.0...v0.7.1

### Added
- Missing `__reduce__` method to `TopHits`.

### Fixed
- Build detection of available platform functions in `setup.py`.


## [v0.7.0] - 2022-12-04
[v0.7.0]: https://github.com/althonos/pyhmmer/compare/v0.6.3...v0.7.0

### Added
- `Bitfield.zeros` and `Bitfield.ones` classmethods for constructing an empty bitfield of known size.
- `Bitfield.copy` method to copy a bitfield object.
- `SequenceBlock` and `OptimizedProfileBlock` classes to store Python objects next to a contiguous array of pointers for iterating with the GIL released.
- `SequenceFile.read_block` method to read a whole sequence block from a file.
- `HMM.sample` class method to generate a HMM at random given a `Randomness` source.
- `hmmscan` function to scan a profile database with sequence queries.
- `deepcopy` implementations to `HMM`, `Profile` and `OptimizedProfile` classes of `plan7`.
- `rewind` method to `HMMFile`, `HMMPressedFile` and `SequenceFile` to reset a file back to its initial position.
- `name` attribute to `HMMFile`, `HMMPressedFile`, `MSAFile` and `SequenceFile` to expose the path of a file (when it was created from path).
- `local` property to `Profile` and `OptimizedProfile`, indicating whether a profile is in local or global mode.
- `multihit` property to `Profile` and `OptimizedProfile`, indicating whether a profile is in unihit or multihit mode, with a setter taking care of the reconfiguration.
- `Domain.included` and `Domain.reported` settable properties to report the inclusion and reporting status of a single domain.
- `TopHits.included` and `TopHits.reported` sized iterator to iterate only on included and reported hits.
- `Domains.included` and `Domains.reported` sized iterator to iterate only on included and reported domains.

### Changed
- `Bitfield`, `Vector` and `Matrix` can now be created from an iterable.
- `Pipeline` search methods now expect a `DigitalSequenceBlock` or a `SequenceFile` for the target sequence database.
- `Pipeline` scan methods now expect an `OptimizedProfileBlock` or a `HMMPressedFile` for the target profile database.
- `TraceAligner` now expect a `DigitalSequenceBlock` for the sequences to align to the HMM.
- `Profile.configure` now uses a default value of 400 for the `L` argument.
- `hmmsearch`, `nhmmer` and `phmmer` support being given a single query instead of requiring an iterable.
- `HMMPressedFile` can now be created, closed and used as a context manager directly without having to manage the source `HMMFile`.
- Renamed `Profile.optimized` method to `Profile.to_optimized`.
- Replaced `Randomness.is_fast` method with the `Randomness.fast` property.
- Rewrite handling of `Hit` flags using settable properties (`Hit.included`, `Hit.reported`, `Hit.new`, `Hit.dropped`, `Hit.duplicate`) instead of methods.

### Fixed
- Memory leak in the `LongTargetsPipeline` search loop.
- PyPy behaviour change of `readinto` methods now expecting `unsigned char*` instead of `char*` memoryview.
- `NULL`-pointer dereference in `Pipeline.search_hmm` when given a query without name.
- `LongTargetsPipeline` not recording the query name and accession.
- Memory leak caused by using a non-default prior scheme when constructing a `Builder`.

### Removed
- `PipelineSearchTargets`, replaced in functionality with `easel.DigitalSequenceBlock`.
- `is_local` and `is_multihit` methods of `Profile` and `OptimizedProfile`, replaced with equivalent properties.
- `Hit.manually_drop` and `Hit.manually_include` methods, replaced with the different `Hit` properties.


## [v0.6.3] - 2022-09-09
[v0.6.3]: https://github.com/althonos/pyhmmer/compare/v0.6.2...v0.6.3

### Fixed
- Error not being raised on alphabet detection failure in `SequenceFile` or `MSAFile`.
- Add check in `DigitalSequence` constructor to make sure encoded characters are in valid range ([#25](https://github.com/althonos/pyhmmer/issues/25)).

### Added
- `SequenceFile.guess_alphabet` and `MSAFile.guess_alphabet` to guess the alphabet from an open file.
- `Alphabet.encode` and `Alphabet.decode` to convert raw sequences between digital and text format.


## [v0.6.2] - 2022-08-12
[v0.6.2]: https://github.com/althonos/pyhmmer/compare/v0.6.1...v0.6.2

### Changed
- `hmmsearch`, `phmmer` and `nhmmer` functions will reduce the requested number of threads to the number of queries, if it can be detected using `operator.length_hint`.

### Added
- Documentation for loading all HMMs from an `HMMFile` object at once ([#23](https://github.com/althonos/pyhmmer/issues/23)).
- List of projects depending on PyHMMER to the `Examples` page of the documentation.


## [v0.6.1] - 2022-06-28
[v0.6.1]: https://github.com/althonos/pyhmmer/compare/v0.6.0...v0.6.1

### Added
- `pickle` protocol support for `TopHits` objects, using the HMMER network serialization.
- `TopHits.write` method to write hits to a file in tabular format.
- `query_name` and `query_accession` properties to `TopHits` objects to access the name and accession of the query that produced the hits.

### Fixed
- Extraction of filename from file-like objects in the `HMMFile` constructor.
- Use `os.cpu_count` instead of `multiprocessing.cpu_count` where applicable to preserve OS scheduling.
- Wrong return type in docstring of `HMM.insert_emissions`.
- `TopHits.searched_nodes` returning the searched number of residues instead of the searched number of model nodes.
- Unsound decoding of pickled `MatrixF` or `VectorF` when data comes from a source of different endianness.

### Changed
- Rewrite `pyhmmer.hmmer` threading code using `Deque` instead of `collections.Queue` to store the queries and results.
- Reduce memory consumption of `pyhmmer.hmmer` by reducing the number of semaphores and event flags used concurrently.
- Make `pyhmmer.hmmer` main threads block on query insertion rather than result retrieval to make sure worker threads are never idling.


## [v0.6.0] - 2022-05-01
[v0.6.0]: https://github.com/althonos/pyhmmer/compare/v0.5.0...v0.6.0

### Added
- `pyhmmer.daemon` module with an client implementation to communicate to a `hmmpgmd` server.
- `Pipeline.arguments` methods to get a list of CLI arguments from the parameters used to initialize the `Pipeline`.
- Setters for `name`, `accession` and `description` properties of `plan7.Hit`.
- Constructor for individual `plan7.Trace` objects outside a `plan7.Traces` list.
- `plan7.Trace.from_sequence` constructor to create a faux trace from a single sequence.
- `manually_include` and `manually_drop` methods to `plan7.Hit` for manually selecting the inclusion status of a `Hit` in a `TopHits` instance.
- `compare_ranking` method to `plan7.TopHits` for comparing the order of the hits compared to a previous run on the same targets stored in an `easel.KeyHash` object.
- `Pipeline.iterate_seq` and `Pipeline.iterate_hmm` to run iterative queries like JackHMMER.
- `repr` implementations for `easel.MSAFile`, `easel.SequenceFile` and `easel.HMMFile` showing the path or file object they were created from.
- `repr` implementation for `easel.Randomness` showing the seed and the RNG algorithm in use.
- `str` implementation for `plan7.Alignment` using HMMER original code to display a domain alignment like in search/scan results.

### Changed
- `plan7.Trace.posterior_probabilities` property may now be `None` in case no memory is allocated for the posteriors in the `P7_TRACE` struct.
- `TopHits.to_msa` can now add additional sequences passed as arguments to the alignment.
- `plan7.HMMPressedFile` now raises an exception on attempts to create a new instance manually.
- `ignore_gaps` argument of `easel.SequenceFile` is now deprecated.
- `repr` implementations for `easel` types now use the fully qualified class name.

### Fixed
- `easel.SequenceFile.readinto` docstring not rendering properly in documentation.
- Type annotations of `hits_included` and `hits_reported` of `plan7.TopHits` marking these properties as `bool` instead of `int`.
- Setters of `name`, `accession`, `description` and `author` properties of `easel.MSA` crashing when given `None` values.
- Exception value raised from Easel code not being properly extracted.
- Plain strings being used in example for `easel.TextSequence` and `easel.TextMSA` constructors where byte strings are expected ([#20](https://github.com/althonos/pyhmmer/pull/20)).


## [v0.5.0] - 2022-03-14
[v0.5.0]: https://github.com/althonos/pyhmmer/compare/v0.4.11...v0.5.0

### Added
- `plan7.PipelineSearchTargets` to reduce the overhead when searching the same sequences several times with different. query profiles.
- `TopHits.copy` method to duplicate a `TopHits` instance.
- `TopHits.merge` method to merge hits obtained with the same query on different targets.
- Buffer protocol implementation for `pyhmmer.easel.Bitfield`.

### Changed
- Renamed `TopHits.included` and `TopHits.reported` properties to `TopHits.hits_included` and `TopHits.hits_included`.
- `MSAFile` and `SequenceFile` are now directly in digital mode if they are instantiated with `digital=True`.
- `SequenceFile.parse` can now return a sequence in digital mode.
- Reorganized tests to make then runnable from a site install.

### Fixed
- Usage of `memcpy` in contexts where it may have had undefined behaviour.
- `VectorF.__eq__` crashing when comparing two empty objects.
- `SequenceFile` and `MSAFile` not closing file handles when raising an error in `__init__`.


## [v0.4.11] - 2021-12-15
[v0.4.11]: https://github.com/althonos/pyhmmer/compare/v0.4.10...v0.4.11

### Added
- `plan7.HMMFile.read` method to read a single `plan7.HMM` from an `plan7.HMMFile` (instead of using `next`).
- `closed` property on `easel.SequenceFile`, `easel.MSAFile` and `plan7.HMMFile` to mark whether a file object is closed.
- `plan7.HMMFile.is_pressed` method to check whether a HMM file has associated pressed data.
- `plan7.HMMFile.optimized_profiles` methods to read the `plan7.OptimizedProfile` entries in an `plan7.HMMFile` is there are associated pressed data available.
- Getters for the `name`, `accession`, `description`, `consensus`, `consensus_structure`, `evalue_parameters` and `cutoffs` properties of a `plan7.OptimizedProfile`.
- `plan7.OptimizedProfile.__eq__` implementation to compare two optimized profiles.
- `__sizeof__` implementations for `plan7.OptimizedProfile` and `plan7.Profile` to get the allocated size of a profile.

### Fixed
- Double-free caused by the Cython cycle breaking feature on several view types (`easel.Randomness`, `easel.Vector`, `easel.Matrix`, `plan7.Cutoffs`, `plan7.EvalueParameters`, `plan7.Offsets`, `plan7.Trace`)
- `plan7.Hit.description` using the pointer to the accession string erroneously, causing occasional NULL dereference.
- `plan7.OptimizedProfile.copy` performing a shallow copy instead of a deep copy as expected.

### Changed
- `pyhmmer.hmmer` type annotations now explicit support for `plan7.Profile` or `plan7.OptimizedProfile` inputs where applicable.


## [v0.4.10] - 2021-12-06
[v0.4.10]: https://github.com/althonos/pyhmmer/compare/v0.4.9...v0.4.10

### Added
- `entropy` and `relative_entropy` methods to `easel.VectorF` to compute the Shannon entropy of a vector and the Kullback-Leibler divergence of two vectors.
- `mean_match_entropy`, `mean_match_information` and `mean_match_relative_entropy` methods to `plan7.HMM` to get information statistics of an HMM model.
- `match_occupancy` method to `plan7.HMM` to compute the occupancy for each match state as an `easel.VectorF`.

### Fixed
- `plan7.Builder.build_msa` using the gap-open and gap-extend probabilities instead of the MSA itself to compute the transition probabilities for the new HMM.

### Changed
- `plan7.Builder.build` will now only load the score system once and reuse it unless a different score system is requested between calls.


## [v0.4.9] - 2021-11-11
[v0.4.9]: https://github.com/althonos/pyhmmer/compare/v0.4.8...v0.4.9

### Added
- `plan7.ScoreData` class to store the substitution scores and maximal extensions for a long target search.
- `plan7.LongTargetsPipeline` to run searches on targets longer than 100,000 residues.
- `Alphabet` methods to check whether an `Alphabet` object is a DNA, RNA, nucleotide or protein alphabet.
- `window_length` and `window_beta` arguments to `plan7.Builder` to set the max length of nucleotide `HMM` created by builder objects.

### Changed
- `pyhmmer.hmmer.nhmmer` now uses a `LongTargetsPipeline` instead of a `Pipeline` to search the target sequences.
- `pyhmmer.hmmer.nhmmer` now supports `HMM` queries in addition to `DigitalSequence` and `DigitalMSA` queries.
- `pyhmmer.hmmer.phmmer` now always assumes protein queries.
- `Z` and `domZ` attributes of `plan7.TopHits` objects is now read-only.

### Fixed
- `nhmmer` now uses DNA as the default alphabet instead of amino acid alphabet like it did before ([#12](https://github.com/althonos/pyhmmer/issues/12)).


## [v0.4.8] - 2021-10-27
[v0.4.8]: https://github.com/althonos/pyhmmer/compare/v0.4.7...v0.4.8

### Added
- Constructor arguments and properties to `plan7.Pipeline` to support bit score thresholds instead to filter top hits.
- Support for creating a `SequenceFile` and an `MSAFile` using a Python file-like object instead of only supporting filenames.
- Support for reading individual sequences from an MSA file with `SequenceFile`.
- `TextMSA.alignment` to access the actual alignment as a tuple of strings.
- Subtraction and division support for `easel.Vector` subclasses

### Changed
- `plan7.Cutoffs` now support setting the bit score cutoffs, but requires both to be set or cleared at the same time.
- `easel.Vector` will always allocate some memory when created manually to avoid having a special empty case in every vector method.
- `pyhmmer.easel.AllocationError` now stores the size it failed to allocate, and the number of elements when allocating an array.

### Fixed
- `TextSequence.digitize` will not raise a `ValueError` when the sequence contains invalid characters for the alphabet (previously was an `UnexpectedError`).


## [v0.4.7] - 2021-09-28
[v0.4.7]: https://github.com/althonos/pyhmmer/compare/v0.4.6...v0.4.7

### Added
- `TraceAligner`, `Trace` and `Traces` classes to `pyhmmer.plan7` to get tracebacks after aligning several sequences against an HMM.
- `pyhmmer.hmmalign` function with the same features as the `hmmalign` binary from HMMER3.
- Support for out-of-band pickling in `easel.Vector` and `easel.Matrix`.

### Changed
- Allow creating an empty `Vector` or `Matrix` by calling their constructor without arguments.

### Fixed
- Potential unreported exceptions in `plan7.OptimizedProfile.write` and several `plan7.SSIWriter` methods.


## [v0.4.6] - 2021-09-10
[v0.4.6]: https://github.com/althonos/pyhmmer/compare/v0.4.5...v0.4.6

### Added
- `pickle` protocol for `easel.Alphabet`, `easel.Bitfield`, `easel.KeyHash`, `easel.Vector`, `easel.Matrix` and `plan7.HMM`.
- `taxonomy_id` and `residue_markups` properties to `easel.Sequence`.
- `sum_score` property to `plan7.Hit`.
- `plan7.EvalueParameters` class to expose the e-value parameters of a `plan7.HMM` or a `plan7.Profile`.
- Equality checks and slicing for `easel.Matrix` and `easel.Vector`.
- Support for creating and manipulating zero-sized `easel` matrices and vectors.
- `plan7.Cutoffs` class to expose the Pfam score cutoffs of a `plan7.HMM` or a `plan7.Profile`.
- Keyword arguments to configure E-value thresholds when creating a `plan7.Pipeline` object.
- Support for using model-specific thresholding options in `plan7.Pipeline`.

### Changed
- Use the *replace* error handler when decoding error messages to skip potential decoding issues when already building an exception.
- Improve `pyhmmer.hmmer` to ensure background threads exit on a `KeyboardInterrupt`.
- `easel.VectorU8.__eq__` accepts any object implementing the buffer protocol.
- `plan7.HMM.creation_time` now takes and returns a `datetime.datetime` object, assuming the field is only ever set with `asctime`.
- Refactor `easel.Vector` and `easel.Matrix` and mark exposed memory as C-contiguous.

### Fixed
- `easel.Alphabet` not reporting potential allocation errors.
- Potential buffer overflow in `easel.Matrix` and `easel.Vector` when calling `__init__` more than once.


## [v0.4.5] - 2021-07-19
[v0.4.5]: https://github.com/althonos/pyhmmer/compare/v0.4.4...v0.4.5

### Added
- `OptimizedProfile.convert` method to configure an optimized profile from a `Profile` without reallocating a new `P7_OPROFILE` struct.

### Changed
- Rewrite the `plan7.Pipeline` search loop to avoid reacquiring the GIL between reference sequences.
- Require the reference sequences to be stored in a collection (instead of an iterable) when passing them to the `search_hmm`, `search_msa` and `search_seq` methods of `plan7.Pipeline`.
- Avoid reallocating a new `OptimizedProfile` every time a new HMM is passed to `Pipeline.search_hmm`.
- Relax the GIL while sorting and thresholding `TopHits` in `Pipeline` search methods.


## [v0.4.4] - 2021-07-07
[v0.4.4]: https://github.com/althonos/pyhmmer/compare/v0.4.3...v0.4.4

### Added
- `ignore_gaps` parameter to `pyhmmer.plan7.SequenceFile`, allowing to skip the gap characters when reading a sequence from an ungapped format.
- `__sizeof__` implementation for some
- Dedicated check for sequence length before running the platform-specific code in `pyhmmer.plan7.Pipeline`.

### Fixed
- Score system not being set in `pyhmmer.plan7.Builder.build_msa`.
- Alphabet not being checked after the first sequence in `Pipeline` search and scan methods.


## [v0.4.3] - 2021-07-03
[v0.4.3]: https://github.com/althonos/pyhmmer/compare/v0.4.2...v0.4.3

### Fixed
- File object wrappers not reporting exceptions raised when seeking on OSX/BSD platforms.


## [v0.4.2] - 2021-06-20
[v0.4.2]: https://github.com/althonos/pyhmmer/compare/v0.4.1...v0.4.2

### Added
- `pyhmmer.easel.Randomness` class exposing a deterministic random number generator.
- `pyhmmer.plan7.Builder.randomness` and `pyhmmer.plan7.Pipeline.randomness` attributes exposing the internal random number generator used by each object.
- `pyhmmer.plan7.Hit.best_domain` property mapping to the highest scoring domain of a hit.
- `pyhmmer.plan7.OptimizedProfile.rbv` property exposing match scores.
- `pyhmmer.plan7.Domain.pvalue` and `pyhmmer.plan7.Hit.pvalue` reporting the p-value for a domain or hit bitscore.

### Fixed
- Dimensions of the `pyhmmer.plan7.OptimizedProfile.sbv` matrix not being properly set.


## [v0.4.1] - 2021-06-06
[v0.4.1]: https://github.com/althonos/pyhmmer/compare/v0.4.0...v0.4.1

### Fixed
- Main buffer not being freed in `MatrixF.__dealloc__` and `MatrixU8.__dealloc__` when created without owner.

### Added
- Additional configuration values for `pyhmmer.plan7.Pipeline` as both constructor arguments and mutable properties.
- `consensus`, `consensus_structure` and `offsets` properties to `pyhmmer.plan7.Profile` objects.

### Changed
- Make `OptimizedProfile.ssv_filter` check the alphabet of the given sequence.


## [v0.4.0] - 2021-06-05 - YANKED
[v0.4.0]: https://github.com/althonos/pyhmmer/compare/v0.3.1...v0.4.0

### Added
- Linear algebra primitives to expose 1D (`Vector`) and 2D (`Matrix`) contiguous buffers containing numerical values to `pyhmmer.easel`.
- Documentation for the `Z` and `domZ` parameters of the `pyhmmer.plan7.Pipeline` constructor.
- `pyhmmer.errors.AlphabetMismatch` exception deriving from `ValueError` to specifically report mismatching Easel alphabets where applicable.
- `scale` and `normalize` methods to `pyhmmer.plan7.HMM` objects.
- Property to access `pyhmmer.plan7.Background` residue frequencies as a `VectorF` object.
- Property to access `pyhmmer.plan7.HMM` mean residue composition as a `VectorF` object.
- Property to access `pyhmmer.plan7.HMM` probabilities and emissions as `MatrixF` objects.
- `ssv_filter` methods to `pyhmmer.plan7.OptimizedProfile` to get the SSV filter score of the profile for a given sequence.
- Several additional properties to access the `pyhmmer.plan7.OptimizedProfile` internals.

### Removed
- Unused `report_e` parameter of `pyhmmer.plan7.Pipeline` constructor.
- `pyhmmer.plan7.TopHits.clear` method which could lead to segfault if it was called while a `Hit` is being held.

### Changed
- Multithreaded loop in `pyhmmer.hmmer` to reduce memory consumption while still yielding hits in order.
- `pyhmmer.easel.DigitalSequence.sequence` property is now a `VectorU8`.

### Fixed
- Type annotations in `pyhmmer.hmmer`.
- Potential double free in `pyhmmer.plan7.HMM.command_line` property setter.
- Minor floating-point precision issues in `pyhmmer.plan7.Builder` constructor.
- Segfault in `TextMSA.digitize` caused by `esl_msa_Copy` not digitizing on-the-fly like `esl_sq_Copy`.
- Exceptions not being raised in some methods of `pyhmmer.plan7.Profile` and `pyhmmer.plan7.TopHits`.


## [v0.3.1] - 2021-05-08
[v0.3.1]: https://github.com/althonos/pyhmmer/compare/v0.3.0...v0.3.1

### Added
- `Pipeline.scan_seq` method to query a database of profiles with one or more sequences.
- `transition_probabilities`, `match_emissions`, `insert_emissions` properties to the `HMM` class, providing access to the numerical parameters of the HMM.
- `consensus_structure` and `consensus_accessibility` properties to the `HMM` class to get consensus lines from the source alignment if the HMM was created from a MSA.
- `nseq` and `nseq_effective` properties to the `HMM` class to get the number of training sequences and effective sequences used to build the HMM.

### Changed
- `HMM.checksum` is now `None` if the `p7H_CHKSUM` flag is not set.
- `Builder` methods will now record `sys.argv` when creating a HMM.

### Fixed
- `HMM.write(..., binary=False)` crashing on HMMs without a consensus line. ([#5](https://github.com/althonos/pyhmmer/issues/5)). Fixed upstream in ([EddyRivasLab/HMMER#236](https://github.com/EddyRivasLab/hmmer/pull/236)).
- `Pipeline.reset` mishandling the `Z` and `domZ` values if those were detected from the number of targets.
- `pyhmmer.hmmer` functions will not block until all results have been collected anymore when run in multithreaded mode.


## [v0.3.0] - 2021-03-11
[v0.3.0]: https://github.com/althonos/pyhmmer/compare/v0.2.2...v0.3.0

### Added
- `easel.MSAFile` to read from a file containing
- `accession`, `author`, `name` and `description` properties to `easel.MSA` objects.
- `plan7.Builder.build_msa` to build a pHMM from a sequence alignment.
- Additional methods to `easel.KeyHash`, allowing to use it as a `dict`/`set` hybrid.
- `Sequence.write` and `MSA.write` methods to format a sequence or an alignment to a file handle.
- `plan7.TopHits.to_msa` method to convert all the top hits of a query against a database into a multiple sequence alignment.
- `easel.MSA.sequences` attribute to access individual sequences of an alignment using the `collections.abc.Sequence` interface.
- `easel.DigitalMSA.textize` method to convert a multiple sequence alignment in digital mode to its text-mode counterpart.
- Read-only `name`, `accession` and `description` properties to `plan7.Profile` showing attributes inherited from the HMM it was configured with.
- `plan7.HMM.consensus` property, allowing to access the consensus sequence of a pHMM.
- `plan7.HMM` equality implementation, using zero tolerance.
- `plan7.Pipeline.search_msa` to query a MSA against a sequence database.
- `easel.Sequence.reverse_complement` method allowing to reverse-complement inplace or to build a copy.
- `errors.AlphabetMismatch` exception for use in cases where an alphabet is expected but not matched by the input.
- `hmmer.nhmmer` function with the same behaviour as `hmmer.phmmer`, except it expects inputs with a DNA alphabet.

### Fixed
- `plan7.Builder.copy` not copying some parameters correctly, causing `pyhmmer.hmmer.phmmer` to give inconsistent results in multithreaded mode.
- `easel.Bitfield` not properly handling index overflows.
- Documentation not rendering for the `__init__` method of all classes.

### Changed
- `plan7.Builder` gap-open and gap-extend probabilities are now set on instantiation and depend on the alphabet type.
- Constructors for `easel.TextMSA` and `easel.DigitalMSA`, which can now be given an iterable of `easel.Sequence` objects to store in the alignment.

### Removed
- Unimplemented `easel.SequenceFile.fetch` and `easel.SequenceFile.fetchinto` methods.


## [v0.2.2] - 2021-03-04
[v0.2.2]: https://github.com/althonos/pyhmmer/compare/v0.2.1...v0.2.2

### Fixed
- Linking issues on OSX caused by aggressive stripping of intermediate libraries.
- `plan7.Builder` RNG not reseeding between different HMMs.


## [v0.2.1] - 2021-01-29
[v0.2.1]: https://github.com/althonos/pyhmmer/compare/v0.2.0...v0.2.1

### Added
- `pyhmmer.plan7.HMM.checksum` property to get the 32-bit checksum of an HMM.


## [v0.2.0] - 2021-01-21
[v0.2.0]: https://github.com/althonos/pyhmmer/compare/v0.1.4...v0.2.0

### Added
- `pyhmmer.plan7.Builder` class to handle building a HMM from a sequence.
- `Pipeline.search_seq` to query a sequence against a sequence database.
- `psutil` dependency to detect the most efficient thread count for `hmmsearch` based on the number of *physical* CPUs.
- `pyhmmer.hmmer.phmmer` function to run a search of query sequences against a sequence database.

### Changed
- `Pipeline.search` was renamed to `Pipeline.search_hmm` for disambiguation.
- `libeasel.random` sequences do not require the GIL anymore.
- Public API now have proper signature annotations.

### Fixed
- Inaccurate exception messages in `Pipeline.search_hmm`.
- Unneeded RNG reallocation, replaced with re-initialisation where possible.
- `SequenceFile.__next__` not working after being set in digital mode.
- `sequences` argument of `hmmsearch` now only requires a `typing.Collection[DigitalSequence]` instead of a `typing.Collection[Sequence]` (not more `__getitem__` needed).

### Removed
- `hits` argument to `Pipeline.search_hmm` to reduce risk of issues with `TopHits` reuse.
- Broken alignment coordinates on `Domain` classes.


## [v0.1.4] - 2021-01-15
[v0.1.4]: https://github.com/althonos/pyhmmer/compare/v0.1.3...v0.1.4

### Added
- `DigitalSequence.textize` to convert a digital sequence to a text sequence.
- `DigitalSequence.__init__` method allowing to create a digital sequence from any object implementing the buffer protocol.
- `Alignment.hmm_accession` property to retrieve the accession of the HMM in an alignment.


## [v0.1.3] - 2021-01-08
[v0.1.3]: https://github.com/althonos/pyhmmer/compare/v0.1.2...v0.1.3

### Fixed
- Compilation issues in OSX-specific Cython code.


## [v0.1.2] - 2021-01-07
[v0.1.2]: https://github.com/althonos/pyhmmer/compare/v0.1.1...v0.1.2

### Fixed
- Required Cython files not being included in source distribution.


## [v0.1.1] - 2020-12-02
[v0.1.1]: https://github.com/althonos/pyhmmer/compare/v0.1.0...v0.1.1

### Fixed
- `HMMFile` calling `file.peek` without arguments, causing it to crash when passed some types, e.g. `gzip.GzipFile`.
- `HMMFile` failing to work with PyPy file objects because of a bug with their implementation of `readinto`.
- C/Python file object implementation using `strcpy` instead of `memcpy`, causing issues when null bytes were read.


## [v0.1.0] - 2020-12-01
[v0.1.0]: https://github.com/althonos/pyhmmer/compare/v0.1.0-a5...v0.1.0

Initial beta release.

### Fixed
- `TextSequence` uses the sequence argument it's given on instantiation.
- Segmentation fault in `Sequence.__eq__` caused by implicit type conversion.
- Segmentation fault on `SequenceFile.read` failure.
- Missing type annotations for the `pyhmmer.easel` module.


## [v0.1.0-a5] - 2020-11-28
[v0.1.0-a5]: https://github.com/althonos/pyhmmer/compare/v0.1.0-a4...v0.1.0-a5

### Added
- `Sequence.__len__` magic method so that `len(seq)` returns the number of letters in `seq`.
- Python file-handle support when opening an `pyhmmer.plan7.HMMFile`.
- Context manager protocol to `pyhmmer.easel.SSIWriter`.
- Type annotations for `pyhmmer.easel.SSIWriter`.
- `add_alias` to `pyhmmer.easel.SSIWriter`.
- `write` method to `pyhmmer.plan7.OptimizedProfile` to write an optimized profile in binary format.
- `offsets` property to interact with the disk offsets of a `pyhmmer.plan7.OptimizedProfile` instance.
- `pyhmmer.hmmer.hmmpress` emulating the `hmmpress` binary from HMMER.
- `M` property to `pyhmmer.plan7.HMM` exposing the number of nodes in the model.

### Changed
- Bumped vendored Easel to `v0.48`.
- Bumped vendored HMMER to `v3.3.2`.
- `pyhmmer.plan7.HMMFile` will raise an `EOFError` when given an empty file.
- Renamed `length` property to `L` in `pyhmmer.plan7.Background`.

### Fixed
- Segmentation fault when `close` method of `pyhmmer.easel.SSIWriter` was called more than once.
- `close` method of `pyhmmer.easel.SSIWriter` not writing the index contents.


## [v0.1.0-a4] - 2020-11-24
[v0.1.0-a4]: https://github.com/althonos/pyhmmer/compare/v0.1.0-a3...v0.1.0-a4

### Added
- `MSA`, `TextMSA` and `DigitalMSA` classes representing a multiple sequence alignment to `pyhmmer.easel`.
- Methods and protocol to copy a `Sequence` and a `MSA`.
- `pyhmmer.plan7.OptimizedProfile` wrapping a platform-specific optimized profile.
- `SSIReader` and `SSIWriter` classes interacting with sequence/subsequence indices to `pyhmmer.easel`.
- Exception handler using Python exceptions to report Easel errors.

### Changed
- `pyhmmer.hmmsearch` returns an iterator of `TopHits`, with one instance per `HMM` in the input.
- `pyhmmer.hmmsearch` properly raises errors happenning in the background threads without deadlock.
- `pyhmmer.plan7.Pipeline` recycles memory between `Pipeline.search` calls.

### Fixed
- Missing type annotations for the `pyhmmer.errors` module.

### Removed
- Unneeded or private methods from `pyhmmer.plan7`.


## [v0.1.0-a3] - 2020-11-19
[v0.1.0-a3]: https://github.com/althonos/pyhmmer/compare/v0.1.0-a2...v0.1.0-a3

### Added
- `TextSequence` and `DigitalSequence` representing a `Sequence` in a given mode.
- E-value properties to `Hit` and `Domain`.
- `TopHits` now stores a reference to the pipeline it was obtained from.
- `Pipeline.Z` and `Pipeline.domZ` properties.
- Experimental pickling support to `Alphabet`.
- Experimental freelist to `Sequence` class to avoid allocation bottlenecks when iterating on a `SequenceFile` without recycling sequence buffers.

### Changed
- Made `Sequence` an abstract base class.
- Additional `Pipeline` parameters can be passed as keyword arguments to `pyhmmer.hmmsearch`.
- `SequenceFile.read` can now be configured to skip reading the metadata or the content of a sequence.

### Removed
- Redundant `SequenceFile` methods.

### Fixed
- `doctest` loader crashing on Python 3.5.
- `TopHits.threshold` segfaulting when being called without prior `Tophits.sort` call
- Unknown `format` argument to `SequenceFile` constructor not raising the right error.


## [v0.1.0-a2] - 2020-11-12
[v0.1.0-a2]: https://github.com/althonos/pyhmmer/compare/v0.1.0-a1...v0.1.0-a2

### Added
- Support for compilation on PowerPC big-endian platforms.
- Type annotations and stub files for Cython modules.

### Changed
- [`distutils`](https://docs.python.org/3/library/distutils.html) is now used to compile the package, instead of calling `autotools` and letting HMMER configure itself.
- `Bitfield.count` now allows passing an argument (for compatibility with [`collections.abc.Sequence`](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)).


## [v0.1.0-a1] - 2020-11-10
[v0.1.0-a1]: https://github.com/althonos/pyhmmer/compare/fe4c279...v0.1.0-a1

Initial alpha release (test deployment to PyPI).
