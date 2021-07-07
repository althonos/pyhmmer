# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyhmmer/compare/v0.4.4...HEAD


## [v0.4.4] - 2021-07-07
[v0.4.3]: https://github.com/althonos/pyhmmer/compare/v0.4.3...v0.4.4

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
- `plan7.TopHist.to_msa` method to convert all the top hits of a query against a database into a multiple sequence alignment.
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
