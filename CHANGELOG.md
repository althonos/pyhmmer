# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyhmmer/compare/v0.1.2...HEAD


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
