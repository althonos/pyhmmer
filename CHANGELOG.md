# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyhmmer/compare/v0.1.0-a2...HEAD


## [v0.1.0-a2] - 2020-11-12
[v0.1.0-a2]: https://github.com/althonos/pyhmmer/compare/v0.1.0-a1...v0.1.0-a2

### Added
- Support for compilation on PowerPC big-endian platforms.
- Type annotations and stub files for Cython modules.

### Changed
- [`distutils`](https://docs.python.org/3/library/distutils.html) is now
  used to compile the package, instead of calling `autotools` and letting
  HMMER configure itself.
- `Bitfield.count` now allows passing an argument (for compatibility with
  [`collections.abc.Sequence`](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)).


## [v0.1.0-a1] - 2020-11-10
[v0.1.0-a1]: https://github.com/althonos/pyhmmer/compare/fe4c279...v0.1.0-a1

Initial alpha release (test deployment to PyPI).
