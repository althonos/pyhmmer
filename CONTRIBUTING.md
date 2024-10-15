# Contributing to PyHMMER

For bug fixes or new features, please file an issue before submitting a
pull request. If the change isn't trivial, it may be best to wait for
feedback.

## Setting up a local repository

Make sure you clone the repository in recursive mode, so you also get the
wrapped code of Easel and HMMER, which are exposed as `git` submodules:

```console
$ git clone --recursive https://github.com/althonos/pyhmmer
```

## Running tests

Tests are written as usual Python unit tests with the `unittest` module of
the standard library. Running them requires the extension to be built
locally:

```console
$ python -m pip install -e . -v --no-build-isolation
$ python -m unittest pyhmmer.tests -vv
```

## Running benchmarks

Benchmarks require the HMMER binaries to be available in the `$PATH`,
as well as the [`hyperfine`](https://github.com/sharkdp/hyperfine) 
binary. Then make sure to have the extensions compiled and installed
in release mode:

```console
$ python -m pip install . -v --no-build-isolation
```

Then run the benchmarks using the UNIX shell script:

```console
$ ./benches/run.sh
```

## Coding guidelines

This project targets Python 3.7 or later.

Python objects should be typed; since it is not supported by Cython,
you must manually declare types in type stubs (`.pyi` files). In Python
files, you can add type annotations to function signatures (supported in
Python 3.5) or in variable assignments (supported from Python 3.6
onward). 

### Interfacing with C

When interfacing with C, and in particular with pointers, use assertions
everywhere you assume the pointer to be non-NULL.
