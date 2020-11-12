# 🐍🟡♦️🟦 pyHMMER [![Stars](https://img.shields.io/github/stars/althonos/pyhmmer.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyhmmer/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [HMMER3](http://hmmer.org/).*

[![TravisCI](https://img.shields.io/travis/com/althonos/pyhmmer/master.svg?logo=travis&maxAge=600&style=flat-square)](https://travis-ci.com/althonos/pyhmmer/branches)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyhmmer?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyhmmer/)
[![PyPI](https://img.shields.io/pypi/v/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer)
[![Wheel](https://img.shields.io/pypi/wheel/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer/#files)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyhmmer/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyhmmer.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyhmmer/issues)
[![Docs](https://img.shields.io/readthedocs/pyhmmer/latest?style=flat-square&maxAge=600)](https://pyhmmer.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyhmmer/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpyhmmer)](https://pepy.tech/project/pyhmmer)

<!-- [![AppVeyor](https://img.shields.io/appveyor/build/althonos/pyhmmer/master.svg?logo=appveyor&maxAge=600&style=flat-square)](https://ci.appveyor.com/project/althonos/pyhmmer/history) -->
<!-- [![Bioconda](https://img.shields.io/conda/vn/bioconda/pyhmmer?style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/pyhmmer) -->


## 🗺️ Overview

HMMER is a biological sequence analysis tool that uses profile hidden Markov
models to search for sequence homologs. HMMER3 is maintained by members of the
the [Eddy/Rivas Laboratory](http://eddylab.org/) at Harvard University.

`pyhmmer` is a Python module, implemented using the [Cython](https://cython.org/)
language, that provides bindings to HMMER3. It directly interacts with the
HMMER internals, which has the following advantages over CLI wrappers
(like [`hmmer-py`](https://pypi.org/project/hmmer/)):

- **single dependency**: If your software or your analysis pipeline is
  distributed as a Python package, you can add `pyhmmer` as a dependency to
  your project, and stop worrying about the HMMER binaries being properly
  setup on the end-user machine.
- **no intermediate files**: Everything happens in memory, in Python objects
  you have control on, making it easier to format your inputs to pass to
  HMMER without needing to write them to a file. Output retrieval is also
  done in memory, through instances of the `pyhmmer.plan7.TopHits` class.
- **no input formatting**: The Easel object model is exposed in the `pyhmmer.easel`
  module, and you have the possibility to build a `Sequence` object yourself
  to pass to the HMMER pipeline. This is useful if your sequences are already
  loaded in memory, for instance because you obtained them from another
  Python library (such as [Pyrodigal](https://github.com/althonos/pyrodigal)
  or [Biopython](https://biopython.org/)).
- **no output formatting**: HMMER3 is notorious for its numerous output files
  and its fixed-width tabular output, which is hard to parse (even
  [`Bio.SearchIO.HmmerIO`](https://biopython.org/docs/dev/api/Bio.SearchIO.HmmerIO.html)
  is struggling on some sequences).
- **no performance penalty**: Using `pyhmmer` to launch `hmmsearch` on sequences
  and HMMs in disk storage is typically not slower than directly using the
  `hmmsearch` binary (see the [Benchmarks section](#%EF%B8%8F-benchmarks)).
  If you already have everything in memory, saving the writing/reading step
  possibly gives you a slight boost!

*This library is still a work-in-progress, and in a very experimental stage,
but it should already pack enough features to run simple biological analyses
involving `hmmsearch`.*


## 🔧 Installing

``pyhmmer`` can be installed from [PyPI](https://pypi.org/project/pyhmmer/),
which hosts some pre-built CPython wheels for x86-64 Linux and OSX,
as well as the code required to compile from source with Cython:
```console
$ pip install pyhmmer
```

Compilation for UNIX PowerPC is not tested in CI, but should work out of the
box. Other architectures (e.g. Arm) and OSes (e.g. Windows) are not
supported by HMMER.


## 📖 Documentation

A complete [API reference](https://pyhmmer.readthedocs.io/en/stable/api/) can
be found in the [online documentation](https://pyhmmer.readthedocs.io/), or
directly from the command line using
[`pydoc`](https://docs.python.org/3/library/pydoc.html):
```console
$ pydoc pyhmmer.easel
$ pydoc pyhmmer.plan7
```


## 💡 Example

Use `pyhmmer` to run `hmmsearch`, and display alignments for hits with
domains longer than 30 letters :

```python
import pyhmmer

# load sequences from the FASTA file
# (HMMER normally rewinds the file containing the sequences and reads it
# again for every HMM, but we can save some time by loading them ahead of
# time, provided the machine has enough memory)
with pyhmmer.easel.SequenceFile("938293.PRJEB85.HG003687.faa") as file:
    sequences = list(file)

# open the HMM file and use it to run `hmmsearch`
# (we don't need to collect ahead of time here, we can just pass an iterator
# to `pyhmmer.hmmsearch`).
with pyhmmer.plan7.HMMFile("Pfam.hmm") as hmms:
    hits = pyhmmer.hmmsearch(hmms, sequences_file)

# find hit domains where the alignment is longer than 30 letters and the
# domain score is higher than 9, then print the aligned sequences
for hit in hits:
    for domain in hit.domains:
        if domain.score > 9 and len(domain.alignment) > 30:
            print(hit.name, domain.alignment.hmm_name)
            print(domain.alignment.hmm_sequence)
            print(domain.alignment.identity_sequence)
            print(domain.alignment.target_sequence)
```


## ⏱️ Benchmarks

Benchmarks were run on a i7-8550U CPU running at 1.80GHz, using a FASTA file
containing 2100 protein sequences (`tests/data/seqs/938293.PRJEB85.HG003687.faa`)
and a subset of the [Pfam](https://pfam.xfam.org/) HMM library containing
2873 domains. Both commands were set to use 8 threads, and were run 20 times.

| Command                     | mean (s) | σ (ms) | min (s) | max (s) |
|-----------------------------|----------|--------|---------|---------|
| `hmmsearch`                 | 40.723   | 223    | 40.402  | 41.321  |
| `python -m hmmer.hmmsearch` | 40.128   | 946    | 37.706  | 42.457  |


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The HMMER3 and Easel code is available under the BSD 3-clause license. See
`vendor/hmmer/LICENSE` and `vendor/easel/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the [original HMMER authors](http://hmmer.org/). It was developed by
[Martin Larralde](https://github.com/althonos/pyhmmer) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
