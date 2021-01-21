# 🐍🟡♦️🟦 pyHMMER [![Stars](https://img.shields.io/github/stars/althonos/pyhmmer.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyhmmer/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [HMMER3](http://hmmer.org/).*

[![GitLabCI](https://img.shields.io/gitlab/pipeline/larralde/pyhmmer/master?gitlab_url=https%3A%2F%2Fgit.embl.de&logo=gitlab&style=flat-square&maxAge=600)](https://git.embl.de/larralde/pyhmmer/-/pipelines)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyhmmer?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyhmmer/)
[![PyPI](https://img.shields.io/pypi/v/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyhmmer?style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/pyhmmer)
[![Wheel](https://img.shields.io/pypi/wheel/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyhmmer.svg?style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/pyhmmer/#files)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyhmmer/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pyhmmer/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyhmmer.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyhmmer/issues)
[![Docs](https://img.shields.io/readthedocs/pyhmmer/latest?style=flat-square&maxAge=600)](https://pyhmmer.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyhmmer/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpyhmmer)](https://pepy.tech/project/pyhmmer)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.4270012-purple?style=flat-square&maxAge=86400)](https://doi.org/10.5281/zenodo.4270012)

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
  you have control on, making it easier to pass your inputs to HMMER without 
  needing to write them to a temporary file. Output retrieval is also done 
  in memory, via instances of the 
  [`pyhmmer.plan7.TopHits`](https://pyhmmer.readthedocs.io/en/stable/api/plan7.html#pyhmmer.plan7.TopHits)
  class.
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
- **efficient**: Using `pyhmmer` to launch `hmmsearch` on sequences
  and HMMs in disk storage is typically faster than directly using the
  `hmmsearch` binary (see the [Benchmarks section](#%EF%B8%8F-benchmarks)).
  `pyhmmer.hmmsearch` uses a different parallelisation strategy compared to
  the `hmmsearch` binary from HMMER, which helps getting the most of
  multiple CPUs.

*This library is still a work-in-progress, and in an experimental stage,
but it should already pack enough features to run biological analyses
involving `hmmsearch` or `phmmer`.*


## 🔧 Installing

`pyhmmer` can be installed from [PyPI](https://pypi.org/project/pyhmmer/),
which hosts some pre-built CPython wheels for x86-64 Linux, as well as the
code required to compile from source with Cython:
```console
$ pip install pyhmmer
```

Compilation for UNIX PowerPC is not tested in CI, but should work out of the
box. Other architectures (e.g. Arm) and OSes (e.g. Windows) are not
supported by HMMER.

A [Bioconda](https://bioconda.github.io/) package is also available,
but only for Linux:
```console
$ conda install -c bioconda pyhmmer
```


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

Use `pyhmmer` to run `hmmsearch`, and obtain an iterable over
[`TopHits`](https://pyhmmer.readthedocs.io/en/stable/api/plan7.html#pyhmmer.plan7.TopHits)
that can be used for further sorting/querying in Python:

```python
import pyhmmer

with pyhmmer.easel.SequenceFile("938293.PRJEB85.HG003687.faa") as file:
    alphabet = file.guess_alphabet()
    sequences = [seq.digitize(alphabet) for seq in file]

with pyhmmer.plan7.HMMFile("Pfam.hmm") as hmms:
    all_hits = list(pyhmmer.hmmsearch(hmms, sequences_file, cpus=4))
```

Processing happens in parallel using Python threads, and a ``TopHits``
object is yielded for every ``HMM`` passed in the input iterable.


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pyhmmer/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See [`CONTRIBUTING.md`](https://github.com/althonos/pyhmmer/blob/master/CONTRIBUTING.md) for more details.

## ⏱️ Benchmarks

Benchmarks were run on a [i7-10710U CPU](https://ark.intel.com/content/www/us/en/ark/products/196448/intel-core-i7-10710u-processor-12m-cache-up-to-4-70-ghz.html) running 1.10GHz with 6 physical / 12
logical cores, using a FASTA file containing 2100 protein sequences extracted
from the genome of *Anaerococcus provencensis*
([`938293.PRJEB85.HG003687.faa`](https://github.com/althonos/pyhmmer/blob/master/tests/data/seqs/938293.PRJEB85.HG003687.faa))
and the version 33.1 of the [Pfam](https://pfam.xfam.org/) HMM library containing
18,259 domains. Commands were run 4 times on a warm SSD. *Plain lines show
the times for pressed HMMs, and dashed-lines the times for HMMs in text format.*

![Benchmarks](https://raw.github.com/althonos/pyhmmer/master/benches/benchmarks.svg)

Raw numbers can be found in the [`benches` folder](https://github.com/althonos/pyhmmer/blob/master/benches/).
They suggest that `phmmer` should be run with the number of *logical* cores,
while `hmmsearch` should be run with the number of *physical* cores (or less).
A possible explanation for this observation would be that HMMER
platform-specific code requires too many [SIMD](https://en.wikipedia.org/wiki/SIMD)
registers per thread to benefit from [simultaneous multi-threading](https://en.wikipedia.org/wiki/Simultaneous_multithreading).


## 🔍 See Also

If despite of all the advantages listed earlier, you would rather use HMMER through its CLI, 
this package will not be of great help. You should then check the 
[`hmmer-py`](https://github.com/EBI-Metagenomics/hmmer-py) package developed 
by [Danilo Horta](https://github.com/horta) at the [EMBL-EBI](https://www.ebi.ac.uk).


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The HMMER3 and Easel code is available under the
[BSD 3-clause](https://choosealicense.com/licenses/bsd-3-clause/) license.
See `vendor/hmmer/LICENSE` and `vendor/easel/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the [original HMMER authors](http://hmmer.org/). It was developed by
[Martin Larralde](https://github.com/althonos/pyhmmer) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
