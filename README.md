# 🐍🟡♦️🟦 PyHMMER [![Stars](https://img.shields.io/github/stars/althonos/pyhmmer.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyhmmer/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [HMMER3](http://hmmer.org/).*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyhmmer/test.yml?branch=master&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyhmmer/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyhmmer?logo=codecov&style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyhmmer/)
[![PyPI](https://img.shields.io/pypi/v/pyhmmer.svg?logo=pypi&style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyhmmer?logo=anaconda&style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/pyhmmer)
[![AUR](https://img.shields.io/aur/version/python-pyhmmer?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyhmmer)
[![Wheel](https://img.shields.io/pypi/wheel/pyhmmer.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyhmmer.svg?logo=python&style=flat-square&maxAge=3600)](https://pypi.org/project/pyhmmer/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyhmmer.svg?logo=python&style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/pyhmmer/#files)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyhmmer/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pyhmmer/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyhmmer.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyhmmer/issues)
[![Docs](https://img.shields.io/readthedocs/pyhmmer/latest?style=flat-square&maxAge=600)](https://pyhmmer.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyhmmer/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyhmmer?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyhmmer)
[![Paper](https://img.shields.io/badge/paper-Bioinformatics-teal.svg?style=flat-square&maxAge=3600)](https://doi.org/10.1093/bioinformatics/btad214)
[![Citations](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fid%2Fpub.1157360482%2Fmetadata.json&query=%24.times_cited&style=flat-square&label=citations&maxAge=86400)](https://badge.dimensions.ai/details/id/pub.1157360482)


## 🗺️ Overview

HMMER is a biological sequence analysis tool that uses profile hidden Markov
models to search for sequence homologs. HMMER3 is developed and maintained by
the [Eddy/Rivas Laboratory](http://eddylab.org/) at Harvard University.

`pyhmmer` is a Python package, implemented using the [Cython](https://cython.org/)
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
- **no input formatting**: The Easel object model is exposed in the
  [`pyhmmer.easel`](https://pyhmmer.readthedocs.io/en/stable/api/easel.html)
  module, and you have the possibility to build a
  [`DigitalSequence`](https://pyhmmer.readthedocs.io/en/stable/api/easel.html#pyhmmer.easel.DigitalSequence)
  object yourself to pass to the HMMER pipeline. This is useful if your sequences are already
  loaded in memory, for instance because you obtained them from another
  Python library (such as [Pyrodigal](https://github.com/althonos/pyrodigal)
  or [Biopython](https://biopython.org/)).
- **no output formatting**: HMMER3 is notorious for its numerous output files
  and its fixed-width tabular output, which is hard to parse (even
  [`Bio.SearchIO.HmmerIO`](https://biopython.org/docs/dev/api/Bio.SearchIO.HmmerIO.html)
  is struggling on some sequences).
- **efficient**: Using `pyhmmer` to launch `hmmsearch` on sequences
  and HMMs in disk storage is typically as fast as directly using the
  `hmmsearch` binary (see the [Benchmarks section](#%EF%B8%8F-benchmarks)).
  [`pyhmmer.hmmer.hmmsearch`](https://pyhmmer.readthedocs.io/en/stable/api/hmmer.html#hmmsearch)
  uses a different parallelisation strategy compared to
  the `hmmsearch` binary from HMMER, which can help getting the most of
  multiple CPUs when annotating smaller sequence databases.

*This library is still a work-in-progress, and in an experimental stage,
but it should already pack enough features to run biological analyses or
workflows involving `hmmsearch`, `hmmscan`, `nhmmer`, `phmmer`, `hmmbuild`
and `hmmalign`.*


## 🔧 Installing

`pyhmmer` can be installed from [PyPI](https://pypi.org/project/pyhmmer/),
which hosts some pre-built CPython wheels for Linux and MacOS on x86-64 and Arm64, as well as the code required to compile from source with Cython:
```console
$ pip install pyhmmer
```

Compilation for UNIX PowerPC is not tested in CI, but should work out of the
box. Note than non-UNIX operating systems (such as Windows) are not
supported by HMMER.

A [Bioconda](https://bioconda.github.io/) package is also available:
```console
$ conda install -c bioconda pyhmmer
```

## 🔖 Citation

PyHMMER is scientific software, with a
[published paper](https://doi.org/10.1093/bioinformatics/btad214)
in the [Bioinformatics](https://academic.oup.com/bioinformatics). Please
cite both [PyHMMER](https://doi.org/10.21105/joss.04296)
and [HMMER](http://hmmer.org) if you are using it in
an academic work, for instance as:

> PyHMMER (Larralde *et al.*, 2023), a Python library binding to HMMER (Eddy, 2011).

Detailed references are available on the [Publications page](https://pyhmmer.readthedocs.io/en/stable/publications.html) of the
[online documentation](https://pyhmmer.readthedocs.io/).


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

Use `pyhmmer` to run `hmmsearch` to search for Type 2 PKS domains
([`t2pks.hmm`](https://raw.githubusercontent.com/althonos/pyhmmer/master/pyhmmer/tests/data/hmms/txt/t2pks.hmm))
inside proteins extracted from the genome of *Anaerococcus provencensis*
([`938293.PRJEB85.HG003687.faa`](https://raw.githubusercontent.com/althonos/pyhmmer/master/pyhmmer/tests/data/seqs/938293.PRJEB85.HG003687.faa)).
This will produce an iterable over
[`TopHits`] that can be used for further sorting/querying in Python.
Processing happens in parallel using Python threads, and a [`TopHits`]
object is yielded for every [`HMM`] passed in the input iterable.

[`HMM`]: https://pyhmmer.readthedocs.io/en/stable/api/plan7.html#pyhmmer.plan7.HMM
[`TopHits`]: https://pyhmmer.readthedocs.io/en/stable/api/plan7.html#pyhmmer.plan7.TopHits

```python
import pyhmmer

with pyhmmer.easel.SequenceFile("pyhmmer/tests/data/seqs/938293.PRJEB85.HG003687.faa", digital=True) as seq_file:
    sequences = list(seq_file)

with pyhmmer.plan7.HMMFile("pyhmmer/tests/data/hmms/txt/t2pks.hmm") as hmm_file:
    for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=4):
      print(f"HMM {hits.query_name.decode()} found {len(hits)} hits in the target sequences")
```

Have a look at more in-depth examples such as [building a HMM from an alignment](https://pyhmmer.readthedocs.io/en/stable/examples/msa_to_hmm.html),
[analysing the active site of a hit](https://pyhmmer.readthedocs.io/en/stable/examples/active_site.html),
or [fetching marker genes from a genome](https://pyhmmer.readthedocs.io/en/stable/examples/fetchmgs.html)
in the [Examples](https://pyhmmer.readthedocs.io/en/stable/examples/index.html)
page of the [online documentation](https://pyhmmer.readthedocs.io/).


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

Benchmarks were run on a [i7-10710U CPU](https://ark.intel.com/content/www/us/en/ark/products/196448/intel-core-i7-10710u-processor-12m-cache-up-to-4-70-ghz.html) running @1.10GHz with 6 physical / 12
logical cores, using a FASTA file containing 4,489 protein sequences extracted
from the genome of *Escherichia coli*
([`562.PRJEB4685`](https://progenomes.embl.de/genome.cgi))
and the version 33.1 of the [Pfam](https://pfam.xfam.org/) HMM library containing
18,259 domains. Commands were run 3 times on a warm SSD. *Plain lines show
the times for pressed HMMs, and dashed-lines the times for HMMs in text format.*

![Benchmarks](https://raw.github.com/althonos/pyhmmer/master/benches/v0.7.0/plot.svg)

Raw numbers can be found in the [`benches` folder](https://github.com/althonos/pyhmmer/blob/master/benches/).
They suggest that `phmmer` should be run with the number of *logical* cores,
while `hmmsearch` should be run with the number of *physical* cores (or less).
A possible explanation for this observation would be that HMMER
platform-specific code requires too many [SIMD](https://en.wikipedia.org/wiki/SIMD)
registers per thread to benefit from [simultaneous multi-threading](https://en.wikipedia.org/wiki/Simultaneous_multithreading).

To read more about how PyHMMER achieves better parallelism than HMMER for
many-to-many searches, have a look at the [Performance page](https://pyhmmer.readthedocs.io/en/stable/performance.html)
of the documentation.


## 🔍 See Also

Building a HMM from scratch? Then you may be interested in the [`pyfamsa`](https://pypi.org/project/pyfamsa/)
package, providing bindings to [FAMSA](https://github.com/refresh-bio/FAMSA),
a very fast multiple sequence aligner. In addition, you may want to trim alignments:
in that case, consider [`pytrimal`](https://pypi.org/project/pytrimal), which
wraps [trimAl 2.0](https://github.com/inab/trimal/tree/2.0_RC).

If despite of all the advantages listed earlier, you would rather use HMMER
through its CLI, this package will not be of great help. You can instead check
the [`hmmer-py`](https://github.com/EBI-Metagenomics/hmmer-py) package developed
by [Danilo Horta](https://github.com/horta) at the [EMBL-EBI](https://www.ebi.ac.uk).


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The HMMER3 and Easel code is available under the
[BSD 3-clause](https://choosealicense.com/licenses/bsd-3-clause/) license.
See `vendor/hmmer/LICENSE` and `vendor/easel/LICENSE` for more information.

*This project is in no way affiliated, sponsored, or otherwise endorsed by
the [original HMMER authors](http://hmmer.org/). It was developed by
[Martin Larralde](https://github.com/althonos/pyhmmer) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
