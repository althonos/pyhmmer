PyHMMER |Stars|
===============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyhmmer.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyhmmer/stargazers


*Cython bindings and Python interface to* `HMMER3 <http://hmmer.org/>`_.

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads| |Paper|


.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyhmmer/test.yml?branch=master&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyhmmer/actions

.. |GitLabCI| image:: https://img.shields.io/gitlab/pipeline/larralde/pyhmmer/master?gitlab_url=https%3A%2F%2Fgit.embl.de&logo=gitlab&style=flat-square&maxAge=600
   :target: https://git.embl.de/larralde/pyhmmer/-/pipelines

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyhmmer?logo=codecov&style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyhmmer/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyhmmer.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyhmmer

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyhmmer?ogo=anaconda&style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyhmmer

.. |AUR| image:: https://img.shields.io/aur/version/python-pyhmmer?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyhmmer

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyhmmer?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyhmmer/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyhmmer.svg?logo=python&style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyhmmer/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyhmmer.svg?logo=python&style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyhmmer/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyhmmer/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pyhmmer/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyhmmer.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyhmmer/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyhmmer?style=flat-square&maxAge=3600
   :target: http://pyhmmer.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyhmmer/blob/master/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpyhmmer
   :target: https://pepy.tech/project/pyhmmer

.. |Paper| image:: https://img.shields.io/badge/paper-Bioinformatics-teal.svg?style=flat-square&maxAge=3600
   :target: https://doi.org/10.1093/bioinformatics/btad214


Overview
--------

HMMER is a biological sequence analysis tool that uses profile hidden Markov
models to search for sequence homologs. HMMER3 is maintained by members of the
the `Eddy/Rivas Laboratory <http://eddylab.org/>`_ at Harvard University.

``pyhmmer`` is a Python module, implemented using the `Cython <https://cython.org/>`_
language, that provides bindings to HMMER3. It directly interacts with the
HMMER internals, which has the following advantages over CLI wrappers:

- **single dependency**:
  If your software or your analysis pipeline is
  distributed as a Python package, you can add `pyhmmer` as a dependency to
  your project, and stop worrying about the HMMER binaries being properly
  setup on the end-user machine.
- **no intermediate files**:
  Everything happens in memory, in Python objects
  you have control on, making it easier to pass your inputs to HMMER without
  needing to write them to a temporary file. Output retrieval is also done
  in memory, via instances of the `pyhmmer.plan7.TopHits` class.
- **no input formatting**:
  The Easel object model is exposed in the `pyhmmer.easel` module, and you
  have the possibility to build a `~pyhmmer.easel.Sequence` object yourself to
  pass to the HMMER pipeline. This is useful if your sequences are already
  loaded in memory, for instance because you obtained them from another
  Python library (such as `Pyrodigal <https://github.com/althonos/pyrodigal>`_
  or `Biopython <https://biopython.org/>`_).
- **no output formatting**:
  HMMER3 is notorious for its numerous output files
  and its fixed-width tabular output, which is hard to parse (even
  `Bio.SearchIO.HmmerIO` is struggling on some sequences).
- **efficient**:
  Using `pyhmmer` to launch ``hmmsearch`` on sequences and HMMs in disk storage
  is typically faster than directly using the ``hmmsearch`` binary.
  `pyhmmer.hmmer.hmmsearch` uses a different parallelisation strategy compared to
  the ``hmmsearch`` binary from HMMER, which helps getting the most of
  multiple CPUs.


Setup
-----

Run ``pip install pyhmmer`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <install>` to find other ways to install ``pyhmmer``.


Citation
--------

PyHMMER is scientific software, with a
`published paper <https://doi.org/10.1093/bioinformatics/btad214>`_
in the `Bioinformatics <https://academic.oup.com/bioinformatics>`_. Check the
:doc:`Publications page <publications>` to see how to cite Pyrodigal PyHMMER.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Examples <examples/index>
   Performance <performance>
   Contributing <contributing>
   Publications <publications>
   Benchmarks <benchmarks>
   API Reference <api/index>
   Changelog <changes>


Related Projects
----------------

Building a HMM from scratch? Then you may be interested in the `PyFAMSA <https://pypi.org/project/pyfamsa/>`_
package, providing bindings to `FAMSA <https://github.com/refresh-bio/FAMSA>`_,
a very fast multiple sequence aligner. In addition, you may want to trim alignments:
in that case, consider `PytrimAl <https://pypi.org/project/pytrimal>`_, which
wraps `trimAl 2.0 <https://github.com/inab/trimal/tree/2.0_RC>`_.

If despite of all the advantages listed earlier, you would rather use HMMER
through its CLI, this package will not be of great help. You can instead check
the `hmmer-py <https://github.com/EBI-Metagenomics/hmmer-py>`_ package developed
by `Danilo Horta <https://github.com/horta>`_ at the `EMBL-EBI <https://www.ebi.ac.uk>`_.


License
-------

This library is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.
The HMMER3 and Easel code is available under
the `BSD 3-clause <https://choosealicense.com/licenses/bsd-3-clause/>`_ license,
which allows redistribution of their sources in the ``pyhmmer`` distribution.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `HMMER <http://hmmer.org>`_ *authors. It was developed by*
`Martin Larralde <https://github.com/althonos/pyhmmer>`_ *during his PhD project
at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
