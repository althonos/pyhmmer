PyHMMER |Stars|
===============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyhmmer.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyhmmer/stargazers


*Cython bindings and Python interface to* `HMMER3 <http://hmmer.org/>`_.

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads| |Paper| |Citations|


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

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyhmmer?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyhmmer

.. |Paper| image:: https://img.shields.io/badge/paper-Bioinformatics-teal.svg?style=flat-square&maxAge=3600
   :target: https://doi.org/10.1093/bioinformatics/btad214

.. |Citations| image:: https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fid%2Fpub.1157360482%2Fmetadata.json&query=%24.times_cited&style=flat-square&label=citations&maxAge=86400
   :target: https://badge.dimensions.ai/details/id/pub.1157360482


Overview
--------

HMMER is a biological sequence analysis tool that uses profile hidden Markov
models to search for sequence homologs. HMMER3 is maintained by members of the
the `Eddy/Rivas Laboratory <http://eddylab.org/>`_ at Harvard University.

``pyhmmer`` is a Python module, implemented using the `Cython <https://cython.org/>`_
language, that provides bindings to HMMER3. It directly interacts with the
HMMER internals, which has the following advantages over CLI wrappers:

.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pyhmmer`` as a ``pip`` or ``conda`` dependency, no need
      for the HMMER binaries or any external dependency.

   .. grid-item-card:: :fas:`screwdriver-wrench` Flexible

      Create input `~pyhmmer.easel.Sequence` and `~pyhmmer.plan7.HMM` objects
      with the :doc:`API <api/index>`, or load them from a file.

   .. grid-item-card:: :fas:`gears` Practical

      Retrieve nested results as dedicated `~pyhmmer.plan7.TopHits` objects,
      write them to a file, or use them for further Python analysis.

   .. grid-item-card:: :fas:`gauge-high` Fast

      Run `hmmsearch` in parallel using an efficient threading model, which
      :doc:`outperforms <guide/benchmarks>` HMMER in some typical usecases.

   .. grid-item-card:: :fas:`dolly` Shareable

      :doc:`Distribute and load <examples/embed_hmms>` `~pyhmmer.plan7.HMM`
      objects from inside a Python package to facilitate sharing analyses.

   .. grid-item-card:: :fas:`eye` Transparent

      Access the internals of a `~pyhmmer.plan7.HMM`, inspect the attributes
      and manually edit transitions or emissions scores.

Setup
-----

Run ``pip install pyhmmer`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <guide/install>` to find other ways to install ``pyhmmer``.


Citation
--------

PyHMMER is scientific software, with a
`published paper <https://doi.org/10.1093/bioinformatics/btad214>`_
in the `Bioinformatics <https://academic.oup.com/bioinformatics>`_. Check the
:doc:`Publications page <guide/publications>` to see how to cite PyHMMER.


Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   Examples <examples/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. grid:: 1 3 5 5
   :gutter: 1

   .. grid-item-card:: :fas:`fire` Pyrodigal
      :link: https://pyrodigal.readthedocs.io

      Gene Finding (with Prodigal).

   .. grid-item-card:: :fas:`virus-covid` Pyrodigal-gv
      :link: https://github.com/althonos/pyrodigal-gv

      Pyrodigal for Giant Viruses.

   .. grid-item-card:: :fas:`align-center` PyFAMSA
      :link: https://pyfamsa.readthedocs.io

      Multiple Sequence Alignment (with FAMSA).

   .. grid-item-card:: :fas:`scissors` PytrimAl
      :link: https://pytrimal.readthedocs.io

      Alignment Trimming (with trimAl).

   .. grid-item-card:: :fas:`music` LightMotif
      :link: https://lightmotif.readthedocs.io

      Platform-accelerated motif scoring.

   .. grid-item-card:: :fas:`knife;fa-custom` Diced
      :link: https://diced.readthedocs.io

      CRISPR Detection (with MinCED).

   .. grid-item-card:: :fas:`table-cells` Scoring Matrices
      :link: https://scoring-matrices.readthedocs.io

      Common scoring matrices for Cython.

   .. grid-item-card:: :fas:`chain` Pyskani
      :link: https://pyskani.readthedocs.io

      Average Nucleotide Identity (with skani).

   .. grid-item-card:: :fas:`forward-fast` PyFastANI
      :link: https://pyfastani.readthedocs.io

      Average Nucleotide Identity (with FastANI).

   .. grid-item-card:: :fas:`magnifying-glass` PyJess
      :link: https://pyjess.readthedocs.io

      Geometric Template Matching (with Jess).

   .. grid-item-card:: :fas:`repeat` PyTantan
      :link: https://pytantan.readthedocs.io

      Tandem Repeat Masking (with Tantan).

   .. grid-item-card:: :fas:`gem` PyOpal
      :link: https://pyopal.readthedocs.io

      Query/Database Aligner (with Opal).

   .. grid-item-card:: :fas:`sword;fa-custom` PySWRD
      :link: https://pyswrd.readthedocs.io

      Database Heuristic Filtering (with SWORD).

   .. grid-item-card:: :fas:`rocket` Mini3di
      :link: https://github.com/althonos/mini3di

      Protein structure to 3di (with FoldSeek).

   .. grid-item-card:: :fas:`calculator` ``peptides.py``
      :link: https://peptides.readthedocs.io

      Peptide descriptors.

   .. grid-item-card:: :fas:`diagram-project` Pronto
      :link: https://pronto.readthedocs.io

      Open Biomedical Ontologies in Python.

   .. grid-item-card:: :fas:`box` NAFcodec
      :link: https://nafcodec.readthedocs.io

      Nucleotide Archival Format in Python.

   .. grid-item-card:: :fas:`bank` ``gb-io.py``
      :link: https://gb-io.readthedocs.io

      Fast GenBank parser written in Rust.


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
