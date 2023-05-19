Create a Python package with embedded HMMs
==========================================


Overview
--------

Let's suppose we are developing a Python package for detecting
and analyzing `thioredoxins <https://en.wikipedia.org/wiki/Thioredoxin>`_,
where PyHMMER is used to identify thioredoxin domains. We want to
distribute the package and all the required HMMs so that the end users
can simply use the functions without having to download additional data.

.. hint::

    If you're new to developing Python packages, you may want to have a look
    at the `Python Packaging User Guide <https://packaging.python.org/en/latest/>`_.

For this example, we will be using pre-made HMMs downloaded from
`InterPro <https://www.ebi.ac.uk/interpro>`_, but this would work
similarly if we were to use custom-made hmms.

.. caution::

    Make sure that you have the permission to redistribute the HMMs along
    with your code! For instance, `Pfam <https://pfam.xfam.org/>`_ is licensed under
    `Creative Commons Zero <https://creativecommons.org/publicdomain/zero/1.0/>`_
    which places the data in the public domain, but `SMART <http://smart.embl-heidelberg.de>`_
    is licensed by `EMBLEM <https://software.embl-em.de/software/18>`_ and
    does not allow redistribution.


Folder structure
----------------

Let's start with a sample project using `setuptools` as the build backend.
Given an example package we would have the following folder structure:

.. code::

    .
    ├── LICENSE
    ├── README.md
    ├── setup.cfg
    ├── setup.py
    ├── redox_detector/
    │   ├── __init__.py
    │   └── search.py
    └── tests/

The easiest way to store the HMMs is to have them right next to
the Python files that will be using them: for instance, suppose we want
to use the `TIGR01068 (thioredoxin) <https://www.ebi.ac.uk/interpro/entry/tigrfams/TIGR01068/>`_
HMM, we can simply download it and put it in the main module folder:

.. code::

    .
    ├── LICENSE
    ├── README.md
    ├── setup.cfg
    ├── setup.py
    ├── redox_detector/
    │   ├── __init__.py
    │   ├── search.py
    │   └── TIGR01068.hmm
    └── tests/


Loading data from Python
------------------------

Since Python 3.7 the standard library contains the `importlib.resources` module
which provides an interface for loading arbitrary package data.

For instance, we can write a function that takes a list of sequences and
return only the sequences that contain a thioredoxin domain, using the
internal HMM for finding hits:

.. code:: Python

    # search.py

    import importlib.resources
    from typing import Iterable

    from pyhmmer.plan7 import HMMFile
    from pyhmmer.easel import Bitfield, TextSequence
    from pyhmmer.hmmer import hmmsearch

    def filter_thioredoxins(sequences: List[str]):
        # turn the input sequences into DigitalSequence objects
        # (we use the index of the sequence as their name)
        digital_sequences = [
            TextSequence(sequence=seq, name=str(i).encode()).digitize(hmm.alphabet).
            for i, seq in enumerate(sequences)
        ]

        # use a bitmap to record which input sequences have had a hit
        is_thioredoxin = Bitfield.zeros(len(sequences))

        # search sequences using the embedded HMMs with `importlib.resources`
        # (using __name__ as the module name tells `open_binary` to
        # look in the same folder as the Python source file)
        with importlib.resources.open_binary(__name__, "TIGR01068.hmm") as src:
            with pyhmmer.plan7.HMMFile(src) as hmm_file:
                # run the search pipeline and get hits with E-value <= 1e-5
                for hits in hmmsearch(hmm_file, digital_sequences, E=1e-5):
                    for hit in hits.included:
                        is_thioredoxin[int(hit.name)] = True

        # return only the sequences that had at least one hit
        return [ seq for i, seq in enumerate(sequences) if is_thioredoxin[i] ]


.. hint::

    In this example we used only a single HMM inside the HMM file, however the
    code above would work even if the HMM file contained more than one HMM.


Distributing data
-----------------

Now that the data is ready and that the Python code knows how to load it,
all that is left is to make sure the data files are actually picked up by
`setuptools` in the distribution files.

Using the appropriate section in the ``setup.cfg`` file, we can instruct
`setuptools` to add any file with the ``.hmm`` extension to the distribution
files:

.. code:: ini

    [options.package_data]
    redox_detector = *.hmm

Now running ``python setup.py sdist bdist_wheel`` will produce source and
wheel distributions that will contain the HMM file. When these archives are
installed, the HMMs will be installed next to the code.

.. caution::

    PyPI has some size limits for uploaded files: by default, a single
    distribution cannot be larger than 100MB. You will be able to use this
    method a relatively small number of HMMs (up to a few thousands), but
    you won't be able to upload a large HMM database such as Pfam in
    its entirety.

