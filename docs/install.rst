Installation
============

.. note::

    Most platforms, such as Linux x86-64, OSX and Windows x86-64 provide
    precompiled wheels, but other less frequent platforms will require building
    the wheel yourself. Building ``pyhmmer`` involves compiling HMMER3 and Easel
    from source, which requires a C compiler to be available on the machine.


PyPi
^^^^

``pyhmmer`` is hosted on GitHub, but the easiest way to install it is to download
the latest release from its `PyPi repository <https://pypi.python.org/pypi/pyhmmer>`_.
It will install all dependencies then install ``pyhmmer`` either from a wheel if
one is available, or from source after compiling the Cython code :

.. code:: console

	$ pip install --user pyhmmer

.. Conda
.. ^^^^^
..
.. Pronto is also available as a `recipe <https://anaconda.org/bioconda/pyhmmer>`_
.. in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
.. use the `conda` installer:
..
.. .. code:: console
..
.. 	 $ conda install -c bioconda pyhmmer


GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

	$ pip install --user https://github.com/althonos/pyhmmer/archive/master.zip

Keep in mind this will install the development version of the library, so not
everything may work as expected compared to a stable versioned release.


GitHub + ``setuptools``
^^^^^^^^^^^^^^^^^^^^^^^

If you do not have ``pip`` installed, you can still clone the repository and
run the ``setup.py`` file manually, although you will need to install the
project dependencies yourself:

.. code:: console

	$ git clone https://github.com/althonos/pyhmmer
	$ cd pyhmmer
	$ python setup.py build_clib build_ext
	# python setup.py install
