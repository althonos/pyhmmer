Installation
============

.. note::

    Wheels are provided for Linux and OSX x86-64 platforms, but other machines
    will have to build the wheel from the source distribution. Building
    ``pyhmmer`` involves compiling HMMER3 and Easel, which requires a C compiler
    to be available.


PyPi
^^^^

``pyhmmer`` is hosted on GitHub, but the easiest way to install it is to download
the latest release from its `PyPi repository <https://pypi.python.org/pypi/pyhmmer>`_.
It will install all dependencies then install ``pyhmmer`` either from a wheel if
one is available, or from source after compiling the Cython code :

.. code:: console

	$ pip install --user pyhmmer

Conda
^^^^^

Pronto is also available as a `recipe <https://anaconda.org/bioconda/pyhmmer>`_
in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
use the ``conda`` installer:

.. code:: console

	 $ conda install -c bioconda pyhmmer


Arch User Repository
^^^^^^^^^^^^^^^^^^^^

A package recipe for Arch Linux can be found in the Arch User Repository
under the name `python-pyhmmer <https://aur.archlinux.org/packages/python-pyhmmer>`_.
It will always match the latest release from PyPI.

Steps to install on ArchLinux depend on your `AUR helper <https://wiki.archlinux.org/title/AUR_helpers>`_
(``yaourt``, ``aura``, ``yay``, etc.). For ``aura``, you'll need to run:

.. code:: console

   $ aura -A python-pyhmmer


.. EMBL Package Registry
.. ^^^^^^^^^^^^^^^^^^^^^
..
.. You can also install ``manylinux`` wheels built from the latest commit that
.. passed the unit tests. Those bleeding-edge releases are available in the GitLab
.. Package Registry hosted on the EMBL ``git`` server. Just instruct ``pip`` to
.. use an extra index URL as follow:
..
.. .. code:: console
..
..   $ pip install --user pyhmmer --extra-index-url https://git.embl.de/api/v4/projects/3638/packages/pypi/simple


GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

	$ pip install --user https://github.com/althonos/pyhmmer/archive/master.zip

.. caution::

    Keep in mind this will install always try to install the latest commit,
    which may not even build, so consider using a versioned release instead.


GitHub + ``setuptools``
^^^^^^^^^^^^^^^^^^^^^^^

If you do not want to use ``pip``, you can still clone the repository and
run the ``setup.py`` file manually, although you will need to install the
build dependencies (mainly `Cython <https://pypi.org/project/cython>`_):

.. code:: console

	$ git clone --recursive https://github.com/althonos/pyhmmer
	$ cd pyhmmer
	$ python setup.py build
	# python setup.py install

.. Danger::

    Installing packages without ``pip`` is strongly discouraged, as they can
    only be uninstalled manually, and may damage your system.
