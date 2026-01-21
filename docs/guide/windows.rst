Windows Support
===============

Since ``v0.12.0``, PyHMMER is provided for Windows on an experimental basis.

.. caution::

    Note that PyHMMER on Windows has only been tested on the internal PyHMMER
    test suite, and not extensively for result consistency across a wide number
    of platforms or build combinations. PyHMMER on Windows is provided for
    convenience and is not officially supported by HMMER and its developers.


By default, HMMER is not compatible with Windows because it uses POSIX
functions extensively: `POSIX threads <https://man7.org/linux/man-pages/man7/pthreads.7.html>`_
for parallelism, `GNU autotools <https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html>`_
for building from source, or POSIX-specific system headers.


To build PyHMMER on Windows, we use `MinGW-w64 <https://www.mingw-w64.org/>`_
which provides a minimum GNU environment to build HMMER with ``gcc`` and
sufficient compatibility to handle most of the core features.


MinGW Build
-----------

PyHMMER was rewritten in ``v0.11.0`` to use ``scikit-build-core`` and CMake
to handle the build of the native extension. As CMake supports Windows
natively, all that is needed for cross-compilation is a working MinGW-w64
install, and pointing CMake to use ``gcc.exe``. This allows us to release
wheels for Windows AMD64 using ``cibuildwheel`` and MinGW-w64 on GitHub
Actions.


Easel and HMMER patches
-----------------------

Some core areas that required patching for Windows/MinGW support include:

    ``esl_buffer``
        The implementation for buffered readers in Easel, which internally
        supports data read from an in-memory location, from a ``FILE*``,
        or from a path using memory-mapping with ``mmap``. Some include
        guards were adjusted to ensure ``mmap`` is not used on a
        ``mmap``-deprived target.

    ``esl_exception``
        On failure, Easel usually attempts to write to the system logger
        (using ``vsyslog``). Since MinGW does not implement that interface
        we simply silence errors reaching ``esl_exception``. In practice,
        PyHMMER should catch them at a higher level and translate them to
        a Python exception with a custom exception handler, although some
        errors may be uncaught and cause a system exit.

    ``ctime_r``
        `~pyhmmer.plan7.HMM.creation_time` code uses ``ctime_r`` to format the
        HMM creation time, which is not a POSIX function. Windows provides
        similar functionality in the ``ctime_s`` function, for which we added
        include guards and conditional compilation.

    ``impl_sse``
        Several parts of the HMMER SSE implementation require manual
        memory alignment for the dynamic programming matrices, and incidentally
        used ``unsigned long`` to cast pointers for unsigned arithmetic.
        Unfortunately, on Windows 64-bit, ``long`` is 32-bit and only ``long long``
        is 64-bit, which caused pointers to be truncated. A patch was added
        to use the ``uintptr_t`` type for casting pointers, which is exposed
        in ``stdint.h`` for that exact purpose.


Memory Management
-----------------

Windows complicates memory management because each DLL handles its own
memory management (see `documentation <https://learn.microsoft.com/en-us/windows/win32/dlls/about-dynamic-link-libraries>`_).
Memory allocated within a DLL must be deallocated within that same DLL
because of virtual address space management. This required some refactoring
since PyHMMER was reallocating routinely in different setters of classes
wrapping a corresponding Easel or HMMER class (such as `~pyhmmer.easel.Sequence`
or `~pyhmmer.plan7.HMM`).

Where applicable, we now use ``esl_strdup`` from Easel to allocate string data
for Easel dataclasses, and corresponding ``esl_free`` when data originating
from Easel must be manually deallocated in PyHMMER. There may still be
locations where "illegal" allocation/deallocation happens, although most of
them should have been tracked and eliminated with the help of the test suite.


File-like object handling
-------------------------

Several I/O classes of PyHMMER, such as `~pyhmmer.easel.SequenceFile`, `~pyhmmer.easel.MSAFile`
or `~pyhmmer.plan7.MSAFile` support reading data from any Python file-like object,
and not just from a file located on the filesystem. This is handled through
OS specific code, where ``fopen_cookie`` is used on Linux and ``funopen`` is
used on BSD and derivatives (MacOS). Unfortunately, such an API is not exposed
by MinGW.

To emulate file-like object support, we use Windows' anonymous pipes, and then
use ``_open_osfhandle`` to get a ``FILE*`` from one end of the handle, while
wrapping the other end inside a background Python thread which manages the
piping accordingly. This has *high* risks of creating deadlocks if an I/O
operation happens in the main thread without having released the GIL.