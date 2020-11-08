import configparser
import glob
import os
import sys
import shlex

import setuptools
from Cython.Build import cythonize
from distutils import log
from distutils.command.clean import clean as _clean
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library


class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    def build_extension(self, ext):
        #
        self.run_command("build_clib")

        # update the extension C flags to use the temporary build folder
        for lib in self.distribution.libraries:
            for include_dir in lib.include_dirs:
                ext.include_dirs.append(os.path.join(self.build_temp, include_dir))
            ext.library_dirs.append(os.path.join(self.build_temp, os.path.dirname(lib.depends[1])))

        # update compile flags if compiling in debug mode
        if self.debug:
            if sys.platform == "linux" or sys.platform == "darwin":
                ext.extra_compile_args.append("-O0")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)


class build_clib(_build_clib):
    """A custom `build_clib` that compiles out of source.
    """

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    def copy(self, libraries):
        # HARDCODED: copy the library sources to the build directory
        src_hmmer = os.path.join(os.path.dirname(__file__), "vendor", "hmmer")
        dst_hmmer = os.path.join(self.build_clib, "vendor", "hmmer")
        if not os.path.exists(dst_hmmer):
            self.copy_tree(src_hmmer, dst_hmmer, preserve_symlinks=1)

        src_easel = os.path.join(os.path.dirname(__file__), "vendor", "easel")
        dst_easel = os.path.join(self.build_clib, "vendor", "hmmer", "easel")
        if not os.path.exists(dst_easel):
            self.copy_tree(src_easel, dst_easel, preserve_symlinks=1)

    def configure(self, libraries):
        for lib in libraries:
            # update environ to invoke `./configure`, which will involve an
            # `os.chdir` eventually.
            _env = os.environ.copy()
            _cflags = " ".join(self.compiler.compiler[1:])
            _cwd = os.getcwd()

            try:
                # pass parameters from `self.compiler` to the `configure` script
                # using environment variables
                os.environ["AR"] = self.compiler.archiver[0]
                os.environ["CC"] = self.compiler.compiler[0]
                os.environ["CFLAGS"] = _cflags

                # chdir to the directory where to run autoconf
                build_dir = os.path.join(self.build_clib, os.path.dirname(lib.depends[0]))
                log.info("entering directory {!r}".format(build_dir))
                os.chdir(build_dir)

                # run autoconf to generate the configure script
                self.make_file(["configure.ac"], "configure", self.spawn, (["autoconf"],))

                # run the configure script
                configure_cmd = ["./configure", "--enable-pic", "--enable-threads"]
                if self.debug:
                    configure_cmd.append("--enable-debugging")
                self.make_file(["configure"], "Makefile", self.spawn, (configure_cmd,))

            finally:
                os.environ.update(_env)
                os.chdir(_cwd)

    def make(self, libraries):
        _env = os.environ.copy()

        try:
            if self.verbose:
                os.environ["V"] = "1"
            for lib in libraries:
                makedir = os.path.join(self.build_temp, os.path.dirname(lib.depends[1]))
                libname = self.compiler.library_filename(lib.name)
                command = ["make", "-C", makedir, libname]
                if self.force:
                    command.append("-B")
                self.spawn(command)
        finally:
            os.environ = _env

    def build_libraries(self, libraries):
        # copy sources and build libraries using autoconf/configure/make
        self.copy(libraries)
        self.configure(libraries)
        self.make(libraries)

        # setup the library dirs and include dirs used by `build_ext`
        self.library_dirs = self.include_dirs = [
            os.path.join(self.build_temp, os.path.dirname(lib.depends[1]))
            for lib in self.libraries
        ]


class clean(_clean):

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pyhmmer")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                log.info("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)


# --- C static libraries -----------------------------------------------------

libraries = [
    Library(
        "easel",
        sources=glob.glob("vendor/hmmer/easel/*.c"),
        include_dirs=["vendor/hmmer/easel"],
        depends=["vendor/hmmer/easel/configure.ac", "vendor/hmmer/easel/Makefile"],
    ),
    Library(
        "hmmer",
        sources=glob.glob("vendor/hmmer/src/*.c"),
        include_dirs=["vendor/hmmer/easel", "vendor/hmmer/src", "vendor/hmmer/libdivsufsort"],
        depends=["vendor/hmmer/configure.ac", "vendor/hmmer/src/Makefile"],
    ),
]


# --- Cython extensions ------------------------------------------------------

extensions = [
    Extension(
        "pyhmmer.errors",
        ["pyhmmer/errors.pyx"],
        libraries=["easel"],
    ),
    Extension(
        "pyhmmer.easel",
        ["pyhmmer/easel.pyx"],
        libraries=["easel"],
    ),
    Extension(
        "pyhmmer.plan7",
        ["pyhmmer/plan7.pyx"],
        libraries=["hmmer", "easel"],
    ),
]


# --- Setup ------------------------------------------------------------------

setuptools.setup(
    ext_modules=cythonize(extensions, annotate=True, include_path=["include"]),
    libraries=libraries,
    cmdclass=dict(sdist=sdist, build_ext=build_ext, clean=clean, build_clib=build_clib),
)
