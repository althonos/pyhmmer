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
        for lib, info in self.distribution.libraries:
            ext.library_dirs.append(os.path.join(self.build_temp, info["makedir"]))
            ext.include_dirs.append(os.path.join(self.build_temp, info["makedir"]))

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

    def copy(self, libraries):
        # copy the library sources to the build directory
        for libname, info in self.libraries:
            src = os.path.join(os.path.dirname(__file__), info["vendordir"])
            dst = os.path.join(self.build_temp, info["builddir"])
            if not os.path.exists(os.path.join(self.build_temp, info["makedir"])):
                self.copy_tree(src, dst, preserve_symlinks=1)

    def configure(self, libraries):
        for lib, info in libraries:
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
                build_dir = os.path.join(self.build_temp, info["builddir"])
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
            for lib, info in libraries:
                makedir = os.path.join(self.build_temp, info["makedir"])
                libname = self.compiler.library_filename(lib)
                self.spawn(["make", "-C", makedir, libname])
        finally:
            os.environ = _env

    def build_libraries(self, libraries):
        # copy sources and build libraries using autoconf/configure/make
        self.copy(libraries)
        self.configure(libraries)
        self.make(libraries)

        # setup the library dirs and include dirs used by `build_ext`
        self.library_dirs = self.include_dirs = [
            os.path.join(self.build_temp, info["makedir"])
            for lib, info in self.libraries
        ]

    def run(self):
        self.build_info = dict(self.libraries)
        _build_clib.run(self)


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


extensions = [
    setuptools.Extension(
        "pyhmmer.errors",
        ["pyhmmer/errors.pyx"],
        libraries=["easel"],
    ),
    setuptools.Extension(
        "pyhmmer.easel",
        ["pyhmmer/easel.pyx"],
        libraries=["easel"],
    ),
    setuptools.Extension(
        "pyhmmer.hmmer",
        ["pyhmmer/hmmer.pyx"],
        libraries=["hmmer", "easel"],
    ),
]

# the libraries dictionary to use with the `build_clib` command
# - sources (required by the base command): name of the autoconf template
# - vendordir: where the vendored sources are located
# - builddir: where to copy the vendored sources (with the temp folder as root)
# - makedir: the path where to invoke `make <staticlib>`
libraries = {
    "easel": {
        "sources": [],
        "vendordir": os.path.join("vendor", "easel"),
        "builddir": os.path.join("hmmer", "easel"),
        "makedir": os.path.join("hmmer", "easel"),
    },
    "hmmer": {
        "sources": [],
        "vendordir": os.path.join("vendor", "hmmer"),
        "builddir": "hmmer",
        "makedir": os.path.join("hmmer", "src"),
    },
}

setuptools.setup(
    ext_modules=cythonize(extensions, annotate=True, include_path=["include"]),
    libraries=list(libraries.items()),
    cmdclass=dict(sdist=sdist, build_ext=build_ext, clean=clean, build_clib=build_clib),
)
