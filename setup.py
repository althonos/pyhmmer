#!/usr/bin/env python
# coding: utf-8

import collections
import configparser
import glob
import os
import platform
import re
import sys
import shlex
import shutil

import setuptools
from distutils import log
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

try:
    import pycparser
except ImportError as err:
    pycparser = err


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

    def run(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {"include_path": ["include"], "compiler_directives": {}}
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False
            cython_args["compiler_directives"]["cdivision"] = True

        # cythonize the extensions
        self.extensions = cythonize(self.extensions, **cython_args)

        # patch the extensions (needed for `_build_ext.run` to work)
        for ext in self.extensions:
            ext._needs_stub = False

        # build the extensions as normal
        _build_ext.run(self)

    def build_extension(self, ext):
        # make sure the C libraries have been built already
        self.run_command("build_clib")
        _clib_cmd = self.get_finalized_command("build_clib")

        # update the extension C flags to use the temporary build folder
        ext.include_dirs.append(_clib_cmd.build_clib)
        ext.library_dirs.append(_clib_cmd.build_clib)

        # update compile flags if compiling in debug mode
        if self.debug:
            if sys.platform == "linux" or sys.platform == "darwin":
                ext.extra_compile_args.append("-O0")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)


class configure(setuptools.Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass


class build_clib(_build_clib):
    """A custom `build_clib` that compiles out of source.
    """

    # --- Autotools-like helpers

    def _has_type(self, typename, **kwargs):
        return True  # FIXME

    def _has_header(self, headername, **kwargs):

        slug = re.sub("[./-]", "_", headername)
        testfile = os.path.join(self.build_temp, "have_{}.c".format(slug))
        objects = []

        with open(testfile, "w") as f:
            f.write('#include "{}"\n'.format(headername))

        try:
            objects = self.compiler.compile(
                [testfile],
                output_dir=self.build_temp,
                debug=self.debug,
            )
            # self.compiler.compile(
            #     library.sources,
            #     output_dir=self.build_temp,
            #     include_dirs=library.include_dirs + [self.build_clib],
            #     debug=self.debug,
            #     depends=library.depends,
            #     extra_preargs=library.extra_compile_args,
            # )

        except CompileError as err:
            log.warn('could not find header "{}"'.format(headername))
            return False
        else:
            log.info('found header "{}"'.format(headername))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

    def _has_function(self, funcname, **kwargs):
        return True  # FIXME

    def _check_sse(self):
        pass

    def _check_vmx(self):
        pass

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Build code ---

    def build_libraries(self, libraries):

        self.mkpath(self.build_clib)
        self._configure_easel()
        self._configure_plan7()

        for library in libraries:
            self.make_file(
                library.sources,
                self.compiler.library_filename(library.name, output_dir=self.build_clib),
                self.build_library,
                (library,)
            )

    def build_library(self, library):
        objects = self.compiler.compile(
            library.sources,
            output_dir=self.build_temp,
            include_dirs=library.include_dirs + [self.build_clib],
            debug=self.debug,
            depends=library.depends,
            extra_preargs=library.extra_compile_args,
        )
        self.compiler.create_static_lib(
            objects,
            library.name,
            output_dir=self.build_clib,
            debug=self.debug,
        )

    def _configure_easel(self):

        defines = collections.OrderedDict()

        defines["EASEL_DATE"] = '"Jul 2020"'
        defines["EASEL_COPYRIGHT"] = '"Copyright (C) 2020 Howard Hughes Medical Institute"'
        defines["EASEL_LICENSE"] = '"Freely distributed under the BSD open source license."'
        defines["EASEL_VERSION"] = '"0.47"'
        defines["EASEL_URL"] = '"http://bioeasel.org/"'

        defines["eslSTOPWATCH_HIGHRES"] = ''
        defines["eslENABLE_SSE"] = ''

        #
        if self._has_function("aligned_alloc", includes=["stdlib.h"]):
            defines["HAVE_ALIGNED_ALLOC"] = 1
        if self._has_function("erfc", includes=["math.h"]):
            defines["HAVE_ERFC"] = 1
        if self._has_function("getpid", includes=["unistd.h"]):
            defines["HAVE_GETPID"] = 1
        # if self._has_function("_mm_malloc", includes=["xmmintrin.h"]):
        #     defines["HAVE__MM_MALLOC"] = 1
        if self._has_function("popen", includes=["stdio.h"]):
            defines["HAVE_POPEN"] = 1
        if self._has_function("posix_memalign", includes=["stdlib.h"]):
            defines["HAVE_POSIX_MEMALIGN"] = 1
        if self._has_function("strcasecmp", includes=["strings.h"]):
            defines["HAVE_STRCASECMP"] = 1
        if self._has_function("strsep", includes=["string.h"]):
            defines["HAVE_STRSEP"] = 1
        if self._has_function("sysconf", includes=["unistd.h"]):
            defines["HAVE_SYSCONF"] = 1
        # if self._has_function("sysctl"):
        #     defines["HAVE_SYSCTL"] = 1
        if self._has_function("times", includes=["sys/times.h"]):
            defines["HAVE_TIMES"] = 1


        #
        if self._has_header("endian.h"):
            defines["HAVE_ENDIAN_H"] = 1
        if self._has_header("inttypes.h"):
            defines["HAVE_INTTYPES_H"] = 1
        if self._has_header("stdint.h"):
            defines["HAVE_STDINT_H"] = 1
        if self._has_header("unistd.h"):
            defines["HAVE_UNISTD_H"] = 1
        if self._has_header("strings.h"):
            defines["HAVE_STRINGS_H"] = 1
        if self._has_header("netinet/in.h"):
            defines["HAVE_NETINET_IN_H"] = 1

        if self._has_header("sys/types.h"):
            defines["HAVE_SYS_TYPES_H"] = 1
        if self._has_header("sys/param.h"):
            defines["HAVE_SYS_PARAM_H"] = 1
        if self._has_header("sys/sysctl.h"):
            defines["HAVE_SYS_SYSCTL_H"] = 1

        with open(os.path.join(self.build_clib, "esl_config.h"), "w") as f:
            f.write("#ifndef eslCONFIG_INCLUDED\n")
            f.write("#define eslCONFIG_INCLUDED\n")
            for k, v in defines.items():
                f.write("#define {} {}\n".format(k, v))
            f.write("#endif  /*eslCONFIG_INCLUDED*/\n")

    def _configure_plan7(self):

        defines = collections.OrderedDict()

        defines["p7_MAXABET"] = 20
        defines["p7_MAXCODE"] = 29
        defines["p7_MAX_SC_TXTLEN"] = 11
        defines["p7_MAXDCHLET"] = 20
        defines["p7_SEQDBENV"] = '"BLASTDB"'
        defines["p7_HMMDBENV"] = '"PFAMDB"'

        defines["p7_ETARGET_AMINO"] = 0.59
        defines["p7_ETARGET_DNA"] = 0.62
        defines["p7_ETARGET_OTHER"] = 1.0

        defines["HMMER_VERSION"] = '"3.3.1"'
        defines["HMMER_DATE"] = '"Jul 2020"'
        defines["HMMER_LICENSE"] = '"Freely distributed under the BSD open source license."'
        defines["HMMER_URL"] = '"http://hmmer.org/"'
        defines["HMMER_COPYRIGHT"] = '"Copyright (C) 2020 Howard Hughes Medical Institute."'

        defines["eslENABLE_SSE"] = ''


        with open(os.path.join(self.build_clib, "p7_config.h"), "w") as f:
            f.write("#ifndef P7_CONFIGH_INCLUDED\n")
            f.write("#define P7_CONFIGH_INCLUDED\n")
            for k, v in defines.items():
                f.write("#define {} {}\n".format(k, v))
            f.write("#endif  /*P7_CONFIGH_INCLUDED*/\n")

        shutil.copy(
            "vendor/hmmer/libdivsufsort/divsufsort.h.in",
            os.path.join(self.build_clib, "divsufsort.h"),
        )


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

if platform.machine().startswith('ppc'):
    hmmer_platform_sources = glob.glob("vendor/hmmer/src/impl_vmx/*.c")
else:
    # just assume we are running on x86-64, but should be checked to be sure
    # (platform.machine doesn't work on Windows though)
    hmmer_platform_sources = glob.glob(os.path.join("vendor", "hmmer", "src", "impl_sse", "*.c"))
    hmmer_platform_sources.remove(os.path.join("vendor", "hmmer", "src", "impl_sse", "vitscore.c"))

libraries = [
    Library(
        "divsufsort",
        sources=["vendor/hmmer/libdivsufsort/divsufsort.c"],
    ),
    Library(
        "easel",
        sources=glob.glob("vendor/easel/*.c"),
        include_dirs=["vendor/easel"],
    ),
    Library(
        "hmmer",
        sources=glob.glob("vendor/hmmer/src/*.c") + hmmer_platform_sources,
        include_dirs=["vendor/easel", "vendor/hmmer/src", "vendor/hmmer/libdivsufsort"],
        extra_compile_args=["-msse3"]
    ),
]


# --- Cython extensions ------------------------------------------------------

extensions = [
    Extension(
        "pyhmmer.errors",
        ["pyhmmer/errors.pyx"],
        libraries=["easel"],
        include_dirs=["vendor/easel"]
    ),
    Extension(
        "pyhmmer.easel",
        ["pyhmmer/easel.pyx"],
        libraries=["easel"],
        include_dirs=["vendor/easel"]
    ),
    Extension(
        "pyhmmer.plan7",
        ["pyhmmer/plan7.pyx"],
        libraries=["hmmer", "easel", "divsufsort"],
        include_dirs=["vendor/easel", "vendor/hmmer/src"]
    ),
]


# --- Setup ------------------------------------------------------------------

setuptools.setup(
    ext_modules=extensions,
    libraries=libraries,
    cmdclass=dict(
        build_ext=build_ext,
        build_clib=build_clib,
        clean=clean,
        configure=configure,
        sdist=sdist,
    ),
)
