#!/usr/bin/env python
# coding: utf-8

import collections
import configparser
import glob
import os
import platform
import re
import sys
import subprocess
from unittest import mock

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

# --- Constants --------------------------------------------------------------

PLATFORM_MACHINE = platform.machine()
SYS_IMPLEMENTATION = sys.implementation.name


# --- Utils ------------------------------------------------------------------

def _split_multiline(value):
    value = value.strip()
    sep = max('\n,;', key=value.count)
    return list(filter(None, map(lambda x: x.strip(), value.split(sep))))


# --- `setup.py` commands ----------------------------------------------------

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

        # check the CPU architecture could be detected
        global machine
        if machine is None:
            raise RuntimeError("Could not detect CPU architecture with `platform.machine`")

        # check a platform-specific implementation of HMMER was selected
        # depending on the detected machine
        global hmmer_impl
        if hmmer_impl is None:
            raise RuntimeError('Could not select implementation for CPU architecture: "{}"'.format(machine))
        else:
            log.info('Building HMMER with {} for CPU architecture: "{}"'.format(hmmer_impl, machine))

        # use debug directives with Cython if building in debug mode
        cython_args = {"include_path": ["include"], "compiler_directives": {}}
        cython_args["compile_time_env"] = {"SYS_IMPLEMENTATION": SYS_IMPLEMENTATION}
        if hmmer_impl is not None:
            cython_args["compile_time_env"]["HMMER_IMPL"] = hmmer_impl
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
        # so that the platform-specific headers and static libs can be found
        ext.include_dirs.append(_clib_cmd.build_clib)
        ext.library_dirs.append(_clib_cmd.build_clib)

        # update compile flags if compiling in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-O0")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Od")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)


class configure(_build_clib):
    """A `./configure`-style command to generate a platform-specific C header.

    Derives from the `build_clib` command from setuptools to be able to use
    the same configuration values (`self.build_temp`, `self.build_clib`, etc).
    """

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Autotools-like helpers ---

    def _silent_spawn(self, cmd):
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as err:
            raise CompileError(err.stderr)

    def _has_header(self, headername):
        slug = re.sub("[./-]", "_", headername)
        testfile = os.path.join(self.build_temp, "have_{}.c".format(slug))
        objects = []

        with open(testfile, "w") as f:
            f.write('#include "{}"\n'.format(headername))
        try:
            with mock.patch.object(self.compiler, "spawn", new=self._silent_spawn):
                objects = self.compiler.compile([testfile], debug=self.debug)
        except CompileError as err:
            log.warn('could not find header "{}"'.format(headername))
            return False
        else:
            log.debug('found header "{}"'.format(headername))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

    def _has_function(self, funcname):
        testfile = os.path.join(self.build_temp, "have_{}.c".format(funcname))
        binfile = os.path.join(self.build_temp, "have_{}.bin".format(funcname))
        objects = []

        with open(testfile, "w") as f:
            f.write('int main() {{return {}(); return 0;}}\n'.format(funcname))
        try:
            with mock.patch.object(self.compiler, "spawn", new=self._silent_spawn):
                objects = self.compiler.compile([testfile], debug=self.debug)
                self.compiler.link_executable(objects, binfile)
        except CompileError:
            log.warn('could not link to function "{}"'.format(funcname))
            return False
        else:
            log.debug('found function "{}"'.format(funcname))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _check_sse(self):
        pass

    def _check_vmx(self):
        pass

    # --- Required interface for `setuptools.Command` ---

    def build_libraries(self, libraries):
        # read `setup.cfg`, which should be next to this file, so we can
        # use the `__file__` magic variable to locate it
        _cfg = configparser.ConfigParser()
        _cfg.read([__file__.replace(".py", ".cfg")])

        # ensure the output directory exists, otherwise create it
        self.mkpath(self.build_clib)

        # run the `configure_library` method sequentially on each library,
        # unless the header already exists and the configuration has not
        # been edited
        for library in self.distribution.libraries:
            section = "configure.{}".format(library.name)
            config = dict(_cfg.items(section))

            headers = _split_multiline(config.get("headers", ""))
            functions = _split_multiline(config.get("functions", ""))
            constants = dict(
                map(str.strip, line.split("="))
                for line in _split_multiline(config.get("constants", ""))
            )

            self.make_file(
                [__file__.replace(".py", ".cfg")],
                os.path.join(self.build_clib, config["filename"]),
                self.configure_library,
                (
                    library,
                    config["filename"],
                    config.get("copy"),
                    constants,
                    headers,
                    functions,
                ),
                exec_msg='generating "{}" for {} library'.format(
                    config["filename"], library.name
                )
            )

    def configure_library(self, library, filename, copy=None, constants=None, headers=None, functions=None):
        # if we only need to copy the header (e.g. for divsufsort) then just
        # do that and then exit
        if copy is not None:
            self.copy_file(copy, os.path.join(self.build_clib, filename))
            return

        # create a mapping to store the defines, and make sure it is ordered
        # (this is to keep compatibility with Python 3.5, a dict would do fine)
        defines = collections.OrderedDict()

        # fill the defines with the constants
        constants = constants or {}
        for k, v in constants.items():
            defines[k] = v

        # check endianness
        if sys.byteorder == "big":
            defines["WORDS_BIGENDIAN"] = 1

        # fill the defines if headers are found
        headers = headers or []
        for header in headers:
            log.info('checking whether or not "{}" can be included'.format(header))
            if self._has_header(header):
                slug = re.sub("[./-]", "_", header).upper()
                defines["HAVE_{}".format(slug)] = 1

        # fill the defines if functions are found
        functions = functions or []
        for func in functions:
            log.info('checking whether or not function "{}" is available'.format(func))
            if self._has_function(func):
                defines["HAVE_{}".format(func.upper())] = 1


        # write the header file
        slug = re.sub("[./-]", "_", filename).upper()
        with open(os.path.join(self.build_clib, filename), "w") as f:
            f.write("#ifndef {}_INCLUDED\n".format(slug))
            f.write("#define {}_INCLUDED\n".format(slug))
            for k, v in defines.items():
                f.write("#define {} {}\n".format(k, '' if v is None else v))
            f.write("#endif\n".format(slug))


class build_clib(_build_clib):
    """A custom `build_clib` that compiles out of source.
    """

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Build code ---

    def run(self):
        # make sure the C headers were generated already
        self.run_command("configure")

        # patch the `p7_hmmfile.c` so that we can use functions it otherwise
        # declares as `static`
        libhmmer = next(lib for lib in self.libraries if lib.name == "hmmer")
        old = next( src for src in libhmmer.sources if src.endswith("p7_hmmfile.c") )
        new = os.path.join(self.build_temp, "p7_hmmfile.c")
        self.make_file([old], new, self.patch_p7_hmmfile, (old, new))
        libhmmer.sources.remove(old)
        libhmmer.sources.append(new)

        # build the libraries normally
        _build_clib.run(self)

    def build_libraries(self, libraries):
        self.mkpath(self.build_clib)
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
            macros=library.define_macros,
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

    def patch_p7_hmmfile(self, old, new):
        with open(old, "r") as f:
            lines = f.readlines()
        with open(new, "w") as f:
            f.writelines(l.replace("static ", "") for l in lines)


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

# fmt: off
hmmer_sources = [
    os.path.join("vendor", "hmmer", "src", basename)
    for basename in [
        "build.c", "cachedb.c", "cachedb_shard.c", "emit.c", "errors.c",
    	"evalues.c", "eweight.c", "generic_decoding.c", "generic_fwdback.c",
    	"generic_fwdback_chk.c", "generic_fwdback_banded.c", "generic_null2.c",
    	"generic_msv.c", "generic_optacc.c", "generic_stotrace.c",
        "generic_viterbi.c", "generic_vtrace.c", "h2_io.c", "heatmap.c",
    	"hmmlogo.c", "hmmdmstr.c", "hmmdmstr_shard.c", "hmmd_search_status.c",
        "hmmdwrkr.c", "hmmdwrkr_shard.c", "hmmdutils.c", "hmmer.c", "logsum.c",
        "modelconfig.c", "modelstats.c", "mpisupport.c", "seqmodel.c",
        "tracealign.c", "p7_alidisplay.c", "p7_bg.c", "p7_builder.c",
        "p7_domain.c", "p7_domaindef.c", "p7_gbands.c", "p7_gmx.c",
        "p7_gmxb.c", "p7_gmxchk.c", "p7_hit.c", "p7_hmm.c", "p7_hmmcache.c",
    	"p7_hmmd_search_stats.c", "p7_hmmfile.c", "p7_hmmwindow.c",
    	"p7_pipeline.c", "p7_prior.c", "p7_profile.c", "p7_spensemble.c",
    	"p7_tophits.c", "p7_trace.c", "p7_scoredata.c", "hmmpgmd2msa.c",
    	"fm_alphabet.c", "fm_general.c", "fm_sse.c", "fm_ssv.c",
    ]
]

# HMMER3 is only supported on x86 CPUs with SSE, and big endian PowerPC
# (see https://github.com/EddyRivasLab/hmmer/issues/142)
machine = platform.machine().lower()
if machine.startswith('ppc') and not machine.endswith('le'):
    hmmer_sources.extend(glob.glob(os.path.join("vendor", "hmmer", "src", "impl_vmx", "*.c")))
    hmmer_sources.remove(os.path.join("vendor", "hmmer", "src", "impl_vmx", "vitscore.c"))
    hmmer_impl = "VMX"
    platform_define_macros = [("eslENABLE_VMX", 1)]
    platform_compile_args = ["-maltivec"]
elif machine.startswith(("x86", "amd", "i386", "i686")):
    hmmer_sources.extend(glob.glob(os.path.join("vendor", "hmmer", "src", "impl_sse", "*.c")))
    hmmer_sources.remove(os.path.join("vendor", "hmmer", "src", "impl_sse", "vitscore.c"))
    hmmer_impl = "SSE"
    platform_define_macros = [("eslENABLE_SSE", 1)]
    platform_compile_args = ["-msse3"]
else:
    log.warn('pyHMMER is not supported on CPU architecture: "{}"'.format(machine))
    platform_define_macros = []
    platform_compile_args = []
    hmmer_impl = None

libraries = [
    Library(
        "divsufsort",
        sources=[os.path.join("vendor", "hmmer", "libdivsufsort", "divsufsort.c")],
    ),
    Library(
        "easel",
        sources=glob.glob(os.path.join("vendor", "easel", "*.c")),
        include_dirs=[os.path.join("vendor", "easel")],
        define_macros=platform_define_macros,
        extra_compile_args=platform_compile_args,
    ),
    Library(
        "hmmer",
        sources=hmmer_sources,
        extra_compile_args=platform_compile_args,
        define_macros=platform_define_macros,
        include_dirs=[
            os.path.join("vendor", "easel"),
            os.path.join("vendor", "hmmer", "src"),
            os.path.join("vendor", "hmmer", "libdivsufsort"),
        ],
    ),
]


# --- Cython extensions ------------------------------------------------------

extensions = [
    Extension(
        "pyhmmer.errors",
        [os.path.join("pyhmmer", "errors.pyx")],
        libraries=["easel"],
        include_dirs=[os.path.join("vendor", "easel")]
    ),
    Extension(
        "pyhmmer.easel",
        [os.path.join("pyhmmer", "easel.pyx")],
        libraries=["easel"],
        include_dirs=[os.path.join("vendor", "easel")],
        define_macros=platform_define_macros,
        extra_compile_args=platform_compile_args,
    ),
    Extension(
        "pyhmmer.plan7",
        [os.path.join("pyhmmer", "plan7.pyx")],
        libraries=["hmmer", "easel", "divsufsort"],
        define_macros=platform_define_macros,
        extra_compile_args=platform_compile_args,
        include_dirs=[
            os.path.join("vendor", "easel"),
            os.path.join("vendor", "hmmer", "src")
        ],
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
