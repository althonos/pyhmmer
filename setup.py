#!/usr/bin/env python
# coding: utf-8

import collections
import configparser
import functools
import glob
import multiprocessing.pool
import os
import platform
import re
import sys
import sysconfig
import subprocess
import tempfile
from pprint import pprint

import setuptools # always import setuptools first
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError, LinkError
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Utils ------------------------------------------------------------------

def _split_multiline(value):
    value = value.strip()
    sep = max('\n,;', key=value.count)
    return list(filter(None, map(lambda x: x.strip(), value.split(sep))))

def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

def _patch_osx_compiler(compiler):
    # On newer OSX, Python has been compiled as a universal binary, so
    # it will attempt to pass universal binary flags when building the
    # extension. This will not work because the code makes use of SSE2.
    for tool in ("compiler", "compiler_so", "linker_so"):
        flags = getattr(compiler, tool)
        i = next((i for i in range(1, len(flags)) if flags[i-1] == "-arch" and flags[i] == "arm64"), None)
        if i is not None:
            flags.pop(i)
            flags.pop(i-1)

def _detect_target_machine(platform):
    if platform == "win32":
        return "x86"
    return platform.rsplit("-", 1)[-1]

def _detect_target_cpu(platform):
    machine = _detect_target_machine(platform)
    if re.match("^mips", machine):
        return "mips"
    elif re.match("^(aarch64|arm64)$", machine):
        return "aarch64"
    elif re.match("^arm", machine):
        return "arm"
    elif re.match("(x86_64)|(x86)|(AMD64|amd64)|(^i.86$)", machine):
        return "x86"
    elif re.match("^(powerpc|ppc)", machine):
        return "ppc"
    return None

def _detect_target_system(platform):
    if platform.startswith("win"):
        return "windows"
    elif platform.startswith("macos"):
        return "macos"
    elif platform.startswith("linux"):
        return "linux_or_android"
    elif platform.startswith("freebsd"):
        return "freebsd"
    return None


# --- Library with platform-specific code ------------------------------------

class Library(setuptools.extension.Library):

    def __init__(self, *args, **kwargs):
        self._needs_stub = False
        self.platform_sources = kwargs.pop("platform_sources", {})
        self.platform_define_macros = kwargs.pop("platform_define_macros", {})
        self.platform_compile_args = kwargs.pop("platform_compile_args", {})
        super().__init__(*args, **kwargs)

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

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # detect if parallel build is enabled
        if self.parallel == 0:
            self.parallel = os.cpu_count()
        # transfer arguments to the build_clib method
        self._clib_cmd = self.distribution.get_command_obj("build_clib", True)
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.parallel = self.parallel
        self._clib_cmd.plat_name = self.plat_name
        self._clib_cmd.finalize_options()

    @property
    def target_machine(self):
        return self._clib_cmd.target_machine

    @property
    def target_system(self):
        return self._clib_cmd.target_system

    @property
    def target_cpu(self):
        return self._clib_cmd.target_cpu

    @property
    def hmmer_impl(self):
        return self._clib_cmd.hmmer_impl

    def _check_getid(self):
        _eprint('checking whether `PyInterpreterState_GetID` is available')

        self.mkpath(self.build_temp)
        fd, testfile = tempfile.mkstemp(prefix="have_getid", dir=self.build_temp, suffix=".c")
        objects = []

        with os.fdopen(fd, "w") as f:
            f.write("""
            #include <stdint.h>
            #include <stdlib.h>
            #include <Python.h>

            int main(int argc, char *argv[]) {{
                PyInterpreterState_GetID(NULL);
                return 0;
            }}
            """)

        if self.compiler.compiler_type == "msvc":
            flags = ["/WX"]
        else:
            flags = ["-Werror=implicit-function-declaration"]

        try:
            objects = self.compiler.compile([testfile], extra_postargs=flags)
        except CompileError:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            if os.path.exists(testfile):
                os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

    # --- Build code ---

    def run(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize
        # check the CPU architecture could be detected
        if self.target_machine is None:
            raise RuntimeError("Could not detect CPU architecture with `platform.machine`")
        # compile the C library if not done already
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # check a platform-specific implementation of HMMER was selected
        # depending on the detected machine
        if self._clib_cmd.hmmer_impl is None:
            raise RuntimeError('Could not select implementation for CPU architecture: "{}"'.format(machine))
        else:
            _eprint('Building HMMER with', self._clib_cmd.hmmer_impl, 'for CPU architecture:', repr(self.target_machine))

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "nthreads": self.parallel,
            "compiler_directives": {
                "binding": True,
                "linetrace": True,
                "embedsignature": True,
                "embedsignature.format": "clinic",
            },
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
                "SYS_BYTEORDER": sys.byteorder,
                "PLATFORM_UNAME_SYSTEM": platform.uname().system,
                "HMMER_IMPL": self._clib_cmd.hmmer_impl,
            }
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["gdb_debug"] = True
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
            cython_args["compiler_directives"]["profile"] = True
            cython_args["compiler_directives"]["nonecheck"] = True
        else:
            cython_args["compiler_directives"]["emit_code_comments"] = False
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False
            cython_args["compiler_directives"]["cdivision"] = True

        # cythonize and patch the extensions
        self.extensions = cythonize(self.extensions, **cython_args)
        for ext in self.extensions:
            ext._needs_stub = False

        # build the extensions as normal
        _build_ext.run(self)

    def build_extensions(self):
        # make sure the PyInterpreterState_GetID() function is available
        if self._check_getid():
            for ext in self.extensions:
                ext.define_macros.append(("HAS_PYINTERPRETERSTATE_GETID", 1))
        # build the extensions as normal
        _build_ext.build_extensions(self)

    def build_extension(self, ext):
        # show the compiler being used
        _eprint("building", ext.name, "for", self.plat_name, "with", self.compiler.compiler_type, "compiler")

        # setup HMMER implementation-specific flags
        self._clib_cmd._setup_impl(ext)
        # update compile flags if compiling in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-Og")
                ext.extra_compile_args.append("--coverage")
                ext.extra_link_args.append("--coverage")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Od")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-Wno-unused-variable")
        # remove universal binary CFLAGS from the compiler if any
        if self.target_system == "macos":
            _patch_osx_compiler(self.compiler)

        # update link and include directories
        ext.include_dirs.append(self._clib_cmd.build_clib)
        ext.library_dirs.append(self._clib_cmd.build_clib)
        for name in ext.libraries:
            lib = self._clib_cmd.get_library(name)
            ext.include_dirs.extend(lib.include_dirs)
            ext.extra_objects.append(self.compiler.library_filename(
                lib.name, output_dir=self._clib_cmd.build_clib
            ))

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)


class configure(_build_clib):
    """A `./configure`-style command to generate a platform-specific C header.

    Derives from the `build_clib` command from setuptools to be able to use
    the same configuration values (`self.build_temp`, `self.build_clib`, etc).
    """

    _FUNCTION_HEADERS = {
        "_mm_malloc": ["malloc.h"],
        "aligned_alloc": ["stdlib.h"],
        "chmod": ["sys/stat.h", "stddef.h"],
        "erfc": ["math.h"],
        "fstat": ["sys/stat.h", "stddef.h"],
        "fseeko": ["stdio.h"],
        "getcwd": ["unistd.h"],
        "getpid": ["unistd.h"],
        "mkstemp": ["stdlib.h"],
        "popen": ["stdio.h"],
        "posix_memalign": ["stdlib.h"],
        "putenv": ["stdlib.h"],
        "strcasecmp": ["strings.h", "stddef.h"],
        "stat": ["sys/stat.h", "stddef.h"],
        "strsep": ["string.h"],
        "sysconf": ["unistd.h"],
        "sysctl": ["sys/types.h", "sys/sysctl.h", "stddef.h"],
        "times": ["sys/times.h", "stddef.h"],
    }

    _FUNCTION_ARGUMENTS = {
        "_mm_malloc": ["0", "0"],
        "aligned_alloc": ["0", "0"],
        "chmod": ["NULL", "0"],
        "erfc": ["0"],
        "fstat": ["1", "NULL"],
        "fseeko": ["NULL", "0", "0"],
        "getcwd": ["NULL", "0"],
        "getpid": [],
        "mkstemp": ["NULL"],
        "popen": ["NULL", "NULL"],
        "posix_memalign": ["NULL", "0", "0"],
        "putenv": ["NULL"],
        "strcasecmp": ["NULL", "NULL"],
        "stat": ["NULL", "NULL"],
        "strsep": ["NULL", "NULL"],
        "sysconf": ["0"],
        "sysctl": ["NULL", "0", "NULL", "NULL", "NULL", "0"],
        "times": ["NULL"],
    }

    def initialize_options(self):
        _build_clib.initialize_options(self)
        self.target_machine = None
        self.target_system = None
        self.target_cpu = None
        self.plat_name = None
        self.hmmer_impl = None

    def finalize_options(self):
        _build_clib.finalize_options(self)
        # detect platform options
        if self.plat_name is None:
            self.plat_name = sysconfig.get_platform()
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
        # detect HMMER implementation
        if self.target_machine.startswith('ppc') and not self.target_machine.endswith('le'):
            self.hmmer_impl = "VMX"
        elif self.target_machine.startswith(("x86", "amd", "i386", "i686")):
            self.hmmer_impl = "SSE"
        elif self.target_machine.lower().startswith(("arm", "aarch")):
            self.hmmer_impl = "NEON"
        else:
            self.hmmer_impl = None

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Autotools-like helpers ---

    def _has_header(self, headername):
        _eprint('checking whether <{}> can be included'.format(headername), end="... ")

        slug = re.sub("[./-]", "_", headername)
        testfile = os.path.join(self.build_temp, "have_{}.c".format(slug))
        objects = []

        with open(testfile, "w") as f:
            f.write('#include "{}"\n'.format(headername))
        try:
            objects = self.compiler.compile([testfile], debug=self.debug)
        except (CompileError, LinkError) as err:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

    def _has_function(self, funcname):
        _eprint('checking whether function', repr(funcname), 'is available', end="... ")

        testfile = os.path.join(self.build_temp, "have_{}.c".format(funcname))
        binfile = os.path.join(self.build_temp, "have_{}.bin".format(funcname))
        objects = []

        with open(testfile, "w") as f:
            for header in self._FUNCTION_HEADERS[funcname]:
                f.write("#include <{}>\n".format(header))
            f.write("int main() {{ {}({}); return 0; }}".format(
                funcname,
                ', '.join(self._FUNCTION_ARGUMENTS[funcname]),
            ))
        try:
            objects = self.compiler.compile([testfile], debug=self.debug, extra_postargs=["-Wno-all"])
            self.compiler.link_executable(objects, binfile)
        except (CompileError, LinkError) as err:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _check_simd_generic(self, name, program, platform_args=()):
        _eprint('checking whether compiler can build', name, 'code', end="... ")

        base = "have_{}".format(name)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write(program)

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_preargs=platform_args)
            self.compiler.link_executable(objects, base, extra_preargs=platform_args, output_dir=self.build_temp)
            subprocess.run([binfile], check=True)
        except CompileError:
            _eprint("no")
            return False
        except (subprocess.SubprocessError, OSError):
            _eprint("yes, but cannot run code")
            return True  # assume we are cross-compiling, and still build
        else:
            if not platform_args:
                _eprint("yes")
            else:
                _eprint("yes, with {}".format(" ".join(platform_args)))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _check_neon(self):
        return self._check_simd_generic(
            "NEON",
            program="""
                #include <arm_neon.h>
                int main(int argc, char *argv[]) {{
                    int16x8_t a = vdupq_n_s16(-1);
                              a = vabsq_s16(a);
                    short     x = vgetq_lane_s16(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """,
            platform_args=[] if "64" in self.target_machine else ["-mfpu=neon"]
        )

    def _check_sse2(self):
        return self._check_simd_generic(
            "SSE2",
            program="""
                #include <emmintrin.h>
                int main(int argc, char *argv[]) {{
                    __m128i a = _mm_set1_epi16(-1);
                            a = _mm_and_si128(a, a);
                    short   x = _mm_extract_epi16(a, 1);
                    return (x == -1) ? 0 : 1;
                }}
            """,
        )

    def _check_vmx(self):
        return self._check_simd_generic(
            "VMX", 
            program="""
                #include <altivec.h>
                int main() {
                    vector float a = vec_splats(1.0);
                    float f;
                    vec_ste(a, 0, &f);
                    return (f == 1) ? 0 : 1;
                }
            """
        )

    # --- Required interface for `setuptools.Command` ---

    def build_libraries(self, libraries):
        # read `setup.cfg`, which should be next to this file, so we can
        # use the `__file__` magic variable to locate it
        _cfg = configparser.ConfigParser()
        _cfg.read([__file__.replace(".py", ".cfg")])

        # ensure the output directory exists, otherwise create it
        self.mkpath(self.build_clib)

        # remove universal binary CFLAGS from the compiler if any
        if self.target_system == "macos":
            _patch_osx_compiler(self.compiler)

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

        # check platform flags
        if self.hmmer_impl is not None:
            if self.hmmer_impl == "SSE":
                supported_feature = self._check_sse2()
            elif self.hmmer_impl == "VMX":
                supported_feature = self._check_vmx()
            elif self.hmmer_impl == "NEON":
                supported_feature = self._check_neon()
            else:
                supported_feature = False
            if not supported_feature:
                raise RuntimeError("failed to compile platform-specific code, aborting.")

        # fill the defines if headers are found
        headers = headers or []
        for header in headers:
            if self._has_header(header):
                slug = re.sub("[./-]", "_", header).upper()
                defines["HAVE_{}".format(slug)] = 1

        # fill the defines if functions are found
        functions = functions or []
        for func in functions:
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

    # --- Distutils command interface ---

    user_options = _build_clib.user_options + [
        ("parallel", "j", "number of parallel build jobs"),
    ]

    def initialize_options(self):
        _build_clib.initialize_options(self)
        self.parallel = None
        self.plat_name = None

    def finalize_options(self):
        _build_clib.finalize_options(self)
        self._configure_cmd = self.distribution.get_command_obj("configure", True)
        self._configure_cmd.debug = self.debug
        self._configure_cmd.verbose = self.verbose
        self._configure_cmd.force = self.force
        self._configure_cmd.plat_name = self.plat_name
        self._configure_cmd.finalize_options()
        # detect if parallel build is enabled
        if self.parallel is not None:
            self.parallel = int(self.parallel)
        if self.parallel == 0:
            self.parallel = os.cpu_count()

    @property
    def target_machine(self):
        return self._configure_cmd.target_machine

    @property
    def target_system(self):
        return self._configure_cmd.target_system

    @property
    def target_cpu(self):
        return self._configure_cmd.target_cpu

    @property
    def hmmer_impl(self):
        return self._configure_cmd.hmmer_impl
                
    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Helpers ---

    def get_library(self, name):
        return next(lib for lib in self.libraries if lib.name == name)

    # --- Build code ---

    def _setup_impl(self, library):
        if self.hmmer_impl == "VMX":
            library.define_macros.append(("eslENABLE_VMX", 1))
            library.extra_compile_args.append("-maltivec")
            library.extra_link_args.append("-maltivec")
        elif self.hmmer_impl == "SSE":
            library.define_macros.append(("eslENABLE_SSE", 1))
        elif self.hmmer_impl == "NEON":
            library.define_macros.append(("eslENABLE_NEON", 1))
            if "64" not in self.target_machine:
                library.extra_compile_args.append("-mfpu=neon")
                library.extra_link_args.append("-mfpu=neon")

    def _patch_easel(self, library):
        # patch the `esl_sqio_ascii.c` so we can use functions it otherwise
        # declares as `static`
        old = next(src for src in library.sources if src.endswith("esl_sqio_ascii.c"))
        new = os.path.join(self.build_temp, "esl_sqio_ascii.c")
        self.make_file([old], new, self.destatic, (old, new))
        library.sources.remove(old)
        library.sources.append(new)
        # add implementation-specific flags and definitions
        self._setup_impl(library)

    def _patch_hmmer(self, library):
        # patch the `p7_hmmfile.c` so that we can use functions it otherwise
        # declares as `static`
        old = next(src for src in library.sources if src.endswith("p7_hmmfile.c"))
        new = os.path.join(self.build_temp, "p7_hmmfile.c")
        self.make_file([old], new, self.destatic, (old, new))
        library.sources.remove(old)
        library.sources.append(new)
        # add implementation-specific sources
        impl_folder = "impl_{}".format(self.hmmer_impl.lower())
        library.sources.extend(glob.glob(os.path.join("vendor", "hmmer", "src", impl_folder, "*.c")))
        library.sources.remove(os.path.join("vendor", "hmmer", "src", impl_folder, "vitscore.c"))
        # add implementation-specific flags and definitions
        self._setup_impl(library)

    def run(self):
        # make sure the C headers were generated already
        if not self.distribution.have_run.get("configure", False):
            self._configure_cmd.run()
        # build the libraries normally
        _build_clib.run(self)

    def build_libraries(self, libraries):
        # remove universal binary CFLAGS from the compiler if any
        if self.target_system == "macos":
            _patch_osx_compiler(self.compiler)
        # build extensions sequentially
        self.mkpath(self.build_clib)
        for library in libraries:
            self.make_file(
                library.sources,
                self.compiler.library_filename(library.name, output_dir=self.build_clib),
                self.build_library,
                (library,)
            )

    def build_library(self, library):
        # update define macros
        if library.name == "easel":
            self._patch_easel(library)
        elif library.name == "hmmer":
            self._patch_hmmer(library)

        # update compile flags if compiling in debug or release mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-Og")
                library.extra_compile_args.append("--coverage")
                library.extra_link_args.append("--coverage")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Od")

        # store compile args
        compile_args = (
            self.build_temp,
            library.define_macros,
            library.include_dirs + [self.build_clib],
            self.debug,
            library.extra_compile_args,
            None,
            library.depends,
        )

        # manually prepare sources and get the names of object files
        sources = library.sources.copy()
        objects = [
            os.path.join(self.build_temp, s.replace(".c", self.compiler.obj_extension))
            for s in sources
        ]

        # compile outdated files in parallel
        with multiprocessing.pool.ThreadPool(self.parallel) as pool:
            pool.starmap(
                functools.partial(self._compile_file, compile_args=compile_args),
                zip(sources, objects)
            )

        # create a static library
        self.compiler.create_static_lib(
            objects,
            library.name,
            output_dir=self.build_clib,
            debug=self.debug,
        )

    def destatic(self, old, new):
        with open(old, "r") as f:
            lines = f.readlines()
        with open(new, "w") as f:
            f.writelines(l.replace("static ", "") for l in lines)

    def _compile_file(self, source, object, compile_args):
        self.make_file(
            [source],
            object,
            self.compiler.compile,
            ([source], *compile_args)
        )


class clean(_clean):

    def remove_file(self, filename):
        if os.path.exists(filename):
            _eprint("removing", repr(filename))
            os.remove(filename)
        else:
            _eprint(repr(filename), "does not exist -- can't clean it")

    def run(self):
        _clean.run(self)

        _build_cmd = self.get_finalized_command("build_ext")
        _build_cmd.inplace = True

        for ext in self.distribution.ext_modules:
            filename = _build_cmd.get_ext_filename(ext.name)
            if self.all:
                self.remove_file(filename)
            basename = _build_cmd.get_ext_fullname(ext.name).replace(".", os.path.sep)
            for ext in ["c", "html"]:
                filename = os.path.extsep.join([basename, ext])
                self.remove_file(filename)


# --- C static libraries -----------------------------------------------------

# fmt: off
libraries = [
    Library(
        "divsufsort",
        sources=[os.path.join("vendor", "hmmer", "libdivsufsort", "divsufsort.c")],
    ),
    Library(
        "easel",
        sources=glob.glob(os.path.join("vendor", "easel", "*.c")),
        include_dirs=[os.path.join("vendor", "easel")],
    ),
    Library(
        "hmmer",
        sources=[
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
        ],
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
    ),
    Extension(
        "pyhmmer.easel",
        [os.path.join("pyhmmer", "easel.pyx")],
        libraries=["easel"],
        depends=glob.glob(os.path.join("pyhmmer", "fileobj", "*.h")),
    ),
    Extension(
        "pyhmmer.plan7",
        [os.path.join("pyhmmer", "plan7.pyx")],
        libraries=["hmmer", "easel", "divsufsort"],
        depends=glob.glob(os.path.join("pyhmmer", "fileobj", "*.h")),
    ),
    Extension(
        "pyhmmer.daemon",
        [os.path.join("pyhmmer", "daemon.pyx")],
        libraries=["hmmer", "easel", "divsufsort"],
    )
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
