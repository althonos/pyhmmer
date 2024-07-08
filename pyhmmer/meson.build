# --- Check if `PyInterpreterState_GetID` is missing 00-------------------------
getid_code = (
  '''
  #include <stdint.h>
  #include <stdlib.h>
  #include <Python.h>

  int main(int argc, char *argv[]) {
      PyInterpreterState_GetID(NULL);
      return 0;
  }
  '''
)
has_getid = cc.compiles(
  getid_code,
  dependencies : py_dep,
  name :  'PyInterpreterState_GetID',
  args : ['-Werror=implicit-function-declaration'],
)

# --- Python extension ---------------------------------------------------------

sys_version_info_major = run_command([python, '-c', 'import sys; print(sys.version_info.major)' ], check : true).stdout().strip()
sys_version_info_minor = run_command([python, '-c', 'import sys; print(sys.version_info.minor)' ], check : true).stdout().strip()
sys_version_info_micro = run_command([python, '-c', 'import sys; print(sys.version_info.micro)' ], check : true).stdout().strip()
sys_implementation_name = run_command([python, '-c', 'import sys; print(sys.implementation.name)' ], check : true).stdout().strip()

cython_args = [
    '-ESYS_VERSION_INFO_MAJOR=' + sys_version_info_major,
    '-ESYS_VERSION_INFO_MINOR=' + sys_version_info_minor,
    '-ESYS_VERSION_INFO_MICRO=' + sys_version_info_micro,
    '-ESYS_IMPLEMENTATION_NAME=' + sys_implementation_name,
    '-ESYS_BYTEORDER=' + target_machine.endian(),
    '-EPLATFORM_UNAME_SYSTEM=' + target_machine.system(),
    '-EHMMER_IMPL=SSE',
    '-Xcdivision=True',
    '-Xnonecheck=False',
]

if get_option('buildtype') == 'debug'
    cython_args += [
        '-Xcdivision_warnings=True',
        '-Xwarn.undeclared=True',
        '-Xwarn.unreachable=True',
        '-Xwarn.maybe_uninitialized=True',
        '-Xwarn.unused=True',
        '-Xwarn.unused_arg=True',
        '-Xwarn.unused_result=True',
        '-Xwarn.multiple_declarators=True',
    ]
else
    cython_args += [
        '-Xboundscheck=False',
        '-Xwraparound=False',
    ]
endif

if has_getid
    c_args = ['-DHAS_PYINTERPRETERSTATE_GETID']
else
    c_args = []
endif

pyhmmer_errors = python.extension_module('errors',
    'errors.pyx',
    dependencies : [py_dep, easel_dep],
    install: true,
    subdir : 'pyhmmer',
    cython_args : cython_args,
    c_args : c_args,
)

pyhmmer_easel = python.extension_module('easel',
    'easel.pyx',
    dependencies : [py_dep, easel_dep],
    install: true,
    subdir : 'pyhmmer',
    cython_args : cython_args,
    c_args : c_args,
)

pyhmmer_plan7 = python.extension_module('plan7',
    'plan7.pyx',
    dependencies : [py_dep, easel_dep, hmmer_dep],
    install: true,
    subdir : 'pyhmmer',
    cython_args : cython_args,
    c_args : c_args,
)

pyhmmer_daemon = python.extension_module('daemon',
    'daemon.pyx',
    dependencies : [py_dep, easel_dep, hmmer_dep],
    install: true,
    subdir : 'pyhmmer',
    cython_args : cython_args,
    c_args : c_args,
)

# --- Python sources -----------------------------------------------------------

python_sources = [
    '__init__.py',
    'daemon.pxd',
    'daemon.pyi',
    'easel.pxd',
    'easel.pyi',
    'errors.pyi',
    'hmmer.py',
    'plan7.pxd',
    'plan7.pyi',
    'utils.py'
]

python.install_sources(
    python_sources, 
    subdir : 'pyhmmer'
)

subdir('tests')