project(
  'pyhmmer',
  'c', 'cython',
  version : '0.10.13',
)

# Load C compiler
cc = meson.get_compiler('c')

# Find optional libraries
math = cc.find_library('m', required : false)
mpi = dependency('mpi', language : 'c', required : false)
threads = dependency('threads', required : false)

# load Python dependencies
py_mod = import('python')
python = py_mod.find_installation(pure: false)
py_dep = python.dependency()

# prepare a script to remove all C `static` keywords
destatic = generator(
  python,
  arguments: [ '@SOURCE_DIR@' / 'tools' / 'destatic.py', '--input', '@INPUT@', '--output', '@OUTPUT@' ],
  output : ['@PLAINNAME@']
)

subdir('lib')
# subdir('include')
subdir('pyhmmer')