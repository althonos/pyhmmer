#!/bin/sh

set -e

export PATH="$HOME/.cargo/bin:$PATH"
export PYBIN="$(echo /opt/python/${1}*/bin)"
export PYTHON_SYS_EXECUTABLE="$PYBIN/python"
export PYTHON_LIB=$(${PYBIN}/python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
export LIBRARY_PATH="$LIBRARY_PATH:$PYTHON_LIB"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PYTHON_LIB"

cd /io

# Clean previous debug files
$PYTHON_SYS_EXECUTABLE setup.py clean --all

# Compile and test in release mode
$PYTHON_SYS_EXECUTABLE setup.py build_clib
$PYTHON_SYS_EXECUTABLE setup.py build_ext --inplace
$PYTHON_SYS_EXECUTABLE -m unittest discover -v

# Build wheels with release library
$PYTHON_SYS_EXECUTABLE setup.py sdist bdist_wheel

# Bundle external shared libraries into the wheels
for whl in /io/dist/*.whl; do
  auditwheel repair "$whl" -w /io/dist && rm $whl
done
