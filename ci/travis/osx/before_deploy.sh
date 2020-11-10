#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh


# --- Using proper Python executable -----------------------------------------

log Activating pyenv
eval "$(pyenv init -)"
pyenv shell $(pyenv versions --bare)


# --- Build and audit wheel --------------------------------------------------

log Cleaning previous build files
$PYTHON setup.py clean --all

log Building extensions in release mode
$PYTHON setup.py build_clib
$PYTHON setup.py build_ext --inplace

log Testing extensions in release mode
$PYTHON -m unittest discover -v

log Building wheel with \`$PYTHON\`
$PYTHON setup.py sdist bdist_wheel

log Verifying distribution files with \`twine\`
twine check dist/*
