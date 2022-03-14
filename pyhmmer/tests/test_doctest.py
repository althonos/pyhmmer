# coding: utf-8
"""Test doctest contained tests in every file of the module.
"""

import configparser
import doctest
import importlib
import os
import pkgutil
import re
import sys
import shutil
import types
import warnings
from unittest import mock

try:
    import numpy
except ImportError:
    numpy = None

import pyhmmer.easel
import pyhmmer.plan7


def _load_tests_from_module(tests, module, globs, setUp=None, tearDown=None):
    """Load tests from module, iterating through submodules.
    """
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith("_")):
        if isinstance(attr, types.ModuleType):
            suite = doctest.DocTestSuite(
                attr,
                globs,
                setUp=setUp,
                tearDown=tearDown,
                optionflags=+doctest.ELLIPSIS,
            )
            tests.addTests(suite)
    return tests


def load_tests(loader, tests, ignore):
    """`load_test` function used by unittest to find the doctests.
    """
    _current_cwd = os.getcwd()

    def setUp(self):
        warnings.simplefilter("ignore")
        os.chdir(os.path.realpath(os.path.join(__file__, "..", "..")))

    def tearDown(self):
        os.chdir(_current_cwd)
        warnings.simplefilter(warnings.defaultaction)

    # doctests are not compatible with `green`, so we may want to bail out
    # early if `green` is running the tests
    if sys.argv[0].endswith("green"):
        return tests

    # doctests require `numpy` to run, which may not be available because
    # it is a pain to get to work out-of-the-box on OSX inside CI
    if numpy is None:
        return tests

    # add a sample HMM and some sequences to use with the globals
    data = os.path.realpath(os.path.join(__file__, os.pardir, "data"))
    hmm_path = os.path.join(data, "hmms", "txt", "Thioesterase.hmm")
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        thioesterase = next(hmm_file)
    seq_path = os.path.join(data, "seqs", "938293.PRJEB85.HG003687.faa")
    with pyhmmer.easel.SequenceFile(seq_path, digital=True, alphabet=thioesterase.alphabet) as seq_file:
        proteins = list(seq_file)
    msa_path = os.path.join(data, "msa", "LuxC.faa")
    with pyhmmer.easel.MSAFile(msa_path, "afa") as msa_file:
        luxc = next(msa_file)

    # recursively traverse all library submodules and load tests from them
    packages = [None, pyhmmer]

    for pkg in iter(packages.pop, None):
        for (_, subpkgname, subispkg) in pkgutil.walk_packages(pkg.__path__):
            # do not import __main__ module to avoid side effects!
            if subpkgname == "__main__" or subpkgname.startswith("tests"):
                continue
            # import the submodule and add it to the tests
            module = importlib.import_module(".".join([pkg.__name__, subpkgname]))
            globs = dict(
                numpy=numpy,
                easel=pyhmmer.easel,
                plan7=pyhmmer.plan7,
                thioesterase=thioesterase,
                proteins=proteins,
                luxc=luxc,
                pyhmmer=pyhmmer,
                **module.__dict__
            )
            tests.addTests(
                doctest.DocTestSuite(
                    module,
                    globs=globs,
                    setUp=setUp,
                    tearDown=tearDown,
                    optionflags=+doctest.ELLIPSIS,
                )
            )
            # if the submodule is a package, we need to process its submodules
            # as well, so we add it to the package queue
            if subispkg and subpkgname != "tests":
                packages.append(module)

    return tests
