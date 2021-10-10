from . import (
    test_doctest,
    test_easel,
    test_errors,
    test_hmmer,
    test_plan7,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_easel))
    suite.addTests(loader.loadTestsFromModule(test_errors))
    suite.addTests(loader.loadTestsFromModule(test_hmmer))
    suite.addTests(loader.loadTestsFromModule(test_plan7))
    return suite
