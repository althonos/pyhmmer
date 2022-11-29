from . import (
    test_background,
    test_block,
    test_builder,
    test_hit,
    test_hmm,
    test_hmmfile,
    test_optimizedprofile,
    test_pipeline,
    test_profile,
    test_tophits,
    test_tracealigner,
    test_traces,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_background))
    suite.addTests(loader.loadTestsFromModule(test_block))
    suite.addTests(loader.loadTestsFromModule(test_builder))
    suite.addTests(loader.loadTestsFromModule(test_hit))
    suite.addTests(loader.loadTestsFromModule(test_hmm))
    suite.addTests(loader.loadTestsFromModule(test_hmmfile))
    suite.addTests(loader.loadTestsFromModule(test_optimizedprofile))
    suite.addTests(loader.loadTestsFromModule(test_pipeline))
    suite.addTests(loader.loadTestsFromModule(test_profile))
    suite.addTests(loader.loadTestsFromModule(test_tophits))
    suite.addTests(loader.loadTestsFromModule(test_tracealigner))
    suite.addTests(loader.loadTestsFromModule(test_traces))
    return suite
