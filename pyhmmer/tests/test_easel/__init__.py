from . import (
    test_alphabet,
    test_bitfield,
    test_block,
    test_keyhash,
    test_matrix,
    test_msa,
    test_msafile,
    test_randomness,
    test_sequence,
    test_sequencefile,
    test_ssi,
    test_vector,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_alphabet))
    suite.addTests(loader.loadTestsFromModule(test_bitfield))
    suite.addTests(loader.loadTestsFromModule(test_block))
    suite.addTests(loader.loadTestsFromModule(test_keyhash))
    suite.addTests(loader.loadTestsFromModule(test_matrix))
    suite.addTests(loader.loadTestsFromModule(test_msa))
    suite.addTests(loader.loadTestsFromModule(test_msafile))
    suite.addTests(loader.loadTestsFromModule(test_randomness))
    suite.addTests(loader.loadTestsFromModule(test_sequence))
    suite.addTests(loader.loadTestsFromModule(test_sequencefile))
    suite.addTests(loader.loadTestsFromModule(test_ssi))
    suite.addTests(loader.loadTestsFromModule(test_vector))
    return suite
