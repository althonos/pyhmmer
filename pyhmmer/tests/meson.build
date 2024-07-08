python_sources = [
    '__init__.py',
    'requirements.txt',
    'test_doctest.py',
    'test_errors.py',
    'test_hmmer.py',
    'utils.py',
]

python.install_sources(
    python_sources,
    subdir : 'pyhmmer' / 'tests',
)

subdir('data')
subdir('test_easel')
subdir('test_plan7')