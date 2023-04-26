HMMER
=====

.. currentmodule:: pyhmmer.hmmer


.. automodule:: pyhmmer.hmmer


hmmsearch
---------

.. autofunction:: pyhmmer.hmmer.hmmsearch(queries, sequences, cpus=0, callback=None, **options)


hmmscan
-------

.. autofunction:: pyhmmer.hmmer.hmmscan(queries, profiles, cpus=0, callback=None, background=None, **options)


phmmer
------

.. autofunction:: pyhmmer.hmmer.phmmer(queries, sequences, cpus=0, callback=None, builder=None, **options)


nhmmer
------

.. autofunction:: pyhmmer.hmmer.nhmmer(queries, sequences, cpus=0, callback=None, builder=None, **options)


hmmpress
--------

.. autofunction:: pyhmmer.hmmer.hmmpress(hmms, output)


hmmalign
--------

.. autofunction:: pyhmmer.hmmer.hmmalign(hmm, sequences, trim=False, digitize=False, all_consensus_cols=True)


jackhmmer
---------

.. autofunction:: pyhmmer.hmmer.jackhmmer(queries, sequences, cpus=0, allback=None, builder=None, max_iterations=5, select_hits=None, checkpoints=False, **options)
