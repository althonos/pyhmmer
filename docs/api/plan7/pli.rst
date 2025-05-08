Pipelines
=========

.. autoclass:: pyhmmer.plan7.Pipeline
   :special-members: __init__
   :members:
   :exclude-members: scan_seq, search_hmm, search_msa, search_seq

   .. automethod:: pyhmmer.plan7.Pipeline.scan_seq

   .. automethod:: pyhmmer.plan7.Pipeline.search_hmm

   .. automethod:: pyhmmer.plan7.Pipeline.search_msa

   .. automethod:: pyhmmer.plan7.Pipeline.search_seq


.. autoclass:: pyhmmer.plan7.LongTargetsPipeline(Pipeline)
   :special-members: __init__
   :members:
   :exclude-members: scan_seq, search_hmm, search_msa, search_seq

   .. automethod:: pyhmmer.plan7.LongTargetsPipeline.scan_seq

   .. automethod:: pyhmmer.plan7.LongTargetsPipeline.search_hmm

   .. automethod:: pyhmmer.plan7.LongTargetsPipeline.search_msa

   .. automethod:: pyhmmer.plan7.LongTargetsPipeline.search_seq


.. autoclass:: pyhmmer.plan7.Builder
   :special-members: __init__
   :members:

.. autoclass:: pyhmmer.plan7.Background
   :special-members: __init__
   :members: