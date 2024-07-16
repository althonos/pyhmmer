Performance
===========

Background
----------

Benchmarks of PyHMMER conducted against the ``hmmsearch`` and ``hmmscan`` binaries
suggest that running a domain search pipeline takes about the same time in
single-threaded mode, and are faster when the right number of CPUs is used.

This comes from several changes in the implementation of the search pipeline
with the PyHMMER API compared to the original HMMER C code, both of which have
absolutely no effect on the final result.


Parallelisation strategy
------------------------

Both PyHMMER and HMMER support searching / scanning several targets with
several queries in parallel using multithreading. However, benchmarks suggest
that PyHMMER takes a greater advantage of the number of available CPUs.

Querying modes
^^^^^^^^^^^^^^

To get a better idea of how this is possible, let's first have a look at how
HMMER parallelises computations on a consumer setup (*i.e.* excluding MPI).
The ``hmmsearch`` binary can process multiple HMM queries against multiple
sequence targets in either of two modes:

- *serial mode*: An outer loop reads the HMMs one after the other and an inner loop
  reads each sequence. This uses a very small amount of memory (there is never
  more than a single sequence and a single HMM loaded).
  **This has the obvious drawback that the file containing the sequences must
  be read several times. The annotation can quickly become I/O-bound in that case.**
- *threaded mode*: An outer loop reads the HMMs one after the other, and clones
  the profiles so that they are available to all worker threads. Sequences are
  read into blocks, which are kept in memory so that they don't have to be read
  again for a subsequent query. Sequence blocks are passed in queues to the worker
  threads, which annotate them in parallel. When all blocks have been processed,
  the results are merged, and the loop advances to the next HMM.
  **However, in HMMER, the block size is set at compile-time (1000 sequences),
  and no load balancing is done between the different worker threads. Consequently,
  parallelism is only achieved for a large number of targets. The few threads
  to receive a full-sized block will be the bottleneck, while the rest of
  the threads will idle**.

Although the threaded mode removes the potential I/O bottleneck, it only works for
a sufficiently large number of targets (:math:`1000 \times n_{cpus}`). To achieve
true parallelism, PyHMMER improves on the threaded mode by switching the worker
thread logic. Target sequences are pre-fetched in memory before looping
over the queries, and are passed by reference to all the worker threads. Each
worker then receives a HMM from the main thread, and process the entirety of
the targets with that HMM. **This allows the pipeline search loop to be truly
run in parallel for several HMMs, and doesn't require threads to wait for each
other before moving on to the next query**.


.. admonition:: Note

    Obviously, the PyHMMER parallelisation strategy will only work for multiple
    queries. But one main motivation to develop PyHMMER was to annotate protein
    sequences with a subset of the `Pfam <http://pfam.xfam.org/>`_ HMM library,
    which is why we benchmark this particular use case.

    The Python API still remains flexible enough for the original HMMER strategy
    to be implemented; if this was of interest, you should consider opening
    a feature request on the `issue tracker <https://github.com/althonos/pyhmmer/issues>`_.


Example
^^^^^^^

To check how well PyHMMER and HMMER3 handle parallelism on a real dataset,
we annotated proteins from representative genomes of the
`proGenomes <https://progenomes.embl.de/>`_ database with the
`Pfam <http://pfam.xfam.org/>`_ collection of HMMs. ``hmmsearch`` runs
are shown in green, ``pyhmmer`` runs in blue.

.. image:: /_images/progenomes-v0.4.5.svg

*The proteins of each genomes were extracted into one FASTA file per genome,
and we measured the runtime of either the* ``hmmsearch`` *binary or the*
``python -m pyhmmer.hmmer hmmsearch`` *module on each file*.


Memory allocation
-----------------

PyHMMER is slightly more conservative with memory: in several places where
the original HMMER binary would reallocate memory within loops, PyHMMER tries
to simply clear the original buffers instead to allow reusing a previous
object.

For instance, the ``hmmsearch`` binary will reallocate a new ``P7_PROFILE`` and
``P7_OPROFILE`` for each HMM in the profile database:

.. code:: C

    // hmmsearch.c

    while (hstatus == eslOK)
      {
          P7_PROFILE      *gm = p7_profile_Create (hmm->M, abc);
          P7_OPROFILE     *om = p7_oprofile_Create(hmm->M, abc)

          /* search all query sequences with the model */
          while (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
          {
            /* ... */
          }

          p7_oprofile_Destroy(om);
          p7_profile_Destroy(gm);

          hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
      }

These ``struct`` are not so large by themselves, but they in turn allocate a
buffer of sufficient size to store the :math:`N` nodes of a HMM.

In PyHMMER, the pipeline will cache memory to be used for the profile and optimized
profiles, and only reallocate if the new HMM to be processed is larger than what the
currently cached ``P7_OPROFILE`` can store.


Cython
------

Cython `freelists <https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#fast-instantiation>`_
help to save allocation cycles for types created in a row (such as `~pyhmmer.plan7.HMM` objects,
which are often yielded one after the other while reading from a `~pyhmmer.plan7.HMMFile`).
