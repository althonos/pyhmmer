Benchmarks
==========

For all the following figures, the X axis (*CPUS*) indicates the number of
background worker threads being spawned, **in addition to the main thread**.

Benchmarks were run on a `i7-10710U CPU <https://ark.intel.com/content/www/us/en/ark/products/196448/intel-core-i7-10710u-processor-12m-cache-up-to-4-70-ghz.html>`_
running @1.10GHz with 6 physical / 12 logical cores, using a FASTA file
containing 2,100 protein sequences extracted from the genome of *Anaerococcus provencensis*
and the version 33.1 of the `Pfam <https://pfam.xfam.org/>`_ HMM library containing
18,259 domains. Commands were run 4 times on a warm SSD. *Plain lines show
the times for pressed HMMs, and dashed-lines the times for HMMs in text format.*


Unreleased
----------


v0.4.0 - 2021-06-05
-------------------

.. image:: _images/bench-v0.4.0.svg

The overhead of pyHMMER has been reduced, and has a much smaller effect when
using a high number of threads.

The main thread has been updated so that it only loads the next HMM from the
queries when a worker needs a new one, instead of saturating the input queue
with all queries like previously. This should avoid running the main thread
and all worker threads at the same time, and help the CPU to run the SIMD code.
The fastest run is now with 4+1 threads (45.8s), which is still smaller than the
total number of physical CPUs (6 cores).

The main thread also now yields result instead of waiting for all HMMs to be
done, thus potentially reducing the size of the result queue, which may in turn
speed up insertions by worker threads.



v0.3.0 - 2021-03-11
-------------------

.. image:: _images/bench-v0.3.0.svg

The small number of proteins renders the HMMER parallelisation useless for
any number of worker threads higher than 2 (because of the block size being
set to 1,000 sequences).

The fastest run is obtained with 3+1 threads (53.8s), which is lower than the
total number of physical CPUs (6 cores). This could be a hint of hindrance
between the different threads.

Loading from a pressed HMM saves a constant time, independently of the number
of threads. pyHMMER also has a constant overhead compared to HMMER for a
higher number of threads.
