Benchmarks
==========

For all the following figures, the X axis (*CPUS*) indicates the number of
background worker threads being spawned, **in addition to the main thread**,
which performs the I/O and dispatches the jobs.

Benchmarks were run on a `i7-10710U CPU <https://ark.intel.com/content/www/us/en/ark/products/196448/intel-core-i7-10710u-processor-12m-cache-up-to-4-70-ghz.html>`_
running @1.10GHz with 6 physical / 12 logical cores, using a FASTA file
containing 2,100 protein sequences extracted from the genome of *Anaerococcus provencensis*
and the version 33.1 of the `Pfam <https://pfam.xfam.org/>`_ HMM library containing
18,259 domains. Commands were run 4 times on a warm SSD. *Plain lines show
the times for pressed HMMs, and dashed-lines the times for HMMs in text format.*


v0.7.0 - 2022-12-04
-------------------

.. image:: /_images/bench-v0.7.0.svg

Collections for storing HMMs and Sequences were updated to allow iterating 
without the GIL. For `hmmscan`, the `OptimizedProfileBlock` store an array 
of semaphores to avoid concurrent reconfiguration by multiple pipeline 
across different threads.


v0.5.0 - 2022-03-14
-------------------

.. image:: /_images/bench-v0.5.0.svg

A new dedicated collection has been added to store the target sequences of a
`~pyhmmer.plan7.Pipeline`, saving some overhead if the same target sequences
are reused with multiple queries.


v0.4.5 - 2021-07-19
-------------------

.. image:: /_images/bench-v0.4.5.svg

By adding an extra requirement on the reference sequences passed to a `~pyhmmer.plan7.Pipeline`,
the Cython code can now evaluate a single HMM against the entirety of the reference
sequences without having to reacquire the GIL between each reference sequence. This now
means all the sequences must be stored in memory and cannot be obtained from an
iterator anymore. However, since this was already required by `pyhmmer.hmmer.hmmsearch`,
the changes in the downstream code should be minimal.

Another consequence of this change is that threads can operate independently for
much longer periods of time. This will greatly improve multi-threaded performance
on workstations with many CPUs. However, the speedup becomes negligible once the
number of threads grows above the number of physical cores because of the SIMD
code making heavy use of CPU registers.

The fastest run with 11+1 threads (40.2s) is barely faster than the run using an
optimal number of 5+1 threads (42.2s).


v0.4.0 - 2021-06-05
-------------------

.. image:: /_images/bench-v0.4.0.svg

The overhead of PyHMMER has been reduced, and has a much smaller effect when
using a high number of threads.

The main thread has been updated so that it only loads the next `pyhmmer.plan7.HMM`
from the queries when a worker needs a new one, instead of saturating the input queue
with all queries like previously. This should avoid running the main thread
and all worker threads at the same time, and help the CPU to run the SIMD code.
The fastest run is now with 4+1 threads (45.8s), which is still smaller than the
total number of physical CPUs (6 cores) for some reason.

The main thread also now yields result instead of waiting for all HMMs to be
done, thus potentially reducing the size of the result queue, which may in turn
speed up insertions by worker threads.


v0.3.0 - 2021-03-11
-------------------

.. image:: /_images/bench-v0.3.0.svg

The small number of proteins renders the HMMER parallelisation useless for
any number of worker threads higher than 2 (because of the block size being
set to 1,000 sequences).

The fastest run is obtained with 3+1 threads (53.8s), which is lower than the
total number of physical CPUs (6 cores). This could be a hint of hindrance
between the different threads.

Loading from a pressed HMM saves a constant time, independently of the number
of threads. PyHMMER also has a constant overhead compared to HMMER for a
higher number of threads.
