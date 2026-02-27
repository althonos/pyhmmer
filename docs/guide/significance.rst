Bitscore, P-values & E-values
=============================

.. currentmodule:: pyhmmer.plan7

When running ``hmmsearch`` or ``hmmscan``, the hits are reported with 
different statistics: bitscores, P-values, and E-values. What exactly 
do these mean, and how should you filter hits based on these values?

.. note:: 

    This page is partially documented from the `HMMER User's Guide <http://eddylab.org/software/hmmer/Userguide.pdf>`_.
    Have a look there too for more detailed explanations, in particular 
    the section titled *The HMMER profile/sequence comparison pipeline*
    which is hastily summarized.


Hit bit-score
-------------

The bit-score (`Hit.score`) is the score given by the alignment of a 
sequence to a HMM. It is expressed in `bits <https://en.wikipedia.org/wiki/Information_theory#Units>`_.
Once a sequence has passed all pipeline filters, the full forward pass is 
computed to obtain the full score. Afterwards, the "domain definition workflow"
is run to identify domains given by multiple HMM matches to the same sequence.
An ad hoc biased composition score correction is calculated for each envelope,
and the corrected bit-score is calculated for each domain and at the sequence 
level.

Hit P-value
-----------

When a HMM is constructed, it runs through a calibration stage which computes
a histogram of hit scores obtained by the HMM on random sequences. The 
latent distribution is modeled after an `Exponential distribution <https://en.wikipedia.org/wiki/Exponential_distribution>`_
and the values obtained with bootstrap are used to determine the two 
parameters of the law, :math:`\tau` and :math:`\lambda`. These parameters
can be accessed through the `HMM.evalue_parameters` property.
Using these parameters, a P-value can be computed from any score given by 
the HMM using the `survival function <https://en.wikipedia.org/wiki/Survival_function>`_
defined by the parameters, i.e. corresponds to :math:`P(X \ge score)`.


Hit E-value
-----------

As the P-value is computed independently for every hit, the 
`multiple testing problem <multiple testing problem>`_ occurs. To address the 
false discovery rate, HMMER performs multiple-testing correction using 
`Bonferroni's method <https://en.wikipedia.org/wiki/Bonferroni_correction>`_.
The number of testing hypotheses, :math:`Z`, can set from the command line or 
from the API, but is usually set automatically by HMMER as the number of targets
(target sequences for ``hmmsearch``, target HMMs for ``hmmscan``). Then, the 
`Hit.evalue` is computed by multiplying `Hit.pvalue` with :math:`Z`.


Domain bit-score
----------------

Once a sequence has passed all pipeline filters, the full forward pass is 
computed to obtain the hit score. Afterwards, the "domain definition workflow"
is run to identify domains given by multiple HMM matches to the same sequence.
The bit score of each domain is then computed based on the local alignment 
of this domain only.


Domain P-value
--------------

The domain P-value is computed identically to the hit P-value, but using the
domain bit-score instead of the hit bit-score (using the same E-value parameters
from the HMM).


Domain E-value(s)
-----------------

For domains, we still have a multiple-testing problem, however the the number
of testing hypotheses is slightly different, because only hits passing the 
pipeline filters (defined with the :math:`F_1`, :math:`F_2` and :math:`F_3`
parameters of a `Pipeline`) actually get through the domain definition workflow
and are tested for significance. Therefore, HMMER reports two ways to correct 
the domain P-values: a *conditional* E-value which takes into account only the 
number of domains passing the first pipeline stages (counted as :math:`domZ`), 
and an *independent* E-value which takes into account all sequence/HMM comparison
independently (counted as :math:`Z`). While `Domain.c_evalue` is more statistically
accurate (and usually reports E-values lower than `Domain.i_evalue`), it is 
harder to replicate the results without setting :math:`domZ` a-priori, which is 
usually not trivial to do accurately.
