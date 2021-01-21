from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.mixdchlet cimport ESL_MIXDCHLET
from libhmmer.p7_hmm cimport P7_HMM

cdef extern from "hmmer.h" nogil:

    ctypedef p7_prior_s P7_PRIOR
    cdef struct p7_prior_s:
        ESL_MIXDCHLET* tm
        ESL_MIXDCHLET* ti
        ESL_MIXDCHLET* td
        ESL_MIXDCHLET* em
        ESL_MIXDCHLET* ei

    P7_PRIOR  *p7_prior_CreateAmino()
    P7_PRIOR  *p7_prior_CreateNucleic()
    P7_PRIOR  *p7_prior_CreateLaplace(const ESL_ALPHABET *abc)
    void       p7_prior_Destroy(P7_PRIOR *pri)

    int        p7_ParameterEstimation(P7_HMM *hmm, const P7_PRIOR* pri)
