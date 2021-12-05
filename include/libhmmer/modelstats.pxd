from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_bg cimport P7_BG


cdef extern from "hmmer.h" nogil:

    cdef double p7_MeanMatchInfo(const P7_HMM* hmm, const P7_BG* bg);
    cdef double p7_MeanMatchEntropy(const P7_HMM *hmm)
    cdef double p7_MeanMatchRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg)
    cdef double p7_MeanForwardScore(const P7_HMM *hmm, const P7_BG *bg)
    cdef int p7_MeanPositionRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg, double *ret_entropy)
    cdef int p7_hmm_CompositionKLD(const P7_HMM *hmm, const P7_BG *bg, float *ret_KL, float **opt_avp)
