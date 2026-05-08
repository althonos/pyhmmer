from libeasel.sq cimport ESL_SQ
from libeasel.random cimport ESL_RANDOMNESS
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_trace cimport P7_TRACE


cdef extern from "hmmer.h" nogil:
    int p7_CoreEmit   (ESL_RANDOMNESS *r, const P7_HMM *hmm,                                        ESL_SQ *sq, P7_TRACE *tr)
    int p7_ProfileEmit(ESL_RANDOMNESS *r, const P7_HMM *hmm, const P7_PROFILE *gm, const P7_BG *bg, ESL_SQ *sq, P7_TRACE *tr)
    int p7_emit_SimpleConsensus(const P7_HMM *hmm, ESL_SQ *sq)
    int p7_emit_FancyConsensus (const P7_HMM *hmm, float min_lower, float min_upper, ESL_SQ *sq)
