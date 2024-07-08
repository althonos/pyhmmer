from .p7_hmm cimport P7_HMM
from .p7_profile cimport P7_PROFILE
from .p7_bg cimport P7_BG


cdef extern from "hmmer.h" nogil:
    int p7_ProfileConfig(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int L, int mode);
    int p7_ReconfigLength  (P7_PROFILE *gm, int L);
    int p7_ReconfigMultihit(P7_PROFILE *gm, int L);
    int p7_ReconfigUnihit  (P7_PROFILE *gm, int L);
