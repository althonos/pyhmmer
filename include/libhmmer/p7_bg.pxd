from libc.stdio cimport FILE

from libeasel cimport ESL_DSQ
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.hmm cimport ESL_HMM


cdef extern from "hmmer.h" nogil:

    ctypedef p7_bg_s P7_BG
    cdef struct p7_bg_s:
        float* f
        float p1
        ESL_HMM* fhmm
        float omega
        const ESL_ALPHABET* abc

    P7_BG *p7_bg_Create(const ESL_ALPHABET *abc)
    P7_BG *p7_bg_CreateUniform(const ESL_ALPHABET *abc)
    P7_BG *p7_bg_Clone(const P7_BG *bg)
    int    p7_bg_Dump(FILE *ofp, const P7_BG *bg)
    void   p7_bg_Destroy(P7_BG *bg)
    int    p7_bg_SetLength(P7_BG *bg, int L)
    int    p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc)

    int    p7_bg_Read(char *bgfile, P7_BG *bg, char *errbuf)
    int    p7_bg_Write(FILE *fp, P7_BG *bg)

    int    p7_bg_SetFilter  (P7_BG *bg, int M, const float *compo)
    int    p7_bg_FilterScore(P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc)
