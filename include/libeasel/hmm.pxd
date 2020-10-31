from libc.stdint cimport uint64_t

from libeasel cimport ESL_DSQ
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.random cimport ESL_RANDOMNESS


cdef extern from "esl_hmm.h" nogil:

    ctypedef struct ESL_HMM:
        int M
        int K
        float* pi
        float** t
        float** e
        float** eo
        const ESL_ALPHABET* abc


    ctypedef struct ESL_HMX:
        float** dp
        float* sc
        int M
        int L
        float* dp_mem
        int allocR
        int validR
        int allocM
        uint64_t ncells


    ESL_HMM *esl_hmm_Create(const ESL_ALPHABET *abc, int M);
    ESL_HMM *esl_hmm_Clone(const ESL_HMM *hmm);
    int      esl_hmm_Configure(ESL_HMM *hmm, float *fq);
    int      esl_hmm_SetDegeneracies(ESL_HMM *hmm);
    void     esl_hmm_Destroy(ESL_HMM *hmm);

    ESL_HMX *esl_hmx_Create(int allocL, int allocM);
    int      esl_hmx_GrowTo (ESL_HMX *mx, int L, int M);
    void     esl_hmx_Destroy(ESL_HMX *mx);

    int      esl_hmm_Emit(ESL_RANDOMNESS *r, const ESL_HMM *hmm, ESL_DSQ **opt_dsq, int **opt_path, int *opt_L);
    int      esl_hmm_Forward(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *fwd, float *opt_sc);
    int      esl_hmm_Backward(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *bck, float *opt_sc);
