from libc.stdint cimport uint32_t
from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel cimport ESL_DSQ
from libeasel.sq cimport ESL_SQ
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.random cimport ESL_RANDOMNESS


cdef extern from "hmmer.h" nogil:

    const size_t p7_NEVPARAM
    const size_t p7_NCUTOFFS
    const size_t p7_NOFFSETS
    const size_t p7_MAXABET

    const size_t p7H_NTRANSITIONS
    cdef enum p7h_transitions_e:
        p7H_MM = 0
        p7H_MI = 1
        p7H_MD = 2
        p7H_IM = 3
        p7H_II = 4
        p7H_DM = 5
        p7H_DD = 6

    cdef enum:
        p7H_HASBITS = (1<<0)
        p7H_DESC    = (1<<1)
        p7H_RF      = (1<<2)
        p7H_CS      = (1<<3)
        p7H_XRAY    = (1<<4)
        p7H_HASPROB = (1<<5)
        p7H_HASDNA  = (1<<6)
        p7H_STATS   = (1<<7)
        p7H_MAP     = (1<<8)
        p7H_ACC     = (1<<9)
        p7H_GA      = (1<<10)
        p7H_TC      = (1<<11)
        p7H_NC      = (1<<12)
        p7H_CA      = (1<<13)
        p7H_COMPO   = (1<<14)
        p7H_CHKSUM  = (1<<15)
        p7H_CONS    = (1<<16)
        p7H_MMASK   = (1<<17)

    ctypedef p7_hmm_s P7_HMM
    cdef struct p7_hmm_s:

        int M
        float** t
        float** mat
        float** ins

        char* name
        char* acc
        char* desc
        char* rf
        char* mm
        char* consensus
        char* cs
        char* ca

        char* comlog
        int nseq
        float eff_nseq
        int max_length
        char* ctime
        int* map
        uint32_t checksum
        float[p7_NEVPARAM] evparam
        float[p7_NCUTOFFS] cutoff
        float[p7_MAXABET] compo

        off_t offset
        ESL_ALPHABET* abc
        int flags


    # The P7_HMM object: allocation, initialization, destruction.
    P7_HMM *p7_hmm_Create(int M, const ESL_ALPHABET *abc)
    P7_HMM *p7_hmm_CreateShell()
    int     p7_hmm_CreateBody(P7_HMM *hmm, int M, const ESL_ALPHABET *abc)
    void    p7_hmm_Destroy(P7_HMM *hmm)
    int     p7_hmm_CopyParameters(const P7_HMM *src, P7_HMM *dest)
    P7_HMM *p7_hmm_Clone(const P7_HMM *hmm)
    int     p7_hmm_Zero(P7_HMM *hmm)
    char    p7_hmm_EncodeStatetype(char *typestring)
    char   *p7_hmm_DecodeStatetype(char st)
    # Convenience routines for setting fields in an HMM.
    int     p7_hmm_SetName       (P7_HMM *hmm, char *name);
    int     p7_hmm_SetAccession  (P7_HMM *hmm, char *acc);
    int     p7_hmm_SetDescription(P7_HMM *hmm, char *desc);
    int     p7_hmm_AppendComlog  (P7_HMM *hmm, int argc, char **argv);
    int     p7_hmm_SetCtime      (P7_HMM *hmm);
    int     p7_hmm_SetComposition(P7_HMM *hmm);
    int     p7_hmm_SetConsensus  (P7_HMM *hmm, ESL_SQ *sq);
    # Renormalization and rescaling counts in core HMMs.
    int     p7_hmm_Scale      (P7_HMM *hmm, double scale);
    int     p7_hmm_ScaleExponential(P7_HMM *hmm, double exp);
    int     p7_hmm_Renormalize(P7_HMM *hmm);
    # Debugging and development code.
    int     p7_hmm_Dump(FILE *fp, P7_HMM *hmm);
    int     p7_hmm_Sample          (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
    int     p7_hmm_SampleUngapped  (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
    int     p7_hmm_SampleEnumerable(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
    int     p7_hmm_SampleUniform   (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,
    				     float tmi, float tii, float tmd, float tdd,  P7_HMM **ret_hmm);
    int     p7_hmm_Compare(P7_HMM *h1, P7_HMM *h2, float tol);
    int     p7_hmm_Validate(P7_HMM *hmm, char *errbuf, float tol);
    # Other routines in the API
    int     p7_hmm_CalculateOccupancy(const P7_HMM *hmm, float *mocc, float *iocc);
