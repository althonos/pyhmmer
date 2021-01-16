from libc.stdio cimport FILE

from libeasel.fileparser cimport ESL_FILEPARSER
from libeasel.random cimport ESL_RANDOMNESS


cdef extern from "esl_mixdchlet.h" nogil:

    ctypedef struct ESL_MIXDCHLET:
        double* q
        double** alpha
        int Q
        int K
        double* postq

    ESL_MIXDCHLET *esl_mixdchlet_Create(int Q, int K)
    void           esl_mixdchlet_Destroy(ESL_MIXDCHLET *dchl)

    double         esl_mixdchlet_logp_c      (ESL_MIXDCHLET *dchl, double *c)
    int            esl_mixdchlet_MPParameters(ESL_MIXDCHLET *dchl, double *c, double *p)

    int            esl_mixdchlet_Fit(double **c, int N, ESL_MIXDCHLET *dchl, double *opt_nll)
    int            esl_mixdchlet_Sample(ESL_RANDOMNESS *rng, ESL_MIXDCHLET *dchl)

    int            esl_mixdchlet_Read(ESL_FILEPARSER *efp, ESL_MIXDCHLET **ret_dchl)
    int            esl_mixdchlet_Write    (FILE *fp, const ESL_MIXDCHLET *dchl)
    int            esl_mixdchlet_WriteJSON(FILE *fp, const ESL_MIXDCHLET *dchl)

    int            esl_mixdchlet_Validate(const ESL_MIXDCHLET *dchl, char *errmsg)
    int            esl_mixdchlet_Compare(const ESL_MIXDCHLET *d1, const ESL_MIXDCHLET *d2, double tol)
    int            esl_mixdchlet_Dump(FILE *fp, const ESL_MIXDCHLET *dchl)
