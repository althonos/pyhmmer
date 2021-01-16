from libc.stdio cimport FILE

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.fileparser cimport ESL_FILEPARSER
from libeasel.dmatrix cimport ESL_DMATRIX


cdef extern from "esl_scorematrix.h" nogil:

    ctypedef struct ESL_SCOREMATRIX:
        int** s
        int   K
        int   Kp

        char* isval
        const ESL_ALPHABET* abc_r

        int nc
        char* outorder

        char* name
        char* path

    ESL_SCOREMATRIX *esl_scorematrix_Create(const ESL_ALPHABET *abc)
    int              esl_scorematrix_Copy(const ESL_SCOREMATRIX *src, ESL_SCOREMATRIX *dest)
    ESL_SCOREMATRIX *esl_scorematrix_Clone(const ESL_SCOREMATRIX *S)
    int              esl_scorematrix_Compare(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2)
    int              esl_scorematrix_CompareCanon(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2)
    int              esl_scorematrix_Max(const ESL_SCOREMATRIX *S)
    int              esl_scorematrix_Min(const ESL_SCOREMATRIX *S)
    int              esl_scorematrix_IsSymmetric(const ESL_SCOREMATRIX *S)
    int              esl_scorematrix_ExpectedScore(ESL_SCOREMATRIX *S, double *fi, double *fj, double *ret_E)
    int              esl_scorematrix_RelEntropy(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, double lambda_, double *ret_D)
    int              esl_scorematrix_JointToConditionalOnQuery(const ESL_ALPHABET *abc, ESL_DMATRIX *P)
    void             esl_scorematrix_Destroy(ESL_SCOREMATRIX *S)

    int              esl_scorematrix_Set(const char *name, ESL_SCOREMATRIX *S)
    int              esl_scorematrix_SetIdentity(ESL_SCOREMATRIX *S)

    int              esl_scorematrix_SetFromProbs(ESL_SCOREMATRIX *S, double lambda_, const ESL_DMATRIX *P, const double *fi, const double *fj)
    int              esl_scorematrix_SetWAG(ESL_SCOREMATRIX *S, double lambda_, double t)

    int  esl_scorematrix_Read(ESL_FILEPARSER *efp, const ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S)
    int  esl_scorematrix_Write(FILE *fp, const ESL_SCOREMATRIX *S)

    int esl_scorematrix_ProbifyGivenBG(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, double *opt_lambda, ESL_DMATRIX **opt_P)

    int esl_scorematrix_Probify(const ESL_SCOREMATRIX *S, ESL_DMATRIX **opt_P, double **opt_fi, double **opt_fj, double *opt_lambda)
