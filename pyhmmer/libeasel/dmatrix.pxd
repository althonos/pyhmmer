from libc.stdio cimport FILE


cdef extern from "esl_dmatrix.h" nogil:

    cdef enum:
        eslGENERAL
        eslUPPER

    ctypedef struct ESL_DMATRIX:
        double** mx
        int n
        int m
        int type
        int ncells

    ctypedef struct ESL_PERMUTATION:
        int* pi
        int n

    # The ESL_DMATRIX object
    ESL_DMATRIX *esl_dmatrix_Create(int n, int m);
    ESL_DMATRIX *esl_dmatrix_CreateUpper(int n);
    int          esl_dmatrix_Destroy(ESL_DMATRIX *A);
    int          esl_dmatrix_Copy       (const ESL_DMATRIX *src, ESL_DMATRIX *dest);
    ESL_DMATRIX *esl_dmatrix_Clone      (const ESL_DMATRIX *old);
    int          esl_dmatrix_Compare    (const ESL_DMATRIX *A, const ESL_DMATRIX *B, double tol);
    int          esl_dmatrix_CompareAbs (const ESL_DMATRIX *A, const ESL_DMATRIX *B, double tol);
    int          esl_dmatrix_Set        (ESL_DMATRIX *A, double x);
    int          esl_dmatrix_SetZero    (ESL_DMATRIX *A);
    int          esl_dmatrix_SetIdentity(ESL_DMATRIX *A);

    # Debugging/validation for ESL_DMATRIX
    int          esl_dmatrix_Dump(FILE *ofp, const ESL_DMATRIX *A,
    				     const char *rowlabel, const char *collabel);

    # Visualization tools
    int          esl_dmatrix_PlotHeatMap(FILE *fp, ESL_DMATRIX *D, double min, double max);

    # The ESL_PERMUTATION object
    ESL_PERMUTATION *esl_permutation_Create(int n);
    int              esl_permutation_Destroy(ESL_PERMUTATION *P);
    int              esl_permutation_Reuse(ESL_PERMUTATION *P);

    # Debugging/validation for ESL_PERMUTATION
    int esl_permutation_Dump(
        FILE *ofp,
        const ESL_PERMUTATION *P,
    		const char *rowlabel,
        const char *collabel
    );

    # The rest of the dmatrix API
    double       esl_dmx_Max    (const ESL_DMATRIX *A);
    double       esl_dmx_Min    (const ESL_DMATRIX *A);
    double       esl_dmx_Sum    (const ESL_DMATRIX *A);
    int          esl_dmx_MinMax(const ESL_DMATRIX *A, double *ret_min, double *ret_max);
    int          esl_dmx_FrobeniusNorm(const ESL_DMATRIX *A, double *ret_fnorm);
    int          esl_dmx_Multiply(const ESL_DMATRIX *A, const ESL_DMATRIX *B, ESL_DMATRIX *C);
    int          esl_dmx_Exp(const ESL_DMATRIX *Q, double t, ESL_DMATRIX *P);
    int          esl_dmx_Transpose(ESL_DMATRIX *A);
    int          esl_dmx_Add(ESL_DMATRIX *A, const ESL_DMATRIX *B);
    int          esl_dmx_Scale(ESL_DMATRIX *A, double k);
    int          esl_dmx_AddScale(ESL_DMATRIX *A, double k, const ESL_DMATRIX *B);
    int          esl_dmx_Permute_PA(const ESL_PERMUTATION *P, const ESL_DMATRIX *A, ESL_DMATRIX *B);
    int          esl_dmx_LUP_decompose(ESL_DMATRIX *A, ESL_PERMUTATION *P);
    int          esl_dmx_LU_separate(const ESL_DMATRIX *LU, ESL_DMATRIX *L, ESL_DMATRIX *U);
    int          esl_dmx_Invert(const ESL_DMATRIX *A, ESL_DMATRIX *Ai);

    # Optional: interoperability with GSL
    #ifdef HAVE_LIBGSL
    #include <gsl/gsl_matrix.h>
    # int          esl_dmx_MorphGSL(const ESL_DMATRIX *E, gsl_matrix **ret_G);
    # int          esl_dmx_UnmorphGSL(const gsl_matrix *G, ESL_DMATRIX **ret_E);
    #endif

    # Optional: interfaces to LAPACK
    #ifdef HAVE_LIBLAPACK
    # int esl_dmx_Diagonalize(const ESL_DMATRIX *A, double **ret_Er, double **ret_Ei, ESL_DMATRIX **ret_UL, ESL_DMATRIX **ret_UR);
    #endif
