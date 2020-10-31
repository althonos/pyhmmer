from libc.stdio cimport FILE


cdef extern from "esl_mininizer.h" nogil:

    cdef int eslMIN_MAXITER
    cdef float eslMIN_CG_RTOL
    cdef float eslMIN_CG_ATOL
    cdef int eslMIN_BRACK_MAXITER
    cdef float eslMIN_BRACK_STEP
    cdef float eslMIN_BRENT_RTOL
    cdef float eslMIN_BRENT_ATOL
    cdef float eslMIN_DERIV_STEP


    ctypedef struct ESL_MIN_CFG:
        int max_iterations
        double cg_rtol
        double cg_atol
        double brent_rtol
        double brent_atol
        int brack_maxiter
        double derive_step
        double* u
        int n


    ctypedef struct ESL_MIN_DAT:
        int niter
        double* fx
        int* brack_n
        double* brack_ax
        double* brack_bx
        double* brack_cx
        double* brack_fa
        double* brack_fb
        double* brack_fc
        int* brent_n
        double* brent_x
        int* nfunc


    int esl_min_ConjugateGradientDescent(ESL_MIN_CFG *cfg, double *x, int n,
					    double (*func)(double *, int, void *),
					    void (*dfunc)(double *, int, void *, double *),
					    void *prm, double *ret_fx, ESL_MIN_DAT *dat)

    ESL_MIN_CFG *esl_min_cfg_Create(int n)
    void         esl_min_cfg_Destroy(ESL_MIN_CFG *cfg)

    ESL_MIN_DAT *esl_min_dat_Create(ESL_MIN_CFG *cfg)
    void         esl_min_dat_Destroy(ESL_MIN_DAT *dat)
    int          esl_min_dat_Dump(FILE *fp, ESL_MIN_DAT *dat)
