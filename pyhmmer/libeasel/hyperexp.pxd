cdef extern from "esl_hyperexp.h" nogil:

    ctypedef struct ESL_HYPEREXP:
        double* q
        double* lambda_
        double* wrk
        double  mu
        int     K
        char*   fixlambda
        int     fixmix

    ESL_HYPEREXP *esl_hyperexp_Create(int K);
    void          esl_hyperexp_Destroy(ESL_HYPEREXP *h);
    int           esl_hyperexp_Copy(ESL_HYPEREXP *src, ESL_HYPEREXP *dest);
    int           esl_hyperexp_FixedUniformMixture(ESL_HYPEREXP *h);
    int           esl_hyperexp_SortComponents(ESL_HYPEREXP *h);
    int           esl_hyperexp_Write(FILE *fp, ESL_HYPEREXP *hxp);
    int           esl_hyperexp_Dump(FILE *fp, ESL_HYPEREXP *hxp);

    int           esl_hyperexp_Read(ESL_FILEPARSER *ef, ESL_HYPEREXP **ret_hxp);
    int           esl_hyperexp_ReadFile(char *filename, ESL_HYPEREXP **ret_hxp);


    double  esl_hxp_pdf    (double x, ESL_HYPEREXP *h);
    double  esl_hxp_logpdf (double x, ESL_HYPEREXP *h);
    double  esl_hxp_cdf    (double x, ESL_HYPEREXP *h);
    double  esl_hxp_logcdf (double x, ESL_HYPEREXP *h);
    double  esl_hxp_surv   (double x, ESL_HYPEREXP *h);
    double  esl_hxp_logsurv(double x, ESL_HYPEREXP *h);
    double  esl_hxp_invcdf (double p, ESL_HYPEREXP *h);

    double  esl_hxp_generic_pdf   (double x, void *params);
    double  esl_hxp_generic_cdf   (double x, void *params);
    double  esl_hxp_generic_surv  (double x, void *params);
    double  esl_hxp_generic_invcdf(double x, void *params);

    int esl_hxp_Plot(FILE *fp, ESL_HYPEREXP *h,
    			double (*func)(double x, ESL_HYPEREXP *h),
    			double xmin, double xmax, double xstep);


    double esl_hxp_Sample(ESL_RANDOMNESS *r, ESL_HYPEREXP *h);

    int esl_hxp_FitGuess   (double *x, int n, ESL_HYPEREXP *h);
    int esl_hxp_FitComplete(double *x, int n, ESL_HYPEREXP *h);

    int esl_hxp_FitGuessBinned   (ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
    int esl_hxp_FitCompleteBinned(ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
