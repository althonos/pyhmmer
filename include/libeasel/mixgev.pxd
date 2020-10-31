cdef extern from "esl_mixgev.h" nogil:

    ctypedef struct ESL_MIXGEV:
        double* q
        double* mu
        double* lambda_
        double* alpha
        double* wrk
        bint* isgumbel
        int K

    ESL_MIXGEV *esl_mixgev_Create(int K)
    void        esl_mixgev_Destroy(ESL_MIXGEV *mg)
    int         esl_mixgev_Copy(ESL_MIXGEV *dest, ESL_MIXGEV *src)
    int         esl_mixgev_ForceGumbel(ESL_MIXGEV *mg, int which)

    double      esl_mixgev_pdf    (double x, ESL_MIXGEV *mg)
    double      esl_mixgev_logpdf (double x, ESL_MIXGEV *mg)
    double      esl_mixgev_cdf    (double x, ESL_MIXGEV *mg)
    double      esl_mixgev_logcdf (double x, ESL_MIXGEV *mg)
    double      esl_mixgev_surv   (double x, ESL_MIXGEV *mg)
    double      esl_mixgev_logsurv(double x, ESL_MIXGEV *mg)
    double      esl_mixgev_invcdf (double p, ESL_MIXGEV *mg)

    double      esl_mixgev_generic_pdf   (double x, void *params)
    double      esl_mixgev_generic_cdf   (double x, void *params)
    double      esl_mixgev_generic_surv  (double x, void *params)
    double      esl_mixgev_generic_invcdf(double p, void *params)

    int         esl_mixgev_Plot(FILE *fp, ESL_MIXGEV *mg,
    				   double (*func)(double x, ESL_MIXGEV *mg),
    				   double xmin, double xmax, double xstep)

    double      esl_mixgev_Sample(ESL_RANDOMNESS *r, ESL_MIXGEV *mg)
    int         esl_mixgev_FitGuess(ESL_RANDOMNESS *r, double *x, int n
    				       ESL_MIXGEV *mg)

    int         esl_mixgev_FitComplete(double *x, int n, ESL_MIXGEV *mg)
