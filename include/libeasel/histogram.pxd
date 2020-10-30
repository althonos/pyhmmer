cdef extern from "esl_histogram.h" nogil:

    cdef enum esl_dataset:
        COMPLETE
        VIRTUAL_CENSORED
        TRUE_CENSORED

    ctypedef struct ESL_HISTOGRAM:
        uint64_t* obs
        int nb
        double w
        double bmin
        double bmax
        int imin
        int imax

        double xmin
        double xmax
        uint64_t n
        double* x
        uint64_t nalloc

        double phi
        int cmin
        uint64_t z
        uint64_t Nc
        uint64_t No

        double* expect
        int emin
        double tailbase
        double tailmass

        bint is_full
        bint is_done
        bint is_sorted
        bint is_tailfit
        bint is_rounded

        esl_dataset dataset_is


    # Creating/destroying histograms and collecting data
    ESL_HISTOGRAM *esl_histogram_Create    (double bmin, double bmax, double w);
    ESL_HISTOGRAM *esl_histogram_CreateFull(double bmin, double bmax, double w);
    void           esl_histogram_Destroy  (ESL_HISTOGRAM *h);
    int            esl_histogram_Score2Bin(ESL_HISTOGRAM *h, double x, int *ret_b);
    int            esl_histogram_Add      (ESL_HISTOGRAM *h, double x);

    # Declarations about the binned data before parameter fitting
    int esl_histogram_DeclareCensoring(ESL_HISTOGRAM *h, int z, double phi);
    int esl_histogram_DeclareRounding (ESL_HISTOGRAM *h);
    int esl_histogram_SetTail         (ESL_HISTOGRAM *h, double phi,
    					  double *ret_newmass);
    int esl_histogram_SetTailByMass   (ESL_HISTOGRAM *h, double pmass,
    					  double *ret_newmass);

    # Accessing data samples in a full histogram:
    int esl_histogram_GetRank(ESL_HISTOGRAM *h, int rank, double *ret_x);
    int esl_histogram_GetData(ESL_HISTOGRAM *h, double **ret_x, int *ret_n);
    int esl_histogram_GetTail(ESL_HISTOGRAM *h, double phi, double **ret_x,
    				 int *ret_n, int *ret_z);
    int esl_histogram_GetTailByMass(ESL_HISTOGRAM *h, double pmass,
    				       double **ret_x, int *ret_n, int *ret_z);


    # Setting expected binned counts:
    int esl_histogram_SetExpect(ESL_HISTOGRAM *h,
    				   double (*cdf)(double x, void *params),
    				   void *params);
    int esl_histogram_SetExpectedTail(ESL_HISTOGRAM *h, double base_val,
    					 double pmass,
    					 double (*cdf)(double x, void *params),
    					 void *params);

    # Output/display of binned data:
    int esl_histogram_Write       (FILE *fp, ESL_HISTOGRAM *h);
    int esl_histogram_Plot        (FILE *fp, ESL_HISTOGRAM *h);
    int esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h);
    int esl_histogram_PlotQQ      (FILE *fp, ESL_HISTOGRAM *h,
    				      double (*invcdf)(double, void *), void *params);

    # Goodness of fit testing
    int esl_histogram_Goodness(ESL_HISTOGRAM *h, int nfitted,
    				  int *ret_nbins,
    				  double *ret_G,  double *ret_Gp,
    				  double *ret_X2, double *ret_X2p);
