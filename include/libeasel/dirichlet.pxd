cdef extern from "esl_dirichlet.h" nogil:

    double esl_dirichlet_logpdf  (double *p, double *alpha, int K)
    double esl_dirichlet_logpdf_c(double *c, double *alpha, int K)

    # Sampling
    int esl_dirichlet_DSample       (ESL_RANDOMNESS *r, double *alpha, int K, double *p)
    int esl_dirichlet_FSample       (ESL_RANDOMNESS *r, float  *alpha, int K, float  *p)
    int esl_dirichlet_DSampleUniform(ESL_RANDOMNESS *r, int K, double *p)
    int esl_dirichlet_FSampleUniform(ESL_RANDOMNESS *r, int K, float  *p)
    int esl_dirichlet_SampleBeta    (ESL_RANDOMNESS *r, double theta1, double theta2, double *ret_answer)
