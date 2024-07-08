from libc.stdio cimport FILE
from libc.stdint cimport uint32_t


cdef extern from "esl_random.h" nogil:

    cdef enum esl_randomness_type:
        eslRND_FAST = 0
        eslRND_MERSENNE = 1

    ctypedef struct ESL_RANDOMNESS:
        esl_randomness_type      type
        int      mti
        uint32_t[624] mt
        uint32_t x
        uint32_t seed

    # The ESL_RANDOMNESS object.
    ESL_RANDOMNESS *esl_randomness_Create    (uint32_t seed);
    ESL_RANDOMNESS *esl_randomness_CreateFast(uint32_t seed);   # DEPRECATED. Use esl_randomness_Create.  The Knuth LCG used to have a speed advantage for us, but MT is fast.
    ESL_RANDOMNESS *esl_randomness_CreateTimeseeded();          # DEPRECATED. Use esl_randomness_Create(0)
    void            esl_randomness_Destroy(ESL_RANDOMNESS *r);
    int             esl_randomness_Init(ESL_RANDOMNESS *r, uint32_t seed);
    uint32_t        esl_randomness_GetSeed(const ESL_RANDOMNESS *r);

    # The generator, esl_random().
    double   esl_random       (ESL_RANDOMNESS *r);
    uint32_t esl_random_uint32(ESL_RANDOMNESS *r);

    uint32_t esl_rnd_mix3(uint32_t a, uint32_t b, uint32_t c);

    # Debugging/development tools.
    int esl_randomness_Dump(FILE *fp, ESL_RANDOMNESS *r);

    # Other fundamental sampling (including Gaussian, gamma).
    double esl_rnd_UniformPositive(ESL_RANDOMNESS *r);
    double esl_rnd_Gaussian (ESL_RANDOMNESS *rng, double mean, double stddev);
    double esl_rnd_Gamma    (ESL_RANDOMNESS *rng, double a);
    int    esl_rnd_Dirichlet(ESL_RANDOMNESS *rng, const double *alpha, int K, double *p);  # Pass alpha=NULL if you just want a uniform draw.
    int    esl_rnd_Deal     (ESL_RANDOMNESS *rng, int m, int n, int *deal);

    # Multinomial sampling from discrete probability n-vectors.
    int    esl_rnd_DChoose   (ESL_RANDOMNESS *r, const double *p,   int N);
    int    esl_rnd_FChoose   (ESL_RANDOMNESS *r, const float  *p,   int N);
    int    esl_rnd_DChooseCDF(ESL_RANDOMNESS *r, const double *cdf, int N);
    int    esl_rnd_FChooseCDF(ESL_RANDOMNESS *r, const float  *cdf, int N);

    # Random data generators (unit testing, etc.)
    int    esl_rnd_mem        (ESL_RANDOMNESS *rng, void *buf, int n);
    int    esl_rnd_floatstring(ESL_RANDOMNESS *rng, char *s);
