from libc.stdint cimport uint64_t, int64_t
from libc.stdio cimport FILE


cdef extern from "esl_rand64.h" nogil:

    ctypedef struct ESL_RAND64:
        int mti
        uint64_t[312] mt
        uint64_t seed


    ESL_RAND64 *esl_rand64_Create (uint64_t seed)
    int         esl_rand64_Init   (ESL_RAND64 *rng, uint64_t seed)
    uint64_t    esl_rand64_GetSeed(ESL_RAND64 *rng)
    void        esl_rand64_Destroy(ESL_RAND64 *rng)

    uint64_t    esl_rand64(ESL_RAND64 *rng)
    uint64_t    esl_rand64_Roll(ESL_RAND64 *rng, uint64_t n);
    double      esl_rand64_double(ESL_RAND64 *rng);
    double      esl_rand64_double_closed(ESL_RAND64 *rng);
    double      esl_rand64_double_open(ESL_RAND64 *rng);

    int         esl_rand64_Deal(ESL_RAND64 *rng, int64_t m, int64_t n, int64_t *deal)

    int         esl_rand64_Dump(FILE *fp, ESL_RAND64 *rng)
