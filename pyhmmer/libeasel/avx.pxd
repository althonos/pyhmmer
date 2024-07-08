cdef extern from "x86intrin.h":

    ctypedef struct __m256i:
        pass

    ctypedef struct __m256:
        pass


cdef extern from "esl_avx.h" nogil:

    void esl_avx_dump_256i_hex4(__m256i v);

    uint8_t esl_avx_hmax_epu8(__m256i a)
    int8_t esl_avx_hmax_epi8(__m256i a)
    int16_t esl_avx_hmax_epi16(__m256i a)
    void esl_avx_hsum_ps(__m256 a, float *ret_sum)
    __m256i esl_avx_rightshift_int8(__m256i v, __m256i neginfmask)
    __m256i esl_avx_rightshift_int16(__m256i v, __m256i neginfmask)
    __m256 esl_avx_rightshiftz_float(__m256 v)
    __m256 esl_avx_leftshiftz_float(__m256 v)
    int esl_avx_any_gt_epi16(__m256i a, __m256i b)
