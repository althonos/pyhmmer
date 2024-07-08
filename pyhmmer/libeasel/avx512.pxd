cdef extern from "x86intrin.h":

    ctypedef struct __m256i:
        pass

    ctypedef struct __m256:
        pass


cdef extern from "esl_avx512.h" nogil:

    void esl_avx512_dump_512i_hex8(__m512i v)
    uint8_t esl_avx512_hmax_epu8(__m512i a)
    int8_t esl_avx512_hmax_epi8(__m512i a)
    int16_t esl_avx512_hmax_epi16(__m512i a)
    void esl_avx512_hsum_ps(__m512 a, float *ret_sum)
    __m512i esl_avx512_rightshift_int8(__m512i v, __m512i neginfmask)
    __m512i esl_avx512_rightshift_int16(__m512i v, __m512i neginfmask)
    __m512 esl_avx512_rightshiftz_float(__m512 v)
    __m512 esl_avx512_leftshiftz_float(__m512 v)
