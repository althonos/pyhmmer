cdef extern from "esl_huffman.h" nogil:

    const size_t eslHUFFMAN_MAXCODE

    ctypedef huffman_s ESL_HUFFMAN
    cdef struct huffman_s:
        int* len
        uint32_t* code
        int K
        int* sorted_at
        int Ku
        int* dt_len
        uint32_t* dt_lcode
        int* dt_rank
        int D
        int Lmax


    int  esl_huffman_Build(const float *fq, int K, ESL_HUFFMAN **ret_hc)
    void esl_huffman_Destroy(ESL_HUFFMAN *hc)

    int  esl_huffman_Encode(const ESL_HUFFMAN *hc, const char     *T, int n,  uint32_t **ret_X, int *ret_nb)
    int  esl_huffman_Decode(const ESL_HUFFMAN *hc, const uint32_t *X, int nb, char     **ret_T, int *ret_n)

    int  esl_huffman_Dump(FILE *fp, ESL_HUFFMAN *hc)
