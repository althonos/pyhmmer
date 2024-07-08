cdef extern from "esl_arr2.h" nogil:

    size_t esl_arr2_SSizeof(char **s, int dim1)
    void   esl_arr2_Destroy(void **p, int dim1)
