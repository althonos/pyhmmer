cdef extern from "esl_arr3.h" nogil:

    size_t esl_arr3_SSizeof(char ***s, int dim1, int dim2)
    void   esl_arr3_Destroy(void ***p, int dim1, int dim2)
