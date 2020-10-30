cdef extern from "esl_alloc.h" nogil:

    void* esl_alloc_aligned(size_t size, size_t alignment)
    void  esl_alloc_free(void *p)
