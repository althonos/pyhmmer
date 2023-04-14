cdef extern from "esl_heap.h" nogil:

    cdef size_t eslHEAP_INITALLOC
    cdef int eslHEAP_MIN
    cdef int eslHEAP_MAX

    ctypedef esl_heap_s ESL_HEAP
    cdef struct esl_heap_s:
        int* idata
        int n
        int nalloc
        int maxormin

    ESL_HEAP *esl_heap_ICreate   (int maxormin);
    int       esl_heap_GetCount  (ESL_HEAP *hp);
    int       esl_heap_IGetTopVal(ESL_HEAP *hp);
    int       esl_heap_Reuse     (ESL_HEAP *hp);
    void      esl_heap_Destroy   (ESL_HEAP *hp);

    int       esl_heap_IInsert(ESL_HEAP *hp, int val);

    int       esl_heap_IExtractTop(ESL_HEAP *hp, int *ret_val);

    int       esl_heap_IGetTop(ESL_HEAP *hp);
