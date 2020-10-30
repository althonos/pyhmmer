cdef extern from "esl_cluster.h" nogil:

    int esl_cluster_SingleLinkage(
        void *base, size_t n, size_t size,
        int (*linkfunc)(const void *, const void *, const void *, int *),
        void *param,
        int *workspace, int *assignments, int *ret_C
    );
