cdef extern from "esl_graph.h" nogil:

    int esl_graph_MaxBipartiteMatch(int **A, int M, int N, int ***opt_G, int *ret_nedges)
