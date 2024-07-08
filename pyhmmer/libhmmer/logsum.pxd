cdef extern from "hmmer.h" nogil:
    int   p7_FLogsumInit()
    float p7_FLogsum(float a, float b)
    float p7_FLogsumError(float a, float b)
    int   p7_ILogsumInit()
    int   p7_ILogsum(int s1, int s2)
