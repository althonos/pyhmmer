cdef extern from "esl_cpu.h" nogil:

    int   esl_cpu_has_sse(void);
    int   esl_cpu_has_sse4(void);
    int   esl_cpu_has_avx(void);
    int   esl_cpu_has_avx512(void);
    char *esl_cpu_Get(void);
