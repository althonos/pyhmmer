cdef extern from "esl_cpu.h" nogil:

    bint   esl_cpu_has_sse(void);
    bint   esl_cpu_has_sse4(void);
    bint   esl_cpu_has_avx(void);
    bint   esl_cpu_has_avx512(void);
    char *esl_cpu_Get(void);
