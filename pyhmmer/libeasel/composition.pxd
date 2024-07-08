cdef extern from "esl_composition.h" nogil:

    int esl_composition_BL62(double *f);
    int esl_composition_WAG (double *f);
    int esl_composition_SW34(double *f);
    int esl_composition_SW50(double *f);
