from libhmmer.p7_tophits cimport P7_TOPHITS

cdef extern from "reexports/p7_tophits.h" nogil:
    int p7_tophits_Reuse(P7_TOPHITS* h)
