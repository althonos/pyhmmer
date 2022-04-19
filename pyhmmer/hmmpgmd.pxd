# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libhmmer.p7_pipeline cimport p7_pipemodes_e


# --- Cython classes ---------------------------------------------------------

cdef class Client:

    cdef str            address
    cdef uint16_t       port
    cdef object         socket
    cdef p7_pipemodes_e mode
