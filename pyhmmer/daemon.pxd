# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport uint16_t, uint64_t
from libhmmer.p7_pipeline cimport p7_pipemodes_e

from pyhmmer.plan7 cimport TopHits, Pipeline


# --- Cython classes ---------------------------------------------------------

cdef class Client:

    cdef readonly str            address
    cdef readonly uint16_t       port
    cdef readonly object         socket
    cdef          p7_pipemodes_e mode

    cdef bytearray _recvall(self, size_t message_size)
    cdef TopHits _client(
        self,
        bytes query,
        uint64_t db,
        list ranges,
        Pipeline pli,
        p7_pipemodes_e mode,
    )
