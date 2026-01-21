# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport FILE


cdef extern from "platform/bsd.h":
    FILE* fileobj_bsd_open(object obj, const char* mode) except NULL


cdef class _BSDReader:
    cdef FILE*  file

    cpdef object close(self)


cdef class _BSDWriter:
    cdef FILE*  file

    cpdef object close(self)

