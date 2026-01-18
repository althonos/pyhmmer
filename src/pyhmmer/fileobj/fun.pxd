# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport FILE


cdef extern from "fun.h":
    FILE* fileobj_bsd_open(object obj, const char* mode) except NULL


cdef class _FunReader:
    cdef FILE*  file

    cpdef object close(self)


cdef class _FunWriter:
    cdef FILE*  file

    cpdef object close(self)

