# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport FILE


cdef FILE* fopen_obj(object obj, const char* mode) except NULL


cdef class _WinSynchronizedReader:
    cdef FILE*  file
    cdef object thread

    cpdef object close(self)


cdef class _WinSynchronizedWriter:
    cdef FILE*  file
    cdef object thread

    cpdef object close(self)


