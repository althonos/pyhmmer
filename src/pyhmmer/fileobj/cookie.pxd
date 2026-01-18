# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport FILE


cdef extern from "cookie.h":
    FILE* fileobj_linux_open(object obj, const char* mode) except NULL


cdef class _CookieReader:
    cdef FILE*  file

    cpdef object close(self)


cdef class _CookieWriter:
    cdef FILE*  file

    cpdef object close(self)

