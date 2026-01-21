# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport FILE


cdef extern from "platform/linux.h":
    FILE* fileobj_linux_open(object obj, const char* mode) except NULL


cdef class _LinuxReader:
    cdef object obj
    cdef FILE*  file

    cpdef object close(self)


cdef class _LinuxWriter:
    cdef object obj
    cdef FILE*  file

    cpdef object close(self)

