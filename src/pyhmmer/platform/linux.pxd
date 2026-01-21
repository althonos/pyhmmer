# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport FILE


cdef extern from "fileobj/linux.h":
    FILE* fileobj_linux_open(object obj, const char* mode) except NULL


cdef class _LinuxReader:
    cdef FILE*  file

    cpdef object close(self)


cdef class _LinuxWriter:
    cdef FILE*  file

    cpdef object close(self)

