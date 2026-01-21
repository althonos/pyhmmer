# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport FILE


cdef class _Win32SynchronizedReader:
    cdef FILE*  file
    cdef object thread

    cpdef object close(self)


cdef class _Win32SynchronizedWriter:
    cdef FILE*  file
    cdef object thread

    cpdef object close(self)


