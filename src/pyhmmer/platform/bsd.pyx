# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport fclose


cdef class _BSDReader:

    def __cinit__(self):
        self.file = NULL

    def __init__(self, object fileobj):
        self.file = fileobj_bsd_open(fileobj, "r")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    cpdef object close(self):
        if self.file != NULL:
            fclose(self.file)


cdef class _BSDWriter:

    def __cinit__(self):
        self.file = NULL

    def __init__(self, object fileobj):
        self.file = fileobj_bsd_open(fileobj, "w")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    cpdef object close(self):
        if self.file != NULL:
            fclose(self.file)

