# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport fclose


cdef class _LinuxReader:

    def __cinit__(self):
        self.file = NULL

    def __init__(self, object fileobj):
        self.file = fileobj_linux_open(fileobj, "r")
        self.obj = fileobj

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    cpdef object close(self):
        if self.file != NULL:
            fclose(self.file)


cdef class _LinuxWriter:

    def __cinit__(self):
        self.file = NULL

    def __init__(self, object fileobj):
        self.file = fileobj_linux_open(fileobj, "w")
        self.obj = fileobj

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    cpdef object close(self):
        if self.file != NULL:
            fclose(self.file)


