from libc.stdio cimport FILE

cdef extern from "fileobj/bsd.h":
    FILE* fileobj_bsd_open(object obj, const char* mode) except NULL