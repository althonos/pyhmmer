from libc.stdio cimport FILE

cdef extern from "fileobj/linux.h":
    FILE* fileobj_linux_open(object obj, const char* mode) except NULL