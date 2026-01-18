from libc.stdio cimport FILE

cdef FILE* fopen_obj(object obj, const char* mode) except NULL
