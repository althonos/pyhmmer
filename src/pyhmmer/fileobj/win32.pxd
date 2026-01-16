from libc.stdio cimport FILE

import io

cdef inline FILE* fopen_obj(object obj, const char* mode) except NULL:
    raise io.UnsupportedOperation("Cannot open a file-like object on Windows")