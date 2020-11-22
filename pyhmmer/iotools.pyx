from cpython.ref cimport PyObject
from cpython.object cimport PyObject_CallMethodObjArgs
from libc.stdio cimport FILE


cdef size_t fwrite_obj(void *cookie, const char *buf, size_t size) except 0:
    cdef PyObject* obj = <PyObject*> cookie
    return (<object> cookie).write(buf[:size])


cdef FILE* fopen_obj(PyObject* obj, str mode = "r"):
    cdef void* cookie = <void*> obj

    cdef cookie_io_functions_t functions
    functions.read = NULL
    functions.write = <cookie_write_function_t*> fwrite_obj
    functions.seek = NULL
    functions.close = NULL

    return fopencookie(cookie, mode.encode("ascii"), functions)
