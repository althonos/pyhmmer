# coding: utf-8
"""Obtain a `FILE*` from a Python object using ``funopen``.
"""

cdef extern from "stdio.h":

    FILE* funopen(
        const void  *cookie,
        int        (*readfn)   (void *,       char *, int),
        int        (*writefn)  (void *, const char *, int),
        fpos_t     (*seekfn)   (void *,       fpos_t, int),
        int        (*closefn)  (void *)
    );

    FILE* fropen(void *cookie, int (*readfn) (void *,       char *, int))
    FILE *fwopen(void *cookie, int (*writefn)(void *, const char *, int));


cdef size_t fwrite_obj(void *cookie, const char *buf, int size) except 0:
    cdef PyObject* obj = <PyObject*> cookie
    return (<object> cookie).write(buf[:size])


cdef FILE* fopen_obj(PyObject* obj, str mode = "r"):

    cdef void* cookie = <void*> obj
    if mode == "w":
        return fwopen(cookie, fwrite_obj)
    else:
        raise NotImplementedError()
