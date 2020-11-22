# coding: utf-8
"""Obtain a `FILE*` from a Python object using ``fopencookie``.
"""

from cpython.ref cimport PyObject
from libc.stdio cimport FILE
from libc.stdint cimport uint64_t

cdef extern from "sys/types.h":
    ctypedef uint64_t off64_t

cdef extern from "stdio.h":

  ctypedef ssize_t (*cookie_read_function_t) (void *cookie,       char *buf, size_t size)
  ctypedef ssize_t (*cookie_write_function_t)(void *cookie, const char *buf, size_t size)
  ctypedef int     (*cookie_seek_function_t) (void *cookie, off64_t *offset, int whence);
  ctypedef int     (*cookie_close_function_t)(void *cookie)

  ctypedef struct cookie_io_functions_t:
      cookie_read_function_t*  read
      cookie_write_function_t* write
      cookie_seek_function_t*  seek
      cookie_close_function_t* close

  FILE *fopencookie(void *cookie, const char *mode, cookie_io_functions_t io_funcs);


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
