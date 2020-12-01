# coding: utf-8
"""Obtain a `FILE*` from a Python object using ``fopencookie``.
"""

# --- C imports --------------------------------------------------------------

from cpython.buffer cimport Py_buffer, PyObject_GetBuffer, PyBuffer_Release, PyBUF_READ, PyBUF_WRITE
from cpython.memoryview cimport PyMemoryView_FromMemory
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF
from libc.stdio cimport EOF, FILE
from libc.stdint cimport uint64_t
from libc.string cimport strcpy, strncpy, memcpy


# --- Linux interface --------------------------------------------------------

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


# --- fwrite implementation --------------------------------------------------

cdef ssize_t fwrite_obj(void *cookie, const char *buf, size_t size) except 0:
    """Zero-copy implementation of `fwrite` for a Python file-handle.
    """
    cdef object obj = <object> cookie
    cdef object mem = PyMemoryView_FromMemory(<char*> buf, size, PyBUF_READ)
    return obj.write(mem)


# --- fread implementations --------------------------------------------------

cdef ssize_t fread_obj_read(void *cookie, char *buf, size_t size) except -1:
    """Copying variant of `fread` for files lacking `readinto`.
    """
    cdef object    obj      = <object> cookie
    cdef object    chunk    = obj.read(size)
    cdef Py_buffer pybuffer

    if PyObject_GetBuffer(chunk, &pybuffer, PyBUF_READ) < 0:
        raise RuntimeError("could not get buffer")
    memcpy(buf, pybuffer.buf, len(chunk))
    PyBuffer_Release(&pybuffer)

    return len(chunk)


cdef ssize_t fread_obj_readinto(void *cookie, char *buf, size_t size) except -1:
    """Zero-copy implementation of `fread` using the `readinto` method.
    """
    cdef object obj   = <object> cookie
    cdef object mem

    IF SYS_IMPLEMENTATION == "pypy":
        # NB: PyPy has a bug in the `readinto` implementation that requires the
        #     memoryview to be read/write and not just write, which is why we
        #     create the memoryview in read/write mode and not just in write mode.
        mem = PyMemoryView_FromMemory(buf, size, PyBUF_READ | PyBUF_WRITE)
    ELSE:
        mem = PyMemoryView_FromMemory(buf, size, PyBUF_WRITE)
        
    return obj.readinto(mem)


# --- fseek implementation ---------------------------------------------------

cdef int fseek_obj(void* cookie, off64_t* offset, int whence) except -1:
    cdef object obj = <object> cookie
    offset[0] = obj.seek(offset[0], whence)
    return 0


# --- fclose implementation --------------------------------------------------

cdef ssize_t fclose_obj(void *cookie) except EOF:
    Py_DECREF(<object> cookie)
    return 0


# --- fopen_obj --------------------------------------------------------------

cdef FILE* fopen_obj(object obj, str mode = "r") except NULL:
    cdef cookie_io_functions_t functions
    functions.close = <cookie_close_function_t*> fclose_obj
    functions.read  = NULL
    functions.seek  = NULL
    functions.write = NULL

    try:
        if obj.readable():
            if hasattr((<object> obj), "readinto"):
                functions.read = <cookie_read_function_t*> fread_obj_readinto
            else:
                functions.read = <cookie_read_function_t*> fread_obj_read
        if obj.writable():
            functions.write = <cookie_write_function_t*> fwrite_obj
        if obj.seekable():
           functions.seek = <cookie_seek_function_t*> fseek_obj
    except AttributeError as err:
        ty = type(obj).__name__
        raise TypeError("expected `io.IOBase` instance, found {}".format(ty)) from err

    Py_INCREF(obj)
    return fopencookie(<void*> obj, mode.encode("ascii"), functions)
