# coding: utf-8
"""Obtain a `FILE*` from a Python object using ``funopen``.
"""

# --- C imports --------------------------------------------------------------

from cpython.buffer cimport Py_buffer, PyObject_GetBuffer, PyBuffer_Release, PyBUF_READ, PyBUF_WRITE
from cpython.memoryview cimport PyMemoryView_FromMemory
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF
from libc.stdio cimport EOF, FILE
from libc.stdint cimport int64_t, uint64_t
from libc.string cimport strcpy, strncpy, memcpy


# --- BSD interface ----------------------------------------------------------

cdef extern from "stdio.h":
    ctypedef long fpos_t

ctypedef int    (*readfn_t) (void *cookie,       char *buf, int size)
ctypedef int    (*writefn_t)(void *cookie, const char *buf, int size)
ctypedef fpos_t (*seekfn_t) (void *cookie,   fpos_t offset, int whence)
ctypedef int    (*closefn_t)(void *cookie)

cdef extern from "stdio.h":
    FILE* funopen(
        const void*  cookie,
        readfn_t*    readfn,
        writefn_t*   writefn,
        seekfn_t*    seekfn,
        closefn_t*   closefn
    )


# --- fwrite implementation --------------------------------------------------

cdef int fwrite_obj(void *cookie, const char *buf, int size) except 0:
    """Zero-copy implementation of `fwrite` for a Python file-handle.
    """
    cdef object obj = <object> cookie
    cdef object mem = PyMemoryView_FromMemory(<char*> buf, size, PyBUF_READ)
    return obj.write(mem)


# --- fread implementations --------------------------------------------------

cdef int fread_obj_read(void *cookie, char *buf, int size) except -1:
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


cdef int fread_obj_readinto(void *cookie, char *buf, int size) except -1:
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

cdef fpos_t fseek_obj(void* cookie, fpos_t offset, int whence):
    cdef object obj = <object> cookie
    return obj.seek(offset, whence)


# --- fclose implementation --------------------------------------------------

cdef int fclose_obj(void *cookie) except EOF:
    Py_DECREF(<object> cookie)
    return 0


# --- fopen_obj --------------------------------------------------------------

cdef FILE* fopen_obj(object obj, str mode = "r") except NULL:
    cdef closefn_t* closefn = <closefn_t*> fclose_obj
    cdef readfn_t*  readfn  = NULL
    cdef writefn_t* writefn = NULL
    cdef seekfn_t*  seekfn  = NULL

    try:
        if obj.readable():
            if hasattr((<object> obj), "readinto"):
                readfn = <readfn_t*> fread_obj_readinto
            else:
                readfn = <readfn_t*> fread_obj_read
        if obj.writable():
            writefn = <writefn_t*> fwrite_obj
        if obj.seekable():
            seekfn = <seekfn_t*> fseek_obj
    except AttributeError as err:
        ty = type(obj).__name__
        raise TypeError("expected `io.IOBase` instance, found {}".format(ty)) from err

    Py_INCREF(obj)
    return funopen(<void*> obj, readfn, writefn, seekfn, closefn)
