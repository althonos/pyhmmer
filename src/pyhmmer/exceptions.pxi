# coding: utf-8
"""Custom exception handler for Easel, that avoids system exits.
"""

# --- C imports --------------------------------------------------------------

from cpython.exc cimport PyErr_Clear, PyErr_Occurred, PyErr_Fetch, PyErr_Restore
from cpython.ref cimport PyObject

cdef extern from "stdarg.h" nogil:
    ctypedef struct va_list:
        pass

cdef extern from "stdio.h" nogil:
    cdef int vsprintf(char* s, const char* format, va_list arg)


# --- Python imports ---------------------------------------------------------

from .errors import EaselError


# --- Error handler ----------------------------------------------------------

cdef object _recover_error():
    cdef PyObject*  type
    cdef PyObject*  value
    cdef PyObject*  traceback
    if PyErr_Occurred():
        PyErr_Fetch(&type, &value, &traceback)
        if isinstance(<object> value, Exception):
            return <object> value
        else:
            return (<object> type)(<object> value)
    else:
        return None

cdef void _reraise_error() except *:
    cdef object error = _recover_error()
    if error is not None:
        raise error

# We register a custom exception handler so that the program should not crash,
# but simply raise an exception, while chaining
cdef void _py_handler(int errcode, int use_errno, char* sourcefile, int sourceline, char* format, va_list argp) except * nogil:
    cdef char[2048] errbuf
    cdef int        errlength

    with gil:
        # recover the internal error, if any
        error = _recover_error()

        # get the error message if possible
        errlength = vsprintf(<char*> &errbuf, format, argp)
        if errlength > 0:
            message = errbuf[:errlength].decode('utf-8', errors="replace")
        else:
            message = None

        # raise the exception and hope it's collected
        raise EaselError(errcode, message) from error

libeasel.esl_exception_SetHandler(<libeasel.esl_exception_handler_f> _py_handler)
