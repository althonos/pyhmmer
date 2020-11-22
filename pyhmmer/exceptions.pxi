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

# We register a custom exception handler so that the program should not crash,
# but simply raise an exception, while chaining
cdef void py_handler(int errcode, int use_errno, char* sourcefile, int sourceline, char* format, va_list argp) nogil except *:
    cdef PyObject*  type
    cdef PyObject*  value
    cdef PyObject*  traceback
    cdef char[2048] errbuf
    cdef int        errlength

    with gil:
        # recover the internal error, if any
        if PyErr_Occurred():
            PyErr_Fetch(&type, &value, &traceback)
            error = (<object> type)(<object> value)
            error.__traceback__ = <object> traceback
        else:
            error = None

        # get the error message if possible
        errlength = vsprintf(<char*> &errbuf, format, argp)
        if errlength > 0:
            message = errbuf[:errlength].decode('utf-8', errors="replace")
        else:
            message = None

        # raise the exception and hope it's collected
        raise EaselError(errcode, message) from error

libeasel.esl_exception_SetHandler(<libeasel.esl_exception_handler_f> py_handler)
