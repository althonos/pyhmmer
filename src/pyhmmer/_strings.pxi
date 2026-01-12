from cpython.unicode cimport PyUnicode_AsUTF8AndSize, PyUnicode_FromString
from cpython.exc cimport PyErr_WarnEx

class _strwrapper(str):

    def __eq__(self, other):
        if isinstance(other, bytes):
            PyErr_WarnEx(
                DeprecationWarning, 
                (
                    "With PyHMMER v0.12.0 most methods and properties returning "
                    "`bytes` for textual data now return `str`, use string comparison "
                    "to make your code forward compatible."
                ), 
                1
            )
            return self.encode() == other
        return super().__eq__(other)

    def __ne__(self, other):
        if isinstance(other, bytes):
            PyErr_WarnEx(
                DeprecationWarning, 
                (
                    "With PyHMMER v0.12.0 most methods and properties returning "
                    "`bytes` for textual data now return `str`, use string comparison "
                    "to make your code forward compatible."
                ), 
                1
            )
            return self.encode().__ne__(other)
        return super().__ne__(other)

    def __hash__(self):
        return super().__hash__()

    def decode(self, encoding='utf-8', errors='strict'):
        PyErr_WarnEx(
            DeprecationWarning, 
            (
                "With PyHMMER v0.12.0 most methods and properties returning "
                "`bytes` for textual data now return `str`, remove this call "
                "to make your code forward compatible."
            ), 
            1
        )
        return self

ctypedef int (*setter_t)(void*, const char*) noexcept nogil 
ctypedef int (*setter_pos_t)(void*, const char*, libeasel.esl_pos_t) noexcept nogil 

cdef object _get_str(const char* ptr):
    return _strwrapper(PyUnicode_FromString(ptr))

cdef int _set_str(
    void* dat, 
    object obj, 
    setter_t setter,
) except -1:
    # FIXME: check for internal NULL bytes
    cdef       int                status
    cdef const unsigned char[::1] view
    cdef       ssize_t            size   = -1
    cdef const char*              text   = NULL

    if isinstance(obj, str):
        text = PyUnicode_AsUTF8AndSize(obj, &size)
    else:
        view = obj
        size = view.shape[0]
        if size > 0:
            text = <const char*> &view[0]

    with nogil:
        status = setter(dat, text)
    return status
