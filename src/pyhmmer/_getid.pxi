# --- Patch for PyPy 3.9 -----------------------------------------------------

cdef extern from *:
    """
    #ifndef HAS_PYINTERPRETERSTATE_GETID
    int64_t PyInterpreterState_GetID(PyInterpreterState *interp) {
        return 0;
    }
    #endif
    """