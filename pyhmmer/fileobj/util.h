#ifndef _PYHMMER_FILEOBJ_UTIL
#define _PYHMMER_FILEOBJ_UTIL

static int is_cpython() {
    PyObject* impl = PySys_GetObject("implementation");
    if (impl == NULL)
        return -1;

    PyObject* name = PyObject_GetAttrString(impl, "name");
    if (name == NULL)
        return -1;
    
    if (!PyUnicode_Check(name)) {
        Py_DECREF(name);
        return -1;
    }

    int cpython = (PyUnicode_CompareWithASCIIString(name, "cpython") == 0);
    Py_DECREF(name);
    return cpython;
}

#endif