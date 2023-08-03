#ifndef _PYHMMER_FILEOBJ_LINUX
#define _PYHMMER_FILEOBJ_LINUX

#include <stdio.h>
#include <sys/types.h>

#include <Python.h>
#include "util.h"

#define _COOKIE_ERROR_CLOSE EOF
#define _COOKIE_ERROR_WRITE 0
#define _COOKIE_ERROR_READ -1
#define _COOKIE_ERROR_SEEK -1


Py_ssize_t fileobj_linux_write(void* cookie, const char* buf, size_t size) {
    PyObject* file = (PyObject*) cookie;
    
    PyObject* out = PyObject_CallMethod(file, "write", "y#", buf, (Py_ssize_t) size);
    if (out == NULL)
        return _COOKIE_ERROR_WRITE;

    if (!PyLong_Check(out)) {
        Py_DECREF(out);
        PyErr_SetString(PyExc_TypeError, "Expected int");
        return _COOKIE_ERROR_WRITE;
    }

    Py_ssize_t n = PyLong_AsSize_t(out);
    Py_DECREF(out);
    return n;
}

Py_ssize_t fileobj_linux_read(void* cookie, char* buf, size_t size) {
    PyObject* file = (PyObject*) cookie;

    PyObject* chunk = PyObject_CallMethod(file, "read", "n", size);
    if (chunk == NULL)
        return _COOKIE_ERROR_READ;

    const char* data = PyBytes_AsString(chunk);
    if (data == NULL) {
        Py_DECREF(chunk);
        return _COOKIE_ERROR_READ;
    }

    Py_ssize_t len = PyBytes_Size(chunk);
    if (len > size) {
        Py_DECREF(chunk);
        PyErr_SetString(PyExc_BufferError, "buffer too small to store `read` result");
        return _COOKIE_ERROR_READ;
    }

    memcpy(buf, data, len);

    Py_DECREF(chunk);
    return len;
}

Py_ssize_t fileobj_linux_readinto(void* cookie, char* buf, size_t size) {
    PyObject* file = (PyObject*) cookie;

    PyObject* mem  = PyMemoryView_FromMemory(buf, (Py_ssize_t) size, PyBUF_WRITE);
    if (mem == NULL) 
        return _COOKIE_ERROR_READ;

    PyObject* out = PyObject_CallMethod(file, "readinto", "O", mem);
    if (out == NULL) {
        Py_DECREF(mem);
        return _COOKIE_ERROR_READ;
    }

    if (!PyLong_Check(out)) {
        Py_DECREF(out);
        Py_DECREF(mem);
        PyErr_SetString(PyExc_TypeError, "Expected int");
        return _COOKIE_ERROR_WRITE;
    }    

    Py_ssize_t len = PyLong_AsSize_t(out);
    Py_DECREF(out);
    Py_DECREF(mem);
    return len;
}

int fileobj_linux_seek(void* cookie, off64_t* offset, int whence) {
    PyObject* file = (PyObject*) cookie;

    PyObject* out = PyObject_CallMethod(file, "seek", "Li", *offset, whence);
    if (out == NULL) 
        return _COOKIE_ERROR_SEEK;

    if (!PyLong_Check(out)) {
        Py_DECREF(out);
        PyErr_SetString(PyExc_TypeError, "Expected int");
        return _COOKIE_ERROR_SEEK;
    }

    *offset = PyLong_AsLongLong(out);
    Py_DECREF(out);
    return 0;
}

int fileobj_linux_close(void* cookie) {
    PyObject* file = (PyObject*) cookie;
    Py_DECREF(file);
    return _COOKIE_ERROR_CLOSE;
}

FILE* fileobj_linux_open(PyObject* obj, const char* mode) {
    Py_INCREF(obj);

    PyTypeObject* ty = Py_TYPE(obj);

    cookie_io_functions_t functions;
    functions.close = fileobj_linux_close;

    PyObject* readable = PyObject_CallMethod(obj, "readable", NULL);
    if (readable == NULL)
        return NULL;
    switch (PyObject_IsTrue(readable)) {
        case 1:
            Py_DECREF(readable);
            functions.read = ((is_cpython() == 1) && PyObject_HasAttrString(obj, "readinto")) ? fileobj_linux_readinto : fileobj_linux_read;
            break;
        case 0:
            Py_DECREF(readable);
            functions.read = NULL;
            break;
        default:
            Py_DECREF(readable);
            PyErr_Format(PyExc_TypeError, "Expected `io.IOBase` instance, found %s", ty->tp_name);
            return NULL;
    }

    PyObject* seekable = PyObject_CallMethod(obj, "seekable", NULL);
    if (seekable == NULL)
        return NULL;
    switch (PyObject_IsTrue(seekable)) {
        case 1: 
            Py_DECREF(seekable);
            functions.seek = fileobj_linux_seek;
            break;
        case 0:
            Py_DECREF(seekable);
            functions.seek = NULL;
            break;
        default:
            Py_DECREF(seekable);
            PyErr_Format(PyExc_TypeError, "Expected `io.IOBase` instance, found %s", ty->tp_name);
            return NULL;
    }
    
    PyObject* writable = PyObject_CallMethod(obj, "writable", NULL);
    if (writable == NULL)
        return NULL;
    switch (PyObject_IsTrue(writable)) {
        case 1: 
            Py_DECREF(writable);
            functions.write = fileobj_linux_write;
            break;
        case 0:
            Py_DECREF(writable);
            functions.write = NULL;
            break;
        default:
            Py_DECREF(writable);
            PyErr_Format(PyExc_TypeError, "Expected `io.IOBase` instance, found %s", ty->tp_name);
            return NULL;
    }

    FILE* file = fopencookie((void*) obj, mode, functions);
    if (file == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to open file-like object");
        Py_DECREF(obj);
    }

    return file;
}



#endif