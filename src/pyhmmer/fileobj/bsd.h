#ifndef _PYHMMER_FILEOBJ_LINUX
#define _PYHMMER_FILEOBJ_LINUX

#include <stdio.h>

#include <Python.h>
#include "util.h"

#define _COOKIE_ERROR_CLOSE -1
#define _COOKIE_ERROR_WRITE -1
#define _COOKIE_ERROR_READ -1
#define _COOKIE_ERROR_SEEK -1


typedef int    (*readfn_t) (void *cookie,       char *buf, int size);
typedef int    (*writefn_t)(void *cookie, const char *buf, int size);
typedef fpos_t (*seekfn_t) (void *cookie,   fpos_t offset, int whence);
typedef int    (*closefn_t)(void *cookie);


int fileobj_bsd_write(void* cookie, const char* buf, int size) {
    PyObject* file = (PyObject*) cookie;
    
    PyObject* out = PyObject_CallMethod(file, "write", "y#", buf, (Py_ssize_t) size);
    if (out == NULL)
        return _COOKIE_ERROR_WRITE;

    if (!PyLong_Check(out)) {
        Py_DECREF(out);
        PyErr_SetString(PyExc_TypeError, "Expected int");
        return _COOKIE_ERROR_WRITE;
    }

    int n = PyLong_AsLongLong(out);
    Py_DECREF(out);
    return n;
}

int fileobj_bsd_read(void* cookie, char* buf, int size) {
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

int fileobj_bsd_readinto(void* cookie, char* buf, int size) {
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

fpos_t fileobj_bsd_seek(void* cookie, fpos_t offset, int whence) {
    PyObject* file = (PyObject*) cookie;

    PyObject* out = PyObject_CallMethod(file, "seek", "Li", offset, whence);
    if (out == NULL) 
        return _COOKIE_ERROR_SEEK;

    if (!PyLong_Check(out)) {
        Py_DECREF(out);
        PyErr_SetString(PyExc_TypeError, "Expected int");
        return _COOKIE_ERROR_SEEK;
    }

    fpos_t offset_new = PyLong_AsLongLong(out);
    Py_DECREF(out);
    return offset_new;
}

int fileobj_bsd_close(void* cookie) {
    PyObject* file = (PyObject*) cookie;
    Py_DECREF(file);
    return _COOKIE_ERROR_CLOSE;
}

FILE* fileobj_bsd_open(PyObject* obj, const char* mode) {
    Py_INCREF(obj);

    PyTypeObject* ty = Py_TYPE(obj);

    readfn_t readfn;
    writefn_t writefn;
    seekfn_t seekfn;

    PyObject* readable = PyObject_CallMethod(obj, "readable", NULL);
    if (readable == NULL)
        return NULL;
    switch (PyObject_IsTrue(readable)) {
        case 1:
            Py_DECREF(readable);
            readfn = ((is_cpython() == 1) && PyObject_HasAttrString(obj, "readinto")) ? fileobj_bsd_readinto : fileobj_bsd_read;
            break;
        case 0:
            Py_DECREF(readable);
            readfn = NULL;
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
            seekfn = fileobj_bsd_seek;
            break;
        case 0:
            Py_DECREF(seekable);
            seekfn = NULL;
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
            writefn = fileobj_bsd_write;
            break;
        case 0:
            Py_DECREF(writable);
            writefn = NULL;
            break;
        default:
            Py_DECREF(writable);
            PyErr_Format(PyExc_TypeError, "Expected `io.IOBase` instance, found %s", ty->tp_name);
            return NULL;
    }

    FILE* file = funopen(
        (void*) obj, 
        readfn,
        writefn,
        seekfn,
        fileobj_bsd_close
    );
    if (file == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to open file-like object");
        Py_DECREF(obj);
    }

    return file;
}



#endif