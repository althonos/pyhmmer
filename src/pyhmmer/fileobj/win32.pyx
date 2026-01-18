# coding: utf-8
# cython: language_level=3
"""File object emulation for Windows.
"""

# --- C imports --------------------------------------------------------------

from cpython.bytearray cimport PyByteArray_AsString
from cpython.bytes cimport PyBytes_AsStringAndSize

cimport libc.errno
from libc.stdio cimport FILE, ferror, feof, fwrite, fread, fclose, printf, fflush, fdopen
from libc.stdint cimport uint32_t, intptr_t
from libc.string cimport strcmp

cimport libeasel

from ..errors import EaselError

# --- Python imports ---------------------------------------------------------

import threading
import os
import io
import shutil
import time

# --- WinApi definitions -----------------------------------------------------

cdef extern from "stdarg.h":
    ctypedef struct va_list

cdef extern from "windows.h" nogil:
    ctypedef size_t   SIZE_T
    ctypedef bint     BOOL
    ctypedef void*    LPVOID
    ctypedef const void* LPCVOID
    ctypedef void*    HANDLE
    ctypedef HANDLE*  PHANDLE
    ctypedef uint32_t DWORD
    ctypedef DWORD*   LPDWORD
    ctypedef char*    LPTSTR

    cdef struct _SECURITY_ATTRIBUTES:
        DWORD  nLength
        void*  lpSecurityDescriptor
        BOOL   bInheritHandle
    ctypedef _SECURITY_ATTRIBUTES  SECURITY_ATTRIBUTES
    ctypedef _SECURITY_ATTRIBUTES* PSECURITY_ATTRIBUTES
    ctypedef _SECURITY_ATTRIBUTES* LPSECURITY_ATTRIBUTES

cdef extern from "errhandlingapi.h" nogil:
    DWORD GetLastError()

cdef extern from "winbase.h" nogil:
    DWORD FormatMessage(DWORD dwFlags, void* lpSource, DWORD dwMessageId, DWORD dwLanguageId, char* lpBuffer, DWORD nSize, va_list* Arguments)

    cdef enum:
        FORMAT_MESSAGE_ALLOCATE_BUFFER
        FORMAT_MESSAGE_ARGUMENT_ARRAY
        FORMAT_MESSAGE_FROM_HMODULE
        FORMAT_MESSAGE_FROM_STRING
        FORMAT_MESSAGE_FROM_SYSTEM
        FORMAT_MESSAGE_IGNORE_INSERTS

cdef extern from "handleapi.h" nogil:
    BOOL CloseHandle(HANDLE hObject)

cdef extern from "namedpipeapi.h" nogil:
    BOOL CreatePipe(PHANDLE hReadPipe, PHANDLE hWritePipe, LPSECURITY_ATTRIBUTES lpPipeAttributes, DWORD nSize)
    BOOL SetNamedPipeHandleState(HANDLE  hNamedPipe, LPDWORD lpMode, LPDWORD lpMaxCollectionCount, LPDWORD lpCollectDataTimeout)

    enum:
        PIPE_WAIT
        PIPE_NOWAIT

cdef extern from "ioapiset.h" nogil:
    enum:
        ERROR_IO_PENDING

cdef extern from "minwinbase.h" nogil:
    cdef struct _OVERLAPPED:
        pass
    ctypedef _OVERLAPPED  OVERLAPPED
    ctypedef _OVERLAPPED* LPOVERLAPPED

cdef extern from "fileapi.h" nogil:
    BOOL WriteFile(HANDLE hFile, LPCVOID lpBuffer, DWORD nNumberOfBytesToWrite, LPDWORD lpNumberOfBytesWritten, LPOVERLAPPED lpOverlapped)

cdef extern from "io.h" nogil:
    int _open_osfhandle(size_t fd, int flags)
    int _close(int fd)


# --- Reader code ------------------------------------------------------------

class _WinReader(threading.Thread):
    daemon = True

    def __init__(self, file, handle):
        super().__init__()
        self.file = file
        self.handle = handle
        self.ready = threading.Event()
        self.ready.clear()

    def run(self):
        cdef int       fdwrite
        cdef bytes     out
        cdef int       error    = 0
        cdef FILE*     f        = NULL
        cdef bytearray data     = bytearray(io.DEFAULT_BUFFER_SIZE)
        cdef bint      readinto = hasattr(self.file, "readinto")
        cdef ssize_t   length   = 1
        cdef char*     b        = PyByteArray_AsString(data)

        # Get a CRT file descriptor from the pipe HANDLE
        fdwrite = _open_osfhandle(<intptr_t> self.handle, 0)
        if fdwrite == -1:
            CloseHandle(<HANDLE> self.handle)
            raise RuntimeError("Failed opening file descriptor")

        # Get a FILE* from the CRF file descriptor
        f = fdopen(fdwrite, "wb")
        if f == NULL:
            _close(fdwrite) # ensure the file descriptor is still closed
            raise RuntimeError("Failed opening pipe")

        try:
            with nogil:
                while True:
                    # Acquire the GIL just long enough to read from the file-like
                    # object
                    with gil:
                        # Use readinto if available, otherwise fallback to using
                        # read which will create a new bytes object
                        if readinto:
                            length = self.file.readinto(data)
                        else:
                            out = self.file.read(io.DEFAULT_BUFFER_SIZE)
                            PyBytes_AsStringAndSize(out, &b, &length)
                        # Having obtained the first chunk, we can now let the consumer
                        # read the pipe, otherwise there is a chance of deadlock
                        # because we could not acquire the GIL to read the data while
                        # the consumer would already be blocking and reading from the pipe
                        # if length > 0:
                        self.ready.set()
                    if length <= 0:
                        break
                    # Write data from the file into the pipe
                    if length > 0:
                        # Write to the pipe while the GIL is not held, so that if
                        # blocking happens because the pipe is full, the code does
                        # not lock
                        w = fwrite(b, sizeof(char), <size_t> length, f)
                        if w < length:
                            fclose(f)
                            # don't error on EPIPE, it means the main thread closed
                            # this file so it shouldn't be treated as an error here.
                            if libc.errno.errno == libc.errno.EPIPE:
                                break
                            else:
                                raise OSError(libc.errno.errno)
                        # Need to flush, otherwise we may try to re-acquire the GIL
                        # before all data has been written to the pipe, and therefore
                        # end up with a deadlock
                        fflush(f)
        finally:
            fclose(f)
            # No need to call _close on the file descriptor, as fclose will also close it
            # (per https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/open-osfhandle?view=msvc-170#remarks)


cdef class _WinSynchronizedReader:
    
    def __init__(self, object fileobj):

        cdef HANDLE hReadPipe
        cdef HANDLE hWritePipe
        cdef int    fdread
        cdef FILE*  fread

        success = CreatePipe(&hReadPipe, &hWritePipe, NULL, 0)
        if not success:
            raise RuntimeError("Failed creating a pipe")

        fdread = _open_osfhandle(<size_t> hReadPipe, 0)
        if fdread == -1:
            CloseHandle(hReadPipe)
            CloseHandle(hWritePipe)
            raise RuntimeError("Failed opening file descriptor")

        fread  = fdopen(fdread, "rb")
        if fread == NULL:
            _close(fdread)
            CloseHandle(hWritePipe)
            raise RuntimeError("Failed opening file")

        self.file = fread
        self.thread = _WinReader(fileobj, <size_t> hWritePipe)
        self.thread.start()
        self.thread.ready.wait()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    cpdef object close(self):
        if self.file != NULL:
            fclose(self.file)
       

cdef FILE* fopen_obj(object obj, const char* mode) except NULL:
    if strcmp(mode, "r") != 0:
        raise ValueError("invalid mode")

    cdef _WinSynchronizedReader r = _WinSynchronizedReader(obj)
    return r.file

# --- Writer code ------------------------------------------------------------

class _WinWriter(threading.Thread):
    daemon = True

    def __init__(self, file, handle):
        super().__init__()
        self.file = file
        self.handle = handle
        self.ready = threading.Barrier(2)
        self.done = threading.Event()
        self.error = None

    def run(self):
        cdef int        fdread
        cdef int        error    = 0
        cdef FILE*      f        = NULL
        cdef bytearray  data     = bytearray(io.DEFAULT_BUFFER_SIZE)
        cdef memoryview mem     = memoryview(data)
        cdef size_t     datalen  = io.DEFAULT_BUFFER_SIZE
        cdef bint       readinto = hasattr(self.file, "readinto")
        cdef ssize_t    length   = -1
        cdef ssize_t    n        = 0
        cdef char*      b        = PyByteArray_AsString(data)

        # Get a CRT file descriptor from the pipe HANDLE
        fdread = _open_osfhandle(<size_t> self.handle, 0)
        if fdread == -1:
            CloseHandle(<HANDLE> self.handle)
            raise RuntimeError("Failed opening file descriptor")

        # Get a FILE* from the CRF file descriptor
        f = fdopen(fdread, "rb")
        if f == NULL:
            _close(fdread) # ensure the file descriptor is still closed
            raise RuntimeError("Failed opening pipe")

        try:
            with nogil:
                while True:
                    length = fread(b, sizeof(char), datalen, f)
                    if length == 0:
                        if ferror(f):
                            raise OSError(libc.errno.errno)
                        break
                    with gil:
                        n = 0
                        while n < length:
                            try:
                                n += self.file.write(mem[n:length])
                            except Exception as err:
                                self.error = err
                                break
        finally:
            fclose(f)
            self.done.set()


cdef class _WinSynchronizedWriter:

    def __init__(self, object fileobj):

        cdef HANDLE hReadPipe
        cdef HANDLE hWritePipe
        cdef int    fdread
        cdef FILE*  fread
        cdef DWORD  mode       = PIPE_NOWAIT

        success = CreatePipe(&hReadPipe, &hWritePipe, NULL, 64)
        if not success:
            raise RuntimeError("Failed creating a pipe")

        # SetNamedPipeHandleState(hReadPipe, &mode, NULL, NULL)
        fdwrite = _open_osfhandle(<intptr_t> hWritePipe, 0)
        if fdwrite == -1:
            CloseHandle(hReadPipe)
            CloseHandle(hWritePipe)
            raise RuntimeError("Failed opening file descriptor")

        fwrite  = fdopen(fdwrite, "wb")
        if fwrite == NULL:
            _close(fdwrite)
            CloseHandle(hReadPipe)
            raise RuntimeError("Failed opening file")

        self.file   = fwrite
        self.thread = _WinWriter(fileobj, <intptr_t> hReadPipe)
        self.thread.start()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    cpdef object close(self):
        if self.file != NULL:
            fclose(self.file)
        self.thread.done.wait()
        if self.thread.error is not None:
            raise EaselError(libeasel.eslEWRITE, "write failed") from self.thread.error
        return None

