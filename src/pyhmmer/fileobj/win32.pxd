from cpython.unicode cimport PyUnicode_FromString

from libc.stdio cimport FILE
from libc.string cimport strcmp
from libc.stdint cimport uint32_t

import io

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

cdef extern from "minwinbase.h" nogil:
    cdef struct _OVERLAPPED:
        pass
    ctypedef _OVERLAPPED  OVERLAPPED
    ctypedef _OVERLAPPED* LPOVERLAPPED

cdef extern from "fileapi.h" nogil:
    BOOL WriteFile(HANDLE hFile, LPCVOID lpBuffer, DWORD nNumberOfBytesToWrite, LPDWORD lpNumberOfBytesWritten, LPOVERLAPPED lpOverlapped)


# cdef extern from "processthreadsapi.h" nogil:
#     HANDLE CreateThread(
#         LPSECURITY_ATTRIBUTES   lpThreadAttributes,
#         SIZE_T                  dwStackSize,
#         LPTHREAD_START_ROUTINE  lpStartAddress,
#         LPVOID                  lpParameter,
#         DWORD                   dwCreationFlags,
#         LPDWORD                 lpThreadId,
#     )



cdef extern from "io.h" nogil:
    int _open_osfhandle(size_t fd, int flags)
    int _close(int fd)

# cdef extern from "stdio.h" nogil:
#     FILE* _fdopen(int filehandle, const char* mode)
from libc.stdio cimport fdopen






# class _WinReader(threading.Thread):
#     def __init__(self, file, fd):
#         self.file = file
#         self.fd = fd

#     def run(self):

#         cdef int fdread = _open_osfhandle(<long> hReadPipe, 0)
#         # cdef FILE* fread = fdopen(fdread, "r")

#         with os.fdopen(fdread) as pipe:
#             shutil.copyfileobj(file, pipe)
       
#         _close(fdread)



cdef inline FILE* fopen_obj(object obj, const char* mode) except NULL

   


    


        # error = GetLastError()
        # success = FormatMessage(
        #     FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        #     NULL,
        #     dw,
        #     MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        #     (LPTSTR) &lpMsgBuf,
        #     0, NULL
        # )
        
        # if not success:
        #     raise RuntimeError("Failed recovering error")
        # else:
        #     raise RuntimeError(PyUnicode_FromString(<char*> lpMsgBuf))


    # raise io.UnsupportedOperation("Cannot open a file-like object on Windows")
