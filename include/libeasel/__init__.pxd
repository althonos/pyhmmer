from libc.stdint cimport uint8_t, int64_t


cdef extern from "stdarg.h" nogil:
    ctypedef struct va_list:
          pass

cdef extern from "easel.h" nogil:

    cdef size_t eslERRBUFSIZE

    cdef double eslCONST_E
    cdef double eslCONST_PI
    cdef double eslCONST_EULER
    cdef double eslCONST_GOLD
    cdef double eslCONST_LOG2
    cdef double eslCONST_LOG2R
    cdef double eslINFINITY

    ctypedef uint8_t ESL_DSQ
    cdef     ESL_DSQ eslDSQ_SENTINEL
    cdef     ESL_DSQ eslDSQ_ILLEGAL
    cdef     ESL_DSQ eslDSQ_IGNORED
    cdef     ESL_DSQ eslDSQ_EOL
    cdef     ESL_DSQ eslDSQ_EOD

    ctypedef int64_t esl_pos_t

    cdef enum:
        eslOK = 0
        eslFAIL = 1
        eslEOL = 2
        eslEOF = 3
        eslEOD = 4
        eslEMEM = 5
        eslENOTFOUND = 6
        eslEFORMAT = 7
        eslEAMBIGUOUS = 8
        eslEDIVZERO = 9
        eslEINCOMPAT = 10
        eslEINVAL = 11
        eslESYS = 12
        eslECORRUPT = 13
        eslEINCONCEIVABLE = 14
        eslESYNTAX = 15
        eslERANGE = 16
        eslEDUP = 17
        eslENOHALT = 18
        eslENORESULT = 19
        eslENODATA = 20
        eslETYPE = 21
        eslEOVERWRITE = 22
        eslENOSPACE = 23
        eslEUNIMPLEMENTED = 24
        eslENOFORMAT = 25
        eslENOALPHABET = 26
        eslEWRITE = 27
        eslEINACCURATE = 28

    # 1. Error handling
    ctypedef void(*esl_exception_handler_f)(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, va_list argp)
    # void esl_fail(char *errbuf, const char *format, ...)
    # void esl_exception(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, ...);
    void esl_exception_SetHandler(esl_exception_handler_f)
    void esl_exception_ResetDefaultHandler()
    # void esl_nonfatal_handler(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, va_list argp);
    # void esl_fatal(const char *format, ...) ESL_ATTRIBUTE_NORETURN;
