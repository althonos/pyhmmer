from libc.stdint cimport uint8_t, int64_t


cdef extern from "stdarg.h" nogil:
    ctypedef struct va_list:
          pass

cdef extern from "easel.h" nogil:

    const size_t eslERRBUFSIZE

    const double eslCONST_E
    const double eslCONST_PI
    const double eslCONST_EULER
    const double eslCONST_GOLD
    const double eslCONST_LOG2
    const double eslCONST_LOG2R
    const double eslINFINITY

    ctypedef uint8_t ESL_DSQ
    const    ESL_DSQ eslDSQ_SENTINEL
    const    ESL_DSQ eslDSQ_ILLEGAL
    const    ESL_DSQ eslDSQ_IGNORED
    const    ESL_DSQ eslDSQ_EOL
    const    ESL_DSQ eslDSQ_EOD

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

    # 4. Improved replacements for some C library functions
    # extern int  esl_fgets(char **buf, int *n, FILE *fp);
    # extern int  esl_fprintf(FILE *fp, const char *format, ...);
    # extern int  esl_printf(const char *format, ...);
    # extern int  esl_strdup(const char *s, int64_t n, char **ret_dup);
    # extern int  esl_strcat(char **dest, int64_t ldest, const char *src, int64_t lsrc);
    # extern int  esl_strmapcat        (const ESL_DSQ *inmap, char **dest, int64_t *ldest, const char *src, esl_pos_t lsrc);
    # extern int  esl_strmapcat_noalloc(const ESL_DSQ *inmap,  char *dest, int64_t *ldest, const char *src, esl_pos_t lsrc);
    # extern int  esl_strtok    (char **s, char *delim, char **ret_tok);
    # extern int  esl_strtok_adv(char **s, char *delim, char **ret_tok, int *opt_toklen, char *opt_endchar);
    # extern int  esl_sprintf (char **ret_s, const char *format, ...);
    # extern int  esl_vsprintf(char **ret_s, const char *format, va_list *ap);
    # extern int  esl_strcmp(const char *s1, const char *s2);
