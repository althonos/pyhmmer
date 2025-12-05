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

    # 5.  Portable drop-in replacements for non-standard C functions
    extern int  esl_strcasecmp(const char *s1, const char *s2)
    extern char *esl_strsep(char **stringp, const char *delim)

    # 6. Additional string functions, esl_str*()
    # int     esl_strchop(char *s, int64_t n)
    # int     esl_strdealign(char *s, const char *aseq, const char *gapchars, int64_t *opt_rlen)
    # int     esl_str_IsBlank(char *s)
    # int     esl_str_IsInteger(char *s)
    # int     esl_str_IsReal(char *s)
    # int64_t esl_str_GetMaxWidth(char **s, int n)

    # 7. File path/name manipulation functions, including tmpfiles
    # int  esl_FileExists(const char *filename)
    # int  esl_FileTail(const char *path, int nosuffix, char **ret_file)
    # int  esl_file_Extension(char *filename, esl_pos_t n_ignore, char **ret_sfx, esl_pos_t *ret_n)
    # int  esl_FileConcat(const char *dir, const char *file, char **ret_path)
    # int  esl_FileNewSuffix(const char *filename, const char *sfx, char **ret_newpath)
    # int  esl_FileEnvOpen(const char *fname, const char *env, FILE **ret_fp, char **ret_path);
    # int  esl_tmpfile(char *basename6X, FILE **ret_fp)
    # int  esl_tmpfile_named(char *basename6X, FILE **ret_fp)
    # int  esl_getcwd(char **ret_cwd)

    # 8. Typed comparison routines.
    int esl_DCompare(double x0, double x, double r_tol, double a_tol)
    int esl_FCompare(float  x0, float  x, float  r_tol, float  a_tol)
    int esl_CCompare(char *s1, char *s2)
    int esl_DCompare_old(double a,  double b, double tol)
    int esl_FCompare_old(float  a,  float  b, float  tol)