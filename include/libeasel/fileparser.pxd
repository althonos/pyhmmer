from libc.stdio cimport FILE

from libeasel cimport eslERRBUFSIZE


cdef extern from "esl_fileparser.h" nogil:

    ctypedef struct ESL_FILEPARSER:
        FILE* fp
        char* buf
        int buflen
        char* s
        char commentchar

        char* filename
        int linenumber
        char[eslERRBUFSIZE] errbuf

        int is_buffer
        char* mem_buffer
        int mem_size
        int mem_pos


    int esl_fileparser_Open(const char *filename, const char *envvar, ESL_FILEPARSER **ret_efp) nogil
    ESL_FILEPARSER *esl_fileparser_Create(FILE *fp) except *
    ESL_FILEPARSER *esl_fileparser_CreateMapped(const void *buffer, int size) nogil
    int  esl_fileparser_SetCommentChar  (ESL_FILEPARSER *efp, char c) nogil
    int  esl_fileparser_GetToken        (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen) except *
    int  esl_fileparser_NextLine        (ESL_FILEPARSER *efp) except *
    int  esl_fileparser_NextLinePeeked  (ESL_FILEPARSER *efp, char *prefix, int plen) except *
    int  esl_fileparser_GetTokenOnLine  (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen) nogil
    int  esl_fileparser_GetRemainingLine(ESL_FILEPARSER *efp, char **ret_s) nogil
    void esl_fileparser_Destroy         (ESL_FILEPARSER *efp) nogil
    void esl_fileparser_Close           (ESL_FILEPARSER *efp) except *
