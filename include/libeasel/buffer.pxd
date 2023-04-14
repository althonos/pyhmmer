from libc.stdio cimport FILE

from libeasel cimport esl_pos_t, eslERRBUFSIZE


cdef extern from "esl_buffer.h" nogil:

    const size_t eslBUFFER_PAGESIZE
    const size_t eslBUFFER_SLURPSIZE

    cdef enum esl_buffer_mode_e:
        eslBUFFER_UNSET   = 0
        eslBUFFER_STREAM  = 1,
        eslBUFFER_CMDPIPE = 2,
        eslBUFFER_FILE    = 3,
        eslBUFFER_ALLFILE = 4,
        eslBUFFER_MMAP    = 5,
        eslBUFFER_STRING  = 6

    ctypedef struct ESL_BUFFER:
        char* mem
        esl_pos_t n
        esl_pos_t balloc
        esl_pos_t pos
        esl_pos_t baseoffset

        esl_pos_t anchor
        int nanchor

        FILE* fp
        char* filename
        char* cmdline

        esl_pos_t pagesize

        char errmsg[eslERRBUFSIZE]
        esl_buffer_mode_e mode_is

    # The ESL_BUFFER object: opening/closing
    int esl_buffer_Open      (const char *filename, const char *envvar, ESL_BUFFER **ret_bf);
    int esl_buffer_OpenFile  (const char *filename,                     ESL_BUFFER **ret_bf);
    int esl_buffer_OpenPipe  (const char *filename, const char *cmdfmt, ESL_BUFFER **ret_bf);
    int esl_buffer_OpenMem   (const char *p,         esl_pos_t  n,      ESL_BUFFER **ret_bf);
    int esl_buffer_OpenStream(FILE *fp,                                 ESL_BUFFER **ret_bf);
    int esl_buffer_Close(ESL_BUFFER *bf);

    # Positioning and anchoring an ESL_BUFFER
    esl_pos_t esl_buffer_GetOffset      (ESL_BUFFER *bf);
    int       esl_buffer_SetOffset      (ESL_BUFFER *bf, esl_pos_t offset);
    int       esl_buffer_SetAnchor      (ESL_BUFFER *bf, esl_pos_t offset);
    int       esl_buffer_SetStableAnchor(ESL_BUFFER *bf, esl_pos_t offset);
    int       esl_buffer_RaiseAnchor    (ESL_BUFFER *bf, esl_pos_t offset);

    # Raw access to the buffer
    int esl_buffer_Get(ESL_BUFFER *bf, char **ret_p, esl_pos_t *ret_n);
    int esl_buffer_Set(ESL_BUFFER *bf, char  *p,     esl_pos_t  nused);

    # Line-based parsing
    int esl_buffer_GetLine       (ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n);
    int esl_buffer_FetchLine     (ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n);
    int esl_buffer_FetchLineAsStr(ESL_BUFFER *bf, char **opt_s, esl_pos_t *opt_n);

    # Token-based parsing
    int esl_buffer_GetToken       (ESL_BUFFER *bf, const char *sep, char **opt_p, esl_pos_t *opt_n);
    int esl_buffer_FetchToken     (ESL_BUFFER *bf, const char *sep, char **opt_p, esl_pos_t *opt_n);
    int esl_buffer_FetchTokenAsStr(ESL_BUFFER *bf, const char *sep, char **opt_p, esl_pos_t *opt_n);

    # Binary (fread-like) parsing
    int esl_buffer_Read(ESL_BUFFER *bf, size_t nbytes, void *p);
