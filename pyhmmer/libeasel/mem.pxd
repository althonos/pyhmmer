cdef extern from "esl_mem.h" nogil:
    int       esl_mem_strtoi32(char *p, esl_pos_t n, int base, int *opt_nc, int32_t *opt_val);
    int       esl_mem_strtoi64(char *p, esl_pos_t n, int base, int *opt_nc, int64_t *opt_val);
    int       esl_mem_strtoi  (char *p, esl_pos_t n, int base, int *opt_nc, int     *opt_val);
    int       esl_mem_strtof  (char *p, esl_pos_t n,           int *opt_nc, float   *opt_val);
    int       esl_memnewline(const char *p, esl_pos_t n, esl_pos_t *ret_nline, int *ret_nterm);
    int       esl_memtok(char **p, esl_pos_t *n, const char *delim, char **ret_tok, esl_pos_t *ret_toklen);
    esl_pos_t esl_memspn (char *p, esl_pos_t n, const char *allow);
    esl_pos_t esl_memcspn(char *p, esl_pos_t n, const char *disallow);
    int       esl_memstrcmp     (const char *p, esl_pos_t n, const char *s);
    int       esl_memstrcmp_case(const char *p, esl_pos_t n, const char *s);
    int       esl_memstrpfx     (const char *p, esl_pos_t n, const char *s);
    int       esl_memstrpfx_case(const char *p, esl_pos_t n, const char *s);
    int       esl_memstrcontains(const char *p, esl_pos_t n, const char *s);
    int       esl_memstrdup(const char *p, esl_pos_t n, char **ret_s);
    int       esl_memstrcpy(const char *p, esl_pos_t n, char *dest);
    int       esl_memtof(const char *p, esl_pos_t n, float  *ret_val);
    int       esl_memtod(const char *p, esl_pos_t n, double *ret_val);
    int       esl_mem_IsReal(const char *p, esl_pos_t n);
