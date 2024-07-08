

cdef extern from "esl_json.h" nogil:

    cdef enum esl_json_type_e:
        eslJSON_UNKNOWN = 0
        eslJSON_OBJECT  = 1
        eslJSON_ARRAY   = 2
        eslJSON_KEY     = 3
        eslJSON_STRING  = 4
        eslJSON_NUMBER  = 5
        eslJSON_BOOLEAN = 6
        eslJSON_NULL    = 7

    cdef enum esl_json_state_e:
        eslJSON_OBJ_NONE      = 0
        eslJSON_OBJ_OPEN      = 1
        eslJSON_OBJ_COLON     = 2
        eslJSON_OBJ_COMMA     = 3
        eslJSON_ARR_OPEN      = 4
        eslJSON_ARR_COMMA     = 5
        eslJSON_STR_OPEN      = 6
        eslJSON_STR_CHAR      = 7
        eslJSON_STR_BACKSLASH = 8
        eslJSON_STR_PROTECTED = 9
        eslJSON_STR_UNICODE   = 10
        eslJSON_KEY_OPEN      = 11
        eslJSON_KEY_CHAR      = 12
        eslJSON_KEY_BACKSLASH = 13
        eslJSON_KEY_PROTECTED = 14
        eslJSON_KEY_UNICODE   = 15
        eslJSON_NUM_SIGN      = 16
        eslJSON_NUM_ZERO      = 17
        eslJSON_NUM_NONZERO   = 18
        eslJSON_NUM_LEADDIGIT = 19
        eslJSON_NUM_POINT     = 20
        eslJSON_NUM_FRACDIGIT = 21
        eslJSON_NUM_EXP       = 22
        eslJSON_NUM_EXPSIGN   = 23
        eslJSON_NUM_EXPDIGIT  = 24
        eslJSON_VAL_TRUE      = 25
        eslJSON_VAL_FALSE     = 26
        eslJSON_VAL_NULL      = 27
        eslJSON_VAL_INOBJ     = 28
        eslJSON_VAL_INARR     = 29
        eslJSON_STR_ASKEY     = 30

    ctypedef struct ESL_JSON_TOK:
        esl_json_type_e type
        esl_pos_t startpos
        esl_pos_t endpos
        int nchild
        int firstchild
        int lastchild
        int nextsib
        int linenum
        int linepos

    ctypedef struct ESL_JSON:
        ESL_JSON_TOK* tok
        int ntok
        int nalloc
        int redline

    ctypedef struct ESL_JSON_PARSER:
        esl_json_state_e state
        ESL_STACK* pda
        int curridx
        int codelen
        esl_pos_t pos
        int linenum
        int linepos

    # Full and incremental JSON parsing
    int esl_json_Parse(ESL_BUFFER *bf, ESL_JSON **ret_pi);
    int esl_json_PartialParse(ESL_JSON_PARSER *parser, ESL_JSON *pi, const char *s, esl_pos_t n, esl_pos_t *ret_nused, char *errbuf);

    # ESL_JSON
    ESL_JSON *esl_json_Create   (void);
    int       esl_json_Grow     (ESL_JSON *pi);
    size_t    esl_json_Sizeof   (ESL_JSON *pi);
    size_t    esl_json_MinSizeof(ESL_JSON *pi);
    int       esl_json_Reuse    (ESL_JSON *pi);
    void      esl_json_Destroy  (ESL_JSON *pi);

    # ESL_JSON_PARSER
    ESL_JSON_PARSER *esl_json_parser_Create(void);
    void             esl_json_parser_Destroy(ESL_JSON_PARSER *parser);

    # Accessing tokenized data
    char      *esl_json_GetMem   (const ESL_JSON *pi, int idx, const ESL_BUFFER *bf);
    esl_pos_t  esl_json_GetLen   (const ESL_JSON *pi, int idx, const ESL_BUFFER *bf);
    int        esl_json_ReadInt  (const ESL_JSON *pi, int idx,       ESL_BUFFER *bf, int   *ret_i);
    int        esl_json_ReadFloat(const ESL_JSON *pi, int idx,       ESL_BUFFER *bf, float *ret_x);

    # Debugging, development 
    int   esl_json_Validate(const ESL_JSON *pi, const ESL_BUFFER *bf, char *errbuf);
    char *esl_json_DecodeType(enum esl_json_type_e type);
    int   esl_json_Dump(FILE *fp, ESL_JSON *pi);
    int   esl_json_SampleDirty(ESL_RANDOMNESS *rng, char **ret_s, int *ret_n);
