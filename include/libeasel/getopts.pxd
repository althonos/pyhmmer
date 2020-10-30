from libc.stdio cimport FILE


cdef extern from "esl_getopts.h" nogil:

    ctypedef struct ESL_OPTIONS:
        char* name
        int type
        char* defval
        char* envvar
        char* range
        char* toggle_opts
        char* required_opts
        char* incompat_opts
        char* help
        int docgrouptag


    ctypedef struct ESL_GETOPTS:
        ESL_OPTIONS* opt
        int nopts
        int argc
        char** argv
        int optind
        int nfiles
        char** val
        int* setby
        int* valloc
        char* optstring
        char* spoof
        char** spoof_argv


    ESL_GETOPTS *esl_getopts_Create(const ESL_OPTIONS *opt);
    ESL_GETOPTS *esl_getopts_CreateDefaultApp(const ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
    int          esl_getopts_Reuse  (ESL_GETOPTS *g);
    void         esl_getopts_Destroy(ESL_GETOPTS *g);
    void         esl_getopts_Dump(FILE *ofp, ESL_GETOPTS *g);

    int esl_opt_ProcessConfigfile (ESL_GETOPTS *g, char *filename, FILE *fp);
    int esl_opt_ProcessEnvironment(ESL_GETOPTS *g);
    int esl_opt_ProcessCmdline    (ESL_GETOPTS *g, int argc, char **argv);
    int esl_opt_ProcessSpoof      (ESL_GETOPTS *g, const char *cmdline);
    int esl_opt_VerifyConfig      (ESL_GETOPTS *g);
    int esl_opt_ArgNumber   (const ESL_GETOPTS *g);
    int esl_opt_SpoofCmdline(const ESL_GETOPTS *g, char **ret_cmdline);

    int esl_opt_GetSetter(const ESL_GETOPTS *g, char *optname);

    int    esl_opt_IsDefault (const ESL_GETOPTS *g, char *optname);
    int    esl_opt_IsOn      (const ESL_GETOPTS *g, char *optname);
    int    esl_opt_IsUsed    (const ESL_GETOPTS *g, char *optname);

    int    esl_opt_GetBoolean(const ESL_GETOPTS *g, char *optname);
    int    esl_opt_GetInteger(const ESL_GETOPTS *g, char *optname);
    double esl_opt_GetReal   (const ESL_GETOPTS *g, char *optname);
    char   esl_opt_GetChar   (const ESL_GETOPTS *g, char *optname);
    char  *esl_opt_GetString (const ESL_GETOPTS *g, char *optname);
    char  *esl_opt_GetArg    (const ESL_GETOPTS *g, int which);

    int esl_opt_DisplayHelp(FILE *ofp, const ESL_GETOPTS *go, int docgroup, int indent, int textwidth);
