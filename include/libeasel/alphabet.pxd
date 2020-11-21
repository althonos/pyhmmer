from libc.stdint cimport int64_t

from libeasel cimport ESL_DSQ, esl_pos_t


cdef extern from "esl_alphabet.h" nogil:

    cdef enum:
        eslUNKNOWN = 0
        eslRNA = 1
        eslDNA = 2
        eslAMINO = 3
        eslCOINS = 4
        eslDICE = 5
        eslNONSTANDARD = 6

    ctypedef struct ESL_ALPHABET:
        int type
        int K
        int Kp
        char* sym
        ESL_DSQ[128] inmap
        char** degen
        int* ndegen
        ESL_DSQ* complement

    # ESL alphabet object
    ESL_ALPHABET *esl_alphabet_Create(int type);
    ESL_ALPHABET *esl_alphabet_CreateCustom(const char *alphabet, int K, int Kp);
    int           esl_alphabet_SetEquiv(ESL_ALPHABET *a, char sym, char c);
    int           esl_alphabet_SetCaseInsensitive(ESL_ALPHABET *a);
    int           esl_alphabet_SetDegeneracy(ESL_ALPHABET *a, char c, char *ds);
    int           esl_alphabet_SetIgnored(ESL_ALPHABET *a, const char *ignoredchars);
    size_t        esl_alphabet_Sizeof(ESL_ALPHABET *a);
    void          esl_alphabet_Destroy(ESL_ALPHABET *a);

    # Digitized sequences
    int     esl_abc_CreateDsq(const ESL_ALPHABET *a, const char    *seq,        ESL_DSQ **ret_dsq);
    int     esl_abc_Digitize (const ESL_ALPHABET *a, const char    *seq,        ESL_DSQ *dsq);
    int     esl_abc_Textize  (const ESL_ALPHABET *a, const ESL_DSQ *dsq,  int64_t L, char   *seq);
    int     esl_abc_TextizeN (const ESL_ALPHABET *a, const ESL_DSQ *dptr, int64_t L, char   *buf);
    int     esl_abc_dsqcpy(const ESL_DSQ *dsq, int64_t L, ESL_DSQ *dcopy);
    int     esl_abc_dsqdup(const ESL_DSQ *dsq, int64_t L, ESL_DSQ **ret_dup);
    int     esl_abc_dsqcat        (const ESL_DSQ *inmap, ESL_DSQ **dsq, int64_t *L, const char *s, esl_pos_t n);
    int     esl_abc_dsqcat_noalloc(const ESL_DSQ *inmap, ESL_DSQ  *dsq, int64_t *L, const char *s, esl_pos_t n);
    int64_t esl_abc_dsqlen(const ESL_DSQ *dsq);
    int64_t esl_abc_dsqrlen(const ESL_ALPHABET *a, const ESL_DSQ *dsq);
    int     esl_abc_CDealign(const ESL_ALPHABET *abc, char    *s, const ESL_DSQ *ref_ax, int64_t *opt_rlen);
    int     esl_abc_XDealign(const ESL_ALPHABET *abc, ESL_DSQ *x, const ESL_DSQ *ref_ax, int64_t *opt_rlen);
    int     esl_abc_ConvertDegen2X(const ESL_ALPHABET *abc, ESL_DSQ *dsq);
    int     esl_abc_revcomp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int n);

    # Other routines in the API
    int    esl_abc_ValidateType(int type);
    int    esl_abc_GuessAlphabet(const int64_t *ct, int *ret_type);
    double esl_abc_Match       (const ESL_ALPHABET *a, ESL_DSQ x, ESL_DSQ y, double *p);
    int    esl_abc_IAvgScore   (const ESL_ALPHABET *a, ESL_DSQ x, const int    *sc);
    float  esl_abc_FAvgScore   (const ESL_ALPHABET *a, ESL_DSQ x, const float  *sc);
    double esl_abc_DAvgScore   (const ESL_ALPHABET *a, ESL_DSQ x, const double *sc);
    int    esl_abc_IExpectScore(const ESL_ALPHABET *a, ESL_DSQ x, const int    *sc, const float  *p);
    float  esl_abc_FExpectScore(const ESL_ALPHABET *a, ESL_DSQ x, const float  *sc, const float  *p);
    double esl_abc_DExpectScore(const ESL_ALPHABET *a, ESL_DSQ x, const double *sc, const double *p);

    int    esl_abc_IAvgScVec   (const ESL_ALPHABET *a, int    *sc);
    int    esl_abc_FAvgScVec   (const ESL_ALPHABET *a, float  *sc);
    int    esl_abc_DAvgScVec   (const ESL_ALPHABET *a, double *sc);
    int    esl_abc_IExpectScVec(const ESL_ALPHABET *a, int    *sc, const float  *p);
    int    esl_abc_FExpectScVec(const ESL_ALPHABET *a, float  *sc, const float  *p);
    int    esl_abc_DExpectScVec(const ESL_ALPHABET *a, double *sc, const double *p);
    int    esl_abc_FCount      (const ESL_ALPHABET *a, float  *ct, ESL_DSQ x, float  wt);
    int    esl_abc_DCount      (const ESL_ALPHABET *a, double *ct, ESL_DSQ x, double wt);
    int    esl_abc_EncodeType   (char *typestring);
    int    esl_abc_EncodeTypeMem(char *type, int n);
    char  *esl_abc_DecodeType   (int type);
    int    esl_abc_ValidateSeq(const ESL_ALPHABET *a, const char *seq, int64_t L, char *errbuf);

    ESL_DSQ esl_abc_DigitizeSymbol (ESL_ALPHABET*, char)

    bint    esl_abc_XIsValid       (ESL_ALPHABET*, int)
    bint    esl_abc_XIsResidue     (ESL_ALPHABET*, int)
    bint    esl_abc_XIsCanonical   (ESL_ALPHABET*, int)
    bint    esl_abc_XIsGap         (ESL_ALPHABET*, int)
    bint    esl_abc_XIsDegenerate  (ESL_ALPHABET*, int)
    bint    esl_abc_XIsUnknown     (ESL_ALPHABET*, int)
    bint    esl_abc_XIsNonresidue  (ESL_ALPHABET*, int)
    bint    esl_abc_XIsMissing     (ESL_ALPHABET*, int)
    int     esl_abc_XGetGap        (ESL_ALPHABET*)
    int     esl_abc_XGetUnknown    (ESL_ALPHABET*)
    int     esl_abc_XGetNonresidue (ESL_ALPHABET*)
    int     esl_abc_XGetMissing    (ESL_ALPHABET*)

    bint    esl_abc_CIsValid       (ESL_ALPHABET*, char)
    bint    esl_abc_CIsResidue     (ESL_ALPHABET*, char)
    bint    esl_abc_CIsCanonical   (ESL_ALPHABET*, char)
    bint    esl_abc_CIsGap         (ESL_ALPHABET*, char)
    bint    esl_abc_CIsDegenerate  (ESL_ALPHABET*, char)
    bint    esl_abc_CIsUnknown     (ESL_ALPHABET*, char)
    bint    esl_abc_CIsNonresidue  (ESL_ALPHABET*, char)
    bint    esl_abc_CIsMissing     (ESL_ALPHABET*, char)
    char    esl_abc_CGetGap        (ESL_ALPHABET*)
    char    esl_abc_CGetUnknown    (ESL_ALPHABET*)
    char    esl_abc_CGetNonresidue (ESL_ALPHABET*)
    char    esl_abc_CGetMissing    (ESL_ALPHABET*)
