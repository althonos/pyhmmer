from libc.stdint cimport int64_t, uint32_t
from posix.types cimport off_t

from libeasel cimport ESL_DSQ, esl_pos_t
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.bitfield cimport ESL_BITFIELD
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.random cimport ESL_RANDOMNESS


cdef extern from "esl_msa.h" nogil:

    cdef enum:
        eslMSA_TC1 = 0
        eslMSA_TC2 = 1
        eslMSA_GA1 = 2
        eslMSA_GA2 = 3
        eslMSA_NC1 = 4
        eslMSA_NC2 = 5
        eslMSA_NCUTS = 6

    cdef enum:
        eslMSA_HASWGTS = 0b01
        eslMSA_DIGITAL = 0b10

    ctypedef struct ESL_MSA:
        char** aseq
        char** sqname
        double* wgt
        int64_t alen
        int nseq
        int flags

        ESL_ALPHABET* abc
        ESL_DSQ** ax

        char* name
        char* desc
        char* acc
        char* au
        char* ss_cons
        char* sa_cons
        char* pp_cons
        char* rf
        char* mm
        char** sqacc
        char** sqdesc
        char** ss
        char** sa
        char** pp
        float[eslMSA_NCUTS] cutoff
        float[eslMSA_NCUTS] cutset

        int sqalloc
        int64_t* sqlen
        int64_t* sslen
        int64_t* salen
        int64_t* pplen
        int lastidx

        char** comment
        int ncomment
        int alloc_ncomment

        char** gf_tag
        char** gf
        int ngf
        int alloc_ngf

        char*** gs_tag
        char*** gs
        int ngs

        char** gr_tag
        char*** gr
        int ngr

        ESL_KEYHASH* index
        ESL_KEYHASH* gs_idx
        ESL_KEYHASH* gc_idx
        ESL_KEYHASH* gr_idx

        off_t offset


    # The ESL_MSA object
    ESL_MSA *esl_msa_Create(int nseq, int64_t alen)
    int      esl_msa_Expand(ESL_MSA *msa)
    int      esl_msa_Copy(const ESL_MSA *msa, ESL_MSA *new_msa)
    ESL_MSA *esl_msa_Clone(const ESL_MSA *msa)
    void     esl_msa_Destroy(ESL_MSA *msa)

    # Digital mode MSA's
    int      esl_msa_GuessAlphabet(const ESL_MSA *msa, int *ret_type)
    ESL_MSA *esl_msa_CreateDigital(const ESL_ALPHABET *abc, int nseq, int64_t alen)
    int      esl_msa_Digitize(const ESL_ALPHABET *abc, ESL_MSA *msa, char *errmsg)
    int      esl_msa_Textize(ESL_MSA *msa)
    int      esl_msa_ConvertDegen2X(ESL_MSA *msa)

    # Setting or checking data fields in an ESL_MSA
    int esl_msa_SetName          (ESL_MSA *msa, const char *s, esl_pos_t n)
    int esl_msa_SetDesc          (ESL_MSA *msa, const char *s, esl_pos_t n)
    int esl_msa_SetAccession     (ESL_MSA *msa, const char *s, esl_pos_t n)
    int esl_msa_SetAuthor        (ESL_MSA *msa, const char *s, esl_pos_t n)
    int esl_msa_SetSeqName       (ESL_MSA *msa, int idx, const char *s, esl_pos_t n)
    int esl_msa_SetSeqAccession  (ESL_MSA *msa, int idx, const char *s, esl_pos_t n)
    int esl_msa_SetSeqDescription(ESL_MSA *msa, int idx, const char *s, esl_pos_t n)
    int esl_msa_SetDefaultWeights(ESL_MSA *msa)

    int esl_msa_FormatName          (ESL_MSA *msa, const char *name,    ...)
    int esl_msa_FormatDesc          (ESL_MSA *msa, const char *desc,    ...)
    int esl_msa_FormatAccession     (ESL_MSA *msa, const char *acc,     ...)
    int esl_msa_FormatAuthor        (ESL_MSA *msa, const char *author,  ...)
    int esl_msa_FormatSeqName       (ESL_MSA *msa, int idx, const char *name, ...)
    int esl_msa_FormatSeqAccession  (ESL_MSA *msa, int idx, const char *acc, ...)
    int esl_msa_FormatSeqDescription(ESL_MSA *msa, int idx, const char *desc, ...)

    int esl_msa_AddComment(ESL_MSA *msa, char *p,   esl_pos_t n)
    int esl_msa_AddGF     (ESL_MSA *msa, char *tag, esl_pos_t taglen,            char *value, esl_pos_t vlen)
    int esl_msa_AddGS     (ESL_MSA *msa, char *tag, esl_pos_t taglen, int sqidx, char *value, esl_pos_t vlen)
    int esl_msa_AppendGC  (ESL_MSA *msa, char *tag, char *value)
    int esl_msa_AppendGR  (ESL_MSA *msa, char *tag, int sqidx, char *value)

    int esl_msa_CheckUniqueNames(const ESL_MSA *msa)

    # Miscellaneous functions for manipulating MSAs
    int esl_msa_ReasonableRF(ESL_MSA *msa, double symfrac, int useconsseq, char *rfline)
    int esl_msa_MarkFragments (const ESL_MSA *msa, float fragthresh, ESL_BITFIELD **ret_fragassign)
    int esl_msa_MarkFragments_old(ESL_MSA *msa, double fragthresh)
    int esl_msa_SequenceSubset(const ESL_MSA *msa, const int *useme, ESL_MSA **ret_new)
    int esl_msa_ColumnSubset (ESL_MSA *msa, char *errbuf, const int *useme)
    int esl_msa_MinimGaps    (ESL_MSA *msa, char *errbuf, const char *gaps, int consider_rf)
    int esl_msa_MinimGapsText(ESL_MSA *msa, char *errbuf, const char *gaps, int consider_rf, int fix_bps)
    int esl_msa_NoGaps       (ESL_MSA *msa, char *errbuf, const char *gaps)
    int esl_msa_NoGapsText   (ESL_MSA *msa, char *errbuf, const char *gaps, int fix_bps)
    int esl_msa_SymConvert(ESL_MSA *msa, const char *oldsyms, const char *newsyms)
    int esl_msa_Checksum(const ESL_MSA *msa, uint32_t *ret_checksum)

    int esl_msa_RemoveBrokenBasepairsFromSS(char *ss, char *errbuf, int len, const int *useme)
    int esl_msa_RemoveBrokenBasepairs(ESL_MSA *msa, char *errbuf, const int *useme)

    int esl_msa_ReverseComplement(ESL_MSA *msa)
    int esl_msa_Hash(ESL_MSA *msa)
    int esl_msa_FlushLeftInserts(ESL_MSA *msa)

    # Debugging, testing, development
    int      esl_msa_Validate(const ESL_MSA *msa, char *errmsg)
    ESL_MSA *esl_msa_CreateFromString(const char *s, int fmt)
    int      esl_msa_Compare         (ESL_MSA *a1, ESL_MSA *a2)
    int      esl_msa_CompareMandatory(ESL_MSA *a1, ESL_MSA *a2)
    int      esl_msa_CompareOptional (ESL_MSA *a1, ESL_MSA *a2)
    int      esl_msa_Sample(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int max_nseq, int max_alen, ESL_MSA **ret_msa)
