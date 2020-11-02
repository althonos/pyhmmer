from libc.stdint cimport int32_t, int64_t, uint32_t
from posix.types cimport off_t

from libeasel cimport ESL_DSQ
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.msa cimport ESL_MSA
from libeasel.random cimport ESL_RANDOMNESS


cdef extern from "esl_sq.h" nogil:

    ctypedef struct ESL_SQ:
        char* name
        char* acc
        char* desc
        int32_t tax_id
        char* seq
        ESL_DSQ* dsq
        char* ss
        int64_t n

        int64_t start
        int64_t end
        int64_t C
        int64_t W
        int64_t L

        char* source
        int nalloc
        int aalloc
        int dalloc
        int64_t salloc
        int srcalloc

        int64_t idx
        off_t roff
        off_t hoff
        off_t doff
        off_t eoff

        char** xr_tag
        char** xr
        int nxr

        ESL_ALPHABET* abc


    ctypedef struct ESL_SQ_BLOCK:
        int count
        int listSize
        int complete
        int64_t first_seqidx
        ESL_SQ* list


    ESL_SQ *esl_sq_Create()
    ESL_SQ *esl_sq_CreateFrom(const char *name, const char *seq,
    				 const char *desc, const char *acc, const char *ss)
    int     esl_sq_Grow  (ESL_SQ *sq, int64_t *ret_nsafe)
    int     esl_sq_GrowTo(ESL_SQ *sq, int64_t  n)
    int     esl_sq_Copy(const ESL_SQ *src, ESL_SQ *dst)
    int     esl_sq_Compare  (ESL_SQ *sq1, ESL_SQ *sq2)
    int     esl_sq_Reuse    (ESL_SQ *sq)
    int     esl_sq_IsDigital(const ESL_SQ *sq)
    int     esl_sq_IsText   (const ESL_SQ *sq)
    void    esl_sq_Destroy  (ESL_SQ *sq)

    int     esl_sq_SetName        (ESL_SQ *sq, const char *name)
    int     esl_sq_SetAccession   (ESL_SQ *sq, const char *acc)
    int     esl_sq_SetDesc        (ESL_SQ *sq, const char *desc)
    int     esl_sq_SetSource      (ESL_SQ *sq, const char *source)
    int     esl_sq_FormatName     (ESL_SQ *sq, const char *name,   ...) #ESL_ATTRIBUTE_FORMAT(printf, 2, 3)
    int     esl_sq_FormatAccession(ESL_SQ *sq, const char *acc,    ...) #ESL_ATTRIBUTE_FORMAT(printf, 2, 3)
    int     esl_sq_FormatDesc     (ESL_SQ *sq, const char *desc,   ...) #ESL_ATTRIBUTE_FORMAT(printf, 2, 3)
    int     esl_sq_FormatSource   (ESL_SQ *sq, const char *source, ...) #ESL_ATTRIBUTE_FORMAT(printf, 2, 3)
    int     esl_sq_AppendDesc     (ESL_SQ *sq, const char *desc)
    int     esl_sq_SetCoordComplete(ESL_SQ *sq, int64_t L)
    int     esl_sq_CAddResidue (ESL_SQ *sq, char c)
    int     esl_sq_ReverseComplement(ESL_SQ *sq)
    int     esl_sq_Checksum(const ESL_SQ *sq, uint32_t *ret_checksum)
    int     esl_sq_CountResidues(const ESL_SQ *sq, int start, int L, float *f)

    ESL_SQ *esl_sq_CreateDigital(const ESL_ALPHABET *abc)
    ESL_SQ *esl_sq_CreateDigitalFrom(const ESL_ALPHABET *abc, const char *name, const ESL_DSQ *dsq,
    					int64_t L, const char *desc, const char *acc,  const char *ss);
    int     esl_sq_Digitize(const ESL_ALPHABET *abc, ESL_SQ *sq)
    int     esl_sq_Textize(ESL_SQ *sq)
    int     esl_sq_GuessAlphabet(ESL_SQ *sq, int *ret_type)
    int     esl_sq_XAddResidue(ESL_SQ *sq, ESL_DSQ x)
    int     esl_sq_ConvertDegen2X(ESL_SQ *sq)

    int     esl_sq_GetFromMSA  (const ESL_MSA *msa, int which, ESL_SQ *sq)
    int     esl_sq_FetchFromMSA(const ESL_MSA *msa, int which, ESL_SQ **ret_sq)

    ESL_SQ_BLOCK *esl_sq_CreateBlock(int count)
    int esl_sq_BlockGrowTo(ESL_SQ_BLOCK *sqblock, int newsize, int do_digital, const ESL_ALPHABET *abc)
    ESL_SQ_BLOCK *esl_sq_CreateDigitalBlock(int count, const ESL_ALPHABET *abc)
    void          esl_sq_DestroyBlock(ESL_SQ_BLOCK *sqBlock)
    int esl_sq_BlockReallocSequences(ESL_SQ_BLOCK *block)
    int esl_sq_Sample(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int maxL, ESL_SQ **ret_sq)
