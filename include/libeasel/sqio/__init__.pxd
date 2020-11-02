from libc.stdint cimport int64_t, uint32_t, uint64_t
from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel cimport ESL_DSQ
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.sq cimport ESL_SQ, ESL_SQ_BLOCK
from libeasel.sqio.ascii cimport ESL_SQASCII_DATA
from libeasel.sqio.ncbi cimport ESL_SQNCBI_DATA


cdef extern from "esl_sqio.h" nogil:

    ctypedef union ESL_SQDATA:
        ESL_SQASCII_DATA ascii
        ESL_SQNCBI_DATA  ncbi


    ctypedef esl_sqio_s ESL_SQFILE
    cdef struct esl_sqio_s:
        char *filename
        bint   do_digital
        const ESL_ALPHABET *abc

        int     format
        ESL_DSQ[128] inmap

        int   (*position)        (esl_sqio_s *sqfp, off_t offset)
        void  (*close)           (esl_sqio_s *sqfp)

        int   (*set_digital)     (esl_sqio_s *sqfp, const ESL_ALPHABET *abc)
        int   (*guess_alphabet)  (esl_sqio_s *sqfp, int *ret_type)

        int   (*read)            (esl_sqio_s *sqfp, ESL_SQ *sq)
        int   (*read_info)       (esl_sqio_s *sqfp, ESL_SQ *sq)
        int   (*read_seq)        (esl_sqio_s *sqfp, ESL_SQ *sq)
        int   (*read_window)     (esl_sqio_s *sqfp, int C, int W, ESL_SQ *sq)
        int   (*echo)            (esl_sqio_s *sqfp, const ESL_SQ *sq, FILE *ofp)

        int   (*read_block)      (esl_sqio_s *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int max_init_window, int long_target)

        int   (*open_ssi)        (esl_sqio_s *sqfp, const char *ssifile_hint)
        int   (*pos_by_key)      (esl_sqio_s *sqfp, const char *key)
        int   (*pos_by_number)   (esl_sqio_s *sqfp, int which)

        int   (*fetch)           (esl_sqio_s *sqfp, const char *key, ESL_SQ *sq)
        int   (*fetch_info)      (esl_sqio_s *sqfp, const char *key, ESL_SQ *sq)
        int   (*fetch_subseq)    (esl_sqio_s *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq)

        int   (*is_rewindable)   (const esl_sqio_s *sqfp)
        const char *(*get_error) (const esl_sqio_s *sqfp)

        ESL_SQDATA data


    ctypedef esl_sqcache_s ESL_SQCACHE
    cdef struct esl_sqcache_s:
        char* filename
        int format
        const ESL_ALPHABET* abc
        uint32_t seq_count
        uint64_t res_count
        uint32_t max_seq
        ESL_SQ* sq_list
        void* residue_mem
        void* header_mem
        uint64_t res_size
        uint64_t hdr_size


    cdef enum:
        eslSQFILE_UNKNOWN = 0
        eslSQFILE_FASTA = 1
        eslSQFILE_EMBL = 2
        eslSQFILE_GENBANK = 3
        eslSQFILE_DDBJ = 4
        eslSQFILE_UNIPROT = 5
        eslSQFILE_NCBI = 6
        eslSQFILE_DAEMON = 7
        eslSQFILE_HMMPGMD = 8
        eslSQFILE_FMINDEX = 9


    cdef int eslREADBUFSIZE


    int  esl_sqfile_Open(const char *seqfile, int fmt, const char *env, ESL_SQFILE **ret_sqfp)
    void esl_sqfile_Close(ESL_SQFILE *sqfp)

    int  esl_sqfile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp)
    int  esl_sqfile_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc)
    int  esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type)

    int   esl_sqio_Read        (ESL_SQFILE *sqfp, ESL_SQ *sq)
    int   esl_sqio_ReadInfo    (ESL_SQFILE *sqfp, ESL_SQ *sq)
    int   esl_sqio_ReadWindow  (ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq)
    int   esl_sqio_ReadSequence(ESL_SQFILE *sqfp, ESL_SQ *sq)
    int   esl_sqio_ReadBlock   (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int max_init_window, int long_target)
    int   esl_sqio_Parse       (char *buffer, int size, ESL_SQ *s, int format)

    int   esl_sqio_Write       (FILE *fp, ESL_SQ *s, int format, int update)
    int   esl_sqio_Echo        (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp)

    const char  *esl_sqfile_GetErrorBuf(const ESL_SQFILE *sqfp)
    int   esl_sqfile_IsRewindable(const ESL_SQFILE *sqfp)
    int   esl_sqio_IsAlignment(int fmt)
    int   esl_sqio_EncodeFormat(char *fmtstring)
    char *esl_sqio_DecodeFormat(int fmt)
    int   esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset)
    int   esl_sqio_Ignore(ESL_SQFILE *sqfp, const char *ignoredchars)
    int   esl_sqio_AcceptAs(ESL_SQFILE *sqfp, char *xchars, char readas)

    int   esl_sqfile_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint)
    int   esl_sqfile_PositionByKey   (ESL_SQFILE *sqfp, const char *key)
    int   esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which)

    int   esl_sqio_Fetch      (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq)
    int   esl_sqio_FetchInfo  (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq)
    int   esl_sqio_FetchSubseq(ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq)

    int   esl_sqfile_Cache(const ESL_ALPHABET *abc, const char *seqfile, int fmt, const char *env, ESL_SQCACHE **ret_sqcache)
    void  esl_sqfile_Free(ESL_SQCACHE *sqcache)
