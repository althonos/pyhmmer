from libc.stdint cimport int64_t
from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel cimport ESL_DSQ
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.sqio cimport ESL_SQFILE
from libeasel.sq cimport ESL_SQ, ESL_SQ_BLOCK


cdef extern from "reexports/esl_sqio_ascii.h":

    int   sqascii_GuessFileFormat(ESL_SQFILE *sqfp, int *ret_fmt)
    int   sqascii_Position       (ESL_SQFILE *sqfp, off_t offset)
    void  sqascii_Close          (ESL_SQFILE *sqfp)
    int   sqascii_SetDigital     (ESL_SQFILE *sqfp, const ESL_ALPHABET *abc)
    int   sqascii_GuessAlphabet  (ESL_SQFILE *sqfp, int *ret_type)
    int   sqascii_Read           (ESL_SQFILE *sqfp, ESL_SQ *sq)
    int   sqascii_ReadInfo       (ESL_SQFILE *sqfp, ESL_SQ *sq)
    int   sqascii_ReadSequence   (ESL_SQFILE *sqfp, ESL_SQ *sq)
    int   sqascii_ReadWindow     (ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq)
    int   sqascii_ReadBlock      (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int max_init_window, int long_target)
    int   sqascii_Echo           (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp)

    int   sqascii_IsRewindable   (const ESL_SQFILE *sqfp)
    const char *sqascii_GetError (const ESL_SQFILE *sqfp)

    int   sqascii_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint)
    int   sqascii_PositionByKey   (ESL_SQFILE *sqfp, const char *key)
    int   sqascii_PositionByNumber(ESL_SQFILE *sqfp, int which)
    int   sqascii_Fetch           (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq)
    int   sqascii_FetchInfo       (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq)
    int   sqascii_FetchSubseq     (ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq)

    # Internal routines shared by parsers.
    int  loadmem  (ESL_SQFILE *sqfp) except -1
    int  loadbuf  (ESL_SQFILE *sqfp) except -1
    int  nextchar (ESL_SQFILE *sqfp, char *ret_c) except -1
    int  seebuf   (ESL_SQFILE *sqfp, int64_t maxn, int64_t *opt_nres, int64_t *opt_endpos) except -1
    void addbuf   (ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nres) except *
    void skipbuf  (ESL_SQFILE *sqfp, int64_t nskip) except *
    int  read_nres(ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nskip, int64_t nres, int64_t *opt_actual_nres) except -1
    int  skip_whitespace(ESL_SQFILE *sqfp) except -1

    # EMBL format also UniProt, TrEMBL
    void config_embl(ESL_SQFILE *sqfp) nogil
    void inmap_embl (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap) nogil
    int  header_embl(ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  skip_embl  (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  end_embl   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # GenBank format also DDBJ
    void config_genbank(ESL_SQFILE *sqfp) nogil
    void inmap_genbank (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap) nogil
    int  header_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  skip_genbank  (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  end_genbank   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # FASTA format
    void config_fasta(ESL_SQFILE *sqfp) nogil
    void inmap_fasta (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap) nogil
    int  header_fasta(ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  skip_fasta  (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  end_fasta   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # daemon format
    void config_daemon(ESL_SQFILE *sqfp) nogil
    void inmap_daemon (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap) nogil
    int  end_daemon   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # HMMPGMD format
    int  fileheader_hmmpgmd(ESL_SQFILE *sqfp) except -1
