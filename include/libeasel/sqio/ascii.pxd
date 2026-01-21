from libc.stdint cimport int32_t, int64_t, uint32_t
from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel cimport ESL_DSQ, eslERRBUFSIZE
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.msa cimport ESL_MSA
from libeasel.msafile cimport ESL_MSAFILE
from libeasel.sq cimport ESL_SQ, ESL_SQ_BLOCK
from libeasel.sqio cimport esl_sqio_s
from libeasel.ssi cimport ESL_SSI


cdef extern from "esl_sqio.h" nogil:
    const size_t MAX_RESIDUE_COUNT

    ctypedef esl_sqascii_s ESL_SQASCII_DATA
    cdef struct esl_sqascii_s:
        FILE *fp
        char[eslERRBUFSIZE]  errbuf

        bint   do_gzip
        bint   do_stdin
        bint   do_buffer

        char    *mem
        int      allocm
        int      mn
        int      mpos
        off_t    moff
        bint      is_recording

        char    *buf
        off_t    boff
        int      balloc
        int      nc
        int      bpos
        int64_t  L
        int64_t  linenumber
        off_t    bookmark_offset
        int64_t  bookmark_linenum

        bint   is_linebased
        bint   eof_is_ok
        int  (*parse_header)(esl_sqio_s *, ESL_SQ *sq)
        int  (*skip_header) (esl_sqio_s *, ESL_SQ *sq)
        int  (*parse_end)   (esl_sqio_s *, ESL_SQ *sq)

        ESL_MSAFILE  *afp
        ESL_MSA      *msa
        int           idx

        char    *ssifile
        int      rpl
        int      bpl
        int      currpl
        int      curbpl
        int      prvrpl
        int      prvbpl
        ESL_SSI *ssi

    int  esl_sqascii_Open(char *seqfile, int format, esl_sqio_s *sqfp)
    int  esl_sqascii_WriteFasta(FILE *fp, ESL_SQ *s, int update)
    int  esl_sqascii_Parse(char *buf, int size, ESL_SQ *s, int format)

    ctypedef esl_sqio_s ESL_SQFILE

    # Private re-exports
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
    void config_embl(ESL_SQFILE *sqfp)
    void inmap_embl (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
    int  header_embl(ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  skip_embl  (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  end_embl   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # GenBank format also DDBJ
    void config_genbank(ESL_SQFILE *sqfp)
    void inmap_genbank (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
    int  header_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  skip_genbank  (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  end_genbank   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # FASTA format
    void config_fasta(ESL_SQFILE *sqfp)
    void inmap_fasta (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
    int  header_fasta(ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  skip_fasta  (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1
    int  end_fasta   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # daemon format
    void config_daemon(ESL_SQFILE *sqfp)
    void inmap_daemon (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
    int  end_daemon   (ESL_SQFILE *sqfp, ESL_SQ *sq) except -1

    # HMMPGMD format
    int  fileheader_hmmpgmd(ESL_SQFILE *sqfp) except -1