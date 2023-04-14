from libc.stdint cimport int32_t, int64_t, uint32_t
from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel cimport ESL_DSQ, eslERRBUFSIZE
from libeasel.msa cimport ESL_MSA
from libeasel.msafile cimport ESL_MSAFILE
from libeasel.sq cimport ESL_SQ
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
