from libc.stdio cimport FILE
from libc.stdint cimport int16_t, int32_t, int64_t, uint16_t, uint32_t, uint64_t
from posix.types cimport off_t

from libeasel cimport eslERRBUFSIZE


cdef extern from "esl_ssi.h" nogil:

    const size_t eslSSI_MAXFILES
    const size_t eslSSI_MAXKEYS
    const size_t eslSSI_MAXRAM


    ctypedef struct ESL_SSI:
        FILE* fp
        uint32_t flags
        uint32_t offsz
        uint16_t nfiles
        uint64_t nprimary
        uint64_t nsecondary
        uint32_t flen
        uint32_t plen
        uint32_t slen
        uint32_t frecsize
        uint32_t precsize
        uint32_t srecsize
        off_t foffset
        off_t poffset
        off_t soffset
        char** filename
        uint32_t* fileformat
        uint32_t* fileflags
        uint32_t* bpl
        uint32_t* rpl


    ctypedef struct ESL_PKEY:
        char* key
        uint16_t fnum
        off_t r_off
        off_t d_off
        int64_t len


    ctypedef struct ESL_SKEY:
        char* key
        char* pkey


    ctypedef struct ESL_NEWSSI:
        char* ssifile
        FILE* ssifp
        int external
        int max_ram

        char** filename
        uint32_t* fileformat
        uint32_t* bpl
        uint32_t* rpl
        uint32_t flen
        uint16_t nfiles

        ESL_PKEY* pkeys
        uint32_t plen
        uint64_t nprimary
        char* ptmpfile
        FILE* ptmp

        ESL_SKEY* skeys
        uint32_t slen
        uint64_t nsecondary
        char* stmpfile
        FILE* stmp

        char errbuf[eslERRBUFSIZE]


    # Using (reading) SSI indices
    int  esl_ssi_Open(const char *filename, ESL_SSI **ret_ssi)
    void esl_ssi_Close(ESL_SSI *ssi)
    int  esl_ssi_FindName(ESL_SSI *ssi, const char *key,
    			     uint16_t *ret_fh, off_t *ret_roff, off_t *opt_doff, int64_t *opt_L)
    int  esl_ssi_FindNumber(ESL_SSI *ssi, int64_t nkey,
    			       uint16_t *opt_fh, off_t *opt_roff, off_t *opt_doff, int64_t *opt_L, char **opt_pkey)
    int  esl_ssi_FindSubseq(ESL_SSI *ssi, const char *key, int64_t requested_start,
    			       uint16_t *ret_fh, off_t *ret_roff, off_t *ret_doff, int64_t *ret_L, int64_t *ret_actual_start)
    int  esl_ssi_FileInfo(ESL_SSI *ssi, uint16_t fh, char **ret_filename, int *ret_format)



    # Creating (writing) SSI indices.
    int  esl_newssi_Open(const char *ssifile, int allow_overwrite, ESL_NEWSSI **ret_newssi)
    int  esl_newssi_AddFile  (ESL_NEWSSI *ns, const char *filename, int fmt, uint16_t *ret_fh)
    int  esl_newssi_SetSubseq(ESL_NEWSSI *ns, uint16_t fh, uint32_t bpl, uint32_t rpl)
    int  esl_newssi_AddKey   (ESL_NEWSSI *ns, const char *key, uint16_t fh, off_t r_off, off_t d_off, int64_t L)
    int  esl_newssi_AddAlias (ESL_NEWSSI *ns, const char *alias, const char *key)
    int  esl_newssi_Write    (ESL_NEWSSI *ns)
    void esl_newssi_Close    (ESL_NEWSSI *ns)


    # Portable binary i/o.
    void     esl_byteswap(char *swap, int nbytes)
    uint16_t esl_ntoh16(uint16_t netshort)
    uint32_t esl_ntoh32(uint32_t netlong)
    uint64_t esl_ntoh64(uint64_t net_int64)
    uint16_t esl_hton16(uint16_t hostshort)
    uint32_t esl_hton32(uint32_t hostlong)
    uint64_t esl_hton64(uint64_t host_int64)
    int      esl_fread_u16(FILE *fp, uint16_t *ret_result)
    int      esl_fread_u32(FILE *fp, uint32_t *ret_result)
    int      esl_fread_u64(FILE *fp, uint64_t *ret_result)
    int      esl_fread_i16(FILE *fp, int16_t  *ret_result)
    int      esl_fread_i32(FILE *fp, int32_t  *ret_result)
    int      esl_fread_i64(FILE *fp, int64_t  *ret_result)
    int      esl_fwrite_u16(FILE *fp, uint16_t n)
    int      esl_fwrite_u32(FILE *fp, uint32_t n)
    int      esl_fwrite_u64(FILE *fp, uint64_t n)
    int      esl_fwrite_i16(FILE *fp, int16_t  n)
    int      esl_fwrite_i32(FILE *fp, int32_t  n)
    int      esl_fwrite_i64(FILE *fp, int64_t  n)
    int	esl_fread_offset(FILE *fp, int mode, off_t *ret_offset)
    int      esl_fwrite_offset(FILE *fp, off_t offset)
