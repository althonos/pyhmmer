from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t
from libc.stdio cimport FILE

from libeasel.sqio cimport esl_sqio_s


cdef extern from "easel.h":
    DEF eslERRBUFSIZE = 128


cdef extern from "esl_sqio.h" nogil:

    DEF MAX_DB_VOLUMES = 100
    DEF MAX_RESIDUE_COUNT = 1024 * 1024

    ctypedef esl_sqncbi_vol_s ESL_SQNCBI_VOLUME
    cdef struct esl_sqncbi_vol_s:
        char* name
        uint32_t start_seq
        uint32_t end_seq
        uint32_t hdr_off
        uint32_t seq_off
        uint32_t amb_off


    ctypedef esl_sqncbi_s ESL_SQNCBI_DATA
    cdef struct esl_sqncbi_s:
        FILE      *fppin
        FILE      *fpphr
        FILE      *fppsq
        char[eslERRBUFSIZE] errbuf

        char *title
        int  version
        char *timestamp

        uint32_t   num_seq
        uint64_t   total_res
        uint32_t   max_seq

        uint32_t   hdr_off
        uint32_t   seq_off
        uint32_t   amb_off

        int        index
        uint32_t   vol_index
        uint32_t   roff
        uint32_t   hoff
        uint32_t   doff
        uint32_t   eoff

        uint32_t   index_start
        uint32_t   index_end
        uint32_t  *hdr_indexes
        uint32_t  *seq_indexes
        uint32_t  *amb_indexes

        uint32_t   volumes
        ESL_SQNCBI_VOLUME[MAX_DB_VOLUMES] vols

        unsigned char *hdr_buf
        unsigned char *hdr_ptr
        int            hdr_alloced

        char          *name_ptr
        int32_t        name_size
        char          *acc_ptr
        int32_t        acc_size
        int32_t        int_id
        char          *str_id_ptr
        int32_t        str_id_size

        uint32_t       seq_apos
        uint32_t       seq_alen
        uint32_t       seq_cpos
        int32_t        seq_L

        int            alphatype
        char          *alphasym


    int  esl_sqncbi_Open(char *seqfile, int format, esl_sqio_s *sqfp)
