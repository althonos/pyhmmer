from libc.stdint cimport uint32_t, int32_t, int64_t

cdef extern from "esl_dsqdata.h" nogil:

    const size_t eslDSQDATA_CHUNK_MAXSEQ
    const size_t eslDSQDATA_CHUNK_MAXPACKET
    const size_t eslDSQDATA_UNPACKERS
    const size_t eslDSQDATA_UMAX

    ctypedef esl_dsqdata_chunk_s ESL_DSQDATA_CHUNK
    cdef struct esl_dsqdata_chunk_s:
        int64_t i0
        int N

        ESL_DSQ** dsq
        char** name
        char** acc
        char** desc
        int32_t* taxid
        int64_t* L

        unsigned char* smem
        uint32_t *psq
        int pn
        char* metadata
        int mdalloc
        esl_dsqdata_chunk_s* nxt


    ctypedef esl_dsqdata_record_s ESL_DSQDATA_RECORD
    cdef struct esl_dsqdata_record_s:
        int64_t metadata_end
        int64_t psq_end


    ctypedef esl_dsqdata_s ESL_DSQDATA
    cdef esl_dsqdata_s:
        char* basename
        FILE* stubfp
        FILE* ifp
        FILE* sfp
        FILE* mfp
        ESL_ALPHABET* abc_r

        uint32_t magic
        uint32_t uniquetag
        uint32_t flags
        uint32_t max_namelen
        uint32_t max_acclen
        uint32_t max_desclen
        uint32_t max_seqlen
        uint32_t nseq
        uint32_t nres

        int chunk_maxseq
        int chunk_maxpacket
        int do_byteswap
        int pack5

        int nconsumers
        int n_unpackers

        ESL_DSQDATA_CHUNK*[eslDSQDATA_UMAX] inbox
        pthread_mutex_t[eslDSQDATA_UMAX] outbox_mutex
        pthread_cond_t[eslDSQDATA_UMAX] outbox_cv
        int[eslDSQDATA_UMAX] outbox_eod

        int64_t nchunk
        pthread_mutex_t nchunk_mutex

        ESL_DSQDATA_CHUNK* recycling
        pthread_mutex_t recycling_mutex
        pthread_cond_t recycling_cv

        int go
        pthread_mutex_t go_mutex
        pthread_cond_t go_cv

        pthread_t loader_t
        pthread_t unpacker_t[eslDSQDATA_UMAX]

        char[eslERRBUFSIZE] errbuf


    int  esl_dsqdata_Open   (ESL_ALPHABET **byp_abc, char *basename, int nconsumers, ESL_DSQDATA **ret_dd)
    int  esl_dsqdata_Read   (ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK **ret_chu)
    int  esl_dsqdata_Recycle(ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK *chu)
    int  esl_dsqdata_Close  (ESL_DSQDATA *dd)
    int  esl_dsqdata_Write  (ESL_SQFILE *sqfp, char *basename, char *errbuf)
