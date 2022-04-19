from libc.stdint cimport uint8_t, uint32_t, uint64_t, int64_t

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.getopts import ESL_GETOPTS
from libeasel.sq cimport ESL_SQ
from libeasel.rand64 cimport ESL_RAND64

from libhmmer.p7_pipeline cimport p7_zsetby_e
from libhmmer.p7_hmm cimport P7_HMM


cdef extern from "hmmpgmd.h" nogil:

    ctypedef struct HMMD_SEARCH_STATUS:
        uint32_t status
        uint64_t msg_size

    ctypedef struct HMMD_SEARCH_STATS:
        double   elapsed
        double   user
        double   sys

        double   Z
        double   domZ

        p7_zsetby_e Z_setby
        p7_zsetby_e domZ_setby

        uint64_t    nmodels
        uint64_t    nseqs
        uint64_t    n_past_msv
        uint64_t    n_past_bias
        uint64_t    n_past_vit
        uint64_t    n_past_fwd

        uint64_t    nhits
        uint64_t    nreported
        uint64_t    nincluded
        uint64_t*   hit_offsets


    const uint32_t HMMD_SEQUENCE
    const uint32_t HMMD_HMM

    const uint32_t HMMD_CMD_SEARCH
    const uint32_t HMMD_CMD_SCAN
    const uint32_t HMMD_CMD_INIT
    const uint32_t HMMD_CMD_SHUTDOWN

    const uint32_t MAX_INIT_DESC


    ctypedef struct HMMD_SEARCH_CMD:
        uint32_t  db_indx
        uint32_t  db_type
        uint32_t  inx
        uint32_t  cnt
        uint32_t  query_type
        uint32_t  query_length
        uint32_t  opts_length
        char      data[1]

    ctypedef struct HMMD_INIT_CMD:
        char      sid[MAX_INIT_DESC]
        char      hid[MAX_INIT_DESC]
        uint32_t  seqdb_off
        uint32_t  hmmdb_off
        uint32_t  db_cnt
        uint32_t  seq_cnt
        uint32_t  hmm_cnt
        uint32_t  model_cnt
        char      data[1]

    ctypedef struct HMMD_INIT_RESET:
        char      ip_addr[1]

    ctypedef struct HMMD_HEADER:
        uint32_t  length
        uint32_t  command
        uint32_t  status

    ctypedef struct HMMD_COMMAND:
        HMMD_HEADER   hdr
        HMMD_INIT_CMD init
        HMMD_SEARCH_CMD srch
        HMMD_INIT_RESET reset

    const size_t HMMD_SEARCH_STATUS_SERIAL_SIZE
    const size_t HMMD_SEARCH_STATS_SERIAL_BASE

    cdef struct queue_data_s:
        uint32_t  cmd_type
        uint32_t  query_type
        P7_HMM*   hmm
        ESL_SQ*   seq
        ESL_ALPHABET* abc
        # ESL_GETOPTS*  opts
        HMMD_COMMAND* cmd

        int  sock
        char ip_addr[64]

        int  dbx
        int  inx
        int  cnt
    ctypedef queue_data_s QUEUE_DATA

    ctypedef struct RANGE_LIST:
        int N
        uint32_t* starts
        uint32_t* ends

    cdef void free_QueueData(QUEUE_DATA* data)
    cdef int  hmmpgmd_IsWithinRanges(int64_t sq_idx, RANGE_LIST* list)
    cdef int  hmmpgmd_GetRanges(RANGE_LIST* list, char* rangestr)

    # cdef int process_searchopts(int fd, char* cmdstr, ESL_GETOPTS** ret_opts)

    # cdef void worker_process(ESL_GETOPTS* go)
    # cdef void master_process(ESL_GETOPTS* go)

    cdef int p7_hmmd_search_stats_Serialize(const HMMD_SEARCH_STATS *obj, uint8_t **buf, uint32_t *n, uint32_t *nalloc)
    cdef int p7_hmmd_search_stats_Deserialize(const uint8_t *buf, uint32_t *pos, HMMD_SEARCH_STATS *ret_obj)

    cdef int hmmd_search_status_Serialize(const HMMD_SEARCH_STATUS *obj, uint8_t **buf, uint32_t *n, uint32_t *nalloc)
    cdef int hmmd_search_status_Deserialize(const uint8_t *buf, uint32_t *n, HMMD_SEARCH_STATUS *ret_obj)
    cdef int hmmd_search_status_TestSample(ESL_RAND64 *rng, HMMD_SEARCH_STATUS **ret_obj)
    cdef int hmmd_search_status_Compare(HMMD_SEARCH_STATUS *first, HMMD_SEARCH_STATUS *second)
