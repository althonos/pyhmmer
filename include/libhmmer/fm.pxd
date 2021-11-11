from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int64_t


cdef extern from "hmmer.h" nogil:

    cdef struct fm_data_s:
        uint64_t N
        uint32_t term_loc
        uint32_t seq_offset
        uint32_t ambig_offset
        uint32_t seq_cnt
        uint32_t ambig_cnt
        uint32_t overlap
        uint8_t* T
        uint8_t* BWT_mem
        uint8_t* BWT
        uint32_t* SA
        int64_t*  C
        uint32_t* occCnts_sb
        uint16_t* occCnts_b
    ctypedef fm_data_s FM_DATA

    ctypedef struct FM_CFG:
        pass
