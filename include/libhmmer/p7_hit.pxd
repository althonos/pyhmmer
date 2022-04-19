from libc.stdint cimport int64_t, uint8_t, uint32_t, uint64_t
from libc.stdio cimport FILE

from libeasel cimport esl_pos_t
from libeasel.rand64 cimport ESL_RAND64
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.getopts cimport ESL_GETOPTS
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.msa cimport ESL_MSA
from libeasel.sq cimport ESL_SQ
from libhmmer.p7_domain cimport P7_DOMAIN
from libhmmer.p7_pipeline cimport p7_pipemodes_e, P7_PIPELINE
from libhmmer.p7_trace cimport P7_TRACE
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY


cdef extern from "hmmer.h" nogil:

    cdef enum p7_hitflags_e:
        p7_HITFLAGS_DEFAULT = 0b00000
        p7_IS_INCLUDED      = 0b00001
        p7_IS_REPORTED      = 0b00010
        p7_IS_NEW           = 0b00100
        p7_IS_DROPPED       = 0b01000
        p7_IS_DUPLICATE     = 0b10000

    ctypedef struct P7_HIT:
        char* name
        char* acc
        char* desc
        int window_length
        double sortkey

        float score
        float pre_score
        float sum_score

        double lnP
        double pre_lnP
        double sum_lnP

        float nexpected
        int nregions
        int nclustered
        int noverlaps
        int nenvelopes
        int ndom

        uint32_t flags
        int nreported
        int nincluded
        int best_domain

        int64_t seqidx
        int64_t subseq_start

        P7_DOMAIN *dcl
        esl_pos_t offset

    cdef P7_HIT *p7_hit_Create_empty()
    cdef void p7_hit_Destroy(P7_HIT *the_hit)
    cdef int p7_hit_Serialize(const P7_HIT *obj, uint8_t **buf, uint32_t *n, uint32_t *nalloc)
    cdef int p7_hit_Deserialize(const uint8_t *buf, uint32_t *n, P7_HIT *ret_obj)
    cdef int p7_hit_TestSample(ESL_RAND64 *rng, P7_HIT **ret_obj)
    cdef int p7_hit_Compare(P7_HIT *first, P7_HIT *second, double atol, double rtol)
