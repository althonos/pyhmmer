from libc.stdint cimport uint8_t, uint32_t, int64_t

from libeasel.rand64 cimport ESL_RAND64


cdef extern from "hmmer.h" nogil:

    # forward declaration of P7_ALIDISPLAY to avoid circular imports
    ctypedef p7_alidisplay_s P7_ALIDISPLAY
    cdef struct p7_alidisplay_s:
        pass

    ctypedef p7_dom_s P7_DOMAIN
    cdef struct p7_dom_s:
        int64_t ienv
        int64_t jenv
        int64_t iali
        int64_t jali
        int64_t iorf
        int64_t jorf
        float envsc
        float domcorrection
        float dombias
        float oasc
        float bitscore
        double lnP
        bint is_reported
        bint is_included
        float* scores_per_pos
        P7_ALIDISPLAY* ad

    P7_DOMAIN *p7_domain_Create_empty();
    void p7_domain_Destroy(P7_DOMAIN *obj);
    int p7_domain_Serialize(const P7_DOMAIN *obj, uint8_t **buf, uint32_t *n, uint32_t *nalloc);
    int p7_domain_Deserialize(const uint8_t *buf, uint32_t *n, P7_DOMAIN *ret_obj);
    int p7_domain_TestSample(ESL_RAND64 *rng, P7_DOMAIN **ret_obj);
    int p7_domain_Compare(P7_DOMAIN *first, P7_DOMAIN *second, double atol, double rtol);
