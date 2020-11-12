from libc.stdint cimport int64_t, uint8_t, uint32_t
from libc.stdio cimport FILE

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.random cimport ESL_RANDOMNESS
from libeasel.sq cimport ESL_SQ
from libhmmer.p7_trace cimport P7_TRACE
from libhmmer.p7_pipeline cimport P7_PIPELINE

IF HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
ELIF HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE


cdef extern from "hmmer.h" nogil:

    ctypedef p7_alidisplay_s P7_ALIDISPLAY
    cdef struct p7_alidisplay_s:
        char* rfline
        char* mmline
        char* csline
        char* model
        char* mline
        char* aseq
        char* ntseq
        char* ppline
        int N

        char* hmmname
        char* hmmacc
        char* hmmdesc
        int hmmfrom
        int hmmto
        int M

        char* sqname
        char* sqacc
        char* sqdesc
        int64_t sqfrom
        int64_t sqto
        int64_t L

        int memsize
        char* mem

    P7_ALIDISPLAY *p7_alidisplay_Create(const P7_TRACE *tr, int which, const P7_OPROFILE *om, const ESL_SQ *sq, const ESL_SQ *ntsq);
    P7_ALIDISPLAY *p7_alidisplay_Create_empty();
    P7_ALIDISPLAY *p7_alidisplay_Clone(const P7_ALIDISPLAY *ad);
    size_t         p7_alidisplay_Sizeof(const P7_ALIDISPLAY *ad);
    int            p7_alidisplay_Serialize(const P7_ALIDISPLAY *obj, uint8_t **buf, uint32_t *n, uint32_t *nalloc);
    int            p7_alidisplay_Deserialize(const uint8_t *buf, uint32_t *n, P7_ALIDISPLAY *ret_obj);
    int            p7_alidisplay_Serialize_old(P7_ALIDISPLAY *ad);
    int            p7_alidisplay_Deserialize_old(P7_ALIDISPLAY *ad);
    void           p7_alidisplay_Destroy(P7_ALIDISPLAY *ad);
    char           p7_alidisplay_EncodePostProb(float p);
    float          p7_alidisplay_DecodePostProb(char pc);
    char           p7_alidisplay_EncodeAliPostProb(float p, float hi, float med, float lo);

    int            p7_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, P7_PIPELINE *pli);
    int            p7_translated_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, P7_PIPELINE *pli);
    int            p7_nontranslated_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, int show_accessions);

    int            p7_alidisplay_Backconvert(const P7_ALIDISPLAY *ad, const ESL_ALPHABET *abc, ESL_SQ **ret_sq, P7_TRACE **ret_tr);
    int            p7_alidisplay_Sample(ESL_RANDOMNESS *rng, int N, P7_ALIDISPLAY **ret_ad);
    int            p7_alidisplay_Dump(FILE *fp, const P7_ALIDISPLAY *ad);
    int            p7_alidisplay_Compare(const P7_ALIDISPLAY *ad1, const P7_ALIDISPLAY *ad2);
