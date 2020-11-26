from libc.stdint cimport int64_t, uint32_t, uint64_t
from libc.stdio cimport FILE

from libeasel cimport esl_pos_t
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.getopts cimport ESL_GETOPTS
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.msa cimport ESL_MSA
from libeasel.sq cimport ESL_SQ
from libhmmer.p7_domain cimport P7_DOMAIN
from libhmmer.p7_pipeline cimport p7_pipemodes_e, P7_PIPELINE
from libhmmer.p7_trace cimport P7_TRACE


cdef extern from "hmmer.h" nogil:

    # forward declaration of P7_ALIDISPLAY to avoid circular imports
    ctypedef p7_alidisplay_s P7_ALIDISPLAY
    cdef struct p7_alidisplay_s:
        pass

    # P7_HIT flags
    cdef enum p7_hitflags_e:
        p7_HITFLAGS_DEFAULT = 0b00000
        p7_IS_INCLUDED      = 0b00001
        p7_IS_REPORTED      = 0b00010
        p7_IS_NEW           = 0b00100
        p7_IS_DROPPED       = 0b01000
        p7_IS_DUPLICATE     = 0b10000

    ctypedef p7_hit_s P7_HIT
    cdef struct p7_hit_s:
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

    ctypedef p7_tophits_s P7_TOPHITS
    cdef struct p7_tophits_s:
        P7_HIT** hit
        P7_HIT* unsrt
        uint64_t Nalloc
        uint64_t N
        uint64_t nreported
        uint64_t nincluded
        bint is_sorted_by_sortkey
        bint is_sorted_by_seqidx

    P7_TOPHITS *p7_tophits_Create()
    int         p7_tophits_Grow(P7_TOPHITS *h)
    int         p7_tophits_CreateNextHit(P7_TOPHITS *h, P7_HIT **ret_hit)
    int         p7_tophits_Add(P7_TOPHITS *h,
    				  char *name, char *acc, char *desc,
    				  double sortkey,
    				  float score,    double lnP,
    				  float mothersc, double mother_lnP,
    				  int sqfrom, int sqto, int sqlen,
    				  int hmmfrom, int hmmto, int hmmlen,
    				  int domidx, int ndom,
    				  P7_ALIDISPLAY *ali)
    int         p7_tophits_SortBySortkey(P7_TOPHITS *h)
    int         p7_tophits_SortBySeqidxAndAlipos(P7_TOPHITS *h)
    int         p7_tophits_SortByModelnameAndAlipos(P7_TOPHITS *h)

    int         p7_tophits_Merge(P7_TOPHITS *h1, P7_TOPHITS *h2)
    int         p7_tophits_GetMaxPositionLength(P7_TOPHITS *h)
    int         p7_tophits_GetMaxNameLength(P7_TOPHITS *h)
    int         p7_tophits_GetMaxAccessionLength(P7_TOPHITS *h)
    int         p7_tophits_GetMaxShownLength(P7_TOPHITS *h)
    void        p7_tophits_Destroy(P7_TOPHITS *h)

    int p7_tophits_ComputeNhmmerEvalues(P7_TOPHITS *th, double N, int W)
    int p7_tophits_RemoveDuplicates(P7_TOPHITS *th, int using_bit_cutoffs)
    int p7_tophits_Threshold(P7_TOPHITS *th, P7_PIPELINE *pli)
    int p7_tophits_CompareRanking(P7_TOPHITS *th, ESL_KEYHASH *kh, int *opt_nnew)
    int p7_tophits_Targets(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw)
    int p7_tophits_Domains(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw)


    int p7_tophits_Alignment(const P7_TOPHITS *th, const ESL_ALPHABET *abc,
    				ESL_SQ **inc_sqarr, P7_TRACE **inc_trarr, int inc_n, int optflags,
    				ESL_MSA **ret_msa)
    int p7_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
    int p7_tophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
    int p7_tophits_TabularXfam(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli)
    int p7_tophits_TabularTail(FILE *ofp, const char *progname, p7_pipemodes_e pipemode,
    				  const char *qfile, const char *tfile, const ESL_GETOPTS *go)
    int p7_tophits_AliScores(FILE *ofp, char *qname, P7_TOPHITS *th )
