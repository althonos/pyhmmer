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
from libhmmer.p7_hit cimport P7_HIT
from libhmmer.p7_pipeline cimport p7_pipemodes_e, P7_PIPELINE
from libhmmer.p7_trace cimport P7_TRACE
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY


cdef extern from "hmmer.h" nogil:

    cdef struct p7_tophits_s:
        P7_HIT** hit
        P7_HIT* unsrt
        uint64_t Nalloc
        uint64_t N
        uint64_t nreported
        uint64_t nincluded
        bint is_sorted_by_sortkey
        bint is_sorted_by_seqidx
    ctypedef p7_tophits_s P7_TOPHITS

    P7_TOPHITS *p7_tophits_Create()
    P7_TOPHITS *p7_tophits_Clone(const P7_TOPHITS* h)
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
    int p7_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header) except *
    int p7_tophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header) except *
    int p7_tophits_TabularXfam(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli)
    int p7_tophits_TabularTail(FILE *ofp, const char *progname, p7_pipemodes_e pipemode,
    				  const char *qfile, const char *tfile, const ESL_GETOPTS *go)
    int p7_tophits_AliScores(FILE *ofp, char *qname, P7_TOPHITS *th )
