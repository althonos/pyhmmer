from libc.stdint cimport uint8_t, uint32_t, int64_t

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.dmatrix cimport ESL_DMATRIX
from libeasel.getopts cimport ESL_GETOPTS
from libeasel.msa cimport ESL_MSA
from libeasel.random cimport ESL_RANDOMNESS
from libeasel.scorematrix cimport ESL_SCOREMATRIX
from libeasel.sq cimport ESL_SQ
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_prior cimport P7_PRIOR
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_trace cimport P7_TRACE

IF HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
ELIF HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE


cdef extern from "hmmer.h" nogil:

    DEF eslERRBUFSIZE = 128

    cdef double p7_DEFAULT_WINDOW_BETA

    cdef enum p7_archchoice_e:
        p7_ARCH_FAST = 0
        p7_ARCH_HAND = 1

    cdef enum p7_wgtchoice_e:
        p7_WGT_NONE   = 0
        p7_WGT_GIVEN  = 1
        p7_WGT_GSC    = 2
        p7_WGT_PB     = 3
        p7_WGT_BLOSUM = 4

    cdef enum p7_effnchoice_e:
        p7_EFFN_NONE        = 0
        p7_EFFN_SET         = 1
        p7_EFFN_CLUST       = 2
        p7_EFFN_ENTROPY     = 3
        p7_EFFN_ENTROPY_EXP = 4

    ctypedef p7_builder_s P7_BUILDER
    cdef struct p7_builder_s:
        p7_archchoice_e arch_strategy
        float symfrac
        float fragthresh

        p7_wgtchoice_e wgt_strategy
        double wid

        p7_effnchoice_e effn_strategy
        double re_target
        double esigma
        double eid
        double eset

        ESL_RANDOMNESS* r
        int do_reseeding

        int EmL
        int EmN
        int EvL
        int EvN
        int EfL
        int EfN
        double Eft

        P7_PRIOR* prior
        int max_insert_len

        ESL_SCOREMATRIX* S
        ESL_DMATRIX* Q
        double popen
        double pextend

        double w_beta
        int w_len

        const ESL_ALPHABET* abc
        char errbuf[eslERRBUFSIZE]

    P7_BUILDER *p7_builder_Create(const ESL_GETOPTS *go, const ESL_ALPHABET *abc)
    int         p7_builder_LoadScoreSystem(P7_BUILDER *bld, const char *matrix,                  double popen, double pextend, P7_BG *bg)
    int         p7_builder_SetScoreSystem (P7_BUILDER *bld, const char *mxfile, const char *env, double popen, double pextend, P7_BG *bg)
    void        p7_builder_Destroy(P7_BUILDER *bld)

    int p7_Builder      (P7_BUILDER *bld, ESL_MSA *msa, P7_BG *bg, P7_HMM **opt_hmm, P7_TRACE ***opt_trarr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om, ESL_MSA **opt_postmsa)
    int p7_SingleBuilder(P7_BUILDER *bld, ESL_SQ *sq,   P7_BG *bg, P7_HMM **opt_hmm, P7_TRACE  **opt_tr,    P7_PROFILE **opt_gm, P7_OPROFILE **opt_om)
    int p7_Builder_MaxLength      (P7_HMM *hmm, double emit_thresh)
