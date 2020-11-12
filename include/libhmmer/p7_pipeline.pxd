from libc.stdint cimport uint64_t

from libeasel.sq cimport ESL_SQ
from libeasel.getopts cimport ESL_GETOPTS
from libeasel.random cimport ESL_RANDOMNESS
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_domaindef cimport P7_DOMAINDEF
from libhmmer.p7_tophits cimport P7_TOPHITS
from libhmmer.p7_hmmfile cimport P7_HMMFILE

IF HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_omx cimport P7_OMX
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
ELIF HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_omx cimport P7_OMX
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE


cdef extern from "easel.h" nogil:

    DEF eslERRBUFSIZE = 128


cdef extern from "hmmer.h" nogil:

    cdef enum p7_pipemodes_e:
        p7_SEARCH_SEQS = 0
        p7_SCAN_MODELS = 1

    cdef enum p7_zsetby_e:
        p7_ZSETBY_NTARGETS = 0
        p7_ZSETBY_OPTION = 1
        p7_ZSETBY_FILEINFO = 2

    cdef enum p7_complementarity_e:
        p7_NOCOMPLEMENT    = 0
        p7_COMPLEMENT   = 1

    # ctypedef p7_pipeline_s P7_PIPELINE
    ctypedef struct P7_PIPELINE:
        P7_OMX* oxf
        P7_OMX* oxb
        P7_OMX* fwd
        P7_OMX* bck

        ESL_RANDOMNESS* r
        bint do_reseeding
        bint do_alignment_score_calc
        P7_DOMAINDEF* ddef

        int by_E
        double E
        double T
        bint dom_by_E
        double domE
        double domT
        int use_bit_cutoffs

        bint inc_by_E
        double incE
        double incT
        bint incdom_by_E
        double incdomE
        double incdomT

        double Z
        double domZ
        p7_zsetby_e Z_setby
        p7_zsetby_e domZ_setby

        bint do_max
        double F1
        double F2
        double F3
        int B1
        int B2
        int B3
        bint do_biasfilter
        bint do_null2

        uint64_t nmodels
        uint64_t nseqs
        uint64_t nres
        uint64_t nnodes
        uint64_t n_past_msv
        uint64_t n_past_bias
        uint64_t n_past_vit
        uint64_t n_past_fwd
        uint64_t n_output
        uint64_t pos_past_msv
        uint64_t pos_past_bias
        uint64_t pos_past_vit
        uint64_t pos_past_fwd
        uint64_t pos_output

        p7_pipemodes_e mode
        bint long_targets
        int strands
        int W
        int block_length

        bint show_accessions
        bint show_alignments

        P7_HMMFILE   *hfp
        char          errbuf[eslERRBUFSIZE];


    P7_PIPELINE *p7_pipeline_Create(const ESL_GETOPTS *go, int M_hint, int L_hint, int long_targets, p7_pipemodes_e mode)
    int          p7_pipeline_Reuse  (P7_PIPELINE *pli)
    void         p7_pipeline_Destroy(P7_PIPELINE *pli)
    int          p7_pipeline_Merge  (P7_PIPELINE *p1, P7_PIPELINE *p2)

    # int p7_pli_ExtendAndMergeWindows (P7_OPROFILE *om, const P7_SCOREDATA *msvdata, P7_HMM_WINDOWLIST *windowlist, float pct_overlap)
    int p7_pli_TargetReportable   (P7_PIPELINE *pli, float score,     double lnP)
    int p7_pli_DomainReportable  (P7_PIPELINE *pli, float dom_score, double lnP)

    int p7_pli_TargetIncludable  (P7_PIPELINE *pli, float score,     double lnP)
    int p7_pli_DomainIncludable  (P7_PIPELINE *pli, float dom_score, double lnP)
    int p7_pli_NewModel          (P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg)
    int p7_pli_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om)
    int p7_pli_NewSeq            (P7_PIPELINE *pli, const ESL_SQ *sq)
    int p7_Pipeline              (P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *th)
    # int p7_Pipeline_LongTarget   (P7_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *data,
    #                                      P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx,
    #                                      const ESL_SQ *sq, int complementarity,
    #                                      const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg)
