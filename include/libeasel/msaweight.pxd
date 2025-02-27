from libc.stdint cimport int64_t, uint32_t, uint64_t
from posix.types cimport off_t

from libeasel.msa cimport ESL_MSA

cdef extern from "esl_msaweight.h" nogil:


    ctypedef struct ESL_MSAWEIGHT_CFG:
        float fragthresh
        float symfrac
        bint ignore_rf
        bint allow_samp
        int sampthresh
        int nsamp
        int maxfrag
        uint64_t seed
        int filterpref

    const float eslMSAWEIGHT_FRAGTHRESH
    const float eslMSAWEIGHT_SYMFRAC
    const bint eslMSAWEIGHT_IGNORE_RF
    const bint eslMSAWEIGHT_ALLOW_SAMP
    const int eslMSAWEIGHT_SAMPTHRESH
    const int eslMSAWEIGHT_NSAMP
    const int eslMSAWEIGHT_MAXFRAG
    const int eslMSAWEIGHT_RNGSEED

    enum:
        eslMSAWEIGHT_FILT_CONSCOVER
        eslMSAWEIGHT_FILT_RANDOM
        eslMSAWEIGHT_FILT_ORIGORDER

    ctypedef struct ESL_MSAWEIGHT_DAT:
        uint64_t seed
        bint cons_by_rf
        bint cons_by_sample
        bint cons_by_all
        bint cons_allcos
        bint rejected_sample
        int ncons
        int* conscols
        int all_nfrag
        int samp_nfrag


    int esl_msaweight_PB(ESL_MSA *msa)
    int esl_msaweight_PB_adv(const ESL_MSAWEIGHT_CFG *cfg, ESL_MSA *msa, ESL_MSAWEIGHT_DAT *dat)

    ESL_MSAWEIGHT_CFG *esl_msaweight_cfg_Create()
    void               esl_msaweight_cfg_Destroy(ESL_MSAWEIGHT_CFG *cfg)
    ESL_MSAWEIGHT_DAT *esl_msaweight_dat_Create()
    int                esl_msaweight_dat_Reuse  (ESL_MSAWEIGHT_DAT *dat)
    void               esl_msaweight_dat_Destroy(ESL_MSAWEIGHT_DAT *dat)

    int esl_msaweight_GSC(ESL_MSA *msa)
    int esl_msaweight_BLOSUM(ESL_MSA *msa, double maxid)

    int esl_msaweight_IDFilter(const ESL_MSA *msa, double maxid, ESL_MSA **ret_newmsa)
    int esl_msaweight_IDFilter_adv(const ESL_MSAWEIGHT_CFG *cfg, const ESL_MSA *msa, double maxid, ESL_MSA **ret_newmsa)
