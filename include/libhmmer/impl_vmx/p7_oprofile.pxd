from libc.stdint cimport uint8_t, int16_t
from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.random cimport ESL_RANDOMNESS
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_profile cimport P7_PROFILE


cdef extern from "p7_config.h" nogil:
    DEF p7_MAXABET = 20


cdef extern from "hmmer.h" nogil:
    DEF p7_NEVPARAM = 6
    DEF p7_NCUTOFFS = 6
    DEF p7_NOFFSETS = 3


cdef extern from "impl_vmx/impl_vmx.h" nogil:

    DEF p7O_NXSTATES = 4
    cdef enum p7o_xstates_e:
        p7O_E = 0
        p7O_N = 1
        p7O_J = 2
        p7O_C = 3

    DEF p7O_NXTRANS = 2
    cdef enum p7o_xtransitions_e:
        p7O_MOVE = 0
        p7O_LOOP = 1

    DEF p7O_NTRANS = 8
    cdef enum p7o_tsc_e:
        p7O_BM   = 0
        p7O_MM   = 1
        p7O_IM   = 2
        p7O_DM   = 3
        p7O_MD   = 4
        p7O_MI   = 5
        p7O_II   = 6
        p7O_DD   = 7

    ctypedef p7_oprofile_s P7_OPROFILE
    cdef struct p7_oprofile_s:
        # vector unsigned char** rbv
        uint8_t   tbm_b
        uint8_t   tec_b
        uint8_t   tjb_b
        float     scale_b
        uint8_t   base_b
        uint8_t   bias_b

        # vector signed short **rwv
        # vector signed short  *twv
        int16_t[p7O_NXSTATES][p7O_NXTRANS]   xw
        float     scale_w
        int16_t   base_w
        int16_t   ddbound_w
        float     ncj_roundoff

        # vector float **rfv
        # vector float  *tfv
        float[p7O_NXSTATES][p7O_NXTRANS]    xf

        # vector unsigned char  *rbv_mem
        # vector signed short   *rwv_mem
        # vector signed short   *twv_mem
        # vector float          *tfv_mem
        # vector float          *rfv_mem

        off_t[p7_NOFFSETS]  offs
        off_t  roff
        off_t  eoff

        char  *name
        char  *acc
        char  *desc
        char  *rf
        char  *mm
        char  *cs
        char  *consensus
        float[p7_NEVPARAM]  evparam
        float[p7_NCUTOFFS]  cutoff
        float[p7_MAXABET]  compo
        const ESL_ALPHABET *abc

        int L
        int M
        int max_length
        int allocM
        int allocQ4
        int allocQ8
        int allocQ16
        int mode
        float nj

        int clone

    ctypedef struct P7_OM_BLOCK:
        int count
        int listSize
        P7_OPROFILE** list


    P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc)
    int          p7_oprofile_IsLocal(const P7_OPROFILE *om)
    void         p7_oprofile_Destroy(P7_OPROFILE *om)
    size_t       p7_oprofile_Sizeof(P7_OPROFILE *om)
    P7_OPROFILE *p7_oprofile_Copy(P7_OPROFILE *om)
    P7_OPROFILE *p7_oprofile_Clone(const P7_OPROFILE *om)
    int          p7_oprofile_UpdateFwdEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
    int          p7_oprofile_UpdateVitEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
    int          p7_oprofile_UpdateMSVEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)

    int          p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om)
    int          p7_oprofile_ReconfigLength    (P7_OPROFILE *om, int L)
    int          p7_oprofile_ReconfigMSVLength (P7_OPROFILE *om, int L)
    int          p7_oprofile_ReconfigRestLength(P7_OPROFILE *om, int L)
    int          p7_oprofile_ReconfigMultihit  (P7_OPROFILE *om, int L)
    int          p7_oprofile_ReconfigUnihit    (P7_OPROFILE *om, int L)

    int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om)
    int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
    				       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om)
    int          p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg)
    int          p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm)
    int          p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm)

    int          p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr )
    int          p7_oprofile_GetSSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr )
    int          p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr )
    int          p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr )
