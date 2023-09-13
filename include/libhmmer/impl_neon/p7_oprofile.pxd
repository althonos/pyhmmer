from libc.stdint cimport uint8_t, int16_t
from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.random cimport ESL_RANDOMNESS
from libhmmer cimport p7_NEVPARAM, p7_NCUTOFFS, p7_NOFFSETS
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_hmm cimport P7_HMM, p7_MAXABET
from libhmmer.p7_profile cimport P7_PROFILE


cdef extern from "<arm_neon.h>":

    ctypedef struct float32x4_t:
        pass

    ctypedef struct uint8x16_t:
        pass
    
    ctypedef struct int16x8_t:
        pass


cdef extern from "impl_neon/impl_neon.h" nogil:

    cdef int p7O_NQB(int)
    cdef int p7O_NQW(int)
    cdef int p7O_NQF(int)

    const size_t p7O_NXSTATES
    cdef enum p7o_xstates_e:
        p7O_E  = 0
        p7O_N  = 1
        p7O_J  = 2
        p7O_C  = 3

    const size_t p7O_NXTRANS
    cdef enum p7o_xtransitions_e:
        p7O_MOVE = 0
        p7O_LOOP = 1

    const size_t p7O_NTRANS
    cdef enum p7o_tsc_e:
        p7O_BM = 0
        p7O_MM = 1
        p7O_IM = 2
        p7O_DM = 3
        p7O_MD = 4
        p7O_MI = 5
        p7O_II = 6
        p7O_DD = 7


    ctypedef p7_oprofile_s P7_OPROFILE
    cdef struct p7_oprofile_s:
        uint8x16_t** rbv
        uint8x16_t** sbv
        uint8_t tbm_b
        uint8_t tec_b
        uint8_t tjb_b
        float scale_b
        uint8_t base_b
        uint8_t bias_b

        int16x8_t** rwv
        int16x8_t*  twv
        int16_t[p7O_NXSTATES][p7O_NXTRANS] xw
        float scale_w
        int16_t base_w
        int16_t ddbound_w
        float ncj_roundoff

        float32x4_t** rfv
        float32x4_t* tfv
        int16_t[p7O_NXSTATES][p7O_NXTRANS] xf

        uint8x16_t  *rbv_mem
        uint8x16_t  *sbv_mem
        int16x8_t   *rwv_mem
        int16x8_t   *twv_mem
        float32x4_t *tfv_mem
        float32x4_t *rfv_mem

        off_t[p7_NOFFSETS]  offs
        off_t  roff
        off_t  eoff

        char* name
        char* acc
        char* desc
        char* rf
        char* mm
        char* cs
        char* consensus
        float[p7_NEVPARAM] evparam
        float[p7_NCUTOFFS] cutoff
        float[p7_MAXABET] compo
        const ESL_ALPHABET* abc

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
