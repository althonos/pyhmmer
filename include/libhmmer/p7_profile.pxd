cdef extern from "hmmer.h" nogil:

    cdef enum p7p_xstates_e:
        p7P_E = 0
        p7P_N = 1
        p7P_J = 2
        p7P_C = 3

    cdef enum p7p_xtransitions_e:
        p7P_LOOP = 0
        p7P_MOVE = 1

    cdef enum p7p_tsc_e:
        p7P_MM = 0
        p7P_IM = 1
        p7P_DM = 2
        p7P_BM = 3
        p7P_MD = 4
        p7P_DD = 5
        p7P_MI = 6
        p7P_II = 7

    cdef enum p7p_rsc_e:
        p7P_MSC = 0
        p7P_ISC = 1
