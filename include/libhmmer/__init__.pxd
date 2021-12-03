from libc.stdint cimport uint8_t


cdef extern from "p7_config.h" nogil:

    const double p7_ETARGET_AMINO
    const double p7_ETARGET_DNA
    const double p7_ETARGET_OTHER


cdef extern from "hmmer.h" nogil:

    # Search modes
    cdef enum:
        p7_NO_MODE   = 0
        p7_LOCAL     = 1
        p7_GLOCAL    = 2
        p7_UNILOCAL  = 3
        p7_UNIGLOCAL = 4

    cdef bint p7_IsLocal(int)
    cdef bint p7_IsMulti(int)


    cdef size_t p7_NEVPARAM = 6
    cdef enum p7_evparams_e:
        p7_MMU     = 0
        p7_MLAMBDA = 1
        p7_VMU     = 2
        p7_VLAMBDA = 3
        p7_FTAU    = 4
        p7_FLAMBDA = 5

    cdef float p7_EVPARAM_UNSET = -99999.0


    cdef size_t p7_NCUTOFFS = 6
    cdef enum p7_cutoffs_e:
        p7_GA1 = 0
        p7_GA2 = 1
        p7_TC1 = 2
        p7_TC2 = 3
        p7_NC1 = 4
        p7_NC2 = 5

    cdef float p7_CUTOFF_UNSET = -99999.0


    cdef size_t p7_NOFFSETS = 3
    cdef enum p7_offsets_e:
        p7_MOFFSET = 0
        p7_FOFFSET = 1
        p7_POFFSET = 2

    cdef int p7_DEFAULT            =     0
    cdef int p7_DIGITIZE           = (1<<0)
    cdef int p7_ALL_CONSENSUS_COLS = (1<<1)
    cdef int p7_TRIM               = (1<<2)
