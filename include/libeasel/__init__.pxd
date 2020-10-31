from libc.stdint cimport uint8_t, int64_t

cdef extern from "easel.h" nogil:

    cdef size_t eslERRBUFSIZE

    ctypedef uint8_t ESL_DSQ
    ctypedef int64_t esl_pos_t

    cdef enum:
        eslOK = 0
        eslFAIL = 1
        eslEOL = 2
        eslEOF = 3
        eslEOD = 4
        eslEMEM = 5
        eslENOTFOUND = 6
        eslEFORMAT = 7
        eslEAMBIGUOUS = 8
        eslEDIVZERO = 9
        eslEINCOMPAT = 10
        eslEINVAL = 11
        eslESYS = 12
        eslECORRUPT = 13
        eslEINCONCEIVABLE = 14
        eslESYNTAX = 15
        eslERANGE = 16
        eslEDUP = 17
        eslENOHALT = 18
        eslENORESULT = 19
        eslENODATA = 20
        eslETYPE = 21
        eslEOVERWRITE = 22
        eslENOSPACE = 23
        eslEUNIMPLEMENTED = 24
        eslENOFORMAT = 25
        eslENOALPHABET = 26
        eslEWRITE = 27
        eslEINACCURATE = 28
