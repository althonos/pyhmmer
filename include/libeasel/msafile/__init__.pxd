from libc.stdint cimport int32_t, int64_t
from libc.stdio cimport FILE

from libeasel cimport ESL_DSQ, esl_pos_t, eslERRBUFSIZE
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.buffer cimport ESL_BUFFER
from libeasel.msa cimport ESL_MSA
from libeasel.ssi cimport ESL_SSI

cdef extern from "esl_msafile.h" nogil:

    ctypedef struct ESL_MSAFILE_FMTDATA:
        int namewidth
        int rpl

    ctypedef struct ESL_MSAFILE:
        ESL_BUFFER* bf
        int32_t format
        ESL_MSAFILE_FMTDATA fmtd

        char* line
        esl_pos_t n
        int64_t linenumber
        esl_pos_t lineoffset

        ESL_DSQ[128] inmap
        const ESL_ALPHABET* abc
        ESL_SSI* ssi
        char[eslERRBUFSIZE] errmsg

    cdef enum:
        eslMSAFILE_UNKNOWN = 0
        eslMSAFILE_STOCKHOLM = 101
        eslMSAFILE_PFAM = 102
        eslMSAFILE_A2M = 103
        eslMSAFILE_PSIBLAST = 104
        eslMSAFILE_SELEX = 105
        eslMSAFILE_AFA = 106
        eslMSAFILE_CLUSTAL = 107
        eslMSAFILE_CLUSTALLIKE = 108
        eslMSAFILE_PHYLIP = 109
        eslMSAFILE_PHYLIPS = 110

    # Opening/closing an ESL_MSAFILE
    extern int   esl_msafile_Open      (ESL_ALPHABET **byp_abc, const char *msafile, const char *env, int format, ESL_MSAFILE_FMTDATA *fmtd, ESL_MSAFILE **ret_afp)
    extern int   esl_msafile_OpenMem   (ESL_ALPHABET **byp_abc, const char *p, esl_pos_t n,           int format, ESL_MSAFILE_FMTDATA *fmtd, ESL_MSAFILE **ret_afp)
    extern int   esl_msafile_OpenBuffer(ESL_ALPHABET **byp_abc, ESL_BUFFER *bf,                       int format, ESL_MSAFILE_FMTDATA *fmtd, ESL_MSAFILE **ret_afp)
    extern void  esl_msafile_OpenFailure(ESL_MSAFILE *afp, int status)
    extern int   esl_msafile_SetDigital (ESL_MSAFILE *afp, const ESL_ALPHABET *abc)
    extern void  esl_msafile_Close(ESL_MSAFILE *afp)

    # ESL_MSAFILE_FMTDATA: optional extra constraints on formats
    extern int   esl_msafile_fmtdata_Init(ESL_MSAFILE_FMTDATA *fmtd)
    extern int   esl_msafile_fmtdata_Copy(ESL_MSAFILE_FMTDATA *src,  ESL_MSAFILE_FMTDATA *dst)

    # Utilities for different file formats
    extern int   esl_msafile_GuessFileFormat(ESL_BUFFER *bf, int *ret_fmtcode, ESL_MSAFILE_FMTDATA *fmtd, char *errbuf)
    extern int   esl_msafile_IsMultiRecord(int fmt)
    extern int   esl_msafile_EncodeFormat(char *fmtstring)
    extern char *esl_msafile_DecodeFormat(int fmt)

    # Utilities for different alphabets
    extern int esl_msafile_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type)

    # Random access in a MSA flatfile database
    extern int esl_msafile_PositionByKey(ESL_MSAFILE *afp, const char *key)

    # Reading an MSA from an ESL_MSAFILE
    extern int  esl_msafile_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
    extern void esl_msafile_ReadFailure(ESL_MSAFILE *afp, int status)

    # Writing an MSA to a stream
    extern int esl_msafile_Write(FILE *fp, ESL_MSA *msa, int fmt)

    # Utilities for specific parsers
    extern int esl_msafile_GetLine(ESL_MSAFILE *afp, char **opt_p, esl_pos_t *opt_n)
    extern int esl_msafile_PutLine(ESL_MSAFILE *afp)


from libeasel.msafile cimport a2m
