cdef extern from "esl_msafile_phylip.h" nogil:
    int esl_msafile_phylip_SetInmap     (ESL_MSAFILE *afp)
    int esl_msafile_phylip_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type)
    int esl_msafile_phylip_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa)
    int esl_msafile_phylip_Write        (FILE *fp, const ESL_MSA *msa, int format, ESL_MSAFILE_FMTDATA *opt_fmtd)

    int esl_msafile_phylip_CheckFileFormat(ESL_BUFFER *bf, int *ret_format, int *ret_namewidth)
