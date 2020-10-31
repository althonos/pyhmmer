cdef extern from "esl_msafile_stockholm.h" nogil:
    int esl_msafile_stockholm_SetInmap     (ESL_MSAFILE *afp)
    int esl_msafile_stockholm_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type)
    int esl_msafile_stockholm_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa)
    int esl_msafile_stockholm_Write        (FILE *fp, const ESL_MSA *msa, int fmt)
