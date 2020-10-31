cdef extern from "esl_msafile_psiblast.h" nogil:
    int esl_msafile_psiblast_SetInmap     (ESL_MSAFILE *afp)
    int esl_msafile_psiblast_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type)
    int esl_msafile_psiblast_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa)
    int esl_msafile_psiblast_Write        (FILE *fp, const ESL_MSA *msa)
