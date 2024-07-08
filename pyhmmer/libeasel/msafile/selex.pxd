cdef extern from "esl_msafile_selex.h" nogil:
    int esl_msafile_selex_SetInmap     (ESL_MSAFILE *afp)
    int esl_msafile_selex_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type)
    int esl_msafile_selex_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa)
    int esl_msafile_selex_Write        (FILE *fp,    const ESL_MSA *msa)
