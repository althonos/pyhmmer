from libc.stdio cimport FILE

cdef extern from "esl_msafile_afa.h" nogil:
    int esl_msafile_afa_SetInmap     (ESL_MSAFILE *afp);
    int esl_msafile_afa_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
    int esl_msafile_afa_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
    int esl_msafile_afa_Write        (FILE *fp, const ESL_MSA *msa);
