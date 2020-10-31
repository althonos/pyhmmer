from libc.stdio cimport FILE

from libeasel.msa cimport ESL_MSA
from libeasel.msafile cimport ESL_MSAFILE


cdef extern from "esl_msafile_a2m.h" nogil:
    int esl_msafile_a2m_SetInmap     (ESL_MSAFILE *afp);
    int esl_msafile_a2m_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
    int esl_msafile_a2m_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
    int esl_msafile_a2m_Write        (FILE *fp,    const ESL_MSA *msa);
