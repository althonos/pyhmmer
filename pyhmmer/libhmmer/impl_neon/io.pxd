from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel.alphabet cimport ESL_ALPHABET
from libhmmer.p7_hmmfile cimport P7_HMMFILE
from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE
from libhmmer.impl_neon.p7_omx cimport P7_OM_BLOCK


cdef extern from "impl_neon/impl_neon.h" nogil:

  int p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om) except *
  int p7_oprofile_ReadMSV (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
  int p7_oprofile_ReadInfoMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
  int p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock)
  int p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om)
  int p7_oprofile_Position(P7_HMMFILE *hfp, off_t offset)
