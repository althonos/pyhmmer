from libc.stdint cimport uint64_t
from libc.stdio cimport FILE

from libeasel cimport ESL_DSQ
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_gmx cimport P7_GMX


cdef extern from "hmmer.h" nogil:

  int p7_GMSV(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float nu, float *opt_sc)
