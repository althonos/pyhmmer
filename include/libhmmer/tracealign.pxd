from libc.stdio cimport FILE

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.msa cimport ESL_MSA
from libeasel.sq cimport ESL_SQ
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_gmx cimport P7_GMX
from libhmmer.p7_trace cimport P7_TRACE


cdef extern from "hmmer.h" nogil:

  int p7_tracealign_Seqs(ESL_SQ **sq, P7_TRACE **tr, int nseq, int M, int optflags, P7_HMM *hmm, ESL_MSA **ret_msa)
  int p7_tracealign_MSA(const ESL_MSA *premsa, P7_TRACE **tr, int M, int optflags, ESL_MSA **ret_postmsa)
  int p7_tracealign_computeTraces(P7_HMM *hmm, ESL_SQ  **sq, int offset, int N, P7_TRACE  **tr)
  int p7_tracealign_getMSAandStats(P7_HMM *hmm, ESL_SQ  **sq, int N, ESL_MSA **ret_msa, float **ret_pp, float **ret_relent, float **ret_scores )
