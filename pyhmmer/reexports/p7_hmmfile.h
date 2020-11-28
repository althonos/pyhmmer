#include "hmmer.h"

extern uint32_t v3a_magic;
extern uint32_t v3b_magic;
extern uint32_t v3c_magic;
extern uint32_t v3d_magic;
extern uint32_t v3e_magic;
extern uint32_t v3f_magic;

extern int read_asc20hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm);
extern int read_asc30hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm);
extern int read_bin30hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm);
