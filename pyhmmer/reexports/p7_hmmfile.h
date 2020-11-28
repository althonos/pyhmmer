#include "hmmer.h"

uint32_t v3a_magic = 0xe8ededb6;
uint32_t v3b_magic = 0xe8ededb7;
uint32_t v3c_magic = 0xe8ededb8;
uint32_t v3d_magic = 0xe8ededb9;
uint32_t v3e_magic = 0xe8ededb0;
uint32_t v3f_magic = 0xe8ededba;

extern int read_asc20hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm);
extern int read_asc30hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm);
extern int read_bin30hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm);
