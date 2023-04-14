cdef extern from "esl_distance.h" nogil:

    # Pairwise distances for aligned text sequences.
    int esl_dst_CPairId(const char *asq1, const char *asq2,
    			   double *opt_pid, int *opt_nid, int *opt_n);
    int esl_dst_CPairmatch(const char *asq1, const char *asq2,
    			      double *opt_pmatch, int *opt_nmatch, int *opt_n);
    int esl_dst_CJukesCantor(int K, const char *as1, const char *as2,
    				double *opt_distance, double *opt_variance);

    # Pairwise distances for aligned digital seqs.
    int esl_dst_XPairId(const ESL_ALPHABET *abc, const ESL_DSQ *ax1, const ESL_DSQ *ax2,
    			   double *opt_pid, int *opt_nid, int *opt_n);
    int esl_dst_XPairMatch(const ESL_ALPHABET *abc, const ESL_DSQ *ax1, const ESL_DSQ *ax2,
    			      double *opt_distance, int *opt_nmatch, int *opt_n);
    int esl_dst_XJukesCantor(const ESL_ALPHABET *abc, const ESL_DSQ *ax, const ESL_DSQ *ay,
    				double *opt_distance, double *opt_variance);

    # Distance matrices for aligned text sequences.
    int esl_dst_CPairIdMx     (char **as, int N, ESL_DMATRIX **ret_S);
    int esl_dst_CDiffMx       (char **as, int N, ESL_DMATRIX **ret_D);
    int esl_dst_CJukesCantorMx(int K, char **as, int N, ESL_DMATRIX **opt_D, ESL_DMATRIX **opt_V);

    # Distance matrices for aligned digital sequences.
    int esl_dst_XPairIdMx(const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_S);
    int esl_dst_XDiffMx  (const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D);

    int esl_dst_XJukesCantorMx(const ESL_ALPHABET *abc, ESL_DSQ **ax, int nseq,
    				  ESL_DMATRIX **opt_D, ESL_DMATRIX **opt_V);

    # Average pairwise identity for multiple alignments.
    int esl_dst_CAverageId   (char **as, int nseq, int max_comparisons, double *ret_id);
    int esl_dst_CAverageMatch(char **as, int N, int max_comparisons, double *ret_match);

    int esl_dst_XAverageId   (const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, int max_comparisons, double *ret_id);
    int esl_dst_XAverageMatch(const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, int max_comparisons, double *ret_match);
