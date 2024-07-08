cdef extern from "esl_msacluster.h":
    int esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid, 
    					int **opt_c, int **opt_nin, int *opt_nc);
