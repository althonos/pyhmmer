cdef extern from "esl_matrixops.h" nogil:
    double **esl_mat_DCreate(int M, int N);
    float  **esl_mat_FCreate(int M, int N);
    int    **esl_mat_ICreate(int M, int N);
    char   **esl_mat_CCreate(int M, int N);

    double **esl_mat_DClone(double **A, int M, int N);
    float  **esl_mat_FClone(float **A,  int M, int N);
    int    **esl_mat_IClone(int **A,    int M, int N);

    int      esl_mat_DGrowTo(double ***ret_A, int M, int N);
    int      esl_mat_FGrowTo(float  ***ret_A, int M, int N);
    int      esl_mat_IGrowTo(int    ***ret_A, int M, int N);
    int      esl_mat_CGrowTo(char   ***ret_A, int M, int N);

    size_t   esl_mat_DSizeof(int M, int N);
    size_t   esl_mat_FSizeof(int M, int N);
    size_t   esl_mat_ISizeof(int M, int N);
    size_t   esl_mat_CSizeof(int M, int N);

    void     esl_mat_DSet(double **A, int M, int N, double value);
    void     esl_mat_FSet(float  **A, int M, int N, float  value);
    void     esl_mat_ISet(int    **A, int M, int N, int    value);

    void     esl_mat_DScale(double **A, int M, int N, double x);
    void     esl_mat_FScale(float **A,  int M, int N, float x);
    void     esl_mat_IScale(int **A, int M, int N, int x);

    void     esl_mat_DCopy(double  **src, int M, int N, double  **dest);
    void     esl_mat_FCopy(float   **src, int M, int N, float   **dest);
    void     esl_mat_ICopy(int     **src, int M, int N, int     **dest);
    void     esl_mat_WCopy(int16_t **src, int M, int N, int16_t **dest);
    void     esl_mat_BCopy(int8_t  **src, int M, int N, int8_t  **dest);

    double   esl_mat_DMax(double **A, int M, int N);
    float    esl_mat_FMax(float  **A, int M, int N);
    int      esl_mat_IMax(int    **A, int M, int N);

    int      esl_mat_DCompare(double **A, double **B, int M, int N, double tol);
    int      esl_mat_FCompare(float  **A, float  **B, int M, int N, float  tol);
    int      esl_mat_ICompare(int    **A, int    **B, int M, int N);

    void     esl_mat_DDestroy(double **A);
    void     esl_mat_FDestroy(float  **A);
    void     esl_mat_IDestroy(int    **A);
    void     esl_mat_CDestroy(char   **A);

    int      esl_mat_DDump(double **A, int M, int N);
    int      esl_mat_FDump( float **A, int M, int N);
    int      esl_mat_IDump(   int **A, int M, int N);
