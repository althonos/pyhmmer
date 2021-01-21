from libc.stdio cimport FILE
from libc.stdint cimport uint32_t

from libeasel cimport esl_pos_t


cdef extern from "esl_keyhash.h" nogil:

    ctypedef struct ESL_KEYHASH:
        int* hashtable
        uint32_t hashsize
        int* key_offset
        int* nxt
        int nkeys
        int kalloc

        char* smem
        int salloc
        int sn

    ESL_KEYHASH *esl_keyhash_Create()
    ESL_KEYHASH *esl_keyhash_CreateCustom(uint32_t hashsize, int kalloc, int salloc)
    ESL_KEYHASH *esl_keyhash_Clone(const ESL_KEYHASH *kh)
    char*        esl_keyhash_Get(const ESL_KEYHASH *kh, int idx)
    int          esl_keyhash_GetNumber(const ESL_KEYHASH *kh)
    size_t       esl_keyhash_Sizeof(const ESL_KEYHASH *kh)
    int          esl_keyhash_Reuse(ESL_KEYHASH *kh)
    void         esl_keyhash_Destroy(ESL_KEYHASH *kh)
    void         esl_keyhash_Dump(FILE *fp, const ESL_KEYHASH *kh)

    int  esl_keyhash_Store (      ESL_KEYHASH *kh, const char *key, esl_pos_t n, int *ret_index)
    int  esl_keyhash_Lookup(const ESL_KEYHASH *kh, const char *key, esl_pos_t n, int *ret_index)
