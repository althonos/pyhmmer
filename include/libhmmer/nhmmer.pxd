from libc.stdint cimport uint64_t
from libc.stdlib cimport calloc, malloc, realloc, free

cimport libeasel
from libhmmer.p7_tophits cimport P7_TOPHITS

from capacity cimport new_capacity


ctypedef struct ID_LENGTH:
    int id
    int length


ctypedef struct ID_LENGTH_LIST:
    ID_LENGTH *id_lengths
    int        count
    int        size


cdef inline ID_LENGTH_LIST* idlen_list_init(int size) nogil:
    cdef ID_LENGTH_LIST* l = <ID_LENGTH_LIST*> malloc(sizeof(ID_LENGTH_LIST))
    if l != NULL:
        l.count = 0
        l.size = size
        l.id_lengths = <ID_LENGTH*> calloc(size, sizeof(ID_LENGTH))
        if l.id_lengths == NULL:
            free(l)
            l = NULL
    return l

cdef inline void idlen_list_destroy(ID_LENGTH_LIST* l) nogil:
    if l != NULL:
        if l.id_lengths != NULL:
            free(l.id_lengths)
        free(l)

cdef inline int idlen_list_add(ID_LENGTH_LIST* l, int id, int L) nogil:
    if l.count > 0 and l.id_lengths[l.count - 1].id == id:
        l.id_lengths[l.count - 1].length = L
    else:
        if l.count == l.size:
            l.size = new_capacity(l.size, l.count + 1)
            l.id_lengths = <ID_LENGTH*> realloc(l.id_lengths, l.size * sizeof(ID_LENGTH))
            if l.id_lengths == NULL:
                return libeasel.eslEMEM
        l.id_lengths[l.count].id = id
        l.id_lengths[l.count].length = L
        l.count += 1
    return libeasel.eslOK

cdef inline int idlen_list_assign(ID_LENGTH_LIST* l, P7_TOPHITS* th) nogil:
    cdef uint64_t i
    cdef int      j = 0
    for i in range(th.N):
        while j < l.count and th.hit[i].seqidx != l.id_lengths[j].id:
            j += 1
        if j == l.count:
            return libeasel.eslENOTFOUND
        th.hit[i].dcl[0].ad.L = l.id_lengths[j].length
    return libeasel.eslOK

cdef inline int idlen_list_clear(ID_LENGTH_LIST* l) nogil:
    l.count = 0
