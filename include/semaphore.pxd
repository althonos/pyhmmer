cdef extern from "semaphore.h" nogil:

    ctypedef struct sem_t:
        pass

    int    sem_close(sem_t *)
    int    sem_destroy(sem_t *)
    int    sem_getvalue(sem_t *restrict, int *restrict)
    int    sem_init(sem_t *, bint, unsigned int)
    int    sem_post(sem_t *)
    int    sem_trywait(sem_t *)
    int    sem_unlink(const char *)
    int    sem_wait(sem_t *)
