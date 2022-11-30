#include <pthread.h>

cdef extern from "pthread.h" nogil:

    ctypedef int pthread_mutex_t

    const pthread_mutex_t PTHREAD_MUTEX_INITIALIZER
    const pthread_mutex_t PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP
    const pthread_mutex_t PTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP

    int pthread_mutex_init(pthread_mutex_t *mutex, const void* mutexattr)
    int pthread_mutex_lock(pthread_mutex_t *mutex)
    int pthread_mutex_trylock(pthread_mutex_t *mutex)
    int pthread_mutex_unlock(pthread_mutex_t *mutex)
    int pthread_mutex_destroy(pthread_mutex_t *mutex)
