cdef inline size_t new_capacity(size_t capacity, size_t length) nogil:
    """new_capacity(capacity, length)\n--

    Compute a new capacity for a buffer reallocation.

    """
    return (<size_t> capacity + (capacity >> 3) + 6) & ~(<size_t> 3)
 
