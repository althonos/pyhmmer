cdef inline size_t new_capacity(size_t capacity, size_t length) nogil:
    """new_capacity(capacity, length)\n--

    Compute a new capacity for a buffer reallocation.

    This is how CPython computes the allocation sizes for the array storing
    the references for a `list` object.

    """
    cdef size_t new_cap = (<size_t> capacity + (capacity >> 3) + 6) & ~(<size_t> 3)
    if (capacity - length) > (new_cap - capacity):
        new_cap = (<size_t> capacity + 3) & ~(<size_t> 3)
    return new_cap