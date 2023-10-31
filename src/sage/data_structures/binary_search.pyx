# sage_setup: distribution = sagemath-categories
# We can probably get away with only having the mpz_binary_searches in here.
# I'm too scared to get rid of it at 2am though.
cdef Py_ssize_t binary_search(Py_ssize_t* v, Py_ssize_t n, Py_ssize_t x, Py_ssize_t* ins) noexcept:
    """
    Find the position of the integer x in the array v, which has length n.

    Return -1 if x is not in the array v, and in this case ins is
    set equal to the position where x should be inserted in order to
    obtain an ordered array.
    """
    if n == 0:
        ins[0] = 0
        return -1

    cdef Py_ssize_t i, j, k
    i = 0
    j = n-1
    while i<=j:
        if i == j:
            if v[i] == x:
                ins[0] = i
                return i
            if v[i] < x:
                ins[0] = i + 1
            else:
                ins[0] = i
            return -1
        k = (i+j)/2
        if v[k] > x:
            j = k-1
        elif v[k] < x:
            i = k+1
        else:   # only possibility is that v[k] == x
            ins[0] = k
            return k
    # end while
    ins[0] = j+1
    return -1


cdef Py_ssize_t binary_search0(Py_ssize_t* v, Py_ssize_t n, Py_ssize_t x) noexcept:
    """
    Find the position of the int x in the array v, which has length n.

    Return -1 if x is not in the array v.
    """
    if n == 0:
        return -1

    cdef Py_ssize_t i, j, k
    i = 0
    j = n-1
    while i<=j:
        if i == j:
            if v[i] == x:
                return i
            return -1
        k = (i+j)/2
        if v[k] > x:
            j = k-1
        elif v[k] < x:
            i = k+1
        else:   # only possibility is that v[k] == x
            return k
    return -1
