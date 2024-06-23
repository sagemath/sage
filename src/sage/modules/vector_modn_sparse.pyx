# sage_setup: distribution = sagemath-modules
#############################################################################
#       Copyright (C) 2004, 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
#############################################################################

from cysignals.memory cimport sig_malloc, sig_free


cdef int allocate_c_vector_modint(c_vector_modint* v, Py_ssize_t num_nonzero) except -1:
    """
    Allocate memory for a c_vector_modint -- most user code won't call this.
    """
    v.entries = <int_fast64_t*>sig_malloc(num_nonzero*sizeof(int_fast64_t))
    if v.entries == NULL:
        raise MemoryError("error allocating memory")
    v.positions = <Py_ssize_t*>sig_malloc(num_nonzero*sizeof(Py_ssize_t))
    if v.positions == NULL:
        sig_free(v.entries)
        raise MemoryError("error allocating memory")
    return 0


cdef int init_c_vector_modint(c_vector_modint* v, int p, Py_ssize_t degree,
                              Py_ssize_t num_nonzero) except -1:
    """
    Initialize a c_vector_modint.
    """
    if (allocate_c_vector_modint(v, num_nonzero) == -1):
        raise MemoryError("error allocating memory for sparse vector")
    if p > 46340:
        clear_c_vector_modint(v)
        raise OverflowError("the prime must be <= 46340")
    v.num_nonzero = num_nonzero
    v.degree = degree
    v.p = p
    return 0


cdef void clear_c_vector_modint(c_vector_modint* v) noexcept:
    sig_free(v.entries)
    sig_free(v.positions)


cdef Py_ssize_t binary_search0_modn(Py_ssize_t* v, Py_ssize_t n, int_fast64_t x) noexcept:
    """
    Find the position of the int x in the array v, which has length n.

    Return -1 if x is not in the array v.
    """
    if n == 0:
        return -1

    cdef Py_ssize_t i, j, k
    i = 0
    j = n-1
    while i <= j:
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


cdef Py_ssize_t binary_search_modn(Py_ssize_t* v, Py_ssize_t n, int_fast64_t x, Py_ssize_t* ins) noexcept:
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
    while i <= j:
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
    ins[0] = j+1
    return -1


cdef int_fast64_t get_entry(c_vector_modint* v, Py_ssize_t n) except -1:
    """
    Return the n-th entry of the sparse vector v.

    This would be v[n] in Python syntax.
    """
    if n >= v.degree or n < 0:
        raise IndexError("index must be between 0 and the degree minus 1")
    cdef Py_ssize_t m
    m = binary_search0_modn(v.positions, v.num_nonzero, n)
    if m == -1:
        return 0
    return v.entries[m]

cdef bint is_entry_zero_unsafe(c_vector_modint* v, Py_ssize_t n) noexcept:
    """
    Return if the ``n``-th entry of the sparse vector ``v`` is zero.

    This is meant for internal use only. If ``n`` is not valid, then
    this might lead to a segfault.
    """
    return binary_search0_modn(v.positions, v.num_nonzero, n) == -1

cdef object to_list(c_vector_modint* v):
    """
    Return a Python list of 2-tuples (i,x), where x=v[i] runs
    through the nonzero elements of x, in order.
    """
    cdef object X
    cdef Py_ssize_t i
    X = []
    for i in range(v.num_nonzero):
        X.append((v.positions[i], v.entries[i]))
    return X


cdef int set_entry(c_vector_modint* v, Py_ssize_t n, int_fast64_t x) except -1:
    """
    Set the n-th component of the sparse vector v equal to x.

    This would be v[n] = x in Python syntax.
    """
    if n < 0 or n >= v.degree:
        raise IndexError("index (=%s) must be between 0 and %s" % (n, v.degree-1))
    cdef Py_ssize_t i, m, ins
    cdef Py_ssize_t *pos
    cdef int_fast64_t *e

    x = x % v.p
    if x < 0:
        x = x + v.p
    m = binary_search_modn(v.positions, v.num_nonzero, n, &ins)

    if m != -1:
        # The position n was found in the array of positions.
        # Now there are two cases:
        #   1. x =/= 0, which is easy, and
        #   2. x = 0, in which case we have to recopy
        #      positions and entries, without the m-th
        #      element, and change num_nonzero.
        if x != 0:   # case 1
            v.entries[m] = x
        else:        # case 2
            e = v.entries
            pos = v.positions
            allocate_c_vector_modint(v, v.num_nonzero - 1)
            for i in range(m):
                v.entries[i] = e[i]
                v.positions[i] = pos[i]
            for i from m < i < v.num_nonzero:
                v.entries[i-1] = e[i]
                v.positions[i-1] = pos[i]
            sig_free(e)
            sig_free(pos)
            v.num_nonzero = v.num_nonzero - 1
    else:
        # Allocate new memory and copy over elements from the
        # old array.  This is similar to case 2 above,
        # except we are inserting a new entry rather than
        # deleting an old one.  The new entry should be inserted
        # at position ins, which was computed using binary search.
        #
        # There is one exception -- if the new entry is 0, we
        # do nothing and return.
        if x == 0:
            return 0
        v.num_nonzero = v.num_nonzero + 1
        e = v.entries
        pos = v.positions
        allocate_c_vector_modint(v, v.num_nonzero)
        for i in range(ins):
            v.entries[i] = e[i]
            v.positions[i] = pos[i]
        v.entries[ins] = x
        v.positions[ins] = n
        for i from ins < i < v.num_nonzero:
            v.entries[i] = e[i-1]
            v.positions[i] = pos[i-1]
        sig_free(e)
        sig_free(pos)


cdef int add_c_vector_modint_init(c_vector_modint* sum, c_vector_modint* v,
                                  c_vector_modint* w, int multiple) except -1:
    """
    Set sum = v + multiple*w.
    """
    if v.p != w.p:
        raise ArithmeticError("the vectors must be modulo the same prime")
    if v.degree != w.degree:
        raise ArithmeticError("the vectors must have the same degree")

    cdef int s
    cdef Py_ssize_t nz, i, j, k
    cdef c_vector_modint* z

    multiple = multiple % v.p    # need this to avoid overflow.
    if multiple < 0:
        multiple = multiple + v.p

    z = sum
    # ALGORITHM:
    # 1. Allocate enough memory to hold the union of the two
    #    lists of positions.  We allocate the sum of the number
    #    of positions of both (up to the degree), to avoid
    #    having to make two passes.  This might be slightly wasteful of
    #    memory, but is faster.
    # 2. Move along the entries of v and w, copying them into the
    #    new position / entry array.  When position are the same,
    #    add modulo p.
    # 3. Set num_nonzero and return success code.

    # 1. Allocate memory:
    nz = v.num_nonzero + w.num_nonzero
    if nz > v.degree:
        nz = v.degree
    init_c_vector_modint(z, v.p, v.degree, nz)
    # 2. Merge entries
    i = 0  # index into entries of v
    j = 0  # index into entries of w
    k = 0  # index into z (the vector we are creating)
    while i < v.num_nonzero or j < w.num_nonzero:
        if i >= v.num_nonzero:   # just copy w in
            z.positions[k] = w.positions[j]
            z.entries[k] = (multiple*w.entries[j])%v.p
            j = j + 1
            k = k + 1
        elif j >= w.num_nonzero:  # just copy v in
            z.positions[k] = v.positions[i]
            z.entries[k] = v.entries[i]
            i = i + 1
            k = k + 1
        elif v.positions[i] < w.positions[j]:  # copy entry from v in
            z.positions[k] = v.positions[i]
            z.entries[k] = v.entries[i]
            i = i + 1
            k = k + 1
        elif v.positions[i] > w.positions[j]: # copy entry from w in
            s = (multiple*w.entries[j])%v.p
            if s != 0:
                z.positions[k] = w.positions[j]
                z.entries[k] = s
                k = k + 1
            j = j + 1
        else:                                 # equal, so add and copy
            s = (v.entries[i] + multiple*w.entries[j]) % v.p
            if s != 0:
                z.positions[k] = v.positions[i]
                z.entries[k] = s
                k = k + 1     # only increment if sum is nonzero!
            i = i + 1
            j = j + 1
    z.num_nonzero = k
    return 0


cdef int scale_c_vector_modint(c_vector_modint* v, int_fast64_t scalar) except -1:
    scalar = scalar % v.p
    if scalar == 0:
        clear_c_vector_modint(v)
        init_c_vector_modint(v, v.p, v.degree, 0)
        return 0
    if scalar < 0:
        scalar = scalar + v.p
    cdef Py_ssize_t i
    for i in range(v.num_nonzero):
        v.entries[i] = (v.entries[i] * scalar) % v.p
    return 0
