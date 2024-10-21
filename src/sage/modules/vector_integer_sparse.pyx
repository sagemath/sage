# sage_setup: distribution = sagemath-modules
#############################################################
#
#    Sparse Vector over mpz_t (the GMP integers)
#
#############################################################

from cysignals.memory cimport sig_malloc, sig_free

from sage.libs.gmp.mpz cimport *
from sage.data_structures.binary_search cimport *
from sage.rings.integer cimport Integer


cdef int allocate_mpz_vector(mpz_vector* v, Py_ssize_t num_nonzero) except -1:
    """
    Allocate memory for a mpz_vector -- most user code won't call this.
    num_nonzero is the number of nonzero entries to allocate memory for.

    This function *does* call mpz_init on each mpz_t that is allocated.
    It does *not* clear the entries of v, if there are any.
    """
    cdef Py_ssize_t i
    v.entries = <mpz_t *>sig_malloc(num_nonzero*sizeof(mpz_t))
    if v.entries == NULL:
        raise MemoryError("Error allocating memory")
    for i from 0 <= i < num_nonzero:
        mpz_init(v.entries[i])
    v.positions = <Py_ssize_t*>sig_malloc(num_nonzero*sizeof(Py_ssize_t))
    if v.positions == NULL:
        for i from 0 <= i < num_nonzero:
            mpz_clear(v.entries[i])
        sig_free(v.entries)
        v.entries = NULL
        raise MemoryError("Error allocating memory")
    return 0

cdef int mpz_vector_init(mpz_vector* v, Py_ssize_t degree, Py_ssize_t num_nonzero) except -1:
    """
    Initialize a mpz_vector -- most user code *will* call this.
    """
    allocate_mpz_vector(v, num_nonzero)
    v.num_nonzero = num_nonzero
    v.degree = degree

cdef void mpz_vector_clear(mpz_vector* v) noexcept:
    cdef Py_ssize_t i
    # Free all mpz objects allocated in creating v
    for i from 0 <= i < v.num_nonzero:
        mpz_clear(v.entries[i])
    # Free entries and positions of those entries.
    # These were allocated from the Python heap.
    # If mpz_vector_init was not called, then this
    # will (of course!) cause a core dump.
    sig_free(v.entries)
    sig_free(v.positions)

cdef Py_ssize_t mpz_binary_search0(mpz_t* v, Py_ssize_t n, mpz_t x) noexcept:
    """
    Find the position of the integers x in the array v, which has length n.
    Returns -1 if x is not in the array v.
    """
    cdef Py_ssize_t i, j, k, c
    if n == 0:
        return -1
    i = 0
    j = n-1
    while i<=j:
        if i == j:
            if mpz_cmp(v[i],x) == 0:
                return i
            return -1
        k = (i+j)/2
        c = mpz_cmp(v[k],x)
        if c > 0:       # v[k] > x
            j = k-1
        elif c < 0:     # v[k] < x
            i = k+1
        else:   # only possibility is that v[k] == x
            return k
    return -1

cdef Py_ssize_t mpz_binary_search(mpz_t* v, Py_ssize_t n, mpz_t x, Py_ssize_t* ins) noexcept:
    """
    Find the position of the integer x in the array v, which has length n.
    Returns -1 if x is not in the array v, and in this case ins is
    set equal to the position where x should be inserted in order to
    obtain an ordered array.

    INPUT:

    - ``v`` -- array of mpz_t  (integer)
    - ``n`` -- integer (length of array v)
    - ``x`` -- mpz_t  (integer)

    OUTPUT:

    position of x (as a Py_ssize_t)
    ins -- (call be pointer), the insertion point if x is not found.
    """
    cdef Py_ssize_t i, j, k, c
    if n == 0:
        ins[0] = 0
        return -1
    i = 0
    j = n-1
    while i<=j:
        if i == j:
            c = mpz_cmp(v[i],x)
            if c == 0:          # v[i] == x
                ins[0] = i
                return i
            if c < 0:           # v[i] < x
                ins[0] = i + 1
            else:
                ins[0] = i
            return -1
        k = (i+j)/2
        c = mpz_cmp(v[k], x)
        if c > 0:               # v[k] > x
            j = k-1
        elif c < 0:             # v[k] < x
            i = k+1
        else:   # only possibility is that v[k] == x
            ins[0] = k
            return k
    # end while
    ins[0] = j+1
    return -1

cdef int mpz_vector_get_entry(mpz_t ans, mpz_vector* v, Py_ssize_t n) except -1:
    """
    Returns the n-th entry of the sparse vector v.  This
    would be ``v[n]`` in Python syntax.

    The return is done using the pointer ans, which is to an mpz_t
    that *must* have been initialized using mpz_init.
    """
    if n >= v.degree:
        raise IndexError("Index (=%s) must be between 0 and %s." % (n, v.degree - 1))
    cdef Py_ssize_t m
    m = binary_search0(v.positions, v.num_nonzero, n)
    if m == -1:
        mpz_set_si(ans, 0)
        return 0
    mpz_set(ans, v.entries[m])
    return 0

cdef bint mpz_vector_is_entry_zero_unsafe(mpz_vector* v, Py_ssize_t n) noexcept:
    """
    Return if the ``n``-th entry of the sparse vector ``v`` is zero.

    This is meant for internal use only. If ``n`` is not valid, then
    this might lead to a segfault.
    """
    return binary_search0(v.positions, v.num_nonzero, n) == -1

cdef object mpz_vector_to_list(mpz_vector* v):
    """
    Return a Python list of 2-tuples (i,x), where ``x=v[i]`` runs
    through the nonzero elements of x, in order.
    """
    cdef object X
    cdef Integer a
    cdef Py_ssize_t i
    X = []
    for i from 0 <= i < v.num_nonzero:
        a = Integer()
        a.set_from_mpz(v.entries[i])
        X.append( (v.positions[i], a) )
    return X


cdef int mpz_vector_set_entry(mpz_vector* v, Py_ssize_t n, mpz_t x) except -1:
    """
    Set the n-th component of the sparse vector v equal to x.
    This would be ``v[n] = x`` in Python syntax.
    """
    if n >= v.degree or n < 0:
        raise IndexError("Index (=%s) must be between 0 and %s." % (n, v.degree - 1))
    cdef Py_ssize_t i, m, ins
    cdef Py_ssize_t *pos
    cdef mpz_t *e

    m = binary_search(v.positions, v.num_nonzero, n, &ins)

    if m != -1:
        # The position n was found in the array of positions.
        # Now there are two cases:
        #   1. x =/= 0, which is easy, and
        #   2. x = 0, in which case we have to recopy
        #      positions and entries, without the m-th
        #      element, and change num_nonzero.
        if mpz_sgn(x) != 0:   # x != 0,  case 1
            # v.entries[m] = x
            mpz_set(v.entries[m], x)
        else:        # case 2
            e = v.entries
            pos = v.positions
            allocate_mpz_vector(v, v.num_nonzero - 1)  # This does *not* change v.num_nonzero
            for i in range(m):
                # v.entries[i] = e[i]
                mpz_set(v.entries[i], e[i])
                mpz_clear(e[i])
                v.positions[i] = pos[i]
            for i from m < i < v.num_nonzero:
                # v.entries[i-1] = e[i]
                mpz_set(v.entries[i-1], e[i])
                mpz_clear(e[i])
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
        if mpz_sgn(x) == 0:
            return 0
        v.num_nonzero = v.num_nonzero + 1
        e = v.entries
        pos = v.positions
        allocate_mpz_vector(v, v.num_nonzero)
        for i in range(ins):
            # v.entries[i] = e[i]
            mpz_set(v.entries[i], e[i])
            mpz_clear(e[i])
            v.positions[i] = pos[i]
        # v.entries[ins] = x
        mpz_set(v.entries[ins], x)
        v.positions[ins] = n
        for i from ins < i < v.num_nonzero:
            mpz_set(v.entries[i], e[i-1])
            mpz_clear(e[i-1])
            v.positions[i] = pos[i-1]
        sig_free(e)
        sig_free(pos)


cdef mpz_t mpz_set_tmp
mpz_init(mpz_set_tmp)
cdef int mpz_vector_set_entry_str(mpz_vector* v, Py_ssize_t n, char *x_str) except -1:
    """
    Set the n-th component of the sparse vector v equal to x.
    This would be v[n] = x in Python syntax.
    """
    mpz_set_str(mpz_set_tmp, x_str, 0)
    mpz_vector_set_entry(v, n, mpz_set_tmp)


cdef int add_mpz_vector_init(mpz_vector* sum,
                             mpz_vector* v,
                             mpz_vector* w,
                             mpz_t multiple) except -1:
    """
    Initialize sum and set sum = v + multiple*w.
    """
    if v.degree != w.degree:
        print("Can't add vectors of degree %s and %s" % (v.degree, w.degree))
        raise ArithmeticError("The vectors must have the same degree.")

    cdef Py_ssize_t nz, i, j, k, do_multiply
    cdef mpz_vector* z
    cdef mpz_t tmp
    if mpz_cmp_si(multiple, 0) == 0:
        mpz_vector_init(sum, v.degree, 0)
        return 0

    mpz_init(tmp)
    # Do not do the multiply if the multiple is 1.
    do_multiply = mpz_cmp_si(multiple, 1)

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
    mpz_vector_init(z, v.degree, nz)
    # 2. Merge entries
    i = 0  # index into entries of v
    j = 0  # index into entries of w
    k = 0  # index into z (the vector we are creating)
    while i < v.num_nonzero or j < w.num_nonzero:
        if i >= v.num_nonzero:   # just copy w in
            z.positions[k] = w.positions[j]

            if do_multiply:
                # This means: z.entries[k] = (multiple*w.entries[j])
                mpz_mul(z.entries[k], multiple, w.entries[j])
            else:
                mpz_set(z.entries[k], w.entries[j])
            j = j + 1
            k = k + 1
        elif j >= w.num_nonzero:  # just copy v in
            z.positions[k] = v.positions[i]
            # This means: z.entries[k] = v.entries[i]
            mpz_set(z.entries[k], v.entries[i])
            i = i + 1
            k = k + 1
        elif v.positions[i] < w.positions[j]:  # copy entry from v in
            z.positions[k] = v.positions[i]
            # This means: z.entries[k] = v.entries[i]
            mpz_set(z.entries[k], v.entries[i])
            i = i + 1
            k = k + 1
        elif v.positions[i] > w.positions[j]: # copy entry from w in
            if do_multiply:
                # This means: tmp = multiple*w.entries[j]
                mpz_mul(tmp, multiple, w.entries[j])
                # This means: z.entries[k] = tmp
                mpz_set(z.entries[k], tmp)
            else:
                mpz_set(z.entries[k], w.entries[j])
            z.positions[k] = w.positions[j]
            k = k + 1
            j = j + 1
        else:                                 # equal, so add and copy
            if do_multiply:
                # This means: tmp = v.entries[i] + multiple*w.entries[j]
                mpz_mul(tmp, multiple, w.entries[j])
                mpz_add(tmp, tmp, v.entries[i])
            else:
                mpz_add(tmp, v.entries[i], w.entries[j])
            if mpz_sgn(tmp) != 0:
                z.positions[k] = v.positions[i]
                # This means: z.entries[k] = tmp
                mpz_set(z.entries[k], tmp)
                k = k + 1     # only increment if sum is nonzero!
            i = i + 1
            j = j + 1
        #end if
    # end while
    for i from k <= i < z.num_nonzero:
        mpz_clear(z.entries[i])
    z.num_nonzero = k
    mpz_clear(tmp)
    return 0

cdef int mpz_vector_scale(mpz_vector* v, mpz_t scalar) except -1:
    if mpz_sgn(scalar) == 0:  # scalar = 0
        mpz_vector_clear(v)
        mpz_vector_init(v, v.degree, 0)
        return 0
    cdef Py_ssize_t i
    for i from 0 <= i < v.num_nonzero:
        # v.entries[i] = scalar * v.entries[i]
        mpz_mul(v.entries[i], v.entries[i], scalar)
    return 0

cdef int mpz_vector_scalar_multiply(mpz_vector* v, mpz_vector* w, mpz_t scalar) except -1:
    """
    v = w * scalar
    """
    cdef Py_ssize_t i
    if v == w:
        # rescale self
        return mpz_vector_scale(v, scalar)
    else:
        mpz_vector_clear(v)
        v.entries = <mpz_t*> sig_malloc(w.num_nonzero * sizeof(mpz_t))
        if v.entries == NULL:
            v.positions = NULL
            raise MemoryError("error allocating rational sparse vector mpz's")
        v.positions = <Py_ssize_t*> sig_malloc(w.num_nonzero * sizeof(Py_ssize_t))
        if v.positions == NULL:
            sig_free(v.entries)
            v.entries = NULL
            raise MemoryError("error allocating rational sparse vector positions")
        v.num_nonzero = w.num_nonzero
        v.degree = w.degree
        for i from 0 <= i < v.num_nonzero:
            mpz_init(v.entries[i])
            mpz_mul(v.entries[i], w.entries[i], scalar)
            v.positions[i] = w.positions[i]
        return 0

cdef int mpz_vector_cmp(mpz_vector* v, mpz_vector* w) noexcept:
    if v.degree < w.degree:
        return -1
    elif v.degree > w.degree:
        return 1
    cdef Py_ssize_t i
    cdef int c
    for i from 0 <= i < v.num_nonzero:
        c = mpz_cmp(v.entries[i], w.entries[i])
        if c < 0:
            return -1
        elif c > 0:
            return 1
    return 0
