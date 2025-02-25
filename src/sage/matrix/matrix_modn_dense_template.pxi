r"""
Dense matrices over `\ZZ/n\ZZ` for `n` small using the LinBox library (FFLAS/FFPACK)

FFLAS/FFPACK are libraries to provide BLAS/LAPACK-style routines for
working with finite fields. Additionally, these routines reduce to
BLAS/LAPACK routines using floating point arithmetic.

EXAMPLES::

    sage: A = matrix(GF(127), 7, 7, range(49))
    sage: A*A
    [  2  23  44  65  86 107   1]
    [ 15  85  28  98  41 111  54]
    [ 28  20  12   4 123 115 107]
    [ 41  82 123  37  78 119  33]
    [ 54  17 107  70  33 123  86]
    [ 67  79  91 103 115   0  12]
    [ 80  14  75   9  70   4  65]
    sage: A.rank()
    2

    sage: A = matrix(GF(127), 4, 4, [106, 98, 24, 84, 108, 7, 94, 71, 96, 100, 15, 42, 80, 56, 72, 35])
    sage: A.rank()
    4
    sage: v = vector(GF(127), 4, (100, 93, 47, 110))
    sage: x = A.solve_right(v)
    sage: A*x == v
    True

AUTHORS:

- William Stein (2004-2006): some functions in this file were copied
  from ``matrix_modn_dense.pyx`` which was mainly written by William
  Stein
- Clement Pernet (2010): LinBox related functions in this file were
  taken from linbox-sage.C by Clement Pernet
- Burcin Erocal (2010-2011): most of the functions present in this file
- Martin Albrecht (2011): some polishing, bug fixes, documentation
- Rob Beezer (2011): documentation

TESTS:

We test corner cases for multiplication::

    sage: v0 = vector(GF(3),[])
    sage: v1 = vector(GF(3),[1])
    sage: m00 = matrix(GF(3),0,0,[])
    sage: m01 = matrix(GF(3),0,1,[])
    sage: m10 = matrix(GF(3),1,0,[])
    sage: m11 = matrix(GF(3),1,1,[1])
    sage: good = [ (v0,m00), (v0,m01), (v1,m10), (v1,m11), (m00,v0), (m10,v0), (m01,v1), (m11,v1), (m00,m00), (m01,m10), (m10,m01), (m11,m11) ]
    sage: for v, m in good:
    ....:     print('{} x {} = {}'.format(v, m, v * m))
    () x [] = ()
    () x [] = (0)
    (1) x [] = ()
    (1) x [1] = (1)
    [] x () = ()
    [] x () = (0)
    [] x (1) = ()
    [1] x (1) = (1)
    [] x [] = []
    [] x [] = []
    [] x [] = [0]
    [1] x [1] = [1]

    sage: bad  = [ (v1,m00), (v1,m01), (v0,m10), (v0,m11), (m00,v1), (m10,v1), (m01,v0), (m11,v0), (m01,m01), (m10,m10), (m11,m01), (m10,m11) ]
    sage: for v, m in bad:
    ....:     try:
    ....:         v*m
    ....:         print('Uncaught dimension mismatch!')
    ....:     except (IndexError, TypeError, ArithmeticError):
    ....:         pass
"""

# ****************************************************************************
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
#       Copyright (C) 2011 Rob Beezer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.stdint cimport uint64_t
from cpython.bytes cimport *

from cysignals.memory cimport check_malloc, check_allocarray, check_calloc, sig_malloc, sig_free
from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.libs.gmp.mpz cimport *
from sage.libs.linbox.fflas cimport FFLAS_TRANSPOSE, FflasNoTrans, FflasTrans, \
    FflasRight, vector, list as std_list
from libcpp cimport bool
from sage.parallel.parallelism import Parallelism

cimport sage.rings.fast_arith
cdef sage.rings.fast_arith.arith_int ArithIntObj
ArithIntObj = sage.rings.fast_arith.arith_int()

# for copying/pickling
from libc.string cimport memcpy
from libc.stdio cimport snprintf

from sage.modules.vector_modn_dense cimport Vector_modn_dense

from sage.arith.misc import is_prime
from sage.structure.element cimport (Element, Vector, Matrix,
                                     ModuleElement, RingElement)
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract
from sage.misc.timing import cputime
from sage.misc.verbose import verbose, get_verbose
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.structure.proof.proof import get_flag as get_proof_flag
from sage.structure.richcmp cimport rich_to_bool
from sage.misc.randstate cimport randstate, current_randstate
import sage.matrix.matrix_space as matrix_space
from sage.matrix.args cimport SparseEntry, MatrixArgs_init


from sage.cpython.string cimport char_to_str

cdef long num = 1
cdef bint little_endian = (<char*>(&num))[0]

cdef inline celement_invert(celement a, celement n):
    """
    Invert the finite field element `a` modulo `n`.
    """
    # This is copied from linbox source linbox/field/modular-float.h
    # The extended Euclidean algorithm
    cdef int x_int, y_int, q, tx, ty, temp
    x_int = <int>n
    y_int = <int>a
    tx = 0
    ty = 1

    while y_int != 0:
        # always: gcd (n,residue) = gcd (x_int,y_int)
        #         sx*n + tx*residue = x_int
        #         sy*n + ty*residue = y_int
        q = x_int / y_int  # integer quotient
        temp = y_int
        y_int = x_int - q * y_int
        x_int = temp
        temp = ty
        ty = tx - q * ty
        tx = temp

    if tx < 0:
        tx += <int>n

    # now x_int = gcd (n,residue)
    return <celement>tx

cdef inline bint linbox_is_zero(celement modulus, celement* entries, Py_ssize_t nrows, Py_ssize_t ncols) except -1:
    """
    Return 1 if all entries of this matrix are zero.
    """
    cdef Py_ssize_t i, j
    for i in range(nrows):
        for j in range(ncols):
            if (entries+i*ncols+j)[0] != 0:
                return 0
    return 1

cdef inline linbox_echelonize(celement modulus, celement* entries, Py_ssize_t nrows, Py_ssize_t ncols):
    """
    Return the reduced row echelon form of this matrix.
    """
    if linbox_is_zero(modulus, entries, nrows, ncols):
        return 0, []

    cdef Py_ssize_t i, j
    cdef ModField *F = new ModField(<long>modulus)
    cdef size_t* P = <size_t*>check_allocarray(nrows, sizeof(size_t))
    cdef size_t* Q = <size_t*>check_allocarray(ncols, sizeof(size_t))

    cdef Py_ssize_t r
    cdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')
    cdef bool transform = False
    if nrows * ncols > 1000:
        sig_on()
    if nbthreads > 1 :
        r = pReducedRowEchelonForm(F[0], nrows, ncols, <ModField.Element*>entries, ncols, P, Q, transform, nbthreads)
    else :
        r = ReducedRowEchelonForm(F[0], nrows, ncols, <ModField.Element*>entries, ncols, P, Q)
    if nrows * ncols > 1000:
        sig_off()

    for i in range(nrows):
        for j in range(r):
            (entries+i*ncols+j)[0] = 0
        if i<r:
            (entries + i*(ncols+1))[0] = 1

    applyP(F[0], FflasRight, FflasNoTrans, nrows, 0, r, <ModField.Element*>entries, ncols, Q)

    cdef list pivots = [int(Q[i]) for i in range(r)]

    sig_free(P)
    sig_free(Q)
    del F
    return r, pivots

cdef inline linbox_echelonize_efd(celement modulus, celement* entries, Py_ssize_t nrows, Py_ssize_t ncols):
    # See trac #13878: This is to avoid sending invalid data to linbox,
    # which would yield a segfault in Sage's debug version. TODO: Fix
    # that bug upstream.
    if nrows == 0 or ncols == 0:
        return 0, []

    cdef ModField *F = new ModField(<long>modulus)
    cdef DenseMatrix *A = new DenseMatrix(F[0], nrows, ncols)

    cdef Py_ssize_t i, j
    for i in range(nrows):
        for j in range(ncols):
            A.setEntry(i, j, entries[i*ncols+j])

    cdef Py_ssize_t r = reducedRowEchelonize(A[0])
    for i in range(nrows):
        for j in range(ncols):
            entries[i*ncols+j] = <celement>A.getEntry(i, j)

    cdef Py_ssize_t ii = 0
    cdef list pivots = []
    for i in range(r):
        for j in range(ii, ncols):
            if entries[i*ncols+j] == 1:
                pivots.append(j)
                ii = j+1
                break

    del F
    return r, pivots

cdef inline celement *linbox_copy(celement modulus, celement *entries,  Py_ssize_t nrows, Py_ssize_t ncols) except? NULL:
    """
    Create a copy of the entries array.
    """
    cdef celement *entries_copy = <celement*>check_allocarray(nrows * ncols, sizeof(celement))
    memcpy(entries_copy, entries, sizeof(celement)*nrows*ncols)
    return entries_copy

cdef inline int linbox_rank(celement modulus, celement* entries, Py_ssize_t nrows, Py_ssize_t ncols) except -1:
    """
    Return the rank of this matrix.
    """
    cdef ModField *F = new ModField(<long>modulus)

    cdef celement *cpy = linbox_copy(modulus, entries, nrows, ncols)

    cdef Py_ssize_t r
    cdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')
    if nrows * ncols > 1000:
        sig_on()
    if nbthreads > 1:
        r = pRank(F[0], nrows, ncols, <ModField.Element*>cpy, ncols, nbthreads)
    else:
        r = Rank(F[0], nrows, ncols, <ModField.Element*>cpy, ncols)
    if nrows * ncols > 1000:
        sig_off()
    sig_free(cpy)
    del F
    return r

cdef inline celement linbox_det(celement modulus, celement* entries, Py_ssize_t n) noexcept:
    """
    Return the determinant of this matrix.
    """
    cdef ModField *F = new ModField(<long>modulus)
    cdef celement *cpy = linbox_copy(modulus, entries, n, n)

    cdef celement d = 0
    cdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')

    if n*n > 1000:
        sig_on()
    if nbthreads > 1 :
        pDet(F[0], d, n, <ModField.Element*>cpy, n, nbthreads)
    else :
        Det(F[0], d, n, <ModField.Element*>cpy, n)
    if n*n > 1000:
        sig_off()
    sig_free(cpy)
    del F
    return d

cdef inline celement linbox_matrix_matrix_multiply(celement modulus, celement* ans, celement* A, celement* B, Py_ssize_t m, Py_ssize_t n, Py_ssize_t k)  noexcept:
    """
    C = A*B
    """
    cdef ModField *F = new ModField(<long>modulus)
    cdef ModField.Element one = 0, zero = 0
    F[0].init(one, <int>1)
    F[0].init(zero, <int>0)

    cdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')

    if m*n*k > 100000:
        sig_on()
    if nbthreads > 1 :
        pfgemm(F[0], FflasNoTrans, FflasNoTrans, m, n, k, one,
               <ModField.Element*>A, k, <ModField.Element*>B, n, zero,
               <ModField.Element*>ans, n, nbthreads)
    else:
        fgemm(F[0], FflasNoTrans, FflasNoTrans, m, n, k, one,
              <ModField.Element*>A, k, <ModField.Element*>B, n, zero,
              <ModField.Element*>ans, n)

    if m*n*k > 100000:
        sig_off()

    del F

cdef inline int linbox_matrix_vector_multiply(celement modulus, celement* C, celement* A, celement* b, Py_ssize_t m, Py_ssize_t n, FFLAS_TRANSPOSE trans) noexcept:
    """
    C = A*v
    """
    cdef ModField *F = new ModField(<long>modulus)
    cdef ModField.Element one = 0, zero = 0
    F.init(one, <int>1)
    F.init(zero, <int>0)

    if m*n > 100000:
        sig_on()

    fgemv(F[0], trans,  m, n, one, <ModField.Element*>A, n, <ModField.Element*>b, 1,
          zero, <ModField.Element*>C, 1)

    if m*n > 100000:
        sig_off()

    del F

cdef inline linbox_minpoly(celement modulus, Py_ssize_t nrows, celement* entries):
    """
    Compute the minimal polynomial.
    """
    cdef Py_ssize_t i
    cdef ModField *F = new ModField(<long>modulus)
    cdef vector[ModField.Element] *minP = new vector[ModField.Element]()

    if nrows*nrows > 1000:
        sig_on()
    MinPoly(F[0], minP[0], nrows, <ModField.Element*>entries, nrows)
    if nrows*nrows > 1000:
        sig_off()

    l = [<celement>minP.at(i) for i in range(minP.size())]

    del F
    return l

cdef inline linbox_charpoly(celement modulus, Py_ssize_t nrows, celement* entries):
    """
    Compute the characteristic  polynomial.
    """
    cdef Py_ssize_t i
    cdef ModField *F = new ModField(<long>modulus)
    cdef ModDensePolyRing * R = new ModDensePolyRing(F[0])
    cdef ModDensePoly  P

    cdef celement *cpy = linbox_copy(modulus, entries, nrows, nrows)

    if nrows * nrows > 1000:
        sig_on()
    CharPoly(R[0], P, nrows, <ModField.Element*>cpy, nrows)
    if nrows * nrows > 1000:
        sig_off()

    sig_free(cpy)

    l = []
    for i in range(P.size()):
        l.append(<celement>P[i])

    del F
    del R
    return l


cpdef __matrix_from_rows_of_matrices(X):
    """
    Return a matrix whose row ``i`` is constructed from the entries of
    matrix ``X[i]``.

    INPUT:

    - ``X`` -- a nonempty list of matrices of the same size mod a
      single modulus `n`

    EXAMPLES::

        sage: X = [random_matrix(GF(17), 4, 4) for _ in range(10)]
        sage: Y = X[0]._matrix_from_rows_of_matrices(X)  # indirect doctest
        sage: all(list(Y[i]) == X[i].list() for i in range(10))
        True

    OUTPUT: a single matrix mod ``p`` whose ``i``-th row is ``X[i].list()``.

    .. NOTE::

         Do not call this function directly but use the static method
         ``Matrix_modn_dense_float/double._matrix_from_rows_of_matrices``
    """
    # The code below is just a fast version of the following:
    #     from constructor import matrix
    #     K = X[0].base_ring()
    #     v = sum([y.list() for y in X],[])
    #     return matrix(K, len(X), X[0].nrows()*X[0].ncols(), v)

    cdef Matrix_modn_dense_template T
    cdef Py_ssize_t i, n, m
    n = len(X)

    T = X[0]
    m = T._nrows * T._ncols
    cdef Matrix_modn_dense_template A = T.new_matrix(nrows=n, ncols=m)

    for i from 0 <= i < n:
        T = X[i]
        memcpy(A._entries + i*m, T._entries, sizeof(celement)*m)
    return A


cdef class Matrix_modn_dense_template(Matrix_dense):
    def __cinit__(self, *args, bint zeroed_alloc=True, **kwds):
        cdef long p = self._base_ring.characteristic()
        self.p = p
        if p >= MAX_MODULUS:
            raise OverflowError("p (=%s) must be < %s." % (p, MAX_MODULUS))

        if zeroed_alloc:
            self._entries = <celement *>check_calloc(self._nrows * self._ncols, sizeof(celement))
        else:
            self._entries = <celement *>check_allocarray(self._nrows * self._ncols, sizeof(celement))

        self._matrix = <celement **>check_allocarray(self._nrows, sizeof(celement*))
        cdef unsigned int k
        cdef Py_ssize_t i
        k = 0
        for i in range(self._nrows):
            self._matrix[i] = self._entries + k
            k = k + self._ncols

    def __dealloc__(self):
        """
        TESTS::

            sage: import gc
            sage: for i in range(10):                                                   # needs sage.libs.linbox sage.rings.finite_rings
            ....:      A = random_matrix(GF(7),1000,1000)
            ....:      B = random_matrix(Integers(10),1000,1000)
            ....:      C = random_matrix(GF(16007),1000,1000)
            ....:      D = random_matrix(Integers(1000),1000,1000)
            ....:      del A
            ....:      del B
            ....:      del C
            ....:      del D
            ....:      _ = gc.collect()
        """
        sig_free(self._entries)
        sig_free(self._matrix)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Create a new matrix.

        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- perform modular reduction first?

        EXAMPLES::

            sage: A = random_matrix(GF(3),1000,1000)
            sage: type(A)
            <class 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
            sage: A = random_matrix(Integers(10),1000,1000)
            sage: type(A)
            <class 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
            sage: A = random_matrix(Integers(2^16),1000,1000)
            sage: type(A)
            <class 'sage.matrix.matrix_modn_dense_double.Matrix_modn_dense_double'>

        TESTS::

            sage: Matrix(GF(7), 2, 2, [-1, int(-2), GF(7)(-3), 1/4])
            [6 5]
            [4 2]

            sage: Matrix(GF(6434383), 2, 2, [-1, int(-2), GF(7)(-3), 1/4])              # needs sage.rings.finite_rings
            [6434382 6434381]
            [      4 1608596]

            sage: Matrix(Integers(4618990), 2, 2, [-1, int(-2), GF(7)(-3), 1/7])        # needs sage.rings.finite_rings
            [4618989 4618988]
            [      4 2639423]

            sage: Matrix(IntegerModRing(200), [[int(2**128+1), int(2**256+1), int(2**1024+1)]])        # needs sage.rings.finite_rings
            [ 57 137  17]
        """
        ma = MatrixArgs_init(parent, entries)
        cdef long i, j
        it = ma.iter(convert=False, sparse=True)
        R = ma.base
        p = R.characteristic()

        for t in it:
            se = <SparseEntry>t
            x = se.entry
            v = self._matrix[se.i]
            if type(x) is IntegerMod_int and (<IntegerMod_int>x)._parent is R:
                v[se.j] = <celement>(<IntegerMod_int>x).ivalue
            elif type(x) is Integer:
                if coerce:
                    v[se.j] = mpz_fdiv_ui((<Integer>x).value, p)
                else:
                    v[se.j] = mpz_get_ui((<Integer>x).value)
            elif coerce:
                v[se.j] = R(x)
            else:
                v[se.j] = <celement>x

    cdef long _hash_(self) except -1:
        """
        EXAMPLES::

            sage: B = random_matrix(GF(127),3,3)
            sage: B.set_immutable()
            sage: _ = {B:0} # indirect doctest

            sage: M = random_matrix(GF(7), 10, 10)
            sage: M.set_immutable()
            sage: _ = hash(M)
            sage: MZ = M.change_ring(ZZ)
            sage: MZ.set_immutable()
            sage: hash(MZ) == hash(M)
            True
            sage: MS = M.sparse_matrix()
            sage: MS.set_immutable()
            sage: hash(MS) == hash(M)
            True

        TESTS::

            sage: A = matrix(GF(2),2,0)
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: hash(A)
            0
        """
        cdef long C[5]
        self.get_hash_constants(C)

        cdef long h = 0, k, l
        cdef Py_ssize_t i, j
        cdef celement* row
        sig_on()
        for i in range(self._nrows):
            k = C[0] if i == 0 else C[1] + C[2] * i
            row = self._matrix[i]
            for j in range(self._ncols):
                l = C[3] * (i - j) * (i ^ j)
                h += (k ^ l) * <long>(row[j])
        h *= C[4]
        sig_off()

        if h == -1:
            return -2
        return h

    def _pickle(self):
        """
        Utility function for pickling.

        If the prime is small enough to fit in a byte, then it is
        stored as a contiguous string of bytes (to save
        space). Otherwise, memcpy is used to copy the raw data in the
        platforms native format. Endianness is dealt with when
        unpickling.

        EXAMPLES::

            sage: m = matrix(Integers(128), 3, 3, [ord(c) for c in "Hi there!"]); m
            [ 72 105  32]
            [116 104 101]
            [114 101  33]
            sage: m._pickle()
            ((1, ..., ...'Hi there!'), 10)

        .. todo::

            The upcoming buffer protocol would be useful to not have
            to do any copying.
        """
        cdef Py_ssize_t i, j
        cdef unsigned char* us
        cdef mod_int *um
        cdef unsigned char* row_us
        cdef mod_int *row_um
        cdef long word_size
        cdef celement *row_self

        if self.p <= 0xFF:
            word_size = sizeof(unsigned char)
        else:
            word_size = sizeof(mod_int)

        cdef void *buf = check_allocarray(self._nrows * self._ncols, word_size)

        sig_on()
        try:
            if word_size == sizeof(unsigned char):
                us = <unsigned char*>buf
                for i in range(self._nrows):
                    row_self = self._matrix[i]
                    row_us = us + i*self._ncols
                    for j in range(self._ncols):
                        row_us[j] = <mod_int>row_self[j]
            else:
                um = <mod_int*>buf
                for i in range(self._nrows):
                    row_self = self._matrix[i]
                    row_um = um + i*self._ncols
                    for j in range(self._ncols):
                        row_um[j] = <mod_int>row_self[j]

            s = PyBytes_FromStringAndSize(<char*>buf, word_size * self._nrows * self._ncols)
        finally:
            sig_free(buf)
            sig_off()
        return (word_size, little_endian, s), 10

    def _unpickle(self, data, int version):
        """
        TESTS:

        Test for char-sized modulus::

            sage: A = random_matrix(GF(7), 5, 9)
            sage: data, version = A._pickle()
            sage: B = A.parent()(0)
            sage: B._unpickle(data, version)
            sage: B == A
            True

        And for larger modulus::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(1009), 51, 5)
            sage: data, version = A._pickle()
            sage: B = A.parent()(0)
            sage: B._unpickle(data, version)
            sage: B == A
            True

        Now test all the bit-packing options::

            sage: A = matrix(Integers(1000), 2, 2)
            sage: A._unpickle((1, True, b'\x01\x02\xFF\x00'), 10)
            sage: A
            [  1   2]
            [255   0]

            sage: A = matrix(Integers(1000), 1, 2)
            sage: A._unpickle((4, True, b'\x02\x01\x00\x00\x01\x00\x00\x00'), 10)
            sage: A
            [258   1]
            sage: A._unpickle((4, False, b'\x00\x00\x02\x01\x00\x00\x01\x03'), 10)
            sage: A
            [513 259]
            sage: A._unpickle((8, True, b'\x03\x01\x00\x00\x00\x00\x00\x00\x05\x00\x00\x00\x00\x00\x00\x00'), 10)
            sage: A
            [259   5]
            sage: A._unpickle((8, False, b'\x00\x00\x00\x00\x00\x00\x02\x08\x00\x00\x00\x00\x00\x00\x01\x04'), 10)
            sage: A
            [520 260]

        Now make sure it works in context::

            sage: A = random_matrix(Integers(33), 31, 31)
            sage: loads(dumps(A)) == A
            True
            sage: A = random_matrix(Integers(3333), 31, 31)
            sage: loads(dumps(A)) == A
            True
        """
        if version < 10:
            return Matrix_dense._unpickle(self, data, version)

        cdef Py_ssize_t i, j
        cdef unsigned char* us
        cdef long word_size
        cdef celement *row_self
        cdef bint little_endian_data
        cdef char* buf
        cdef Py_ssize_t buflen
        cdef Py_ssize_t expectedlen
        cdef mod_int v

        if version == 10:
            word_size, little_endian_data, s = data
            expectedlen = word_size * self._nrows * self._ncols

            PyBytes_AsStringAndSize(s, &buf, &buflen)
            if buflen != expectedlen:
                raise ValueError("incorrect size in matrix pickle (expected %d, got %d)" % (expectedlen, buflen))

            sig_on()
            try:
                if word_size == 1:
                    us = <unsigned char*>buf
                    for i from 0 <= i < self._nrows:
                        row_self = self._matrix[i]
                        for j from 0 <= j < self._ncols:
                            row_self[j] = <celement>(us[0])
                            us += word_size

                elif word_size >= 4 and little_endian_data:
                    us = <unsigned char*>buf
                    for i from 0 <= i < self._nrows:
                        row_self = self._matrix[i]
                        for j from 0 <= j < self._ncols:
                            v = <mod_int>(us[0])
                            v += <mod_int>(us[1]) << 8
                            v += <mod_int>(us[2]) << 16
                            v += <mod_int>(us[3]) << 24
                            row_self[j] = <celement>v
                            us += word_size

                elif word_size >= 4 and not little_endian_data:
                    us = <unsigned char*>buf
                    for i from 0 <= i < self._nrows:
                        row_self = self._matrix[i]
                        for j from 0 <= j < self._ncols:
                            v = <mod_int>(us[word_size-1])
                            v += <mod_int>(us[word_size-2]) << 8
                            v += <mod_int>(us[word_size-3]) << 16
                            v += <mod_int>(us[word_size-4]) << 24
                            row_self[j] = <celement>v
                            us += word_size

                else:
                    raise ValueError("unknown matrix pickle format")
            finally:
                sig_off()
        else:
            raise ValueError("unknown matrix pickle version")

    def __neg__(self):
        """
        EXAMPLES::

            sage: A = matrix(GF(19), 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]

            sage: -A
            [ 0 18 17]
            [16 15 14]
            [13 12 11]
        """
        cdef Py_ssize_t i, j
        cdef Matrix_modn_dense_template M
        cdef celement p = self.p

        M = self.__class__.__new__(self.__class__, self._parent,
                                   None, None, None, zeroed_alloc=False)

        sig_on()
        for i in range(self._nrows*self._ncols):
            if self._entries[i]:
                M._entries[i] = p - self._entries[i]
            else:
                M._entries[i] = 0
        sig_off()
        return M

    cpdef _lmul_(self, Element left):
        """
        EXAMPLES::

            sage: A = matrix(GF(101), 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A * 5
            [ 0  5 10]
            [15 20 25]
            [30 35 40]
            sage: A * 50
            [  0  50 100]
            [ 49  99  48]
            [ 98  47  97]

        ::

            sage: A = random_matrix(Integers(60), 400, 500)
            sage: 3*A + 9*A == 12*A
            True
        """
        cdef Py_ssize_t i, j
        cdef Matrix_modn_dense_template M
        cdef celement p = self.p
        cdef celement a = left

        M = self.__class__.__new__(self.__class__, self._parent,
                                   None, None, None, zeroed_alloc=False)

        sig_on()
        for i in range(self._nrows*self._ncols):
            M._entries[i] = (a*self._entries[i]) % p
        sig_off()
        return M

    def __copy__(self):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(127), 100, 100)
            sage: copy(A) == A
            True
            sage: copy(A) is A
            False
        """
        cdef Matrix_modn_dense_template A
        A = self.__class__.__new__(self.__class__, self._parent,
                                   None, None, None, zeroed_alloc=False)
        memcpy(A._entries, self._entries, sizeof(celement)*self._nrows*self._ncols)
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    cpdef _add_(self, right):
        r"""
        Add two dense matrices over `\Z/n\Z`.

        INPUT:

        - ``right`` -- a matrix

        EXAMPLES::

            sage: A = MatrixSpace(GF(19),3)(range(9))
            sage: A+A
            [ 0  2  4]
            [ 6  8 10]
            [12 14 16]

            sage: B = MatrixSpace(GF(19),3)(range(9))
            sage: B.swap_rows(1,2)
            sage: A+B
            [ 0  2  4]
            [ 9 11 13]
            [ 9 11 13]

            sage: B+A
            [ 0  2  4]
            [ 9 11 13]
            [ 9 11 13]
        """
        cdef Py_ssize_t i
        cdef celement k, p
        cdef Matrix_modn_dense_template M

        M = self.__class__.__new__(self.__class__, self._parent,
                                   None, None, None, zeroed_alloc=False)
        p = self.p
        cdef celement* other_ent = (<Matrix_modn_dense_template>right)._entries

        sig_on()
        for i in range(self._nrows*self._ncols):
            k = self._entries[i] + other_ent[i]
            M._entries[i] = k - (k >= p) * p
        sig_off()
        return M

    cpdef _sub_(self, right):
        r"""
        Subtract two dense matrices over `\Z/n\Z`.

        EXAMPLES::

            sage: A = matrix(GF(11), 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]

            sage: A - 4
            [7 1 2]
            [3 0 5]
            [6 7 4]

            sage: A - matrix(GF(11), 3, 3, range(1, 19, 2))
            [10  9  8]
            [ 7  6  5]
            [ 4  3  2]
        """
        cdef Py_ssize_t i
        cdef celement k, p
        cdef Matrix_modn_dense_template M

        M = self.__class__.__new__(self.__class__, self._parent, None, None, None, zeroed_alloc=False)
        p = self.p
        cdef celement* other_ent = (<Matrix_modn_dense_template>right)._entries

        sig_on()
        for i in range(self._nrows*self._ncols):
            k = p + self._entries[i] - other_ent[i]
            M._entries[i] = k - (k >= p) * p
        sig_off()
        return M

    cpdef _richcmp_(self, right, int op):
        r"""
        Compare two dense matrices over `\Z/n\Z`.

        EXAMPLES::

            sage: A = matrix(GF(17), 4, range(3, 83, 5)); A
            [ 3  8 13  1]
            [ 6 11 16  4]
            [ 9 14  2  7]
            [12  0  5 10]
            sage: A == A
            True
            sage: B = A - 3; B
            [ 0  8 13  1]
            [ 6  8 16  4]
            [ 9 14 16  7]
            [12  0  5  7]
            sage: B < A
            True
            sage: B > A
            False
            sage: B == A
            False
            sage: B + 3 == A
            True

        ::

            sage: A = matrix(ZZ, 10, 10, range(1000, 1100))
            sage: A.change_ring(GF(17)) == A.change_ring(GF(17))
            True
            sage: A.change_ring(GF(17)) == A.change_ring(GF(19))                        # needs sage.rings.finite_rings
            False
            sage: A.change_ring(GF(17)) == A.change_ring(Integers(2000))                # needs sage.rings.finite_rings
            False
            sage: A.change_ring(GF(17)) == A.change_ring(Integers(2000))                # needs sage.rings.finite_rings
            False
        """
        cdef Py_ssize_t i
        cdef celement* other_ent = (<Matrix_modn_dense_template>right)._entries
        sig_on()
        for i in range(self._nrows * self._ncols):
            if self._entries[i] < other_ent[i]:
                sig_off()
                return rich_to_bool(op, -1)
            elif self._entries[i] > other_ent[i]:
                sig_off()
                return rich_to_bool(op, 1)
        sig_off()
        return rich_to_bool(op, 0)

    cdef _matrix_times_matrix_(self, Matrix right):
        """
        Return ``self*right``.

        INPUT:

        - ``right``-  a matrix

        EXAMPLES::

            sage: A = random_matrix(GF(7),2,2)
            sage: B = random_matrix(GF(7),2,2)
            sage: C = A*B
            sage: all(C[i, j] == sum(A[i, k]*B[k, j] for k in range(2)) for i in range(2) for j in range(2))
            True

            sage: MS = parent(A)
            sage: MS(3) * A == 3*A
            True

        ::

            sage: A = random_matrix(GF(17), 201, 117)
            sage: B = random_matrix(GF(17), 117, 195)
            sage: C = random_matrix(GF(17), 201, 117)
            sage: D = random_matrix(GF(17), 117, 195)

            sage: E = (A+C)*(B+D)

            sage: F = A*B + A*D + C*B + C*D

            sage: E == F
            True

            sage: A = random_matrix(GF(17), 200, 200)
            sage: MS = parent(A)
            sage: (MS(0) * A) == 0
            True

            sage: (MS(1) * A) == A
            True

        ::

            sage: A = random_matrix(Integers(8),2,2)
            sage: B = random_matrix(Integers(8),2,2)
            sage: C = A*B
            sage: all(C[i, j] == sum(A[i, k]*B[k, j] for k in range(2)) for i in range(2) for j in range(2))
            True

            sage: MS = parent(A)
            sage: MS(3) * A == 3*A
            True

        ::

            sage: A = random_matrix(Integers(16), 201, 117)
            sage: B = random_matrix(Integers(16), 117, 195)
            sage: C = random_matrix(Integers(16), 201, 117)
            sage: D = random_matrix(Integers(16), 117, 195)

            sage: E = (A+C)*(B+D)

            sage: F = A*B + A*D + C*B + C*D

            sage: E == F
            True

            sage: A = random_matrix(Integers(16), 200, 200)
            sage: MS = parent(A)
            sage: (MS(0) * A) == 0
            True

            sage: (MS(1) * A) == A
            True

        ::

            sage: A = random_matrix(GF(16007),2,2)                                      # needs sage.rings.finite_rings
            sage: B = random_matrix(GF(16007),2,2)                                      # needs sage.rings.finite_rings
            sage: C = A*B
            sage: all(C[i, j] == sum(A[i, k]*B[k, j] for k in range(2)) for i in range(2) for j in range(2))
            True

            sage: MS = parent(A)
            sage: MS(3) * A == 3*A
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(15991), 201, 117)
            sage: B = random_matrix(GF(15991), 117, 195)
            sage: C = random_matrix(GF(15991), 201, 117)
            sage: D = random_matrix(GF(15991), 117, 195)

            sage: E = (A+C)*(B+D)

            sage: F = A*B + A*D + C*B + C*D

            sage: E == F
            True

        ::

            sage: A = random_matrix(GF(16007), 200, 200)                                # needs sage.rings.finite_rings
            sage: MS = parent(A)
            sage: (MS(0) * A) == 0
            True

            sage: (MS(1) * A) == A
            True

        ::

            sage: A = random_matrix(Integers(1008),2,2)
            sage: B = random_matrix(Integers(1008),2,2)
            sage: C = A*B
            sage: all(C[i, j] == sum(A[i, k]*B[k, j] for k in range(2)) for i in range(2) for j in range(2))
            True

            sage: MS = parent(A)
            sage: MS(3) * A == 3*A
            True

        ::

            sage: A = random_matrix(Integers(1600), 201, 117)
            sage: B = random_matrix(Integers(1600), 117, 195)
            sage: C = random_matrix(Integers(1600), 201, 117)
            sage: D = random_matrix(Integers(1600), 117, 195)

            sage: E = (A+C)*(B+D)

            sage: F = A*B + A*D + C*B + C*D

            sage: E == F
            True
        """
        if get_verbose() >= 2:
            verbose('mod-p multiply of %s x %s matrix by %s x %s matrix modulo %s' % (
                    self._nrows, self._ncols, right._nrows, right._ncols, self.p))

        if self._ncols != right._nrows:
            raise ArithmeticError("right's number of rows must match self's number of columns")

        cdef int e
        cdef Matrix_modn_dense_template ans, B

        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())

        B = right

        linbox_matrix_matrix_multiply(self.p, ans._entries, self._entries,
                                      B._entries, self._nrows, B._ncols, B._nrows)

        return ans

    cdef _vector_times_matrix_(self, Vector v):
        """
        Return ``v*self``.

        INPUT:

        - ``v`` -- a vector

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 20)
            sage: v = random_vector(GF(17), 10)
            sage: matrix(v*A) == matrix(v)*A
            True

            sage: A = random_matrix(Integers(126), 10, 20)
            sage: v = random_vector(Integers(126), 10)
            sage: matrix(v*A) == matrix(v)*A
            True

            sage: A = random_matrix(GF(4796509), 10, 20)                                # needs sage.rings.finite_rings
            sage: v = random_vector(GF(4796509), 10)                                    # needs sage.rings.finite_rings
            sage: matrix(v*A) == matrix(v)*A
            True

            sage: A = random_matrix(Integers(16337), 10, 20)
            sage: v = random_vector(Integers(16337), 10)                                # needs sage.rings.finite_rings
            sage: matrix(v*A) == matrix(v)*A
            True
        """
        if not isinstance(v, Vector_modn_dense):
            return (self.new_matrix(1, self._nrows, entries=v.list()) * self)[0]

        M = self.row_ambient_module()
        cdef Vector_modn_dense c = M.zero_vector()

        if self._ncols == 0 or self._nrows == 0:
            return c

        cdef Py_ssize_t i
        cdef Vector_modn_dense b = v

        cdef celement *_b = <celement*>check_allocarray(self._nrows, sizeof(celement))
        cdef celement *_c = <celement*>check_allocarray(self._ncols, sizeof(celement))

        for i in range(self._nrows):
            _b[i] = <celement>b._entries[i]

        linbox_matrix_vector_multiply(self.p, _c, self._entries, _b, self._nrows, self._ncols, FflasTrans)

        for i in range(self._ncols):
            c._entries[i] = <mod_int>_c[i]
        sig_free(_b)
        sig_free(_c)
        return c

    cdef _matrix_times_vector_(self, Vector v):
        """
        Return ``self*v``.

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 20)
            sage: v = random_vector(GF(17), 20)
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True

            sage: A = random_matrix(Integers(126), 10, 20)
            sage: v = random_vector(Integers(126), 20)
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True

            sage: A = random_matrix(GF(4796509), 10, 20)                                # needs sage.rings.finite_rings
            sage: v = random_vector(GF(4796509), 20)                                    # needs sage.rings.finite_rings
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True

            sage: A = random_matrix(Integers(16337), 10, 20)
            sage: v = random_vector(Integers(16337), 20)                                # needs sage.rings.finite_rings
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True
        """
        if not isinstance(v, Vector_modn_dense):
            r = (self * self.new_matrix(nrows=len(v), ncols=1, entries=v.list()))
            from sage.modules.free_module_element import vector
            return vector(r.list())

        M = self.column_ambient_module()
        cdef Vector_modn_dense c = M.zero_vector()

        if self._ncols == 0 or self._nrows == 0:
            return c

        cdef Py_ssize_t i
        cdef Vector_modn_dense b = v

        cdef celement *_b = <celement*>check_allocarray(self._ncols, sizeof(celement))
        cdef celement *_c = <celement*>check_allocarray(self._nrows, sizeof(celement))

        for i in range(self._ncols):
            _b[i] = <celement>b._entries[i]

        linbox_matrix_vector_multiply(self.p, _c, self._entries, _b, self._nrows, self._ncols, FflasNoTrans)

        for i in range(self._nrows):
            c._entries[i] = <mod_int>_c[i]
        sig_free(_b)
        sig_free(_c)
        return c

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #  x * cdef _sub_
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #        - all row/column operations, but optimized
    # x      - echelon form in place
    #        - Hessenberg forms of matrices
    ########################################################################

    def charpoly(self, var='x', algorithm='linbox'):
        """
        Return the characteristic polynomial of ``self``.

        INPUT:

        - ``var`` -- a variable name

        - ``algorithm`` -- 'generic', 'linbox' or 'all' (default: linbox)

        EXAMPLES::

            sage: A = random_matrix(GF(19), 10, 10)
            sage: B = copy(A)
            sage: char_p = A.characteristic_polynomial()
            sage: char_p(A) == 0
            True
            sage: B == A              # A is not modified
            True

            sage: min_p = A.minimal_polynomial(proof=True)
            sage: min_p.divides(char_p)
            True

        ::

            sage: A = random_matrix(GF(2916337), 7, 7)                                  # needs sage.rings.finite_rings
            sage: B = copy(A)
            sage: char_p = A.characteristic_polynomial()
            sage: char_p(A) == 0
            True
            sage: B == A               # A is not modified
            True

            sage: min_p = A.minimal_polynomial(proof=True)
            sage: min_p.divides(char_p)
            True

            sage: A = Mat(Integers(6),3,3)(range(9))
            sage: A.charpoly()
            x^3

        TESTS::

            sage: for i in range(10):
            ....:     A = random_matrix(GF(17), 50, 50, density=0.1)
            ....:     _ = A.characteristic_polynomial(algorithm='all')

            sage: A = random_matrix(GF(19), 0, 0)
            sage: A.minimal_polynomial()
            1

            sage: A = random_matrix(GF(19), 0, 1)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = random_matrix(GF(19), 1, 0)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = matrix(GF(19), 10, 10)
            sage: A.minimal_polynomial()                                                # needs sage.libs.pari
            x

            sage: A = random_matrix(GF(4198973), 0, 0)                                  # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            1

            sage: A = random_matrix(GF(4198973), 0, 1)                                  # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = random_matrix(GF(4198973), 1, 0)                                  # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = matrix(GF(4198973), 10, 10)                                       # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            x

            sage: A = Mat(GF(7),3,3)([0, 1, 2] * 3)
            sage: A.charpoly()
            x^3 + 4*x^2

        ALGORITHM: Uses LinBox if ``self.base_ring()`` is a field,
        otherwise use Hessenberg form algorithm.

        TESTS:

        The cached polynomial should be independent of the ``var``
        argument (:issue:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached. The polynomial here is not unique, so we only
        check the polynomial's variable.

            sage: M = MatrixSpace(Integers(37), 2)
            sage: A = M(range(0, 2^2))
            sage: type(A)
            <class 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
            sage: A.charpoly('x').variables()
            (x,)
            sage: A.charpoly('y').variables()
            (y,)
            sage: A._cache['charpoly_linbox'].variables()
            (x,)
        """
        cache_key = 'charpoly_%s' % algorithm
        g = self.fetch(cache_key)
        if g is not None:
            return g.change_variable_name(var)

        if algorithm == 'linbox' and (self.p == 2 or not self.base_ring().is_field()):
            algorithm = 'generic'  # LinBox only supports Z/pZ (p prime)

        if algorithm == 'linbox':
            g = self._charpoly_linbox(var)
        elif algorithm == 'generic':
            g = Matrix_dense.charpoly(self, var)
        elif algorithm == 'all':
            g = self._charpoly_linbox(var)
            h = Matrix_dense.charpoly(self, var)
            if g != h:
                raise ArithmeticError("Characteristic polynomials do not match.")
        else:
            raise ValueError("no algorithm '%s'" % algorithm)

        self.cache(cache_key, g)
        return g

    def minpoly(self, var='x', algorithm='linbox', proof=None):
        """
        Return the minimal polynomial of ``self``.

        INPUT:

        - ``var`` -- a variable name

        - ``algorithm`` -- ``generic`` or ``linbox`` (default:
          ``linbox``)

        - ``proof`` -- (default: ``True``) whether to provably return
          the true minimal polynomial; if ``False``, we only guarantee
          to return a divisor of the minimal polynomial.  There are
          also certainly cases where the computed results is
          frequently not exactly equal to the minimal polynomial (but
          is instead merely a divisor of it).

         .. warning::

             If ``proof=True``, minpoly is insanely slow compared to
             ``proof=False``. This matters since proof=True is the
             default, unless you first type
             ``proof.linear_algebra(False)``.

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 10)
            sage: B = copy(A)
            sage: min_p = A.minimal_polynomial(proof=True)
            sage: min_p(A) == 0
            True
            sage: B == A
            True

            sage: char_p = A.characteristic_polynomial()
            sage: min_p.divides(char_p)
            True

        ::

            sage: A = random_matrix(GF(1214471), 10, 10)                                # needs sage.rings.finite_rings
            sage: B = copy(A)
            sage: min_p = A.minimal_polynomial(proof=True)
            sage: min_p(A) == 0
            True
            sage: B == A
            True

            sage: char_p = A.characteristic_polynomial()
            sage: min_p.divides(char_p)
            True

        TESTS::

            sage: A = random_matrix(GF(17), 0, 0)
            sage: A.minimal_polynomial()
            1

            sage: A = random_matrix(GF(17), 0, 1)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = random_matrix(GF(17), 1, 0)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = matrix(GF(17), 10, 10)
            sage: A.minimal_polynomial()                                                # needs sage.libs.pari
            x

        ::

            sage: A = random_matrix(GF(2535919), 0, 0)                                  # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            1

            sage: A = random_matrix(GF(2535919), 0, 1)                                  # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = random_matrix(GF(2535919), 1, 0)                                  # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = matrix(GF(2535919), 10, 10)                                       # needs sage.rings.finite_rings
            sage: A.minimal_polynomial()                                                # needs sage.rings.finite_rings
            x

        EXAMPLES::

            sage: R.<x>=GF(3)[]
            sage: A = matrix(GF(3),2,[0,0,1,2])
            sage: A.minpoly()
            x^2 + x

            sage: A.minpoly(proof=False) in [x, x+1, x^2+x]
            True
        """
        proof = get_proof_flag(proof, "linear_algebra")

        if algorithm == 'linbox' and (self.p == 2 or not self.base_ring().is_field()):
            algorithm='generic'  # LinBox only supports fields

        if algorithm == 'linbox':
            if self._nrows != self._ncols:
                raise ValueError("matrix must be square")

            if self._nrows <= 1:
                return Matrix_dense.minpoly(self, var)

            R = self._base_ring[var]
            v = linbox_minpoly(self.p, self._nrows, self._entries)
            g = R(v)

            if proof:
                while g(self):  # insanely toy slow (!)
                    g = g.lcm(R(linbox_minpoly(self.p, self._nrows, self._entries)))

        elif algorithm == 'generic':
            raise NotImplementedError("Minimal polynomials are not implemented for Z/nZ.")

        else:
            raise ValueError("no algorithm '%s'" % algorithm)

        self.cache('minpoly_%s_%s' % (algorithm, var), g)
        return g

    def _charpoly_linbox(self, var='x'):
        """
        Compute the characteristic polynomial using LinBox. No checks
        are performed.

        This function is called internally by ``charpoly``.

        INPUT:

        - ``var`` -- a variable name

        EXAMPLES::

            sage: A = random_matrix(GF(19), 10, 10)
            sage: B = copy(A)
            sage: char_p = A._charpoly_linbox()
            sage: char_p(A) == 0
            True
            sage: B == A              # A is not modified
            True

            sage: min_p = A.minimal_polynomial(proof=True)
            sage: min_p.divides(char_p)
            True
        """
        verbose('_charpoly_linbox...')

        if self._nrows != self._ncols:
            raise ValueError("matrix must be square")
        R = self._base_ring[var]
        # call linbox for charpoly
        v = linbox_charpoly(self.p, self._nrows, self._entries)
        r = R(v)
        return r

    def echelonize(self, algorithm='linbox_noefd', **kwds):
        """
        Put ``self`` in reduced row echelon form.

        INPUT:

        - ``self`` -- a mutable matrix

        - ``algorithm``

          - ``linbox`` -- uses the LinBox library (wrapping fflas-ffpack)

          - ``linbox_noefd`` -- uses the FFPACK directly, less memory and faster (default)

          - ``gauss`` -- uses a custom slower `O(n^3)` Gauss
            elimination implemented in Sage

          - ``all`` -- compute using both algorithms and verify that
            the results are the same

        - ``**kwds`` -- these are all ignored

        OUTPUT: ``self`` is put in reduced row echelon form

        - the rank of ``self`` is computed and cached

        - the pivot columns of ``self`` are computed and cached

        - the fact that ``self`` is now in echelon form is recorded and
          cached so future calls to echelonize return immediately

        EXAMPLES::

            sage: A = random_matrix(GF(7), 10, 20)
            sage: E = A.echelon_form()
            sage: A.row_space() == E.row_space()
            True
            sage: all(r[r.nonzero_positions()[0]] == 1 for r in E.rows() if r)
            True

        ::

            sage: A = random_matrix(GF(13), 10, 10)
            sage: while A.rank() != 10:
            ....:     A = random_matrix(GF(13), 10, 10)
            sage: MS = parent(A)
            sage: B = A.augment(MS(1))
            sage: B.echelonize()
            sage: A.rank()
            10
            sage: C = B.submatrix(0,10,10,10)
            sage: ~A == C
            True

        ::

            sage: A = random_matrix(Integers(10), 10, 20)
            sage: A.echelon_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 10'.

        ::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(16007), 10, 20)
            sage: E = A.echelon_form()
            sage: A.row_space() == E.row_space()
            True
            sage: all(r[r.nonzero_positions()[0]] == 1 for r in E.rows() if r)
            True

        ::

            sage: A = random_matrix(Integers(10000), 10, 20)
            sage: A.echelon_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 10000'.

        Parallel computation::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(65521),100,200)
            sage: Parallelism().set('linbox', nproc=2)
            sage: E = A.echelon_form()
            sage: Parallelism().set('linbox', nproc=1) # switch off parallelization
            sage: F = A.echelon_form()
            sage: E==F
            True

        TESTS::

            sage: A = random_matrix(GF(7),  0, 10)
            sage: A.echelon_form()
            []
            sage: A = random_matrix(GF(7), 10,  0)
            sage: A.echelon_form()
            []
            sage: A = random_matrix(GF(7),  0,  0)
            sage: A.echelon_form()
            []
            sage: A = matrix(GF(7),  10,  10)
            sage: A.echelon_form()
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(16007),  0, 10)
            sage: A.echelon_form()
            []
            sage: A = random_matrix(GF(16007), 10,  0)
            sage: A.echelon_form()
            []
            sage: A = random_matrix(GF(16007),  0,  0)
            sage: A.echelon_form()
            []
            sage: A = matrix(GF(16007),  10,  10)
            sage: A.echelon_form()
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]

            sage: A = matrix(GF(97),3,4,range(12))
            sage: A.echelonize(); A
            [ 1  0 96 95]
            [ 0  1  2  3]
            [ 0  0  0  0]
            sage: A.pivots()
            (0, 1)

            sage: for p in (3,17,97,127,1048573):
            ....:    for i in range(10):
            ....:        A = random_matrix(GF(3), 100, 100)
            ....:        A.echelonize(algorithm='all')
        """
        x = self.fetch('in_echelon_form')
        if x is not None:
            return  # already known to be in echelon form

        if not self.base_ring().is_field():
            raise NotImplementedError("Echelon form not implemented over '%s'." % self.base_ring())

        if algorithm == 'linbox':
            self._echelonize_linbox(efd=True)
        elif algorithm == 'linbox_noefd':
            self._echelonize_linbox(efd=False)
        elif algorithm == 'gauss':
            self._echelon_in_place_classical()

        elif algorithm == 'all':
            A = self.__copy__()
            B = self.__copy__()
            self._echelonize_linbox(efd=True)
            A._echelon_in_place_classical()
            B._echelonize_linbox(efd=False)
            if A != self or A != B:
                raise ArithmeticError("Bug in echelon form.")
        else:
            raise ValueError("Algorithm '%s' not known" % algorithm)

    def _echelonize_linbox(self, efd=True):
        """
        Puts ``self`` in row echelon form using LinBox.

        This function is called by echelonize if
        ``algorithm='linbox'``.

        INPUT:

        - ``efd`` -- if ``True`` LinBox's ``EchelonFormDomain``
          implementation is used, which is faster than the direct
          ``LinBox::FFPACK`` implementation, since the latter also
          computes the transformation matrix (which we
          ignore). However, ``efd=True`` uses more memory than FFLAS
          directly (default=``True``)

        EXAMPLES::

            sage: A = random_matrix(GF(7), 10, 20)
            sage: B = copy(A)
            sage: A._echelonize_linbox()
            sage: A.row_space() == B.row_space()
            True
            sage: all(r[r.nonzero_positions()[0]] == 1 for r in A.rows() if r)
            True
        """
        self.check_mutability()
        self.clear_cache()

        t = verbose('Calling echelonize mod %d.' % self.p)
        if efd:
            r, pivots = linbox_echelonize_efd(self.p, self._entries,
                                              self._nrows, self._ncols)
        else:
            r, pivots = linbox_echelonize(self.p, self._entries,
                                          self._nrows, self._ncols)
        verbose('done with echelonize', t)
        self.cache('in_echelon_form', True)
        self.cache('rank', r)
        self.cache('pivots', tuple(pivots))

    def _echelon_in_place_classical(self):
        """
        Puts ``self`` in row echelon form using LinBox.

        This function is called by echelonize if
        ``algorithm='gauss'``.

        EXAMPLES::

            sage: A = random_matrix(GF(7), 10, 20)
            sage: B = copy(A)
            sage: A._echelon_in_place_classical()
            sage: A.row_space() == B.row_space()
            True
            sage: all(r[r.nonzero_positions()[0]] == 1 for r in A.rows() if r)
            True
        """
        self.check_mutability()
        self.clear_cache()

        cdef Py_ssize_t start_row, c, r, nr, nc, i
        cdef celement p, a, s, t, b
        cdef celement **m

        start_row = 0
        p = self.p
        m = self._matrix
        nr = self._nrows
        nc = self._ncols
        pivots = []
        fifth = self._ncols / 10 + 1
        do_verb = (get_verbose() >= 2)
        for c from 0 <= c < nc:
            if do_verb and (c % fifth == 0 and c > 0):
                tm = verbose('on column %s of %s' % (c, self._ncols),
                             level = 2,
                             caller_name = 'matrix_modn_dense echelon')
            # end if
            sig_check()
            for r from start_row <= r < nr:
                a = m[r][c]
                if a:
                    pivots.append(c)
                    a_inverse = celement_invert(a, p)
                    self.rescale_row_c(r, a_inverse, c)
                    self.swap_rows_c(r, start_row)
                    for i from 0 <= i < nr:
                        if i != start_row:
                            b = m[i][c]
                            if b != 0:
                                self.add_multiple_of_row_c(i, start_row, p-b, c)
                    start_row = start_row + 1
                    break
        self.cache('pivots', tuple(pivots))
        self.cache('in_echelon_form', True)

    def right_kernel_matrix(self, algorithm='linbox', basis='echelon'):
        r"""
        Return a matrix whose rows form a basis for the right kernel
        of ``self``.

        If the base ring is the ring of integers modulo a composite,
        the keyword arguments are ignored and the computation is
        delegated to :meth:`Matrix_dense.right_kernel_matrix`.

        INPUT:

        - ``algorithm`` -- (default: ``'linbox'``) a parameter that is
          passed on to ``self.echelon_form``, if computation of an echelon
          form is required; see that routine for allowable values

        - ``basis`` -- (default: ``'echelon'``) a keyword that describes the
          format of the basis returned, allowable values are:

          - ``'echelon'``: the basis matrix is in echelon form
          - ``'pivot'``: the basis matrix is such that the submatrix obtained
             by taking the columns that in ``self`` contain no pivots, is the
             identity matrix
          - ``'computed'``: no work is done to transform the basis; in
             the current implementation the result is the negative of
             that returned by ``'pivot'``

        OUTPUT:

        A matrix ``X`` whose rows are a basis for the right kernel of
        ``self``. This means that ``self * X.transpose()`` is a zero matrix.

        The result is not cached, but the routine benefits when ``self`` is
        known to be in echelon form already.

        EXAMPLES::

            sage: M = matrix(GF(5),6,6,range(36))
            sage: M.right_kernel_matrix(basis='computed')
            [4 2 4 0 0 0]
            [3 3 0 4 0 0]
            [2 4 0 0 4 0]
            [1 0 0 0 0 4]
            sage: M.right_kernel_matrix(basis='pivot')
            [1 3 1 0 0 0]
            [2 2 0 1 0 0]
            [3 1 0 0 1 0]
            [4 0 0 0 0 1]
            sage: M.right_kernel_matrix()
            [1 0 0 0 0 4]
            [0 1 0 0 1 3]
            [0 0 1 0 2 2]
            [0 0 0 1 3 1]
            sage: M * M.right_kernel_matrix().transpose()
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
        """
        if self.fetch('in_echelon_form') is None:
            try:
                self = self.echelon_form(algorithm=algorithm)
            except NotImplementedError:  # composite modulus
                return Matrix_dense.right_kernel_matrix(self)

        cdef Py_ssize_t r = self.rank()
        cdef Py_ssize_t nrows = self._nrows
        cdef Py_ssize_t ncols = self._ncols
        cdef Py_ssize_t i, j, k

        cdef Py_ssize_t* nonpivots = <Py_ssize_t*>sig_malloc(sizeof(Py_ssize_t)*(ncols-r))
        cdef Py_ssize_t* pivots = <Py_ssize_t*>sig_malloc(sizeof(Py_ssize_t)*(r))
        cdef tuple pivot_tuple = self.pivots()

        for i in range(r):
            pivots[i] = pivot_tuple[i]
        j = 0
        k = 0
        for i in range(ncols):
            if j < r and i == pivots[j]:
                j += 1
            else:
                nonpivots[k] = i
                k += 1

        cdef Matrix_modn_dense_template M = self.new_matrix(nrows=ncols-r, ncols=ncols)
        cdef celement pm1 = self.p - 1

        k = 0
        for i in range(ncols-r):
            for j in range(ncols-r):
                M._entries[nonpivots[i]+j*ncols] = 0
            M._entries[nonpivots[i]+k*ncols] = pm1
            k += 1
            for j in range(r):
                M._entries[i*ncols+pivots[j]] = self._entries[nonpivots[i]+j*ncols]

        sig_free(pivots)
        sig_free(nonpivots)
        if basis == 'computed':
            return M
        elif basis == 'pivot':
            return -M
        elif basis != 'echelon':
            raise ValueError("matrix kernel basis format not recognized")
        M.echelonize(algorithm=algorithm)
        return M

    def hessenbergize(self):
        """
        Transform ``self`` in place to its Hessenberg form.

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 10, density=0.1)
            sage: B = copy(A)
            sage: A.hessenbergize()
            sage: all(A[i,j] == 0 for j in range(10) for i in range(j+2, 10))
            True
            sage: A.charpoly() == B.charpoly()
            True
        """
        self.check_mutability()
        x = self.fetch('in_hessenberg_form')
        if x is not None and x:
            return  # already known to be in Hessenberg form
        if self._nrows != self._ncols:
            raise ArithmeticError("Matrix must be square to compute Hessenberg form.")

        cdef Py_ssize_t n
        n = self._nrows

        cdef celement **h
        h = self._matrix

        cdef celement p, t, t_inv, u
        cdef Py_ssize_t i, j, m
        p = self.p

        sig_on()
        for m from 1 <= m < n-1:
            # Search for a nonzero entry in column m-1
            i = -1
            for r from m+1 <= r < n:
                if h[r][m-1]:
                    i = r
                    break

            if i != -1:
                # Found a nonzero entry in column m-1 that is strictly
                # below row m.  Now set i to be the first nonzero position >=
                # m in column m-1.
                if h[m][m-1]:
                    i = m
                t = h[i][m-1]
                t_inv = celement_invert(t, p)
                if i > m:
                    self.swap_rows_c(i, m)
                    self.swap_columns_c(i, m)

                # Now the nonzero entry in position (m,m-1) is t.
                # Use t to clear the entries in column m-1 below m.
                for j from m+1 <= j < n:
                    if h[j][m-1]:
                        u = (h[j][m-1] * t_inv) % p
                        self.add_multiple_of_row_c(j, m, p - u, 0)  # h[j] -= u*h[m]
                        # To maintain charpoly, do the corresponding
                        # column operation, which doesn't mess up the
                        # matrix, since it only changes column m, and
                        # we're only worried about column m-1 right
                        # now.  Add u*column_j to column_m.
                        self.add_multiple_of_column_c(m, j, u, 0)
                # end for
            # end if
        # end for
        sig_off()
        self.cache('in_hessenberg_form', True)

    def _charpoly_hessenberg(self, var):
        """
        Transform ``self`` in place to its Hessenberg form then computes
        and returns the coefficients of the characteristic polynomial
        of this matrix.

        INPUT:

        - ``var`` -- name of the indeterminate of the charpoly

        OUTPUT:

           The characteristic polynomial is represented as a vector of
           ints, where the constant term of the characteristic
           polynomial is the `0`-th coefficient of the vector.

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 10, density=0.1)
            sage: B = copy(A)
            sage: P.<x> = GF(17)[]
            sage: A._charpoly_hessenberg('x') == B.charpoly()
            True
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("charpoly not defined for non-square matrix.")

        cdef Py_ssize_t i, m, n,
        n = self._nrows

        cdef celement p, t
        p = self.p

        # Replace self by its Hessenberg form, and set H to this form
        # for notation below.
        cdef Matrix_modn_dense_template H
        H = self.__copy__()
        H.hessenbergize()

        # We represent the intermediate polynomials that come up in
        # the calculations as rows of an (n+1)x(n+1) matrix, since
        # we've implemented basic arithmetic with such a matrix.
        # Please see the generic implementation of charpoly in
        # matrix.py to see more clearly how the following algorithm
        # actually works.  (The implementation is clearer (but slower)
        # if one uses polynomials to represent polynomials instead of
        # using the rows of a matrix.)  Also see Cohen's first GTM,
        # Algorithm 2.2.9.

        cdef Matrix_modn_dense_template c
        c = self.new_matrix(nrows=n+1, ncols=n+1)    # the 0 matrix
        c._matrix[0][0] = 1
        for m from 1 <= m <= n:
            # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
            # We do this by hand by setting the m-th row to c[m-1]
            # shifted to the right by one.  We then add
            # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
            for i from 1 <= i <= n:
                c._matrix[m][i] = c._matrix[m-1][i-1]
            # the p-.. below is to keep scalar normalized between 0 and p.
            c.add_multiple_of_row_c(m, m-1, p - H._matrix[m-1][m-1], 0)
            t = 1
            for i from 1 <= i < m:
                t = (t*H._matrix[m-i][m-i-1]) % p
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c.add_multiple_of_row_c(m, m-i-1, p - (t*H._matrix[m-i-1][m-1]) % p, 0)

        # The answer is now the n-th row of c.
        v = []
        for i from 0 <= i <= n:
            v.append(int(c._matrix[n][i]))
        R = self._base_ring[var]    # polynomial ring over the base ring
        return R(v)

    def rank(self):
        """
        Return the rank of this matrix.

        EXAMPLES::

            sage: A = random_matrix(GF(3), 100, 100)
            sage: B = copy(A)
            sage: _ = A.rank()
            sage: B == A
            True

            sage: A = random_matrix(GF(3), 100, 100, density=0.01)
            sage: A.transpose().rank() == A.rank()
            True

            sage: A = matrix(GF(3), 100, 100)
            sage: A.rank()
            0

        Rank is not implemented over the integers modulo a composite
        yet.::

            sage: M = matrix(Integers(4), 2, [2,2,2,2])
            sage: M.rank()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 4'.

        ::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(16007), 100, 100)
            sage: B = copy(A)
            sage: A.rank()
            100
            sage: B == A
            True
            sage: MS = A.parent()
            sage: MS(1) == ~A*A
            True

        TESTS::

            sage: A = random_matrix(GF(7), 0, 0)
            sage: A.rank()
            0
            sage: A = random_matrix(GF(7), 1, 0)
            sage: A.rank()
            0
            sage: A = random_matrix(GF(7), 0, 1)
            sage: A.rank()
            0

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(16007), 0, 0)
            sage: A.rank()
            0
            sage: A = random_matrix(GF(16007), 1, 0)
            sage: A.rank()
            0
            sage: A = random_matrix(GF(16007), 0, 1)
            sage: A.rank()
            0
        """
        cdef Matrix_modn_dense_template A
        if self.p > 2 and is_prime(self.p):
            x = self.fetch('rank')
            if x is not None:
                return x
            r = Integer(linbox_rank(self.p, self._entries, self._nrows, self._ncols))
            self.cache('rank', r)
            return r
        else:
            # linbox is very buggy for p=2, but this code should never
            # be called since p=2 is handled via M4RI
            return Matrix_dense.rank(self)

    def determinant(self):
        """
        Return the determinant of this matrix.

        EXAMPLES::

            sage: s = set()
            sage: while s != set(GF(7)):
            ....:     A = random_matrix(GF(7), 10, 10)
            ....:     s.add(A.determinant())

        ::

            sage: A = random_matrix(GF(7), 100, 100)
            sage: A.determinant() == A.transpose().determinant()
            True

            sage: B = random_matrix(GF(7), 100, 100)
            sage: (A*B).determinant() == A.determinant() * B.determinant()
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(16007), 10, 10)
            sage: A.determinant().parent() is GF(16007)
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(16007), 100, 100)
            sage: A.determinant().parent() is GF(16007)
            True
            sage: A.determinant() == A.transpose().determinant()
            True
            sage: B = random_matrix(GF(16007), 100, 100)
            sage: (A*B).determinant() == A.determinant() * B.determinant()
            True

        Parallel computation::

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(65521),200)
            sage: B = copy(A)
            sage: Parallelism().set('linbox', nproc=2)
            sage: d = A.determinant()
            sage: Parallelism().set('linbox', nproc=1) # switch off parallelization
            sage: e = B.determinant()
            sage: d==e
            True

        TESTS::

            sage: A = random_matrix(GF(7), 0, 0); A.det()
            1

            sage: A = random_matrix(GF(7), 0, 1); A.det()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix

            sage: A = random_matrix(GF(7), 1, 0); A.det()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix

            sage: A = matrix(GF(7), 5, 5); A.det()                                      # needs sage.libs.pari
            0

            sage: # needs sage.rings.finite_rings
            sage: A = random_matrix(GF(16007), 0, 0); A.det()
            1
            sage: A = random_matrix(GF(16007), 0, 1); A.det()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
            sage: A = random_matrix(GF(16007), 1, 0); A.det()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
            sage: A = matrix(GF(16007), 5, 5); A.det()
            0
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        if self._nrows == 0:
            return self._coerce_element(1)

        if self.p > 2 and is_prime(self.p):
            x = self.fetch('det')
            if x is not None:
                return x
            d = linbox_det(self.p, self._entries, self._nrows)
            d2 = self._coerce_element(d)
            self.cache('det', d2)
            return d2
        else:
            return Matrix_dense.determinant(self)

    cdef xgcd_eliminate(self, celement * row1, celement* row2, Py_ssize_t start_col):
        r"""
        Reduces ``row1`` and ``row2`` by a unimodular transformation
        using the xgcd relation between their first coefficients ``a`` and
        ``b``.

        INPUT:

        - ``row1``, ``row2`` -- the two rows to be transformed (within
          ``self``)

        - ``start_col`` -- the column of the pivots in ``row1`` and
          ``row2``. It is assumed that all entries before ``start_col``
          in ``row1`` and ``row2`` are zero.

        OUTPUT:

        - g: the gcd of the first elements of row1 and
          row2 at column start_col

        - put row1 = s \* row1 + t \* row2 row2 = w \*
          row1 + t \* row2 where g = sa + tb
        """
        cdef int p = <int>self.p
        cdef celement *row1_p
        cdef celement *row2_p
        cdef celement tmp
        cdef int g, s, t, v, w
        cdef Py_ssize_t nc, i
        cdef int a = <int>row1[start_col]
        cdef int b = <int>row2[start_col]
        g = ArithIntObj.c_xgcd_int(a, b, <int*>&s, <int*>&t)
        v = a/g
        w = -<int>b/g
        nc = self.ncols()

        for i from start_col <= i < nc:
            tmp = (s * <int>row1[i] + t * <int>row2[i]) % p
            row2[i] = (w* <int>row1[i] + v*<int>row2[i]) % p
            row1[i] = tmp
        return g

    cdef rescale_row_c(self, Py_ssize_t row, multiple, Py_ssize_t start_col):
        """
        Rescale ``self[row]`` by ``multiple`` but only start at column
        index ``start_col``.

        INPUT:

        - ``row`` -- integer
        - ``multiple`` -- finite field element
        - ``start_col`` -- integer

        EXAMPLES::

            sage: A = matrix(GF(19), 4, 4, range(16)); A
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [12 13 14 15]

            sage: A.rescale_row(1, 17, 2); A
            [ 0  1  2  3]
            [ 4  5  7  5]
            [ 8  9 10 11]
            [12 13 14 15]

            sage: 6*17 % 19 == A[1,2]
            True

            sage: A = matrix(Integers(2^4), 4, 4, range(16)); A
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [12 13 14 15]

            sage: A.rescale_row(1, 3, 2); A
            [ 0  1  2  3]
            [ 4  5  2  5]
            [ 8  9 10 11]
            [12 13 14 15]

            sage: 6*3 % 16 == A[1,2]
            True
        """
        cdef celement r, p
        cdef celement* v
        cdef Py_ssize_t i
        p = self.p
        v = self._matrix[row]
        for i from start_col <= i < self._ncols:
            v[i] = (v[i]*<celement>multiple) % p

    cdef rescale_col_c(self, Py_ssize_t col, multiple, Py_ssize_t start_row):
        """
        EXAMPLES::

            sage: B = MatrixSpace(Integers(37),3,3)([1]*9)
            sage: B
            [1 1 1]
            [1 1 1]
            [1 1 1]
            sage: B.rescale_col(1,5)
            sage: B
            [1 5 1]
            [1 5 1]
            [1 5 1]

        Recaling need not include the entire row.::

            sage: B.rescale_col(0,2,1); B
            [1 5 1]
            [2 5 1]
            [2 5 1]

        Bounds are checked.::

            sage: B.rescale_col(3,2)
            Traceback (most recent call last):
            ...
            IndexError: matrix column index out of range

        Rescaling by a negative number::

            sage: B.rescale_col(2,-3); B
            [ 1  5 34]
            [ 2  5 34]
            [ 2  5 34]
        """
        cdef celement r, p
        cdef celement* v
        cdef Py_ssize_t i
        p = self.p
        for i from start_row <= i < self._nrows:
            self._matrix[i][col] = (self._matrix[i][col]*<celement>multiple) % p

    cdef add_multiple_of_row_c(self, Py_ssize_t row_to, Py_ssize_t row_from, multiple, Py_ssize_t start_col):
        """
        Add ``multiple`` times ``self[row_from]`` to ``self[row_to]``
        statting in column ``start_col``.

        EXAMPLES::

            sage: A = random_matrix(GF(37), 10, 10)
            sage: B = copy(A)

            sage: A.add_multiple_of_row(2, 3, 10)
            sage: all(A[i] == B[i] for i in range(10) if not i == 2)
            True
            sage: A[2] == B[2] + 10*B[3]
            True

            sage: A.add_multiple_of_row(2, 3, 10, 4)
            sage: all(A[i] == B[i] for i in range(10) if not i == 2)
            True
            sage: A[2][:4] == B[2][:4] + 10*B[3][:4]
            True
            sage: A[2][4:] == B[2][4:] + 20*B[3][4:]
            True
        """
        cdef celement p
        cdef celement *v_from
        cdef celement *v_to

        p = self.p
        v_from = self._matrix[row_from]
        v_to = self._matrix[row_to]

        cdef Py_ssize_t i, nc
        nc = self._ncols
        for i from start_col <= i < nc:
            v_to[i] = ((<celement>multiple) * v_from[i] + v_to[i]) % p

    cdef add_multiple_of_column_c(self, Py_ssize_t col_to, Py_ssize_t col_from, multiple, Py_ssize_t start_row):
        """
        Add ``multiple`` times ``self[row_from]`` to ``self[row_to]``
        statting in column ``start_col``.

        EXAMPLES::

            sage: A = random_matrix(GF(37), 10, 10)
            sage: B = copy(A)

            sage: A.add_multiple_of_column(2, 3, 10)
            sage: all(A.column(i) == B.column(i) for i in range(10) if not i == 2)
            True
            sage: A.column(2) == B.column(2) + 10*B.column(3)
            True

            sage: A.add_multiple_of_column(2, 3, 10, 4)
            sage: all(A.column(i) == B.column(i) for i in range(10) if not i == 2)
            True
            sage: A.column(2)[:4] == B.column(2)[:4] + 10*B.column(3)[:4]
            True
            sage: A.column(2)[4:] == B.column(2)[4:] + 20*B.column(3)[4:]
            True
        """
        cdef celement  p
        cdef celement **m

        m = self._matrix
        p = self.p

        cdef Py_ssize_t i, nr
        nr = self._nrows
        for i from start_row <= i < self._nrows:
            m[i][col_to] = (m[i][col_to] + (<celement>multiple) * m[i][col_from]) %p

    cdef swap_rows_c(self, Py_ssize_t row1, Py_ssize_t row2):
        """
        EXAMPLES::

            sage: A = matrix(Integers(8), 2,[1,2,3,4])
            sage: A.swap_rows(0,1)
            sage: A
            [3 4]
            [1 2]
        """
        cdef celement* r1 = self._matrix[row1]
        cdef celement* r2 = self._matrix[row2]
        cdef celement temp
        for i in range(self._ncols):
            temp = r1[i]
            r1[i] = r2[i]
            r2[i] = temp

    cdef swap_columns_c(self, Py_ssize_t col1, Py_ssize_t col2):
        """
        EXAMPLES::

            sage: A = matrix(Integers(8), 2,[1,2,3,4])
            sage: A.swap_columns(0,1)
            sage: A
            [2 1]
            [4 3]
        """
        cdef Py_ssize_t i, nr
        cdef celement t
        cdef celement **m
        m = self._matrix
        nr = self._nrows
        for i from 0 <= i < nr:
            t = m[i][col1]
            m[i][col1] = m[i][col2]
            m[i][col2] = t

    def randomize(self, density=1, nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this
        matrix, leaving the rest unchanged.

        INPUT:

        - ``density`` -- integer; proportion (roughly) to be considered
          for changes
        - ``nonzero`` -- boolean (default: ``False``); whether the new
          entries are forced to be nonzero

        OUTPUT: none, the matrix is modified in-space

        EXAMPLES::

            sage: A = matrix(GF(5), 5, 5, 0)
            sage: total_count = 0
            sage: from collections import defaultdict
            sage: dic = defaultdict(Integer)
            sage: def add_samples(density):
            ....:     global dic, total_count
            ....:     for _ in range(100):
            ....:         A = Matrix(GF(5), 5, 5, 0)
            ....:         A.randomize(density)
            ....:         for a in A.list():
            ....:             dic[a] += 1
            ....:             total_count += 1.0

            sage: add_samples(1.0)
            sage: while not all(abs(dic[a]/total_count - 1/5) < 0.01 for a in dic):
            ....:     add_samples(1.0)

            sage: def add_sample(density):
            ....:     global density_sum, total_count
            ....:     total_count += 1.0
            ....:     density_sum += random_matrix(GF(5), 1000, 1000, density=density).density()

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.5)
            sage: expected_density = 1.0 - (999/1000)^500
            sage: expected_density
            0.3936...
            sage: while abs(density_sum/total_count - expected_density) > 0.001:
            ....:     add_sample(0.5)

        The matrix is updated instead of overwritten::

            sage: def add_sample(density):
            ....:     global density_sum, total_count
            ....:     total_count += 1.0
            ....:     A = random_matrix(GF(5), 1000, 1000, density=density)
            ....:     A.randomize(density=density, nonzero=True)
            ....:     density_sum += A.density()

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.5)
            sage: expected_density = 1.0 - (999/1000)^1000
            sage: expected_density
            0.6323...
            sage: while abs(density_sum/total_count - expected_density) > 0.001:
            ....:     add_sample(0.5)

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.1)
            sage: expected_density = 1.0 - (999/1000)^200
            sage: expected_density
            0.1813...
            sage: while abs(density_sum/total_count - expected_density) > 0.001:
            ....:     add_sample(0.1)
        """
        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = float(1)

        self.check_mutability()
        self.clear_cache()

        cdef randstate rstate = current_randstate()

        cdef int nc, p = <int>self.p
        cdef long pm1

        if not nonzero:
            # Original code, before adding the ``nonzero`` option.
            if density == 1:
                for i from 0 <= i < self._nrows*self._ncols:
                    self._entries[i] = rstate.c_random() % p
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                sig_on()
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random() % nc
                        self._matrix[i][k] = rstate.c_random() % p
                sig_off()
        else:
            # New code, to implement the ``nonzero`` option.
            pm1 = p - 1
            if density == 1:
                for i from 0 <= i < self._nrows*self._ncols:
                    self._entries[i] = (rstate.c_random() % pm1) + 1
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                sig_on()
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random() % nc
                        self._matrix[i][k] = (rstate.c_random() % pm1) + 1
                sig_off()

    def _magma_init_(self, magma):
        """
        Return a string representation of ``self`` in Magma form.

        INPUT:

        - ``magma`` -- a Magma session

        OUTPUT: string

        EXAMPLES::

            sage: a = matrix(GF(389),2,2,[1..4])
            sage: magma(a)                                         # optional - magma
            [  1   2]
            [  3   4]
            sage: a._magma_init_(magma)                            # optional - magma
            'Matrix(GF(389),2,2,StringToIntegerSequence("1 2 3 4"))'

        A consistency check::

            sage: a = random_matrix(GF(13),50); b = random_matrix(GF(13),50)
            sage: magma(a*b) == magma(a)*magma(b)                  # optional - magma
            True
        """
        s = self.base_ring()._magma_init_(magma)
        return 'Matrix(%s,%s,%s,StringToIntegerSequence("%s"))' % (
            s, self._nrows, self._ncols, self._export_as_string())

    cpdef _export_as_string(self):
        """
        Return space separated string of the entries in this matrix.

        EXAMPLES::

            sage: A = matrix(GF(997),2,3,[1,2,5,-3,8,2]); A
            [  1   2   5]
            [994   8   2]
            sage: A._export_as_string()
            '1 2 5 994 8 2'
        """
        cdef int ndigits = len(str(self.p))

        cdef Py_ssize_t i, n
        cdef char *s
        cdef char *t

        if self._nrows == 0 or self._ncols == 0:
            data = ''
        else:
            n = self._nrows*self._ncols*(ndigits + 1) + 2  # spaces between each number plus trailing null
            s = <char*>check_malloc(n * sizeof(char))
            t = s
            sig_on()
            for i in range(self._nrows * self._ncols):
                t += snprintf(t, ndigits+2, "%ld ", <long>self._entries[i])

            sig_off()
            data = char_to_str(s)[:-1]
            sig_free(s)
        return data

    def _list(self):
        """
        Return list of elements of ``self``.  This method is called by ``self.list()``.

        EXAMPLES::

            sage: w = matrix(GF(19), 2, 3, [1..6])
            sage: w.list()
            [1, 2, 3, 4, 5, 6]
            sage: w._list()
            [1, 2, 3, 4, 5, 6]
            sage: w.list()[0].parent()
            Finite Field of size 19

        TESTS::

            sage: w = random_matrix(GF(3),100)
            sage: w.parent()(w.list()) == w
            True
        """
        cdef Py_ssize_t i
        F = self.base_ring()
        return [F(<int>self._entries[i]) for i in range(self._nrows*self._ncols)]

    def lift(self):
        """
        Return the lift of this matrix to the integers.

        EXAMPLES::

            sage: A = matrix(GF(7),2,3,[1..6])
            sage: A.lift()
            [1 2 3]
            [4 5 6]
            sage: A.lift().parent()
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

            sage: # needs sage.rings.finite_rings
            sage: A = matrix(GF(16007),2,3,[1..6])
            sage: A.lift()
            [1 2 3]
            [4 5 6]
            sage: A.lift().parent()
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

        Subdivisions are preserved when lifting::

            sage: A.subdivide([], [1,1]); A
            [1||2 3]
            [4||5 6]
            sage: A.lift()
            [1||2 3]
            [4||5 6]
        """
        cdef Py_ssize_t i, j

        cdef Matrix_integer_dense L
        cdef object P = matrix_space.MatrixSpace(ZZ, self._nrows,
                                                 self._ncols, sparse=False)
        L = Matrix_integer_dense(P, ZZ(0), False, False)
        cdef celement* A_row
        for i in range(self._nrows):
            A_row = self._matrix[i]
            for j in range(self._ncols):
                L.set_unsafe_double(i, j, A_row[j])
        if self._subdivisions is not None:
            L.subdivide(*self.subdivisions())
        return L

    def transpose(self):
        """
        Return the transpose of ``self``, without changing ``self``.

        EXAMPLES:

        We create a matrix, compute its transpose, and note that
        the original matrix is not changed. ::

            sage: M = MatrixSpace(GF(41),  2)
            sage: A = M([1,2,3,4])
            sage: B = A.transpose()
            sage: B
            [1 3]
            [2 4]
            sage: A
            [1 2]
            [3 4]

        ``.T`` is a convenient shortcut for the transpose::

           sage: A.T
           [1 3]
           [2 4]

        ::

            sage: A.subdivide(None, 1); A
            [1|2]
            [3|4]
            sage: A.transpose()
            [1 3]
            [---]
            [2 4]
        """
        cdef Py_ssize_t nrows = self._nrows
        cdef Py_ssize_t ncols = self._ncols

        cdef Matrix_modn_dense_template M = self.new_matrix(nrows=ncols, ncols=nrows)
        cdef Py_ssize_t i, j

        for i from 0 <= i < ncols:
            for j from 0 <= j < nrows:
                M._entries[j+i*nrows] = self._entries[i+j*ncols]

        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            M.subdivide(col_divs, row_divs)

        return M

    cdef _stack_impl(self, bottom):
        r"""
        Implementation of :meth:`stack` by returning a new matrix
        formed by appending the matrix ``bottom`` beneath ``self``.

        Assume that ``self`` and ``other`` are compatible in the sense
        that they have the same base ring and that both are dense.

        INPUT:

        - ``bottom`` -- a matrix compatible with ``self``

        EXAMPLES:

        Stacking with a matrix::

            sage: A = matrix(GF(41), 4, 3, range(12))
            sage: B = matrix(GF(41), 3, 3, range(9))
            sage: A.stack(B)
            [ 0  1  2]
            [ 3  4  5]
            [ 6  7  8]
            [ 9 10 11]
            [ 0  1  2]
            [ 3  4  5]
            [ 6  7  8]

        Stacking with a vector::

            sage: A = matrix(GF(41), 3, 2, [0, 2, 4, 6, 8, 10])
            sage: v = vector(GF(41), 2, [100, 200])
            sage: A.stack(v)
            [ 0  2]
            [ 4  6]
            [ 8 10]
            [18 36]

        Errors are raised if the sizes are incompatible::

            sage: A = matrix(GF(41), [[1, 2],[3, 4]])
            sage: B = matrix(GF(41), [[10, 20, 30], [40, 50, 60]])
            sage: A.stack(B)
            Traceback (most recent call last):
            ...
            TypeError: number of columns must be the same, not 2 and 3

            sage: v = vector(GF(41), [100, 200, 300])
            sage: A.stack(v)
            Traceback (most recent call last):
            ...
            TypeError: number of columns must be the same, not 2 and 3

        Setting ``subdivide`` to ``True`` will, in its simplest form,
        add a subdivision between ``self`` and ``bottom``::

            sage: A = matrix(GF(41), 2, 5, range(10))
            sage: B = matrix(GF(41), 3, 5, range(15))
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [--------------]
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [10 11 12 13 14]

        Row subdivisions are preserved by stacking, and enriched,
        if subdivisions are requested.  (So multiple stackings can
        be recorded.) ::

            sage: A = matrix(GF(41), 2, 4, range(8))
            sage: A.subdivide([1], None)
            sage: B = matrix(GF(41), 3, 4, range(12))
            sage: B.subdivide([2], None)
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3]
            [-----------]
            [ 4  5  6  7]
            [-----------]
            [ 0  1  2  3]
            [ 4  5  6  7]
            [-----------]
            [ 8  9 10 11]

        Column subdivisions can be preserved, but only if they are identical.
        Otherwise, this information is discarded and must be managed
        separately. ::

            sage: A = matrix(GF(41), 2, 5, range(10))
            sage: A.subdivide(None, [2,4])
            sage: B = matrix(GF(41), 3, 5, range(15))
            sage: B.subdivide(None, [2,4])
            sage: A.stack(B, subdivide=True)
            [ 0  1| 2  3| 4]
            [ 5  6| 7  8| 9]
            [-----+-----+--]
            [ 0  1| 2  3| 4]
            [ 5  6| 7  8| 9]
            [10 11|12 13|14]

            sage: A.subdivide(None, [1,2])
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [--------------]
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [10 11 12 13 14]

        The result retains the base ring of ``self`` by coercing
        the elements of ``bottom`` into the base ring of ``self``::

            sage: A = matrix(GF(41), 1, 2, [1,2])
            sage: B = matrix(ZZ, 1, 2, [100, 100])
            sage: C = A.stack(B); C
            [ 1  2]
            [18 18]

            sage: C.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 41

            sage: D = B.stack(A); D
            [18 18]
            [ 1  2]

            sage: D.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 41
        """
        cdef Matrix_modn_dense_template other = <Matrix_modn_dense_template> bottom
        cdef Matrix_modn_dense_template M = self.new_matrix(nrows=self._nrows+other._nrows,
                                                            ncols=self._ncols)
        cdef Py_ssize_t selfsize = self._ncols * self._nrows
        memcpy(M._entries, self._entries, sizeof(celement)*selfsize)
        memcpy(M._entries+selfsize, other._entries, sizeof(celement)*other._ncols*other._nrows)
        return M

    def submatrix(self,
                  Py_ssize_t row=0, Py_ssize_t col=0,
                  Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        r"""
        Return the matrix constructed from ``self`` using the specified
        range of rows and columns.

        INPUT:

        - ``row``, ``col`` -- index of the starting row and column;
          indices start at zero

        - ``nrows``, ``ncols`` -- (optional) number of rows and columns to
          take; if not provided, take all rows below and all columns to
          the right of the starting entry

        .. SEEALSO::

            The functions :func:`matrix_from_rows`,
            :func:`matrix_from_columns`, and
            :func:`matrix_from_rows_and_columns` allow one to select
            arbitrary subsets of rows and/or columns.

        EXAMPLES:

        Take the `3 \times 3` submatrix starting from entry `(1,1)` in a
        `4 \times 4` matrix::

            sage: m = matrix(GF(17),4, [1..16])
            sage: m.submatrix(1, 1)
            [ 6  7  8]
            [10 11 12]
            [14 15 16]

        Same thing, except take only two rows::

            sage: m.submatrix(1, 1, 2)
            [ 6  7  8]
            [10 11 12]

        And now take only one column::

            sage: m.submatrix(1, 1, 2, 1)
            [ 6]
            [10]

        You can take zero rows or columns if you want::

            sage: m.submatrix(0, 0, 0)
            []
            sage: parent(m.submatrix(0, 0, 0))
            Full MatrixSpace of 0 by 4 dense matrices over Finite Field of size 17
        """
        if ncols == -1:
            ncols = self._ncols - col

        if nrows == -1:
            nrows = self._nrows - row

        if nrows < 0 or row < 0 or row + nrows > self._nrows:
            raise IndexError("rows out of range")
        if ncols < 0 or col < 0 or col + ncols > self._ncols:
            raise IndexError("columns out of range")

        cdef Matrix_modn_dense_template M = self.new_matrix(nrows=nrows, ncols=ncols)

        if col == 0 and ncols == self._ncols:
            memcpy(M._entries, self._matrix[row], sizeof(celement)*ncols*nrows)
            return M

        cdef Py_ssize_t i, r
        for i, r in enumerate(range(row, row+nrows)) :
            memcpy(M._matrix[i], self._matrix[r]+col, sizeof(celement)*ncols)

        return M

    def _matrices_from_rows(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Make a list of matrix from the rows of this matrix.  This is a
        fairly technical function which is used internally, e.g., by
        the cyclotomic field linear algebra code.

        INPUT:

        - ``nrows`` -- integer

        - ``ncols`` -- integer

        EXAMPLES::

            sage: A = matrix(GF(127), 4, 4, range(16))
            sage: A
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [12 13 14 15]
            sage: A._matrices_from_rows(2,2)
            [
            [0 1]  [4 5]  [ 8  9]  [12 13]
            [2 3], [6 7], [10 11], [14 15]
            ]

        OUTPUT: list of matrices
        """
        if nrows * ncols != self._ncols:
            raise ValueError("nrows * ncols must equal self's number of columns")

        cdef Matrix_modn_dense_template M
        cdef Py_ssize_t i
        cdef Py_ssize_t n = nrows * ncols
        ans = []
        for i from 0 <= i < self._nrows:
            M = self.new_matrix(nrows = nrows, ncols = ncols)
            memcpy(M._entries, self._entries+i*n, sizeof(celement)*n)
            ans.append(M)
        return ans

    def matrix_from_columns(self, columns):
        """
        Return the matrix constructed from ``self`` using columns with indices
        in the columns list.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_columns([2,1])
            [2 1]
            [5 4]
            [0 7]
        """
        cdef Py_ssize_t ncols = len(columns)

        # Construct new matrix
        cdef Matrix_modn_dense_template A = self.new_matrix(ncols=ncols)
        cdef Py_ssize_t i, j, col
        for j, col in enumerate(columns):
            if col < 0 or col >= self._ncols:
                raise IndexError("column index out of range")
            for i in range(self._nrows):
                A._matrix[i][j] = self._matrix[i][col]

        return A

    def matrix_from_rows(self, rows):
        """
        Return the matrix constructed from ``self`` using rows with indices in
        the rows list.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows([2,1])
            [6 7 0]
            [3 4 5]
        """
        cdef Py_ssize_t nrows = len(rows)

        # Construct new matrix
        cdef Matrix_modn_dense_template A = self.new_matrix(nrows=nrows)

        cdef Py_ssize_t i, row
        for i, row in enumerate(rows):
            if row < 0 or row >= self._nrows:
                raise IndexError("row index out of range")
            memcpy(A._matrix[i], self._matrix[row], sizeof(celement)*self._ncols)

        return A

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from ``self`` from the given rows and
        columns.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows_and_columns([1], [0,2])
            [3 5]
            sage: A.matrix_from_rows_and_columns([1,2], [1,2])
            [4 5]
            [7 0]

        Note that row and column indices can be reordered or repeated::

            sage: A.matrix_from_rows_and_columns([2,1], [2,1])
            [0 7]
            [5 4]

        For example here we take from row 1 columns 2 then 0 twice, and do
        this 3 times::

            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0,0])
            [5 3 3]
            [5 3 3]
            [5 3 3]

        AUTHORS:

        - Jaap Spies (2006-02-18)

        - Didier Deshommes: some Pyrex speedups implemented
        """
        cdef Py_ssize_t ncols = len(columns)
        cdef Py_ssize_t nrows = len(rows)

        # Check whether column indices are valid
        cdef Py_ssize_t i, j, row, col
        for col in columns:
            if col < 0 or col >= self._ncols:
                raise IndexError("column index out of range")

        # Construct new matrix
        cdef Matrix_modn_dense_template A = self.new_matrix(nrows=nrows, ncols=ncols)
        for i, row in enumerate(rows):
            if row < 0 or row >= self._nrows:
                raise IndexError("row index out of range")
            for j, col in enumerate(columns):
                A._matrix[i][j] = self._matrix[row][col]

        return A

    def __bool__(self):
        """
        Test whether this matrix is zero.

        EXAMPLES::

            sage: A = matrix(GF(7), 10, 10, range(100))
            sage: A == 0 # indirect doctest
            False
            sage: A.is_zero()
            False

            sage: A = matrix(Integers(10), 10, 10)
            sage: bool(A)
            False

            sage: # needs sage.rings.finite_rings
            sage: A = matrix(GF(16007), 0, 0)
            sage: A.is_zero()
            True
            sage: A = matrix(GF(16007), 1, 0)
            sage: A.is_zero()
            True
            sage: A = matrix(GF(16007), 0, 1)
            sage: A.is_zero()
            True
        """
        return not linbox_is_zero(self.p, self._entries, self._nrows, self._ncols)

    _matrix_from_rows_of_matrices = staticmethod(__matrix_from_rows_of_matrices)

    cdef int _copy_row_to_mod_int_array(self, mod_int *to, Py_ssize_t i) noexcept:
        cdef Py_ssize_t j
        cdef celement *_from = self._entries+(i*self._ncols)
        for j in range(self._ncols):
            to[j] = <mod_int>_from[j]

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        r"""
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: M = Matrix(GF(49), 2, [1,2,-2,0])                                     # needs sage.rings.finite_rings
            sage: M.zero_pattern_matrix()  # indirect doctest                           # needs sage.rings.finite_rings
            [0 0]
            [0 1]

            sage: M = Matrix(Integers(10), 2, [1,2,-2,0])
            sage: M.zero_pattern_matrix()  # indirect doctest
            [0 0]
            [0 1]
        """
        return self._entries[j+i*self._ncols] == 0
