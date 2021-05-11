"""
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
    sage: x = A\v
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

#*****************************************************************************
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
#       Copyright (C) 2011 Rob Beezer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.stdint cimport uint64_t
from cpython.bytes cimport *

from cysignals.memory cimport check_malloc, check_allocarray, sig_malloc, sig_free
from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.libs.gmp.mpz cimport *
from sage.libs.linbox.fflas cimport FFLAS_TRANSPOSE, FflasNoTrans, FflasTrans, \
    FflasRight, vector, list as std_list
from libcpp cimport bool
from sage.parallel.parallelism import Parallelism

cimport sage.rings.fast_arith
cdef sage.rings.fast_arith.arith_int ArithIntObj
ArithIntObj  = sage.rings.fast_arith.arith_int()

# for copying/pickling
from libc.string cimport memcpy
from libc.stdio cimport snprintf

from sage.modules.vector_modn_dense cimport Vector_modn_dense

from sage.arith.all import is_prime
from sage.structure.element cimport (Element, Vector, Matrix,
        ModuleElement, RingElement)
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract
from sage.misc.misc import cputime
from sage.misc.verbose import verbose, get_verbose
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.structure.proof.proof import get_flag as get_proof_flag
from sage.structure.richcmp cimport rich_to_bool
from sage.misc.randstate cimport randstate, current_randstate
import sage.matrix.matrix_space as matrix_space
from .args cimport MatrixArgs_init


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
        q = x_int / y_int # integer quotient
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
        return 0,[]

    cdef Py_ssize_t i, j
    cdef ModField *F = new ModField(<long>modulus)
    cdef size_t* P = <size_t*>check_allocarray(nrows, sizeof(size_t))
    cdef size_t* Q = <size_t*>check_allocarray(ncols, sizeof(size_t))

    cdef Py_ssize_t r
    cpdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')
    cdef bool transform = False
    if nrows*ncols > 1000: sig_on()
    if nbthreads > 1 :
        r = pReducedRowEchelonForm(F[0], nrows, ncols, <ModField.Element*>entries, ncols, P, Q, transform, nbthreads)
    else :
        r = ReducedRowEchelonForm(F[0], nrows, ncols, <ModField.Element*>entries, ncols, P, Q)
    if nrows*ncols > 1000: sig_off()

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
        return 0,[]

    cdef ModField *F = new ModField(<long>modulus)
    cdef DenseMatrix *A = new DenseMatrix(F[0], <ModField.Element*>entries,<Py_ssize_t>nrows, <Py_ssize_t>ncols)
    cdef Py_ssize_t r = reducedRowEchelonize(A[0])
    cdef Py_ssize_t i,j
    for i in range(nrows):
        for j in range(ncols):
            entries[i*ncols+j] = <celement>A.getEntry(i,j)

    cdef Py_ssize_t ii = 0
    cdef list pivots = []
    for i in range(r):
        for j in range(ii,ncols):
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

    cpdef Py_ssize_t r
    cpdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')
    if nrows*ncols > 1000: sig_on()
    if nbthreads > 1:
        r = pRank(F[0], nrows, ncols, <ModField.Element*>cpy, ncols, nbthreads)
    else:
        r = Rank(F[0], nrows, ncols, <ModField.Element*>cpy, ncols)
    if nrows*ncols > 1000: sig_off()
    sig_free(cpy)
    del F
    return r

cdef inline celement linbox_det(celement modulus, celement* entries, Py_ssize_t n):
    """
    Return the determinant of this matrix.
    """
    cdef ModField *F = new ModField(<long>modulus)
    cdef celement *cpy = linbox_copy(modulus, entries, n, n)

    cdef celement d
    cpdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')

    if n*n > 1000: sig_on()
    if nbthreads > 1 :
        pDet(F[0], d, n, <ModField.Element*>cpy, n, nbthreads)
    else :
        Det(F[0], d, n, <ModField.Element*>cpy, n)
    if n*n > 1000: sig_off()
    sig_free(cpy)
    del F
    return d

cdef inline celement linbox_matrix_matrix_multiply(celement modulus, celement* ans, celement* A, celement* B, Py_ssize_t m, Py_ssize_t n, Py_ssize_t k) :
    """
    C = A*B
    """
    cdef ModField *F = new ModField(<long>modulus)
    cdef ModField.Element one, zero
    F[0].init(one, <int>1)
    F[0].init(zero, <int>0)

    cpdef size_t nbthreads
    nbthreads = Parallelism().get('linbox')

    if m*n*k > 100000: sig_on()
    if nbthreads > 1 :
        pfgemm(F[0], FflasNoTrans, FflasNoTrans, m, n, k, one,
               <ModField.Element*>A, k, <ModField.Element*>B, n, zero,
               <ModField.Element*>ans, n, nbthreads)
    else :
        fgemm(F[0], FflasNoTrans, FflasNoTrans, m, n, k, one,
               <ModField.Element*>A, k, <ModField.Element*>B, n, zero,
               <ModField.Element*>ans, n)

    if m*n*k > 100000: sig_off()

    del F

cdef inline int linbox_matrix_vector_multiply(celement modulus, celement* C, celement* A, celement* b, Py_ssize_t m, Py_ssize_t n, FFLAS_TRANSPOSE trans):
    """
    C = A*v
    """
    cdef ModField *F = new ModField(<long>modulus)
    cdef ModField.Element one, zero
    F.init(one, <int>1)
    F.init(zero, <int>0)

    if m*n > 100000: sig_on()

    fgemv(F[0], trans,  m, n, one, <ModField.Element*>A, n, <ModField.Element*>b, 1,
               zero, <ModField.Element*>C, 1)

    if m*n > 100000: sig_off()

    del F

cdef inline linbox_minpoly(celement modulus, Py_ssize_t nrows, celement* entries):
    """
    Compute the minimal polynomial.
    """
    cdef Py_ssize_t i
    cdef ModField *F = new ModField(<long>modulus)
    cdef vector[ModField.Element] *minP = new vector[ModField.Element]()

    if nrows*nrows > 1000: sig_on()
    MinPoly(F[0], minP[0], nrows, <ModField.Element*>entries, nrows)
    if nrows*nrows > 1000: sig_off()

    l = []
    for i in range(minP.size()):
        l.append( <celement>minP.at(i) )

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

    if nrows*nrows > 1000: sig_on()
    CharPoly(R[0], P, nrows, <ModField.Element*>cpy, nrows)
    if nrows*nrows > 1000: sig_off()

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

    - ``X`` - a nonempty list of matrices of the same size mod a
       single modulus `n`

    EXAMPLES::

        sage: X = [random_matrix(GF(17), 4, 4) for _ in range(10)]; X
        [
        [ 2 14  0 15]  [12 14  3 13]  [ 9 15  8  1]  [ 2 12  6 10]
        [11 10 16  2]  [10  1 14  6]  [ 5  8 10 11]  [12  0  6  9]
        [ 9  4 10 14]  [ 2 14 13  7]  [ 5 12  4  9]  [ 7  7  3  8]
        [ 1 14  3 14], [ 6 14 10  3], [15  2  6 11], [ 2  9  1  5],
        <BLANKLINE>
        [12 13  7 16]  [ 5  3 16  2]  [14 15 16  4]  [ 1 15 11  0]
        [ 7 11 11  1]  [11 10 12 14]  [14  1 12 13]  [16 13  8 14]
        [ 0  2  0  4]  [ 0  7 16  4]  [ 5  5 16 13]  [13 14 16  4]
        [ 7  9  8 15], [ 6  5  2  3], [10 12  1  7], [15  6  6  6],
        <BLANKLINE>
        [ 4 10 11 15]  [13 12  5  1]
        [11  2  9 14]  [16 13 16  7]
        [12  5  4  4]  [12  2  0 11]
        [ 2  0 12  8], [13 11  6 15]
        ]
        sage: X[0]._matrix_from_rows_of_matrices(X) # indirect doctest
        [ 2 14  0 15 11 10 16  2  9  4 10 14  1 14  3 14]
        [12 14  3 13 10  1 14  6  2 14 13  7  6 14 10  3]
        [ 9 15  8  1  5  8 10 11  5 12  4  9 15  2  6 11]
        [ 2 12  6 10 12  0  6  9  7  7  3  8  2  9  1  5]
        [12 13  7 16  7 11 11  1  0  2  0  4  7  9  8 15]
        [ 5  3 16  2 11 10 12 14  0  7 16  4  6  5  2  3]
        [14 15 16  4 14  1 12 13  5  5 16 13 10 12  1  7]
        [ 1 15 11  0 16 13  8 14 13 14 16  4 15  6  6  6]
        [ 4 10 11 15 11  2  9 14 12  5  4  4  2  0 12  8]
        [13 12  5  1 16 13 16  7 12  2  0 11 13 11  6 15]

    OUTPUT: A single matrix mod ``p`` whose ``i``-th row is ``X[i].list()``.

    .. note::

         Do not call this function directly but use the static method
         ``Matrix_modn_dense_float/double._matrix_from_rows_of_matrices``
    """
    # The code below is just a fast version of the following:
    ##     from constructor import matrix
    ##     K = X[0].base_ring()
    ##     v = sum([y.list() for y in X],[])
    ##     return matrix(K, len(X), X[0].nrows()*X[0].ncols(), v)

    cdef Matrix_modn_dense_template T
    cdef Py_ssize_t i, n, m
    n = len(X)

    T = X[0]
    m = T._nrows * T._ncols
    cdef Matrix_modn_dense_template A = T.new_matrix(nrows = n, ncols = m)

    for i from 0 <= i < n:
        T = X[i]
        memcpy(A._entries + i*m, T._entries, sizeof(celement)*m)
    return A


cdef class Matrix_modn_dense_template(Matrix_dense):
    def __cinit__(self):
        cdef long p = self._base_ring.characteristic()
        self.p = p
        if p >= MAX_MODULUS:
            raise OverflowError("p (=%s) must be < %s."%(p, MAX_MODULUS))

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
            sage: for i in range(10):
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

        - ``coerce`` - perform modular reduction first?

        EXAMPLES::

            sage: A = random_matrix(GF(3),1000,1000)
            sage: type(A)
            <type 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
            sage: A = random_matrix(Integers(10),1000,1000)
            sage: type(A)
            <type 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
            sage: A = random_matrix(Integers(2^16),1000,1000)
            sage: type(A)
            <type 'sage.matrix.matrix_modn_dense_double.Matrix_modn_dense_double'>

        TESTS::

            sage: Matrix(GF(7), 2, 2, [-1, int(-2), GF(7)(-3), 1/4])
            [6 5]
            [4 2]

            sage: Matrix(GF(6434383), 2, 2, [-1, int(-2), GF(7)(-3), 1/4])
            [6434382 6434381]
            [      4 1608596]

            sage: Matrix(Integers(4618990), 2, 2, [-1, int(-2), GF(7)(-3), 1/7])
            [4618989 4618988]
            [      4 2639423]
        """
        ma = MatrixArgs_init(parent, entries)
        cdef long i, j
        it = ma.iter(False)
        R = ma.base
        p = R.characteristic()
        for i in range(ma.nrows):
            v = self._matrix[i]
            for j in range(ma.ncols):
                x = next(it)
                if type(x) is int:
                    tmp = (<long>x) % p
                    v[j] = tmp + (tmp<0)*p
                elif type(x) is IntegerMod_int and (<IntegerMod_int>x)._parent is R:
                    v[j] = <celement>(<IntegerMod_int>x).ivalue
                elif type(x) is Integer:
                    if coerce:
                        v[j] = mpz_fdiv_ui((<Integer>x).value, p)
                    else:
                        v[j] = mpz_get_ui((<Integer>x).value)
                elif coerce:
                    v[j] = R(x)
                else:
                    v[j] = <celement>x

    cdef long _hash_(self) except -1:
        """
        EXAMPLES::

            sage: B = random_matrix(GF(127),3,3)
            sage: B.set_immutable()
            sage: {B:0} # indirect doctest
            {[  9  75  94]
            [  4  57 112]
            [ 59  85  45]: 0}

            sage: M = random_matrix(GF(7), 10, 10)
            sage: M.set_immutable()
            sage: hash(M)
            -5724333594806680561  # 64-bit
            -1581874161           # 32-bit
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
                raise ValueError("incorrect size in matrix pickle (expected %d, got %d)"%(expectedlen, buflen))

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
                            v  = <mod_int>(us[0])
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
                            v  = <mod_int>(us[word_size-1])
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

        M = self.__class__.__new__(self.__class__, self._parent,None,None,None)

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
        cdef Py_ssize_t i,j
        cdef Matrix_modn_dense_template M
        cdef celement p = self.p
        cdef celement a = left

        M = self.__class__.__new__(self.__class__, self._parent,None,None,None)

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
        A = self.__class__.__new__(self.__class__, self._parent, 0, 0, 0)
        memcpy(A._entries, self._entries, sizeof(celement)*self._nrows*self._ncols)
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A


    cpdef _add_(self, right):
        """
        Add two dense matrices over `\Z/n\Z`

        INPUT:

        - ``right`` - a matrix

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

        M = self.__class__.__new__(self.__class__, self._parent,None,None,None)
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
        Subtract two dense matrices over `\Z/n\Z`

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

        M = self.__class__.__new__(self.__class__, self._parent, None, None, None)
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
            sage: A.change_ring(GF(17)) == A.change_ring(GF(19))
            False
            sage: A.change_ring(GF(17)) == A.change_ring(Integers(2000))
            False
            sage: A.change_ring(GF(17)) == A.change_ring(Integers(2000))
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
        return ``self*right``

        INPUT:

        - ``right``-  a matrix

        EXAMPLES::

            sage: A = random_matrix(GF(7),2,2); A
            [3 1]
            [6 6]

            sage: B = random_matrix(GF(7),2,2); B
            [4 4]
            [2 2]

            sage: A*B
            [0 0]
            [1 1]

            sage: 3*A
            [2 3]
            [4 4]

            sage: MS = parent(A)
            sage: MS(3) * A
            [2 3]
            [4 4]

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

            sage: A = random_matrix(Integers(8),2,2); A
            [0 5]
            [6 4]

            sage: B = random_matrix(Integers(8),2,2); B
            [4 4]
            [5 6]

            sage: A*B
            [1 6]
            [4 0]

            sage: 3*A
            [0 7]
            [2 4]

            sage: MS = parent(A)
            sage: MS(3) * A
            [0 7]
            [2 4]

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

            sage: A = random_matrix(GF(16007),2,2); A
            [ 6194 13327]
            [ 5985  5926]

            sage: B = random_matrix(GF(16007),2,2); B
            [ 6901  1242]
            [13032   859]

            sage: A*B
            [ 7618 12476]
            [14289  6330]

            sage: 3*A
            [2575 7967]
            [1948 1771]

            sage: MS = parent(A)
            sage: MS(3) * A
            [2575 7967]
            [1948 1771]

        ::

            sage: A = random_matrix(GF(15991), 201, 117)
            sage: B = random_matrix(GF(15991), 117, 195)
            sage: C = random_matrix(GF(15991), 201, 117)
            sage: D = random_matrix(GF(15991), 117, 195)

            sage: E = (A+C)*(B+D)

            sage: F = A*B + A*D + C*B + C*D

            sage: E == F
            True

        ::

            sage: A = random_matrix(GF(16007), 200, 200)
            sage: MS = parent(A)
            sage: (MS(0) * A) == 0
            True

            sage: (MS(1) * A) == A
            True

        ::

            sage: A = random_matrix(Integers(1008),2,2); A
            [ 41 973]
            [851 876]

            sage: B = random_matrix(Integers(1008),2,2); B
            [180 234]
            [680 640]

            sage: A*B
            [716 298]
            [924 750]

            sage: 3*A
            [123 903]
            [537 612]

            sage: MS = parent(A)
            sage: MS(3) * A
            [123 903]
            [537 612]

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
            verbose('mod-p multiply of %s x %s matrix by %s x %s matrix modulo %s'%(
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
        ``v*self``

        INPUT:

        - ``v`` - a vector

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 20)
            sage: v = random_vector(GF(17), 10)
            sage: matrix(v*A) == matrix(v)*A
            True

            sage: A = random_matrix(Integers(126), 10, 20)
            sage: v = random_vector(Integers(126), 10)
            sage: matrix(v*A) == matrix(v)*A
            True

            sage: A = random_matrix(GF(4796509), 10, 20)
            sage: v = random_vector(GF(4796509), 10)
            sage: matrix(v*A) == matrix(v)*A
            True

            sage: A = random_matrix(Integers(16337), 10, 20)
            sage: v = random_vector(Integers(16337), 10)
            sage: matrix(v*A) == matrix(v)*A
            True

        """
        if not isinstance(v, Vector_modn_dense):
            return (self.new_matrix(1,self._nrows, entries=v.list()) * self)[0]

        M = self._row_ambient_module()
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
        ``self*v``

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 20)
            sage: v = random_vector(GF(17), 20)
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True

            sage: A = random_matrix(Integers(126), 10, 20)
            sage: v = random_vector(Integers(126), 20)
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True

            sage: A = random_matrix(GF(4796509), 10, 20)
            sage: v = random_vector(GF(4796509), 20)
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True

            sage: A = random_matrix(Integers(16337), 10, 20)
            sage: v = random_vector(Integers(16337), 20)
            sage: matrix(A*v).transpose() == A*matrix(v).transpose()
            True
        """
        if not isinstance(v, Vector_modn_dense):
            r = (self * self.new_matrix(nrows=len(v), ncols=1, entries=v.list()))
            from sage.modules.free_module_element import vector
            return vector(r.list())

        M = self._column_ambient_module()
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

        - ``var`` - a variable name

        - ``algorithm`` - 'generic', 'linbox' or 'all' (default: linbox)

        EXAMPLES::

            sage: A = random_matrix(GF(19), 10, 10); A
            [ 3  1  8 10  5 16 18  9  6  1]
            [ 5 14  4  4 14 15  5 11  3  0]
            [ 4  1  0  7 11  6 17  8  5  6]
            [ 4  6  9  4  8  1 18 17  8 18]
            [11  2  0  6 13  7  4 11 16 10]
            [12  6 12  3 15 10  5 11  3  8]
            [15  1 16  2 18 15 14  7  2 11]
            [16 16 17  7 14 12  7  7  0  5]
            [13 15  9  2 12 16  1 15 18  7]
            [10  8 16 18  9 18  2 13  5 10]

            sage: B = copy(A)
            sage: char_p = A.characteristic_polynomial(); char_p
            x^10 + 2*x^9 + 18*x^8 + 4*x^7 + 13*x^6 + 11*x^5 + 2*x^4 + 5*x^3 + 7*x^2 + 16*x + 6
            sage: char_p(A) == 0
            True
            sage: B == A              # A is not modified
            True

            sage: min_p = A.minimal_polynomial(proof=True); min_p
            x^10 + 2*x^9 + 18*x^8 + 4*x^7 + 13*x^6 + 11*x^5 + 2*x^4 + 5*x^3 + 7*x^2 + 16*x + 6
            sage: min_p.divides(char_p)
            True

        ::

            sage: A = random_matrix(GF(2916337), 7, 7); A
            [ 514193 1196222 1242955 1040744   99523 2447069   40527]
            [ 930282 2685786 2892660 1347146 1126775 2131459  869381]
            [1853546 2266414 2897342 1342067 1054026  373002   84731]
            [1270068 2421818  569466  537440  572533  297105 1415002]
            [2079710  355705 2546914 2299052 2883413 1558788 1494309]
            [1027319 1572148  250822  522367 2516720  585897 2296292]
            [1797050 2128203 1161160  562535 2875615 1165768  286972]

            sage: B = copy(A)
            sage: char_p = A.characteristic_polynomial(); char_p
            x^7 + 1274305*x^6 + 1497602*x^5 + 12362*x^4 + 875330*x^3 + 31311*x^2 + 1858466*x + 700510
            sage: char_p(A) == 0
            True
            sage: B == A               # A is not modified
            True

            sage: min_p = A.minimal_polynomial(proof=True); min_p
            x^7 + 1274305*x^6 + 1497602*x^5 + 12362*x^4 + 875330*x^3 + 31311*x^2 + 1858466*x + 700510
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
            sage: A.minimal_polynomial()
            x

            sage: A = random_matrix(GF(4198973), 0, 0)
            sage: A.minimal_polynomial()
            1

            sage: A = random_matrix(GF(4198973), 0, 1)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = random_matrix(GF(4198973), 1, 0)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = matrix(GF(4198973), 10, 10)
            sage: A.minimal_polynomial()
            x

            sage: A = Mat(GF(7),3,3)([0, 1, 2] * 3)
            sage: A.charpoly()
            x^3 + 4*x^2

        ALGORITHM: Uses LinBox if ``self.base_ring()`` is a field,
        otherwise use Hessenberg form algorithm.

        TESTS:

        The cached polynomial should be independent of the ``var``
        argument (:trac:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached. The polynomial here is not unique, so we only
        check the polynomial's variable.

            sage: M = MatrixSpace(Integers(37), 2)
            sage: A = M(range(0, 2^2))
            sage: type(A)
            <type 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
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
            algorithm = 'generic' # LinBox only supports Z/pZ (p prime)

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
        Returns the minimal polynomial of`` self``.

        INPUT:

        - ``var`` - a variable name

        - ``algorithm`` - ``generic`` or ``linbox`` (default:
          ``linbox``)

        - ``proof`` -- (default: ``True``); whether to provably return
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

            sage: A = random_matrix(GF(17), 10, 10); A
            [ 2 14  0 15 11 10 16  2  9  4]
            [10 14  1 14  3 14 12 14  3 13]
            [10  1 14  6  2 14 13  7  6 14]
            [10  3  9 15  8  1  5  8 10 11]
            [ 5 12  4  9 15  2  6 11  2 12]
            [ 6 10 12  0  6  9  7  7  3  8]
            [ 2  9  1  5 12 13  7 16  7 11]
            [11  1  0  2  0  4  7  9  8 15]
            [ 5  3 16  2 11 10 12 14  0  7]
            [16  4  6  5  2  3 14 15 16  4]

            sage: B = copy(A)
            sage: min_p = A.minimal_polynomial(proof=True); min_p
            x^10 + 13*x^9 + 10*x^8 + 9*x^7 + 10*x^6 + 4*x^5 + 10*x^4 + 10*x^3 + 12*x^2 + 14*x + 7
            sage: min_p(A) == 0
            True
            sage: B == A
            True

            sage: char_p = A.characteristic_polynomial(); char_p
            x^10 + 13*x^9 + 10*x^8 + 9*x^7 + 10*x^6 + 4*x^5 + 10*x^4 + 10*x^3 + 12*x^2 + 14*x + 7
            sage: min_p.divides(char_p)
            True

        ::

            sage: A = random_matrix(GF(1214471), 10, 10); A
            [ 160562  831940   65852  173001  515930  714380  778254  844537  584888  392730]
            [ 502193  959391  614352  775603  240043 1156372  104118 1175992  612032 1049083]
            [ 660489 1066446  809624   15010 1002045  470722  314480 1155149 1173111   14213]
            [1190467 1079166  786442  429883  563611  625490 1015074  888047 1090092  892387]
            [   4724  244901  696350  384684  254561  898612   44844   83752 1091581  349242]
            [ 130212  580087  253296  472569  913613  919150   38603  710029  438461  736442]
            [ 943501  792110  110470  850040  713428  668799 1122064  325250 1084368  520553]
            [1179743  791517   34060 1183757 1118938  642169   47513   73428 1076788  216479]
            [ 626571  105273  400489 1041378 1186801  158611  888598 1138220 1089631   56266]
            [1092400  890773 1060810  211135  719636 1011640  631366  427711  547497 1084281]

            sage: B = copy(A)
            sage: min_p = A.minimal_polynomial(proof=True); min_p
            x^10 + 384251*x^9 + 702437*x^8 + 960299*x^7 + 202699*x^6 + 409368*x^5 + 1109249*x^4 + 1163061*x^3 + 333802*x^2 + 273775*x + 55190

            sage: min_p(A) == 0
            True
            sage: B == A
            True

            sage: char_p = A.characteristic_polynomial(); char_p
            x^10 + 384251*x^9 + 702437*x^8 + 960299*x^7 + 202699*x^6 + 409368*x^5 + 1109249*x^4 + 1163061*x^3 + 333802*x^2 + 273775*x + 55190

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
            sage: A.minimal_polynomial()
            x

        ::

            sage: A = random_matrix(GF(2535919), 0, 0)
            sage: A.minimal_polynomial()
            1

            sage: A = random_matrix(GF(2535919), 0, 1)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = random_matrix(GF(2535919), 1, 0)
            sage: A.minimal_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square

            sage: A = matrix(GF(2535919), 10, 10)
            sage: A.minimal_polynomial()
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
            algorithm='generic' # LinBox only supports fields

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
            raise ValueError("no algorithm '%s'"%algorithm)

        self.cache('minpoly_%s_%s'%(algorithm, var), g)
        return g

    def _charpoly_linbox(self, var='x'):
        """
        Computes the characteristic polynomial using LinBox. No checks
        are performed.

        This function is called internally by ``charpoly``.

        INPUT:

        - ``var`` - a variable name

        EXAMPLES::

            sage: A = random_matrix(GF(19), 10, 10); A
            [ 3  1  8 10  5 16 18  9  6  1]
            [ 5 14  4  4 14 15  5 11  3  0]
            [ 4  1  0  7 11  6 17  8  5  6]
            [ 4  6  9  4  8  1 18 17  8 18]
            [11  2  0  6 13  7  4 11 16 10]
            [12  6 12  3 15 10  5 11  3  8]
            [15  1 16  2 18 15 14  7  2 11]
            [16 16 17  7 14 12  7  7  0  5]
            [13 15  9  2 12 16  1 15 18  7]
            [10  8 16 18  9 18  2 13  5 10]

            sage: B = copy(A)
            sage: char_p = A._charpoly_linbox(); char_p
            x^10 + 2*x^9 + 18*x^8 + 4*x^7 + 13*x^6 + 11*x^5 + 2*x^4 + 5*x^3 + 7*x^2 + 16*x + 6
            sage: char_p(A) == 0
            True
            sage: B == A              # A is not modified
            True

            sage: min_p = A.minimal_polynomial(proof=True); min_p
            x^10 + 2*x^9 + 18*x^8 + 4*x^7 + 13*x^6 + 11*x^5 + 2*x^4 + 5*x^3 + 7*x^2 + 16*x + 6
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

    def echelonize(self, algorithm="linbox_noefd", **kwds):
        """
        Put ``self`` in reduced row echelon form.

        INPUT:

        - ``self`` - a mutable matrix

        - ``algorithm``

          - ``linbox`` - uses the LinBox library (wrapping fflas-ffpack)

          - ``linbox_noefd`` - uses the FFPACK directly, less memory and faster (default)

          - ``gauss`` - uses a custom slower `O(n^3)` Gauss
            elimination implemented in Sage.

          - ``all`` - compute using both algorithms and verify that
            the results are the same.

        - ``**kwds`` - these are all ignored

        OUTPUT:

        - ``self`` is put in reduced row echelon form.

        - the rank of self is computed and cached

        - the pivot columns of self are computed and cached.

        - the fact that self is now in echelon form is recorded and
          cached so future calls to echelonize return immediately.

        EXAMPLES::

            sage: A = random_matrix(GF(7), 10, 20); A
            [3 1 6 6 4 4 2 2 3 5 4 5 6 2 2 1 2 5 0 5]
            [3 2 0 5 0 1 5 4 2 3 6 4 5 0 2 4 2 0 6 3]
            [2 2 4 2 4 5 3 4 4 4 2 5 2 5 4 5 1 1 1 1]
            [0 6 3 4 2 2 3 5 1 1 4 2 6 5 6 3 4 5 5 3]
            [5 2 4 3 6 2 3 6 2 1 3 3 5 3 4 2 2 1 6 2]
            [0 5 6 3 2 5 6 6 3 2 1 4 5 0 2 6 5 2 5 1]
            [4 0 4 2 6 3 3 5 3 0 0 1 2 5 5 1 6 0 0 3]
            [2 0 1 0 0 3 0 2 4 2 2 4 4 4 5 4 1 2 3 4]
            [2 4 1 4 3 0 6 2 2 5 2 5 3 6 4 2 2 6 4 4]
            [0 0 2 2 1 6 2 0 5 0 4 3 1 6 0 6 0 4 6 5]

            sage: A.echelon_form()
            [1 0 0 0 0 0 0 0 0 0 6 2 6 0 1 1 2 5 6 2]
            [0 1 0 0 0 0 0 0 0 0 0 4 5 4 3 4 2 5 1 2]
            [0 0 1 0 0 0 0 0 0 0 6 3 4 6 1 0 3 6 5 6]
            [0 0 0 1 0 0 0 0 0 0 0 3 5 2 3 4 0 6 5 3]
            [0 0 0 0 1 0 0 0 0 0 0 6 3 4 5 3 0 4 3 2]
            [0 0 0 0 0 1 0 0 0 0 1 1 0 2 4 2 5 5 5 0]
            [0 0 0 0 0 0 1 0 0 0 1 0 1 3 2 0 0 0 5 3]
            [0 0 0 0 0 0 0 1 0 0 4 4 2 6 5 4 3 4 1 0]
            [0 0 0 0 0 0 0 0 1 0 1 0 4 2 3 5 4 6 4 0]
            [0 0 0 0 0 0 0 0 0 1 2 0 5 0 5 5 3 1 1 4]

        ::

            sage: A = random_matrix(GF(13), 10, 10); A
            [ 8  3 11 11  9  4  8  7  9  9]
            [ 2  9  6  5  7 12  3  4 11  5]
            [12  6 11 12  4  3  3  8  9  5]
            [ 4  2 10  5 10  1  1  1  6  9]
            [12  8  5  5 11  4  1  2  8 11]
            [ 2  6  9 11  4  7  1  0 12  2]
            [ 8  9  0  7  7  7 10  4  1  4]
            [ 0  8  2  6  7  5  7 12  2  3]
            [ 2 11 12  3  4  7  2  9  6  1]
            [ 0 11  5  9  4  5  5  8  7 10]

            sage: MS = parent(A)
            sage: B = A.augment(MS(1))
            sage: B.echelonize()
            sage: A.rank()
            10
            sage: C = B.submatrix(0,10,10,10); C
            [ 4  9  4  4  0  4  7 11  9 11]
            [11  7  6  8  2  8  6 11  9  5]
            [ 3  9  9  2  4  8  9  2  9  4]
            [ 7  0 11  4  0  9  6 11  8  1]
            [12 12  4 12  3 12  6  1  7 12]
            [12  2 11  6  6  6  7  0 10  6]
            [ 0  7  3  4  7 11 10 12  4  6]
            [ 5 11  0  5  3 11  4 12  5 12]
            [ 6  7  3  5  1  4 11  7  4  1]
            [ 4  9  6  7 11  1  2 12  6  7]

            sage: ~A == C
            True

        ::

            sage: A = random_matrix(Integers(10), 10, 20)
            sage: A.echelon_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 10'.

        ::

            sage: A = random_matrix(GF(16007), 10, 20); A
            [15455  1177 10072  4693  3887  4102 10746 15265  6684 14559  4535 13921  9757  9525  9301  8566  2460  9609  3887  6205]
            [ 8602 10035  1242  9776   162  7893 12619  6660 13250  1988 14263 11377  2216  1247  7261  8446 15081 14412  7371  7948]
            [12634  7602   905  9617 13557  2694 13039  4936 12208 15480  3787 11229   593 12462  5123 14167  6460  3649  5821  6736]
            [10554  2511 11685 12325 12287  6534 11636  5004  6468  3180  3607 11627 13436  5106  3138 13376  8641  9093  2297  5893]
            [ 1025 11376 10288   609 12330  3021   908 13012  2112 11505    56  5971   338  2317  2396  8561  5593  3782  7986 13173]
            [ 7607   588  6099 12749 10378   111  2852 10375  8996  7969   774 13498 12720  4378  6817  6707  5299  9406 13318  2863]
            [15545   538  4840  1885  8471  1303 11086 14168  1853 14263  3995 12104  1294  7184  1188 11901 15971  2899  4632   711]
            [  584 11745  7540 15826 15027  5953  7097 14329 10889 12532 13309 15041  6211  1749 10481  9999  2751 11068    21  2795]
            [  761 11453  3435 10596  2173  7752 15941 14610  1072  8012  9458  5440   612 10581 10400   101 11472 13068  7758  7898]
            [10658  4035  6662   655  7546  4107  6987  1877  4072  4221  7679 14579  2474  8693  8127 12999 11141   605  9404 10003]
            sage: A.echelon_form()
            [    1     0     0     0     0     0     0     0     0     0  8416  8364 10318  1782 13872  4566 14855  7678 11899  2652]
            [    0     1     0     0     0     0     0     0     0     0  4782 15571  3133 10964  5581 10435  9989 14303  5951  8048]
            [    0     0     1     0     0     0     0     0     0     0 15688  6716 13819  4144   257  5743 14865 15680  4179 10478]
            [    0     0     0     1     0     0     0     0     0     0  4307  9488  2992  9925 13984 15754  8185 11598 14701 10784]
            [    0     0     0     0     1     0     0     0     0     0   927  3404 15076  1040  2827  9317 14041 10566  5117  7452]
            [    0     0     0     0     0     1     0     0     0     0  1144 10861  5241  6288  9282  5748  3715 13482  7258  9401]
            [    0     0     0     0     0     0     1     0     0     0   769  1804  1879  4624  6170  7500 11883  9047   874   597]
            [    0     0     0     0     0     0     0     1     0     0 15591 13686  5729 11259 10219 13222 15177 15727  5082 11211]
            [    0     0     0     0     0     0     0     0     1     0  8375 14939 13471 12221  8103  4212 11744 10182  2492 11068]
            [    0     0     0     0     0     0     0     0     0     1  6534   396  6780 14734  1206  3848  7712  9770 10755   410]

        ::

            sage: A = random_matrix(Integers(10000), 10, 20)
            sage: A.echelon_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 10000'.

        Parallel computation::

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
        if not x is None:
            return  # already known to be in echelon form

        if not self.base_ring().is_field():
            raise NotImplementedError("Echelon form not implemented over '%s'."%self.base_ring())

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
            raise ValueError("Algorithm '%s' not known"%algorithm)

    def _echelonize_linbox(self, efd=True):
        """
        Puts ``self`` in row echelon form using LinBox.

        This function is called by echelonize if
        ``algorithm='linbox'``.

        INPUT:

        - ``efd`` - if ``True`` LinBox's ``EchelonFormDomain``
          implementation is used, which is faster than the direct
          ``LinBox::FFPACK`` implementation, since the latter also
          computes the transformation matrix (which we
          ignore). However, ``efd=True`` uses more memory than FFLAS
          directly (default=``True``)

        EXAMPLES::

            sage: A = random_matrix(GF(7), 10, 20); A
            [3 1 6 6 4 4 2 2 3 5 4 5 6 2 2 1 2 5 0 5]
            [3 2 0 5 0 1 5 4 2 3 6 4 5 0 2 4 2 0 6 3]
            [2 2 4 2 4 5 3 4 4 4 2 5 2 5 4 5 1 1 1 1]
            [0 6 3 4 2 2 3 5 1 1 4 2 6 5 6 3 4 5 5 3]
            [5 2 4 3 6 2 3 6 2 1 3 3 5 3 4 2 2 1 6 2]
            [0 5 6 3 2 5 6 6 3 2 1 4 5 0 2 6 5 2 5 1]
            [4 0 4 2 6 3 3 5 3 0 0 1 2 5 5 1 6 0 0 3]
            [2 0 1 0 0 3 0 2 4 2 2 4 4 4 5 4 1 2 3 4]
            [2 4 1 4 3 0 6 2 2 5 2 5 3 6 4 2 2 6 4 4]
            [0 0 2 2 1 6 2 0 5 0 4 3 1 6 0 6 0 4 6 5]

            sage: A._echelonize_linbox(); A
            [1 0 0 0 0 0 0 0 0 0 6 2 6 0 1 1 2 5 6 2]
            [0 1 0 0 0 0 0 0 0 0 0 4 5 4 3 4 2 5 1 2]
            [0 0 1 0 0 0 0 0 0 0 6 3 4 6 1 0 3 6 5 6]
            [0 0 0 1 0 0 0 0 0 0 0 3 5 2 3 4 0 6 5 3]
            [0 0 0 0 1 0 0 0 0 0 0 6 3 4 5 3 0 4 3 2]
            [0 0 0 0 0 1 0 0 0 0 1 1 0 2 4 2 5 5 5 0]
            [0 0 0 0 0 0 1 0 0 0 1 0 1 3 2 0 0 0 5 3]
            [0 0 0 0 0 0 0 1 0 0 4 4 2 6 5 4 3 4 1 0]
            [0 0 0 0 0 0 0 0 1 0 1 0 4 2 3 5 4 6 4 0]
            [0 0 0 0 0 0 0 0 0 1 2 0 5 0 5 5 3 1 1 4]
        """
        self.check_mutability()
        self.clear_cache()

        t = verbose('Calling echelonize mod %d.'%self.p)
        if efd:
            r, pivots = linbox_echelonize_efd(self.p, self._entries, self._nrows, self._ncols)
        else:
            r, pivots = linbox_echelonize(self.p, self._entries, self._nrows, self._ncols)
        verbose('done with echelonize',t)
        self.cache('in_echelon_form',True)
        self.cache('rank', r)
        self.cache('pivots', tuple(pivots))

    def _echelon_in_place_classical(self):
        """
        Puts ``self`` in row echelon form using LinBox.

        This function is called by echelonize if
        ``algorithm='gauss'``.

        EXAMPLES::

            sage: A = random_matrix(GF(7), 10, 20); A
            [3 1 6 6 4 4 2 2 3 5 4 5 6 2 2 1 2 5 0 5]
            [3 2 0 5 0 1 5 4 2 3 6 4 5 0 2 4 2 0 6 3]
            [2 2 4 2 4 5 3 4 4 4 2 5 2 5 4 5 1 1 1 1]
            [0 6 3 4 2 2 3 5 1 1 4 2 6 5 6 3 4 5 5 3]
            [5 2 4 3 6 2 3 6 2 1 3 3 5 3 4 2 2 1 6 2]
            [0 5 6 3 2 5 6 6 3 2 1 4 5 0 2 6 5 2 5 1]
            [4 0 4 2 6 3 3 5 3 0 0 1 2 5 5 1 6 0 0 3]
            [2 0 1 0 0 3 0 2 4 2 2 4 4 4 5 4 1 2 3 4]
            [2 4 1 4 3 0 6 2 2 5 2 5 3 6 4 2 2 6 4 4]
            [0 0 2 2 1 6 2 0 5 0 4 3 1 6 0 6 0 4 6 5]

            sage: A._echelon_in_place_classical(); A
            [1 0 0 0 0 0 0 0 0 0 6 2 6 0 1 1 2 5 6 2]
            [0 1 0 0 0 0 0 0 0 0 0 4 5 4 3 4 2 5 1 2]
            [0 0 1 0 0 0 0 0 0 0 6 3 4 6 1 0 3 6 5 6]
            [0 0 0 1 0 0 0 0 0 0 0 3 5 2 3 4 0 6 5 3]
            [0 0 0 0 1 0 0 0 0 0 0 6 3 4 5 3 0 4 3 2]
            [0 0 0 0 0 1 0 0 0 0 1 1 0 2 4 2 5 5 5 0]
            [0 0 0 0 0 0 1 0 0 0 1 0 1 3 2 0 0 0 5 3]
            [0 0 0 0 0 0 0 1 0 0 4 4 2 6 5 4 3 4 1 0]
            [0 0 0 0 0 0 0 0 1 0 1 0 4 2 3 5 4 6 4 0]
            [0 0 0 0 0 0 0 0 0 1 2 0 5 0 5 5 3 1 1 4]
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
            if do_verb and (c % fifth == 0 and c>0):
                tm = verbose('on column %s of %s'%(c, self._ncols),
                             level = 2,
                             caller_name = 'matrix_modn_dense echelon')
            #end if
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
        self.cache('in_echelon_form',True)

    def right_kernel_matrix(self, algorithm='linbox', basis='echelon'):
        r"""
        Returns a matrix whose rows form a basis for the right kernel
        of ``self``, where ``self`` is a matrix over a (small) finite field.

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
            self = self.echelon_form(algorithm=algorithm)

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
        Transforms self in place to its Hessenberg form.

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 10, density=0.1); A
            [ 0  0  0  0 12  0  0  0  0  0]
            [ 0  0  0  4  0  0  0  0  0  0]
            [ 0  0  0  0  2  0  0  0  0  0]
            [ 0 14  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0 10  0  0  0  0]
            [ 0  0  0  0  0 16  0  0  0  0]
            [ 0  0  0  0  0  0  6  0  0  0]
            [15  0  0  0  0  0  0  0  0  0]
            [ 0  0  0 16  0  0  0  0  0  0]
            [ 0  5  0  0  0  0  0  0  0  0]
            sage: A.hessenbergize(); A
            [ 0  0  0  0  0  0  0 12  0  0]
            [15  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  2  0  0]
            [ 0  0  0  0 14  0  0  0  0  0]
            [ 0  0  0  4  0  0  0  0  0  0]
            [ 0  0  0  0  5  0  0  0  0  0]
            [ 0  0  0  0  0  0  6  0  0  0]
            [ 0  0  0  0  0  0  0  0  0 10]
            [ 0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0 16]
        """
        self.check_mutability()
        x = self.fetch('in_hessenberg_form')
        if not x is None and x: return  # already known to be in Hessenberg form

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
                 t_inv = celement_invert(t,p)
                 if i > m:
                     self.swap_rows_c(i,m)
                     self.swap_columns_c(i,m)

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
        self.cache('in_hessenberg_form',True)

    def _charpoly_hessenberg(self, var):
        """
        Transforms self in place to its Hessenberg form then computes
        and returns the coefficients of the characteristic polynomial
        of this matrix.

        INPUT:

        - ``var`` - name of the indeterminate of the charpoly.

        OUTPUT:

           The characteristic polynomial is represented as a vector of
           ints, where the constant term of the characteristic
           polynomial is the 0th coefficient of the vector.

        EXAMPLES::

            sage: A = random_matrix(GF(17), 10, 10, density=0.1); A
            [ 0  0  0  0 12  0  0  0  0  0]
            [ 0  0  0  4  0  0  0  0  0  0]
            [ 0  0  0  0  2  0  0  0  0  0]
            [ 0 14  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0 10  0  0  0  0]
            [ 0  0  0  0  0 16  0  0  0  0]
            [ 0  0  0  0  0  0  6  0  0  0]
            [15  0  0  0  0  0  0  0  0  0]
            [ 0  0  0 16  0  0  0  0  0  0]
            [ 0  5  0  0  0  0  0  0  0  0]
            sage: A.characteristic_polynomial()
            x^10 + 12*x^9 + 6*x^8 + 8*x^7 + 13*x^6
            sage: P.<x> = GF(17)[]
            sage: A._charpoly_hessenberg('x')
            x^10 + 12*x^9 + 6*x^8 + 8*x^7 + 13*x^6
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
        c = self.new_matrix(nrows=n+1,ncols=n+1)    # the 0 matrix
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
                c.add_multiple_of_row_c(m, m-i-1, p - (t*H._matrix[m-i-1][m-1])%p, 0)

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
            sage: A.rank()
            99
            sage: B == A
            True

            sage: A = random_matrix(GF(3), 100, 100, density=0.01)
            sage: A.rank()
            63

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
            if not x is None:
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

            sage: A = random_matrix(GF(7), 10, 10); A
            [3 1 6 6 4 4 2 2 3 5]
            [4 5 6 2 2 1 2 5 0 5]
            [3 2 0 5 0 1 5 4 2 3]
            [6 4 5 0 2 4 2 0 6 3]
            [2 2 4 2 4 5 3 4 4 4]
            [2 5 2 5 4 5 1 1 1 1]
            [0 6 3 4 2 2 3 5 1 1]
            [4 2 6 5 6 3 4 5 5 3]
            [5 2 4 3 6 2 3 6 2 1]
            [3 3 5 3 4 2 2 1 6 2]

            sage: A.determinant()
            6

        ::

            sage: A = random_matrix(GF(7), 100, 100)
            sage: A.determinant()
            2

            sage: A.transpose().determinant()
            2

            sage: B = random_matrix(GF(7), 100, 100)
            sage: B.determinant()
            4

            sage: (A*B).determinant() == A.determinant() * B.determinant()
            True

        ::

            sage: A = random_matrix(GF(16007), 10, 10); A
            [ 5037  2388  4150  1400   345  5945  4240 14022 10514   700]
            [15552  8539  1927  3870  9867  3263 11637   609 15424  2443]
            [ 3761 15836 12246 15577 10178 13602 13183 15918 13942  2958]
            [ 4526 10817  6887  6678  1764  9964  6107  1705  5619  5811]
            [13537 15004  8307 11846 14779   550 14113  5477  7271  7091]
            [13338  4927 11406 13065  5437 12431  6318  5119 14198   496]
            [ 1044   179 12881   353 12975 12567  1092 10433 12304   954]
            [10072  8821 14118 13895  6543 13484 10685 14363  2612 11070]
            [15113   237  2612 14127 11589  5808   117  9656 15957 14118]
            [15233 11080  5716  9029 11402  9380 13045 13986 14544  5771]

            sage: A.determinant()
            10207

        ::

            sage: A = random_matrix(GF(16007), 100, 100)
            sage: A.determinant()
            3576


            sage: A.transpose().determinant()
            3576

            sage: B = random_matrix(GF(16007), 100, 100)
            sage: B.determinant()
            4075

            sage: (A*B).determinant() == A.determinant() * B.determinant()
            True

        Parallel computation::

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

            sage: A = matrix(GF(7), 5, 5); A.det()
            0

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
            if not x is None:
                return x
            d = linbox_det(self.p, self._entries, self._nrows)
            d2 = self._coerce_element(d)
            self.cache('det', d2)
            return d2
        else:
            return Matrix_dense.determinant(self)

    cdef xgcd_eliminate(self, celement * row1, celement* row2, Py_ssize_t start_col):
        """
        Reduces ``row1`` and ``row2`` by a unimodular transformation
        using the xgcd relation between their first coefficients ``a`` and
        ``b``.

        INPUT:

        - ``row1, row2`` - the two rows to be transformed (within
          self)

        -``start_col`` - the column of the pivots in ``row1`` and
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
        g = ArithIntObj.c_xgcd_int (a,b,<int*>&s,<int*>&t)
        v = a/g
        w = -<int>b/g
        nc = self.ncols()

        for i from start_col <= i < nc:
            tmp = ( s * <int>row1[i] + t * <int>row2[i]) % p
            row2[i] = (w* <int>row1[i] + v*<int>row2[i]) % p
            row1[i] = tmp
        return g

    cdef rescale_row_c(self, Py_ssize_t row, multiple, Py_ssize_t start_col):
        """
        Rescale ``self[row]`` by ``multiple`` but only start at column
        index ``start_col``.

        INPUT:

        - ``row`` - integer
        - ``multiple`` - finite field element
        - ``start_col`` - integer

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

            sage: A = random_matrix(GF(37), 10, 10); A
            [24 15  7 27 32 34 16 32 25 23]
            [11  3 22 13 35 33  1 10 33 25]
            [33  9 25  3 15 27 30 30  7 12]
            [10  0 35  4 12 34 32 16 19 17]
            [36  4 21 17  3 34 11 10 10 17]
            [32 15 23  2 23 32  5  8 18 11]
            [24  5 28 13 21 22 29 18 33 30]
            [26 18 10 26 17 31 35 18 25 30]
            [21  1  4 14 11 17 29 16 18 12]
            [34 19 14 11 35 30 35 34 25 33]

            sage: A[2] + 10*A[3]
            (22, 9, 5, 6, 24, 34, 17, 5, 12, 34)

            sage: A.add_multiple_of_row(2, 3, 10)
            sage: A
            [24 15  7 27 32 34 16 32 25 23]
            [11  3 22 13 35 33  1 10 33 25]
            [22  9  5  6 24 34 17  5 12 34]
            [10  0 35  4 12 34 32 16 19 17]
            [36  4 21 17  3 34 11 10 10 17]
            [32 15 23  2 23 32  5  8 18 11]
            [24  5 28 13 21 22 29 18 33 30]
            [26 18 10 26 17 31 35 18 25 30]
            [21  1  4 14 11 17 29 16 18 12]
            [34 19 14 11 35 30 35 34 25 33]

            sage: A.add_multiple_of_row(2, 3, 10, 4)
            sage: A
            [24 15  7 27 32 34 16 32 25 23]
            [11  3 22 13 35 33  1 10 33 25]
            [22  9  5  6 33  4  4 17 17 19]
            [10  0 35  4 12 34 32 16 19 17]
            [36  4 21 17  3 34 11 10 10 17]
            [32 15 23  2 23 32  5  8 18 11]
            [24  5 28 13 21 22 29 18 33 30]
            [26 18 10 26 17 31 35 18 25 30]
            [21  1  4 14 11 17 29 16 18 12]
            [34 19 14 11 35 30 35 34 25 33]
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
            v_to[i] = ((<celement>multiple) * v_from[i] +  v_to[i]) % p

    cdef add_multiple_of_column_c(self, Py_ssize_t col_to, Py_ssize_t col_from, multiple, Py_ssize_t start_row):
        """
        Add ``multiple`` times ``self[row_from]`` to ``self[row_to]``
        statting in column ``start_col``.

        EXAMPLES::

            sage: A = random_matrix(GF(37), 10, 10); A
            [24 15  7 27 32 34 16 32 25 23]
            [11  3 22 13 35 33  1 10 33 25]
            [33  9 25  3 15 27 30 30  7 12]
            [10  0 35  4 12 34 32 16 19 17]
            [36  4 21 17  3 34 11 10 10 17]
            [32 15 23  2 23 32  5  8 18 11]
            [24  5 28 13 21 22 29 18 33 30]
            [26 18 10 26 17 31 35 18 25 30]
            [21  1  4 14 11 17 29 16 18 12]
            [34 19 14 11 35 30 35 34 25 33]

            sage: A.column(2) + 10*A.column(3)
            (18, 4, 18, 1, 6, 6, 10, 11, 33, 13)

            sage: A.add_multiple_of_column(2, 3, 10)
            sage: A
            [24 15 18 27 32 34 16 32 25 23]
            [11  3  4 13 35 33  1 10 33 25]
            [33  9 18  3 15 27 30 30  7 12]
            [10  0  1  4 12 34 32 16 19 17]
            [36  4  6 17  3 34 11 10 10 17]
            [32 15  6  2 23 32  5  8 18 11]
            [24  5 10 13 21 22 29 18 33 30]
            [26 18 11 26 17 31 35 18 25 30]
            [21  1 33 14 11 17 29 16 18 12]
            [34 19 13 11 35 30 35 34 25 33]

            sage: A.add_multiple_of_column(2, 3, 10, 4)
            sage: A
            [24 15 18 27 32 34 16 32 25 23]
            [11  3  4 13 35 33  1 10 33 25]
            [33  9 18  3 15 27 30 30  7 12]
            [10  0  1  4 12 34 32 16 19 17]
            [36  4 28 17  3 34 11 10 10 17]
            [32 15 26  2 23 32  5  8 18 11]
            [24  5 29 13 21 22 29 18 33 30]
            [26 18 12 26 17 31 35 18 25 30]
            [21  1 25 14 11 17 29 16 18 12]
            [34 19 12 11 35 30 35 34 25 33]
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

        - ``density`` - Integer; proportion (roughly) to be considered
           for changes
        - ``nonzero`` - Bool (default: ``False``); whether the new
           entries are forced to be non-zero

        OUTPUT:

        -  None, the matrix is modified in-space

        EXAMPLES::

            sage: A = matrix(GF(5), 5, 5, 0)
            sage: A.randomize(0.5); A
            [0 0 0 2 0]
            [0 3 0 0 2]
            [4 0 0 0 0]
            [4 0 0 0 0]
            [0 1 0 0 0]

            sage: A.randomize(); A
            [3 3 2 1 2]
            [4 3 3 2 2]
            [0 3 3 3 3]
            [3 3 2 2 4]
            [2 2 2 1 4]

        The matrix is updated instead of overwritten::

            sage: A = random_matrix(GF(5), 100, 100, density=0.1)
            sage: A.density()
            961/10000

            sage: A.randomize(density=0.1)
            sage: A.density()
            801/5000
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
        Returns a string representation of ``self`` in Magma form.

        INPUT:

        -  ``magma`` - a Magma session

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
        return 'Matrix(%s,%s,%s,StringToIntegerSequence("%s"))'%(
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
        cdef object P =  matrix_space.MatrixSpace(ZZ, self._nrows, self._ncols, sparse=False)
        L = Matrix_integer_dense(P,ZZ(0),False,False)
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

        cdef Matrix_modn_dense_template M = self.new_matrix(nrows = ncols, ncols = nrows)
        cdef Py_ssize_t i,j

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

    def submatrix(self, Py_ssize_t row=0, Py_ssize_t col=0,
                        Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        r"""
        Return the matrix constructed from self using the specified
        range of rows and columns.

        INPUT:

        - ``row``, ``col`` -- index of the starting row and column.
          Indices start at zero

        - ``nrows``, ``ncols`` -- (optional) number of rows and columns to
          take. If not provided, take all rows below and all columns to
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

        if col != 0 or ncols != self._ncols:
            return self.matrix_from_rows_and_columns(range(row, row+nrows), range(col, col+ncols))

        if nrows < 0 or row < 0 or row + nrows > self._nrows:
            raise IndexError("rows out of range")

        cdef Matrix_modn_dense_template M = self.new_matrix(nrows=nrows, ncols=self._ncols)
        memcpy(M._entries, self._entries+row*ncols, sizeof(celement)*ncols*nrows)
        return M

    def _matrices_from_rows(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Make a list of matrix from the rows of this matrix.  This is a
        fairly technical function which is used internally, e.g., by
        the cyclotomic field linear algebra code.

        INPUT:

        - ``nrows`` - integer

        - ``ncols`` - integer

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

        OUTPUT:

        - ``list`` - a list of matrices
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

    def __nonzero__(self):
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

    cdef int _copy_row_to_mod_int_array(self, mod_int *to, Py_ssize_t i):
        cdef Py_ssize_t j
        cdef celement *_from = self._entries+(i*self._ncols)
        for j in range(self._ncols):
            to[j] = <mod_int>_from[j]

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        r"""
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: M = Matrix(GF(49), 2, [1,2,-2,0])
            sage: M.zero_pattern_matrix()  # indirect doctest
            [0 0]
            [0 1]

            sage: M = Matrix(Integers(10), 2, [1,2,-2,0])
            sage: M.zero_pattern_matrix()  # indirect doctest
            [0 0]
            [0 1]
        """
        return self._entries[j+i*self._ncols] == 0

