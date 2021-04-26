# distutils: extra_compile_args = -D_XPG6
# flags chosen from libs/flint/fmpz_poly.pyx
r"""
FLINT nmod_mat class wrapper

This file implements matrices over `\ZZ/N\ZZ` for `N < 2^{63}`.
It adds some capabilities for composite `N` that are not present in FLINT.

AUTHORS:

- Edgar Costa, David Roe (2021) Initial version.
"""

from cpython.sequence cimport *
from cysignals.signals cimport sig_on, sig_str, sig_off
from libc.string cimport memcpy

from sage.structure.element cimport Element, Matrix
from sage.structure.element import is_Vector
from sage.categories.cartesian_product import cartesian_product
from sage.structure.richcmp cimport rich_to_bool, Py_EQ, Py_NE
from sage.structure.factorization import Factorization
from sage.structure.proof.all import linear_algebra as linalg_proof
from sage.misc.prandom import randrange
from sage.arith.all import gcd
from sage.arith.power cimport generic_power
from sage.arith.long cimport integer_check_long_py
from sage.rings.polynomial.polynomial_zmod_flint cimport Polynomial_zmod_flint
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.integer_mod cimport (
    IntegerMod_abstract,
    IntegerMod_int,
    IntegerMod_int64,
    IntegerMod_gmp
)
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from .args cimport SparseEntry, MatrixArgs_init

from sage.libs.flint.nmod_mat cimport *
from sage.libs.flint.nmod_poly cimport (
    nmod_poly_set,
    nmod_poly_set_coeff_ui,
    nmod_poly_get_coeff_ui,
    nmod_poly_fit_length
)
from sage.libs.flint.ulong_extras cimport (
    n_precompute_inverse,
    n_preinvert_limb,
    n_pow,
    n_addmod,
    n_submod,
    n_negmod,
    n_invmod,
    n_mod2_preinv,
    n_mulmod2_preinv,
    n_remove2_precomp,
    n_CRT,
    n_gcd,
    n_div2_preinv,
    n_divrem2_precomp,
)
from sage.libs.gmp.mpz cimport mpz_sgn,  mpz_fits_ulong_p, mpz_get_ui, mpz_get_si


cdef class Matrix_modn_dense_flint(Matrix_dense):
    r"""
    Matrices modulo `N` for `N < 2^{63}`

    EXAMPLES::

        sage: A = matrix(Zmod(36), 3, 3, range(9))
        sage: type(A)
        <class 'sage.matrix.matrix_modn_dense_flint.Matrix_modn_dense_flint'>
    """
    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    def __cinit__(self, parent, *args, **kwds):
        """
        Memory initialization

        EXAMPLES::

            sage: M = MatrixSpace(Zmod(36), 2, 2)
            sage: M()
            [0 0]
            [0 0]
        """
        self._modulus = parent._base._pyx_order
        sig_str("FLINT exception")
        nmod_mat_init(self._matrix, self._nrows, self._ncols, mpz_get_ui(self._modulus.sageInteger.value))
        sig_off()

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        """
        Initialization

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: TestSuite(A).run()
        """
        self._parent = parent
        ma = MatrixArgs_init(parent, entries)
        cdef SparseEntry se
        for se in ma.iter(coerce, sparse=True):
            nmod_mat_set_entry(self._matrix, se.i, se.j, se.entry)

    def __dealloc__(self):
        """
        Memory deallocation

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: del A
        """
        sig_on()
        nmod_mat_clear(self._matrix)
        sig_off()

    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        nmod_mat_set_entry(self._matrix, i, j, value)

    cdef void set_unsafe_ui(self, Py_ssize_t i, Py_ssize_t j, unsigned long value):
        nmod_mat_set_entry(self._matrix, i, j, value)

    cdef unsigned long get_unsafe_ui(self, Py_ssize_t i, Py_ssize_t j):
        return nmod_mat_get_entry(self._matrix, i, j)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Low level interface for setting entries, used by generic matrix code.
        """
        e = self._modulus.element_class()
        cdef mp_limb_t ivalue
        if e is IntegerMod_int:
            ivalue = (<IntegerMod_int?>x).ivalue
        elif e is IntegerMod_int64:
            ivalue = (<IntegerMod_int64?>x).ivalue
        else:
            ivalue = mpz_get_ui((<IntegerMod_gmp?>x).value)
        nmod_mat_set_entry(self._matrix, i, j, ivalue)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Low level interface for getting entries as an element of the base ring,
        used by generic matrix code.
        """
        cdef type t = self._modulus.element_class()
        cdef IntegerMod_abstract x = t.__new__(t)
        x._parent = self._parent._base
        x.__modulus = self._modulus
        x.set_from_ulong_fast(nmod_mat_get_entry(self._matrix, i, j))
        return x

    cdef Matrix_modn_dense_flint _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Return a new matrix over the parent from given parent
        All memory is allocated for this matrix, but its
        entries have not yet been filled in.
        """
        if nrows == self._nrows and ncols == self._ncols:
            P = self._parent
        else:
            from sage.matrix.matrix_space import MatrixSpace
            P = MatrixSpace(self._parent._base, nrows, ncols, sparse=False)
        cdef Matrix_modn_dense_flint ans = Matrix_modn_dense_flint.__new__(Matrix_modn_dense_flint, P)
        ans._parent = P
        return ans

    def _change_implementation(self, implementation):
        """
        Return this matrix with a different underlying implementation.

        INPUT::

        - ``implementation`` -- an argument accepted as an implementation argument to :class:`sage.matrix.matrix_space.MatrixSpace`.

        EXAMPLES::

            sage: A = matrix(Zmod(5), 2, range(4))
            sage: type(A)
            <class 'sage.matrix.matrix_modn_dense_flint.Matrix_modn_dense_flint'>
            sage: B = A._change_implementation("linbox")
            sage: type(B)
            <class 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
        """
        from sage.matrix.matrix_space import MatrixSpace
        P = MatrixSpace(self._parent._base, self._nrows, self._ncols, sparse=False, implementation=implementation)
        if P is self.parent():
            return self
        cdef Matrix_dense mat = P()
        cdef Py_ssize_t i, j
        for i in range(self._nrows):
            for j in range(self._ncols):
                mat.set_unsafe_ui(i, j, self.get_unsafe_ui(i, j))
        if self._subdivisions is not None:
            mat.subdivide(self.subdivisions())
        return mat

    ########################################################################
    # LEVEL 2 helpers:
    #   These function support the implementation of the level 2 functionality.
    ########################################################################

    cpdef _add_(self, _right):
        """
        Addition

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: B = matrix(Zmod(36), 3, 3, range(2, 20, 2))
            sage: A + B
            [ 2  5  8]
            [11 14 17]
            [20 23 26]
        """
        cdef Matrix_modn_dense_flint right = _right
        cdef Matrix_modn_dense_flint M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_add(M._matrix, self._matrix, right._matrix)
        sig_off()
        return M

    cdef Matrix _matrix_times_matrix_(left, Matrix _right):
        """
        Multiplication

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 3, range(6))
            sage: B = matrix(Zmod(36), 3, 2, range(6))
            sage: A * B
            [10 13]
            [28  4]
        """
        if left._ncols != _right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right.")
        cdef Matrix_modn_dense_flint right = _right
        cdef Py_ssize_t i, j
        cdef Matrix_modn_dense_flint M = left._new(left._nrows, right._ncols)
        sig_on()
        nmod_mat_mul(M._matrix, left._matrix, right._matrix)
        sig_off()
        return M

    cpdef _lmul_(self, Element right):
        """
        Scalar multiplication

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(0, 18, 2))
            sage: 9*A
            [ 0 18  0]
            [18  0 18]
            [ 0 18  0]
        """
        cdef Matrix_modn_dense_flint M = self._new(self._nrows, self._ncols)
        e = self._modulus.element_class()
        cdef mp_limb_t ivalue
        if e is IntegerMod_int:
            ivalue = (<IntegerMod_int>right).ivalue
        elif e is IntegerMod_int64:
            ivalue = (<IntegerMod_int64>right).ivalue
        else:
            ivalue = mpz_get_ui((<IntegerMod_gmp>right).value)
        sig_on()
        nmod_mat_scalar_mul(M._matrix, self._matrix, ivalue)
        sig_off()
        return M


    def __copy__(self):
        """
        Return a copy of this matrix. Changing the entries of the copy will
        not change the entries of this matrix.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 2, range(4))
            sage: B = copy(A)
            sage: B[1,1] = 35
            sage: A
            [0 1]
            [2 3]
        """
        cdef Matrix_modn_dense_flint M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_set(M._matrix, self._matrix)
        sig_off()
        return M

    def __neg__(self):
        r"""
        Return the negative of this matrix.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 2, range(4))
            sage: -A
            [ 0 35]
            [34 33]
        """
        cdef Matrix_modn_dense_flint M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_neg(M._matrix, self._matrix)
        sig_off()
        return M

    cpdef _richcmp_(self, right, int op):
        r"""
        Compare ``self`` with ``right``, examining entries in
        lexicographic (row major) ordering.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 2, range(4))
            sage: B = matrix(Zmod(36), 2, 2, [0, 1, 3, 2])
            sage: A == A
            True
            sage: A != A
            False
            sage: A == B
            False
            sage: A != B
            True
            sage: A < B
            True
        """
        cdef Py_ssize_t i, j
        cdef int k

        if op == Py_EQ:
            return bool(nmod_mat_equal(self._matrix, (<Matrix_modn_dense_flint>right)._matrix))
        elif op == Py_NE:
            return not bool(nmod_mat_equal(self._matrix, (<Matrix_modn_dense_flint>right)._matrix))
        else:
            sig_on()
            for i in range(self._nrows):
                for j in range(self._ncols):
                    k = nmod_mat_get_entry(self._matrix,i,j) - nmod_mat_get_entry((<Matrix_modn_dense_flint>right)._matrix,i,j)
                    if k:
                        sig_off()
                        if k < 0:
                            return rich_to_bool(op, -1)
                        else:
                            return rich_to_bool(op, 1)
            sig_off()
            return rich_to_bool(op, 0)

    ########################################################################
    # LEVEL 3 helpers:
    #   These function support the implementation of the level 2 functionality.
    ########################################################################

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: bool(matrix(Zmod(36), 2, 2, 1))
            True
            sage: bool(matrix(Zmod(36), 0, 0))
            False
        """
        return not nmod_mat_is_zero(self._matrix)

    cpdef _sub_(self, _right):
        """
        Subtraction

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: B = matrix(Zmod(36), 3, 3, range(2, 20, 2))
            sage: A - B
            [34 33 32]
            [31 30 29]
            [28 27 26]
        """
        cdef Matrix_modn_dense_flint right = _right
        cdef Matrix_modn_dense_flint M = self._new(self._nrows, self._ncols)
        nmod_mat_sub(M._matrix, self._matrix, right._matrix)
        return M

    def __invert__(self):
        r"""
        Return the inverse of this matrix.

        EXAMPLES::

            sage: A = matrix(Zmod(36), [[28, 32, 19], [25, 24, 2], [15, 11, 30]])
            sage: B = ~A
            sage: A * B == B * A == identity_matrix(Zmod(36), 3, 3)
            True
        """
        return self.inverse_of_unit()

    def inverse_of_unit(self, algorithm=None):
        """
        Return the inverse of this matrix.

        Raises a ``ZeroDivisionError`` if the determinant is not a unit,
        and raises an ``ArithmeticError`` if the
        inverse doesn't exist because the matrix is nonsquare.

        INPUT:

        - ``algorithm`` -- either ``"crt"`` or ``"lift"``.
            In the first case, CRT is used to assemble results from prime powers
            dividing the modulus, with quadratic Hensel lifting for prime powers.
            In the second case, the inverse is computed over the integers and
            then reduced.

        EXAMPLES::

            sage: A = matrix(Zmod(36), [[28, 32, 19], [25, 24, 2], [15, 11, 30]])
            sage: A.inverse_of_unit()
            [14  5  4]
            [ 0 15 23]
            [23 28 16]

            sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]:
            ....:     algs = ["flint", "linbox", "lift"] if N == 5 else ["crt", "lift"]
            ....:     for n in [3, 6]:
            ....:         A = random_matrix(Zmod(N), n, n)
            ....:         while not A.det().is_unit():
            ....:             A = random_matrix(Zmod(N), n, n)
            ....:         Is = [A.inverse_of_unit(alg) for alg in algs]
            ....:         for I in Is: I.set_immutable()
            ....:         assert len(set(Is)) == 1

            sage: matrix(Zmod(36), [[2, 0], [0, 2]]).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
            sage: matrix(Zmod(36), [[2], [2]]).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: inverse only defined for square matrix
        """
        if not self.is_square():
            raise ArithmeticError("inverse only defined for square matrix")
        if not self.nrows():
            return self

        cdef Matrix_modn_dense_flint M
        R = self._parent._base
        cdef Py_ssize_t i, j, n = self._nrows
        cdef long k, e, b, nlifts
        cdef mp_limb_t p, N = 1
        cdef nmod_mat_t inv, A, B, tmp
        cdef bint lift_required, ok
        if algorithm is None:
            if R.is_field():
                if n > 1000:
                    from .matrix_modn_dense_double import MAX_MODULUS
                    if R.order() < MAX_MODULUS:
                        algorithm = "linbox"
                    else:
                        algorithm = "flint"
                else:
                    algorithm = "flint"
            else:
                algorithm = "crt"
        if algorithm == "flint":
            if not R.is_field():
                raise ValueError("Flint inversion only supported over fields; use algorithm=crt")
            M = self._new(n, n)
            ok = nmod_mat_inv(M._matrix, self._matrix)
            if not ok:
                raise ZeroDivisionError("input matrix must be nonsingular")
            return M
        elif algorithm == "linbox":
            if not R.is_field():
                raise ValueError("Linbox inversion only supported over fields; use algorithm=crt")
            return self._change_implementation("linbox").inverse_of_unit()._change_implementation("flint")
        elif algorithm == "lift":
            return (~self.lift_centered()).change_ring(R)
        elif algorithm == "crt":
            # The modulus is small, so factoring is feasible.
            # We find inverses modulo each prime dividing the modulus, lift p-adically, then CRT the results together
            M = self._new(n, n)
            F = R.factored_order()
            maxe = max(fac[1] for fac in F)
            lift_required = (maxe > 1)
            if lift_required:
                nlifts = Integer(maxe - 1).nbits()
                lift_mods = [1 for _ in range(nlifts)]
            try:
                nmod_mat_init(A, n, n, 1)
                nmod_mat_init(inv, n, n, 1)
                if lift_required:
                    nmod_mat_init(tmp, n, n, 1)
                for pz, ez in F:
                    p = pz
                    _nmod_mat_set_mod(A, p)
                    for i in range(n):
                        for j in range(n):
                            nmod_mat_set_entry(A, i, j, nmod_mat_get_entry(self._matrix, i, j) % p)
                    ok = nmod_mat_inv(A, A)
                    if not ok:
                        raise ZeroDivisionError("input matrix must be nonsingular")
                    _nmod_mat_set_mod(inv, N*p)
                    for i in range(n):
                        for j in range(n):
                            nmod_mat_set_entry(inv, i, j, n_CRT(nmod_mat_get_entry(inv, i, j), N, nmod_mat_get_entry(A, i, j), p))
                    N *= p
                    # Now inv is accurate modulo N
                    if lift_required:
                        e = ez
                        # Update the moduli for lifting so that the exponent of p at most doubles
                        for k in reversed(range(nlifts)):
                            lift_mods[k] *= p**e
                            e = (e + 1) // 2
                if lift_required:
                    for k in range(nlifts):
                        N = lift_mods[k]
                        _nmod_mat_set_mod(tmp, N)
                        _nmod_mat_set_mod(inv, N)
                        _nmod_mat_set_mod(A, N)
                        for i in range(n):
                            for j in range(n):
                                nmod_mat_set_entry(A, i, j, nmod_mat_get_entry(self._matrix, i, j) % N)
                        # Suppose B is an approximation to the inverse of A.
                        # I - AB = p^k E
                        # (I - AB)(I - AB) = I - A(2B - BAB) = p^(2k)E^2
                        # So 2B - BAB is a better approximation
                        nmod_mat_mul(tmp, A, inv)
                        nmod_mat_mul(tmp, inv, tmp)
                        nmod_mat_scalar_mul(inv, inv, 2)
                        nmod_mat_sub(inv, inv, tmp)
                        # Now inv is accurate mod N
                nmod_mat_set(M._matrix, inv)
                return M
            finally:
                nmod_mat_clear(A)
                nmod_mat_clear(inv)
                if lift_required:
                    nmod_mat_clear(tmp)
        raise ValueError("Unrecognized algorithm '%s'" % algorithm)

    def hessenberg_form(self):
        """
        The Hessenberg form of a matrix is almost upper triangular,
        and allows for efficient computation of the characteristic polynomial.

        Requires that the base ring have prime power modulus.

        EXAMPLES::

            sage: A = matrix(Zmod(125), 4, [2,1,1,-2,2,2,-1,-1,-1,1,2,3,4,5,6,7])
            sage: H = A.hessenberg_form(); H
            [  2  59  88   1]
            [  2  63  34 124]
            [  0  15 104   8]
            [  0   0  90  94]
            sage: H.charpoly() == A.charpoly()
            True
        """
        B = self.__copy__()
        B.hessenbergize()
        return B

    def hessenbergize(self):
        """
        Transform this matrix into Hessenberg form.

        The hessenberg form of a matrix `A` is a matrix that is
        similar to `A`, so has the same characteristic polynomial
        as `A`, and is upper triangular except possible for entries
        right below the diagonal.

        Requires that the base ring have prime power modulus.

        The algorithm is a modification of the standard one over fields,
        where rather than swapping in an arbitrary nonzero element in each column
        one instead picks an element of minimal valuation.

        EXAMPLES::

            sage: A = matrix(Zmod(125), 4, [2,1,1,-2,2,2,-1,-1,-1,1,2,3,4,5,6,7])
            sage: A.hessenbergize(); A
            [  2  59  88   1]
            [  2  63  34 124]
            [  0  15 104   8]
            [  0   0  90  94]

        You cannot Hessenbergize an immutable matrix::

            sage: A = matrix(Zmod(125), 3, range(9))
            sage: A.set_immutable()
            sage: A.hessenbergize()
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
        """
        R = self._parent._base
        F = R.factored_order()
        if len(F) != 1:
            raise ValueError("Hessenbergize only supported for prime power modulus")
        if not self.is_square():
            raise TypeError("self must be square")
        self.check_mutability()
        self.clear_cache()

        pz, ez = F[0]
        cdef Py_ssize_t i, j, k, m, n, r
        n = self._nrows

        cdef bint found_nonzero
        cdef mp_limb_t u, minu, inv, q, x, y, prepe, prepe1, preq, pe, p = pz
        cdef int v, minv, e = ez
        cdef double pinv = n_precompute_inverse(p)
        pe = n_pow(p, e)
        prepe = n_preinvert_limb(pe)
        prepe1 = n_preinvert_limb(n_pow(p, e-1))

        for m in range(1, n-1):
            i = -1
            minv = e
            for r in range(m, n):
                u = nmod_mat_get_entry(self._matrix, r, m-1)
                if u != 0:
                    v = n_remove2_precomp(&u, p, pinv)
                    if v < minv:
                        i = r
                        minv = v
                        minu = u
            if i != -1:
                # i is the first row below m with an entry in colunmn m-1 with minimal valuation
                q = n_pow(p, e-minv)
                if minv == 0:
                    preq = prepe
                elif minv == 1:
                    preq = prepe1
                else:
                    preq = n_preinvert_limb(q)
                inv = n_invmod(minu, q)
                if i > m:
                    self.swap_rows_c(i, m)
                    # We must do the corresponding column swap to
                    # maintain the characteristic polynomial (which is
                    # an invariant of Hessenberg form)
                    self.swap_columns_c(i, m)
                # Now the nonzero entry in position (m,m-1) is minu * p^minv.
                # Use it to clear the entries in column m-1 below m.
                for j in range(m+1, n):
                    u = nmod_mat_get_entry(self._matrix, j, m-1)
                    if u != 0:
                        # multiply by the inverse of the (m,m-1) entry
                        if minv == 0:
                            # Simpler case if minv=0
                            u = n_mulmod2_preinv(u, inv, q, preq)
                        else:
                            # extract powers of p from u, multiply by the unit
                            # and put the right number of powers of p back.
                            v = n_remove2_precomp(&u, p, pinv)
                            u = n_pow(p, v-minv) * n_mulmod2_preinv(u, inv, q, preq)
                        # self.add_multiple_of_row(j, m, -u, 0)
                        for k in range(n):
                            x = nmod_mat_get_entry(self._matrix, j, k)
                            y = nmod_mat_get_entry(self._matrix, m, k)
                            y = n_mulmod2_preinv(y, u, pe, prepe)
                            nmod_mat_set_entry(self._matrix, j, k, n_submod(x, y, pe))
                        # To maintain charpoly, do the corresponding column operation,
                        # which doesn't mess up the matrix, since it only changes
                        # column m, and we're only worried about column m-1 right now.
                        # Add u*column_j to column_m.
                        for k in range(n):
                            x = nmod_mat_get_entry(self._matrix, k, m)
                            y = nmod_mat_get_entry(self._matrix, k, j)
                            y = n_mulmod2_preinv(y, u, pe, prepe)
                            nmod_mat_set_entry(self._matrix, k, m, n_addmod(x, y, pe))

    def charpoly(self, var='x', algorithm=None):
        """
        Return the characteristic polynomial of this matrix, as a polynomial over the base ring.

        INPUT:

        - ``var`` -- a string, the variable name for the polynomial returned.

        - ``algorithm`` -- either ``"flint"``, ``"lift"``, ``"crt"``, or ``"linbox"``.
            In the first case, use FLINT's characteristic polynomial function,
            in the second case, lift to the integers and compute there,
            in the third form, compute the hessenberg form for each prime power
            dividing the modulus.
            In the last, converts to a linbox matrix and computes the charpoly there.
            If not given, defaults to ``"crt"`` when not over a field, to ``"linbox"`` when
            the dimension is large, and ``"flint"``  otherwise.

        EXAMPLES::

            sage: A = matrix(Zmod(36), [[28, 32, 19], [25, 24, 2], [15, 11, 30]])
            sage: A.charpoly()
            x^3 + 26*x^2 + 9*x + 35
            sage: A = matrix(Zmod(43^10), 6, [0, 5664461354126771, 12212357361910300, 15947020959157478, 0, 16792952041597449, 14690359073749623, 11237259451999884, 5117434014428142, 15157488677243483, 9004103062307752, 20761679499270441, 4620722392655416, 5445142895231681, 6605357538252496, 7608812697273777, 18542817615638637, 18194689690271501, 0, 20341333098836812, 12117922812876054, 1270149214447437, 0, 10999401748338075, 4620722392655416, 10891113386038365, 956055025271903, 2162842206467093, 18542817615638637, 1143972982339214, 13128267973348003, 15817056104759912, 20531311511260484, 13598045280630823, 7585782589268305, 14053895308766769])
            sage: A.charpoly()
            x^6 + 13124967810747524*x^5 + 20067912494391006*x^4 + 11204731077775359*x^3

            sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]:
            ....:     for n in [3, 6]:
            ....:         A = random_matrix(Zmod(N), n, n)
            ....:         assert A.charpoly(algorithm='crt') == A.charpoly(algorithm='lift')

            sage: A = random_matrix(GF(73), 6)
            sage: fs = [A.charpoly(algorithm=alg) for alg in ['flint', 'linbox', 'crt', 'lift', 'hessenberg', 'df']]
            sage: len(set(fs))
            1
            sage: A = random_matrix(Zmod(1728), 6)
            sage: fs = [A.charpoly(algorithm=alg) for alg in ['crt', 'lift', 'df']]
            sage: len(set(fs))
            1
        """
        if not self.is_square():
            raise TypeError("self must be square")
        ans = self.fetch('charpoly')
        if ans is not None:
            return ans.change_variable_name(var)
        R = self._parent._base
        F = R.factored_order()
        if algorithm is None:
            if len(F) == 1 and F[0][1] == 1: # N prime
                if self._nrows > 200:
                    from .matrix_modn_dense_double import MAX_MODULUS
                    if R.order() < MAX_MODULUS:
                        algorithm = "linbox"
                    else:
                        algorithm = "flint"
                else:
                    algorithm = "flint"
            else:
                algorithm = "crt"
        cdef Polynomial_zmod_flint f
        cdef Py_ssize_t i, j, m, jstart, n = self._nrows
        cdef mp_limb_t pe, prepe, preN, x, y, scalar, t, N = R.order()
        cdef Matrix_modn_dense_flint A
        cdef nmod_mat_t H, c
        S = PolynomialRing(R, var, implementation="FLINT")
        if algorithm == "linbox":
            ans = self._change_implementation("linbox").charpoly(var)
        elif algorithm == "flint":
            f = S()
            sig_on()
            nmod_mat_charpoly(&f.x, self._matrix)
            sig_off()
            ans = f
        elif algorithm == "lift":
            ans = self.lift_centered().charpoly(var).change_ring(R)
        elif algorithm == "crt":
            preN = n_preinvert_limb(N)
            # We hessenbergize at each prime, CRT, then compute the charpoly from the Hessenberg form
            try:
                nmod_mat_init(H, n, n, 1)
                nmod_mat_init(c, n+1, n+1, N)
                # Now we use N cummulatively as we work through the factorization
                N = 1
                for pz, ez in F:
                    pe = n_pow(pz, ez)
                    prepe = n_preinvert_limb(pe)
                    R = Zmod(pe)
                    R.factored_order.set_cache(Factorization([(pz, ez)]))
                    from sage.matrix.matrix_space import MatrixSpace
                    A = MatrixSpace(R, self._nrows, self._ncols, sparse=False, implementation="flint")()
                    for i in range(n):
                        for j in range(n):
                            x = nmod_mat_get_entry(self._matrix, i, j)
                            x = n_mod2_preinv(x, pe, prepe)
                            nmod_mat_set_entry(A._matrix, i, j, x)
                    A.hessenbergize()
                    _nmod_mat_set_mod(H, N*pe)
                    for i in range(n):
                        jstart = i - 1
                        if i == 0:
                            jstart = 0
                        for j in range(jstart, n):
                            nmod_mat_set_entry(H, i, j, n_CRT(nmod_mat_get_entry(H, i, j), N, nmod_mat_get_entry(A._matrix, i, j), pe))
                    N *= pe
                # We represent the intermediate polynomials that come up in
                # the calculations as rows of an (n+1)x(n+1) matrix, since
                # we've implemented basic arithmetic with such a matrix.
                # Also see Cohen's first GTM, Algorithm 2.2.9.
                nmod_mat_set_entry(c, 0, 0, 1)
                for m in range(1, n+1):
                    # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
                    # We do this by hand by setting the m-th row to c[m-1]
                    # shifted to the right by one.  We then add
                    # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
                    for i in range(1, n+1):
                        nmod_mat_set_entry(c, m, i, nmod_mat_get_entry(c, m-1, i-1))
                    # self.sub_multiple_of_row(m, m-1, -H[m-1,m-1], 0)
                    scalar = nmod_mat_get_entry(H, m-1, m-1)
                    for j in range(n+1):
                        x = nmod_mat_get_entry(c, m, j)
                        y = nmod_mat_get_entry(c, m-1, j)
                        y = n_mulmod2_preinv(y, scalar, N, preN)
                        nmod_mat_set_entry(c, m, j, n_submod(x, y, N))
                    t = 1
                    for i in range(1, m):
                        t = n_mulmod2_preinv(t, nmod_mat_get_entry(H, m-i, m-i-1), N, preN)
                        scalar = n_mulmod2_preinv(t, nmod_mat_get_entry(H, m-i-1, m-1), N, preN)
                        for j in range(n+1):
                            x = nmod_mat_get_entry(c, m, j)
                            y = nmod_mat_get_entry(c, m-i-1, j)
                            y = n_mulmod2_preinv(y, scalar, N, preN)
                            nmod_mat_set_entry(c, m, j, n_submod(x, y, N))
                # the answer is now the n-th row of c.
                f = S()
                for j in range(n+1):
                    nmod_poly_set_coeff_ui(&f.x, j, nmod_mat_get_entry(c, n, j))
                return f
            finally:
                nmod_mat_clear(H)
                nmod_mat_clear(c)
        elif algorithm == "hessenberg":
            ans = self._charpoly_hessenberg(var)
        elif algorithm == "df":
            ans = self._charpoly_df(var)
        else:
            raise ValueError("Unknown algorithm '%s'" % algorithm)
        self.cache('charpoly', ans)
        return ans

    cpdef _shift_mod(self, mp_limb_t modulus, mp_limb_t shift=1, bint mul=True, bint domod=True):
        """
        A fast method for returning a copy of this matrix with
        different modulus or scaled by an integer.

        It is the caller's responsibility to ensure that the combination of
        shift and modulus yields a result that is mathematically meaningful.

        INPUT:

        - ``modulus`` -- an integer, the modulus for the result.

        - ``shift`` -- an integer to multiply or divide by (default 1)

        - ``mul`` -- boolean.  If true, multiply by ``shift``, otherwise divide.
            In the case of division, floor division is used; the intended purpose
            is for when all entries are divisible by ``shift`` and the modulus
            will be divided as well.

        - ``domod`` -- boolean.  If false, assumes that the entries
            don't need to be reduced after scaling.

        EXAMPLES::

            sage: A = matrix(Zmod(12), 2, [2,4,6,8])
            sage: A._shift_mod(6, 2, mul=False, domod=False)
            [1 2]
            [3 4]
            sage: A._shift_mod(24, 2, mul=True, domod=False)
            [ 4  8]
            [12 16]
            sage: A._shift_mod(12, 3, mul=True)
            [6 0]
            [6 0]

        If not all entries are divisible by shift, the result may not be sensible::

            sage: A._shift_mod(4, 3, mul=False)
            [0 1]
            [2 2]

        It's possible to get unreduced matrices::

            sage: A._shift_mod(12, 2, domod=False)
            [ 4  8]
            [12 16]

        Note that you also need to be careful about not reducing if there's a
        possibility of overflow::

            sage: N = 2^63 - 5
            sage: A = matrix(Zmod(N), 1, [-1])
            sage: A._shift_mod(N, 3, domod=False)
            [9223372036854775790]
            sage: (3*(N-1)) % (2^64)
            9223372036854775790
            sage: A.lift() * 3
            [27670116110564327406]
            sage: (A.lift() * 3).change_ring(Zmod(N))
            [9223372036854775800]

        This works correctly if you use reduction::

            sage: A._shift_mod(N, 3, domod=True)
            [9223372036854775800]
        """
        cdef Py_ssize_t i, j, m = self._nrows, n = self._ncols
        cdef mp_limb_t x, preN, preshift, N = modulus
        preN = n_preinvert_limb(N)
        if not mul:
            preshift = n_preinvert_limb(shift)
        R = Zmod(N)
        from sage.matrix.matrix_space import MatrixSpace
        P = MatrixSpace(R, m, n, sparse=False)
        cdef Matrix_modn_dense_flint ans = Matrix_modn_dense_flint.__new__(Matrix_modn_dense_flint, P)
        ans._parent = P
        for i in range(m):
            for j in range(n):
                x = nmod_mat_get_entry(self._matrix, i, j)
                if shift == 1:
                    if domod:
                        x = n_mod2_preinv(x, N, preN)
                elif mul:
                    if domod:
                        x = n_mulmod2_preinv(x, shift, N, preN)
                    else:
                        x *= shift
                else:
                    x = n_div2_preinv(x, shift, preshift)
                    if domod:
                        x = n_mod2_preinv(x, N, preN)
                nmod_mat_set_entry(ans._matrix, i, j, x)
        return ans

    def minpoly_ideal(self, var='x', proof=None, **kwds):
        """
        The ideal of polynomials over the base ring that vanish on this matrix.

        When the base ring is not a field, this ideal is not necessarily principal.

        INPUT:

        - ``var`` -- the variable name for the polynomial ring

        - ``proof`` -- boolean, whether to check that the answer computed by using a minimal polynomial is correct.  If not specificied, uses the linear alegbra default proof state.  Note that if the modulus is composite and divisible by small primes the probability of an incorrect result is substantial.  The result is not cached when proof is ``False``, so this function can be called multiple times to get a desired level of certainty.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, [6, 0, 0, -6])
            sage: A.minpoly_ideal()
            Ideal (x^2, 3*x + 18) of Univariate Polynomial Ring in x over Ring of integers modulo 36

            sage: A = matrix(Zmod(43^10), 6, [0, 5664461354126771, 12212357361910300, 15947020959157478, 0, 16792952041597449, 14690359073749623, 11237259451999884, 5117434014428142, 15157488677243483, 9004103062307752, 20761679499270441, 4620722392655416, 5445142895231681, 6605357538252496, 7608812697273777, 18542817615638637, 18194689690271501, 0, 20341333098836812, 12117922812876054, 1270149214447437, 0, 10999401748338075, 4620722392655416, 10891113386038365, 956055025271903, 2162842206467093, 18542817615638637, 1143972982339214, 13128267973348003, 15817056104759912, 20531311511260484, 13598045280630823, 7585782589268305, 14053895308766769])
            sage: I = A.minpoly_ideal(); I
            Ideal (x^4 + 3889828461*x^3 + 5524445805550*x^2 + 182758215997616*x, 6321363049*x^3 + 1630911666642*x^2 + 140258403331212*x, 11688200277601*x^2, 502592611936843*x) of Univariate Polynomial Ring in x over Ring of integers modulo 21611482313284249

        We check that all generators vanish on A::

            sage: all(f(A) == 0 for f in I.gens())
            True

        We check that the ideal is preserved by conjugation::

            sage: B = random_matrix(Zmod(43^10), 6)
            sage: while not B.is_unit():
            ....:     B = random_matrix(Zmod(43^10), 6)
            sage: J = (~B * A * B).minpoly_ideal()
            sage: J == I
            True

        We check that the polynomials for random matrices vanish for various moduli::

            sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]:
            ....:     for n in [3, 6]:
            ....:         A = random_matrix(Zmod(N), n, n)
            ....:         I = A.minpoly_ideal()
            ....:         assert all(f(A) == 0 for f in I.gens())
            ....:         assert I.gens()[0].is_monic()
        """
        if not self.is_square():
            raise ValueError("Minimal polynomial ideal not defined for non-square matrices")
        if proof is None:
            proof = linalg_proof()
        R = self.base_ring()
        Rx = PolynomialRing(R, var, implementation="FLINT")
        key = "minpoly_ideal"
        ans = self.fetch(key)
        if ans is not None:
            if ans.ring().variable_name() != var:
                ans = Rx.ideal([Rx(g) for g in ans.gens()])
            return ans
        n = self.nrows()
        F = R.factored_order()
        P = self.parent()
        cdef mp_limb_t piv, c, N = self.base_ring().order()
        cdef Matrix_modn_dense_flint C
        cdef Polynomial_zmod_flint f, mpoly
        cdef Py_ssize_t i, j, jj, k, d
        if len(F) == 1:
            p, e = F[0]
            # We start with the minimal polynomial mod p
            S = Zmod(p)
            Sx = PolynomialRing(S, var, implementation="FLINT")
            from sage.matrix.matrix_space import MatrixSpace
            SM = MatrixSpace(S, self._nrows, self._ncols, implementation="flint")
            C = SM(self)
            mpoly = Sx()
            sig_on()
            nmod_mat_minpoly(&mpoly.x, C._matrix)
            sig_off()
            if e == 1:
                ans = Rx.ideal(mpoly)
                self.cache(key, ans)
                return ans
            d = mpoly.degree()
            if d == n:
                # The minimal polynomial is the same as the characteristic polynomial
                ans = Rx.ideal(self.charpoly(var=var))
                self.cache(key, ans)
                return ans
        # We pick a random vector and iteratively multiply by the matrix n times
        cdef Matrix_modn_dense_flint v = self._new(n, 1), B = self._new(n+1, n+1)
        pows = [self.parent().identity_matrix()]
        def check(polys, pows):
            d = polys[0].degree()
            while d >= len(pows):
                pows.append(pows[-1] * self)
            for f in polys:
                C = self._new(n, n)
                for k in range(f.degree() + 1):
                    C += f[k] * pows[k]
                if C:
                    return False
            return True
        while True:
            for i in range(n):
                c = randrange(N)
                nmod_mat_set_entry(v._matrix, i, 0, c)
                nmod_mat_set_entry(B._matrix, i, n, c)
            # It would be nice to do the kernel computation in parallel with computing self*v so that we know if we can stop early
            for j in range(n-1, -1, -1):
                v = self * v
                for i in range(n):
                    nmod_mat_set_entry(B._matrix, i, j, nmod_mat_get_entry(v._matrix, i, 0))
            # _right_kernel_matrix caches echelon form, and we're changing B
            B.clear_cache()
            _, C = B._right_kernel_matrix()
            C.howellize()
            polys = []
            # scan for last nonzero row
            for i in range(n, -1, -1):
                for j in range(n):
                    piv = nmod_mat_get_entry(C._matrix, i, j)
                    if piv:
                        f = Rx()
                        nmod_poly_fit_length(&f.x, n+1-j)
                        for jj in range(n+1-j):
                            nmod_poly_set_coeff_ui(&f.x, jj, nmod_mat_get_entry(C._matrix, i, n-jj))
                        polys.append(f)
                        break
                if piv == 1:
                    # reached the monic generator, so try to return an answer
                    polys.reverse()
                    # We only know that the polynomials vanish on a random vector.
                    # We need to check that they actually vanish on the matrix
                    if proof and check(polys, pows) or not proof:
                        ans = Rx.ideal(polys)
                        if proof: self.cache(key, ans)
                        return ans
                    break

    def minpoly(self, var='x', proof=None, **kwds):
        """
        Returns the minimal polynomial of this matrix.

        Note that, in general, the minimal polynomial of a matrix over a ring with composite modulus is not well defined.

        .. SEEALSO::

            :meth:`minpoly_ideal`

        INPUT:

        - ``var`` -- the variable name for the polynomial ring

        - ``proof`` -- boolean, whether to check that the answer computed by using a minimal polynomial is correct.  If not specificied, uses the linear alegbra default proof state.  Note that if the modulus is composite and divisible by small primes the probability of an incorrect result is substantial.  The result is not cached when proof is ``False``, so this function can be called multiple times to get a desired level of certainty.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 6, 6, 1)
            sage: A.minpoly()
            x + 35

            sage: A = matrix(Zmod(43^4), 4, 4, [1547380, 1211593, 1900136, 2872531, 3010350, 3130824, 2590437, 420249, 1236104, 1628956, 1300526, 145927, 1817994, 1412525, 1305313, 2775020])
            sage: f = A.minpoly(); f
            x^4 + 1502653*x^3 + 2556099*x^2 + 595635*x + 865202
            sage: f(A) == 0
            True

            sage: A = matrix(Zmod(36), 2, [6, 0, 0, -6])
            sage: A.minpoly()
            Traceback (most recent call last):
            ...
            ValueError: Matrix does not have a minimal polynomial; try minpoly_ideal

        We build matrices with a given minimal polynomial, then check that the correct result is returned::

           sage: R.<x> = Zmod(2^63-1)[]
           sage: f = x^3 + R.random_element(2)
           sage: g = x^2 + R.random_element(1)
           sage: A = block_diagonal_matrix([companion_matrix(f), companion_matrix(g), companion_matrix(f*g)])
           sage: B = random_matrix(Zmod(2^63-1), 10, 10)
           sage: while not B.is_unit():
           ....:     B = random_matrix(Zmod(2^63-1), 10, 10)
           sage: C = ~B * A * B
           sage: C.minpoly() == f*g
           True
        """
        I = self.minpoly_ideal(var, proof=proof)
        gens = I.gens()
        if len(gens) != 1:
            raise ValueError("Matrix does not have a minimal polynomial; try minpoly_ideal")
        return gens[0]

    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2):
        """
        Swap two columns using FLINT.
        """
        nmod_mat_swap_cols(self._matrix, NULL, c1, c2)

    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2):
        """
        Swap two rows using FLINT.
        """
        nmod_mat_swap_rows(self._matrix, NULL, r1, r2)

    def _solve_right_modn(self, B, check=True):
        """
        Solve the matrix equation `A X = B`.

        sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]: # indirect doctest
        ....:     for n in [3, 6]:
        ....:         for m in [4, 5, 6]:
        ....:             A = random_matrix(Zmod(N), m, n)
        ....:             for k in [1, 2, 4]:
        ....:                 X0 = random_matrix(Zmod(N), n, k)
        ....:                 B = A * X0
        ....:                 X = A.solve_right(B)
        ....:                 assert A * X == B
        """
        b_is_vec = is_Vector(B)
        cdef Py_ssize_t i, j, ii, jj, n, Cm, An = self.ncols()
        cdef mp_limb_t piv, q, r, tmp, entry, N = self.base_ring().order(), Ninv = n_preinvert_limb(N)
        cdef double pivinv
        if b_is_vec:
            n = 1
        else:
            n = B.ncols()
        cdef Matrix_modn_dense_flint C, X
        R = self.base_ring()
        X = self._new(An, n)
        if R.is_field():
            C = B.column() if b_is_vec else B
            sig_on()
            nmod_mat_can_solve(X._matrix, self._matrix, C._matrix)
            sig_off()
        else:
            # We use Howell form to solve the system
            C = self.augment(B)
            if C.nrows() < C.ncols():
                C = C.stack(C.new_matrix(nrows=C.ncols() - C.nrows()))
            C.echelonize()
            # save the echelon form if not already cached
            if self.fetch('echelon_form') is None:
                E = C.submatrix(max(self.ncols(), self.nrows()), self.ncols())
                self.cache('echelon_form', E)
            Cm = C.nrows()
            for i in range(Cm - 1, -1, -1):
                for j in range(An):
                    piv = nmod_mat_get_entry(C._matrix, i, j)
                    if not piv:
                        continue
                    elif piv != 1:
                        pivinv = n_precompute_inverse(piv)
                    for jj in range(n):
                        if piv == 1:
                            q = nmod_mat_get_entry(C._matrix, i, An + jj)
                        else:
                            entry = nmod_mat_get_entry(C._matrix, i, An + jj)
                            r = n_divrem2_precomp(&q, entry, piv, pivinv)
                            if r:
                                raise ValueError("matrix equation has no solution")
                        nmod_mat_set_entry(X._matrix, j, jj, q)
                        # We now update the right augmented block based on the entries above the pivot.
                        for ii in range(i - 1, -1, -1):
                            entry = n_mulmod2_preinv(q, nmod_mat_get_entry(C._matrix, ii, j), N, Ninv)
                            entry = n_submod(nmod_mat_get_entry(C._matrix, ii, An + jj), entry, N)
                            nmod_mat_set_entry(C._matrix, ii, An + jj, entry)
                    break
                else:
                    # full zero row on the left augmented block; check that it's also zero on the right
                    for jj in range(n):
                        if nmod_mat_get_entry(C._matrix, i, An + jj):
                            raise ValueError("matrix equation has no solution")
        if b_is_vec:
            # Convert back to a vector
            return X.column(0)
        else:
            return X

    def transpose(self):
        """
        Return the transpose of this matrix, without changing this matrix.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 2, [0, 1, 0, 0])
            sage: A.transpose()
            [0 0]
            [1 0]
        """
        cdef Matrix_modn_dense_flint M = self._new(self._ncols, self._nrows)
        sig_on()
        nmod_mat_transpose(M._matrix, self._matrix)
        sig_off()
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            M.subdivide(col_divs, row_divs)
        return M


    def echelonize(self, algorithm="default", **kwds):
        """
        Transform this matrix into echelon form in place, over the same base ring.

        .. NOTE::

            Over a field, transforms to standard reduced row echelon form.
            Otherwise, uses Howell form, which may increase the number of nonzero rows.  In this case, we require that the number of rows is at least the number of columns.

        INPUT:

        - ``algorithm`` -- a string, either "default", "flint" or "linbox"

        - ``transformation`` -- boolean. Whether to return a matrix `T` so that left multiplication by `T` transforms the original matrix into echelon form.

        OUTPUT:

        This matrix is put into echelon form.  Nothing is returned unless
        the keyword option ``transformation=True`` is specified, in
        which case the transformation matrix is returned.

        EXAMPLES::

            sage: A = matrix(Zmod(625), 4, 3, [[404, 355, 133], [375, 482, 448], [506, 115,  77], [370, 384, 66]])
            sage: A
            [404 355 133]
            [375 482 448]
            [506 115  77]
            [370 384  66]
            sage: A.echelonize()
            sage: A
            [1 0 2]
            [0 1 4]
            [0 0 5]
            [0 0 0]
        """
        self.check_mutability()
        self.clear_cache()
        cdef bint transformation = 'transformation' in kwds and kwds['transformation']
        cdef Matrix_modn_dense_flint aug, trans, E
        cdef Py_ssize_t i, j, m, n
        R = self._parent._base
        if algorithm == "default":
            if R.is_field() and min(self._nrows, self._ncols) > 200:
                from .matrix_modn_dense_double import MAX_MODULUS
                if R.order() < MAX_MODULUS:
                    algorithm = "linbox"
                else:
                    algorithm = "flint"
            else:
                algorithm = "flint"
        if algorithm not in ["flint", "linbox"]:
            raise ValueError("Unknown algorithm '%s'" % algorithm)
        if algorithm == "linbox":
            transformation = False
            E = self._change_implementation("linbox").echelon_form(algorithm="linbox_noefd")._change_implementation("flint")
            nmod_mat_set(self._matrix, E._matrix)
        elif self._parent._base.is_field():
            if transformation:
                m = self.nrows()
                n = self.ncols()
                from sage.matrix.special import identity_matrix
                aug = self.augment(identity_matrix(self.base_ring(), m))
                trans = self._new(m, m)
                sig_on()
                nmod_mat_rref(aug._matrix)
                sig_off()
                for i in range(m):
                    for j in range(n):
                        nmod_mat_set_entry(self._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, j))
                    for j in range(m):
                        nmod_mat_set_entry(trans._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, n+j))
            else:
                sig_on()
                rank = nmod_mat_rref(self._matrix)
                sig_off()
                self.cache('rank', Integer(rank))
        else:
            trans = self.howellize(transformation)
        self.cache('in_echelon_form', True)
        if transformation:
            return trans

    def _echelon_copy(self):
        """
        Copies the matrix and adds zero rows at the bottom if necessary.

        EXAMPLES::

            sage: A = matrix(Zmod(5), [1, 2])
            sage: A._echelon_copy()
            [1 2]
            sage: B = matrix(Zmod(6), [1, 2])
            sage: B._echelon_copy()
            [1 2]
            [0 0]
        """
        if self._nrows < self._ncols and not self.base_ring().is_field():
            M = self._new(self._ncols - self._nrows, self._ncols)
            return self.stack(M)
        else:
            return self.__copy__()

    def _pivots(self):
        """
        When not over a field, it is possible to have leading entries
        that are zero divisors.  This function distinguishes
        between such pivots and pivots that are 1.

        OUTPUT:

        - a list of columns containing a leading 1

        - a list of columns containing a leading nonzero zero divisor

        EXAMPLES::

            sage: R = Zmod(625)
            sage: A = matrix(Zmod(625), 3, 4, [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25])
            sage: A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            sage: A._pivots()
            ((0,), (1, 3))
        """
        key = '_pivots'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        cdef Matrix_modn_dense_flint E
        cdef Py_ssize_t i, j, k = 0
        # howell form has all the zero rows at the bottom
        if self.fetch('in_echelon_form'):
            E = self
        else:
            E = self.echelon_form()
        p = []
        zdp = []
        for i in range(E._nrows):
            for j in range(k, E._ncols):
                if nmod_mat_get_entry(E._matrix, i, j) != 0:  # nonzero position
                    if nmod_mat_get_entry(E._matrix, i, j) == 1:
                        p.append(j)
                    else:
                        zdp.append(j) # is a zero divisor
                    k = j+1  # so start at next position next time
                    break
        ans = (tuple(p), tuple(zdp))
        self.cache(key, ans)
        return ans

    def pivots(self):
        """
        Return the columns containing a leading 1 in the echelon form of this matrix.

        When the base ring is not a field, there may be other rows with leading entries
        a zero divisor.  These columns are available using the :meth:`_pivots` method.

        EXAMPLES::

            sage: R = Zmod(625)
            sage: A = matrix(Zmod(625), 3, 4, [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25])
            sage: A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            sage: A.pivots()
            (0,)
        """
        return self._pivots()[0]

    def rank(self):
        """
        Return the rank of this matrix.

        When the base ring is not a field, this is defined as the number of leading 1 pivots.
        This choice has the benefit that a square matrix will be invertible exactly when
        the rank is the same as the number of rows.

        .. SEEALSO:

        - :meth:`pivots`

        - :meth:`_pivots`

        EXAMPlES::

            sage: entries = [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25]
            sage: A = matrix(GF(17), 3, 4, entries)
            sage: A
            [1 2 3 4]
            [0 5 5 6]
            [0 0 0 8]
            sage: A.rank()
            3
            sage: A = matrix(Zmod(625), 3, 4, entries)
            sage: A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            sage: A.rank()
            1
        """
        key = 'rank'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        if not self._parent._base.is_field():
            ans = len(self.pivots())
        else:
            p = self.fetch('pivots')
            if ans is not None:
                ans = len(p[0])
            else:
                sig_on()
                ans = nmod_mat_rank(self._matrix)
                sig_off()
        self.cache(key, ans)
        return ans


    def determinant(self, algorithm=None):
        """
        Return the determinant of this matrix.

        INPUT:

        - ``algorithm`` -- a string, one of ``flint``, ``charpoly``, ``linbox`` or ``lift``

        EXAMPLES::

            sage: A = matrix(ZZ, 3, [1, 6, 11, -2, 3, 1, 3, 1, 0])
            sage: all(A.change_ring(Zmod(N)).det() == A.det() for N in range(2, 100))
            True

            sage: A = random_matrix(Zmod(36), 6)
            sage: all(A.det()^n == (A^n).det() for n in range(2, 10))
            True
            sage: A.det("charpoly") == A.det("lift")
            True
            sage: B = random_matrix(Zmod(36), 6)
            sage: while not B.is_unit():
            ....:     B = random_matrix(Zmod(36), 6)
            sage: A.det() == (~B * A * B).det()
            True

            sage: A = random_matrix(Zmod(73), 6)
            sage: a = A.det("flint")
            sage: b = A.det("linbox")
            sage: c = A.det("lift")
            sage: d = A.det("charpoly")
            sage: a == b == c == d
            True
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        key = 'det'
        d = self.fetch(key)
        if d is not None:
            return d
        R = self._parent._base
        if algorithm is None:
            if R.is_field():
                if self._nrows > 400:
                    from .matrix_modn_dense_double import MAX_MODULUS
                    if R.order() < MAX_MODULUS:
                        algorithm = "linbox"
                    else:
                        algorithm = "flint"
                else:
                    algorithm = "flint"
            else:
                algorithm = "charpoly"

        # If charpoly known, then det is easy.
        f = self.fetch('charpoly')
        if f is not None:
            d = f[0]
            if self._ncols % 2:
                d = -d
        elif algorithm == "flint":
            sig_on()
            d = nmod_mat_det(self._matrix)
            sig_off()
        elif algorithm == "linbox":
            d = self._change_implementation("linbox").det()
        elif algorithm == "charpoly":
            f = self.charpoly()
            d = f[0]
            if self._ncols % 2:
                d = -d
        elif algorithm == "lift":
            d = self.lift_centered().det()
        else:
            raise ValueError("Unknown algorithm '%s'" % algorithm)

        d = R(d)
        self.cache(key, d)
        return d

    def trace(self):
        """
        Return the trace of this matrix, which is the sum of the diagonal entries.

        The input must be square.

        EXAMPLES::

            sage: A = matrix(ZZ, 3, [1, 6, 11, -2, 3, 1, 3, 1, 0])
            sage: all(A.change_ring(Zmod(N)).trace() == A.trace() for N in range(2, 100))
            True

        If the characteristic polynomial is known, negating the second coefficient
        is used to recover the trace::

            sage: A = matrix(Zmod(36), 3, [12, 12, 10, 6, 23, 24, 35, 24, 30])
            sage: A.charpoly()
            x^3 + 7*x^2 + 4*x + 22
            sage: A.trace()
            29
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        key = 'trace'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        f = self.fetch('charpoly')
        if f is not None:
            ans = -f[self._nrows - 1]
        else:
            sig_on()
            ans = nmod_mat_trace(self._matrix)
            sig_off()
            ans = self.base_ring()(ans)
        self.cache(key, ans)
        return ans

    def strong_echelonize(self):
        """
        Transform this matrix in place into its strong echelon form.

        The strong echelon form obtained from the Howell form by permuting rows:
        rather than having all zero rows at the bottom,
        they are placed so that the pivots occur on the diagonal.

        The input must have at least as many rows as columns.

        EXAMPLES::

            sage: entries = [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25, 0, 0, 0, 0]
            sage: A = matrix(Zmod(625), 4, entries); A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            [ 0  0  0  0]
            sage: A.strong_echelonize(); A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0  0]
            [ 0  0  0 25]
        """
        if self._nrows < self._ncols:
            raise ValueError("self must must have at least as many rows as columns.")
        self.check_mutability()
        self.clear_cache()
        sig_on()
        nmod_mat_strong_echelon_form(self._matrix)
        sig_off()

    def strong_echelon_form(self):
        """
        Strong echelon form of this matrix.

        This is obtained from the Howell form (the form used
        as echelon form when not over a field) by permuting rows:
        rather than having all zero rows at the bottom,
        they are placed so that the pivots occur on the diagonal.

        The result will always have at least as many rows as columns, with zero
        rows added if necessary to accomplish this.

        EXAMPLES::

            sage: entries = [624, 353, 352, 497, 182, 374, 556, 365, 271, 557, 203, 327]
            sage: A = matrix(Zmod(625), 3, 4, entries)
            sage: A.strong_echelon_form()
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0  0]
            [ 0  0  0 25]
        """
        key = 'strong_echelon_form'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        if self._nrows < self._ncols:
            M = self._new(self._ncols - self._nrows, self._ncols)
            ans = self.stack(M)
        else:
            ans = self.__copy__()
        ans.strong_echelonize()
        self.cache(key, ans)
        return ans

    def howellize(self, transformation=False):
        r"""
        Transform this matrix in place into its Howell form.

        See :meth:`howell_form` for the definition.

        Since the Howell form may have more rows than `A`, to transform `A` in place
        we require that it has at least as many rows as columns.

        EXAMPLES::

            sage: entries = [624, 353, 352, 497, 182, 374, 556, 365, 271, 557, 203, 327, 181, 102, 283, 237]
            sage: A = matrix(Zmod(625), 4, entries)
            sage: A.howellize(); A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            [ 0  0  0  0]
        """
        if self._nrows < self._ncols:
            raise ValueError("matrix must have at least as many rows as columns.")
        self.check_mutability()
        cdef Matrix_modn_dense_flint aug, trans
        cdef Py_ssize_t i, j, m, n
        self.clear_cache()
        if transformation:
            m = self.nrows()
            n = self.ncols()
            from sage.matrix.special import identity_matrix
            aug = self.augment(identity_matrix(self.base_ring(), m)).stack(self._new(n, m + n))
            trans = self._new(m, m)
            sig_on()
            nmod_mat_howell_form(aug._matrix)
            sig_off()
            for i in range(m):
                for j in range(n):
                    nmod_mat_set_entry(self._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, j))
                for j in range(m):
                    nmod_mat_set_entry(trans._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, n+j))
            return trans
        else:
            sig_on()
            nmod_mat_howell_form(self._matrix)
            sig_off()

    def howell_form(self):
        """
        Return the Howell form of this matrix.

        The Howell form is an echelon form with a few additional properties
        that make it unique and useful for solving equations modulo `N`.

        For a matrix `A`, let `S(A)` be the module spanned by the rows, and
        `S_j(A)` the submodule consisting of vectors with the first `j` entries
        zero.  Recall that elements `a, b` of a ring `R` are associates if there
        is a unit `u \in R` with `a = ub`; this is an equivalence relation.
        We can choose a set `X` of representatives in `\ZZ/N\ZZ` by taking
        all products of the primes dividing `N`, raised to arbitrary powers.

        The Howell form of `A` is a matrix `H` with the same number of columns
        as `A` so that

        - it is an echelon form: the zero rows are at the bottom, and the leading
          entries (pivots) in each row each occur down and to the right of the previous.

        - Each pivot lies in `X`, and the entries above a pivot are smaller
          (as in Hermite normal form).

        - If `(i, j_i)` is the location of a pivot, then rows `i+1`, ... generate
          `S_{j_i}(A)`.

        if `r` is the number of nonzero rows of `H`, then the first `r` rows
          of `H` are nonzero

        .. NOTE::

            The Howell form will always have at least as many rows as columns.

        EXAMPLES::

            sage: entries = [624, 353, 352, 497, 182, 374, 556, 365, 271, 557, 203, 327]
            sage: A = matrix(Zmod(625), 3, 4, entries)
            sage: A.howell_form()
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            [ 0  0  0  0]

        The following example illustrates why the number of rows needs to be increased::

            sage: A = matrix(Zmod(625), [125, 25, 5, 1])
            sage: A.howell_form()
            [125  25   5   1]
            [  0 125  25   5]
            [  0   0 125  25]
            [  0   0   0 125]
        """
        key='howell_form'
        ans = self.fetch(key)
        if ans is not None:
            return ans

        if self._nrows < self._ncols:
            M = self._new(self._ncols - self._nrows, self._ncols)
            ans = self.stack(M)
        else:
            ans = self.__copy__()
        ans.howellize()
        self.cache(key, ans)
        return ans

    def __pow__(Matrix_modn_dense_flint self, n, dummy):
        """
        Exponentiation.

        EXAMPLES::

            sage: A = matrix(Zmod(625), 2, [1, 1, 0, 1])
            sage: A^625
            [1 0]
            [0 1]

            sage: N = 2^63 - 1
            sage: A = random_matrix(Zmod(N), 6)
            sage: B = random_matrix(Zmod(N), 6)
            sage: while not B.is_unit():
            ....:     B = random_matrix(Zmod(N), 6)
            sage: all((~B * A * B)^k == ~B * A^k * B for k in range(2, 10))
            True
        """
        #cdef Matrix_modn_dense_flint self = <Matrix_modn_dense_flint?>sself

        if dummy is not None:
            raise ValueError
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        cdef unsigned long e
        cdef long e_sgn
        cdef int err

        if integer_check_long_py(n, &e_sgn, &err) and not err:
            if e_sgn < 0:
                return (~self) ** (-n)
            e = <unsigned long>e_sgn
        else:
            if not isinstance(n, Integer):
                n = Integer(n)
            if mpz_sgn((<Integer>n).value) < 0:
                return (~self) ** (-n)

            if mpz_fits_ulong_p((<Integer>n).value):
                e = mpz_get_ui((<Integer>n).value)
            else:
                return generic_power(self, n)

        if e == 0:
            return self._parent.identity_matrix()
        if e == 1:
            return self

        cdef Matrix_modn_dense_flint M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_pow(M._matrix, self._matrix, e)
        sig_off()
        return M

    def _right_kernel_matrix(self, algorithm=None, proof=None, zero_divisors_are_pivots=False):
        """
        A matrix whose rows form a basis for the right kernel of this matrix.

        INPUT:

        - ``algorithm`` -- either ``flint`` or ``linbox``.

        - ``proof`` -- ignored, for compatibility with :meth:`sage.matrix.matrix2.Matrix.right_kernel_matrix`.

        - ``zero_divisors_are_pivots`` -- 

        OUTPUT:

        - ``format`` -- a string indicating what kind of echelonization has been done on the result

        - ``K`` -- a matrix so that right multiplication by ``K.transpose()`` yields zero.

        EXAMPLES::

        When the base ring is not a field, the kernel is computed using the Howell form::

            sage: A = matrix(Zmod(43^10), 6, [0, 5664461354126771, 12212357361910300, 15947020959157478, 0, 16792952041597449, 14690359073749623, 11237259451999884, 5117434014428142, 15157488677243483, 9004103062307752, 20761679499270441, 4620722392655416, 5445142895231681, 6605357538252496, 7608812697273777, 18542817615638637, 18194689690271501, 0, 20341333098836812, 12117922812876054, 1270149214447437, 0, 10999401748338075, 4620722392655416, 10891113386038365, 956055025271903, 2162842206467093, 18542817615638637, 1143972982339214, 13128267973348003, 15817056104759912, 20531311511260484, 13598045280630823, 7585782589268305, 14053895308766769])

            sage: A.howell_form()
            [               43                 0               430   463321664584385  6642784060940877  1272143397612677]
            [                0                43              1290    13718919206186  9222987387198236  4718167694579285]
            [                0                 0              1849   175909299267694  5349721549210673 15609543848590133]
            [                0                 0                 0   502592611936843 17088148805852662 20606297089410563]
            [                0                 0                 0                 0                 0                 0]
            [                0                 0                 0                 0                 0                 0]

            sage: K0 = A.right_kernel_matrix(basis='computed'); K0 # indirect doctest
            [502592611936843               0               0               0               0               0]
            [              0 502592611936843               0               0               0               0]
            [385710609160833 151946603608813  11688200277601               0               0               0]
            [465890695877871 260955102496367   7597286341143              43               0               0]
            [171748388290949  47073416965373   7938657150170               9               1               0]
            [420900330750694 300556574130784   3055770614472               2               0               1]
            sage: A * K0.T == 0
            True

        The default is to echelonize the kernel::

            sage: K = A.right_kernel_matrix(); K
            [              1               0               7  11391562385871 163874381071062   4349304245336]
            [              0               1              22    218579974005 169987310309906  32297642657385]
            [              0               0              43   3510762799844 121697346420828 148529822771001]
            [              0               0               0  11688200277601  35064600832803 350646008328030]
            [              0               0               0               0 502592611936843               0]
            [              0               0               0               0               0 502592611936843]
            sage: A * K.T == 0
            True
        """
        R = self._parent._base
        if algorithm == "default":
            if R.is_field() and min(self._nrows, self._ncols) > 200:
                from .matrix_modn_dense_double import MAX_MODULUS
                if R.order() < MAX_MODULUS:
                    algorithm = "linbox"
                else:
                    algorithm = "flint"
            else:
                algorithm = "flint"
        if algorithm == "linbox":
            K = self._change_implementation("linbox").right_kernel_matrix("linbox-noefd", basis="pivot")
            return "pivot-linboxed", K._change_implementation("flint")
        cdef Py_ssize_t i, j, k, l, cur_row, pivl
        cdef mp_limb_t s, x, y, N, xinv, yinv, Ninv
        cdef Matrix_modn_dense_flint X, ans, E = self.echelon_form()
        if R.is_field():
            # nmod_mut_nullspace will do this regardless
            # so we are better off to start in echelon form to have the rank
            X = self._new(self._ncols, self._ncols - self.rank())
            ans = self._new(self._ncols - self.rank(), self._ncols)
            sig_on()
            nmod_mat_nullspace(X._matrix, E._matrix) # columns of X form a basis
            nmod_mat_transpose(ans._matrix, X._matrix)
            sig_off()
            return "pivot-nmod-field", ans
        else:
            # We need to have a square matrix so that it's possible to echelonize.
            ans = self._new(self._ncols, self._ncols)
            N = mpz_get_si(self._modulus.sageInteger.value)
            Ninv = n_preinvert_limb(N)

            p, zdp = map(set, self._pivots())
            set_pivots = p.union(zdp)
            pivot = sorted(set_pivots)

            # In the comments below, we write v for the current row of ans being constructed
            cur_row = 0
            # k tracks the number of pivots we've passed
            k = 0
            for j in range(self._ncols):
                if j in p:
                    k += 1
                    continue
                i = k
                if j in zdp:
                    k += 1
                    if zero_divisors_are_pivots:
                        continue
                    #v[j] = N // E[i,j]
                    nmod_mat_set_entry(ans._matrix, cur_row, j, N // nmod_mat_get_entry(E._matrix, i, j))
                else:
                    #v[j] = 1
                    nmod_mat_set_entry(ans._matrix, cur_row, j, 1)

                # figure out the remaining coefficients of v
                # note that v might need to be rescaled several times
                for l in reversed(range(i)):
                    pivl = pivot[l]
                    x = nmod_mat_get_entry(E._matrix, l, pivl)
                    xinv = n_preinvert_limb(x)

                    # solve for T
                    # v[pivot[l]] * E[l, pivot[l]] + T*sum(v[m] * E[l,m] for m in range(pivot[l] + 1, j + 1)) = 0
                    # Start by computing the sum
                    #s = sum(v[m] * E[l, m] for m in range(pivot[l] + 1, j + 1))
                    s = 0
                    for m in range(pivl + 1, j + 1):
                        y = n_mulmod2_preinv(nmod_mat_get_entry(ans._matrix, cur_row, m),
                                             nmod_mat_get_entry(E._matrix, l, m), N, Ninv)
                        s = n_addmod(s, y, N)
                    # Make s divisible by x by working mod N/x
                    if n_mod2_preinv(s, x, xinv): # make sure we can work mod N/x
                    #if s % x != 0: # make sure we can work mod N/x
                        # set y = x // gcd(s, x), then multiply s by y
                        y = n_gcd(s, x)
                        yinv = n_preinvert_limb(y)
                        y = n_div2_preinv(x, y, yinv)
                        s = n_mulmod2_preinv(s, y, N, Ninv)
                        # now s is divisible by x, but we also have to multiply
                        # the other entries of the row by y
                        for m in range(pivl + 1, j + 1):
                            # v[m] *= y
                            nmod_mat_set_entry(
                                ans._matrix, cur_row, m,
                                n_mulmod2_preinv(
                                    nmod_mat_get_entry(ans._matrix, cur_row, m),
                                    y, N, Ninv))
                        # assert v[j] % N != 0
                    # This is correct modulo N/x, and since the pivot is x
                    # the kernel is now well defined modulo N.
                    #v[pivot[l]] = -s / x
                    nmod_mat_set_entry(
                        ans._matrix, cur_row, pivl,
                        n_div2_preinv(n_negmod(s, N), x, xinv))
                cur_row += 1
            return "pivot-nmod-ring", ans

    cdef int _copy_row_to_mod_int_array(self, mod_int *to, Py_ssize_t i):
        cdef Py_ssize_t j
        for j in range(self._ncols):
            to[j] = nmod_mat_get_entry(self._matrix, i, j)

    def _matrices_from_rows(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Make a list of matrices from the rows of this matrix.  This is a
        fairly technical function which is used internally, e.g., by
        the cyclotomic field linear algebra code.

        INPUT:

        - ``nrows`` - integer

        - ``ncols`` - integer

        OUTPUT:

        - ``list`` - a list of matrices, one for each row of the input, of size ``nrows`` by ``ncols``.

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

        """
        if nrows * ncols != self._ncols:
            raise ValueError("nrows * ncols must equal self's number of columns")

        cdef Matrix_modn_dense_flint M
        cdef Py_ssize_t i
        cdef Py_ssize_t n = self._ncols
        ans = []
        for i in range(self._nrows):
            M = self._new(nrows, ncols)
            if n:
                memcpy(nmod_mat_entry_ptr(M._matrix, 0, 0), nmod_mat_entry_ptr(self._matrix, i, 0), sizeof(mp_limb_t)*n)
            ans.append(M)
        return ans

    @staticmethod
    def _matrix_from_rows_of_matrices(X):
        """
        Return a matrix whose row ``i`` is constructed from the entries of
        matrix ``X[i]``.

        INPUT:

        - ``X`` - a nonempty list of matrices of the same size mod a
            single modulus `n`

        EXAMPLES::

            sage: X = [random_matrix(GF(17), 4, 4) for _ in range(10)]; X
            [
            [ 1  8  0  5]  [12  8  6 14]  [ 4  1  4  3]  [11  6  1 15]
            [ 8  3 11  1]  [ 7 15 10 11]  [13  1 10 16]  [ 3  9 14  6]
            [15  9 15 13]  [11  5  9 15]  [15 13 15  2]  [ 9 10  3 16]
            [ 2 14 14  8], [13 16 10 16], [13 13 14  1], [ 7  3  8  7],
            <BLANKLINE>
            [12 14  7 14]  [ 1 15 11 11]  [13  8  8  0]  [ 2 16 14 12]
            [ 6  1  7 11]  [ 2  2 11  0]  [ 5  9  6  6]  [ 7 11 15  7]
            [16 12  5  9]  [ 7 15  0  0]  [12 14  2 15]  [ 3  2 14  6]
            [ 0  2 14  0], [13  1 15  6], [ 1  2  5  3], [ 7  9 11  0],
            <BLANKLINE>
            [ 5  1  6  5]  [15  1 15 14]
            [15  1  9  2]  [ 6 14  9 10]
            [ 0 12  2 13]  [ 9  3  5  9]
            [ 1 13  8 14], [ 6  3  1 16]
            ]
            sage: X[0]._matrix_from_rows_of_matrices(X) # indirect doctest
            [ 1  8  0  5  8  3 11  1 15  9 15 13  2 14 14  8]
            [12  8  6 14  7 15 10 11 11  5  9 15 13 16 10 16]
            [ 4  1  4  3 13  1 10 16 15 13 15  2 13 13 14  1]
            [11  6  1 15  3  9 14  6  9 10  3 16  7  3  8  7]
            [12 14  7 14  6  1  7 11 16 12  5  9  0  2 14  0]
            [ 1 15 11 11  2  2 11  0  7 15  0  0 13  1 15  6]
            [13  8  8  0  5  9  6  6 12 14  2 15  1  2  5  3]
            [ 2 16 14 12  7 11 15  7  3  2 14  6  7  9 11  0]
            [ 5  1  6  5 15  1  9  2  0 12  2 13  1 13  8 14]
            [15  1 15 14  6 14  9 10  9  3  5  9  6  3  1 16]

        OUTPUT: A single matrix mod ``p`` whose ``i``-th row is ``X[i].list()``.
        """
        # The code below is just a fast version of the following:
        ##     from constructor import matrix
        ##     K = X[0].base_ring()
        ##     v = sum([y.list() for y in X],[])
        ##     return matrix(K, len(X), X[0].nrows()*X[0].ncols(), v)

        cdef Matrix_modn_dense_flint T = X[0]
        cdef Py_ssize_t i, j, copysize, n = len(X), m = T._nrows * T._ncols

        cdef Matrix_modn_dense_flint A = T._new(n, m)
        for i in range(n):
            T = X[i]
            copysize = T._ncols*sizeof(mp_limb_t)
            # rows could have been swapped around
            for j in range(T._nrows):
                memcpy(nmod_mat_entry_ptr(A._matrix, i, j*T._ncols), nmod_mat_entry_ptr(T._matrix, j, 0), copysize)
        return A
