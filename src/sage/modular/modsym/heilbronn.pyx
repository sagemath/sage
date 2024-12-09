# distutils: extra_compile_args = -D_XPG6
"""
Heilbronn matrix computation
"""
# ****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_allocarray, sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off, sig_check

from sage.misc.verbose import verbose
from sage.arith.misc import is_prime

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.flint.fmpz cimport fmpz_add
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpq_mat cimport fmpq_mat_entry_num

cdef extern from "<math.h>":
    float roundf(float x)

cimport sage.modular.modsym.p1list as p1list
from sage.modular.modsym import p1list
cdef p1list.export export
export = p1list.export()

from sage.modular.modsym.apply cimport Apply
cdef Apply PolyApply= Apply()

from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense
from sage.matrix.matrix_cyclo_dense cimport Matrix_cyclo_dense

ctypedef long long llong

cdef int llong_prod_mod(int a, int b, int N) noexcept:
    cdef int c
    c = <int> (((<llong> a) * (<llong> b)) % (<llong> N))
    if c < 0:
        c = c + N
    return c

cdef struct list:
    int *v
    int i   # how many positions of list are filled
    int n   # how much memory has been allocated

cdef int* expand(int *v, int n, int new_length) except NULL:
    cdef int *w
    cdef int i
    w = <int*>check_allocarray(new_length, sizeof(int))
    if v:
        for i in range(n):
            w[i] = v[i]
        sig_free(v)
    return w

cdef int list_append(list* L, int a) except -1:
    cdef int j
    if L.i >= L.n:
        j = 10 + 2*L.n
        L.v = expand(L.v, L.n, j)
        L.n = j
    L.v[L.i] = a
    L.i = L.i + 1

cdef int list_append4(list* L, int a, int b, int c, int d) except -1:
    list_append(L, a)
    list_append(L, b)
    list_append(L, c)
    list_append(L, d)

cdef void list_clear(list L) noexcept:
    sig_free(L.v)

cdef void list_init(list* L) noexcept:
    L.n = 16
    L.i = 0
    L.v = expand(<int*>0, 0, L.n)


cdef class Heilbronn:
    cdef int length
    cdef list list

    def __dealloc__(self):
        list_clear(self.list)

    def _initialize_list(self):
        """
        Initialize the list of matrices corresponding to ``self``.

        (This function is automatically called during initialization.)

        .. NOTE::

           This function must be overridden by all derived classes!

        EXAMPLES::

            sage: from sage.modular.modsym.heilbronn import Heilbronn
            sage: H = Heilbronn()
            sage: H._initialize_list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __getitem__(self, int n):
        r"""
        Return the `n`-th matrix in ``self``.

        EXAMPLES::

            sage: H = HeilbronnCremona(11)
            sage: H[17]
            [-1, 0, 1, -11]
            sage: H[98234]
            Traceback (most recent call last):
            ...
            IndexError
        """
        if n < 0 or n >= self.length:
            raise IndexError
        return [self.list.v[4*n], self.list.v[4*n+1],
                self.list.v[4*n+2], self.list.v[4*n+3]]

    def __len__(self):
        """
        Return the number of matrices in ``self``.

        EXAMPLES::

            sage: HeilbronnCremona(2).__len__()
            4
        """
        return self.length

    def to_list(self):
        """
        Return the list of Heilbronn matrices corresponding to ``self``.

        Each matrix is given as a list of four integers.

        EXAMPLES::

            sage: H = HeilbronnCremona(2); H
            The Cremona-Heilbronn matrices of determinant 2
            sage: H.to_list()
            [[1, 0, 0, 2], [2, 0, 0, 1], [2, 1, 0, 1], [1, 0, 1, 2]]
        """
        cdef int i
        L = []
        for i in range(self.length):
            L.append([self.list.v[4*i], self.list.v[4*i+1],
                      self.list.v[4*i+2], self.list.v[4*i+3]])
        return L

    cdef apply_only(self, int u, int v, int N, int* a, int* b):
        r"""
        INPUT:

        - ``u``, ``v``, ``N`` -- integers

        - ``a``, ``b`` -- preallocated integer arrays of the length of ``self``

        OUTPUT: sets the entries of `a`, `b`

        EXAMPLES::

            sage: M = ModularSymbols(33,2,1)       # indirect test
            sage: a, b = sage.modular.modsym.heilbronn.hecke_images_gamma0_weight2(1,0,33,[2,3],M.manin_gens_to_basis())
            sage: z = M((1,0))
            sage: [M.T(n)(z).element() for n in [2,2]] == [a, a]
            True
        """
        cdef Py_ssize_t i
        sig_on()
        if N == 1:   # easy special case
            for i in range(self.length):
                a[i] = b[i] = 0
        if N < 32768:   # use ints with no reduction modulo N
            for i in range(self.length):
                a[i] = u * self.list.v[4*i] + v * self.list.v[4*i+2]
                b[i] = u * self.list.v[4*i+1] + v * self.list.v[4*i+3]
        elif N < 46340:    # use ints but reduce mod N so can add two
            for i in range(self.length):
                a[i] = (u * self.list.v[4*i]) % N + (v * self.list.v[4*i+2]) % N
                b[i] = (u * self.list.v[4*i+1]) % N + (v * self.list.v[4*i+3]) % N
        else:
            for i in range(self.length):
                a[i] = llong_prod_mod(u, self.list.v[4*i], N) + llong_prod_mod(v, self.list.v[4*i+2], N)
                b[i] = llong_prod_mod(u, self.list.v[4*i+1], N) + llong_prod_mod(v, self.list.v[4*i+3], N)
        sig_off()

    cdef apply_to_polypart(self, fmpz_poly_t* ans, int i, int k):
        r"""
        INPUT:

        - ``ans`` -- ``fmpz_poly_t*``; pre-allocated an
          initialized array of ``self.length`` ``fmpz_poly_t``s
        - ``i`` -- integer
        - ``k`` -- integer

        OUTPUT: sets entries of ``ans``
        """
        cdef int j, m = k - 2
        for j in range(self.length):
            PolyApply.apply_to_monomial_flint(ans[j], i, m,
                                              self.list.v[4 * j],
                                              self.list.v[4 * j + 1],
                                              self.list.v[4 * j + 2],
                                              self.list.v[4 * j + 3])

    def apply(self, int u, int v, int N):
        r"""
        Return a list of pairs `((c,d),m)`, which is obtained as follows:

        1) Compute the images `(a,b)` of the vector `(u,v) \pmod N` acted on by
        each of the HeilbronnCremona matrices in ``self``.

        2) Reduce each `(a,b)` to canonical form `(c,d)` using ``p1normalize``.

        3) Sort.

        4) Create the list `((c,d),m)`, where `m` is the number of times
        that `(c,d)` appears in the list created in steps 1-3
        above. Note that the pairs `((c,d),m)` are sorted
        lexicographically by `(c,d)`.

        INPUT:

        - ``u``, ``v``, ``N`` -- integers

        OUTPUT: list

        EXAMPLES::

            sage: H = sage.modular.modsym.heilbronn.HeilbronnCremona(2); H
            The Cremona-Heilbronn matrices of determinant 2
            sage: H.apply(1,2,7)
            [((1, 1), 1), ((1, 4), 1), ((1, 5), 1), ((1, 6), 1)]
        """
        cdef int i, a, b, c, d, s
        cdef object X
        cdef dict M = {}
        t = verbose("start making list M.", level=5)
        if N < 32768:   # use ints with no reduction modulo N
            for i in range(self.length):
                sig_check()
                a = u*self.list.v[4*i] + v*self.list.v[4*i+2]
                b = u*self.list.v[4*i+1] + v*self.list.v[4*i+3]
                export.c_p1_normalize_int(N, a, b, &c, &d, &s, 0)
                X = (c, d)
                if X in M:
                    M[X] += 1
                else:
                    M[X] = 1
        elif N < 46340:    # use ints but reduce mod N so can add two
            for i in range(self.length):
                sig_check()
                a = (u * self.list.v[4*i]) % N + (v * self.list.v[4*i+2]) % N
                b = (u * self.list.v[4*i+1]) % N + (v * self.list.v[4*i+3]) % N
                export.c_p1_normalize_int(N, a, b, &c, &d, &s, 0)
                X = (c, d)
                if X in M:
                    M[X] += 1
                else:
                    M[X] = 1
        else:
            for i in range(self.length):
                sig_check()
                a = llong_prod_mod(u, self.list.v[4*i], N) + llong_prod_mod(v, self.list.v[4*i+2], N)
                b = llong_prod_mod(u, self.list.v[4*i+1], N) + llong_prod_mod(v, self.list.v[4*i+3], N)
                export.c_p1_normalize_llong(N, a, b, &c, &d, &s, 0)
                X = (c, d)
                if X in M:
                    M[X] += 1
                else:
                    M[X] = 1
        t = verbose("finished making list M.", t, level=5)
        mul = sorted(M.items())
        t = verbose("finished making mul list.", t, level=5)
        return mul


cdef class HeilbronnCremona(Heilbronn):
    cdef public int p

    def __init__(self, int p):
        r"""
        Create the list of Heilbronn-Cremona matrices of determinant `p`.

        EXAMPLES::

            sage: H = HeilbronnCremona(3) ; H
            The Cremona-Heilbronn matrices of determinant 3
            sage: H.to_list()
            [[1, 0, 0, 3],
            [3, 1, 0, 1],
            [1, 0, 1, 3],
            [3, 0, 0, 1],
            [3, -1, 0, 1],
            [-1, 0, 1, -3]]
        """
        if p <= 1 or not is_prime(p):
            raise ValueError("p must be >= 2 and prime")
        self.p = p
        self._initialize_list()

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: HeilbronnCremona(691).__repr__()
            'The Cremona-Heilbronn matrices of determinant 691'
        """
        return f"The Cremona-Heilbronn matrices of determinant {self.p}"

    def _initialize_list(self):
        """
        Initialize the list of matrices corresponding to ``self``.

        (This function is automatically called during initialization.)

        EXAMPLES::

            sage: from sage.modular.modsym.heilbronn import HeilbronnCremona
            sage: H = HeilbronnCremona.__new__(HeilbronnCremona)
            sage: H.p = 5
            sage: H
            The Cremona-Heilbronn matrices of determinant 5
            sage: H.to_list()
            []
            sage: H._initialize_list()
            sage: H.to_list()
            [[1, 0, 0, 5],
            [5, 2, 0, 1],
            [2, 1, 1, 3],
            [1, 0, 3, 5],
            [5, 1, 0, 1],
            [1, 0, 1, 5],
            [5, 0, 0, 1],
            [5, -1, 0, 1],
            [-1, 0, 1, -5],
            [5, -2, 0, 1],
            [-2, 1, 1, -3],
            [1, 0, -3, 5]]
        """
        cdef int r, x1, x2, y1, y2, a, b, c, x3, y3, q, p
        cdef list *L
        list_init(&self.list)
        L = &self.list
        p = self.p

        list_append4(L, 1, 0, 0, p)

        # When p==2, then Heilbronn matrices are
        #    [[1,0,0,2], [2,0,0,1], [2,1,0,1], [1,0,1,2]]
        # which are not given by the algorithm below.
        if p == 2:
            list_append4(L, 2, 0, 0, 1)
            list_append4(L, 2, 1, 0, 1)
            list_append4(L, 1, 0, 1, 2)
            self.length = 4
            return

        # NOTE: In C, -p/2 means "round toward 0", but in Python it
        # means "round down."
        sig_on()
        for r in range(-p // 2, p // 2 + 1):
            x1 = p
            x2 = -r
            y1 = 0
            y2 = 1
            a = -p
            b = r
            c = 0
            x3 = 0
            y3 = 0
            q = 0
            list_append4(L, x1, x2, y1, y2)
            while b:
                q = <int>roundf(<float>a / <float> b)
                c = a - b*q
                a = -b
                b = c
                x3 = q*x2 - x1
                x1 = x2
                x2 = x3
                y3 = q*y2 - y1
                y1 = y2
                y2 = y3
                list_append4(L, x1, x2, y1, y2)
        self.length = L.i // 4
        sig_off()


cdef class HeilbronnMerel(Heilbronn):
    cdef public int n

    def __init__(self, int n):
        r"""
        Initialize the list of Merel-Heilbronn matrices of determinant `n`.

        EXAMPLES::

            sage: H = HeilbronnMerel(3) ; H
            The Merel-Heilbronn matrices of determinant 3
            sage: H.to_list()
            [[1, 0, 0, 3],
            [1, 0, 1, 3],
            [1, 0, 2, 3],
            [2, 1, 1, 2],
            [3, 0, 0, 1],
            [3, 1, 0, 1],
            [3, 2, 0, 1]]
        """
        if n <= 0:
            raise ValueError("n (=%s) must be >= 1" % n)
        self.n = n
        self._initialize_list()

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: HeilbronnMerel(8).__repr__()
            'The Merel-Heilbronn matrices of determinant 8'
        """
        return "The Merel-Heilbronn matrices of determinant %s" % self.n

    def _initialize_list(self):
        """
        Initialize the list of matrices corresponding to ``self``.

        (This function is automatically called during initialization.)

        EXAMPLES::

            sage: from sage.modular.modsym.heilbronn import HeilbronnMerel
            sage: H = HeilbronnMerel.__new__(HeilbronnMerel)
            sage: H.n = 5
            sage: H
            The Merel-Heilbronn matrices of determinant 5
            sage: H.to_list()
            []
            sage: H._initialize_list()
            sage: H.to_list()
            [[1, 0, 0, 5],
            [1, 0, 1, 5],
            [1, 0, 2, 5],
            [1, 0, 3, 5],
            [1, 0, 4, 5],
            [2, 1, 1, 3],
            [2, 1, 3, 4],
            [3, 1, 1, 2],
            [3, 2, 2, 3],
            [4, 3, 1, 2],
            [5, 0, 0, 1],
            [5, 1, 0, 1],
            [5, 2, 0, 1],
            [5, 3, 0, 1],
            [5, 4, 0, 1]]
        """
        cdef int a, q, d, b, c, n
        cdef llong bc
        cdef list *L
        list_init(&self.list)
        L = &self.list
        n = self.n

        sig_on()
        for a in range(1, n+1):
            # We have ad-bc=n so c=0 and ad=n, or b=(ad-n)/c
            # Must have ad - n >= 0, so d must be >= Ceiling(n/a).
            q = n // a
            if q*a == n:
                d = q
                for b in range(a):
                    list_append4(L, a, b, 0, d)
                for c in range(1, d):
                    list_append4(L, a, 0, c, d)
            for d in range(q+1, n+1):
                bc = (<llong>a) * (<llong>d) - (<llong>n)
                # Divisor c of bc must satisfy Floor(bc/c) lt a and c lt d.
                # c ge (bc div a + 1)  <=>  Floor(bc/c) lt a  (for integers)
                # c le d - 1           <=>  c lt d
                for c in range(bc // a + 1, d):
                    if bc % c == 0:
                        list_append4(L, a, bc // c, c, d)
        self.length = L.i // 4
        sig_off()


############################################################################
# Fast application of Heilbronn operators to computation of
# systems of Hecke eigenvalues.
#   GAMMA0 trivial character weight 2 case
############################################################################


def hecke_images_gamma0_weight2(int u, int v, int N, indices, R):
    r"""
    INPUT:

    - ``u``, ``v``, ``N`` -- integers so that `\gcd(u,v,N) = 1`
    - ``indices`` -- list of positive integers
    - ``R`` -- matrix over `\QQ` that writes each elements of
      `\textnormal{P1} = \textnormal{P1List}(N)` in terms of a
      subset of `\textnormal{P1}`

    OUTPUT: a dense matrix whose columns are the images `T_n(x)`
    for `n` in indices and `x` the Manin symbol `(u,v)`, expressed
    in terms of the basis.

    EXAMPLES::

        sage: M = ModularSymbols(23,2,1)
        sage: A = sage.modular.modsym.heilbronn.hecke_images_gamma0_weight2(1,0,23,[1..6],M.manin_gens_to_basis())
        sage: rowsA = A.rows()
        sage: z = M((1,0))
        sage: all(M.T(n)(z).element() == rowsA[n-1] for n in [1..6])
        True

    TESTS::

        sage: M = ModularSymbols(389,2,1,GF(7))
        sage: C = M.cuspidal_subspace()
        sage: N = C.new_subspace()
        sage: D = N.decomposition()
        sage: D[1].q_eigenform(10, 'a') # indirect doctest
        q + 4*q^2 + 2*q^3 + 6*q^5 + q^6 + 5*q^7 + 6*q^8 + q^9 + O(q^10)
    """
    cdef p1list.P1List P1 = p1list.P1List(N)

    # Create a zero dense matrix over `\QQ` with len(indices) rows
    # and #P^1(N) columns.
    cdef Matrix_rational_dense T
    from sage.matrix.constructor import matrix
    from sage.rings.rational_field import QQ
    T = matrix(QQ, len(indices), len(P1), sparse=False)
    original_base_ring = R.base_ring()
    if original_base_ring != QQ:
        R = R.change_ring(QQ)

    cdef Py_ssize_t i, j
    cdef int *a
    cdef int *b
    cdef int k

    cdef Heilbronn H

    t = verbose("computing non-reduced images of symbol under Hecke operators",
                level=1, caller_name='hecke_images_gamma0_weight2')
    for i, n in enumerate(indices):
        # List the Heilbronn matrices of determinant n defined by Cremona or Merel
        H = HeilbronnCremona(n) if is_prime(n) else HeilbronnMerel(n)

        # Allocate memory to hold images of (u,v) under all Heilbronn matrices
        a = <int*> sig_malloc(sizeof(int)*H.length)
        if not a:
            raise MemoryError
        b = <int*> sig_malloc(sizeof(int)*H.length)
        if not b:
            raise MemoryError

        # Compute images of (u,v) under all Heilbronn matrices
        H.apply_only(u, v, N, a, b)

        # Compute the indexes of these images.
        # We just store them in the array a for simplicity.
        for j in range(H.length):
            # Compute index of the symbol a[j], b[j] in the standard list.
            k = P1.index(a[j], b[j])

            # Now fill in row i of the matrix T.
            if k != -1:
                # The following line is just a dangerous direct way to do: T[i,k] += 1
                T._add_ui_unsafe_assuming_int(i, k, 1)

        # Free a and b
        sig_free(a)
        sig_free(b)

    t = verbose("finished computing non-reduced images",
                t, level=1, caller_name='hecke_images_gamma0_weight2')

    t = verbose("Now reducing images of symbol",
                level=1, caller_name='hecke_images_gamma0_weight2')

    # Return the product T * R, whose rows are the image of (u,v) under
    # the Hecke operators T_n for n in indices.
    if max(indices) <= 30:   # In this case T tends to be very sparse
        ans = T.sparse_matrix()._matrix_times_matrix_dense(R)
        verbose("did reduction using sparse multiplication",
                t, level=1, caller_name='hecke_images_gamma0_weight2')
    elif R.is_sparse():
        ans = T * R.dense_matrix()
        verbose("did reduction using dense multiplication",
                t, level=1, caller_name='hecke_images_gamma0_weight2')
    else:
        ans = T * R
        verbose("did reduction using dense multiplication",
                t, level=1, caller_name='hecke_images_gamma0_weight2')

    if original_base_ring != QQ:
        ans = ans.change_ring(original_base_ring)

    return ans


############################################################################
# Fast application of Heilbronn operators to computation of
# systems of Hecke eigenvalues.
#   Nontrivial character but weight 2.
############################################################################


def hecke_images_nonquad_character_weight2(int u, int v, int N, indices, chi, R):
    r"""
    Return images of the Hecke operators `T_n` for `n`
    in the list indices, where `\chi` must be a quadratic Dirichlet
    character with values in `\QQ`.

    `R` is assumed to be the relation matrix of a weight modular symbols
    space over `\QQ` with character `\chi`.

    INPUT:

    - ``u``, ``v``, ``N`` -- integers so that `\gcd(u,v,N) = 1`
    - ``indices`` -- list of positive integers
    - ``chi`` -- a Dirichlet character that takes values
      in a nontrivial extension of `\QQ`
    - ``R`` -- matrix over `\QQ` that writes each elements of
      `\textnormal{P1} = \textnormal{P1List}(N)` in terms of a
      subset of `\textnormal{P1}`

    OUTPUT: a dense matrix with entries in the field `\QQ(\chi)` (the
    values of `\chi`) whose columns are the images `T_n(x)` for `n` in
    indices and `x` the Manin symbol `(u,v)`, expressed in terms of the
    basis.

    EXAMPLES::

        sage: chi = DirichletGroup(13).0^2
        sage: M = ModularSymbols(chi)
        sage: eps = M.character()
        sage: R = M.manin_gens_to_basis()
        sage: sage.modular.modsym.heilbronn.hecke_images_nonquad_character_weight2(1,0,13,[1,2,6],eps,R)
        [           1            0            0            0]
        [   zeta6 + 2            0            0           -1]
        [           7 -2*zeta6 + 1   -zeta6 - 1     -2*zeta6]
        sage: x = M((1,0)); x.element()
        (1, 0, 0, 0)
        sage: M.T(2)(x).element()
        (zeta6 + 2, 0, 0, -1)
        sage: M.T(6)(x).element()
        (7, -2*zeta6 + 1, -zeta6 - 1, -2*zeta6)
    """
    cdef p1list.P1List P1 = p1list.P1List(N)

    from sage.rings.rational_field import QQ
    K = chi.base_ring()

    if K == QQ:
        raise TypeError("character must not be trivial or quadratic")

    if R.base_ring() != K:
        R = R.change_ring(K)

    # Create a zero dense matrix over K with len(indices) rows
    # and #P^1(N) columns.
    cdef Matrix_cyclo_dense T
    from sage.matrix.constructor import matrix
    T = matrix(K, len(indices), len(P1), sparse=False)

    cdef Py_ssize_t i, j
    cdef int *a
    cdef int *b
    cdef int k, scalar

    cdef Heilbronn H

    t = verbose("computing non-reduced images of symbol under Hecke operators",
                level=1, caller_name='hecke_images_character_weight2')

    # Make a matrix over the rational numbers each of whose columns
    # are the values of the character chi.
    cdef Matrix_rational_dense chi_vals
    z = [t.list() for t in chi.values()]
    chi_vals = matrix(QQ, z).transpose()

    for i, n in enumerate(indices):
        H = HeilbronnCremona(n) if is_prime(n) else HeilbronnMerel(n)

        # Allocate memory to hold images of (u,v) under all Heilbronn matrices
        a = <int*> sig_malloc(sizeof(int)*H.length)
        if not a:
            raise MemoryError
        b = <int*> sig_malloc(sizeof(int)*H.length)
        if not b:
            raise MemoryError

        # Compute images of (u,v) under all Heilbronn matrices
        H.apply_only(u, v, N, a, b)

        for j in range(H.length):
            # Compute index of the symbol a[j], b[j] in the standard list.
            P1.index_and_scalar(a[j], b[j], &k, &scalar)
            # Now fill in row i of the matrix T.
            if k != -1:
                # The following line is just a dangerous direct way to do: T[i,k] += chi(scalar)
                # T[i,k] += chi(scalar)
                # This code makes assumptions about the internal structure
                # of matrices over cyclotomic fields.  It's nasty, but it
                # is exactly what is needed to get a solid 100 or more
                # times speedup.
                scalar %= N
                if scalar < 0:
                    scalar += N
                # Note that the next line totally dominates the runtime of this whole function.
                T._matrix._add_col_j_of_A_to_col_i_of_self(i * T._ncols + k, chi_vals, scalar)

        # Free a and b
        sig_free(a)
        sig_free(b)

    return T * R


def hecke_images_quad_character_weight2(int u, int v, int N, indices, chi, R):
    r"""
    INPUT:

    - ``u``, ``v``, ``N`` -- integers so that `\gcd(u,v,N) = 1`
    - ``indices`` -- list of positive integers
    - ``chi`` -- a Dirichlet character that takes values in `\QQ`
    - ``R`` -- matrix over `\QQ(\chi)` that writes each elements of
      `\textnormal{P1} = \textnormal{P1List}(N)` in terms of a subset
      of `\textnormal{P1}`

    OUTPUT: a dense matrix with entries in the rational field `\QQ` (the
    values of `\chi`) whose columns are the images `T_n(x)` for `n` in
    indices and `x` the Manin symbol `(u,v)`, expressed in terms of the
    basis.

    EXAMPLES::

        sage: chi = DirichletGroup(29,QQ).0
        sage: M = ModularSymbols(chi)
        sage: R = M.manin_gens_to_basis()
        sage: sage.modular.modsym.heilbronn.hecke_images_quad_character_weight2(2,1,29,[1,3,4],chi,R)
        [ 0  0  0  0  0 -1]
        [ 0  1  0  1  1  1]
        [ 0 -2  0  2 -2 -1]
        sage: x = M((2,1)) ; x.element()
        (0, 0, 0, 0, 0, -1)
        sage: M.T(3)(x).element()
        (0, 1, 0, 1, 1, 1)
        sage: M.T(4)(x).element()
        (0, -2, 0, 2, -2, -1)
    """
    cdef p1list.P1List P1 = p1list.P1List(N)
    from sage.rings.rational_field import QQ
    if chi.base_ring() != QQ:
        raise TypeError("character must takes values in QQ")

    # Create a zero dense matrix over `\QQ` with len(indices) rows
    # and #P^1(N) columns.
    cdef Matrix_rational_dense T
    from sage.matrix.constructor import matrix
    T = matrix(QQ, len(indices), len(P1), sparse=False)

    if R.base_ring() != QQ:
        R = R.change_ring(QQ)

    cdef Py_ssize_t i, j
    cdef int *a
    cdef int *b
    cdef int k, scalar
    cdef Heilbronn H

    _ = verbose("computing non-reduced images of symbol under Hecke operators",
                level=1, caller_name='hecke_images_quad_character_weight2')

    # Make a matrix over the rational numbers each of whose columns
    # are the values of the character chi.
    _chivals = chi.values()
    cdef int *chi_vals = <int*>sig_malloc(sizeof(int)*len(_chivals))
    if not chi_vals:
        raise MemoryError
    for i in range(len(_chivals)):
        chi_vals[i] = _chivals[i]

    for i, n in enumerate(indices):
        H = HeilbronnCremona(n) if is_prime(n) else HeilbronnMerel(n)
        a = <int*> sig_malloc(sizeof(int)*H.length)
        if not a:
            raise MemoryError
        b = <int*> sig_malloc(sizeof(int)*H.length)
        if not b:
            raise MemoryError

        H.apply_only(u, v, N, a, b)
        for j in range(H.length):
            P1.index_and_scalar(a[j], b[j], &k, &scalar)
            if k != -1:
                # This is just T[i,k] += chi(scalar)
                scalar %= N
                if scalar < 0:
                    scalar += N
                if chi_vals[scalar] > 0:
                    T._add_ui_unsafe_assuming_int(i, k, 1)
                elif chi_vals[scalar] < 0:
                    T._sub_ui_unsafe_assuming_int(i, k, 1)
        sig_free(a)
        sig_free(b)

    sig_free(chi_vals)
    return T * R


# ##########################################################################
# Fast application of Heilbronn operators to computation of
# systems of Hecke eigenvalues.
#   Trivial character and weight > 2.
# ##########################################################################

def hecke_images_gamma0_weight_k(int u, int v, int i, int N, int k, indices, R):
    r"""
    INPUT:

    - ``u``, ``v``, ``N`` -- integers so that `\gcd(u,v,N) = 1`
    - ``i`` -- integer with `0 \le i \le k-2`
    - ``k`` -- weight
    - ``indices`` -- list of positive integers
    - ``R`` -- matrix over `\QQ` that writes each elements of
      `\textnormal{P1} = \textnormal{P1List}(N)` in terms of a
      subset of `\textnormal{P1}`

    OUTPUT: a dense matrix with rational entries whose columns are the
    images `T_n(x)` for `n` in indices and `x` the Manin symbol
    [`X^i*Y^{(k-2-i)}`, `(u,v)`], expressed in terms of the basis.

    EXAMPLES::

        sage: M = ModularSymbols(15,6,sign=-1)
        sage: R = M.manin_gens_to_basis()
        sage: a,b,c = sage.modular.modsym.heilbronn.hecke_images_gamma0_weight_k(4,1,3,15,6,[1,11,12], R)
        sage: x = M((3,4,1)) ; x.element() == a
        True
        sage: M.T(11)(x).element() == b
        True
        sage: M.T(12)(x).element() == c
        True
    """
    cdef p1list.P1List P1 = p1list.P1List(N)

    # The Manin symbols are enumerated as
    #   all [0, (u,v)] for (u,v) in P^1(N) then
    #   all [1, (u,v)] for (u,v) in P^1(N) etc.
    # So we create a zero dense matrix over `\QQ` with len(indices) rows
    # and #P^1(N) * (k-1) columns.
    cdef Matrix_rational_dense T
    from sage.matrix.constructor import matrix
    from sage.rings.rational_field import QQ
    T = matrix(QQ, len(indices), len(P1)*(k-1), sparse=False)

    if R.base_ring() != QQ:
        R = R.change_ring(QQ)

    cdef Py_ssize_t j, m, z, w, n, p
    cdef int *a
    cdef int *b

    n = len(P1)

    cdef Heilbronn H
    cdef fmpz_poly_t* poly

    for z, m in enumerate(indices):
        H = HeilbronnCremona(m) if is_prime(m) else HeilbronnMerel(m)

        # Allocate memory to hold images of (u,v) under all Heilbronn matrices
        a = <int*> sig_malloc(sizeof(int)*H.length)
        if not a:
            raise MemoryError
        b = <int*> sig_malloc(sizeof(int)*H.length)
        if not b:
            raise MemoryError

        # Compute images of (u,v) under all Heilbronn matrices
        H.apply_only(u, v, N, a, b)

        # Compute images of X^i Y^(2-k-i) under each Heilbronn matrix
        poly = <fmpz_poly_t*> sig_malloc(sizeof(fmpz_poly_t)*H.length)
        for j in range(H.length):
            fmpz_poly_init(poly[j])

        # The following line dominates the runtime of this entire function:
        H.apply_to_polypart(poly, i, k)

        # Compute the indexes of these images.
        # We just store them in the array a for simplicity.
        for j in range(H.length):
            # Compute index of the symbol a[j], b[j] in the standard list.
            p = P1.index(a[j], b[j])
            # Now fill in row z of the matrix T.
            if p != -1:
                for w in range(min(fmpz_poly_length(poly[j]), k-1)):
                    # The following two lines are just a vastly faster version of:
                    #           T[z, n*w + p] += poly[j][w]
                    # They use that we know that T has only integer entries.
                    fmpz_add(fmpq_mat_entry_num(T._matrix, z, n*w+p),
                             fmpq_mat_entry_num(T._matrix, z, n*w+p),
                             fmpz_poly_get_coeff_ptr(poly[j], w))

        # Free a and b
        sig_free(a)
        sig_free(b)

        # Free poly part
        for j in range(H.length):
            fmpz_poly_clear(poly[j])
        sig_free(poly)

    # Return the product T * R, whose rows are the image of (u,v) under
    # the Hecke operators T_n for n in indices.
    return T * R.dense_matrix()


############################################################################
# Fast application of Heilbronn operators to computation of
# systems of Hecke eigenvalues.
#   Nontrivial character of order > 2 and weight > 2
############################################################################

# TODO

############################################################################
# Fast application of Heilbronn operators to computation of
# systems of Hecke eigenvalues.
#   Nontrivial character of order = 2 and weight > 2
############################################################################

# TODO
