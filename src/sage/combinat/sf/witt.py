"""
Witt symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#                     2013 Darij Grinberg <darijgrinberg@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from . import multiplicative
from sage.matrix.constructor import matrix

class SymmetricFunctionAlgebra_witt(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    r"""
    The Witt symmetric function basis (or Witt basis, to be short).

    The Witt basis of the ring of symmetric functions is
    denoted by `(x_{\lambda})` in [HazWitt1]_, section 9.63, and by
    `(q_{\lambda})` in [DoranIV1996]_. We will denote this basis by
    `(w_{\lambda})` (which is precisely how it is denoted in
    [GriRei18]_, Exercise 2.9.3(d)). It is a multiplicative basis
    (meaning that `w_{\emptyset} = 1` and that every partition
    `\lambda` satisfies
    `w_{\lambda} = w_{\lambda_1} w_{\lambda_2} w_{\lambda_3} \cdots`,
    where `w_i` means `w_{(i)}` for every nonnegative integer `i`).

    This basis can be defined in various ways. Probably the most
    well-known one is using the equation

    .. MATH::

        \prod_{d=1}^{\infty} (1 - w_d t^d)^{-1} = \sum_{n=0}^{\infty} h_n t^n

    where `t` is a formal variable and `h_n` are the complete
    homogeneous symmetric functions, extended to `0` by `h_0 = 1`.
    This equation allows one to uniquely determine the functions
    `w_1, w_2, w_3, \ldots` by recursion; one consequently extends the
    definition to all `w_{\lambda}` by requiring multiplicativity.

    A way to rewrite the above equation without power series is:

    .. MATH::

        h_n = \sum_{\lambda \vdash n} w_{\lambda}

    for all nonnegative integers `n`, where `\lambda \vdash n` means
    that `\lambda` is a partition of `n`.

    A similar equation (which is easily seen to be equivalent to the
    former) is

    .. MATH::

        e_n = \sum_{\lambda} (-1)^{n - \ell(\lambda)} w_{\lambda},

    with the sum running only over *strict* partitions `\lambda` of
    `n` this time. This equation can also be used to recursively
    define the `w_n`. Furthermore, every positive integer `n`
    satisfies

    .. MATH::

        p_n = \sum_{d\mid n} d w_d^{n/d},

    and this can be used to define the `w_n` recursively over any
    ring which is torsion-free as a `\ZZ`-module. While these
    equations all yield easy formulas for classical bases of the
    ring of symmetric functions in terms of the Witt symmetric
    functions, it seems difficult to obtain explicit formulas in
    the other direction.

    The Witt symmetric functions owe their name to the fact that
    the ring of symmetric functions can be viewed as the coordinate
    ring of the group scheme of Witt vectors, and the Witt
    symmetric functions are the functions that send a Witt vector
    to its components (whereas the powersum symmetric functions
    send a Witt vector to its ghost components). Details can be
    found in [HazWitt1]_ or section 3.2 of [BorWi2004]_.

    INPUT:

    - ``Sym`` -- an instance of the ring of the symmetric functions.
    - ``coerce_h`` -- (default: ``True``) a boolean that determines
      whether the transition maps between the Witt basis and the
      complete homogeneous basis will be cached and registered as
      coercions.
    - ``coerce_e`` -- (default: ``False``) a boolean that determines
      whether the transition maps between the Witt basis and the
      elementary symmetric basis will be cached and registered as
      coercions.
    - ``coerce_p`` -- (default: ``False``) a boolean that determines
      whether the transition maps between the Witt basis and the
      powersum basis will be cached and registered as coercions (or
      conversions, if the base ring is not a `\QQ`-algebra).

    REFERENCES:

    .. [HazWitt1] Michiel Hazewinkel. *Witt vectors. Part 1*.
       :arxiv:`0804.3888v1`

    .. [DoranIV1996] William F. Doran IV.
       *A Proof of Reutenauer's `-q_{(n)}` Conjecture*.
       Journal of combinatorial theory, Series A 74, pp. 342-344 (1996),
       article no. 0056. :doi:`10.1006/jcta.1996.0056`

    .. [BorWi2004] James Borger, Ben Wieland.
       *Plethystic algebra*.
       :arxiv:`math/0407227v1`

    .. [GriRei18]_

    EXAMPLES:

    Here are the first few Witt symmetric functions, in various bases::

        sage: Sym = SymmetricFunctions(QQ)
        sage: w = Sym.w()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: p = Sym.p()
        sage: s = Sym.s()
        sage: m = Sym.m()

        sage: p(w([1]))
        p[1]
        sage: m(w([1]))
        m[1]
        sage: e(w([1]))
        e[1]
        sage: h(w([1]))
        h[1]
        sage: s(w([1]))
        s[1]

        sage: p(w([2]))
        -1/2*p[1, 1] + 1/2*p[2]
        sage: m(w([2]))
        -m[1, 1]
        sage: e(w([2]))
        -e[2]
        sage: h(w([2]))
        -h[1, 1] + h[2]
        sage: s(w([2]))
        -s[1, 1]

        sage: p(w([3]))
        -1/3*p[1, 1, 1] + 1/3*p[3]
        sage: m(w([3]))
        -2*m[1, 1, 1] - m[2, 1]
        sage: e(w([3]))
        -e[2, 1] + e[3]
        sage: h(w([3]))
        -h[2, 1] + h[3]
        sage: s(w([3]))
        -s[2, 1]

        sage: Sym = SymmetricFunctions(ZZ)
        sage: w = Sym.w()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: s = Sym.s()
        sage: m = Sym.m()
        sage: p = Sym.p()
        sage: m(w([4]))
        -9*m[1, 1, 1, 1] - 4*m[2, 1, 1] - 2*m[2, 2] - m[3, 1]
        sage: e(w([4]))
        -e[2, 1, 1] + e[3, 1] - e[4]
        sage: h(w([4]))
        -h[1, 1, 1, 1] + 2*h[2, 1, 1] - h[2, 2] - h[3, 1] + h[4]
        sage: s(w([4]))
        -s[1, 1, 1, 1] - s[2, 1, 1] - s[2, 2] - s[3, 1]

    Some examples of conversions the other way::

        sage: w(h[3])
        w[1, 1, 1] + w[2, 1] + w[3]
        sage: w(e[3])
        -w[2, 1] + w[3]
        sage: w(m[2,1])
        2*w[2, 1] - 3*w[3]
        sage: w(p[3])
        w[1, 1, 1] + 3*w[3]

    Antipodes::

        sage: w([1]).antipode()
        -w[1]
        sage: w([2]).antipode()
        -w[1, 1] - w[2]

    The following holds for all odd `i` and is easily proven by
    induction::

        sage: all( w([i]).antipode() == -w([i]) for i in range(1, 10, 2) )
        True

    The Witt basis does not allow for simple expressions for
    comultiplication and antipode in general (this is related to the
    fact that the sum of two Witt vectors isn't easily described in
    terms of the components). Therefore, most computations with Witt
    symmetric functions, as well as conversions and coercions, pass
    through the complete homogeneous symmetric functions by default.
    However, one can also use the elementary symmetric functions
    instead, or (if the base ring is a `\QQ`-algebra) the powersum
    symmetric functions. This is what the optional keyword variables
    ``coerce_e``, ``coerce_h`` and ``coerce_p`` are for. These
    variables do not affect the results of the (non-underscored)
    methods of ``self``, but they affect the speed of the computations
    (the more of these variables are set to ``True``, the
    faster these are) and the size of the cache (the more of
    these variables are set to ``True``, the bigger the cache). Let us
    check that the results are the same no matter to what the
    variables are set::

        sage: Sym = SymmetricFunctions(QQ)
        sage: p = Sym.p()
        sage: wh = Sym.w()
        sage: we = Sym.w(coerce_h=False, coerce_e=True)
        sage: wp = Sym.w(coerce_h=False, coerce_p=True)
        sage: all( p(wh(lam)) == p(we(lam)) == p(wp(lam)) for lam in Partitions(4) )
        True
        sage: all ( wh(p(lam)).monomial_coefficients()
        ....:       == we(p(lam)).monomial_coefficients()
        ....:       == wp(p(lam)).monomial_coefficients() for lam in Partitions(4) )
        True

    TESTS:

    Let us check that all the above computations work with a
    non-default setting as well::

        sage: Sym = SymmetricFunctions(QQ)
        sage: w = Sym.w(coerce_h=False, coerce_p=True)
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: p = Sym.p()
        sage: s = Sym.s()
        sage: m = Sym.m()

        sage: p(w([1]))
        p[1]
        sage: m(w([1]))
        m[1]
        sage: e(w([1]))
        e[1]
        sage: h(w([1]))
        h[1]
        sage: s(w([1]))
        s[1]

        sage: p(w([2]))
        -1/2*p[1, 1] + 1/2*p[2]
        sage: m(w([2]))
        -m[1, 1]
        sage: e(w([2]))
        -e[2]
        sage: h(w([2]))
        -h[1, 1] + h[2]
        sage: s(w([2]))
        -s[1, 1]

        sage: p(w([3]))
        -1/3*p[1, 1, 1] + 1/3*p[3]
        sage: m(w([3]))
        -2*m[1, 1, 1] - m[2, 1]
        sage: e(w([3]))
        -e[2, 1] + e[3]
        sage: h(w([3]))
        -h[2, 1] + h[3]
        sage: s(w([3]))
        -s[2, 1]

        sage: Sym = SymmetricFunctions(ZZ)
        sage: w = Sym.w()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: s = Sym.s()
        sage: m = Sym.m()
        sage: p = Sym.p()
        sage: m(w([4]))
        -9*m[1, 1, 1, 1] - 4*m[2, 1, 1] - 2*m[2, 2] - m[3, 1]
        sage: e(w([4]))
        -e[2, 1, 1] + e[3, 1] - e[4]
        sage: h(w([4]))
        -h[1, 1, 1, 1] + 2*h[2, 1, 1] - h[2, 2] - h[3, 1] + h[4]
        sage: s(w([4]))
        -s[1, 1, 1, 1] - s[2, 1, 1] - s[2, 2] - s[3, 1]

        sage: w(h[3])
        w[1, 1, 1] + w[2, 1] + w[3]
        sage: w(e[3])
        -w[2, 1] + w[3]
        sage: w(m[2,1])
        2*w[2, 1] - 3*w[3]
        sage: w(p[3])
        w[1, 1, 1] + 3*w[3]

        sage: w([1]).antipode()
        -w[1]
        sage: w([2]).antipode()
        -w[1, 1] - w[2]
        sage: all( w([i]).antipode() == -w([i]) for i in range(1, 10, 2) )
        True

    Another non-default setting::

        sage: Sym = SymmetricFunctions(QQ)
        sage: w = Sym.w(coerce_h=False, coerce_e=True)
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: p = Sym.p()
        sage: s = Sym.s()
        sage: m = Sym.m()

        sage: p(w([1]))
        p[1]
        sage: m(w([1]))
        m[1]
        sage: e(w([1]))
        e[1]
        sage: h(w([1]))
        h[1]
        sage: s(w([1]))
        s[1]

        sage: p(w([2]))
        -1/2*p[1, 1] + 1/2*p[2]
        sage: m(w([2]))
        -m[1, 1]
        sage: e(w([2]))
        -e[2]
        sage: h(w([2]))
        -h[1, 1] + h[2]
        sage: s(w([2]))
        -s[1, 1]

        sage: p(w([3]))
        -1/3*p[1, 1, 1] + 1/3*p[3]
        sage: m(w([3]))
        -2*m[1, 1, 1] - m[2, 1]
        sage: e(w([3]))
        -e[2, 1] + e[3]
        sage: h(w([3]))
        -h[2, 1] + h[3]
        sage: s(w([3]))
        -s[2, 1]

        sage: Sym = SymmetricFunctions(ZZ)
        sage: w = Sym.w()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: s = Sym.s()
        sage: m = Sym.m()
        sage: p = Sym.p()
        sage: m(w([4]))
        -9*m[1, 1, 1, 1] - 4*m[2, 1, 1] - 2*m[2, 2] - m[3, 1]
        sage: e(w([4]))
        -e[2, 1, 1] + e[3, 1] - e[4]
        sage: h(w([4]))
        -h[1, 1, 1, 1] + 2*h[2, 1, 1] - h[2, 2] - h[3, 1] + h[4]
        sage: s(w([4]))
        -s[1, 1, 1, 1] - s[2, 1, 1] - s[2, 2] - s[3, 1]
        sage: [type(coeff) for a, coeff in h(w([4]))]
        [<class 'sage.rings.integer.Integer'>,
         <class 'sage.rings.integer.Integer'>,
         <class 'sage.rings.integer.Integer'>,
         <class 'sage.rings.integer.Integer'>,
         <class 'sage.rings.integer.Integer'>]

        sage: w(h[3])
        w[1, 1, 1] + w[2, 1] + w[3]
        sage: w(e[3])
        -w[2, 1] + w[3]
        sage: w(m[2,1])
        2*w[2, 1] - 3*w[3]
        sage: w(p[3])
        w[1, 1, 1] + 3*w[3]

        sage: w([1]).antipode()
        -w[1]
        sage: w([2]).antipode()
        -w[1, 1] - w[2]
        sage: all( w([i]).antipode() == -w([i]) for i in range(1, 10, 2) )
        ....:      #this holds for all odd i and is easily proven by induction
        True
    """
    def __init__(self, Sym, coerce_h=True, coerce_e=False, coerce_p=False):
        """
        Initialize ``self``.

        TESTS::

            sage: w = SymmetricFunctions(QQ).w()
            sage: TestSuite(w).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(w).run(elements = [w[1,1]+w[2], w[1]+2*w[1,1]])
        """
        self._coerce_h = coerce_h
        self._coerce_e = coerce_e
        self._coerce_p = coerce_p
        multiplicative.SymmetricFunctionAlgebra_multiplicative.__init__(self, Sym, "Witt", 'w')

    def _precompute_cache(self, n, to_self_cache, from_self_cache, transition_matrices, inverse_transition_matrices, to_self_gen_function):
        """
        Compute the transition matrices between ``self`` and another
        multiplicative homogeneous basis in the homogeneous components of
        degree `n`.

        The results are not returned, but rather stored in the caches.

        This assumes that the transition matrices in all degrees smaller
        than `n` have already been computed and cached!

        INPUT:

        - ``n`` -- nonnegative integer
        - ``to_self_cache`` -- a cache which stores the coordinates of
          the elements of the other basis with respect to the
          basis ``self``
        - ``from_self_cache`` -- a cache which stores the coordinates
          of the elements of ``self`` with respect to the other
          basis
        - ``transition_matrices`` -- a cache for transition matrices
          which contain the coordinates of the elements of the other
          basis with respect to ``self``
        - ``inverse_transition_matrices`` -- a cache for transition
          matrices which contain the coordinates of the elements of
          ``self`` with respect to the other basis
        - ``to_self_gen_function`` -- a function which takes a
          positive integer `n` and returns the element of the other
          basis corresponding to the partition `[n]` expanded with
          respect to the Witt basis ``self`` (as an element of
          ``self``, not as a dictionary)

        Examples for usage of this function are the ``_precompute_h``,
        ``_precompute_e`` and ``_precompute_p`` methods of this class.

        EXAMPLES:

        The examples below demonstrate how the caches are built
        step by step using the ``_precompute_cache`` method. In order
        not to influence the outcome of other doctests, we make sure
        not to use the caches internally used by this class, but
        rather to create new caches::

            sage: Sym = SymmetricFunctions(QQ)
            sage: w = Sym.w()
            sage: toy_to_self_cache = {}
            sage: toy_from_self_cache = {}
            sage: toy_transition_matrices = {}
            sage: toy_inverse_transition_matrices = {}
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(toy_to_self_cache)
            []
            sage: def toy_gen_function(n):
            ....:     if n > 1:
            ....:         return w(Partition([n])) + n * w(Partition([n-1,1]))
            ....:     return w(Partition([n]))
            sage: w._precompute_cache(0, toy_to_self_cache,
            ....:                        toy_from_self_cache,
            ....:                        toy_transition_matrices,
            ....:                        toy_inverse_transition_matrices,
            ....:                        toy_gen_function)
            sage: l(toy_to_self_cache)
            [([], [([], 1)])]
            sage: w._precompute_cache(1, toy_to_self_cache,
            ....:                        toy_from_self_cache,
            ....:                        toy_transition_matrices,
            ....:                        toy_inverse_transition_matrices,
            ....:                        toy_gen_function)
            sage: l(toy_to_self_cache)
            [([], [([], 1)]), ([1], [([1], 1)])]
            sage: w._precompute_cache(2, toy_to_self_cache,
            ....:                        toy_from_self_cache,
            ....:                        toy_transition_matrices,
            ....:                        toy_inverse_transition_matrices,
            ....:                        toy_gen_function)
            sage: l(toy_to_self_cache)
            [([], [([], 1)]),
             ([1], [([1], 1)]),
             ([1, 1], [([1, 1], 1)]),
             ([2], [([1, 1], 2), ([2], 1)])]
            sage: toy_transition_matrices[2]
            [1 2]
            [0 1]
            sage: toy_inverse_transition_matrices[2]
            [ 1 -2]
            [ 0  1]
            sage: sorted(toy_transition_matrices)
            [0, 1, 2]
        """
        # Much of this code is adapted from dual.py
        base_ring = self.base_ring()
        zero = base_ring.zero()

        from sage.combinat.partition import Partition, Partitions_n

        # Handle the n == 0 case separately
        if n == 0:
            part = Partition([])
            one = base_ring.one()
            to_self_cache[ part ] = { part: one }
            from_self_cache[ part ] = { part: one }
            transition_matrices[n] = matrix(base_ring, [[one]])
            inverse_transition_matrices[n] = matrix(base_ring, [[one]])
            return

        partitions_n = Partitions_n(n).list()

        # The other basis will be called B from now on.

        # This contains the data for the transition matrix from the
        # basis B to the Witt basis self.
        transition_matrix_n = matrix(base_ring, len(partitions_n), len(partitions_n))

        # This first section calculates how the basis elements of the
        # basis B are expressed in terms of the Witt basis ``self``.

        # For every partition p of size n, expand B[p] in terms of
        # the Witt basis self using multiplicativity and
        # to_self_gen_function.
        i = 0
        for s_part in partitions_n:
            # s_mcs will be self(B[s_part])._monomial_coefficients
            s_mcs = {}

            # We need to compute the coordinates of B[s_part] in the Witt basis.
            hsp_in_w_basis = self.one()
            for p in s_part:
                hsp_in_w_basis *= to_self_gen_function(p)
            # Now, hsp_in_w_basis is B[s_part] expanded in the Witt
            # basis self (this is the same as the coercion self(B[s_part]).
            j = 0
            for p_part in partitions_n:

                if p_part in hsp_in_w_basis._monomial_coefficients:
                    sp = hsp_in_w_basis._monomial_coefficients[p_part]
                    s_mcs[p_part] = sp
                    transition_matrix_n[i,j] = sp

                j += 1

            to_self_cache[ s_part ] = s_mcs
            i += 1

        # Save the transition matrix
        transition_matrices[n] = transition_matrix_n

        # This second section calculates how the basis elements of
        # self expand in terms of the basis B.  We do this by
        # computing the inverse of the matrix transition_matrix_n
        # obtained above.
        # TODO: Possibly this can be sped up by using properties
        # of this matrix (e. g., it being triangular in most standard cases).
        # Are there significantly faster ways to invert a triangular
        # matrix (compared to the usual matrix inversion algorithms)?
        inverse_transition = (~transition_matrix_n).change_ring(base_ring)
        # Note that we don't simply write
        # "inverse_transition = ~transition_matrix_n" because that
        # tends to cast the entries of the matrix into a quotient
        # field even if this is unnecessary.

        # TODO: This still looks fragile when the base ring is weird!
        # Possibly work over ZZ in this method?

        for i in range(len(partitions_n)):
            d_mcs = {}
            for j in range(len(partitions_n)):
                if inverse_transition[i,j] != zero:
                    d_mcs[ partitions_n[j] ] = inverse_transition[i,j]

            from_self_cache[ partitions_n[i] ] = d_mcs

        inverse_transition_matrices[n] = inverse_transition

    def _precompute_h(self, n):
        """
        Compute the transition matrices between ``self`` and the complete
        homogeneous basis in the homogeneous components of degree `n`
        (and in those of smaller degree, if not already computed).
        The result is not returned, but rather stored in the cache.

        This assumes that the ``coerce_h`` keyword has been set to
        ``True`` in the initialization of ``self`` (otherwise the cache
        does not exist).

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES:

        The examples below demonstrate how the caches of ``w`` are built
        step by step using the ``_precompute_h`` method. Thus they rely on
        an untouched Witt symmetric basis that hasn't already seen some
        of its cache filled by other computations. We obtain such a basis
        by choosing a ground ring unlikely to appear elsewhere::

            sage: Sym = SymmetricFunctions(ZZ['hell', 'yeah'])
            sage: w = Sym.Witt()
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(w._h_to_self_cache)
            []
            sage: w._precompute_h(0)
            sage: l(w._h_to_self_cache)
            [([], [([], 1)])]
            sage: w._precompute_h(1)
            sage: l(w._h_to_self_cache)
            [([], [([], 1)]), ([1], [([1], 1)])]
            sage: w._precompute_h(2)
            sage: l(w._h_to_self_cache)
            [([], [([], 1)]),
             ([1], [([1], 1)]),
             ([1, 1], [([1, 1], 1)]),
             ([2], [([1, 1], 1), ([2], 1)])]
            sage: w._h_transition_matrices[2]
            [1 1]
            [0 1]
            sage: w._h_inverse_transition_matrices[2]
            [ 1 -1]
            [ 0  1]
            sage: sorted(w._h_transition_matrices)
            [0, 1, 2]
        """
        l = len(self._h_transition_matrices)
        if l <= n:
            from sage.combinat.partition import Partitions_n
            from sage.misc.cachefunc import cached_function
            @cached_function
            def wsum(m):     # expansion of h_m in w-basis, for m > 0
                return self._from_dict({lam: 1 for lam in Partitions_n(m)})
            for i in range(l, n + 1):
                self._precompute_cache(i, self._h_to_self_cache,
                                       self._h_from_self_cache,
                                       self._h_transition_matrices,
                                       self._h_inverse_transition_matrices,
                                       wsum)

    def _precompute_e(self, n):
        """
        Compute the transition matrices between ``self`` and the elementary
        symmetric basis in the homogeneous components of degree `n`
        (and in those of smaller degree, if not already computed).
        The result is not returned, but rather stored in the cache.

        This assumes that the ``coerce_e`` keyword has been set to
        ``True`` in the initialization of ``self`` (otherwise the cache
        does not exist).

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES:

        The examples below demonstrate how the caches of ``w`` are built
        step by step using the ``_precompute_e`` method. Thus they rely on
        an untouched Witt symmetric basis that hasn't already seen some
        of its cache filled by other computations. We obtain such a basis
        by choosing a ground ring unlikely to appear elsewhere::

            sage: Sym = SymmetricFunctions(ZZ['hell', 'yeah'])
            sage: w = Sym.Witt(coerce_e=True)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(w._e_to_self_cache)
            []
            sage: w._precompute_e(0)
            sage: l(w._e_to_self_cache)
            [([], [([], 1)])]
            sage: w._precompute_e(1)
            sage: l(w._e_to_self_cache)
            [([], [([], 1)]), ([1], [([1], 1)])]
            sage: w._precompute_e(2)
            sage: l(w._e_to_self_cache)
            [([], [([], 1)]),
             ([1], [([1], 1)]),
             ([1, 1], [([1, 1], 1)]),
             ([2], [([2], -1)])]
            sage: w._e_transition_matrices[2]
            [-1  0]
            [ 0  1]
            sage: w._e_inverse_transition_matrices[2]
            [-1  0]
            [ 0  1]
        """
        l = len(self._e_transition_matrices)
        if l <= n:
            from sage.combinat.partition import Partitions
            from sage.misc.cachefunc import cached_function
            @cached_function
            def wsum_e(m):     # expansion of e_m in w-basis, for m > 0
                return self._from_dict({lam: (-1 if (m + len(lam)) % 2 == 1 else 1)
                                        for lam in Partitions(m, max_slope=-1)})
            for i in range(l, n + 1):
                self._precompute_cache(i, self._e_to_self_cache,
                                       self._e_from_self_cache,
                                       self._e_transition_matrices,
                                       self._e_inverse_transition_matrices,
                                       wsum_e)

    def _precompute_p(self, n):
        """
        Compute the transition matrices between ``self`` and the powersum
        basis in the homogeneous components of degree `n`
        (and in those of smaller degree, if not already computed).
        The result is not returned, but rather stored in the cache.

        This assumes that the ``coerce_p`` keyword has been set to
        ``True`` in the initialization of ``self`` (otherwise the cache
        does not exist).

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES:

        The examples below demonstrate how the caches of ``w`` are built
        step by step using the ``_precompute_p`` method. Thus they rely on
        an untouched Witt symmetric basis that hasn't already seen some
        of its cache filled by other computations. We obtain such a basis
        by choosing a ground ring unlikely to appear elsewhere::

            sage: Sym = SymmetricFunctions(QQ['hell', 'yeah'])
            sage: w = Sym.Witt(coerce_h=False, coerce_e=True, coerce_p=True)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(w._p_to_self_cache)
            []
            sage: w._precompute_p(0)
            sage: l(w._p_to_self_cache)
            [([], [([], 1)])]
            sage: w._precompute_p(1)
            sage: l(w._p_to_self_cache)
            [([], [([], 1)]), ([1], [([1], 1)])]
            sage: w._precompute_p(2)
            sage: l(w._p_to_self_cache)
            [([], [([], 1)]), ([1], [([1], 1)]), ([1, 1], [([1, 1], 1)]), ([2], [([1, 1], 1), ([2], 2)])]
            sage: w._p_transition_matrices[2]
            [2 1]
            [0 1]
            sage: w._p_inverse_transition_matrices[2]
            [ 1/2 -1/2]
            [   0    1]
        """
        l = len(self._p_transition_matrices)
        if l <= n:
            from sage.arith.all import divisors
            from sage.combinat.partition import Partition
            from sage.misc.cachefunc import cached_function
            @cached_function
            def wsum_p(m):     # expansion of p_m in w-basis, for m > 0
                return self._from_dict({Partition([d] * (m // d)): d
                                        for d in divisors(m)})
            for i in range(l, n + 1):
                self._precompute_cache(i, self._p_to_self_cache,
                                       self._p_from_self_cache,
                                       self._p_transition_matrices,
                                       self._p_inverse_transition_matrices,
                                       wsum_p)

    def _h_to_w_on_basis(self, lam):
        r"""
        Return the complete homogeneous symmetric function ``h[lam]``
        expanded in the Witt basis, where ``lam`` is a partition.

        This assumes that the ``coerce_h`` keyword has been set to ``True`` in
        the initialization of ``self`` (otherwise the cache does not exist).

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - the expansion of ``h[lam]`` in the Witt basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: h = Sym.homogeneous()
            sage: w = Sym.w()
            sage: w._h_to_w_on_basis(Partition([]))
            w[]
            sage: w._h_to_w_on_basis(Partition([4,2,1]))
            w[1, 1, 1, 1, 1, 1, 1] + 2*w[2, 1, 1, 1, 1, 1] + 2*w[2, 2, 1, 1, 1] + w[2, 2, 2, 1] + w[3, 1, 1, 1, 1] + w[3, 2, 1, 1] + w[4, 1, 1, 1] + w[4, 2, 1]
            sage: h(w._h_to_w_on_basis(Partition([3,1]))) == h[3,1]
            True
        """
        n = sum(lam)
        self._precompute_h(n)
        return self._from_dict(self._h_to_self_cache[lam])

    def _w_to_h_on_basis(self, lam):
        r"""
        Return the Witt symmetric function ``w[lam]``  expanded in the
        complete homogeneous basis, where ``lam`` is a partition.

        This assumes that the ``coerce_h`` keyword has been set to ``True`` in
        the initialization of ``self`` (otherwise the cache does not exist).

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - the expansion of ``w[lam]`` in the complete
          homogeneous basis of ``self.realization_of()``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: h = Sym.homogeneous()
            sage: w = Sym.w()
            sage: w._w_to_h_on_basis(Partition([]))
            h[]
            sage: w._w_to_h_on_basis(Partition([4,2,1]))
            h[1, 1, 1, 1, 1, 1, 1] - 3*h[2, 1, 1, 1, 1, 1] + 3*h[2, 2, 1, 1, 1] - h[2, 2, 2, 1] + h[3, 1, 1, 1, 1] - h[3, 2, 1, 1] - h[4, 1, 1, 1] + h[4, 2, 1]
            sage: w(w._w_to_h_on_basis(Partition([3,1]))) == w[3,1]
            True
        """
        n = sum(lam)
        self._precompute_h(n)
        return self._h._from_dict(self._h_from_self_cache[lam])

    def _e_to_w_on_basis(self, lam):
        r"""
        Return the elementary symmetric function ``e[lam]`` expanded in
        the Witt basis, where ``lam`` is a partition.

        This assumes that the ``coerce_e`` keyword has been set to ``True`` in
        the initialization of ``self`` (otherwise the cache does not exist).

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - the expansion of ``e[lam]`` in the Witt basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: e = Sym.elementary()
            sage: w = Sym.w(coerce_e=True)
            sage: w._e_to_w_on_basis(Partition([]))
            w[]
            sage: w._e_to_w_on_basis(Partition([4,2,1]))
            -w[3, 2, 1, 1] + w[4, 2, 1]
            sage: e(w._e_to_w_on_basis(Partition([3,1]))) == e[3,1]
            True
        """
        n = sum(lam)
        self._precompute_e(n)
        return self._from_dict(self._e_to_self_cache[lam])

    def _w_to_e_on_basis(self, lam):
        r"""
        Return the Witt symmetric function ``w[lam]``
        expanded in the elementary symmetric basis, where
        ``lam`` is a partition.

        This assumes that the ``coerce_e`` keyword has been set to ``True`` in
        the initialization of ``self`` (otherwise the cache does not exist).

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - the expansion of ``w[lam]`` in the elementary
          symmetric basis of ``self.realization_of()``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: e = Sym.elementary()
            sage: w = Sym.w(coerce_e=True)
            sage: w._w_to_e_on_basis(Partition([]))
            e[]
            sage: w._w_to_e_on_basis(Partition([4,2,1]))
            e[2, 2, 1, 1, 1] - e[3, 2, 1, 1] + e[4, 2, 1]
            sage: w(w._w_to_e_on_basis(Partition([3,1]))) == w[3,1]
            True
        """
        n = sum(lam)
        self._precompute_e(n)
        return self._e._from_dict(self._e_from_self_cache[lam])

    def _p_to_w_on_basis(self, lam):
        r"""
        Return the powersum symmetric function ``p[lam]`` expanded in
        the Witt basis, where ``lam`` is a partition.

        This assumes that the ``coerce_p`` keyword has been set to ``True`` in
        the initialization of ``self`` (otherwise the cache does not exist).

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - the expansion of ``p[lam]`` in the Witt basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.power()
            sage: w = Sym.w(coerce_p=True)
            sage: w._p_to_w_on_basis(Partition([]))
            w[]
            sage: w._p_to_w_on_basis(Partition([4,2,1]))
            w[1, 1, 1, 1, 1, 1, 1] + 2*w[2, 1, 1, 1, 1, 1] + 2*w[2, 2, 1, 1, 1] + 4*w[2, 2, 2, 1] + 4*w[4, 1, 1, 1] + 8*w[4, 2, 1]
            sage: p(w._p_to_w_on_basis(Partition([3,1]))) == p[3,1]
            True
        """
        n = sum(lam)
        self._precompute_p(n)
        return self._from_dict(self._p_to_self_cache[lam])

    def _w_to_p_on_basis(self, lam):
        r"""
        Return the Witt symmetric function ``w[lam]`` expanded in the
        powersum basis, where ``lam`` is a  partition.

        This assumes that the ``coerce_p`` keyword has been set to ``True`` in
        the initialization of ``self`` (otherwise the cache does not exist).

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - the expansion of ``w[lam]`` in the powersum
          basis of ``self.realization_of()``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.power()
            sage: w = Sym.w(coerce_p=True)
            sage: w._w_to_p_on_basis(Partition([]))
            p[]
            sage: w._w_to_p_on_basis(Partition([4,2,1]))
            3/16*p[1, 1, 1, 1, 1, 1, 1] - 5/16*p[2, 1, 1, 1, 1, 1] + 3/16*p[2, 2, 1, 1, 1] - 1/16*p[2, 2, 2, 1] - 1/8*p[4, 1, 1, 1] + 1/8*p[4, 2, 1]
            sage: w(w._w_to_p_on_basis(Partition([3,1]))) == w[3,1]
            True
        """
        n = sum(lam)
        self._precompute_p(n)
        return self._p._from_dict(self._p_from_self_cache[lam])

    def __init_extra__(self):
        """
        Sets up caches for the transition maps to other bases, and registers
        them as coercions.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ) # indirect doctest
            sage: h = Sym.h(); w = Sym.w()

            sage: phi = h.coerce_map_from(w); phi
            Generic morphism:
              From: Symmetric Functions over Rational Field in the Witt basis
              To:   Symmetric Functions over Rational Field in the homogeneous basis
            sage: phi(w.an_element()) == h(w.an_element())
            True
            sage: e = Sym.e(); w2 = Sym.w(coerce_e=True)
            sage: psi = e.coerce_map_from(w2); psi
            Generic morphism:
              From: Symmetric Functions over Rational Field in the Witt basis
              To:   Symmetric Functions over Rational Field in the elementary basis
            sage: psi(w2.an_element()) == e(w2.an_element())
            True
        """

        #category = sage.categories.all.ModulesWithBasis(self.base_ring())

        # Set up coercions and conversions with appropriate other bases.
        # self._p, self._e and self._h will be the powersum basis, the elementary
        # symmetric basis and the complete homogeneous basis (over the same base
        # ring as self), respectively (but they are only set if the respective
        # arguments ``coerce_p``, ``coerce_e`` and ``coerce_h`` are True).
        # self._friendly will be the one available basis which makes computations
        # the easiest.

        self._friendly = None

        if self._coerce_p:
            self._p = self.realization_of().p()
            # Set up the cache for conversion from the Witt basis
            # to the powersum basis.

            # cache for the coordinates of the elements
            # of the powersum basis with respect to the Witt basis
            self._p_to_self_cache = {}
            # cache for the coordinates of the elements
            # of the Witt basis with respect to the powersum basis
            self._p_from_self_cache = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the powersum basis with respect to the Witt basis
            self._p_transition_matrices = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the Witt basis with respect to the powersum basis
            self._p_inverse_transition_matrices = {}

            self   .register_coercion(self._p._module_morphism(self._p_to_w_on_basis, codomain = self))
            from sage.rings.rational_field import RationalField
            if self.base_ring().has_coerce_map_from(RationalField):
                self._p.register_coercion(self._module_morphism(self._w_to_p_on_basis, codomain = self._p))
                self._friendly = self._p
            else:
                # self._w_to_p_on_basis is a partial map at best
                self._p.register_conversion(self._module_morphism(self._w_to_p_on_basis, codomain = self._p))
                if (not self._coerce_e) and (not self._coerce_h):
                    # ensure that self has coercion at least to one other basis,
                    # or else coercion-based computations will fail
                    self._coerce_h = True
        elif (not self._coerce_e) and (not self._coerce_h):
            self._coerce_h = True     # at least one coercion is needed!

        if self._coerce_h:
            self._h = self.realization_of().h()
            # Set up the cache for conversion from the Witt basis to the complete
            # homogeneous basis. (This is the conversion that is used by default.)

            # cache for the coordinates of the elements
            # of the homogeneous basis with respect to the Witt basis
            self._h_to_self_cache = {}
            # cache for the coordinates of the elements
            # of the Witt basis with respect to the homogeneous basis
            self._h_from_self_cache = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the homogeneous basis with respect to the Witt basis
            self._h_transition_matrices = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the Witt basis with respect to the homogeneous basis
            self._h_inverse_transition_matrices = {}
            self   .register_coercion(self._h._module_morphism(self._h_to_w_on_basis, codomain = self))
            self._h.register_coercion(self._module_morphism(self._w_to_h_on_basis, codomain = self._h))
            if self._friendly is None:
                self._friendly = self._h

        if self._coerce_e:
            self._e = self.realization_of().e()
            # Set up the cache for conversion from the Witt basis to the elementary
            # symmetric basis.

            # cache for the coordinates of the elements
            # of the elementary basis with respect to the Witt basis
            self._e_to_self_cache = {}
            # cache for the coordinates of the elements
            # of the Witt basis with respect to the elementary basis
            self._e_from_self_cache = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the elementary basis with respect to the Witt basis
            self._e_transition_matrices = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the Witt basis with respect to the elementary basis
            self._e_inverse_transition_matrices = {}
            self   .register_coercion(self._e._module_morphism(self._e_to_w_on_basis, codomain = self))
            self._e.register_coercion(self._module_morphism(self._w_to_e_on_basis, codomain = self._e))
            if self._friendly is None:
                self._friendly = self._e

    def from_other_uncached(self, u):
        r"""
        Return an element ``u`` of another basis of the ring of
        symmetric functions, expanded in the Witt basis ``self``.
        The result is the same as ``self(u)``, but the
        ``from_other_uncached`` method does not precompute a
        cache with transition matrices. Thus,
        ``from_other_uncached`` is faster when ``u`` is sparse.

        INPUT:

        - ``u`` -- an element of ``self.realization_of()``

        OUTPUT:

        - the expansion of ``u`` in the Witt basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.p()
            sage: w = Sym.w()
            sage: a = p([3,2]) - p([4,1]) + 27 * p([3])
            sage: w.from_other_uncached(a) == w(a)
            True

        Here's a verification of an obvious fact that would take
        long with regular coercion::

            sage: fouc = w.from_other_uncached
            sage: fouc(p([15]))
            w[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] + 3*w[3, 3, 3, 3, 3] + 5*w[5, 5, 5] + 15*w[15]
            sage: fouc(p([15])) * fouc(p([14])) == fouc(p([15, 14]))
            True

        Other bases::

            sage: e = Sym.e()
            sage: h = Sym.h()
            sage: s = Sym.s()
            sage: all( fouc(e(lam)) == w(e(lam)) for lam in Partitions(5) )
            True
            sage: all( fouc(h(lam)) == w(h(lam)) for lam in Partitions(5) )
            True
            sage: all( fouc(p(lam)) == w(p(lam)) for lam in Partitions(5) )
            True
            sage: all( fouc(s(lam)) == w(s(lam)) for lam in Partitions(5) )
            True
        """
        parent_name = u.parent().basis_name()
        from sage.misc.cachefunc import cached_function

        if parent_name == "homogeneous":
            from sage.combinat.partition import Partitions_n
            @cached_function
            def wsum(m):     # expansion of h_m in w-basis, for m > 0
                return self._from_dict({lam: 1 for lam in Partitions_n(m)})
            result = self.zero()
            for lam, a in u.monomial_coefficients().items():
                product = self.one()
                for i in lam:
                    product *= wsum(i)
                result += a * product
            return result

        if parent_name == "powersum":
            from sage.arith.all import divisors
            from sage.combinat.partition import Partition
            @cached_function
            def wsum_p(m):     # expansion of p_m in w-basis, for m > 0
                return self._from_dict({Partition([d] * (m // d)): d
                                        for d in divisors(m)})
            result = self.zero()
            for lam, a in u.monomial_coefficients().items():
                product = self.one()
                for i in lam:
                    product *= wsum_p(i)
                result += a * product
            return result

        # Coerce u into elementary symmetric basis.
        if parent_name != "elementary":
            u = u.parent().realization_of().elementary()(u)

        from sage.combinat.partition import Partitions
        @cached_function
        def wsum_e(m):     # expansion of e_m in w-basis, for m > 0
            return self._from_dict({lam: (-1 if (m + len(lam)) % 2 == 1 else 1)
                                    for lam in Partitions(m, max_slope=-1)})
        result = self.zero()
        for lam, a in u.monomial_coefficients().items():
            product = self.one()
            for i in lam:
                product *= wsum_e(i)
            result += a * product
        return result

    def coproduct(self, elt):
        r"""
        Return the coproduct of the element ``elt``.

        INPUT:

        - ``elt`` -- a symmetric function written in this basis

        OUTPUT:

        - The coproduct acting on ``elt``; the result is an element of the
          tensor squared of the basis ``self``

        EXAMPLES::

            sage: w = SymmetricFunctions(QQ).w()
            sage: w[2].coproduct()
            w[] # w[2] - w[1] # w[1] + w[2] # w[]
            sage: w.coproduct(w[2])
            w[] # w[2] - w[1] # w[1] + w[2] # w[]
            sage: w[2,1].coproduct()
            w[] # w[2, 1] - w[1] # w[1, 1] + w[1] # w[2] - w[1, 1] # w[1] + w[2] # w[1] + w[2, 1] # w[]
            sage: w.coproduct(w[2,1])
            w[] # w[2, 1] - w[1] # w[1, 1] + w[1] # w[2] - w[1, 1] # w[1] + w[2] # w[1] + w[2, 1] # w[]

        TESTS:

        The same, but with other settings::

            sage: w = SymmetricFunctions(QQ).w(coerce_h=False, coerce_e=True)
            sage: w[2].coproduct()
            w[] # w[2] - w[1] # w[1] + w[2] # w[]
            sage: w.coproduct(w[2])
            w[] # w[2] - w[1] # w[1] + w[2] # w[]
            sage: w[2,1].coproduct()
            w[] # w[2, 1] - w[1] # w[1, 1] + w[1] # w[2] - w[1, 1] # w[1] + w[2] # w[1] + w[2, 1] # w[]
            sage: w.coproduct(w[2,1])
            w[] # w[2, 1] - w[1] # w[1, 1] + w[1] # w[2] - w[1, 1] # w[1] + w[2] # w[1] + w[2, 1] # w[]

            sage: w = SymmetricFunctions(QQ).w(coerce_h=False, coerce_p=True)
            sage: w[2].coproduct()
            w[] # w[2] - w[1] # w[1] + w[2] # w[]
            sage: w.coproduct(w[2])
            w[] # w[2] - w[1] # w[1] + w[2] # w[]
            sage: w[2,1].coproduct()
            w[] # w[2, 1] - w[1] # w[1, 1] + w[1] # w[2] - w[1, 1] # w[1] + w[2] # w[1] + w[2, 1] # w[]
            sage: w.coproduct(w[2,1])
            w[] # w[2, 1] - w[1] # w[1, 1] + w[1] # w[2] - w[1, 1] # w[1] + w[2] # w[1] + w[2, 1] # w[]
        """
        from sage.categories.tensor import tensor
        friendly = self._friendly
        return self.tensor_square().sum(coeff * tensor([self(friendly[x]), self(friendly[y])])
                                        for ((x,y), coeff) in friendly(elt).coproduct())

    def verschiebung(self, n):
        r"""
        Return the image of the symmetric function ``self`` under the
        `n`-th Verschiebung operator.

        The `n`-th Verschiebung operator `\mathbf{V}_n` is defined to be
        the unique algebra endomorphism `V` of the ring of symmetric
        functions that satisfies `V(h_r) = h_{r/n}` for every positive
        integer `r` divisible by `n`, and satisfies `V(h_r) = 0` for
        every positive integer `r` not divisible by `n`. This operator
        `\mathbf{V}_n` is a Hopf algebra endomorphism. For every
        nonnegative integer `r` with `n \mid r`, it satisfies

        .. MATH::

            \mathbf{V}_n(h_r) = h_{r/n},
            \quad \mathbf{V}_n(p_r) = n p_{r/n},
            \quad \mathbf{V}_n(e_r) = (-1)^{r - r/n} e_{r/n},
            \quad \mathbf{V}_n(w_r) = w_{r/n},

        (where `h` is the complete homogeneous basis, `p` is the
        powersum basis, `e` is the elementary basis, and `w` is the
        Witt basis). For every nonnegative integer `r` with `n \nmid r`,
        it satisfes

        .. MATH::

            \mathbf{V}_n(h_r) = \mathbf{V}_n(p_r) = \mathbf{V}_n(e_r)
            = \mathbf{V}_n(w_r) = 0.

        The `n`-th Verschiebung operator is also called the `n`-th
        Verschiebung endomorphism. Its name derives from the Verschiebung
        (German for "shift") endomorphism of the Witt vectors.

        The `n`-th Verschiebung operator is adjoint to the `n`-th
        Frobenius operator (see :meth:`frobenius` for its definition)
        with respect to the Hall scalar product (:meth:`scalar`).

        The action of the `n`-th Verschiebung operator on the Schur basis
        can also be computed explicitly. The following (probably clumsier
        than necessary) description can be obtained by solving exercise
        7.61 in Stanley's [STA]_.

        Let `\lambda` be a partition. Let `n` be a positive integer. If
        the `n`-core of `\lambda` is nonempty, then
        `\mathbf{V}_n(s_\lambda) = 0`. Otherwise, the following method
        computes `\mathbf{V}_n(s_\lambda)`: Write the partition `\lambda`
        in the form `(\lambda_1, \lambda_2, \ldots, \lambda_{ns})` for some
        nonnegative integer `s`. (If `n` does not divide the length of
        `\lambda`, then this is achieved by adding trailing zeroes to
        `\lambda`.) Set `\beta_i = \lambda_i + ns - i` for every
        `s \in \{ 1, 2, \ldots, ns \}`. Then,
        `(\beta_1, \beta_2, \ldots, \beta_{ns})` is a strictly decreasing
        sequence of nonnegative integers. Stably sort the list
        `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
        `-1 - \beta_i` modulo `n`. Let `\xi` be the sign of the
        permutation that is used for this sorting. Let `\psi` be the sign
        of the permutation that is used to stably sort the list
        `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
        `i - 1` modulo `n`. (Notice that `\psi = (-1)^{n(n-1)s(s-1)/4}`.)
        Then, `\mathbf{V}_n(s_\lambda) = \xi \psi \prod_{i = 0}^{n - 1}
        s_{\lambda^{(i)}}`, where
        `(\lambda^{(0)}, \lambda^{(1)}, \ldots, \lambda^{(n - 1)})`
        is the `n`-quotient of `\lambda`.

        INPUT:

        - ``n`` -- a positive integer

        OUTPUT:

        The result of applying the `n`-th Verschiebung operator (on the ring of
        symmetric functions) to ``self``.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: w = Sym.w()
            sage: w[3].verschiebung(2)
            0
            sage: w[4].verschiebung(4)
            w[1]

        TESTS:

        Let us check that this method on the Witt basis gives the
        same result as the implementation in sfa.py on the complete
        homogeneous basis::

            sage: Sym = SymmetricFunctions(QQ)
            sage: w = Sym.w(); h = Sym.h()
            sage: all( w(h(lam)).verschiebung(3) == w(h(lam).verschiebung(3))
            ....:      for lam in Partitions(6) )
            True
            sage: all( h(w(lam)).verschiebung(2) == h(w(lam).verschiebung(2))
            ....:      for lam in Partitions(4) )
            True
        """
        parent = self.parent()
        w_coords_of_self = self.monomial_coefficients().items()
        from sage.combinat.partition import Partition
        dct = {Partition([i // n for i in lam]): coeff
               for (lam, coeff) in w_coords_of_self
               if all( i % n == 0 for i in lam )}
        result_in_w_basis = parent._from_dict(dct)
        return parent(result_in_w_basis)


