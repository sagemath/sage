# sage.doctest: needs sage.modules
"""
Yokonuma-Hecke Algebras

AUTHORS:

- Travis Scrimshaw (2015-11): initial version
- Travis Scrimshaw (2025-03): general type version
"""

# ****************************************************************************
#  Copyright (C) 2015-2025 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.rational_field import QQ
from sage.categories.algebras import Algebras
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations
from sage.sets.family import Family


class YokonumaHeckeAlgebra(CombinatorialFreeModule):
    r"""
    Abstract base class for Yokonuma-Hecke algebras that
    implements common features.

    .. TODO::

        Factor out the near-common features.
    """
    @staticmethod
    def __classcall_private__(cls, d, n, q=None, R=None):
        r"""
        Standardize input to ensure a unique representation and dispatch
        to the correct implementation.

        TESTS::

            sage: Y1 = algebras.YokonumaHecke(5, 3)
            sage: q = LaurentPolynomialRing(QQ, 'q').gen()
            sage: Y2 = algebras.YokonumaHecke(5, 3, q)
            sage: Y3 = algebras.YokonumaHecke(5, 3, q, q.parent())
            sage: Y1 is Y2 and Y2 is Y3
            True
        """
        if q is None:
            q = LaurentPolynomialRing(QQ, 'q').gen()
        if R is None:
            R = q.parent()
        q = R(q)
        if R not in Rings().Commutative():
            raise TypeError("base ring must be a commutative ring")
        if n not in ZZ:
            from sage.combinat.root_system.cartan_type import CartanType
            n = CartanType(n)
            return YokonumaHeckeAlgebraWeyl(d, n, q, R)
        return YokonumaHeckeAlgebraGL(d, n, q, R)

    def __init__(self, d, W, q, R, indices, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: elts = Y.some_elements() + list(Y.algebra_generators())
            sage: TestSuite(Y).run(elements=elts)
        """
        self._d = d
        self._W = W
        self._cartan_type = W.cartan_type()
        self._q = q
        cat = Algebras(R).WithBasis().or_subcategory(category)
        CombinatorialFreeModule.__init__(self, R, indices, prefix='Y',
                                         category=cat)
        self._assign_names(self.algebra_generators().keys())

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['F',4])
            sage: Y.cartan_type()
            ['F', 4]
        """
        return self._cartan_type

    def index_set(self):
        r"""
        Return the index set of ``self``, which is the index set of
        the Cartan type of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['F',4])
            sage: Y.index_set() == Y.cartan_type().index_set()
            True
        """
        return self._cartan_type.index_set()

    def q(self):
        r"""
        Return the parameter `q` of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['F',4])
            sage: Y.q()
            q
            sage: Y.q().parent() is Y.base_ring()
            True
        """
        return self._q

    def g(self, i=None):
        """
        Return the generator(s) `g_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `g_i` or if ``None``,
          then the family of all generators `g_i`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(8, 3)
            sage: Y.g(1)
            g[1]
            sage: list(Y.g())
            [g[1], g[2]]

            sage: Y = algebras.YokonumaHecke(8, ['G',2])
            sage: Y.g(1)
            g[1]
            sage: Y.g()
            Finite family {1: g[1], 2: g[2]}
        """
        G = self.algebra_generators()
        if i is None:
            I = self._W.index_set()
            d = {i: G['g%s' % i] for i in I}
            return Family(I, d.__getitem__)
        return G['g%s' % i]

    @cached_method
    def gens(self) -> tuple:
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: Y.gens()
            (g[1], g[2], t1, t2, t3)

            sage: Y = algebras.YokonumaHecke(5, ['B',2])
            sage: Y.gens()
            (g[1], g[2], h1, h2)
        """
        return tuple(self.algebra_generators())


class YokonumaHeckeAlgebraGL(YokonumaHeckeAlgebra):
    r"""
    The Yokonuma-Hecke algebra `Y_{d,n}(q)` for `GL_n(\GF{d})`.

    Let `R` be a commutative ring and `q` be a unit in `R`. The
    *Yokonuma-Hecke algebra* `Y_{d,n}(q)` is the associative, unital
    `R`-algebra generated by `t_1, t_2, \ldots, t_n, g_1, g_2, \ldots,
    g_{n-1}` and subject to the relations:

    - `g_i g_j = g_j g_i` for all `|i - j| > 1`,
    - `g_i g_{i+1} g_i = g_{i+1} g_i g_{i+1}`,
    - `t_i t_j = t_j t_i`,
    - `t_j g_i = g_i t_{j s_i}`, and
    - `t_j^d = 1`,

    where `s_i` is the simple transposition `(i, i+1)`, along with
    the quadratic relation

    .. MATH::

        g_i^2 = 1 + \frac{(q - q^{-1})}{d} \left( \sum_{s=0}^{d-1}
        t_i^s t_{i+1}^{-s} \right) g_i.

    Thus the Yokonuma-Hecke algebra can be considered a quotient of
    the framed braid group `(\ZZ / d\ZZ) \wr B_n`, where `B_n` is the
    classical braid group on `n` strands, by the quadratic relations.
    Moreover, all of the algebra generators are invertible. In
    particular, we have

    .. MATH::

        g_i^{-1} = g_i - (q - q^{-1}) e_i.

    When we specialize `q = \pm 1`, we obtain the group algebra of
    the complex reflection group `G(d, 1, n) = (\ZZ / d\ZZ) \wr S_n`.
    Moreover for `d = 1`, the Yokonuma-Hecke algebra is equal to the
    :class:`Iwahori-Hecke <IwahoriHeckeAlgebra>` of type `A_{n-1}`.

    This was considered for more general Chevalley groups (Lie groups
    over finite fields); see :class:`YokonumaHeckeAlgebraWeyl`.

    INPUT:

    - ``d`` -- the maximum power of `t`
    - ``n`` -- the number of generators or a Cartan type
    - ``q`` -- (optional) an invertible element in a commutative ring;
      the default is `q \in \QQ[q,q^{-1}]`
    - ``R`` -- (optional) a commutative ring containing ``q``; the
      default is the parent of `q`

    EXAMPLES:

    We construct `Y_{4,3}` and do some computations::

        sage: Y = algebras.YokonumaHecke(4, 3)
        sage: g1, g2, t1, t2, t3 = Y.algebra_generators()
        sage: g1 * g2
        g[1,2]
        sage: t1 * g1
        t1*g[1]
        sage: g2 * t2
        t3*g[2]
        sage: g2 * t3
        t2*g[2]
        sage: (g2 + t1) * (g1 + t2*t3)
        g[2,1] + t2*t3*g[2] + t1*g[1] + t1*t2*t3
        sage: g1 * g1
        1 - (1/4*q^-1-1/4*q)*g[1] - (1/4*q^-1-1/4*q)*t1*t2^3*g[1]
         - (1/4*q^-1-1/4*q)*t1^2*t2^2*g[1] - (1/4*q^-1-1/4*q)*t1^3*t2*g[1]
        sage: g2 * g1 * t1
        t3*g[2,1]

    We construct the elements `e_i` and show that they are idempotents::

        sage: e1 = Y.e(1); e1
        1/4 + 1/4*t1*t2^3 + 1/4*t1^2*t2^2 + 1/4*t1^3*t2
        sage: e1 * e1 == e1
        True
        sage: e2 = Y.e(2); e2
        1/4 + 1/4*t2*t3^3 + 1/4*t2^2*t3^2 + 1/4*t2^3*t3
        sage: e2 * e2 == e2
        True

    REFERENCES:

    - [CL2013]_
    - [CPdA2014]_
    - [ERH2015]_
    - [JPdA15]_
    """
    def __init__(self, d, n, q, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: elts = Y.some_elements() + list(Y.algebra_generators())
            sage: TestSuite(Y).run(elements=elts)
        """
        self._n = n
        W = Permutations(n)
        import itertools
        C = itertools.product(*([range(d)] * n))
        indices = list(itertools.product(C, W))
        YokonumaHeckeAlgebra.__init__(self, d, W, q, R, indices)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.YokonumaHecke(5, 2)
            Yokonuma-Hecke algebra of rank 5 and order 2 with q=q
             over Univariate Laurent Polynomial Ring in q over Rational Field
        """
        return "Yokonuma-Hecke algebra of rank {} and order {} with q={} over {}".format(
            self._d, self._n, self._q, self.base_ring())

    def _latex_(self) -> str:
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 2)
            sage: latex(Y)
            \mathcal{Y}_{5,2}(q)
        """
        return "\\mathcal{Y}_{%s,%s}(%s)" % (self._d, self._n, self._q)

    def _repr_term(self, m) -> str:
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: Y._repr_term( ((1, 0, 2), Permutation([3,2,1])) )
            't1*t3^2*g[2,1,2]'
        """
        def gen_str(e):
            return '' if e == 1 else '^%s' % e
        lhs = '*'.join('t%s' % (j + 1) + gen_str(i)
                       for j, i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        rhs = 'g[{}]'.format(','.join(str(i) for i in redword))
        if not lhs:
            return rhs
        return lhs + '*' + rhs

    def _latex_term(self, m) -> str:
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: Y._latex_term( ((1, 0, 2), Permutation([3,2,1])) )
            't_{1} t_{3}^{2} g_{2} g_{1} g_{2}'
        """
        def gen_str(e):
            return '' if e == 1 else '^{%s}' % e
        lhs = ' '.join('t_{%s}' % (j + 1) + gen_str(i)
                       for j, i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        return lhs + ' ' + ' '.join("g_{%d}" % i for i in redword)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: dict(Y.algebra_generators())
            {'g1': g[1], 'g2': g[2], 't1': t1, 't2': t2, 't3': t3}
        """
        one = self._W.one()
        zero = [0] * self._n
        d = {}
        for i in range(self._n):
            r = list(zero) # Make a copy
            r[i] = 1
            d['t%s' % (i+1)] = self.monomial((tuple(r), one))
        G = self._W.group_generators()
        for i in range(1, self._n):
            d['g%s' % i] = self.monomial((tuple(zero), G[i]))
        return Family(sorted(d), lambda i: d[i])

    @cached_method
    def one_basis(self):
        """
        Return the index of the basis element of `1`.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: Y.one_basis()
            ((0, 0, 0), [1, 2, 3])
        """
        one = self._W.one()
        zero = (0,) * self._n
        return (zero, one)

    @cached_method
    def e(self, i):
        """
        Return the element `e_i`.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: Y.e(1)
            1/4 + 1/4*t1*t2^3 + 1/4*t1^2*t2^2 + 1/4*t1^3*t2
            sage: Y.e(2)
            1/4 + 1/4*t2*t3^3 + 1/4*t2^2*t3^2 + 1/4*t2^3*t3
        """
        if i < 1 or i >= self._n:
            raise ValueError("invalid index")
        c = ~self.base_ring()(self._d)
        zero = [0]*self._n
        one = self._W.one()
        d = {}
        for s in range(self._d):
            r = list(zero) # Make a copy
            r[i-1] = s
            if s != 0:
                r[i] = self._d - s
            d[(tuple(r), one)] = c
        return self._from_dict(d, remove_zeros=False)

    def t(self, i=None):
        """
        Return the generator(s) `t_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `t_i` or if ``None``,
          then the family of all generators `t_i`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(8, 3)
            sage: Y.t(2)
            t2
            sage: list(Y.t())
            [t1, t2, t3]
        """
        G = self.algebra_generators()
        if i is None:
            I = tuple(range(1, self._n+1))
            d = {i: G['t%s' % i] for i in I}
            return Family(I, d.__getitem__)
        return G['t%s' % i]

    def product_on_basis(self, m1, m2):
        """
        Return the product of the basis elements indexed by ``m1`` and ``m2``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: 4 * Y.product_on_basis(m, m)
            -(q^-1-q)*t2^2*g[1] + 4*t1*t2 - (q^-1-q)*t1*t2*g[1]
             - (q^-1-q)*t1^2*g[1] - (q^-1-q)*t1^3*t2^3*g[1]

        Check that we apply the permutation correctly on `t_i`::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: g1, g2, t1, t2, t3 = Y.algebra_generators()
            sage: g21 = g2 * g1
            sage: g21 * t1
            t3*g[2,1]
        """
        t1,g1 = m1
        t2,g2 = m2
        # Commute g1 and t2, then multiply t1 and t2
        # ig1 = g1
        t = [(t1[i] + t2[g1.index(i+1)]) % self._d for i in range(self._n)]
        one = self._W.one()
        if g1 == one:
            return self.monomial((tuple(t), g2))
        ret = self.monomial((tuple(t), g1))
        # We have to reverse the reduced word due to Sage's convention
        #   for permutation multiplication
        for i in g2.reduced_word():
            ret = self.linear_combination((self._product_by_basis_gen(m, i), c)
                                          for m,c in ret)
        return ret

    def _product_by_basis_gen(self, m, i):
        r"""
        Return the product `t g_w g_i`.

        If the quadratic relation is `g_i^2 = 1 + (q - q^{-1}) e_i g_i`,
        then we have

        .. MATH::

            g_w g_i = \begin{cases}
            g_{ws_i} & \text{if } \ell(ws_i) = \ell(w) + 1, \\
            g_{ws_i} + (q - q^{-1}) g_w e_i & \text{if }
            \ell(w s_i) = \ell(w) - 1.
            \end{cases}

        INPUT:

        - ``m`` -- a pair ``[t, w]``, where ``t`` encodes the monomial
          and ``w``  is an element of the permutation group
        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: 4 * Y._product_by_basis_gen(m, 1)
            -(q^-1-q)*t2*t3^2*g[1] + 4*t1*t3^2 - (q^-1-q)*t1*t3^2*g[1]
             - (q^-1-q)*t1^2*t2^3*t3^2*g[1] - (q^-1-q)*t1^3*t2^2*t3^2*g[1]
        """
        t, w = m
        wi = w.apply_simple_reflection(i, side='right')
        if not w.has_descent(i, side='right'):
            return self.monomial((t, wi))

        R = self.base_ring()
        c = (self._q - ~self._q) * ~R(self._d)
        d = {(t, wi): R.one()}
        # We commute g_w and e_i and then multiply by t
        for s in range(self._d):
            r = list(t)
            r[w[i-1]-1] = (r[w[i-1]-1] + s) % self._d
            if s != 0:
                r[w[i]-1] = (r[w[i]-1] + self._d - s) % self._d
            d[(tuple(r), w)] = c
        return self._from_dict(d, remove_zeros=False)

    @cached_method
    def inverse_g(self, i):
        r"""
        Return the inverse of the generator `g_i`.

        From the quadratic relation, we have

        .. MATH::

            g_i^{-1} = g_i - (q - q^{-1}) e_i.

        INPUT:

        - ``i`` -- (default: ``None``) the inverse generator `g_i^{-1}` or
          if ``None``, then the family of all inverse generators `g_i^{-1}`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(2, 4)
            sage: [2*Y.inverse_g(i) for i in range(1, 4)]
            [(q^-1-q) + 2*g[1] + (q^-1-q)*t1*t2,
             (q^-1-q) + 2*g[2] + (q^-1-q)*t2*t3,
             (q^-1-q) + 2*g[3] + (q^-1-q)*t3*t4]
            sage: all(Y.inverse_g(i) * Y.g(i) == Y.one() for i in range(1, 4))
            True
            sage: all(Y.g(i) * Y.inverse_g(i) == Y.one() for i in range(1, 4))
            True
        """
        if i is None:
            I = self._W.index_set()
            d = {i: self.inverse_g(i) for i in I}
            return Family(I, d.__getitem__)
        if i < 1 or i >= self._n:
            raise ValueError("invalid index")
        return self.g(i) + (~self._q - self._q) * self.e(i)

    class Element(CombinatorialFreeModule.Element):
        def __invert__(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: Y = algebras.YokonumaHecke(3, 3)
                sage: t = prod(Y.t()); t
                t1*t2*t3
                sage: t.inverse()   # indirect doctest
                t1^2*t2^2*t3^2
                sage: [3*~(t*g) for g in Y.g()]
                [(q^-1-q)*t2*t3^2 + (q^-1-q)*t1*t3^2
                   + (q^-1-q)*t1^2*t2^2*t3^2 + 3*t1^2*t2^2*t3^2*g[1],
                 (q^-1-q)*t1^2*t3 + (q^-1-q)*t1^2*t2
                   + (q^-1-q)*t1^2*t2^2*t3^2 + 3*t1^2*t2^2*t3^2*g[2]]
                sage: g = prod(Y.g())
                sage: ~g * g == Y.one()
                True
                sage: g * ~g == Y.one()
                True

                sage: tp = t * Y.t(2)
                sage: all(tp*g * ~(tp*g) == Y.one() for g in Y.g())
                True

            TESTS:

            Check that :issue:`26424` is fixed::

                sage: Y = algebras.YokonumaHecke(3, 3)
                sage: t = 3 * prod(Y.t())
                sage: ~t
                1/3*t1^2*t2^2*t3^2

                sage: ~Y.zero()
                Traceback (most recent call last):
                ...
                ZeroDivisionError
            """
            if not self:
                raise ZeroDivisionError
            if len(self) != 1:
                raise NotImplementedError("inverse only implemented for basis elements (monomials in the generators)" % self)
            H = self.parent()
            t, w = self.support_of_term()
            c = ~self.coefficients()[0]
            telt = H.monomial((tuple((H._d - e) % H._d for e in t), H._W.one()))
            return c * H.prod(H.inverse_g(i) for i in reversed(w.reduced_word())) * telt


class YokonumaHeckeAlgebraWeyl(YokonumaHeckeAlgebra):
    r"""
    The Yokonuma-Hecke algebra associated to a Cartan type.

    Let `R` be a commutative ring and `q` be a unit in `R`. Let
    `W` be the Weyl group acting on a root lattice `Q`. The
    *Yokonuma-Hecke algebra* `Y_{d,W}(q)` is the associative, unital
    `R`-algebra generated by `\{h_i, g_i \mid i \in I\}`, where `I` is
    the index set of simple roots of `Q`, and subject to the relations:

    - `g_i` and `g_j` satisfy the braid relations of the corresponding
      simple reflections `s_i` and `s_j` in `W`,
    - `h_i h_j = h_j h_i`,
    - `h_j g_i = g_i (s_i \cdot h_j)` with considering `h_j` as the simple
      root `\alpha_j \in Q`, and
    - `h_j^d = 1`,

    along with the quadratic relation

    .. MATH::

        g_i^2 = 1 + (q - 1) e_i (1 + g_i),
        \qquad\qquad
        e_i := \frac{1}{d} \sum_{s=0}^{d-1} h_i^s.

    In particular, we can identify the subalgebra generated by `\{h_i \mid
    i \in I\}` with `(\ZZ / d \ZZ) \otimes_{\ZZ} Q`. The Yokonuma-Hecke
    algebra, when `d = p^m - 1` for a prime `p` and some `m \geq 1`, can
    be identified with functions invariant under the left *and* right actions
    of the unipotent group `U` on `G(\GF{d})`, the semisimple Chevalley
    (or Lie) group associated with `W`. Moreover, all of the algebra
    generators are invertible. In particular, we have

    .. MATH::

        g_i^{-1} = g_i + (q^{-1} - 1) e_i (1 + g_i).

    For `d = 1`, the Yokonuma-Hecke algebra is equal to the
    :class:`Iwahori-Hecke <IwahoriHeckeAlgebra>` of `W`.

    INPUT:

    - ``d`` -- the maximum power of `t`
    - ``ct`` -- the Cartan type
    - ``q`` -- (optional) an invertible element in a commutative ring;
      the default is `q \in \QQ[q,q^{-1}]`
    - ``R`` -- (optional) a commutative ring containing ``q``; the
      default is the parent of `q`

    .. WARNING::

        For type `A_n`, this returns the Yokonuma-Hecke algebra associated
        to the Lie (or Chevalley) group `SL_n(\GF{d})`. For the Yokonuma-Hecke
        algebra corresponding to the (reductive) Lie group `GL_n(\GF{d})`, use
        :class:`YokonumaHeckeAlgebraGL`. Additionally, this uses a different
        quadratic relation.

    REFERENCES:

    - [Marin2018]_
    """
    def __init__(self, d, ct, q, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(2, ['F',4])
            sage: TestSuite(Y).run()
            sage: Y = algebras.YokonumaHecke(3, ['G',2])
            sage: elts = list(Y.gens()) + [Y.an_element()] + [sum(Y.gens())]
            sage: TestSuite(Y).run(elements=elts)  # long time
        """
        from sage.categories.sets_cat import cartesian_product
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

        self._Q = ct.root_system().root_space(IntegerModRing(d))
        self._Qp = ct.root_system().root_lattice()
        W = self._Qp.weyl_group(prefix='s')
        indices = cartesian_product([self._Q, W])
        YokonumaHeckeAlgebra.__init__(self, d, W, q, R, indices)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.YokonumaHecke(5, ['E',6])
            Yokonuma-Hecke algebra of rank 5 for ['E', 6] with q=q
             over Univariate Laurent Polynomial Ring in q over Rational Field
        """
        return "Yokonuma-Hecke algebra of rank {} for {} with q={} over {}".format(
            self._d, self._cartan_type, self._q, self.base_ring())

    def _latex_(self) -> str:
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, ['E',6])
            sage: latex(Y)
            \mathcal{Y}_{5,E_6}(q)
        """
        from sage.misc.latex import latex
        return "\\mathcal{Y}_{%s,%s}(%s)" % (self._d, latex(self._cartan_type), self._q)

    def _repr_term(self, m) -> str:
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['E',6])
            sage: al = Y._Q.simple_root(1) + 3*Y._Q.simple_root(5)
            sage: Y._repr_term((al, prod(Y._W.gens())))
            'h1*h5^3*g[1,3,2,4,5,6]'
        """
        def gen_str(e):
            return '' if e == 1 else '^%s' % e

        I = self._cartan_type.index_set()
        lhs = '*'.join('h%s' % j + gen_str(m[0][j]) for j in I if m[0][j])
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        rhs = 'g[{}]'.format(','.join(str(i) for i in redword))
        if not lhs:
            return rhs
        return lhs + '*' + rhs

    def _latex_term(self, m) -> str:
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['E',6])
            sage: al = Y._Q.simple_root(1) + 3*Y._Q.simple_root(5)
            sage: Y._latex_term((al, prod(Y._W.gens())))
            'h_{1} h_{5}^{3} g_{1} g_{3} g_{2} g_{4} g_{5} g_{6}'
        """
        def gen_str(e):
            return '' if e == 1 else '^{%s}' % e

        I = self._cartan_type.index_set()
        lhs = ' '.join('h_{%s}' % j + gen_str(m[0][j]) for j in I if m[0][j])
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        return lhs + ' ' + ' '.join("g_{%d}" % i for i in redword)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, ['G',2])
            sage: dict(Y.algebra_generators())
            {'g1': g[1], 'g2': g[2], 'h1': h1, 'h2': h2}
        """
        one = self._W.one()
        zero = self._Q.zero()
        d = {}
        for i, al in self._Q.simple_roots().items():
            d['h%s' % i] = self.monomial((al, one))
        for i, g in self._W.simple_reflections().items():
            d['g%s' % i] = self.monomial((zero, g))
        return Family(sorted(d), d.__getitem__)

    @cached_method
    def one_basis(self):
        """
        Return the index of the basis element of `1`.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, ['D',6])
            sage: Y.one_basis()
            (0, 1)
        """
        return (self._Q.zero(), self._W.one())

    @cached_method
    def e(self, i=None):
        r"""
        Return the element(s) `e_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the element `e_i` or if ``None``,
          then the family of all idempotents `e_i`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['B',3])
            sage: Y.e(1)
            1/4 + 1/4*h1 + 1/4*h1^2 + 1/4*h1^3
            sage: Y.e(2)
            1/4 + 1/4*h2 + 1/4*h2^2 + 1/4*h2^3

        We verify that they are idempotents::

            sage: all(Y.e(i)^2 == Y.e(i) for i in Y.index_set())
            True

        Another example::

            sage: Y = algebras.YokonumaHecke(3, ['G',2])
            sage: e = Y.e()
            sage: all(e[i]^2 == e[i] for i in Y.index_set())
            True
        """
        if i is None:
            I = self._W.index_set()
            d = {i: self.e(i) for i in I}
            return Family(I, d.__getitem__)
        if i not in self._W.index_set():
            raise ValueError("invalid index")
        c = ~self.base_ring()(self._d)
        al = self._Q.simple_root(i)
        one = self._W.one()
        d = {(k*al, one): c for k in self._Q.base_ring()}
        return self._from_dict(d, remove_zeros=False)

    def h(self, i=None):
        r"""
        Return the generator(s) `h_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `h_i` or if ``None``,
          then the family of all generators `h_i`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(8, ['B',3])
            sage: Y.h(2)
            h2
            sage: Y.h()
            Finite family {1: h1, 2: h2, 3: h3}
        """
        G = self.algebra_generators()
        if i is None:
            I = self._W.index_set()
            d = {i: G['h%s' % i] for i in I}
            return Family(I, d.__getitem__)
        return G['h%s' % i]

    @cached_method
    def inverse_g(self, i=None):
        r"""
        Return the inverse of the generator(s) `g_i`.

        From the quadratic relation, we have

        .. MATH::

            g_i^{-1} = g_i + (q^{-1} - 1) e_i (1 + g_i).

        INPUT:

        - ``i`` -- (default: ``None``) the inverse generator `g_i^{-1}` or
          if ``None``, then the family of all inverse generators `g_i^{-1}`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(2, ['B',3])
            sage: [2*Y.inverse_g(i) for i in Y.index_set()]
            [(q^-1+1)*g[1] + (q^-1-1) + (q^-1-1)*h1*g[1] + (q^-1-1)*h1,
             (q^-1-1) + (q^-1+1)*g[2] + (q^-1-1)*h2 + (q^-1-1)*h2*g[2],
             (q^-1-1) + (q^-1+1)*g[3] + (q^-1-1)*h3 + (q^-1-1)*h3*g[3]]
            sage: all(Y.inverse_g(i) * Y.g(i) == Y.one() for i in range(1, 4))
            True
            sage: all(Y.g(i) * Y.inverse_g(i) == Y.one() for i in range(1, 4))
            True

            sage: Y = algebras.YokonumaHecke(3, ['G',2])
            sage: ginv = Y.inverse_g()
            sage: all(Y.g(i) * ginv[i] == Y.one() for i in Y.index_set())
            True
            sage: all(ginv[i] * Y.g(i) == Y.one() for i in Y.index_set())
            True
        """
        if i is None:
            I = self._W.index_set()
            d = {i: self.inverse_g(i) for i in I}
            return Family(I, d.__getitem__)
        if i not in self._W.index_set():
            raise ValueError("invalid index")
        return self.g(i) + (~self._q - 1) * self.e(i) * (self.one() + self.g(i))

    def product_on_basis(self, m1, m2):
        r"""
        Return the product of the basis elements indexed by ``m1`` and ``m2``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['C',3])
            sage: al = Y._Q.simple_root(1) + 2*Y._Q.simple_root(3)
            sage: w = Y._W.from_reduced_word([3,2,1,2]); w.length()
            4
            sage: Y.product_on_basis((Y._Q.zero(), w), (al, Y._W.one()))
            h2^3*h3*g[3,1,2,1]
            sage: Y.product_on_basis((al, w), (al, Y._W.one()))
            h1*h2^3*h3^3*g[3,1,2,1]
            sage: Y.product_on_basis((al, Y._W.one()), (al, w))
            h1^2*g[3,1,2,1]
            sage: 4 * Y.product_on_basis((al, w), (al, w))
            -(1-q)*h1*g[3,1,2,3,2,1] - (1-q)*h1*g[3,1,2,3,1,2,1]
             - (1-q)*h1*h2*h3*g[3,1,2,3,2,1] - (1-q)*h1*h2*h3*g[3,1,2,3,1,2,1]
             - (1-q)*h1*h2^2*h3^2*g[3,1,2,3,2,1] - (1-q)*h1*h2^2*h3^2*g[3,1,2,3,1,2,1]
             + (3+q)*h1*h2^3*h3^3*g[3,1,2,3,2,1] - (1-q)*h1*h2^3*h3^3*g[3,1,2,3,1,2,1]

        Check that we apply the permutation correctly on `h_i`::

            sage: Y = algebras.YokonumaHecke(4, ['B',3])
            sage: g1, g2, g3, h1, h2, h3 = Y.algebra_generators()
            sage: (g2 * g1) * h1
            h1^3*h2^3*g[2,1]
            sage: g2 * (g1 * h1)
            h1^3*h2^3*g[2,1]
        """
        h1, g1 = m1
        h2, g2 = m2
        # Commute g1 and t2, then multiply h1 and h2
        # ig1 = g1
        h = h1 + h2.weyl_action(g1)
        one = self._W.one()
        if g1 == one:
            return self.monomial((h, g2))
        ret = self.monomial((h, g1))
        for i in g2.reduced_word():
            ret = self.linear_combination((self._product_by_basis_gen(m, i), c)
                                          for m, c in ret)
        return ret

    def _product_by_basis_gen(self, m, i):
        r"""
        Return the product `t g_w g_i`.

        If the quadratic relation is `g_i^2 = 1 + (q - 1) (1 + g_i) e_i`,
        then we have

        .. MATH::

            g_w g_i = \begin{cases}
            g_{ws_i} & \text{if } \ell(ws_i) = \ell(w) + 1, \\
            g_{ws_i} + (q-1) (g_{ws_i} + g_w) e_i & \text{if }
              \ell(w s_i) = \ell(w) - 1.
            \end{cases}

        INPUT:

        - ``m`` -- a pair ``[h, w]``, where ``h`` encodes the monomial
          and ``w``  is an element of the Weyl group
        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, ['D',4])
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: 4 * Y._product_by_basis_gen(m, 1)  # not tested
            -(q^-1-q)*t2*t3^2*g[1] + 4*t1*t3^2 - (q^-1-q)*t1*t3^2*g[1]
             - (q^-1-q)*t1^2*t2^3*t3^2*g[1] - (q^-1-q)*t1^3*t2^2*t3^2*g[1]
        """
        h, w = m
        wi = w.apply_simple_reflection(i, side='right')
        if not w.has_descent(i, side='right'):
            return self.monomial((h, wi))

        q = self._q
        mon = self.monomial((h, wi))
        # TODO: Optimize this by computing an explicit expression
        #   for the commutation of w with ei.
        one = self.base_ring().one()
        binomial = self.element_class(self, {(h,wi): one, (h,w): one})
        return mon + (q-1) * binomial * self.e(i)

    class Element(CombinatorialFreeModule.Element):
        def __invert__(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: Y = algebras.YokonumaHecke(3, ['B',3])
                sage: all(g * ~g == Y.one() for g in Y.g())
                True
                sage: h = prod(Y.h()) * Y.h(2)
                sage: all(h*g * ~(h*g) == Y.one() for g in Y.g())
                True
            """
            if not self:
                raise ZeroDivisionError
            if len(self) != 1:
                raise NotImplementedError("inverse only implemented for basis elements (monomials in the generators)" % self)
            H = self.parent()
            t, w = self.support_of_term()
            c = ~self.coefficients()[0]
            telt = H.monomial((-t, H._W.one()))
            return c * H.prod(H.inverse_g(i) for i in reversed(w.reduced_word())) * telt
