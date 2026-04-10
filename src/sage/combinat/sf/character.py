# sage.doctest: needs sage.combinat sage.modules
r"""
Characters of the symmetric group as bases of the symmetric functions

Just as the Schur functions are the irreducible characters of `Gl_n`
and form a basis of the symmetric functions, the irreducible
symmetric group character basis are the irreducible characters of
of `S_n` when the group is realized as the permutation matrices.
"""

# ****************************************************************************
#       Copyright (C) 2015 Mike Zabrocki <zabrocki@mathstat.yorku.ca
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

from sage.arith.misc import binomial, divisors, moebius
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.combinat.sf.sfa import SymmetricFunctionAlgebra_generic as SFA_generic
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer


class Character_generic(SFA_generic):

    def _my_key(self, la):
        r"""
        A rank function for partitions.

        The leading term of a homogeneous expression will
        be the partition with the largest key value.

        This key value is `|\lambda|^2 + \lambda_0` and
        using the ``max`` function on a list of Partitions.

        Of course it is possible that this rank function
        is equal for some partitions, but the leading
        term should be any one of these partitions.

        INPUT:

        - ``la`` -- a partition

        OUTPUT: integer

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: ht._my_key(Partition([2,1,1]))
            18
            sage: ht._my_key(Partition([2,2]))
            18
            sage: ht._my_key(Partition([3,1]))
            19
            sage: ht._my_key(Partition([1,1,1,1]))
            17
        """
        if la:
            return la.size()**2 + la[0]
        return 0

    def _other_to_self(self, sexpr):
        r"""
        Convert an expression the target basis to the character-basis.

        We use triangularity to determine the expansion
        by subtracting off the leading term.  The target basis
        is specified by the method ``self._other``.

        INPUT:

        - ``sexpr`` -- an element of ``self._other`` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: h = Sym.h()
            sage: ht._other_to_self(h[2] + h([]))
            ht[] + ht[1] + ht[2]
            sage: st = SymmetricFunctions(QQ).st()
            sage: s = Sym.s()
            sage: st._other_to_self(s[1] + s([]))
            2*st[] + st[1]
            sage: 7 * st[[]] * st[[]]
            7*st[]
        """
        if sexpr == 0:
            return self(0)
        if list(sexpr.support()) == [[]]:
            return self._from_dict({self.one_basis(): sexpr.coefficient([])},
                                   remove_zeros=False)
        out = self.zero()
        while sexpr:
            mup = max(sexpr.support(), key=self._my_key)
            out += sexpr.coefficient(mup) * self(mup)
            sexpr -= sexpr.coefficient(mup) * self._self_to_other_on_basis(mup)
        return out

    def _b_power_k(self, k):
        r"""
        An expression involving Moebius inversion in the powersum generators.

        For a positive value of ``k``, this expression is

        .. MATH::

            \frac{1}{k} \sum_{d|k} \mu(d/k) p_d.

        INPUT:

        - ``k`` -- positive integer

        OUTPUT: an expression in the powersum basis of the symmetric functions

        EXAMPLES::

            sage: st = SymmetricFunctions(QQ).st()
            sage: st._b_power_k(1)
            p[1]
            sage: st._b_power_k(2)
            -1/2*p[1] + 1/2*p[2]
            sage: st._b_power_k(6)
            1/6*p[1] - 1/6*p[2] - 1/6*p[3] + 1/6*p[6]
        """
        if k == 1:
            return self._p([1])
        if k > 0:
            return ~k * self._p.linear_combination((self._p([d]),moebius(k//d))
                                    for d in divisors(k))


class InducedCharacterBases(Character_generic):
    r"""
    Character basis with Frobenius image of other times complete.
    """
    def __init__(self, Sym, other_basis, **kwds):
        r"""
        Initialize the basis and register coercions.

        INPUT:

        - ``Sym`` -- an instance of the symmetric function algebra
        - ``other_basis`` -- the target basis of the symmetric functions
        """
        SFA_generic.__init__(self, Sym, **kwds)
        self._other = other_basis
        self._p = Sym.powersum()
        self.module_morphism(self._self_to_power_on_basis,
                             codomain=Sym.powersum()).register_as_coercion()
        self.register_coercion(SetMorphism(Hom(self._other, self),
                                           self._other_to_self))

    def _b_bar_power_k_r(self, k, r):
        r"""
        An expression involving Moebius inversion in the powersum generators.

        For a positive value of ``k``, this expression is

        .. MATH::

            \sum_{j=0}^r (-1)^{r-j}k^j\binom{r,j}
            \left( \frac{1}{k} \sum_{d|k} \mu(d/k) p_d \right)_k.

        INPUT:

        - ``k``, ``r`` -- positive integers

        OUTPUT: an expression in the powersum basis of the symmetric functions

        EXAMPLES::

            sage: ht = SymmetricFunctions(QQ).ht()
            sage: ht._b_bar_power_k_r(1,1)
            p[1]
            sage: ht._b_bar_power_k_r(2,2)
            2*p[1] + p[1, 1] - 2*p[2] - 2*p[2, 1] + p[2, 2]
            sage: ht._b_bar_power_k_r(3,2)
            3*p[1] + p[1, 1] - 3*p[3] - 2*p[3, 1] + p[3, 3]
        """
        return k**r * self._p.prod( self._b_power_k(k)-j for j in range(r) )

    def _b_bar_power_gamma(self, gamma):
        r"""
        An expression involving Moebius inversion in the powersum generators.

        For a partition `\gamma = (1^{m_1}, 2^{m_2}, \ldots, r^{m_r})`,
        this expression is

        .. MATH::

            {\mathbf p}_{\ga} = \sum_{k \geq 1} {\mathbf p}_{k^{m_k}},

        where

        .. MATH::

            {\mathbf p}_{k^r} = \sum_{j=0}^r (-1)^{r-j}k^j\binom{r,j}
            \left( \frac{1}{k} \sum_{d|k} \mu(d/k) p_d \right)_k.

        INPUT:

        - ``gamma`` -- a partition

        OUTPUT: an expression in the powersum basis of the symmetric functions

        EXAMPLES::

            sage: ht = SymmetricFunctions(QQ).ht()
            sage: ht._b_bar_power_gamma(Partition([2,2,1]))
            2*p[1, 1] + p[1, 1, 1] - 2*p[2, 1] - 2*p[2, 1, 1] + p[2, 2, 1]
            sage: ht._b_bar_power_gamma(Partition([1,1,1]))
            2*p[1] - 3*p[1, 1] + p[1, 1, 1]
            sage: ht._b_bar_power_gamma(Partition([3,3,1]))
            3*p[1, 1] + p[1, 1, 1] - 3*p[3, 1] - 2*p[3, 1, 1] + p[3, 3, 1]
        """
        return self._p.prod(self._b_bar_power_k_r(Integer(k), Integer(r))
                            for k, r in gamma.to_exp_dict().items())

    def _self_to_power_on_basis(self, lam):
        r"""
        An expansion of the character basis in the powersum basis.

        The formula for the induced trivial character basis indexed by the
        partition ``lam`` is given by the formula

        .. MATH::

            \sum_{\gamma} \left\langle h_\lambda, p_\gamma \right\rangle
            \frac{{\overline {\mathbf p}}_\gamma}{z_\gamma},

        and for the irreducible rook character basis is given by

        .. MATH::

            \sum_{\gamma} \left\langle s_\lambda, p_\gamma \right\rangle
            \frac{{\overline {\mathbf p}}_\gamma}{z_\gamma},

        where if `\gamma = (1^{m_1}2^{m_2}\cdots \ell^{m_\ell})` then

        .. MATH::

            {\overline {\mathbf p}}_\gamma =
            \prod_{i \ge 1} i^{m_i} \prod_{n = 0}^{m_i-1} \left(
            \Big( \frac{1}{i} \sum_{d|i} \mu(i/d) p_d \Big) - n \right) .

        INPUT:

        - ``lam`` -- a partition

        OUTPUT: an expression in the power sum basis

        EXAMPLES::

            sage: ht = SymmetricFunctions(QQ).ht()
            sage: ht._self_to_power_on_basis([2,1])
            p[1] - 2*p[1, 1] + 1/2*p[1, 1, 1] + 1/2*p[2, 1]
            sage: ht._self_to_power_on_basis([1,1,1])
            2*p[1] - 3*p[1, 1] + p[1, 1, 1]

            sage: Sym = SymmetricFunctions(QQ)
            sage: xt = Sym.xt()
            sage: xt._self_to_power_on_basis([2,1])
            p[1] - p[1, 1] + 1/3*p[1, 1, 1] - 1/3*p[3]
            sage: xt._self_to_power_on_basis([1,1,1])
            1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3]
        """
        return self._p.sum( c*self._b_bar_power_gamma(ga)
                            for (ga, c) in self._p(self._other(lam)) )

    @cached_method
    def _self_to_other_on_basis(self, lam):
        r"""
        An expansion of the induced character basis in the ``other`` basis.

        Compute the ``self._other`` expansion by first computing it in the
        powersum basis and the coercing to the ``self._other`` basis.

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - an expression in the other basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: ht._self_to_other_on_basis(Partition([2,1]))
            h[1] - 2*h[1, 1] + h[2, 1]

            sage: Sym = SymmetricFunctions(QQ)
            sage: xt = Sym.xt()
            sage: xt._self_to_other_on_basis(Partition([2,1]))
            s[1] - s[1, 1] - s[2] + s[2, 1]

        TESTS::

            sage: h = SymmetricFunctions(QQ).h()
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: st = SymmetricFunctions(QQ).st()
            sage: all(ht(h(ht[la])) == ht[la] for i in range(5) for la in Partitions(i))
            True
            sage: all(h(ht(h[la])) == h[la] for i in range(5) for la in Partitions(i))
            True
            sage: all(st(h(st[la])) == st[la] for i in range(5) for la in Partitions(i))
            True
            sage: all(h(st(h[la])) == h[la] for i in range(5) for la in Partitions(i))
            True
        """
        return self._other(self._self_to_power_on_basis(lam))


class InducedTrivialCharacterBasis(InducedCharacterBases):
    r"""
    The induced trivial symmetric group character basis of
    the symmetric functions.

    This is a basis of the symmetric functions that has the
    property that ``self(la).character_to_frobenius_image(n)``
    is equal to ``h([n-sum(la)]+la)``.

    It has the property that the (outer) structure
    constants are the analogue of the stable Kronecker
    coefficients on the complete basis.

    This basis is introduced in [OZ2015]_.

    EXAMPLES::

        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.s()
        sage: h = Sym.h()
        sage: ht = SymmetricFunctions(QQ).ht()
        sage: st = SymmetricFunctions(QQ).st()
        sage: ht(s[2,1])
        ht[1, 1] + ht[2, 1] - ht[3]
        sage: s(ht[2,1])
        s[1] - 2*s[1, 1] - 2*s[2] + s[2, 1] + s[3]
        sage: ht(h[2,1])
        ht[1] + 2*ht[1, 1] + ht[2, 1]
        sage: h(ht[2,1])
        h[1] - 2*h[1, 1] + h[2, 1]
        sage: st(ht[2,1])
        st[] + 2*st[1] + st[1, 1] + 2*st[2] + st[2, 1] + st[3]
        sage: ht(st[2,1])
        ht[1] - ht[1, 1] + ht[2, 1] - ht[3]
        sage: ht[2]*ht[1,1]
        ht[1, 1] + 2*ht[1, 1, 1] + ht[2, 1, 1]
        sage: h[4,2].kronecker_product(h[4,1,1])
        h[2, 2, 1, 1] + 2*h[3, 1, 1, 1] + h[4, 1, 1]
        sage: s(st[2,1])
        3*s[1] - 2*s[1, 1] - 2*s[2] + s[2, 1]
        sage: st(s[2,1])
        st[] + 3*st[1] + 2*st[1, 1] + 2*st[2] + st[2, 1]
        sage: st[2]*st[1]
        st[1] + st[1, 1] + st[2] + st[2, 1] + st[3]
        sage: s[4,2].kronecker_product(s[5,1])
        s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2] + s[5, 1]
    """
    def __init__(self, Sym):
        r"""
        Initialize the basis and register coercions.

        The coercions are set up between the ``other_basis``.

        INPUT:

        - ``Sym`` -- an instance of the symmetric function algebra

        EXAMPLES::

            sage: ht = SymmetricFunctions(QQ).ht()
            sage: TestSuite(ht).run()
        """
        InducedCharacterBases.__init__(self, Sym, Sym.complete(),
            basis_name="induced trivial symmetric group character",
            prefix='ht', graded=False)


class RookIrreducibleCharacterBasis(InducedCharacterBases):
    r"""
    The irreducible character basis of the rook monoid as
    symmetric functions.

    It might also be called the induced irreducible character
    basis because it is the character of an irreducible
    `S_k` module induced to `S_n` where `n > k`.

    This basis appears implicitly in the paper by Assaf and Speyer [AS2020]_.
    The basis appears explicitly in [OZ2019]_.

    EXAMPLES::

        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.s()
        sage: h = Sym.h()
        sage: ht = Sym.ht()
        sage: st = Sym.st()
        sage: xt = Sym.xt(); xt
        Symmetric Functions over Rational Field in the irreducible rook monoid character basis
        sage: xt(s[2,1])
        xt[1, 1] + xt[2] + xt[2, 1]
        sage: s(xt[2,1])
        s[1] - s[1, 1] - s[2] + s[2, 1]
        sage: xt(h[2,1])
        xt[1] + 2*xt[1, 1] + 2*xt[2] + xt[2, 1] + xt[3]
        sage: h(xt[2,1])
        h[1] - h[1, 1] + h[2, 1] - h[3]
        sage: xt(ht[2,1])
        xt[2, 1] + xt[3]
        sage: st(xt[2,1])
        st[1] + st[1, 1] + st[2] + st[2, 1]
        sage: xt(st[2,1])
        xt[1] - xt[1, 1] - xt[2] + xt[2, 1]
        sage: xt[2]*xt[1,1]
        xt[1, 1] + xt[1, 1, 1] + 2*xt[2, 1] + xt[2, 1, 1] + xt[3] + xt[3, 1]
        sage: xt(s[2,1])
        xt[1, 1] + xt[2] + xt[2, 1]
        sage: xt[2]*xt[1]
        xt[1, 1] + xt[2] + xt[2, 1] + xt[3]

    This is a basis of the symmetric functions that has the
    property that ``self(la).character_to_frobenius_image(n)``
    is equal to ``s([n-sum(la)]) * s(la)``.::

        sage: s(xt[2,1].character_to_frobenius_image(3))
        s[2, 1]
        sage: s(xt[2,1].character_to_frobenius_image(9)) == s[6] * s[2,1]
        True
    """
    def __init__(self, Sym):
        r"""
        Initialize the basis and register coercions.

        The coercions are set up between the ``other_basis``.
        This code is almost exactly the same as the
        induced trivial character basis, but the other basis
        is the complete symmetric functions instead of the
        Schur basis.

        INPUT:

        - ``Sym`` -- an instance of the symmetric function algebra

        EXAMPLES::

            sage: xt = SymmetricFunctions(QQ).xt()
            sage: TestSuite(xt).run()
        """
        InducedCharacterBases.__init__(self, Sym, Sym.Schur(),
            basis_name="irreducible rook monoid character",
            prefix='xt', graded=False)


class IrreducibleCharacterBasis(Character_generic):
    r"""
    The irreducible symmetric group character basis of
    the symmetric functions.

    This is a basis of the symmetric functions that has the
    property that ``self(la).character_to_frobenius_image(n)``
    is equal to ``s([n-sum(la)]+la)``.

    It should also have the property that the (outer) structure
    constants are the analogue of the stable Kronecker
    coefficients on the Schur basis.

    This basis is introduced in [OZ2015]_.

    EXAMPLES::

        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.s()
        sage: h = Sym.h()
        sage: ht = SymmetricFunctions(QQ).ht()
        sage: st = SymmetricFunctions(QQ).st()
        sage: st(ht[2,1])
        st[] + 2*st[1] + st[1, 1] + 2*st[2] + st[2, 1] + st[3]
        sage: ht(st[2,1])
        ht[1] - ht[1, 1] + ht[2, 1] - ht[3]
        sage: s(st[2,1])
        3*s[1] - 2*s[1, 1] - 2*s[2] + s[2, 1]
        sage: st(s[2,1])
        st[] + 3*st[1] + 2*st[1, 1] + 2*st[2] + st[2, 1]
        sage: st[2]*st[1]
        st[1] + st[1, 1] + st[2] + st[2, 1] + st[3]
        sage: s[4,2].kronecker_product(s[5,1])
        s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2] + s[5, 1]
        sage: st[1,1,1].counit()
        -1
        sage: all(sum(c*st(la)*st(mu).antipode() for
        ....:    ((la,mu),c) in st(ga).coproduct())==st(st(ga).counit())
        ....:    for ga in Partitions(3))
        True

    TESTS::

        sage: TestSuite(st).run()
    """

    def __init__(self, Sym):
        r"""
        Initialize the basis and register coercions.

        The coercions are set up between the ``other_basis``

        INPUT:

        - ``Sym`` -- an instance of the symmetric function algebra
        - ``pfix`` -- a prefix to use for the basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht(); ht
            Symmetric Functions over Rational Field in the induced trivial
             symmetric group character basis
            sage: st = SymmetricFunctions(QQ).st(); st
            Symmetric Functions over Rational Field in the irreducible
             symmetric group character basis
        """
        SFA_generic.__init__(self, Sym,
                             basis_name="irreducible symmetric group character",
                             prefix='st', graded=False)
        self._other = Sym.Schur()
        self._p = Sym.powersum()

        self.module_morphism(self._self_to_power_on_basis,
                             codomain=Sym.powersum()).register_as_coercion()
        self.register_coercion(SetMorphism(Hom(self._other, self),
                                           self._other_to_self))

    def _b_power_k_r(self, k, r):
        r"""
        An expression involving Moebius inversion in the powersum generators.

        For a positive value of ``k``, this expression is

        .. MATH::

            \sum_{j=0}^r (-1)^{r-j}k^j\binom{r,j} \left(
            \frac{1}{k} \sum_{d|k} \mu(d/k) p_d \right)_k.

        INPUT:

        - ``k``, ``r`` -- positive integers

        OUTPUT: an expression in the powersum basis of the symmetric functions

        EXAMPLES::

            sage: st = SymmetricFunctions(QQ).st()
            sage: st._b_power_k_r(1,1)
            -p[] + p[1]
            sage: st._b_power_k_r(2,2)
            p[] + 4*p[1] + p[1, 1] - 4*p[2] - 2*p[2, 1] + p[2, 2]
            sage: st._b_power_k_r(3,2)
            p[] + 5*p[1] + p[1, 1] - 5*p[3] - 2*p[3, 1] + p[3, 3]
        """
        p = self._p
        return p.sum( (-1)**(r-j) * k**j * binomial(r,j)
                      * p.prod(self._b_power_k(k) - i*p.one() for i in range(j))
                      for j in range(r+1) )

    def _b_power_gamma(self, gamma):
        r"""
        An expression involving Moebius inversion in the powersum generators.

        For a partition `\gamma = (1^{m_1}, 2^{m_2}, \ldots, r^{m_r})`,
        this expression is

        .. MATH::

            {\mathbf p}_{\ga} = \sum_{k \geq 1} {\mathbf p}_{k^{m_k}},

        where

        .. MATH::

            {\mathbf p}_{k^r} = \sum_{j=0}^r (-1)^{r-j}k^j\binom{r,j}
            \left( \frac{1}{k} \sum_{d|k} \mu(d/k) p_d \right)_k.

        INPUT:

        - ``gamma`` -- a partition

        OUTPUT: an expression in the powersum basis of the symmetric functions

        EXAMPLES::

            sage: st = SymmetricFunctions(QQ).st()
            sage: st._b_power_gamma(Partition([2,2]))
            p[] + 4*p[1] + p[1, 1] - 4*p[2] - 2*p[2, 1] + p[2, 2]
            sage: st._b_power_gamma(Partition([1,1,1]))
            -p[] + 8*p[1] - 6*p[1, 1] + p[1, 1, 1]
            sage: st._b_power_gamma(Partition([3,1]))
            p[] - p[1, 1] - p[3] + p[3, 1]
        """
        return self._p.prod(self._b_power_k_r(Integer(k), Integer(r))
                            for k, r in gamma.to_exp_dict().items())

    def _self_to_power_on_basis(self, lam):
        r"""
        An expansion of the irreducible character in the powersum basis.

        The formula for the irreducible character basis indexed by the
        partition ``lam`` is given by the formula

        .. MATH::

            \sum_{\gamma} \chi^{\lambda}(\gamma)
            \frac{{\mathbf p}_\gamma}{z_\gamma},

        where `\chi^{\lambda}(\gamma)` is the irreducible character
        indexed by the partition `\lambda` and evaluated at an element
        of cycle structure `\gamma` and `{\mathbf p}_\gamma` is the
        power sum expression calculated in the method
        :meth:`_b_power_gamma`.

        INPUT:

        - ``lam`` -- a partition

        OUTPUT: an expression in the power sum basis

        EXAMPLES::

            sage: st = SymmetricFunctions(QQ).st()
            sage: st._self_to_power_on_basis([2,1])
            3*p[1] - 2*p[1, 1] + 1/3*p[1, 1, 1] - 1/3*p[3]
            sage: st._self_to_power_on_basis([1,1])
            p[] - p[1] + 1/2*p[1, 1] - 1/2*p[2]
        """
        return self._p.sum( c*self._b_power_gamma(ga)
                            for (ga, c) in self._p(self._other(lam)) )

    @cached_method
    def _self_to_other_on_basis(self, lam):
        r"""
        An expansion of the irreducible character basis in the Schur basis.

        Compute the Schur expansion by first computing it in the
        powersum basis and the coercing to the Schur basis.

        INPUT:

        - ``lam`` -- a partition

        OUTPUT: an expression in the Schur basis

        EXAMPLES::

            sage: st = SymmetricFunctions(QQ).st()
            sage: st._self_to_other_on_basis(Partition([1,1]))
            s[] - s[1] + s[1, 1]
            sage: st._self_to_other_on_basis(Partition([2,1]))
            3*s[1] - 2*s[1, 1] - 2*s[2] + s[2, 1]
        """
        return self._other(self._self_to_power_on_basis(lam))
