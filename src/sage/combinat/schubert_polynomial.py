# sage.doctest: needs sage.combinat sage.modules
r"""
Schubert Polynomials


See :wikipedia:`Schubert_polynomial` and
`SymmetricFunctions.com <https://www.symmetricfunctions.com/schubert.htm#schubert>`_.
Schubert polynomials are representatives of cohomology classes in flag varieties.
In `n` variables, they are indexed by permutations `w \in S_n`. They also form
a basis for the coinvariant ring of the `S_n` action on
`\ZZ[x_1, x_2, \ldots, x_n]`.

EXAMPLES::

    sage: X = SchubertPolynomialRing(ZZ)
    sage: w = [1,2,5,4,3];  # a list representing an element of `S_5`
    sage: X(w)
    X[1, 2, 5, 4, 3]

This can be expanded in terms of polynomial variables::

    sage: X(w).expand()
    x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2
     + x0*x2^2 + x1*x2^2 + x0^2*x3 + x0*x1*x3 + x1^2*x3
     + x0*x2*x3 + x1*x2*x3 + x2^2*x3

We can also convert back from polynomial variables. For example,
the longest permutation is a single term. In `S_5`, this is the
element (in one line notation) `w_0 = 54321`::

    sage: w0 = [5,4,3,2,1]
    sage: R.<x0, x1, x2, x3, x4> = PolynomialRing(ZZ)
    sage: Sw0 = X(x0^4*x1^3*x2^2*x3);  Sw0
    X[5, 4, 3, 2, 1]

The polynomials also have the property that if the indexing permutation `w` is
multiplied by a simple transposition `s_i = (i, i+1)` such that the length of
`w` is more than the length of `ws_i`, then the Schubert
polynomial of the permutation `ws_i` is computed by applying the divided
difference operator :meth:`~SchubertPolynomial_class.divided_difference` to
the polynomial indexed by `w`. For example, applying the divided difference
operator `\partial_2` to the Schubert polynomial `\mathfrak{S}_{w_0}`::

    sage: Sw0.divided_difference(2)
    X[5, 3, 4, 2, 1]

We can also check the properties listed in :wikipedia:`Schubert_polynomial`::

    sage: X([1,2,3,4,5])  # the identity in one-line notation
    X[1]
    sage: X([1,3,2,4,5]).expand()  # the transposition swapping 2 and 3
    x0 + x1
    sage: X([2,4,5,3,1]).expand()
    x0^2*x1^2*x2*x3 + x0^2*x1*x2^2*x3 + x0*x1^2*x2^2*x3

    sage: w = [4,5,1,2,3]
    sage: s = SymmetricFunctions(QQ).schur()
    sage: s[3,3].expand(2)
    x0^3*x1^3
    sage: X(w).expand()
    x0^3*x1^3
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.key_polynomial import KeyPolynomial
from sage.combinat.permutation import Permutations, Permutation
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

lazy_import('sage.libs.symmetrica', 'all', as_='symmetrica')


def SchubertPolynomialRing(R):
    """
    Return the Schubert polynomial ring over ``R`` on the X basis.

    This is the basis made of the Schubert polynomials.

    EXAMPLES::

        sage: X = SchubertPolynomialRing(ZZ); X
        Schubert polynomial ring with X basis over Integer Ring
        sage: TestSuite(X).run()
        sage: X(1)
        X[1]
        sage: X([1,2,3])*X([2,1,3])
        X[2, 1]
        sage: X([2,1,3])*X([2,1,3])
        X[3, 1, 2]
        sage: X([2,1,3])+X([3,1,2,4])
        X[2, 1] + X[3, 1, 2]
        sage: a = X([2,1,3])+X([3,1,2,4])
        sage: a^2
        X[3, 1, 2] + 2*X[4, 1, 2, 3] + X[5, 1, 2, 3, 4]
    """
    return SchubertPolynomialRing_xbasis(R)


class SchubertPolynomial_class(CombinatorialFreeModule.Element):
    def expand(self):
        """
        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: X([2,1,3]).expand()
            x0
            sage: [X(p).expand() for p in Permutations(3)]
            [1, x0 + x1, x0, x0*x1, x0^2, x0^2*x1]

        TESTS:

        Calling .expand() should always return an element of an
        MPolynomialRing::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: f = X([1]); f
            X[1]
            sage: type(f.expand())
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f.expand()
            1
            sage: f = X([1,2])
            sage: type(f.expand())
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f = X([1,3,2,4])
            sage: type(f.expand())
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>

        Now we check for correct handling of the empty
        permutation (:issue:`23443`)::

            sage: X([1]).expand() * X([2,1]).expand()
            x0
        """
        p = symmetrica.t_SCHUBERT_POLYNOM(self)
        if not isinstance(p, MPolynomial):
            R = PolynomialRing(self.parent().base_ring(), 1, 'x0')
            p = R(p)
        return p

    def divided_difference(self, i, algorithm='sage'):
        r"""
        Return the ``i``-th divided difference operator, applied to ``self``.

        Here, ``i`` can be either a permutation or a positive integer.

        INPUT:

        - ``i`` -- permutation or positive integer

        - ``algorithm`` -- (default: ``'sage'``) either ``'sage'``
          or ``'symmetrica'``; this determines which software is
          called for the computation

        OUTPUT:

        The result of applying the ``i``-th divided difference
        operator to ``self``.

        If `i` is a positive integer, then the `i`-th divided
        difference operator `\delta_i` is the linear operator sending
        each polynomial `f = f(x_1, x_2, \ldots, x_n)` (in
        `n \geq i+1` variables) to the polynomial

        .. MATH::

            \frac{f - f_i}{x_i - x_{i+1}}, \qquad \text{ where }
            f_i = f(x_1, x_2, ..., x_{i-1}, x_{i+1}, x_i,
            x_{i+1}, ..., x_n) .

        If `\sigma` is a permutation in the `n`-th symmetric group,
        then the `\sigma`-th divided difference operator `\delta_\sigma`
        is the composition
        `\delta_{i_1} \delta_{i_2} \cdots \delta_{i_k}`, where
        `\sigma = s_{i_1} \circ s_{i_2} \circ \cdots \circ s_{i_k}` is
        any reduced expression for `\sigma` (the precise choice of
        reduced expression is immaterial).

        .. NOTE::

            The :meth:`expand` method results in a polynomial
            in `n` variables named ``x0, x1, ..., x(n-1)`` rather than
            `x_1, x_2, \ldots, x_n`.
            The variable named ``xi`` corresponds to `x_{i+1}`.
            Thus, ``self.divided_difference(i)`` involves the variables
            ``x(i-1)`` and ``xi`` getting switched (in the numerator).

        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: a = X([3,2,1])
            sage: a.divided_difference(1)
            X[2, 3, 1]
            sage: a.divided_difference([3,2,1])
            X[1]
            sage: a.divided_difference(5)
            0

        Any divided difference of `0` is `0`::

            sage: X.zero().divided_difference(2)
            0

        This is compatible when a permutation is given as input::

            sage: a = X([3,2,4,1])
            sage: a.divided_difference([2,3,1])
            0
            sage: a.divided_difference(1).divided_difference(2)
            0

        ::

            sage: a = X([4,3,2,1])
            sage: a.divided_difference([2,3,1])
            X[3, 2, 4, 1]
            sage: a.divided_difference(1).divided_difference(2)
            X[3, 2, 4, 1]
            sage: a.divided_difference([4,1,3,2])
            X[1, 4, 2, 3]
            sage: b = X([4, 1, 3, 2])
            sage: b.divided_difference(1).divided_difference(2)
            X[1, 3, 4, 2]
            sage: b.divided_difference(1).divided_difference(2).divided_difference(3)
            X[1, 3, 2]
            sage: b.divided_difference(1).divided_difference(2).divided_difference(3).divided_difference(2)
            X[1]
            sage: b.divided_difference(1).divided_difference(2).divided_difference(3).divided_difference(3)
            0
            sage: b.divided_difference(1).divided_difference(2).divided_difference(1)
            0

        TESTS:

        Check that :issue:`23403` is fixed::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: a = X([3,2,4,1])
            sage: a.divided_difference(2)
            0
            sage: a.divided_difference([3,2,1])
            0
            sage: a.divided_difference(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot apply \delta_{0} to a (= X[3, 2, 4, 1])
        """
        if not self:  # if self is 0
            return self
        Perms = Permutations()
        if i in ZZ:
            if algorithm == "sage":
                if i <= 0:
                    raise ValueError(r"cannot apply \delta_{%s} to a (= %s)" % (i, self))
                # The operator `\delta_i` sends the Schubert
                # polynomial `X_\pi` (where `\pi` is a finitely supported
                # permutation of `\{1, 2, 3, \ldots\}`) to:
                # - the Schubert polynomial X_\sigma`, where `\sigma` is
                #   obtained from `\pi` by switching the values at `i` and `i+1`,
                #   if `i` is a descent of `\pi` (that is, `\pi(i) > \pi(i+1)`);
                # - `0` otherwise.
                # Notice that distinct `\pi`s lead to distinct `\sigma`s,
                # so we can use `_from_dict` here.
                res_dict = {}
                for pi, coeff in self:
                    pi = pi[:]
                    n = len(pi)
                    if n <= i:
                        continue
                    if pi[i - 1] < pi[i]:
                        continue
                    pi[i - 1], pi[i] = pi[i], pi[i - 1]
                    pi = Perms(pi).remove_extra_fixed_points()
                    res_dict[pi] = coeff
                return self.parent()._from_dict(res_dict)
            else:  # if algorithm == "symmetrica":
                return symmetrica.divdiff_schubert(i, self)
        elif i in Perms:
            if algorithm == "sage":
                i = Permutation(i)
                redw = i.reduced_word()
                res_dict = {}
                for pi, coeff in self:
                    next_pi = False
                    pi = pi[:]
                    n = len(pi)
                    for j in redw:
                        if n <= j:
                            next_pi = True
                            break
                        if pi[j - 1] < pi[j]:
                            next_pi = True
                            break
                        pi[j - 1], pi[j] = pi[j], pi[j - 1]
                    if next_pi:
                        continue
                    pi = Perms(pi).remove_extra_fixed_points()
                    res_dict[pi] = coeff
                return self.parent()._from_dict(res_dict)
            else:  # if algorithm == "symmetrica":
                return symmetrica.divdiff_perm_schubert(i, self)
        else:
            raise TypeError("i must either be an integer or permutation")

    def scalar_product(self, x):
        """
        Return the standard scalar product of ``self`` and ``x``.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: a = X([3,2,4,1])
            sage: a.scalar_product(a)
            0
            sage: b = X([4,3,2,1])
            sage: b.scalar_product(a)
            X[1, 3, 4, 6, 2, 5]
            sage: Permutation([1, 3, 4, 6, 2, 5, 7]).to_lehmer_code()
            [0, 1, 1, 2, 0, 0, 0]
            sage: s = SymmetricFunctions(ZZ).schur()
            sage: c = s([2,1,1])
            sage: b.scalar_product(a).expand()
            x0^2*x1*x2 + x0*x1^2*x2 + x0*x1*x2^2 + x0^2*x1*x3 + x0*x1^2*x3
             + x0^2*x2*x3 + 3*x0*x1*x2*x3 + x1^2*x2*x3 + x0*x2^2*x3 + x1*x2^2*x3
             + x0*x1*x3^2 + x0*x2*x3^2 + x1*x2*x3^2
            sage: c.expand(4)
            x0^2*x1*x2 + x0*x1^2*x2 + x0*x1*x2^2 + x0^2*x1*x3 + x0*x1^2*x3
             + x0^2*x2*x3 + 3*x0*x1*x2*x3 + x1^2*x2*x3 + x0*x2^2*x3 + x1*x2^2*x3
             + x0*x1*x3^2 + x0*x2*x3^2 + x1*x2*x3^2
        """
        if isinstance(x, SchubertPolynomial_class):
            return symmetrica.scalarproduct_schubert(self, x)
        else:
            raise TypeError("x must be a Schubert polynomial")

    def multiply_variable(self, i):
        """
        Return the Schubert polynomial obtained by multiplying ``self``
        by the variable `x_i`.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: a = X([3,2,4,1])
            sage: a.multiply_variable(0)
            X[4, 2, 3, 1]
            sage: a.multiply_variable(1)
            X[3, 4, 2, 1]
            sage: a.multiply_variable(2)
            X[3, 2, 5, 1, 4] - X[3, 4, 2, 1] - X[4, 2, 3, 1]
            sage: a.multiply_variable(3)
            X[3, 2, 4, 5, 1]
        """
        if isinstance(i, Integer):
            return symmetrica.mult_schubert_variable(self, i)
        else:
            raise TypeError("i must be an integer")


class SchubertPolynomialRing_xbasis(CombinatorialFreeModule):

    Element = SchubertPolynomial_class

    def __init__(self, R):
        """
        EXAMPLES::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X == loads(dumps(X))
            True
        """
        self._name = "Schubert polynomial ring with X basis"
        self._repr_option_bracket = False
        CombinatorialFreeModule.__init__(self, R, Permutations(),
                                         category=GradedAlgebrasWithBasis(R),
                                         prefix='X')

    @cached_method
    def one_basis(self):
        """
        Return the index of the unit of this algebra.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X.one()  # indirect doctest
            X[1]
        """
        return self._indices([1])

    def _element_constructor_(self, x):
        """
        Coerce x into ``self``.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X._element_constructor_([2,1,3])
            X[2, 1]
            sage: X._element_constructor_(Permutation([2,1,3]))
            X[2, 1]

            sage: R.<x1, x2, x3> = QQ[]
            sage: X(x1^2*x2)
            X[3, 2, 1]

            sage: S.<x> = InfinitePolynomialRing(QQ)
            sage: X(x[0]^2*x[1])
            X[3, 2, 1]
            sage: X(x[0]*x[1]^2*x[2]^2*x[3] + x[0]^2*x[1]^2*x[2]*x[3] + x[0]^2*x[1]*x[2]^2*x[3])
            X[2, 4, 5, 3, 1]

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: X(k([3,2,1]))
            X[4, 3, 2, 1]

        TESTS:

        We check that :issue:`12924` is fixed::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X._element_constructor_([1,2,1])
            Traceback (most recent call last):
            ...
            ValueError: the input [1, 2, 1] is not a valid permutation

        Now we check for correct handling of the empty
        permutation (:issue:`23443`)::

            sage: X([])
            X[1]

        Check the round trip from key polynomials::

            sage: k = KeyPolynomials(ZZ)
            sage: X = SchubertPolynomialRing(ZZ)
            sage: it = iter(Permutations())
            sage: for _ in range(50):
            ....:     P = next(it)
            ....:     assert X(k(X(P))) == X(P), P
        """
        if isinstance(x, list):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            if x not in Permutations():
                raise ValueError(f"the input {x} is not a valid permutation")
            perm = Permutation(x).remove_extra_fixed_points()
            return self._from_dict({perm: self.base_ring().one()})
        elif isinstance(x, Permutation):
            perm = x.remove_extra_fixed_points()
            return self._from_dict({perm: self.base_ring().one()})
        elif isinstance(x, MPolynomial):
            return symmetrica.t_POLYNOM_SCHUBERT(x)
        elif isinstance(x, InfinitePolynomial):
            R = x.polynomial().parent()
            # massage the term order to be what symmetrica expects
            S = PolynomialRing(R.base_ring(),
                               names=list(map(repr, reversed(R.gens()))))
            return symmetrica.t_POLYNOM_SCHUBERT(S(x.polynomial()))
        elif isinstance(x, KeyPolynomial):
            return self(x.expand())
        else:
            raise TypeError

    def some_elements(self):
        """
        Return some elements.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X.some_elements()
            [X[1], X[1] + 2*X[2, 1], -X[3, 2, 1] + X[4, 2, 1, 3]]
        """
        return [self.one(), self([1]) + 2 * self([2, 1]),
                self([4, 2, 1, 3]) - self([3, 2, 1])]

    def product_on_basis(self, left, right):
        """
        EXAMPLES::

            sage: p1 = Permutation([3,2,1])
            sage: p2 = Permutation([2,1,3])
            sage: X = SchubertPolynomialRing(QQ)
            sage: X.product_on_basis(p1,p2)
            X[4, 2, 1, 3]
        """
        return symmetrica.mult_schubert_schubert(left, right)
