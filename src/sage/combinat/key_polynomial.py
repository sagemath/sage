r"""
Key polynomials

Key polynomials (a.k.a. Demazure characters) are defined by applying an
operator `\pi_\sigma` to a monomial corresponding to an integer partition `\mu
\vdash n`. For a full definition, see
`SymmetricFunctions.com<https://www.symmetricfunctions.com/key.htm>`_. Key
polynomials are indexed by weak compositions with no trailing zeros, and
`\sigma` is the permutation of shortest length which sorts the indexing
composition into a partition.

.. SEEALSO::

    For Demazure characters in other types, see
    - :meth:`sage.combinat.root_system.weyl_characters.WeylCharacterRing.demazure_character`
    - :meth:`sage.categories.classical_crystals.ClassicalCrystals.ParentMethods.demazure_character`

EXAMPLES:

    Key polynomials are a basis, indexed by (weak) compositions with no
    trailing zeros, for infinite variable polynomial rings over a
    characteristic-`0` field::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
        sage: k = KeyPolynomialBasis(QQ)
        sage: k([3,0,1,2])
        k[3, 0, 1, 2]
        sage: R = k.polynomial_ring(); R
        Infinite polynomial ring in z over Rational Field

    We can expand them in the standard monomial basis::

        sage: k([3,0,1,2]).expand()
        z_3^2*z_2*z_0^3 + z_3^2*z_1*z_0^3 + z_3*z_2^2*z_0^3
         + 2*z_3*z_2*z_1*z_0^3 + z_3*z_1^2*z_0^3 + z_2^2*z_1*z_0^3
         + z_2*z_1^2*z_0^3

        sage: k([0,0,2]).expand()
        z_2^2 + z_2*z_1 + z_2*z_0 + z_1^2 + z_1*z_0 + z_0^2

    If we have a polynomial, we can express it in the key basis::

        sage: z, = R.gens()
        sage: k.from_polynomial(z[2]^2*z[1]*z[0])
        k[1, 1, 2] - k[1, 2, 1]

        sage: f = z[3]^2*z[2]*z[0]^3 + z[3]^2*z[1]*z[0]^3 + z[3]*z[2]^2*z[0]^3 + \
        ....: 2*z[3]*z[2]*z[1]*z[0]^3 + z[3]*z[1]^2*z[0]^3 + z[2]^2*z[1]*z[0]^3 + \
        ....: z[2]*z[1]^2*z[0]^3
        sage: k.from_polynomial(f)
        k[3, 0, 1, 2]

    Since the ring of key polynomials may be regarded as a different choice of
    basis for a polynomial ring, it forms an algebra, so we have
    multiplication::

        sage: k([10,5,2])*k([1,1,1])
        k[11, 6, 3]

    When the sorting permutation is the longest element, the key polynomial
    agrees with the Schur polynomial::

        sage: s = SymmetricFunctions(QQ).schur()
        sage: k([1,2,3]).expand()
        z_2^3*z_1^2*z_0 + z_2^3*z_1*z_0^2 + z_2^2*z_1^3*z_0
         + 2*z_2^2*z_1^2*z_0^2 + z_2^2*z_1*z_0^3 + z_2*z_1^3*z_0^2
         + z_2*z_1^2*z_0^3
        sage: s[3,2,1].expand(3)
        x0^3*x1^2*x2 + x0^2*x1^3*x2 + x0^3*x1*x2^2 + 2*x0^2*x1^2*x2^2
         + x0*x1^3*x2^2 + x0^2*x1*x2^3 + x0*x1^2*x2^3

    The polynomial expansions can be computed using crystals and expressed in
    terms of the key basis::

        sage: T = crystals.Tableaux(['A',3],shape=[2,1])
        sage: f = T.demazure_character([3,2,1])
        sage: k.from_polynomial(f)
        k[1, 0, 0, 2]

AUTHORS:

- Trevor K. Karn (2022-08-17): initial version

"""

# ****************************************************************************
#       Copyright (C) 2022 Trevor K. Karn <karnx018 (at) umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.composition import Composition
from sage.combinat.integer_vector import IntegerVector, IntegerVectors
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial_sparse
from sage.rings.rational_field import QQ


class KeyPolynomial(CombinatorialFreeModule.Element):

    def expand(self):
        r"""
        Return ``self`` written in the monomial basis.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: f = k([4,3,2,1])
            sage: f.expand()
            z_3*z_2^2*z_1^3*z_0^4
            sage: f = k([1,2,3])
            sage: f.expand()
            z_2^3*z_1^2*z_0 + z_2^3*z_1*z_0^2 + z_2^2*z_1^3*z_0 + 2*z_2^2*z_1^2*z_0^2 + z_2^2*z_1*z_0^3 + z_2*z_1^3*z_0^2 + z_2*z_1^2*z_0^3
        """
        R = self.parent().polynomial_ring()
        out = R.zero()
        z, = R.gens()

        for m, c in self.monomial_coefficients().items():
            # get the partition itself
            mu = sorted(m, reverse=True)

            # create the monomial to apply
            monom = R.prod(z[i] ** mi for i, mi in enumerate(mu) if mi)
            
            # find the permutation sorting m into a mu
            w = _sorting_word(m)

            out += c * _pi(self.parent(), w, monom)

        return out

    def pi(self, w):
        r"""
        Apply the operator `\pi_w` to ``self``.

        ``w`` may be either a ``Permutation`` or a list of indices of simple
        transpositions. In the case that it is a list of indices of simple
        transpositions, the convention follows that of the output of
        :meth:`sage.combinat.permutations.Permutation.reduced_word`
        namely that the transpositions are applied from left to right
        rather than from right to left.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k([3,2,1]).pi(1)
            k[3, 1, 2]

            sage: k([3,2,1]).pi([1,0])
            k[1, 3, 2]

            sage: k([3,2,1]).pi(Permutation([3,2,1]))
            k[1, 2, 3]
        """
        P = self.parent()
        f = self.expand()

        if isinstance(w, Permutation):
            w = w.reduced_word()
            return P.from_polynomial(_pi(P, [i-1 for i in w], f))

        # This can be done without punting to the polynomial ring I think.

        return P.from_polynomial(_pi(P, w, f))

    def divided_difference(self, i):
        r"""
        Apply the divided difference operator `\partial_i` to ``self``.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k([3,2,1]).divided_difference(1)
            k[3, 1, 1]
        """
        P = self.parent()
        f = self.expand()
        return P.from_polynomial(_divided_difference(P, i, f))

class KeyPolynomialBasis(CombinatorialFreeModule):

    Element = KeyPolynomial

    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k == loads(dumps(k))
            True
            sage: TestSuite(k).run()
        """
        self._name = "Ring of key polynomials"
        self._repr_option_bracket = False
        CombinatorialFreeModule.__init__(self, R, _IntegerVectors,
                                         category=GradedAlgebrasWithBasis(R),
                                         prefix='k')

    def _element_constructor_(self, alpha):
        r"""
        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k([10,0,1,2])
            k[10, 0, 1, 2]
            sage: IV = IntegerVectors()
            sage: k(IV([9,0,3,2,4]))
            k[9, 0, 3, 2, 4]
            sage: k(Composition([9,1,2,3]))
            k[9, 1, 2, 3]

            sage: R.<z> = InfinitePolynomialRing(QQ)
            sage: k(z[4]*z[3]*z[2]*z[1]^2*z[0]^3)
            k[3, 2, 1, 1, 1]

            sage: z0, z1, z2 = var('z_0, z_1, z_2')
            sage: f = z2^4*z1^5*z0^9
            sage: k(R(f))
            k[9, 5, 4]
        """
        if isinstance(alpha, (list, IntegerVector, Composition)):
            alpha = _IntegerVectors(alpha).trim()
            return self._from_dict({alpha: self.base_ring().one()})
        if isinstance(alpha, InfinitePolynomial_sparse):
            return self.from_polynomial(alpha)
        raise ValueError("alpha can not be interpreted as a weak composition or polynomial")

    def one_basis(self):
        r"""
        Return the basis element indexing the identity.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k.one_basis()
            []
        """
        return _IntegerVectors([])

    def polynomial_ring(self):
        r"""
        Return the polynomial ring associated to ``self``.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k.polynomial_ring()
            Infinite polynomial ring in z over Rational Field
        """
        return InfinitePolynomialRing(self.base_ring(), 'z')

    def from_polynomial(self, f):
        r"""
        Expand a polynomial in terms of the key basis.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: K = KeyPolynomialBasis(QQ)
            sage: R = K.polynomial_ring(); z, = R.gens()
            sage: p = z[0]^4*z[1]^2*z[2]*z[3] + z[0]^4*z[1]*z[2]^2*z[3]
            sage: K.from_polynomial(p)
            k[4, 1, 2, 1]

            sage: all(K(c) == K.from_polynomial(K(c).expand()) for c in Compositions(5))
            True

            sage: T = crystals.Tableaux(['A', 4], shape=[4,2,1,1])
            sage: K.from_polynomial(T.demazure_character([2]))
            k[4, 1, 2, 1]

        """
        if f not in self.polynomial_ring():
            try:  # to accept elements of SymbolicRing
                from sage.calculus.var import var
                f = f.substitute(list(d == var(f'z_{i}')
                                 for i, d in enumerate(f.variables())))
                f = self.polynomial_ring()(f)
            except AttributeError:
                raise ValueError(f"f must be an element of {self.polynomial_ring()}")

        out = self.zero()
        counter = 0

        while f.monomials():
            M = f.monomials()[0]
            c = f.monomial_coefficient(M)
            m = list(reversed(f.exponents()[0]))

            new_term = self._from_dict({_IntegerVectors(m).trim(): c})

            f -= new_term.expand()
            out += new_term

        return out

    def product(self, a, b):
        r"""
        Multiply the basis elements indexed by ``a`` and ``b``.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: K = KeyPolynomialBasis(QQ)
            sage: K([4,3,2]) * K([1,1,1])
            k[5, 4, 3]
        """
        p = a.expand() * b.expand()
        return self.from_polynomial(p)

def _divided_difference(P, i, f):
    r"""
    Apply the `i`th divided difference operator to the polynomial `f`.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis, _divided_difference
        sage: K = KeyPolynomialBasis(QQ)
        sage: R = K.polynomial_ring(); R
        Infinite polynomial ring in z over Rational Field
        sage: z, = R.gens()
        sage: f = z[1]*z[2]^3 + z[1]*z[2]*z[3]
        sage: _divided_difference(K, 2, f)
        z_3^2*z_1 + z_3*z_2*z_1 + z_2^2*z_1

    """
    out = P.polynomial_ring().zero()
    # linearly extend the divided difference on monomials
    for m in f.monomials():
        out += f.monomial_coefficient(m) * _divided_difference_on_monomial(P, i, m)

    return out

def _divided_difference_on_monomial(P, i, m):
    r"""
    Apply the `i`th divided difference operator to `m`.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis, _divided_difference_on_monomial
        sage: K = KeyPolynomialBasis(QQ)
        sage: R = K.polynomial_ring(); R
        Infinite polynomial ring in z over Rational Field
        sage: z, = R.gens()
        sage: _divided_difference_on_monomial(K, 1, z[1]*z[2]^3)
        -z_2^2*z_1 - z_2*z_1^2
        sage: _divided_difference_on_monomial(K, 2, z[1]*z[2]*z[3])
        0
        sage: _divided_difference_on_monomial(K, 3, z[1]*z[2]*z[3])
        z_2*z_1
        sage: _divided_difference_on_monomial(K, 3, z[1]*z[2]*z[4])
        -z_2*z_1
    """
    # The polynomial ring and generators
    R = P.polynomial_ring()
    x, = R.gens()

    # The exponent vector for the monomial
    exp = list(reversed(m.exponents()[0]))

    if i >= len(exp) - 1:
        if i == len(exp) - 1:
            exp = exp + [0]
        else:
            # if the transposition acts on two varibles which aren't
            # present, then the numerator is f - f == 0.
            return R.zero()

    # Apply the transposition s_i to the exponent vector
    exp[i+1], exp[i] = exp[i], exp[i+1]

    # Create the corresponding list of variables in the monomial
    terms = list(x[i]**j for i, j in enumerate(exp) if j)

    # Create the monomial from the list of variables
    si_m = R.prod(terms)

    # Create the numerator of the divided difference operator
    f = si_m - m

    # Division using the / operator is not implemented for
    # InfinitePolynomialRing, so we want to remove the factor of
    # x[i+1] - x[i]. If x[i+1] - x[i] is not a factor, it is because
    # the numerator is zero.
    try:
        factors = f.factor()
        factors_dict = dict(factors)
    except ArithmeticError:  # if f is zero already
        return R.zero()

    factors_dict[x[i + 1] - x[i]] = factors_dict[x[i + 1] - x[i]] - 1

    return R.prod(k**v for k, v in factors_dict.items()) * factors.unit()

def _pi(P, w, f):
    r"""
    Apply the operator `\pi_w` to the polynomial `f`.
    ``w`` may be either a single index or a list of
    indices of simple transpositions.

    .. WARNING::

        The simple transpositions should be applied from left to right. 

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis, _pi
        sage: K = KeyPolynomialBasis(QQ)
        sage: R = K.polynomial_ring()
        sage: z, = R.gens()
        sage: w = [2, 3]
        sage: _pi(K, w, z[1]^4*z[2]^2*z[3]*z[4])
        z_4^2*z_3*z_2*z_1^4 + z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4

        sage: w = [3, 2]
        sage: _pi(K, w, z[1]^4*z[2]^2*z[3]*z[4])
        z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4

        sage: w = [2, 1]
        sage: _pi(K, w, z[1]^2*z[2])
        z_3*z_2^2 + z_3*z_2*z_1 + z_3*z_1^2 + z_2^2*z_1 + z_2*z_1^2
    """
    if not hasattr(w, "__iter__"):  # this allows us to pass i instead of a word
        return _pi_i(P, w, f)
    for i in w:
        f = _pi_i(P, i, f)
    return f

def _pi_i(P, i, f):
    r"""
    Apply `\pi_i` for a single simple transposition `s_i = (i, i+1)`.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis, _pi_i
        sage: K = KeyPolynomialBasis(QQ)
        sage: R = K.polynomial_ring()
        sage: z, = R.gens()
        sage: _pi_i(K, 3, z[1]^4*z[2]^2*z[4])
        0

        sage: _pi_i(K, 2, z[1]^4*z[2]^2*z[3]*z[4])
        z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4
    """
    R = P.polynomial_ring()
    z, = R.gens()
    return _divided_difference(P, i, z[i] * f)

def _sorting_word(alpha):
    r"""
    Get a reduced word for the permutation which sorts ``alpha``
    into a partition. 

    The result is a list ``l = [i0, i1, i2, ...]`` where each ``ij``
    is a nonnegative integer such that it applies the simple
    transposition `(i_j, i_j+1)`. 

    # TODO:: reword the next paragraph
    The convention is that they are applied from zero index
    to ``len(l)``. This is the oposite of how they would be
    applied as function composition, where we would start with
    the rightmost simple transposition and work to the left.
    For instance, the result of the action of ``[0, 1]`` on 
    the integer
    vector ``[9,2,3]`` is ``[2,3,9]`` (**not** ``[3,9,2]``)

    EXAMPLES::

        sage: IV = IntegerVectors()
        sage: from sage.combinat.key_polynomial import _sorting_word
        sage: list(_sorting_word(IV([2,3,2])))
        [0]
        sage: list(_sorting_word(IV([5,6,7])))
        [0, 1, 0]
        sage: list(_sorting_word(IV([0,3,2])))
        [1, 0]
        sage: list(_sorting_word(IV([0,3,0,2])))
        [1, 2, 0]
        sage: list(_sorting_word(IV([3,2,1])))
        []
        sage: list(_sorting_word(IV([2,3,3])))
        [1, 0]
    """
    w = []
    n = len(alpha)

    # bubble sort to get the shortest sorting word
    with alpha.clone() as ac:
        for i in range(n-1):
            for j in range(n-i-1):
                if ac[j] < ac[j + 1]:
                    w.append(j)
                    ac[j], ac[j + 1] = ac[j + 1], ac[j]
    return reversed(w)

_IntegerVectors = IntegerVectors()
