# sage.doctest: needs sage.combinat
r"""
F-algebra for motivic multiple zeta values.

This is a commutative algebra, defined as the tensor product of the
polynomial algebra over one generator `f_2` by the shuffle algebra in
infinitely many generators indexed by odd integers and denoted by
`f_3`, `f_5`, ... It serves as an auxiliary algebra in the study of
the ring of motivic multiple zeta values.

Here we provide a basic direct implementation, endowed with the
motivic coproduct.

The similar algebra where the shuffle algebra has generators
`f_1, f_3, f_5, \ldots` is now also available. The implementation is even more
general, allowing any positive odd integer as start index.

AUTHORS:

- Frédéric Chapoton (2022-09): Initial version
"""
# ****************************************************************************
#  Copyright (C) 2022 Frédéric Chapoton <chapoton-unistra-fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from collections.abc import Iterator

from sage.arith.misc import bernoulli
from sage.categories.rings import Rings
from sage.categories.bialgebras_with_basis import BialgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
from sage.combinat.words.finite_word import FiniteWord_class
from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2 as shuffle
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.integer_range import IntegerRange
from sage.rings.integer_ring import ZZ
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.rings.infinity import Infinity
from sage.modules.free_module_element import vector


def W_Odds(start=3):
    r"""
    Indexing set for the odd generators.

    This is the set of pairs
    (integer power of `f_2`, word in `s, s+2, s+4, \ldots`)
    where `s` is the chosen odd start index.

    INPUT:

    - ``start`` -- (default: ``3``) odd start index for odd generators

    EXAMPLES::

        sage: from sage.modular.multiple_zeta_F_algebra import W_Odds
        sage: W_Odds(3)
        Finite words over {3, 5, ...}
    """
    return Words(IntegerRange(start, Infinity, 2), infinite=False)


def str_to_index(x: str) -> tuple:
    r"""
    Convert a string to an index.

    Every letter ``'2'`` contributes to the power of `f_2`. Other letters
    are odd and define a word in `f_1, f_3, f_5, \ldots`

    Usually the letters ``'2'`` form a prefix of the input.

    EXAMPLES::

        sage: from sage.modular.multiple_zeta_F_algebra import str_to_index
        sage: str_to_index("22357")
        (2, [3, 5, 7])

        sage: str_to_index("22317")
        (2, [3, 1, 7])
    """
    p = x.count("2")
    w = [int(i) for i in x if i != '2']
    return (p, w)


def basis_f_odd_iterator(n, start=3) -> Iterator[tuple]:
    r"""
    Return an iterator over compositions of `n` with odd parts.

    Let `s` be the chosen odd start index. The allowed parts are the
    odd integers at least equal to `s`, in the set `s,s+2,s+4,s+6,\ldots`.

    This set of compositions is used to index a basis.

    INPUT:

    - ``n`` -- integer

    - ``start`` -- odd integer (default: `3`); start index for odd generators

    EXAMPLES::

        sage: from sage.modular.multiple_zeta_F_algebra import basis_f_odd_iterator
        sage: [list(basis_f_odd_iterator(i)) for i in range(2,9)]
        [[], [(3,)], [], [(5,)], [(3, 3)], [(7,)], [(5, 3), (3, 5)]]
        sage: list(basis_f_odd_iterator(14))
        [(11, 3),
         (5, 3, 3, 3),
         (3, 5, 3, 3),
         (3, 3, 5, 3),
         (9, 5),
         (3, 3, 3, 5),
         (7, 7),
         (5, 9),
         (3, 11)]
    """
    if n == 0:
        yield ()
        return
    if n % 2 and n >= start:
        yield (n,)
    for k in range(start, n, 2):
        for word in basis_f_odd_iterator(n - k, start=start):
            yield word + (k, )


def basis_f_iterator(n, start=3) -> Iterator[tuple]:
    r"""
    Return an iterator for decompositions of `n` using `2` and odd integers.

    Let `s` be the chosen odd start index. The allowed odd parts are the
    odd integers at least equal to `s`, in the set `s,s+2,s+4,s+6,\ldots`.

    The means that each term is made of a power of 2 and a composition
    of the remaining integer with parts in `(s,s+2,s+4,\ldots)`.

    This set is indexing a basis of the homogeneous component of weight ``n``.

    INPUT:

    - ``n`` -- integer

    - ``start`` -- (default: `3`) odd start index for odd generators

    Each term is returned as a pair (integer, word) where
    the integer is the exponent of 2.

    EXAMPLES::

        sage: from sage.modular.multiple_zeta_F_algebra import basis_f_iterator
        sage: [list(basis_f_iterator(i)) for i in range(2,9)]
        [[(1, word: )],
         [(0, word: 3)],
         [(2, word: )],
         [(0, word: 5), (1, word: 3)],
         [(0, word: 33), (3, word: )],
         [(0, word: 7), (1, word: 5), (2, word: 3)],
         [(0, word: 53), (0, word: 35), (1, word: 33), (4, word: )]]
        sage: list(basis_f_iterator(11))
        [(0, word: 11),
         (0, word: 533),
         (0, word: 353),
         (0, word: 335),
         (1, word: 9),
         (1, word: 333),
         (2, word: 7),
         (3, word: 5),
         (4, word: 3)]

    TESTS::

        sage: list(basis_f_iterator(0))
        [(0, word: )]
        sage: list(basis_f_iterator(3, start=1))
        [(0, word: 3), (0, word: 111), (1, word: 1)]
    """
    wodds = W_Odds(start)
    for k in range(n // 2 + 1):
        for word in basis_f_odd_iterator(n - 2 * k, start):
            yield (k, wodds(word, check=False))


def morphism_constructor(data: dict, start=3):
    r"""
    Build a morphism from the F-algebra to some codomain.

    Let `s` be the chosen odd start index.

    INPUT:

    - ``data`` -- dictionary with integer keys containing the images of
      `f_2, f_s, f_{s+2}, f_{s+4}, \ldots`

    - ``start`` -- (default: 3) start index for odd generators

    OUTPUT: the unique morphism defined by the dictionary ``data``

    The codomain must be a zinbiel algebra, namely have both a
    commutative associative product ``*`` and a zinbiel product
    available as ``half_product``.

    EXAMPLES::

        sage: from sage.modular.multiple_zeta_F_algebra import F_algebra, morphism_constructor
        sage: Z = Multizeta
        sage: D = {2: Z(2), 3: Z(3)}
        sage: rho = morphism_constructor(D)
        sage: F = rho.domain()
        sage: rho(F("2"))
        ζ(2)
        sage: rho(F("3"))
        ζ(3)
        sage: rho(F("33"))
        6*ζ(1,5) + 3*ζ(2,4) + ζ(3,3)
        sage: rho(F("23"))
        6*ζ(1,4) + 3*ζ(2,3) + ζ(3,2)
    """
    im_f2 = data[2]
    codomain = im_f2.parent()
    domain = F_algebra(codomain.base_ring(), start=start)

    def morphism_on_basis(pw):
        p, w = pw
        if not w:
            return im_f2**p
        v = im_f2**p * data[w[-1]]
        for letter in w[-2::-1]:
            v = codomain.half_product(data[letter], v)
        return v

    morphism = domain._module_morphism(morphism_on_basis, codomain=codomain)

    return morphism


class F_algebra(CombinatorialFreeModule):
    r"""
    Auxiliary algebra for the study of motivic multiple zeta values.

    INPUT:

    - ``R`` -- ring

    - ``start`` -- (default: ``3``) odd start index for odd generators

    EXAMPLES::

        sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
        sage: F = F_algebra(QQ); F
        F-ring over Rational Field
        sage: F.base_ring()
        Rational Field
        sage: F.is_commutative()
        True
        sage: TestSuite(F).run()

        sage: f2 = F("2")
        sage: f3 = F("3")
        sage: f5 = F("5")

        sage: s = f2*f3+f5; s
        f5 + f2*f3
    """
    def __init__(self, R, start=3):
        r"""
        Initialize ``self``.

        INPUT:

        - ``R`` -- base ring

        - ``start`` -- (default: ``3``) odd start index for odd generators

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(QQ); F
            F-ring over Rational Field

        TESTS::

            sage: F_algebra(24)
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        if not start % 2 and start > 0:
            raise ValueError("argument start must be odd and positive")
        self._start = start
        Indices = NonNegativeIntegers().cartesian_product(W_Odds(start))
        cat = BialgebrasWithBasis(R).Commutative().Graded()
        CombinatorialFreeModule.__init__(self, R, Indices,
                                         latex_prefix='', prefix='f',
                                         category=cat)

    def _repr_term(self, pw) -> str:
        r"""
        Return the custom representation of terms.

        Each monomial is written as a power of `f_2` times a word
        in `f_1, f_3, f_5, \ldots`.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(QQ)
            sage: f2 = F.gen(2)
            sage: f3 = F.gen(3)
            sage: f5 = F.gen(5)
            sage: f2*f3+f5+f2**2  # indirect doctest
            f5 + f2*f3 + f2^2
        """
        p, w = pw
        if not p:
            if not w:
                return "1"
            resu = ""
        elif p == 1:
            resu = "f2"
        else:
            resu = f"f2^{p}"
        if p and w:
            resu += "*"
        return resu + "".join(f"f{i}" for i in w)

    def _repr_(self) -> str:
        r"""
        Text representation of this algebra.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(ZZ)
            sage: F  # indirect doctest
            F-ring over Integer Ring
        """
        return f"F-ring over {self.base_ring()}"

    @cached_method
    def one_basis(self):
        r"""
        Return the pair (0, empty word), which index of `1` of this algebra.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: A = F_algebra(QQ)
            sage: A.one_basis()
            (0, word: )
        """
        return self.basis().keys()((0, []))

    def product_on_basis(self, pw1, pw2):
        r"""
        Return the product of basis elements ``pw1`` and ``pw2``.

        INPUT:

        - ``pw1``, ``pw2`` -- basis elements

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: A = F_algebra(QQ)
            sage: W = A.basis().keys()
            sage: A.product(A("23"), A("25"))  # indirect doctest
            f2^2*f3f5 + f2^2*f5f3
        """
        p1, w1 = pw1
        p2, w2 = pw2
        p = p1 + p2
        return self.sum_of_monomials((p, u) for u in w1.shuffle(w2))

    def half_product_on_basis(self, pw1, pw2):
        r"""
        Return the half product of basis elements ``pw1`` and ``pw2``.

        This is an extension of the zinbiel product of the shuffle algebra.

        INPUT:

        - ``pw1``, ``pw2`` -- basis elements

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: A = F_algebra(QQ)
            sage: W = A.basis().keys()
            sage: t = A.half_product(A("23"), A("25")); t  # indirect doctest
            f2^2*f3f5

        TESTS::

            sage: [list(pw[1]) for pw, cf in t]
            [[3, 5]]
        """
        p1, w1 = pw1
        p2, w2 = pw2
        p = p1 + p2
        if not w1:
            return self.basis()[(p, w2)]
        letter = w1[:1]
        return self.sum_of_monomials((p, letter + u)
                                     for u in w1[1:].shuffle(w2))

    @lazy_attribute
    def half_product(self):
        r"""
        Return the `<` product.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: A = F_algebra(QQ)
            sage: W = A.basis().keys()
            sage: A.half_product(A("235"), A("227"))
            f2^3*f3f5f7 + f2^3*f3f7f5
        """
        half = self.half_product_on_basis
        return self._module_morphism(self._module_morphism(half, position=0,
                                                           codomain=self),
                                     position=1)

    def gen(self, i):
        r"""
        Return the generator of the F ring over `\QQ`.

        INPUT:

        - ``i`` -- nonnegative integer (at least 2)

        If ``i`` is odd, this returns a single generator `f_i` of the free
        shuffle algebra.

        Otherwise, it returns an appropriate multiple of a power of `f_2`.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: A = F_algebra(QQ)
            sage: [A.gen(i) for i in range(2,8)]
            [f2, f3, 2/5*f2^2, f5, 8/35*f2^3, f7]
        """
        f2 = self.monomial(self._indices((1, [])))
        if i == 2:
            return f2
        # now i odd >= start
        if i % 2:
            return self.monomial(self._indices((0, [i])))
        # now powers of f2
        i = i // 2
        B = bernoulli(2 * i) * (-1)**(i - 1)
        B *= ZZ(2)**(3 * i - 1) * ZZ(3)**i / ZZ(2 * i).factorial()
        return B * f2**i

    def an_element(self):
        """
        Return a typical element.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(ZZ)
            sage: F.an_element()
            3*f2*f3f5 + f2*f5f3
        """
        return self("253") + 3 * self("235")

    def some_elements(self):
        """
        Return some typical elements.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(ZZ)
            sage: F.some_elements()
            [0, 1, f2, f3 + f5]
        """
        return [self.zero(), self.one(), self.gen(2),
                self.gen(3) + self.gen(5)]

    def coproduct_on_basis(self, pw):
        r"""
        Return the coproduct of the basis element indexed by the pair ``pw``.

        The coproduct is given by deconcatenation on the shuffle part,
        and extended by the value

        .. MATH::

            \Delta(f_2) = 1 \otimes f_2.

        INPUT:

        - ``pw`` -- an index

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(QQ)
            sage: W = F.basis().keys()
            sage: F.coproduct_on_basis(W((1,[])))
            1 # f2
            sage: F.coproduct_on_basis(W((0,[3])))
            1 # f3 + f3 # 1
            sage: F.coproduct_on_basis(W((1,[3])))
            1 # f2*f3 + f3 # f2
            sage: F.coproduct_on_basis(W((0,[3,5])))
            1 # f3f5 + f3 # f5 + f3f5 # 1
            sage: F.coproduct_on_basis(W((0,[])))
            1 # 1

        TESTS::

            sage: F = F_algebra(QQ)
            sage: S = F.an_element(); S
            3*f2*f3f5 + f2*f5f3
            sage: F.coproduct(S)
            3*1 # f2*f3f5 + 1 # f2*f5f3 + 3*f3 # f2*f5 + 3*f3f5 # f2
            + f5 # f2*f3 + f5f3 # f2
        """
        p, w = pw
        TS = self.tensor_square()
        return TS.sum_of_monomials(((0, w[:i]), (p, w[i:]))
                                   for i in range(len(w) + 1))

    def degree_on_basis(self, pw):
        """
        Return the degree of the element ``w``.

        This is the sum of the power of `f_2` and the indices in the word.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: A = F_algebra(QQ)
            sage: [A.degree_on_basis(x.leading_support()) for x in A.some_elements() if x != 0]
            [0, 1, 5]
        """
        p, w = pw
        return ZZ(p + sum(w))

    def homogeneous_from_vector(self, vec, N):
        """
        Convert back a vector to an element of the F-algebra.

        INPUT:

        - ``vec`` -- a vector with coefficients in some base ring

        - ``N`` -- integer; the homogeneous weight

        OUTPUT: a homogeneous element of :func:`F_ring` over this base ring

        .. SEEALSO:: :meth:`F_algebra.homogeneous_to_vector`

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(QQ)
            sage: F.homogeneous_from_vector((4,5),6)
            4*f3f3 + 5*f2^3
            sage: _.homogeneous_to_vector()
            (4, 5)
        """
        if isinstance(vec, (list, tuple)):
            vec = vector(vec)
        return self.sum(cf * self.monomial(bi)
                        for cf, bi in zip(vec, basis_f_iterator(N, self._start)))

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: R = F_algebra(QQ)
            sage: R("3")
            f3
            sage: R("2")
            f2
            sage: R("2235")
            f2^2*f3f5
        """
        if isinstance(x, (str, FiniteWord_class)):
            return self.monomial(self._indices(str_to_index(x)))

        P = x.parent()
        if isinstance(P, F_algebra):
            if P is self:
                return x
            if P is not self.base_ring():
                return self.element_class(self, x.monomial_coefficients())

        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self, {})
        return self.from_base_ring_from_one_basis(x)

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - an ``F_algebra`` over a base with a coercion
          map into ``self.base_ring()``.

        - Anything with a coercion into ``self.base_ring()``.

        EXAMPLES::

            sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
            sage: F = F_algebra(GF(7)); F
            F-ring over Finite Field of size 7

        Elements of the algebra itself canonically coerce in::

            sage: F.coerce(F("2")*F("3")) # indirect doctest
            f2*f3

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            1

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field
            to F-ring over Finite Field of size 7

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5

        The algebra over `\ZZ` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = F_algebra(ZZ)
            sage: Gx,Gy = G.gen(2), G.gen(3)
            sage: z = F.coerce(Gx**2 * Gy);z
            f2^2*f3
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(F("2"))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from F-ring over Finite Field
            of size 7 to F-ring over Integer Ring

        TESTS::

            sage: F = F_algebra(ZZ)
            sage: G = F_algebra(QQ)

            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True

            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True

            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        if isinstance(R, F_algebra):
            return self.base_ring().has_coerce_map_from(R.base_ring())

        return self.base_ring().has_coerce_map_from(R)

    class Element(CombinatorialFreeModule.Element):
        def coefficient(self, w):
            """
            Return the coefficient of the given index.

            EXAMPLES::

                sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
                sage: F = F_algebra(QQ)
                sage: S = F.an_element(); S
                3*f2*f3f5 + f2*f5f3
                sage: S.coefficient("235")
                3
                sage: S.coefficient((1,[5,3]))
                1
            """
            if isinstance(w, str):
                w = str_to_index(w)
            w = self.parent()._indices(w)
            return super().coefficient(w)

        def homogeneous_to_vector(self):
            """
            Convert an homogeneous element to a vector.

            This is using a fixed enumeration of the basis.

            OUTPUT: a vector with coefficients in the base ring

            .. SEEALSO:: :meth:`F_algebra.homogeneous_from_vector`

            EXAMPLES::

                sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
                sage: F = F_algebra(QQ)
                sage: f2 = F("2")
                sage: x = f2**4 + 34 * F("233")
                sage: x.homogeneous_to_vector()
                (0, 0, 34, 1)
                sage: x.coefficients()
                [34, 1]

            TESTS::

                sage: x = F.monomial(F._indices((0,[11]))); x
                f11
                sage: x.homogeneous_to_vector()
                (1, 0, 0, 0, 0, 0, 0, 0, 0)
            """
            F = self.parent()
            BR = F.base_ring()
            if not self:
                return vector(BR, [])
            a, b = next(iter(self))[0]
            N = 2 * a + sum(int(x) for x in b)
            return vector(BR, [self.coefficient(b)
                               for b in basis_f_iterator(N, F._start)])

        def without_f2(self):
            """
            Remove all terms containing a power of `f_2`.

            EXAMPLES::

                sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
                sage: F = F_algebra(QQ)
                sage: t = 4 * F("35") + F("27")
                sage: t.without_f2()
                4*f3f5
            """
            F = self.parent()
            return F._from_dict({(0, w): cf for (p, w), cf in self if not p})

        def single_valued(self):
            """
            Return the single-valued version of ``self``.

            EXAMPLES::

                sage: from sage.modular.multiple_zeta_F_algebra import F_algebra
                sage: F = F_algebra(QQ)
                sage: t = 4 * F("2") + F("3")
                sage: t.single_valued()
                2*f3
                sage: t = 4 * F("35") + F("27")
                sage: t.single_valued()
                8*f3f5 + 8*f5f3
            """
            F = self.parent()
            no_f2 = self.without_f2()
            return F.sum_of_terms(((0, w), cf)
                                  for (a, b), cf in no_f2.coproduct()
                                  for w in shuffle(a[1], b[1].reversal(), False))
