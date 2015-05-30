# -*- coding: utf-8 -*-
r"""
Shuffle algebras

AUTHORS:

- Frédéric Chapoton (2013-03): Initial version
- Matthieu Deneufchatel (2013-07): Implemented dual PBW basis
"""

#*****************************************************************************
#  Copyright (C) 2013 Frédéric Chapoton <chapoton-math-univ-lyon1-fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.rings import Rings
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.coalgebras_with_basis import CoalgebrasWithBasis
from sage.categories.tensor import TensorProductsCategory, tensor
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word
from sage.combinat.words.finite_word import FiniteWord_class
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.sets.family import Family

class ShuffleAlgebra(CombinatorialFreeModule):
    r"""
    The shuffle algebra on some generators over a base ring.

    Shuffle algebras are commutative and associative algebras, with a
    basis indexed by words. The product of two words `w_1 \cdot w_2` is given
    by the sum over the shuffle product of `w_1` and `w_2`.

    .. SEEALSO::

        For more on shuffle products, see
        :mod:`~sage.combinat.words.shuffle_product` and
        :meth:`~sage.combinat.words.finite_word.FiniteWord_class.shuffle()`.

    REFERENCES:

    - :wikipedia:`Shuffle algebra`

    INPUT:

    -  ``R`` -- ring

    -  ``names`` -- generator names (string or an alphabet)

    EXAMPLES::

        sage: F = ShuffleAlgebra(QQ, 'xyz'); F
        Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Rational Field

        sage: mul(F.gens())
        B[word: xyz] + B[word: xzy] + B[word: yxz] + B[word: yzx] + B[word: zxy] + B[word: zyx]

        sage: mul([ F.gen(i) for i in range(2) ]) + mul([ F.gen(i+1) for i in range(2) ])
        B[word: xy] + B[word: yx] + B[word: yz] + B[word: zy]

        sage: S = ShuffleAlgebra(ZZ, 'abcabc'); S
        Shuffle Algebra on 3 generators ['a', 'b', 'c'] over Integer Ring
        sage: S.base_ring()
        Integer Ring

        sage: G = ShuffleAlgebra(S, 'mn'); G
        Shuffle Algebra on 2 generators ['m', 'n'] over Shuffle Algebra on 3 generators ['a', 'b', 'c'] over Integer Ring
        sage: G.base_ring()
        Shuffle Algebra on 3 generators ['a', 'b', 'c'] over Integer Ring

    Shuffle algebras commute with their base ring::

        sage: K = ShuffleAlgebra(QQ,'ab')
        sage: a,b = K.gens()
        sage: K.is_commutative()
        True
        sage: L = ShuffleAlgebra(K,'cd')
        sage: c,d = L.gens()
        sage: L.is_commutative()
        True
        sage: s = a*b^2 * c^3; s
        (12*B[word:abb]+12*B[word:bab]+12*B[word:bba])*B[word: ccc]
        sage: parent(s)
        Shuffle Algebra on 2 generators ['c', 'd'] over Shuffle Algebra on 2 generators ['a', 'b'] over Rational Field
        sage: c^3 * a * b^2
        (12*B[word:abb]+12*B[word:bab]+12*B[word:bba])*B[word: ccc]

    Shuffle algebras are commutative::

        sage: c^3 * b * a * b == c * a * c * b^2 * c
        True

    We can also manipulate elements in the basis and coerce elements from our
    base field::

        sage: F = ShuffleAlgebra(QQ, 'abc')
        sage: B = F.basis()
        sage: B[Word('bb')] * B[Word('ca')]
        B[word: bbca] + B[word: bcab] + B[word: bcba] + B[word: cabb] + B[word: cbab] + B[word: cbba]
        sage: 1 - B[Word('bb')] * B[Word('ca')] / 2
        B[word: ] - 1/2*B[word: bbca] - 1/2*B[word: bcab] - 1/2*B[word: bcba] - 1/2*B[word: cabb] - 1/2*B[word: cbab] - 1/2*B[word: cbba]
    """
    @staticmethod
    def __classcall_private__(cls, R, names):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = ShuffleAlgebra(QQ, 'xyz')
            sage: F2 = ShuffleAlgebra(QQ, ['x','y','z'])
            sage: F3 = ShuffleAlgebra(QQ, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        return super(ShuffleAlgebra, cls).__classcall__(cls, R, Alphabet(names))

    def __init__(self, R, names):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = ShuffleAlgebra(QQ, 'xyz'); F
            Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Rational Field
            sage: TestSuite(F).run()

        TESTS::

            sage: ShuffleAlgebra(24, 'toto')
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        self._alphabet = names
        self.__ngens = self._alphabet.cardinality()
        CombinatorialFreeModule.__init__(self, R, Words(names, infinite=False),
            latex_prefix="",
            category=(AlgebrasWithBasis(R), CommutativeAlgebras(R), CoalgebrasWithBasis(R)))

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R = ShuffleAlgebra(QQ,'xy')
            sage: R.variable_names()
            {'x', 'y'}
        """
        return self._alphabet

    def is_commutative(self):
        r"""
        Return ``True`` as the shuffle algebra is commutative.

        EXAMPLES::

            sage: R = ShuffleAlgebra(QQ,'x')
            sage: R.is_commutative()
            True
            sage: R = ShuffleAlgebra(QQ,'xy')
            sage: R.is_commutative()
            True
        """
        return True

    def _repr_(self):
        r"""
        Text representation of this shuffle algebra.

        EXAMPLES::

            sage: F = ShuffleAlgebra(QQ,'xyz')
            sage: F  # indirect doctest
            Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Rational Field

            sage: ShuffleAlgebra(ZZ,'a')
            Shuffle Algebra on one generator ['a'] over Integer Ring
        """
        if self.__ngens == 1:
            gen = "one generator"
        else:
            gen = "%s generators" %self.__ngens
        return "Shuffle Algebra on "+ gen +" %s over %s"%(
            self._alphabet.list(), self.base_ring())

    @cached_method
    def one_basis(self):
        r"""
        Return the empty word, which index of `1` of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = ShuffleAlgebra(QQ,'a')
            sage: A.one_basis()
            word:
            sage: A.one()
            B[word: ]
        """
        return self.basis().keys()([])

    def product_on_basis(self, w1, w2):
        r"""
        Return the product of basis elements ``w1`` and ``w2``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``w1``, ``w2`` -- Basis elements

        EXAMPLES::

            sage: A = ShuffleAlgebra(QQ,'abc')
            sage: W = A.basis().keys()
            sage: A.product_on_basis(W("acb"), W("cba"))
            B[word: acbacb] + B[word: acbcab] + 2*B[word: acbcba] + 2*B[word: accbab] + 4*B[word: accbba] + B[word: cabacb] + B[word: cabcab] + B[word: cabcba] + B[word: cacbab] + 2*B[word: cacbba] + 2*B[word: cbaacb] + B[word: cbacab] + B[word: cbacba]

            sage: (a,b,c) = A.algebra_generators()
            sage: a * (1-b)^2 * c
            2*B[word: abbc] - 2*B[word: abc] + 2*B[word: abcb] + B[word: ac] - 2*B[word: acb] + 2*B[word: acbb] + 2*B[word: babc] - 2*B[word: bac] + 2*B[word: bacb] + 2*B[word: bbac] + 2*B[word: bbca] - 2*B[word: bca] + 2*B[word: bcab] + 2*B[word: bcba] + B[word: ca] - 2*B[word: cab] + 2*B[word: cabb] - 2*B[word: cba] + 2*B[word: cbab] + 2*B[word: cbba]
        """
        return sum(self.basis()[u] for u in w1.shuffle(w2))

    def gen(self,i):
        r"""
        The ``i``-th generator of the algebra.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: F = ShuffleAlgebra(ZZ,'xyz')
            sage: F.gen(0)
            B[word: x]

            sage: F.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError("argument i (= %s) must be between 0 and %s"%(i, n-1))
        return self.algebra_generators()[i]

    def coproduct_on_basis(self, w):
        """
        Return the coproduct of the element of the basis indexed by
        the word ``w``.

        INPUT:

        - ``w`` -- a word

        EXAMPLES::

            sage: F = ShuffleAlgebra(QQ,'ab')
            sage: F.coproduct_on_basis(Word('a'))
            B[word: ] # B[word: a] + B[word: a] # B[word: ]
            sage: F.coproduct_on_basis(Word('aba'))
            B[word: ] # B[word: aba] + B[word: a] # B[word: ab] + B[word: a] # B[word: ba]
             + B[word: aa] # B[word: b] + B[word: ab] # B[word: a] + B[word: aba] # B[word: ]
             + B[word: b] # B[word: aa] + B[word: ba] # B[word: a]
            sage: F.coproduct_on_basis(Word())
            B[word: ] # B[word: ]
        """
        if len(w) == 0:
            return self.tensor_square().monomial((self.one_basis(), self.one_basis()))
        if len(w) == 1:
            return self.tensor_square().sum_of_terms([
                        ((w, self.one_basis()), 1),
                        ((self.one_basis(), w), 1) ], distinct=True)

        B = self.basis()
        result = self.coproduct_on_basis(Word([w[0]]))
        for i in w[1:]:
            result = self.tensor_square().sum_of_terms([
                    ((Word(v1)*Word(u1), Word(v2)*Word(u2)), coeff1 * coeff2)
                    for ((u1,u2),coeff1) in self.coproduct_on_basis(Word([i]))
                    for ((v1,v2),coeff2) in result ])
        return result

    def coproduct(self, S):
        """
        Return the coproduct of the series ``S``.

        EXAMPLES::

            sage: F = ShuffleAlgebra(QQ,'ab')
            sage: S = F.an_element(); S
            2*B[word: ] + 2*B[word: a] + 3*B[word: b]
            sage: F.coproduct(S)
            2*B[word: ] # B[word: ] + 2*B[word: ] # B[word: a] + 3*B[word: ] # B[word: b]
             + 2*B[word: a] # B[word: ] + 3*B[word: b] # B[word: ]
            sage: F.coproduct(F.one())
            B[word: ] # B[word: ]
        """
        return sum([c * self.coproduct_on_basis(i)
                    for i,c in S.monomial_coefficients().items()])

    def counit(self,S):
        """
        Return the counit of ``S``.

        EXAMPLES::

            sage: F = ShuffleAlgebra(QQ,'ab')
            sage: S = F.an_element(); S
            2*B[word: ] + 2*B[word: a] + 3*B[word: b]
            sage: F.counit(S)
            2
        """
        return S.monomial_coefficients().get(Word(), 0)

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        EXAMPLES::

            sage: A = ShuffleAlgebra(ZZ,'fgh'); A
            Shuffle Algebra on 3 generators ['f', 'g', 'h'] over Integer Ring
            sage: A.algebra_generators()
            Family (B[word: f], B[word: g], B[word: h])

            sage: A = ShuffleAlgebra(QQ, ['x1','x2'])
            sage: A.algebra_generators()
            Family (B[word: x1], B[word: x2])
        """
        Words = self.basis().keys()
        return Family( [self.monomial(Words([a])) for a in self._alphabet] )
        # FIXME: use this once the keys argument of FiniteFamily will be honoured
        # for the specifying the order of the elements in the family
        #return Family(self._alphabet, lambda a: self.term(self.basis().keys()(a)))

    gens = algebra_generators

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = ShuffleAlgebra(QQ,'xy')
            sage: x, y = R.gens()
            sage: R(3)
            3*B[word: ]
            sage: R(x)
            B[word: x]
            sage: R('xyy')
            B[word: xyy]
            sage: R(Word('xyx'))
            B[word: xyx]
        """
        if isinstance(x, (str, FiniteWord_class)):
            W = self.basis().keys()
            return self.monomial(W(x))

        P = x.parent()
        if isinstance(P, ShuffleAlgebra):
            if P is self:
                return x
            if not (P is self.base_ring()):
                return self.element_class(self, x.monomial_coefficients())
        if isinstance(P, DualPBWBasis):
            return self(P.expansion(x))
        # ok, not a shuffle algebra element (or should not be viewed as one).
        if isinstance(x, basestring):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x,locals=self.gens_dict())
        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self,{})
        else:
            return self.from_base_ring_from_one_basis(x)

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here is what canonically coerces to ``self``:

        - this shuffle algebra,

        - anything that coerces to the base ring of this shuffle algebra,

        - any shuffle algebra on the same variables, whose base ring
          coerces to the base ring of this shuffle algebra.

        EXAMPLES::

            sage: F = ShuffleAlgebra(GF(7), 'xyz'); F
            Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Finite Field of size 7

        Elements of the shuffle algebra canonically coerce in::

            sage: x, y, z = F.gens()
            sage: F.coerce(x*y) # indirect doctest
            B[word: xy] + B[word: yx]

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            B[word: ]

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Finite Field of size 7

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5*B[word: ]

        The shuffle algebra over `\ZZ` on `x, y, z` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = ShuffleAlgebra(ZZ,'xyz')
            sage: Gx,Gy,Gz = G.gens()
            sage: z = F.coerce(Gx**2 * Gy);z
            2*B[word: xxy] + 2*B[word: xyx] + 2*B[word: yxx]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the shuffle
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(x^3*y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Shuffle Algebra on 3 generators
            ['x', 'y', 'z'] over Finite Field of size 7 to Shuffle Algebra on 3
            generators ['x', 'y', 'z'] over Integer Ring
        """
        try:
            R = x.parent()

            # shuffle algebras in the same variables over any base
            # that coerces in:
            if isinstance(R,ShuffleAlgebra):
                if R.variable_names() == self.variable_names():
                    if self.has_coerce_map_from(R.base_ring()):
                        return self(x)
                    else:
                        raise TypeError("no natural map between bases of shuffle algebras")

            if isinstance(R, DualPBWBasis):
                return self(R.expansion(x))

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this shuffle algebra.
        return self._coerce_try(x, [self.base_ring()])

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - Shuffle Algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``.

        - Anything with a coercion into ``self.base_ring()``.

        TESTS::

            sage: F = ShuffleAlgebra(ZZ, 'xyz')
            sage: G = ShuffleAlgebra(QQ, 'xyz')
            sage: H = ShuffleAlgebra(ZZ, 'y')
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            False
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
            sage: F._coerce_map_from_(F.dual_pbw_basis())
            True
        """
        # shuffle algebras in the same variable over any base that coerces in:
        if isinstance(R, ShuffleAlgebra):
            if R.variable_names() == self.variable_names():
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
                else:
                    return False

        if isinstance(R, DualPBWBasis):
            return self.has_coerce_map_from(R._alg)

        return self.base_ring().has_coerce_map_from(R)

    def dual_pbw_basis(self):
        """
        Return the dual PBW of ``self``.

        EXAMPLES::

            sage: A = ShuffleAlgebra(QQ, 'ab')
            sage: A.dual_pbw_basis()
            The dual Poincare-Birkhoff-Witt basis of Shuffle Algebra on 2 generators ['a', 'b'] over Rational Field
        """
        return DualPBWBasis(self.base_ring(), self._alphabet)

    def to_dual_pbw_element(self, w):
        """
        Return the element `w` of ``self`` expressed in the dual PBW basis.

        INPUT:

        - ``w`` -- an element of the shuffle algebra

        EXAMPLES::

            sage: A = ShuffleAlgebra(QQ, 'ab')
            sage: f = 2 * A(Word()) + A(Word('ab')); f
            2*B[word: ] + B[word: ab]
            sage: A.to_dual_pbw_element(f)
            2*S[word: ] + S[word: ab]
            sage: A.to_dual_pbw_element(A.one())
            S[word: ]
            sage: S = A.dual_pbw_basis()
            sage: elt = S.expansion_on_basis(Word('abba')); elt
            2*B[word: aabb] + B[word: abab] + B[word: abba]
            sage: A.to_dual_pbw_element(elt)
            S[word: abba]
            sage: A.to_dual_pbw_element(2*A(Word('aabb')) + A(Word('abab')))
            S[word: abab]
            sage: S.expansion(S('abab'))
            2*B[word: aabb] + B[word: abab]
        """
        D = self.dual_pbw_basis()
        l = {}
        W = self.basis().keys()
        while w != self.zero():
            support = [W(i[0]) for i in list(w)]
            min_elt = W(support[0])
            if len(support) > 1:
                for word in support[1:len(support)-1]:
                    if min_elt.lex_less(word):
                        min_elt = W(word)
            coeff = list(w)[support.index(min_elt)][1]
            l[min_elt] = l.get(min_elt, 0) + coeff
            w = w - coeff * D.expansion_on_basis(W(min_elt))

        return D.sum_of_terms([(m, c) for m,c in l.items() if c != 0])

class DualPBWBasis(CombinatorialFreeModule):
    r"""
    The basis dual to the Poincare-Birkhoff-Witt basis of the free algebra.

    We recursively define the dual PBW basis as the basis of the
    shuffle algebra given by

    .. MATH::

        S_w = \begin{cases}
        w & |w| = 1, \\
        x S_u & w = xu \text{ and } w \in \mathrm{Lyn}(X), \\
        \displaystyle \frac{S_{\ell_{i_1}}^{\ast \alpha_1} \ast \cdots
        \ast S_{\ell_{i_k}}^{\ast \alpha_k}}{\alpha_1! \cdots \alpha_k!} &
        w = \ell_{i_1}^{\alpha_1} \cdots \ell_{i_k}^{\alpha_k} \text{ with }
        \ell_1 > \cdots > \ell_k \in \mathrm{Lyn}(X).
        \end{cases}

    where `S \ast T` denotes the shuffle product of `S` and `T` and
    `\mathrm{Lyn}(X)` is the set of Lyndon words in the alphabet `X`.

    The definition may be found in Theorem 5.3 of [Reuten1993]_.

    INPUT:

    - ``R`` -- ring

    - ``names`` -- names of the generators (string or an alphabet)

    REFERENCES:

    .. [Reuten1993] C. Reutenauer. *Free Lie Algebras*. Number 7 in
       London Math. Soc. Monogr. (N.S.). Oxford University Press. (1993).

    EXAMPLES::

        sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
        sage: S
        The dual Poincare-Birkhoff-Witt basis of Shuffle Algebra on 2 generators ['a', 'b'] over Rational Field
        sage: S.one()
        S[word: ]
        sage: S.one_basis()
        word:
        sage: T = ShuffleAlgebra(QQ, 'abcd').dual_pbw_basis(); T
        The dual Poincare-Birkhoff-Witt basis of Shuffle Algebra on 4 generators ['a', 'b', 'c', 'd'] over Rational Field
        sage: T.algebra_generators()
        (S[word: a], S[word: b], S[word: c], S[word: d])

    TESTS:

    We check conversion between the bases::

        sage: A = ShuffleAlgebra(QQ, 'ab')
        sage: S = A.dual_pbw_basis()
        sage: W = Words('ab', 5)
        sage: all(S(A(S(w))) == S(w) for w in W)
        True
        sage: all(A(S(A(w))) == A(w) for w in W)
        True
    """
    @staticmethod
    def __classcall_private__(cls, R, names):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.algebras.shuffle_algebra import DualPBWBasis
            sage: D1 = DualPBWBasis(QQ, 'ab')
            sage: D2 = DualPBWBasis(QQ, Alphabet('ab'))
            sage: D1 is D2
            True
        """
        return super(DualPBWBasis, cls).__classcall__(cls, R, Alphabet(names))

    def __init__(self, R, names):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: D = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: TestSuite(D).run()
        """
        self._alphabet = names
        self._alg = ShuffleAlgebra(R, names)
        CombinatorialFreeModule.__init__(self, R, Words(names), prefix='S',
            category=(AlgebrasWithBasis(R), CommutativeAlgebras(R), CoalgebrasWithBasis(R)))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            The dual Poincare-Birkhoff-Witt basis of Shuffle Algebra on 2 generators ['a', 'b'] over Rational Field
        """
        return "The dual Poincare-Birkhoff-Witt basis of {}".format(self._alg)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: A = ShuffleAlgebra(QQ, 'ab')
            sage: S = A.dual_pbw_basis()
            sage: S('abaab')
            S[word: abaab]
            sage: S(Word('aba'))
            S[word: aba]
            sage: S(A('ab'))
            S[word: ab]
        """
        if isinstance(x, (str, FiniteWord_class)):
            W = self.basis().keys()
            x = W(x)
        elif isinstance(x.parent(), ShuffleAlgebra):
            return self._alg.to_dual_pbw_element(self._alg(x))
        return super(DualPBWBasis, self)._element_constructor_(x)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion from ``R`` into ``self`` and
        ``False`` otherwise. The things that coerce into ``self`` are:

        - Anything that coerces into the associated shuffle algebra of ``self``

        TESTS::

            sage: F = ShuffleAlgebra(ZZ, 'xyz').dual_pbw_basis()
            sage: G = ShuffleAlgebra(QQ, 'xyz').dual_pbw_basis()
            sage: H = ShuffleAlgebra(ZZ, 'y').dual_pbw_basis()
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            False
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
            sage: F._coerce_map_from_(F._alg)
            True
        """
        return self._alg.has_coerce_map_from(R)

    def one_basis(self):
        """
        Return the indexing element of the basis element `1`.

        EXAMPLES::

            sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: S.one_basis()
            word:
        """
        W = self.basis().keys()
        return W([])

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: S.algebra_generators()
            (S[word: a], S[word: b])
        """
        W = self.basis().keys()
        return tuple(self.monomial(W(a)) for a in self._alphabet)

    gens = algebra_generators

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: S.gen(0)
            S[word: a]
            sage: S.gen(1)
            S[word: b]
        """
        return self.algebra_generators()[i]

    def shuffle_algebra(self):
        """
        Return the associated shuffle algebra of ``self``.

        EXAMPLES::

            sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: S.shuffle_algebra()
            Shuffle Algebra on 2 generators ['a', 'b'] over Rational Field
        """
        return self._alg

    def product(self, u, v):
        """
        Return the product of two elements ``u`` and ``v``.

        EXAMPLES::

            sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: a,b = S.gens()
            sage: S.product(a, b)
            S[word: ba]
            sage: S.product(b, a)
            S[word: ba]
            sage: S.product(b^2*a, a*b*a)
            36*S[word: bbbaaa]

        TESTS:

        Check that multiplication agrees with the multiplication in the
        shuffle algebra::

            sage: A = ShuffleAlgebra(QQ, 'ab')
            sage: S = A.dual_pbw_basis()
            sage: a,b = S.gens()
            sage: A(a*b)
            B[word: ab] + B[word: ba]
            sage: A(a*b*a)
            2*B[word: aab] + 2*B[word: aba] + 2*B[word: baa]
            sage: S(A(a)*A(b)*A(a)) == a*b*a
            True
        """
        return self(self.expansion(u) * self.expansion(v))

    @lazy_attribute
    def expansion(self):
        """
        Return the morphism corresponding to the expansion into words of
        the shuffle algebra.

        EXAMPLES::

            sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: f = S('ab') + S('aba')
            sage: S.expansion(f)
            2*B[word: aab] + B[word: ab] + B[word: aba]
        """
        return self.module_morphism(self.expansion_on_basis, codomain=self._alg)

    @cached_method
    def expansion_on_basis(self, w):
        r"""
        Return the expansion of `S_w` in words of the shuffle algebra.

        INPUT:

        - ``w`` -- a word

        EXAMPLES::

            sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
            sage: S.expansion_on_basis(Word())
            B[word: ]
            sage: S.expansion_on_basis(Word()).parent()
            Shuffle Algebra on 2 generators ['a', 'b'] over Rational Field
            sage: S.expansion_on_basis(Word('abba'))
            2*B[word: aabb] + B[word: abab] + B[word: abba]
            sage: S.expansion_on_basis(Word())
            B[word: ]
            sage: S.expansion_on_basis(Word('abab'))
            2*B[word: aabb] + B[word: abab]
        """
        from sage.functions.other import factorial
        if len(w) == 0:
            return self._alg.one()
        if len(w) == 1:
           return self._alg.monomial(w)

        if w.is_lyndon():
            W = self.basis().keys()
            letter = W(w[0])
            expansion = self.expansion_on_basis(W(w[1:]))
            return self._alg.sum_of_terms([(letter * i, c) for i,c in expansion])

        lf = w.lyndon_factorization()
        powers = {}
        for i in lf:
            powers[i] = powers.get(i, 0) + 1
        denom = prod(factorial(p) for p in powers.values())
        result = self._alg.prod(self.expansion_on_basis(i) for i in lf)
        return self._alg(result / denom)

    class Element(CombinatorialFreeModule.Element):
        """
        An element in the dual PBW basis.
        """
        def expand(self):
            """
            Expand ``self`` in words of the shuffle algebra.

            EXAMPLES::

                sage: S = ShuffleAlgebra(QQ, 'ab').dual_pbw_basis()
                sage: f = S('ab') + S('bab')
                sage: f.expand()
                B[word: ab] + 2*B[word: abb] + B[word: bab]
            """
            return self.parent().expansion(self)

