# sage.doctest: needs sage.combinat sage.modules
r"""
Free dendriform algebras

AUTHORS:

- Frédéric Chapoton (2017)
"""
# ****************************************************************************
#       Copyright (C) 2010-2015 Frédéric Chapoton <chapoton@unistra.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.hopf_algebras import HopfAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.binary_tree import (BinaryTrees, BinaryTree,
                                       LabelledBinaryTrees,
                                       LabelledBinaryTree)
from sage.categories.pushout import (ConstructionFunctor,
                                     CompositeConstructionFunctor,
                                     IdentityConstructionFunctor)
from sage.categories.rings import Rings
from sage.categories.functor import Functor
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.structure.coerce_exceptions import CoercionException
from sage.rings.infinity import Infinity


class FreeDendriformAlgebra(CombinatorialFreeModule):
    r"""
    The free dendriform algebra.

    Dendriform algebras are associative algebras, where the associative
    product `*` is decomposed as a sum of two binary operations

    .. MATH::

        x * y = x \succ y + x \prec y

    that satisfy the axioms:

    .. MATH::

        (x \succ y) \prec z = x \succ (y \prec z),

    .. MATH::

        (x \prec y) \prec z = x \prec (y * z).

    .. MATH::

        (x * y) \succ z = x \succ (y \succ z).

    The free Dendriform algebra on a given set `E` has an explicit
    description using (planar) binary trees, just as the free
    associative algebra can be described using words. The underlying
    vector space has a basis indexed by finite binary trees endowed
    with a map from their vertices to `E`. In this basis, the
    associative product of two (decorated) binary trees `S * T` is the
    sum over all possible ways of identifying (glueing) the rightmost path in
    `S` and the leftmost path in `T`.

    The decomposition of the associative product as the sum of two
    binary operations `\succ` and
    `\prec` is made by separating the terms according to the origin of
    the root vertex. For `x \succ y`, one keeps the terms where the root
    vertex comes from `y`, whereas for `x \prec y` one keeps the terms
    where the root vertex comes from `x`.

    The free dendriform algebra can also be considered as the free
    algebra over the Dendriform operad.

    .. NOTE::

        The usual binary operator `*` is used for the
        associative product.

    EXAMPLES::

        sage: F = algebras.FreeDendriform(ZZ, 'xyz')
        sage: x,y,z = F.gens()
        sage: (x * y) * z
        B[x[., y[., z[., .]]]] + B[x[., z[y[., .], .]]] + B[y[x[., .], z[., .]]]
         + B[z[x[., y[., .]], .]] + B[z[y[x[., .], .], .]]

    The free dendriform algebra is associative::

        sage: x * (y * z) == (x * y) * z
        True

    The associative product decomposes in two parts::

        sage: x * y == F.prec(x, y) + F.succ(x, y)
        True

    The axioms hold::

        sage: F.prec(F.succ(x, y), z) == F.succ(x, F.prec(y, z))
        True
        sage: F.prec(F.prec(x, y), z) == F.prec(x, y * z)
        True
        sage: F.succ(x * y, z) == F.succ(x, F.succ(y, z))
        True

    When there is only one generator, unlabelled trees are used instead::

        sage: F1 = algebras.FreeDendriform(QQ)
        sage: w = F1.gen(0); w
        B[[., .]]
        sage: w * w * w
        B[[., [., [., .]]]] + B[[., [[., .], .]]] + B[[[., .], [., .]]]
         + B[[[., [., .]], .]] + B[[[[., .], .], .]]

    The set `E` can be infinite::

        sage: F = algebras.FreeDendriform(QQ, ZZ)
        sage: w = F.gen(1); w
        B[1[., .]]
        sage: x = F.gen(2); x
        B[-1[., .]]
        sage: w*x
        B[-1[1[., .], .]] + B[1[., -1[., .]]]

    REFERENCES:

    - [LR1998]_
    """
    @staticmethod
    def __classcall_private__(cls, R, names=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = algebras.FreeDendriform(QQ, 'xyz')
            sage: F2 = algebras.FreeDendriform(QQ, ['x','y','z'])
            sage: F3 = algebras.FreeDendriform(QQ, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        if names is not None:
            if ',' in names:
                names = [u for u in names if u != ',']
            names = Alphabet(names)

        if R not in Rings():
            raise TypeError("argument R must be a ring")
        return super().__classcall__(cls, R, names)

    def __init__(self, R, names=None):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.FreeDendriform(QQ, '@'); A
            Free Dendriform algebra on one generator ['@'] over Rational Field
            sage: TestSuite(A).run()  # long time (3s)

            sage: F = algebras.FreeDendriform(QQ, 'xy')
            sage: TestSuite(F).run() # long time (3s)
        """
        if names is None:
            Trees = BinaryTrees()
            key = BinaryTree._sort_key
            self._alphabet = Alphabet(['o'])
        else:
            Trees = LabelledBinaryTrees()
            key = LabelledBinaryTree._sort_key
            self._alphabet = names
        # Here one would need LabelledBinaryTrees(names)
        # so that one can restrict the labels to some fixed set

        cat = HopfAlgebras(R).WithBasis().Graded().Connected()
        CombinatorialFreeModule.__init__(self, R, Trees,
                                         latex_prefix='',
                                         sorting_key=key,
                                         category=cat)

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R = algebras.FreeDendriform(QQ, 'xy')
            sage: R.variable_names()
            {'x', 'y'}
        """
        return self._alphabet

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.FreeDendriform(QQ, '@')  # indirect doctest
            Free Dendriform algebra on one generator ['@'] over Rational Field
        """
        n = self.algebra_generators().cardinality()
        finite = bool(n < Infinity)
        if not finite:
            gen = "generators indexed by"
        elif n == 1:
            gen = "one generator"
        else:
            gen = "{} generators".format(n)
        s = "Free Dendriform algebra on {} {} over {}"
        if finite:
            try:
                return s.format(gen, self._alphabet.list(), self.base_ring())
            except NotImplementedError:
                return s.format(gen, self._alphabet, self.base_ring())
        else:
            return s.format(gen, self._alphabet, self.base_ring())

    def gen(self, i):
        r"""
        Return the `i`-th generator of the algebra.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: F = algebras.FreeDendriform(ZZ, 'xyz')
            sage: F.gen(0)
            B[x[., .]]

            sage: F.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        G = self.algebra_generators()
        n = G.cardinality()
        if i < 0 or not i < n:
            m = "argument i (= {}) must be between 0 and {}".format(i, n - 1)
            raise IndexError(m)
        return G[G.keys().unrank(i)]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        These are the binary trees with just one vertex.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(ZZ, 'fgh'); A
            Free Dendriform algebra on 3 generators ['f', 'g', 'h']
             over Integer Ring
            sage: list(A.algebra_generators())
            [B[f[., .]], B[g[., .]], B[h[., .]]]

            sage: A = algebras.FreeDendriform(QQ, ['x1','x2'])
            sage: list(A.algebra_generators())
            [B[x1[., .]], B[x2[., .]]]
        """
        Trees = self.basis().keys()
        return Family(self._alphabet, lambda a: self.monomial(Trees([], a)))

    def change_ring(self, R):
        """
        Return the free dendriform algebra in the same variables over `R`.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: A = algebras.FreeDendriform(ZZ, 'fgh')
            sage: A.change_ring(QQ)
            Free Dendriform algebra on 3 generators ['f', 'g', 'h'] over
            Rational Field
        """
        return FreeDendriformAlgebra(R, names=self.variable_names())

    def gens(self) -> tuple:
        """
        Return the generators of ``self`` (as an algebra).

        EXAMPLES::

            sage: A = algebras.FreeDendriform(ZZ, 'fgh')
            sage: A.gens()
            (B[f[., .]], B[g[., .]], B[h[., .]])
        """
        return tuple(self.algebra_generators())

    def degree_on_basis(self, t):
        """
        Return the degree of a binary tree in the free Dendriform algebra.

        This is the number of vertices.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ,'@')
            sage: RT = A.basis().keys()
            sage: u = RT([], '@')
            sage: A.degree_on_basis(u.over(u))
            2
        """
        return t.node_number()

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, 'xy')
            sage: A.an_element()
            B[x[., .]] + 2*B[x[., x[., .]]] + 2*B[x[x[., .], .]]
        """
        o = self.gen(0)
        return o + 2 * o * o

    def some_elements(self):
        """
        Return some elements of the free dendriform algebra.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: A.some_elements()
            [B[.],
             B[[., .]],
             B[[., [., .]]] + B[[[., .], .]],
             B[.] + B[[., [., .]]] + B[[[., .], .]]]

        With several generators::

            sage: A = algebras.FreeDendriform(QQ, 'xy')
            sage: A.some_elements()
            [B[.],
             B[x[., .]],
             B[x[., x[., .]]] + B[x[x[., .], .]],
             B[.] + B[x[., x[., .]]] + B[x[x[., .], .]]]
        """
        u = self.one()
        o = self.gen(0)
        x = o * o
        y = u + x
        return [u, o, x, y]

    def one_basis(self):
        """
        Return the index of the unit.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: A.one_basis()
            .
            sage: A = algebras.FreeDendriform(QQ, 'xy')
            sage: A.one_basis()
            .
        """
        Trees = self.basis().keys()
        return Trees(None)

    def product_on_basis(self, x, y):
        r"""
        Return the `*` associative dendriform product of two trees.

        This is the sum over all possible ways of identifying the
        rightmost path in `x` and the leftmost path in `y`. Every term
        corresponds to a shuffle of the vertices on the rightmost path
        in `x` and the vertices on the leftmost path in `y`.

        .. SEEALSO::

            - :meth:`succ_product_on_basis`, :meth:`prec_product_on_basis`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: RT = A.basis().keys()
            sage: x = RT([])
            sage: A.product_on_basis(x, x)
            B[[., [., .]]] + B[[[., .], .]]
        """
        return self.sum(self.basis()[u] for u in x.dendriform_shuffle(y))

    def succ_product_on_basis(self, x, y):
        r"""
        Return the `\succ` dendriform product of two trees.

        This is the sum over all possible ways to identify the rightmost path
        in `x` and the leftmost path in `y`, with the additional condition
        that the root vertex of the result comes from `y`.

        The usual symbol for this operation is `\succ`.

        .. SEEALSO::

            - :meth:`product_on_basis`, :meth:`prec_product_on_basis`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: RT = A.basis().keys()
            sage: x = RT([])
            sage: A.succ_product_on_basis(x, x)
            B[[[., .], .]]

        TESTS::

            sage: u = A.one().support()[0]
            sage: A.succ_product_on_basis(u, u)
            Traceback (most recent call last):
            ...
            ValueError: dendriform products | < | and | > | are not defined
        """
        if y.is_empty():
            if x.is_empty():
                raise ValueError("dendriform products | < | and | > | are "
                                 "not defined")
            else:
                return []
        if x.is_empty():
            return [y]
        K = self.basis().keys()
        if hasattr(y, 'label'):
            return self.sum(self.basis()[K([u, y[1]], y.label())]
                            for u in x.dendriform_shuffle(y[0]))

        return self.sum(self.basis()[K([u, y[1]])]
                        for u in x.dendriform_shuffle(y[0]))

    @lazy_attribute
    def succ(self):
        r"""
        Return the `\succ` dendriform product.

        This is the sum over all possible ways of identifying the
        rightmost path in `x` and the leftmost path in `y`, with the
        additional condition that the root vertex of the result comes
        from `y`.

        The usual symbol for this operation is `\succ`.

        .. SEEALSO::

            :meth:`product`, :meth:`prec`, :meth:`over`, :meth:`under`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.succ(x, x)
            B[[[., .], .]]
        """
        suc = self.succ_product_on_basis
        return self._module_morphism(self._module_morphism(suc, position=0,
                                                           codomain=self),
                                     position=1)

    def prec_product_on_basis(self, x, y):
        r"""
        Return the `\prec` dendriform product of two trees.

        This is the sum over all possible ways of identifying the
        rightmost path in `x` and the leftmost path in `y`, with the
        additional condition that the root vertex of the result comes
        from `x`.

        The usual symbol for this operation is `\prec`.

        .. SEEALSO::

            - :meth:`product_on_basis`, :meth:`succ_product_on_basis`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: RT = A.basis().keys()
            sage: x = RT([])
            sage: A.prec_product_on_basis(x, x)
            B[[., [., .]]]

        TESTS::

            sage: u = A.one().support()[0]
            sage: A.prec_product_on_basis(u, u)
            Traceback (most recent call last):
            ...
            ValueError: dendriform products | < | and | > | are not defined
        """
        if x.is_empty() and y.is_empty():
            raise ValueError("dendriform products | < | and | > | are "
                             "not defined")
        if x.is_empty():
            return []
        if y.is_empty():
            return [x]
        K = self.basis().keys()
        if hasattr(y, 'label'):
            return self.sum(self.basis()[K([x[0], u], x.label())]
                            for u in x[1].dendriform_shuffle(y))

        return self.sum(self.basis()[K([x[0], u])]
                        for u in x[1].dendriform_shuffle(y))

    @lazy_attribute
    def prec(self):
        r"""
        Return the `\prec` dendriform product.

        This is the sum over all possible ways to identify the rightmost path
        in `x` and the leftmost path in `y`, with the additional condition
        that the root vertex of the result comes from `x`.

        The usual symbol for this operation is `\prec`.

        .. SEEALSO::

            :meth:`product`, :meth:`succ`, :meth:`over`, :meth:`under`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.prec(x, x)
            B[[., [., .]]]
        """
        pre = self.prec_product_on_basis
        return self._module_morphism(self._module_morphism(pre, position=0,
                                                           codomain=self),
                                     position=1)

    @lazy_attribute
    def over(self):
        r"""
        Return the over product.

        The over product `x/y` is the binary tree obtained by
        grafting the root of `y` at the rightmost leaf of `x`.

        The usual symbol for this operation is `/`.

        .. SEEALSO::

            :meth:`product`, :meth:`succ`, :meth:`prec`, :meth:`under`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.over(x, x)
            B[[., [., .]]]
        """
        def ov(x, y):
            return self._monomial(x.over(y))
        return self._module_morphism(self._module_morphism(ov, position=0,
                                                           codomain=self),
                                     position=1)

    @lazy_attribute
    def under(self):
        r"""
        Return the under product.

        The over product `x \backslash y` is the binary tree obtained by
        grafting the root of `x` at the leftmost leaf of `y`.

        The usual symbol for this operation is `\backslash`.

        .. SEEALSO::

            :meth:`product`, :meth:`succ`, :meth:`prec`, :meth:`over`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.under(x, x)
            B[[[., .], .]]
        """
        def und(x, y):
            return self._monomial(x.under(y))
        return self._module_morphism(self._module_morphism(und, position=0,
                                                           codomain=self),
                                     position=1)

    def coproduct_on_basis(self, x):
        """
        Return the coproduct of a binary tree.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ)
            sage: x = A.gen(0)
            sage: ascii_art(A.coproduct(A.one()))  # indirect doctest
            1 # 1

            sage: ascii_art(A.coproduct(x))  # indirect doctest
            1 # B  + B  # 1
                 o    o

            sage: A = algebras.FreeDendriform(QQ, 'xyz')
            sage: x, y, z = A.gens()
            sage: w = A.under(z,A.over(x,y))
            sage: A.coproduct(z)
            B[.] # B[z[., .]] + B[z[., .]] # B[.]
            sage: A.coproduct(w)
            B[.] # B[x[z[., .], y[., .]]] + B[x[., .]] # B[z[., y[., .]]] +
            B[x[., .]] # B[y[z[., .], .]] + B[x[., y[., .]]] # B[z[., .]] +
            B[x[z[., .], .]] # B[y[., .]] + B[x[z[., .], y[., .]]] # B[.]
        """
        B = self.basis()
        Trees = B.keys()
        if not x.node_number():
            return self.one().tensor(self.one())
        L, R = list(x)
        try:
            root = x.label()
        except AttributeError:
            root = '@'
        resu = self.one().tensor(self.monomial(x))
        resu += sum(cL * cR *
                    self.monomial(Trees([LL[0], RR[0]], root)).tensor(
                        self.monomial(LL[1]) * self.monomial(RR[1]))
                    for LL, cL in self.coproduct_on_basis(L)
                    for RR, cR in self.coproduct_on_basis(R))
        return resu

    # after this line : coercion
    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = algebras.FreeDendriform(QQ, 'xy')
            sage: x, y = R.gens()
            sage: R(x)
            B[x[., .]]
            sage: R(x+4*y)
            B[x[., .]] + 4*B[y[., .]]

            sage: Trees = R.basis().keys()
            sage: R(Trees([],'x'))
            B[x[., .]]
            sage: D = algebras.FreeDendriform(ZZ, 'xy')
            sage: X, Y = D.gens()
            sage: R(X-Y).parent()
            Free Dendriform algebra on 2 generators ['x', 'y'] over Rational Field
        """
        if x in self.basis().keys():
            return self.monomial(x)
        try:
            P = x.parent()
            if isinstance(P, FreeDendriformAlgebra):
                if P is self:
                    return x
                return self.element_class(self, x.monomial_coefficients())
        except AttributeError:
            raise TypeError('not able to coerce this in this algebra')
        # Ok, not a dendriform algebra element (or should not be viewed as one).

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - free dendriform algebras in a subset of variables of ``self``
          over a base with a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: F = algebras.FreeDendriform(GF(7), 'xyz'); F
            Free Dendriform algebra on 3 generators ['x', 'y', 'z']
             over Finite Field of size 7

        Elements of the free dendriform algebra canonically coerce in::

            sage: x, y, z = F.gens()
            sage: F.coerce(x+y) == x+y
            True

        The free dendriform algebra over `\ZZ` on `x, y, z` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = algebras.FreeDendriform(ZZ, 'xyz')
            sage: Gx,Gy,Gz = G.gens()
            sage: z = F.coerce(Gx+Gy); z
            B[x[., .]] + B[y[., .]]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the free dendriform
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Dendriform algebra
             on 3 generators ['x', 'y', 'z'] over Finite Field of size
             7 to Free Dendriform algebra on 3 generators ['x', 'y', 'z']
             over Integer Ring

        TESTS::

            sage: F = algebras.FreeDendriform(ZZ, 'xyz')
            sage: G = algebras.FreeDendriform(QQ, 'xyz')
            sage: H = algebras.FreeDendriform(ZZ, 'y')
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            True
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            False
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # free dendriform algebras in a subset of variables
        # over any base that coerces in:
        if isinstance(R, FreeDendriformAlgebra):
            if all(x in self.variable_names() for x in R.variable_names()):
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
        return False

    def construction(self):
        """
        Return a pair ``(F, R)``, where ``F`` is a :class:`DendriformFunctor`
        and `R` is a ring, such that ``F(R)`` returns ``self``.

        EXAMPLES::

            sage: P = algebras.FreeDendriform(ZZ, 'x,y')
            sage: x,y = P.gens()
            sage: F, R = P.construction()
            sage: F
            Dendriform[x,y]
            sage: R
            Integer Ring
            sage: F(ZZ) is P
            True
            sage: F(QQ)
            Free Dendriform algebra on 2 generators ['x', 'y'] over Rational Field
        """
        return DendriformFunctor(self.variable_names()), self.base_ring()


class DendriformFunctor(ConstructionFunctor):
    """
    A constructor for dendriform algebras.

    EXAMPLES::

        sage: P = algebras.FreeDendriform(ZZ, 'x,y')
        sage: x,y = P.gens()
        sage: F = P.construction()[0]; F
        Dendriform[x,y]

        sage: A = GF(5)['a,b']
        sage: a, b = A.gens()
        sage: F(A)
        Free Dendriform algebra on 2 generators ['x', 'y']
         over Multivariate Polynomial Ring in a, b over Finite Field of size 5

        sage: f = A.hom([a+b,a-b],A)
        sage: F(f)
        Generic endomorphism of Free Dendriform algebra on 2 generators ['x', 'y']
         over Multivariate Polynomial Ring in a, b over Finite Field of size 5

        sage: F(f)(a * F(A)(x))
        (a+b)*B[x[., .]]
    """
    rank = 9

    def __init__(self, vars):
        """
        EXAMPLES::

            sage: F = sage.combinat.free_dendriform_algebra.DendriformFunctor(['x','y'])
            sage: F
            Dendriform[x,y]
            sage: F(ZZ)
            Free Dendriform algebra on 2 generators ['x', 'y']  over Integer Ring
        """
        Functor.__init__(self, Rings(), Rings())
        self.vars = vars

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        EXAMPLES::

            sage: R = algebras.FreeDendriform(ZZ, 'x,y,z')
            sage: F = R.construction()[0]; F
            Dendriform[x,y,z]
            sage: type(F)
            <class 'sage.combinat.free_dendriform_algebra.DendriformFunctor'>
            sage: F(ZZ)          # indirect doctest
            Free Dendriform algebra on 3 generators ['x', 'y', 'z'] over Integer Ring
        """
        return FreeDendriformAlgebra(R, self.vars)

    def _apply_functor_to_morphism(self, f):
        """
        Apply the functor ``self`` to the ring morphism `f`.

        TESTS::

            sage: R = algebras.FreeDendriform(ZZ, 'x').construction()[0]
            sage: R(ZZ.hom(GF(3)))  # indirect doctest
            Generic morphism:
              From: Free Dendriform algebra on one generator ['x'] over Integer Ring
              To:   Free Dendriform algebra on one generator ['x'] over Finite Field of size 3
        """
        dom = self(f.domain())
        codom = self(f.codomain())

        def action(x):
            return codom._from_dict({a: f(b)
                                     for a, b in
                                     x.monomial_coefficients().items()})
        return dom.module_morphism(function=action, codomain=codom)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: F = algebras.FreeDendriform(ZZ, 'x,y,z').construction()[0]
            sage: G = algebras.FreeDendriform(QQ, 'x,y,z').construction()[0]
            sage: F == G
            True
            sage: G == loads(dumps(G))
            True
            sage: G = algebras.FreeDendriform(QQ, 'x,y').construction()[0]
            sage: F == G
            False
        """
        if not isinstance(other, DendriformFunctor):
            return False
        return self.vars == other.vars

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: F = algebras.FreeDendriform(ZZ, 'x,y,z').construction()[0]
            sage: G = algebras.FreeDendriform(QQ, 'x,y,z').construction()[0]
            sage: F != G
            False
            sage: G != loads(dumps(G))
            False
            sage: G = algebras.FreeDendriform(QQ, 'x,y').construction()[0]
            sage: F != G
            True
        """
        return not (self == other)

    def __mul__(self, other):
        """
        If two Dendriform functors are given in a row, form a single Dendriform functor
        with all of the variables.

        EXAMPLES::

            sage: F = sage.combinat.free_dendriform_algebra.DendriformFunctor(['x','y'])
            sage: G = sage.combinat.free_dendriform_algebra.DendriformFunctor(['t'])
            sage: G * F
            Dendriform[x,y,t]
        """
        if isinstance(other, IdentityConstructionFunctor):
            return self
        if isinstance(other, DendriformFunctor):
            if set(self.vars).intersection(other.vars):
                raise CoercionException("Overlapping variables (%s,%s)" %
                                        (self.vars, other.vars))
            return DendriformFunctor(other.vars + self.vars)
        elif (isinstance(other, CompositeConstructionFunctor) and
              isinstance(other.all[-1], DendriformFunctor)):
            return CompositeConstructionFunctor(other.all[:-1],
                                                self * other.all[-1])
        else:
            return CompositeConstructionFunctor(other, self)

    def merge(self, other):
        """
        Merge ``self`` with another construction functor, or return ``None``.

        EXAMPLES::

            sage: F = sage.combinat.free_dendriform_algebra.DendriformFunctor(['x','y'])
            sage: G = sage.combinat.free_dendriform_algebra.DendriformFunctor(['t'])
            sage: F.merge(G)
            Dendriform[x,y,t]
            sage: F.merge(F)
            Dendriform[x,y]

        Now some actual use cases::

            sage: R = algebras.FreeDendriform(ZZ, 'x,y,z')
            sage: x,y,z = R.gens()
            sage: 1/2 * x
            1/2*B[x[., .]]
            sage: parent(1/2 * x)
            Free Dendriform algebra on 3 generators ['x', 'y', 'z'] over Rational Field

            sage: S = algebras.FreeDendriform(QQ, 'zt')
            sage: z,t = S.gens()
            sage: x + t
            B[t[., .]] + B[x[., .]]
            sage: parent(x + t)
            Free Dendriform algebra on 4 generators ['z', 't', 'x', 'y'] over Rational Field
        """
        if isinstance(other, DendriformFunctor):
            if self.vars == other.vars:
                return self
            ret = list(self.vars)
            cur_vars = set(ret)
            ret.extend(v for v in other.vars if v not in cur_vars)
            return DendriformFunctor(Alphabet(ret))

        return None

    def _repr_(self) -> str:
        """
        TESTS::

            sage: algebras.FreeDendriform(QQ,'x,y,z,t').construction()[0]
            Dendriform[x,y,z,t]
        """
        return "Dendriform[%s]" % ','.join(self.vars)
