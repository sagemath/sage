# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
"""
Free Zinbiel Algebras

AUTHORS:

- Travis Scrimshaw (2015-09): initial version
"""

# ****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.categories.magmatic_algebras import MagmaticAlgebras
from sage.categories.magmas import Magmas
from sage.categories.pushout import (ConstructionFunctor,
                                     CompositeConstructionFunctor,
                                     IdentityConstructionFunctor)
from sage.categories.coalgebras_with_basis import CoalgebrasWithBasis
from sage.categories.rings import Rings
from sage.categories.functor import Functor
from sage.categories.sets_cat import Sets
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
from sage.combinat.words.alphabet import Alphabet
from sage.sets.family import Family
from sage.structure.coerce_exceptions import CoercionException


class FreeZinbielAlgebra(CombinatorialFreeModule):
    r"""
    The free Zinbiel algebra on `n` generators.

    Let `R` be a ring. A *Zinbiel algebra* is a non-associative
    algebra with multiplication `\circ` that satisfies

    .. MATH::

        (a \circ b) \circ c = a \circ (b \circ c) + a \circ (c \circ b).

    Zinbiel algebras were first introduced by Loday (see [Lod1995]_ and
    [LV2012]_) as the Koszul dual to Leibniz algebras (hence the name
    coined by Lemaire).

    By default, the convention above is used. The opposite product,
    which satisfy the opposite axiom, can be used instead by setting
    the ``side`` parameter to ``'>'`` instead of the default value ``'<'``.

    Zinbiel algebras are divided power algebras, in that for

    .. MATH::

        x^{\circ n} = \bigl(x \circ (x \circ \cdots \circ( x \circ x) \cdots
        ) \bigr)

    we have

    .. MATH::

        x^{\circ m} \circ x^{\circ n} = \binom{n+m-1}{m} x^{n+m}

    and

    .. MATH::

        \underbrace{\bigl( ( x \circ \cdots \circ x \circ (x \circ x) \cdots
        ) \bigr)}_{n+1 \text{ times}} = n! x^n.

    .. NOTE::

        This implies that Zinbiel algebras are not power associative.

    To every Zinbiel algebra, we can construct a corresponding commutative
    associative algebra by using the symmetrized product:

    .. MATH::

        a * b = a \circ b + b \circ a.

    The free Zinbiel algebra on `n` generators is isomorphic as `R`-modules
    to the reduced tensor algebra `\bar{T}(R^n)` with the product

    .. MATH::

        (x_0 x_1 \cdots x_p) \circ (x_{p+1} x_{p+2} \cdots x_{p+q})
        = \sum_{\sigma \in S_{p,q}} x_0 (x_{\sigma(1)} x_{\sigma(2)}
        \cdots x_{\sigma(p+q)},

    where `S_{p,q}` is the set of `(p,q)`-shuffles.

    The free Zinbiel algebra is free as a divided power algebra. Moreover,
    the corresponding commutative algebra is isomorphic to the (non-unital)
    shuffle algebra.

    INPUT:

    - ``R`` -- a ring
    - ``n`` -- (optional) the number of generators
    - ``names`` -- the generator names

    .. WARNING::

        Currently the basis is indexed by all finite words over the variables,
        including the empty word. This is a slight abuse as it is supposed
        to be indexed by all non-empty words.

    EXAMPLES:

    We create the free Zinbiel algebra and check the defining relation::

        sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
        sage: (x*y)*z
        Z[xyz] + Z[xzy]
        sage: x*(y*z) + x*(z*y)
        Z[xyz] + Z[xzy]

    We see that the Zinbiel algebra is not associative, not even
    power associative::

        sage: x*(y*z)
        Z[xyz]
        sage: x*(x*x)
        Z[xxx]
        sage: (x*x)*x
        2*Z[xxx]

    We verify that it is a divided power algebra::

        sage: (x*(x*x)) * (x*(x*(x*x)))
        15*Z[xxxxxxx]
        sage: binomial(3+4-1,4)
        15
        sage: (x*(x*(x*x))) * (x*(x*x))
        20*Z[xxxxxxx]
        sage: binomial(3+4-1,3)
        20
        sage: ((x*x)*x)*x
        6*Z[xxxx]
        sage: (((x*x)*x)*x)*x
        24*Z[xxxxx]

    A few tests with the opposite convention for the product::

        sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ, side='>')
        sage: (x*y)*z
        Z[xyz]
        sage: x*(y*z)
        Z[xyz] + Z[yxz]

    TESTS::

        sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
        sage: Z.basis().keys()
        Finite words over {'x', 'y', 'z'}

        sage: A = algebras.FreeZinbiel(QQ,'z2,z3')
        sage: x, y = A.gens()
        sage: x*y
        Z[z2,z3]

    REFERENCES:

    - :wikipedia:`Zinbiel_algebra`

    - [Lod1995]_

    - [LV2012]_
    """
    @staticmethod
    def __classcall_private__(cls, R, n=None, names=None,
                              prefix=None, side=None):
        """
        Standardize input to ensure a unique representation.

        TESTS::

            sage: Z1.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: Z2.<x,y,z> = algebras.FreeZinbiel(QQ, 3)
            sage: Z3 = algebras.FreeZinbiel(QQ, 3, 'x,y,z')
            sage: Z4.<x,y,z> = algebras.FreeZinbiel(QQ, 'x,y,z')
            sage: Z1 is Z2 and Z1 is Z3 and Z1 is Z4
            True

            sage: algebras.FreeZinbiel(QQ, ['x', 'y'])
            Free Zinbiel algebra on generators (Z[x], Z[y]) over Rational Field
            sage: algebras.FreeZinbiel(QQ, ('x', 'y'))
            Free Zinbiel algebra on generators (Z[x], Z[y]) over Rational Field

            sage: Z = algebras.FreeZinbiel(QQ, ZZ)
        """
        if isinstance(n, (list, tuple)):
            names = n
            n = len(names)
        elif isinstance(n, str):
            names = n.split(',')
            n = len(names)
        elif isinstance(names, str):
            names = names.split(',')
        elif n is None:
            n = len(names)
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        if prefix is None:
            prefix = 'Z'
        if side is None:
            side = '<'
        if side not in ['<', '>']:
            raise ValueError("side must be either '<' or '>'")
        if names is None:
            return super().__classcall__(cls, R, n, None, prefix, side)
        return super().__classcall__(cls, R, n, tuple(names), prefix, side)

    def __init__(self, R, n, names, prefix, side):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: TestSuite(Z).run()

            sage: Z = algebras.FreeZinbiel(QQ, ZZ)
            sage: G = Z.algebra_generators()
            sage: TestSuite(Z).run(elements=[Z.an_element(), G[1], G[1]*G[2]*G[0]])

        TESTS::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(5)
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring

            sage: algebras.FreeZinbiel(QQ, ['x', 'y'], prefix='f')
            Free Zinbiel algebra on generators (f[x], f[y]) over Rational Field
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        if names is None:
            indices = Words(Alphabet(n), infinite=False)
            self._n = None
        else:
            indices = Words(Alphabet(n, names=names), infinite=False)
            self._n = n
        self._side = side
        if side == '<':
            self.product_on_basis = self.product_on_basis_left
        else:
            self.product_on_basis = self.product_on_basis_right
        cat = MagmaticAlgebras(R).WithBasis().Graded()
        cat &= CoalgebrasWithBasis(R)
        CombinatorialFreeModule.__init__(self, R, indices, prefix=prefix,
                                         category=cat)
        if self._n is not None:
            self._assign_names(names)

    def _repr_term(self, t):
        """
        Return a string representation of the basis element indexed by ``t``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: Z._repr_term(Z._indices('xyzxxy'))
            'Z[xyzxxy]'
        """
        return "{!s}[{!s}]".format(self._print_options['prefix'], repr(t)[6:])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Z.<x,y> = algebras.FreeZinbiel(QQ)
            sage: Z
            Free Zinbiel algebra on generators (Z[x], Z[y]) over Rational Field

            sage: Z = algebras.FreeZinbiel(QQ, ZZ)
            sage: Z
            Free Zinbiel algebra on generators indexed by Integer Ring over Rational Field
        """
        if self._n is None:
            return "Free Zinbiel algebra on generators indexed by {} over {}".format(
                self._indices.alphabet(), self.base_ring())
        return "Free Zinbiel algebra on generators {} over {}".format(
            self.gens(), self.base_ring())

    def side(self):
        """
        Return the choice of side for the product.

        This is either ``'<'`` or ``'>'``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: Z.side()
            '<'
        """
        return self._side

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: list(Z.algebra_generators())
            [Z[x], Z[y], Z[z]]
        """
        if self._n is None:
            A = self._indices.alphabet()
        else:
            A = self.variable_names()
        return Family(A, lambda g: self.monomial(self._indices([g])))

    def change_ring(self, R):
        """
        Return the free Zinbiel algebra in the same variables over ``R``.

        INPUT:

        - ``R`` -- a ring

        The same side convention is used for the product.

        EXAMPLES::

            sage: A = algebras.FreeZinbiel(ZZ, 'f,g,h')
            sage: A.change_ring(QQ)
            Free Zinbiel algebra on generators (Z[f], Z[g], Z[h])
            over Rational Field
        """
        A = self.variable_names()
        return FreeZinbielAlgebra(R, n=len(A), names=A, side=self._side)

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: Z.gens()
            (Z[x], Z[y], Z[z])
        """
        if self._n is None:
            return self.algebra_generators()
        return tuple(self.algebra_generators())

    def degree_on_basis(self, t):
        """
        Return the degree of a word in the free Zinbiel algebra.

        This is the length.

        EXAMPLES::

            sage: A = algebras.FreeZinbiel(QQ, 'x,y')
            sage: W = A.basis().keys()
            sage: A.degree_on_basis(W('xy'))
            2
        """
        return len(t)

    def product_on_basis_left(self, x, y):
        """
        Return the product < of the basis elements indexed by ``x`` and ``y``.

        This is one half of the shuffle product, where the first letter
        comes from the first letter of the first argument.

        INPUT:

        - ``x``, ``y`` -- two words

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: (x*y)*z  # indirect doctest
            Z[xyz] + Z[xzy]

        TESTS::

            sage: Z.<x,y> = algebras.FreeZinbiel(QQ)
            sage: Z.product_on_basis(Word(), Word('y'))
            Z[y]
        """
        if not x:
            return self.monomial(y)
        x0 = self._indices([x[0]])
        return self.sum_of_monomials(x0 + sh for sh in x[1:].shuffle(y))

    def product_on_basis_right(self, x, y):
        """
        Return the product > of the basis elements indexed by ``x`` and ``y``.

        This is one half of the shuffle product, where the last letter
        comes from the last letter of the second argument.

        INPUT:

        - ``x``, ``y`` -- two words

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ, side='>')
            sage: (x*y)*z  # indirect doctest
            Z[xyz]

        TESTS::

            sage: Z.<x,y> = algebras.FreeZinbiel(QQ, side='>')
            sage: Z.product_on_basis(Word('x'), Word())
            Z[x]
        """
        if not y:
            return self.monomial(x)
        yf = self._indices([y[-1]])
        return self.sum_of_monomials(sh + yf for sh in x.shuffle(y[:-1]))

    def coproduct_on_basis(self, w):
        """
        Return the coproduct of the element of the basis indexed by
        the word ``w``.

        The coproduct is given by deconcatenation.

        INPUT:

        - ``w`` -- a word

        EXAMPLES::

            sage: F = algebras.FreeZinbiel(QQ,['a','b'])
            sage: F.coproduct_on_basis(Word('a'))
            Z[] # Z[a] + Z[a] # Z[]
            sage: F.coproduct_on_basis(Word('aba'))
            Z[] # Z[aba] + Z[a] # Z[ba] + Z[ab] # Z[a] + Z[aba] # Z[]
            sage: F.coproduct_on_basis(Word())
            Z[] # Z[]

        TESTS::

            sage: F = algebras.FreeZinbiel(QQ,['a','b'])
            sage: S = F.an_element(); S
            Z[] + 2*Z[a] + 3*Z[b] + Z[bab]
            sage: F.coproduct(S)
            Z[] # Z[] + 2*Z[] # Z[a] + 3*Z[] # Z[b] + Z[] # Z[bab] +
            2*Z[a] # Z[] + 3*Z[b] # Z[] + Z[b] # Z[ab] + Z[ba] # Z[b] +
            Z[bab] # Z[]
        """
        TS = self.tensor_square()
        return TS.sum_of_monomials((w[:i], w[i:]) for i in range(len(w) + 1))

    def counit(self, S):
        """
        Return the counit of ``S``.

        EXAMPLES::

            sage: F = algebras.FreeZinbiel(QQ,['a','b'])
            sage: S = F.an_element(); S
            Z[] + 2*Z[a] + 3*Z[b] + Z[bab]
            sage: F.counit(S)
            1
        """
        W = self.basis().keys()
        return S.coefficient(W())

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = algebras.FreeZinbiel(QQ, 'x,y')
            sage: x, y = R.gens()
            sage: R(x)
            Z[x]
            sage: R(x+4*y)
            Z[x] + 4*Z[y]

            sage: W = R.basis().keys()
            sage: R(W('x'))
            Z[x]
            sage: D = algebras.FreeZinbiel(ZZ, 'x,y')
            sage: X, Y = D.gens()
            sage: R(X-Y).parent()
            Free Zinbiel algebra on generators (Z[x], Z[y]) over Rational Field

        TESTS::

            sage: R.<x,y> = algebras.FreeZinbiel(QQ)
            sage: S.<z> = algebras.FreeZinbiel(GF(3))
            sage: R(z)
            Traceback (most recent call last):
            ...
            TypeError: not able to convert this to this algebra
        """
        if x in self.basis().keys():
            return self.monomial(x)
        try:
            P = x.parent()
        except AttributeError:
            raise TypeError('not able to convert this to this algebra')
        if isinstance(P, FreeZinbielAlgebra) and self._coerce_map_from_(P):
            if self._side == P._side:
                return self.element_class(self,
                                          x.monomial_coefficients(copy=False))
            else:
                dic = x.monomial_coefficients(copy=False)
                # canonical isomorphism when switching side
                return self.element_class(self,
                                          {w.reversal(): cf
                                           for w, cf in dic.items()})
        else:
            raise TypeError('not able to convert this to this algebra')
        # Ok, not a Zinbiel algebra element (or should not be viewed as one).

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - free Zinbiel algebras whose set `E` of labels is
          a subset of the corresponding ``self`` of ``set``, and whose base
          ring has a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: F = algebras.FreeZinbiel(GF(7), 'x,y,z'); F
            Free Zinbiel algebra on generators (Z[x], Z[y], Z[z])
            over Finite Field of size 7

        Elements of the free Zinbiel algebra canonically coerce in::

            sage: x, y, z = F.gens()
            sage: F.coerce(x+y) == x+y
            True

        The free Zinbiel algebra over `\ZZ` on `x, y, z` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = algebras.FreeZinbiel(ZZ, 'x,y,z')
            sage: Gx,Gy,Gz = G.gens()
            sage: z = F.coerce(Gx+Gy); z
            Z[x] + Z[y]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the free Zinbiel
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Zinbiel algebra on
            generators (Z[x], Z[y], Z[z]) over Finite Field of size 7 to
            Free Zinbiel algebra on generators (Z[x], Z[y], Z[z])
            over Integer Ring

        TESTS::

            sage: F = algebras.FreeZinbiel(ZZ, 'x,y,z')
            sage: G = algebras.FreeZinbiel(QQ, 'x,y,z')
            sage: H = algebras.FreeZinbiel(ZZ, 'y')
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            True
            sage: F._coerce_map_from_(QQ) is None
            True
            sage: G._coerce_map_from_(QQ) is None
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False

            sage: I = algebras.FreeZinbiel(ZZ, ZZ)
            sage: F._coerce_map_from_(I)
            False
            sage: I._coerce_map_from_(F)
            False
        """
        # free Zinbiel algebras in a subset of variables
        # over any base that coerces in:
        if isinstance(R, FreeZinbielAlgebra):
            if self._n is None or R._n is None:
                return False
            return (all(x in self.variable_names()
                        for x in R.variable_names()) and
                    self.base_ring().has_coerce_map_from(R.base_ring()))
        return super()._coerce_map_from_(R)

    def construction(self):
        """
        Return a pair ``(F, R)``, where ``F`` is a :class:`ZinbielFunctor`
        and ``R`` is a ring, such that ``F(R)`` returns ``self``.

        EXAMPLES::

            sage: P = algebras.FreeZinbiel(ZZ, 'x,y')
            sage: x,y = P.gens()
            sage: F, R = P.construction()
            sage: F
            Zinbiel[x,y]
            sage: R
            Integer Ring
            sage: F(ZZ) is P
            True
            sage: F(QQ)
            Free Zinbiel algebra on generators (Z[x], Z[y]) over Rational Field
        """
        if self._n is None:
            A = self._indices.alphabet()
        else:
            A = self.variable_names()
        return ZinbielFunctor(A, side=self._side), self.base_ring()


class ZinbielFunctor(ConstructionFunctor):
    """
    A constructor for free Zinbiel algebras.

    EXAMPLES::

        sage: P = algebras.FreeZinbiel(ZZ, 'x,y')
        sage: x,y = P.gens()
        sage: F = P.construction()[0]; F
        Zinbiel[x,y]

        sage: A = GF(5)['a,b']
        sage: a, b = A.gens()
        sage: F(A)
        Free Zinbiel algebra on generators (Z[x], Z[y])
        over Multivariate Polynomial Ring in a, b over Finite Field of size 5

        sage: f = A.hom([a+b,a-b],A)
        sage: F(f)
        Generic endomorphism of Free Zinbiel algebra on generators (Z[x], Z[y])
        over Multivariate Polynomial Ring in a, b over Finite Field of size 5

        sage: F(f)(a * F(A)(x))
        (a+b)*Z[x]
    """
    rank = 9

    def __init__(self, variables, side):
        """
        EXAMPLES::

            sage: functor = sage.algebras.free_zinbiel_algebra.ZinbielFunctor
            sage: F = functor(['x','y'], '<'); F
            Zinbiel[x,y]
            sage: F(ZZ)
            Free Zinbiel algebra on generators (Z[x], Z[y]) over Integer Ring
        """
        Functor.__init__(self, Rings(), Magmas())
        self.vars = variables
        self._side = side
        self._finite_vars = bool(isinstance(variables, (list, tuple)) or variables in Sets().Finite())

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        EXAMPLES::

            sage: R = algebras.FreeZinbiel(ZZ, 'x,y,z')
            sage: F = R.construction()[0]; F
            Zinbiel[x,y,z]
            sage: type(F)
            <class 'sage.algebras.free_zinbiel_algebra.ZinbielFunctor'>
            sage: F(ZZ)          # indirect doctest
            Free Zinbiel algebra on generators (Z[x], Z[y], Z[z])
             over Integer Ring

            sage: R = algebras.FreeZinbiel(QQ, ZZ)
            sage: F = R.construction()[0]; F
            Zinbiel[Integer Ring]
            sage: F(ZZ)  # indirect doctest
            Free Zinbiel algebra on generators indexed by Integer Ring over Integer Ring
        """
        if self._finite_vars:
            return FreeZinbielAlgebra(R, len(self.vars), self.vars,
                                      side=self._side)
        return FreeZinbielAlgebra(R, self.vars, side=self._side)

    def _apply_functor_to_morphism(self, f):
        """
        Apply the functor ``self`` to the ring morphism `f`.

        TESTS::

            sage: R = algebras.FreeZinbiel(ZZ, 'x').construction()[0]
            sage: R(ZZ.hom(GF(3)))  # indirect doctest
            Generic morphism:
              From: Free Zinbiel algebra on generators (Z[x],)
                    over Integer Ring
              To:   Free Zinbiel algebra on generators (Z[x],)
                    over Finite Field of size 3
        """
        dom = self(f.domain())
        codom = self(f.codomain())

        def action(x):
            return codom._from_dict({a: f(b)
                                     for a, b in x.monomial_coefficients(copy=False).items()})
        return dom.module_morphism(function=action, codomain=codom)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: F = algebras.FreeZinbiel(ZZ, 'x,y,z').construction()[0]
            sage: G = algebras.FreeZinbiel(QQ, 'x,y,z').construction()[0]
            sage: F == G
            True
            sage: G == loads(dumps(G))
            True
            sage: G = algebras.FreeZinbiel(QQ, 'x,y').construction()[0]
            sage: F == G
            False
        """
        if not isinstance(other, ZinbielFunctor):
            return False
        return self.vars == other.vars and self._side == other._side

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: F = algebras.FreeZinbiel(ZZ, 'x,y,z').construction()[0]
            sage: G = algebras.FreeZinbiel(QQ, 'x,y,z').construction()[0]
            sage: hash(F) == hash(G)
            True
        """
        return hash(repr(self))

    def __mul__(self, other):
        """
        If two Zinbiel functors are given in a row, form a single
        Zinbiel functor with all of the variables.

        EXAMPLES::

            sage: from sage.algebras.free_zinbiel_algebra import ZinbielFunctor as functor
            sage: F = functor(['x','y'], '<')
            sage: G = functor(['t'], '<')
            sage: G * F
            Zinbiel[x,y,t]

        With an infinite generating set::

            sage: H = functor(ZZ, '<')
            sage: H * G
            Traceback (most recent call last):
            ...
            CoercionException: Unable to determine overlap for infinite sets
            sage: G * H
            Traceback (most recent call last):
            ...
            CoercionException: Unable to determine overlap for infinite sets
        """
        if isinstance(other, IdentityConstructionFunctor):
            return self
        if isinstance(other, ZinbielFunctor):
            if self._side != other._side:
                raise CoercionException("detection of distinct sides")
            if not self._finite_vars or not other._finite_vars:
                raise CoercionException("Unable to determine overlap for infinite sets")
            if set(self.vars).intersection(other.vars):
                raise CoercionException("Overlapping variables (%s,%s)" %
                                        (self.vars, other.vars))
            return ZinbielFunctor(other.vars + self.vars, self._side)
        elif (isinstance(other, CompositeConstructionFunctor) and
              isinstance(other.all[-1], ZinbielFunctor)):
            return CompositeConstructionFunctor(other.all[:-1],
                                                self * other.all[-1])
        else:
            return CompositeConstructionFunctor(other, self)

    def merge(self, other):
        """
        Merge ``self`` with another construction functor, or return ``None``.

        EXAMPLES::

            sage: functor = sage.algebras.free_zinbiel_algebra.ZinbielFunctor
            sage: F = functor(['x','y'], '<')
            sage: G = functor(['t'], '<')
            sage: F.merge(G)
            Zinbiel[x,y,t]
            sage: F.merge(F)
            Zinbiel[x,y]

        With an infinite generating set::

            sage: H = functor(ZZ, '<')
            sage: H.merge(H) is H
            True
            sage: H.merge(F) is None
            True
            sage: F.merge(H) is None
            True

        Now some actual use cases::

            sage: R = algebras.FreeZinbiel(ZZ, 'x,y,z')
            sage: x,y,z = R.gens()
            sage: 1/2 * x
            1/2*Z[x]
            sage: parent(1/2 * x)
            Free Zinbiel algebra on generators (Z[x], Z[y], Z[z])
            over Rational Field

            sage: S = algebras.FreeZinbiel(QQ, 'z,t')
            sage: z,t = S.gens()
            sage: x * t
            Z[xt]
            sage: parent(x * t)
            Free Zinbiel algebra on generators (Z[z], Z[t], Z[x], Z[y])
            over Rational Field

        TESTS:

        Using the other side convention::

            sage: F = functor(['x','y'], '>')
            sage: G = functor(['t'], '>')
            sage: H = functor(['t'], '<')
            sage: F.merge(G)
            Zinbiel[x,y,t]
            sage: F.merge(H)
            Traceback (most recent call last):
            ...
            TypeError: cannot merge free Zinbiel algebras with distinct sides
        """
        if isinstance(other, ZinbielFunctor):
            if self._side != other._side:
                raise TypeError('cannot merge free Zinbiel algebras '
                                'with distinct sides')
            if self.vars == other.vars:
                return self

            def check(x):
                return isinstance(x, (list, tuple)) or x in Sets().Finite()
            if not check(self.vars) or not check(other.vars):
                return None
            ret = list(self.vars)
            cur_vars = set(ret)
            ret.extend(v for v in other.vars if v not in cur_vars)
            return ZinbielFunctor(ret, self._side)

        return None

    def _repr_(self) -> str:
        """
        TESTS::

            sage: algebras.FreeZinbiel(QQ,'x,y,z,t').construction()[0]
            Zinbiel[x,y,z,t]

            sage: algebras.FreeZinbiel(QQ, ZZ).construction()[0]
            Zinbiel[Integer Ring]
        """
        if self._finite_vars:
            return "Zinbiel[%s]" % ','.join(self.vars)
        return "Zinbiel[{}]".format(self.vars)
