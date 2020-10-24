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

    TESTS::

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
    def __classcall_private__(cls, R, n=None, names=None, prefix=None):
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
        superclass = super(FreeZinbielAlgebra, cls)
        if names is None:
            return superclass.__classcall__(cls, R, n, None, prefix)
        return superclass.__classcall__(cls, R, n, tuple(names), prefix)

    def __init__(self, R, n, names, prefix):
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
        cat = MagmaticAlgebras(R).WithBasis().Graded()
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

        EXAMPLES::

            sage: A = algebras.FreeZinbiel(ZZ, 'f,g,h')
            sage: A.change_ring(QQ)
            Free Zinbiel algebra on generators (Z[f], Z[g], Z[h])
            over Rational Field
        """
        A = self.variable_names()
        return FreeZinbielAlgebra(R, n=len(A), names=A)

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

    def product_on_basis(self, x, y):
        """
        Return the product of the basis elements indexed by ``x`` and ``y``.

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
            return self.element_class(self,
                                      x.monomial_coefficients(copy=False))
        else:
            raise TypeError('not able to convert this to this algebra')
        # Ok, not a Zinbiel algebra element (or should not be viewed as one).

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - free Zinbiel algebras whose set `E` of labels is
          a subset of the corresponding self of ``set`, and whose base
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
            return (all(x in self.variable_names() for x in R.variable_names())
                    and self.base_ring().has_coerce_map_from(R.base_ring()))
        return super(FreeZinbielAlgebra, self)._coerce_map_from_(R)

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
        return ZinbielFunctor(A), self.base_ring()


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

    def __init__(self, variables):
        """
        EXAMPLES::

            sage: functor = sage.algebras.free_zinbiel_algebra.ZinbielFunctor
            sage: F = functor(['x','y']); F
            Zinbiel[x,y]
            sage: F(ZZ)
            Free Zinbiel algebra on generators (Z[x], Z[y]) over Integer Ring
        """
        Functor.__init__(self, Rings(), Magmas())
        self.vars = variables
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
            return FreeZinbielAlgebra(R, len(self.vars), self.vars)
        return FreeZinbielAlgebra(R, self.vars)

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
        return self.vars == other.vars

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
            sage: F = functor(['x','y'])
            sage: G = functor(['t'])
            sage: G * F
            Zinbiel[x,y,t]

        With an infinite generating set::

            sage: H = functor(ZZ)
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
            if not self._finite_vars or not other._finite_vars:
                raise CoercionException("Unable to determine overlap for infinite sets")
            if set(self.vars).intersection(other.vars):
                raise CoercionException("Overlapping variables (%s,%s)" %
                                        (self.vars, other.vars))
            return ZinbielFunctor(other.vars + self.vars)
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
            sage: F = functor(['x','y'])
            sage: G = functor(['t'])
            sage: F.merge(G)
            Zinbiel[x,y,t]
            sage: F.merge(F)
            Zinbiel[x,y]

        With an infinite generating set::

            sage: H = functor(ZZ)
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
        """
        if isinstance(other, ZinbielFunctor):
            if self.vars == other.vars:
                return self
            def check(x):
                return isinstance(x, (list, tuple)) or x in Sets().Finite()
            if not check(self.vars) or not check(other.vars):
                return None
            ret = list(self.vars)
            cur_vars = set(ret)
            for v in other.vars:
                if v not in cur_vars:
                    ret.append(v)
            return ZinbielFunctor(ret)
        else:
            return None

    def _repr_(self):
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

