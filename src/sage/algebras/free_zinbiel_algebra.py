"""
Free Zinbiel Algebras

AUTHORS:

- Travis Scrimshaw (2015-09): initial version
"""

#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.magmatic_algebras import MagmaticAlgebras
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
from sage.combinat.words.alphabet import Alphabet
from sage.sets.family import Family

class FreeZinbielAlgebra(CombinatorialFreeModule):
    r"""
    The free Zinbiel algebra on `n` generators.

    Let `R` be a ring. A *Zinbiel algebra* is a non-associative
    algebra with multiplication `\circ` that satisfies

    .. MATH::

        a \circ (b \circ c) = a \circ (b \circ c) + a \circ (c \circ b).

    Zinbiel algebras were first introduced by Loday as the Koszul
    dual to Leibniz algebras (hence the name coined by Lemaire).

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

        Currently the basis is indexed by all words over the variables,
        incuding the empty word. This is a slight abuse as it is suppose
        to be the indexed by all non-empty words.

    EXAMPLES:

    We create the free Zinbiel algebra and check the defining relation::

        sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
        sage: (x*y)*z
        Z[xyz] + Z[xzy]
        sage: x*(y*z) + x*(z*y)
        Z[xyz] + Z[xzy]

    We see that the Zinbiel algebra is not associative, nor even
    power associative::

        sage: x*(y*z)
        Z[xyz]
        sage: x*(x*x)
        Z[xxx]
        sage: (x*x)*x
        2*Z[xxx]

    We verify that it is a divided powers algebra::

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

    REFERENCES:

    - :wikipedia:`Zinbiel_algebra`

    .. [Loday95] Jean-Louis Loday.
       *Cup-product for Leibniz cohomology and dual Leibniz algebras*.
       Math. Scand., pp. 189--196 (1995).
       http://www.math.uiuc.edu/K-theory/0015/cup_product.pdf
    .. [LV12] Jean-Louis Loday and Bruno Vallette. *Algebraic Operads*.
       Springer-Verlag Berlin Heidelberg (2012).
       :doi:`10.1007/978-3-642-30362-3`.
    """
    @staticmethod
    def __classcall_private__(cls, R, n=None, names=None):
        """
        Standardize input to ensure a unique representation.

        TESTS::

            sage: Z1.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: Z2.<x,y,z> = algebras.FreeZinbiel(QQ, 3)
            sage: Z3 = algebras.FreeZinbiel(QQ, 3, 'x,y,z')
            sage: Z4.<x,y,z> = algebras.FreeZinbiel(QQ, 'x,y,z')
            sage: Z1 is Z2 and Z1 is Z3 and Z1 is Z4
            True
        """
        if isinstance(n, (list,tuple)):
            names = n
            n = len(names)
        elif isinstance(n, str):
            names = n.split(',')
            n = len(names)
        elif isinstance(names, str):
            names = names.split(',')
        elif n is None:
            n = len(names)
        return super(FreeZinbielAlgebra, cls).__classcall__(cls, R, n, tuple(names))

    def __init__(self, R, n, names):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: TestSuite(Z).run()
        """
        if R not in Rings:
            raise TypeError("argument R must be a ring")
        indices = Words(Alphabet(n, names=names))
        cat = MagmaticAlgebras(R).WithBasis()
        self._n = n
        CombinatorialFreeModule.__init__(self, R, indices, prefix='Z',
                                         category=cat)
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

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: Z
            Free Zinbiel algebra on generators (Z[x], Z[y], Z[z]) over Rational Field
        """
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
        A = self.variable_names()
        return Family( A, lambda g: self.monomial(self._indices(g)) )

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: Z.gens()
            (Z[x], Z[y], Z[z])
        """
        return tuple(self.algebra_generators())

    def product_on_basis(self, x, y):
        """
        Return the product of the basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: Z.<x,y,z> = algebras.FreeZinbiel(QQ)
            sage: (x*y)*z  # indirect doctest
            Z[xyz] + Z[xzy]
        """
        if not x:
            return self.monomial(y)
        x0 = self._indices(x[0])
        return self.sum_of_monomials(x0 + sh for sh in x[1:].shuffle(y))

