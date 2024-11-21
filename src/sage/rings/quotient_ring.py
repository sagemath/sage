r"""
Quotient Rings

AUTHORS:

- William Stein
- Simon King (2011-04): Put it into the category framework, use the
  new coercion model.
- Simon King (2011-04): Quotients of non-commutative rings by
  twosided ideals.

TESTS::

    sage: R.<x> = PolynomialRing(ZZ)
    sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
    sage: S = R.quotient_ring(I)

.. TODO::

    The following skipped tests should be removed once :issue:`13999` is fixed::

        sage: TestSuite(S).run(skip=['_test_nonzero_equal', '_test_elements', '_test_zero'])

In :issue:`11068`, non-commutative quotient rings `R/I` were
implemented.  The only requirement is that the two-sided ideal `I`
provides a ``reduce`` method so that ``I.reduce(x)`` is the normal
form of an element `x` with respect to `I` (i.e., we have
``I.reduce(x) == I.reduce(y)`` if `x-y \in I`, and
``x - I.reduce(x) in I``). Here is a toy example::

    sage: from sage.rings.noncommutative_ideals import Ideal_nc
    sage: from itertools import product
    sage: class PowerIdeal(Ideal_nc):
    ....:     def __init__(self, R, n):
    ....:         self._power = n
    ....:         self._power = n
    ....:         Ideal_nc.__init__(self, R, [R.prod(m) for m in product(R.gens(), repeat=n)])
    ....:     def reduce(self, x):
    ....:         R = self.ring()
    ....:         return add([c*R(m) for m,c in x if len(m)<self._power],R(0))
    sage: F.<x,y,z> = FreeAlgebra(QQ, 3)                                                # needs sage.combinat sage.modules
    sage: I3 = PowerIdeal(F,3); I3                                                      # needs sage.combinat sage.modules
    Twosided Ideal (x^3, x^2*y, x^2*z, x*y*x, x*y^2, x*y*z, x*z*x, x*z*y,
    x*z^2, y*x^2, y*x*y, y*x*z, y^2*x, y^3, y^2*z, y*z*x, y*z*y, y*z^2,
    z*x^2, z*x*y, z*x*z, z*y*x, z*y^2, z*y*z, z^2*x, z^2*y, z^3) of
    Free Algebra on 3 generators (x, y, z) over Rational Field

Free algebras have a custom quotient method that serves at creating
finite dimensional quotients defined by multiplication matrices. We
are bypassing it, so that we obtain the default quotient::

    sage: # needs sage.combinat sage.modules
    sage: Q3.<a,b,c> = F.quotient(I3)
    sage: Q3
    Quotient of Free Algebra on 3 generators (x, y, z) over Rational Field by
    the ideal (x^3, x^2*y, x^2*z, x*y*x, x*y^2, x*y*z, x*z*x, x*z*y, x*z^2,
    y*x^2, y*x*y, y*x*z, y^2*x, y^3, y^2*z, y*z*x, y*z*y, y*z^2, z*x^2, z*x*y,
    z*x*z, z*y*x, z*y^2, z*y*z, z^2*x, z^2*y, z^3)
    sage: (a+b+2)^4
    16 + 32*a + 32*b + 24*a^2 + 24*a*b + 24*b*a + 24*b^2
    sage: Q3.is_commutative()
    False

Even though `Q_3` is not commutative, there is commutativity for
products of degree three::

    sage: a*(b*c)-(b*c)*a==F.zero()                                                     # needs sage.combinat sage.modules
    True

If we quotient out all terms of degree two then of course the resulting
quotient ring is commutative::

    sage: # needs sage.combinat sage.modules
    sage: I2 = PowerIdeal(F,2); I2
    Twosided Ideal (x^2, x*y, x*z, y*x, y^2, y*z, z*x, z*y, z^2) of Free Algebra
    on 3 generators (x, y, z) over Rational Field
    sage: Q2.<a,b,c> = F.quotient(I2)
    sage: Q2.is_commutative()
    True
    sage: (a+b+2)^4
    16 + 32*a + 32*b

Since :issue:`7797`, there is an implementation of free algebras
based on Singular's implementation of the Letterplace Algebra. Our
letterplace wrapper allows to provide the above toy example more
easily::

    sage: # needs sage.combinat sage.libs.singular sage.modules
    sage: from itertools import product
    sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
    sage: Q3 = F.quo(F*[F.prod(m) for m in product(F.gens(), repeat=3)]*F)
    sage: Q3
    Quotient of Free Associative Unital Algebra on 3 generators (x, y, z)
     over Rational Field by the ideal (x*x*x, x*x*y, x*x*z, x*y*x, x*y*y, x*y*z,
      x*z*x, x*z*y, x*z*z, y*x*x, y*x*y, y*x*z, y*y*x, y*y*y, y*y*z, y*z*x, y*z*y,
      y*z*z, z*x*x, z*x*y, z*x*z, z*y*x, z*y*y, z*y*z, z*z*x, z*z*y, z*z*z)
    sage: Q3.0*Q3.1 - Q3.1*Q3.0
    xbar*ybar - ybar*xbar
    sage: Q3.0*(Q3.1*Q3.2) - (Q3.1*Q3.2)*Q3.0
    0
    sage: Q2 = F.quo(F*[F.prod(m) for m in product(F.gens(), repeat=2)]*F)
    sage: Q2.is_commutative()
    True
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************
import sage.interfaces.abc
import sage.misc.latex as latex
import sage.structure.parent_gens
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.rings import Rings
from sage.misc.cachefunc import cached_method
from sage.rings import ideal, quotient_ring_element, ring
from sage.structure.category_object import normalize_names
from sage.structure.richcmp import richcmp, richcmp_method

_Rings = Rings()
_CommRings = CommutativeRings()


MPolynomialIdeal_quotient = None


def QuotientRing(R, I, names=None, **kwds):
    r"""
    Create a quotient ring of the ring `R` by the twosided ideal `I`.

    Variables are labeled by ``names`` (if the quotient ring is a quotient
    of a polynomial ring).  If ``names`` isn't given, 'bar' will be appended
    to the variable names in `R`.

    INPUT:

    - ``R`` -- a ring

    - ``I`` -- a twosided ideal of `R`

    - ``names`` -- (optional) a list of strings to be used as names for
      the variables in the quotient ring `R/I`

    - further named arguments that will be passed to the constructor
      of the quotient ring instance

    OUTPUT: `R/I` - the quotient ring `R` mod the ideal `I`

    ASSUMPTION:

    ``I`` has a method ``I.reduce(x)`` returning the normal form
    of elements `x\in R`. In other words, it is required that
    ``I.reduce(x)==I.reduce(y)`` `\iff x-y \in I`, and
    ``x-I.reduce(x) in I``, for all `x,y\in R`.

    EXAMPLES:

    Some simple quotient rings with the integers::

        sage: R = QuotientRing(ZZ, 7*ZZ); R
        Quotient of Integer Ring by the ideal (7)
        sage: R.gens()
        (1,)
        sage: 1*R(3); 6*R(3); 7*R(3)
        3
        4
        0

    ::

        sage: S = QuotientRing(ZZ,ZZ.ideal(8)); S
        Quotient of Integer Ring by the ideal (8)
        sage: 2*S(4)
        0

    With polynomial rings (note that the variable name of the quotient
    ring can be specified as shown below)::

        sage: # needs sage.libs.pari
        sage: P.<x> = QQ[]
        sage: R.<xx> = QuotientRing(P, P.ideal(x^2 + 1))
        sage: R
        Univariate Quotient Polynomial Ring in xx over Rational Field
         with modulus x^2 + 1
        sage: R.gens(); R.gen()
        (xx,)
        xx
        sage: for n in range(4): xx^n
        1
        xx
        -1
        -xx

    ::

        sage: # needs sage.libs.pari
        sage: P.<x> = QQ[]
        sage: S = QuotientRing(P, P.ideal(x^2 - 2))
        sage: S
        Univariate Quotient Polynomial Ring in xbar over Rational Field
         with modulus x^2 - 2
        sage: xbar = S.gen(); S.gen()
        xbar
        sage: for n in range(3): xbar^n
        1
        xbar
        2

    Sage coerces objects into ideals when possible::

        sage: P.<x> = QQ[]
        sage: R = QuotientRing(P, x^2 + 1); R                                           # needs sage.libs.pari
        Univariate Quotient Polynomial Ring in xbar over Rational Field
         with modulus x^2 + 1

    By Noether's homomorphism theorems, the quotient of a quotient ring
    of `R` is just the quotient of `R` by the sum of the ideals. In this
    example, we end up modding out the ideal `(x)` from the ring
    `\QQ[x,y]`::

        sage: # needs sage.libs.pari sage.libs.singular
        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = QuotientRing(R, R.ideal(1 + y^2))
        sage: T.<c,d> = QuotientRing(S, S.ideal(a))
        sage: T
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field
         by the ideal (x, y^2 + 1)
        sage: R.gens(); S.gens(); T.gens()
        (x, y)
        (a, b)
        (0, d)
        sage: for n in range(4): d^n
        1
        d
        -1
        -d

    TESTS:

    By :issue:`11068`, the following does not return a generic
    quotient ring but a usual quotient of the integer ring::

        sage: R = Integers(8)
        sage: I = R.ideal(2)
        sage: R.quotient(I)
        Ring of integers modulo 2

    Here is an example of the quotient of a free algebra by a
    twosided homogeneous ideal (see :issue:`7797`)::

        sage: # needs sage.combinat sage.libs.singular sage.modules
        sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
        sage: I = F * [x*y + y*z, x^2 + x*y - y*x - y^2] * F
        sage: Q.<a,b,c> = F.quo(I); Q
        Quotient of Free Associative Unital Algebra on 3 generators (x, y, z)
         over Rational Field by the ideal (x*y + y*z, x*x + x*y - y*x - y*y)
        sage: a*b
        -b*c
        sage: a^3
        -b*c*a - b*c*b - b*c*c
        sage: J = Q * [a^3 - b^3] * Q
        sage: R.<i,j,k> = Q.quo(J); R
        Quotient of Free Associative Unital Algebra on 3 generators (x, y, z)
         over Rational Field by the ideal
         (-y*y*z - y*z*x - 2*y*z*z, x*y + y*z, x*x + x*y - y*x - y*y)
        sage: i^3
        -j*k*i - j*k*j - j*k*k
        sage: j^3
        -j*k*i - j*k*j - j*k*k

    Check that :issue:`5978` is fixed by if we quotient by the zero ideal `(0)`
    then we just return ``R``::

        sage: R = QQ['x']
        sage: R.quotient(R.zero_ideal())
        Univariate Polynomial Ring in x over Rational Field
        sage: R.<x> = PolynomialRing(ZZ)
        sage: R is R.quotient(R.zero_ideal())
        True
        sage: I = R.ideal(0)
        sage: R is R.quotient(I)
        True
    """
    # 1. Not all rings inherit from the base class of rings.
    # 2. We want to support quotients of free algebras by homogeneous two-sided ideals.
    from sage.rings.finite_rings.integer_mod_ring import Integers
    from sage.rings.integer_ring import ZZ
    if R not in _Rings:
        raise TypeError("R must be a ring")
    is_commutative = R in _CommRings
    if names is None:
        try:
            names = tuple([x + 'bar' for x in R.variable_names()])
        except ValueError:  # no names are assigned
            pass
    else:
        names = normalize_names(R.ngens(), names)
    if kwds.get('implementation') == 'pbori':
        from sage.rings.polynomial.polynomial_ring_constructor import (
            BooleanPolynomialRing_constructor as BooleanPolynomialRing,
        )
        kwds.pop('implementation')
        return BooleanPolynomialRing(R.ngens(), names=names, **kwds)
    # workaround to silence warning from #34806
    from sage.rings.abc import Order
    if isinstance(R, Order):
        if not R.is_maximal():
            raise NotImplementedError('only implemented for maximal orders')
        I = R.number_field().ideal(I)
    elif not isinstance(I, ideal.Ideal_generic) or I.ring() != R:
        I = R.ideal(I)
    if I.is_zero():
        return R
    try:
        if I.is_principal():
            return R.quotient_by_principal_ideal(I.gen(), names)
    except (AttributeError, NotImplementedError):
        pass
    if not is_commutative:
        try:
            if I.side() != 'twosided':
                raise AttributeError
        except AttributeError:
            raise TypeError("A twosided ideal is required.")
    if isinstance(R, QuotientRing_nc):
        pi = R.cover()
        S = pi.domain()
        G = [pi.lift(x) for x in I.gens()]
        I_lift = S.ideal(G)
        J = R.defining_ideal()
        if S == ZZ:
            return Integers((I_lift+J).gen(), **kwds)
        return R.__class__(S, I_lift + J, names=names)
    if R in _CommRings:
        return QuotientRing_generic(R, I, names, **kwds)
    return QuotientRing_nc(R, I, names, **kwds)


def is_QuotientRing(x):
    """
    Test whether or not ``x`` inherits from :class:`QuotientRing_nc`.

    EXAMPLES::

        sage: from sage.rings.quotient_ring import is_QuotientRing
        sage: R.<x> = PolynomialRing(ZZ,'x')
        sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
        sage: S = R.quotient_ring(I)
        sage: is_QuotientRing(S)
        doctest:warning...
        DeprecationWarning: The function is_QuotientRing is deprecated;
        use 'isinstance(..., QuotientRing_nc)' instead.
        See https://github.com/sagemath/sage/issues/38266 for details.
        True
        sage: is_QuotientRing(R)
        False

    ::

        sage: # needs sage.combinat sage.libs.singular sage.modules
        sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
        sage: I = F * [x*y + y*z, x^2 + x*y - y*x - y^2] * F
        sage: Q = F.quo(I)
        sage: is_QuotientRing(Q)
        True
        sage: is_QuotientRing(F)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38266,
                "The function is_QuotientRing is deprecated; "
                "use 'isinstance(..., QuotientRing_nc)' instead.")
    return isinstance(x, QuotientRing_nc)


_RingsQuotients = _Rings.Quotients()
_CommutativeRingsQuotients = _CommRings.Quotients()
from sage.structure.category_object import check_default_category


@richcmp_method
class QuotientRing_nc(ring.Ring, sage.structure.parent_gens.ParentWithGens):
    """
    The quotient ring of `R` by a twosided ideal `I`.

    This class is for rings that do not inherit from
    :class:`~sage.rings.ring.CommutativeRing`.

    EXAMPLES:

    Here is a quotient of a free algebra by a twosided homogeneous ideal::

        sage: # needs sage.combinat sage.libs.singular sage.modules
        sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
        sage: I = F * [x*y + y*z, x^2 + x*y - y*x - y^2]*F
        sage: Q.<a,b,c> = F.quo(I); Q
        Quotient of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
         by the ideal (x*y + y*z, x*x + x*y - y*x - y*y)
        sage: a*b
        -b*c
        sage: a^3
        -b*c*a - b*c*b - b*c*c

    A quotient of a quotient is just the quotient of the original top
    ring by the sum of two ideals::

        sage: # needs sage.combinat sage.libs.singular sage.modules
        sage: J = Q * [a^3 - b^3] * Q
        sage: R.<i,j,k> = Q.quo(J); R
        Quotient of
         Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
         by the ideal (-y*y*z - y*z*x - 2*y*z*z, x*y + y*z, x*x + x*y - y*x - y*y)
        sage: i^3
        -j*k*i - j*k*j - j*k*k
        sage: j^3
        -j*k*i - j*k*j - j*k*k

    For rings that *do* inherit from :class:`~sage.rings.ring.CommutativeRing`,
    we provide a subclass :class:`QuotientRing_generic`, for backwards
    compatibility.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(ZZ,'x')
        sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
        sage: S = R.quotient_ring(I); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring
         by the ideal (x^2 + 3*x + 4, x^2 + 1)

    ::

        sage: R.<x,y> = PolynomialRing(QQ)
        sage: S.<a,b> = R.quo(x^2 + y^2)                                                # needs sage.libs.singular
        sage: a^2 + b^2 == 0                                                            # needs sage.libs.singular
        True
        sage: S(0) == a^2 + b^2                                                         # needs sage.libs.singular
        True

    Again, a quotient of a quotient is just the quotient of the original top
    ring by the sum of two ideals.

    ::

        sage: # needs sage.libs.singular
        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quo(1 + y^2)
        sage: T.<c,d> = S.quo(a)
        sage: T
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field
         by the ideal (x, y^2 + 1)
        sage: T.gens()
        (0, d)
    """
    Element = quotient_ring_element.QuotientRingElement

    def __init__(self, R, I, names, category=None):
        """
        Create the quotient ring of `R` by the twosided ideal `I`.

        INPUT:

        - ``R`` -- a ring

        - ``I`` -- a twosided ideal of `R`

        - ``names`` -- list of generator names

        EXAMPLES::

            sage: # needs sage.combinat sage.libs.singular sage.modules
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F * [x*y + y*z, x^2 + x*y - y*x - y^2] * F
            sage: Q.<a,b,c> = F.quo(I); Q
            Quotient of
             Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
             by the ideal (x*y + y*z, x*x + x*y - y*x - y*y)
            sage: a*b
            -b*c
            sage: a^3
            -b*c*a - b*c*b - b*c*c
        """
        if R not in _Rings:
            raise TypeError("The first argument must be a ring, but %s is not" % R)
        # workaround to silence warning from #34806
        from sage.rings.abc import Order
        if isinstance(R, Order):
            M = R.number_field().ideal_monoid()
        else:
            M = R.ideal_monoid()
        if I not in M:
            raise TypeError("The second argument must be an ideal of the given ring, but %s is not" % I)
        self.__R = R
        self.__I = I
        #sage.structure.parent_gens.ParentWithGens.__init__(self, R.base_ring(), names)
        ##
        # Unfortunately, computing the join of categories, which is done in
        # check_default_category, is very expensive.
        # However, we don't just want to use the given category without mixing in
        # some quotient stuff - unless Parent.__init__ was called
        # previously, in which case the quotient ring stuff is just
        # a vaste of time. This is the case for FiniteField_prime_modn.
        if not self._is_category_initialized():
            if category is None:
                try:
                    commutative = R.is_commutative()
                except (AttributeError, NotImplementedError):
                    commutative = False
                if commutative:
                    category = check_default_category(_CommutativeRingsQuotients,category)
                else:
                    category = check_default_category(_RingsQuotients,category)
            ring.Ring.__init__(self, R.base_ring(), names=names, category=category)
        # self._populate_coercion_lists_([R]) # we don't want to do this, since subclasses will often implement improved coercion maps.

    def construction(self):
        """
        Return the functorial construction of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ,'x')
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: R.quotient_ring(I).construction()
            (QuotientFunctor, Univariate Polynomial Ring in x over Integer Ring)

            sage: # needs sage.combinat sage.libs.singular sage.modules
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F * [x*y + y*z, x^2 + x*y - y*x - y^2] * F
            sage: Q = F.quo(I)
            sage: Q.construction()
            (QuotientFunctor,
             Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field)

        TESTS::

            sage: F, R = Integers(5).construction()
            sage: F(R)
            Ring of integers modulo 5
            sage: F, R = GF(5).construction()
            sage: F(R)
            Finite Field of size 5
        """
        from sage.categories.pushout import QuotientFunctor

        # Is there a better generic way to distinguish between things like Z/pZ as a field and Z/pZ as a ring?
        from sage.rings.ring import Field
        try:
            names = self.variable_names()
        except ValueError:
            try:
                names = self.cover_ring().variable_names()
            except ValueError:
                names = None
        if self in _CommRings:
            return QuotientFunctor(self.__I, names=names, domain=_CommRings,
                                   codomain=_CommRings,
                                   as_field=isinstance(self, Field)), self.__R
        else:
            return QuotientFunctor(self.__I, names=names,
                                   as_field=isinstance(self, Field)), self.__R

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ,'x')
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: R.quotient_ring(I)._repr_()
            'Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 3*x + 4, x^2 + 1)'
        """
        return "Quotient of %s by the ideal %s" % (self.cover_ring(), self.defining_ideal()._repr_short())

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ,'x')
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: R.quotient_ring(I)._latex_()
            '\\Bold{Z}[x]/\\left(x^{2} + 3x + 4, x^{2} + 1\\right)\\Bold{Z}[x]'
        """
        return "%s/%s" % (latex.latex(self.cover_ring()), latex.latex(self.defining_ideal()))

    def is_commutative(self):
        """
        Tell whether this quotient ring is commutative.

        .. NOTE::

            This is certainly the case if the cover ring is commutative.
            Otherwise, if this ring has a finite number of generators, it
            is tested whether they commute. If the number of generators is
            infinite, a :exc:`NotImplementedError` is raised.

        AUTHOR:

        - Simon King (2011-03-23): See :issue:`7797`.

        EXAMPLES:

        Any quotient of a commutative ring is commutative::

            sage: P.<a,b,c> = QQ[]
            sage: P.quo(P.random_element()).is_commutative()
            True

        The non-commutative case is more interesting::

            sage: # needs sage.combinat sage.libs.singular sage.modules
            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F * [x*y + y*z, x^2 + x*y - y*x - y^2] * F
            sage: Q = F.quo(I)
            sage: Q.is_commutative()
            False
            sage: Q.1*Q.2 == Q.2*Q.1
            False

        In the next example, the generators apparently commute::

            sage: # needs sage.combinat sage.libs.singular sage.modules
            sage: J = F * [x*y - y*x, x*z - z*x, y*z - z*y, x^3 - y^3] * F
            sage: R = F.quo(J)
            sage: R.is_commutative()
            True
        """
        try:
            if self.__R.is_commutative():
                return True
        except (AttributeError, NotImplementedError):
            pass
        from sage.rings.infinity import Infinity
        if self.ngens() == Infinity:
            raise NotImplementedError("This quotient ring has an infinite number of generators.")
        for i in range(self.ngens()):
            gi = self.gen(i)
            for j in range(i + 1, self.ngens()):
                gj = self.gen(j)
                if gi * gj != gj * gi:
                    return False
        return True

    @cached_method
    def cover(self):
        r"""
        The covering ring homomorphism `R \to R/I`, equipped with a section.

        EXAMPLES::

            sage: R = ZZ.quo(3 * ZZ)
            sage: pi = R.cover()
            sage: pi
            Ring morphism:
              From: Integer Ring
              To:   Ring of integers modulo 3
              Defn: Natural quotient map
            sage: pi(5)
            2
            sage: l = pi.lift()

        ::

            sage: # needs sage.libs.singular
            sage: R.<x,y>  = PolynomialRing(QQ)
            sage: Q = R.quo((x^2, y^2))
            sage: pi = Q.cover()
            sage: pi(x^3 + y)
            ybar
            sage: l = pi.lift(x + y^3)
            sage: l
            x
            sage: l = pi.lift(); l
            Set-theoretic ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field
                    by the ideal (x^2, y^2)
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn: Choice of lifting map
            sage: l(x + y^3)
            x
        """
        try:
            return self.__cover
        except AttributeError:
            from . import morphism
            pi = morphism.RingHomomorphism_cover(self.__R.Hom(self))
            lift = self.lifting_map()
            pi._set_lift(lift)
            self.__cover = pi
            return self.__cover

    @cached_method
    def lifting_map(self):
        """
        Return the lifting map to the cover.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S = R.quotient(x^2 + y^2)
            sage: pi = S.cover(); pi
            Ring morphism:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field
                    by the ideal (x^2 + y^2)
              Defn: Natural quotient map
            sage: L = S.lifting_map(); L
            Set-theoretic ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field
                    by the ideal (x^2 + y^2)
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn: Choice of lifting map
            sage: L(S.0)
            x
            sage: L(S.1)
            y

        Note that some reduction may be applied so that the lift of a
        reduction need not equal the original element::

            sage: z = pi(x^3 + 2*y^2); z                                                # needs sage.libs.singular
            -xbar*ybar^2 + 2*ybar^2
            sage: L(z)                                                                  # needs sage.libs.singular
            -x*y^2 + 2*y^2
            sage: L(z) == x^3 + 2*y^2                                                   # needs sage.libs.singular
            False

        Test that there also is a lift for rings that are no
        instances of :class:`~sage.rings.ring.Ring` (see :issue:`11068`)::

            sage: # needs sage.modules
            sage: MS = MatrixSpace(GF(5), 2, 2)
            sage: I = MS * [MS.0*MS.1, MS.2 + MS.3] * MS
            sage: Q = MS.quo(I)
            sage: Q.lift()
            Set-theoretic ring morphism:
              From: Quotient of Full MatrixSpace of 2 by 2 dense matrices
                    over Finite Field of size 5 by the ideal
            (
              [0 1]
              [0 0],
            <BLANKLINE>
              [0 0]
              [1 1]
            )
            <BLANKLINE>
              To:   Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 5
              Defn: Choice of lifting map
        """
        try:
            return self.__lift
        except AttributeError:
            pass
        from .morphism import RingMap_lift
        m = RingMap_lift(self, self.__R)
        self.__lift = m
        return m

    # The following is to make the category framework happy.
    def lift(self, x=None):
        """
        Return the lifting map to the cover, or the image
        of an element under the lifting map.

        .. NOTE::

            The category framework imposes that ``Q.lift(x)`` returns
            the image of an element `x` under the lifting map. For
            backwards compatibility, we let ``Q.lift()`` return the
            lifting map.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S = R.quotient(x^2 + y^2)
            sage: S.lift()                                                              # needs sage.libs.singular
            Set-theoretic ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field
                    by the ideal (x^2 + y^2)
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn: Choice of lifting map
            sage: S.lift(S.0) == x                                                      # needs sage.libs.singular
            True
        """
        if x is None:
            return self.lifting_map()
        return self.lifting_map()(x)

    def retract(self, x):
        """
        The image of an element of the cover ring under the quotient map.

        INPUT:

        - ``x`` -- an element of the cover ring

        OUTPUT: the image of the given element in ``self``

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S = R.quotient(x^2 + y^2)
            sage: S.retract((x+y)^2)                                                    # needs sage.libs.singular
            2*xbar*ybar
        """
        return self.cover()(x)

    def characteristic(self):
        r"""
        Return the characteristic of the quotient ring.

        .. TODO::

            Not yet implemented!

        EXAMPLES::

            sage: Q = QuotientRing(ZZ,7*ZZ)
            sage: Q.characteristic()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def defining_ideal(self):
        r"""
        Return the ideal generating this quotient ring.

        EXAMPLES:

        In the integers::

            sage: Q = QuotientRing(ZZ,7*ZZ)
            sage: Q.defining_ideal()
            Principal ideal (7) of Integer Ring

        An example involving a quotient of a quotient. By Noether's
        homomorphism theorems, this is actually a quotient by a sum of two
        ideals::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = QuotientRing(R, R.ideal(1 + y^2))
            sage: T.<c,d> = QuotientRing(S, S.ideal(a))
            sage: S.defining_ideal()
            Ideal (y^2 + 1) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: T.defining_ideal()
            Ideal (x, y^2 + 1) of Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self.__I

    @cached_method
    def is_field(self, proof=True):
        r"""
        Return ``True`` if the quotient ring is a field. Checks to see if the
        defining ideal is maximal.

        TESTS::

            sage: Q = QuotientRing(ZZ, 7*ZZ)
            sage: Q.is_field()
            True

        Requires the ``is_maximal`` method of the defining ideal to be
        implemented::

            sage: R.<x, y> = ZZ[]
            sage: R.quotient_ring(R.ideal([2, 4 + x])).is_field()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if proof:
            return self.defining_ideal().is_maximal()
        else:
            try:
                return self.defining_ideal().is_maximal()
            except NotImplementedError:
                return False

    @cached_method
    def is_integral_domain(self, proof=True):
        r"""
        With ``proof`` equal to ``True``  (the default), this function may
        raise a :exc:`NotImplementedError`.

        When ``proof`` is ``False``, if ``True`` is returned, then ``self`` is
        definitely an integral domain.  If the function returns ``False``,
        then either ``self`` is not an integral domain or it was unable to
        determine whether or not ``self`` is an integral domain.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: R.quo(x^2 - y).is_integral_domain()                                   # needs sage.libs.singular
            True
            sage: R.quo(x^2 - y^2).is_integral_domain()                                 # needs sage.libs.singular
            False
            sage: R.quo(x^2 - y^2).is_integral_domain(proof=False)                      # needs sage.libs.singular
            False
            sage: R.<a,b,c> = ZZ[]
            sage: Q = R.quotient_ring([a, b])
            sage: Q.is_integral_domain()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: Q.is_integral_domain(proof=False)
            False
        """
        if proof:
            return self.defining_ideal().is_prime()
        else:
            try:
                return self.defining_ideal().is_prime()
            except NotImplementedError:
                return False

    def is_noetherian(self):
        r"""
        Return ``True`` if this ring is Noetherian.

        EXAMPLES::

            sage: R = QuotientRing(ZZ, 102 * ZZ)
            sage: R.is_noetherian()
            True

            sage: P.<x> = QQ[]
            sage: R = QuotientRing(P, x^2 + 1)                                          # needs sage.libs.pari
            sage: R.is_noetherian()
            True

        If the cover ring of ``self`` is not Noetherian, we currently
        have no way of testing whether ``self`` is Noetherian, so we
        raise an error::

            sage: R.<x> = InfinitePolynomialRing(QQ)
            sage: R.is_noetherian()
            False
            sage: I = R.ideal([x[1]^2, x[2]])
            sage: S = R.quotient(I)
            sage: S.is_noetherian()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # Naive test: if this is the quotient of a Noetherian ring,
        # then it is Noetherian.  Otherwise we give up.
        if self.cover_ring().is_noetherian():
            return True

        raise NotImplementedError

    def cover_ring(self):
        r"""
        Return the cover ring of the quotient ring: that is, the original
        ring `R` from which we modded out an ideal, `I`.

        EXAMPLES::

            sage: Q = QuotientRing(ZZ, 7 * ZZ)
            sage: Q.cover_ring()
            Integer Ring

        ::

            sage: P.<x> = QQ[]
            sage: Q = QuotientRing(P, x^2 + 1)                                          # needs sage.libs.pari
            sage: Q.cover_ring()                                                        # needs sage.libs.pari
            Univariate Polynomial Ring in x over Rational Field
        """
        return self.__R

    # This is to make the category framework happy
    ambient = cover_ring

    def ideal(self, *gens, **kwds):
        """
        Return the ideal of ``self`` with the given generators.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2 + y^2)
            sage: S.ideal()                                                             # needs sage.libs.singular
            Ideal (0) of Quotient of Multivariate Polynomial Ring in x, y
             over Rational Field by the ideal (x^2 + y^2)
            sage: S.ideal(x + y + 1)                                                    # needs sage.libs.singular
            Ideal (xbar + ybar + 1) of Quotient of Multivariate Polynomial Ring in x, y
             over Rational Field by the ideal (x^2 + y^2)

        TESTS:

        We create an ideal of a fairly generic integer ring (see
        :issue:`5666`)::

            sage: R = Integers(10)
            sage: R.ideal(1)
            Principal ideal (1) of Ring of integers modulo 10
        """
        if len(gens) == 1:
            gens = gens[0]
        from sage.rings.polynomial.multi_polynomial_ring_base import (
            MPolynomialRing_base,
        )
        if not (isinstance(self.__R, MPolynomialRing_base) and self.__R._has_singular):
            # pass through
            return super().ideal(gens, **kwds)
        if isinstance(gens, sage.interfaces.abc.SingularElement):
            gens = list(gens)
        elif not isinstance(gens, (list, tuple)):
            gens = [gens]
        if 'coerce' in kwds and kwds['coerce']:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!

        global MPolynomialIdeal_quotient
        if MPolynomialIdeal_quotient is None:
            from sage.rings.polynomial.multi_polynomial_ideal import (
                MPolynomialIdeal_quotient,
            )
        return MPolynomialIdeal_quotient(self, gens, **kwds)

    def _element_constructor_(self, x, coerce=True):
        """
        Construct an element with ``self`` as the parent.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2 + y^2)
            sage: S(x)  # indirect doctest                                              # needs sage.libs.singular
            xbar
            sage: S(x^2 + y^2)                                                          # needs sage.libs.singular
            0

        The rings that coerce into the quotient ring canonically, are:

        - this ring

        - anything that coerces into the ring of which this is the
          quotient

        ::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quotient(x^2 + y^2)
            sage: S.coerce(0)
            0
            sage: S.coerce(2/3)
            2/3
            sage: S.coerce(a^2 - b)
            -b^2 - b
            sage: S.coerce(GF(7)(3))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field of size 7
            to Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)

        TESTS::

            sage: S(x, coerce=False)                                                    # needs sage.libs.singular
            a
        """
        if isinstance(x, quotient_ring_element.QuotientRingElement):
            if x.parent() is self:
                return x
            x = x.lift()
        if isinstance(x, sage.interfaces.abc.SingularElement):
            # self._singular_().set_ring()
            x = self.element_class(self, x.sage_poly(self.cover_ring()))
            return x
        if coerce:
            R = self.cover_ring()
            x = R(x)
        return self.element_class(self, x)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion map from ``R`` to ``self``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2 + y^2)
            sage: S.has_coerce_map_from(R)  # indirect doctest
            True
            sage: S.has_coerce_map_from(QQ)
            True
            sage: T = S.quotient_ring(x^3 - y)                                          # needs sage.libs.singular
            sage: S.has_coerce_map_from(T)                                              # needs sage.libs.singular
            False
            sage: T.has_coerce_map_from(R)                                              # needs sage.libs.singular
            True

        TESTS:

        We check that :issue:`13682` is fixed::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: I = R.ideal(x^2 + y^2)
            sage: J = R.ideal(x^2 + y^2, x^3 - y)
            sage: S = R.quotient(I)
            sage: T = R.quotient(J)

            sage: # needs sage.libs.singular
            sage: I < J
            True
            sage: T.has_coerce_map_from(S)
            True
            sage: S.quotient_ring(x^4 - x*y + 1).has_coerce_map_from(S)
            True
            sage: S.has_coerce_map_from(T)
            False

        We also allow coercions with the cover rings::

            sage: Rp.<x,y> = PolynomialRing(ZZ)
            sage: Ip = Rp.ideal(x^2 + y^2)
            sage: Jp = Rp.ideal(x^2 + y^2, x^3 - y)
            sage: Sp = Rp.quotient(Ip)
            sage: Tp = Rp.quotient(Jp)
            sage: R.has_coerce_map_from(Rp)
            True
            sage: Sp.has_coerce_map_from(Sp)
            True
            sage: T.has_coerce_map_from(Sp)                                             # needs sage.libs.singular
            True
            sage: Sp.has_coerce_map_from(T)                                             # needs sage.libs.singular
            False
        """
        C = self.cover_ring()
        if isinstance(R, QuotientRing_nc):
            if C == R.cover_ring():
                if R.defining_ideal() <= self.defining_ideal():
                    return True
            elif C.has_coerce_map_from(R.cover_ring()):
                try:
                    if R.defining_ideal().change_ring(C) <= self.defining_ideal():
                        return True
                except AttributeError: # Not all ideals have a change_ring
                    pass
        return C.has_coerce_map_from(R)

    def __richcmp__(self, other, op):
        r"""
        Only quotients by the *same* ring and same ideal (with the same
        generators!!) are considered equal.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2 + y^2)
            sage: S == R.quotient_ring(x^2 + y^2)
            True

        The ideals `(x^2 + y^2)` and `(-x^2-y^2)` are
        equal, but since the generators are different, the corresponding
        quotient rings are not equal::

            sage: R.ideal(x^2 + y^2) == R.ideal(-x^2 - y^2)                             # needs sage.libs.singular
            True
            sage: R.quotient_ring(x^2 + y^2) == R.quotient_ring(-x^2 - y^2)
            False
        """
        if not isinstance(other, QuotientRing_nc):
            return NotImplemented
        return richcmp((self.cover_ring(), self.defining_ideal().gens()),
                       (other.cover_ring(), other.defining_ideal().gens()), op)

    def ngens(self):
        r"""
        Return the number of generators for this quotient ring.

        .. TODO::

            Note that ``ngens`` counts 0 as a generator. Does
            this make sense? That is, since 0 only generates itself and the
            fact that this is true for all rings, is there a way to "knock it
            off" of the generators list if a generator of some original ring is
            modded out?

        EXAMPLES::

            sage: R = QuotientRing(ZZ, 7*ZZ)
            sage: R.gens(); R.ngens()
            (1,)
            1

        ::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = QuotientRing(R, R.ideal(1 + y^2))
            sage: T.<c,d> = QuotientRing(S, S.ideal(a))
            sage: T
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field
             by the ideal (x, y^2 + 1)
            sage: R.gens(); S.gens(); T.gens()
            (x, y)
            (a, b)
            (0, d)
            sage: R.ngens(); S.ngens(); T.ngens()
            2
            2
            2
        """
        return self.cover_ring().ngens()

    def gen(self, i=0):
        r"""
        Return the `i`-th generator for this quotient ring.

        EXAMPLES::

            sage: R = QuotientRing(ZZ, 7*ZZ)
            sage: R.gen(0)
            1

        ::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = QuotientRing(R, R.ideal(1 + y^2))
            sage: T.<c,d> = QuotientRing(S, S.ideal(a))
            sage: T
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field
             by the ideal (x, y^2 + 1)
            sage: R.gen(0); R.gen(1)
            x
            y
            sage: S.gen(0); S.gen(1)
            a
            b
            sage: T.gen(0); T.gen(1)
            0
            d
        """
        return self(self.__R.gen(i))

    def _singular_(self, singular=None):
        """
        Return the Singular quotient ring of ``self`` if the base ring is
        coercible to Singular.

        If a valid Singular representation is found it is used otherwise a
        new 'qring' is created.

        INPUT:

        - ``singular`` -- Singular instance (default: the
          default Singular instance)

        .. NOTE::

           This method also sets the current ring in Singular to ``self``

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2 + y^2)
            sage: S._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            // quotient ring from ideal
            _[1]=x2+y2
        """
        if singular is None:
            from sage.interfaces.singular import singular

        try:
            Q = self.__singular
            if Q.parent() is not singular:
                raise ValueError
            Q._check_valid()
            return Q
        except (AttributeError, ValueError):
            self.__singular = self._singular_init_(singular)
            return self.__singular

    def _singular_init_(self, singular=None):
        """
        Return a newly created Singular quotient ring matching ``self`` if
        the base ring is coercible to Singular.

        See ``_singular_``

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2 + y^2)
            sage: T = S._singular_init_()                                               # needs sage.libs.singular
            sage: parent(S)
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category'>
            sage: parent(T)                                                             # needs sage.libs.singular
            Singular
        """
        if singular is None:
            from sage.interfaces.singular import singular
        self.__R._singular_().set_ring()
        self.__singular = singular("%s" % self.__I._singular_().name(),"qring")
        return self.__singular

    def _magma_init_(self, magma):
        r"""
        Return string that evaluates to Magma version of this quotient
        ring. This is called implicitly when doing conversions to Magma.

        INPUT:

        - ``magma`` -- a Magma instance

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(GF(2))
            sage: Q = P.quotient(sage.rings.ideal.FieldIdeal(P))
            sage: magma(Q)  # indirect doctest                      # optional - magma
            Affine Algebra of rank 2 over GF(2)
            Graded Reverse Lexicographical Order
            Variables: x, y
            Quotient relations:
            [
            x^2 + x,
            y^2 + y
            ]
        """
        R = magma(self.__R)
        I = magma(self.__I.gens())
        return "quo<%s|%s>" % (R.name(), I._ref())

    def term_order(self):
        """
        Return the term order of this ring.

        EXAMPLES::

            sage: P.<a,b,c> = PolynomialRing(QQ)
            sage: I = Ideal([a^2 - a, b^2 - b, c^2 - c])
            sage: Q = P.quotient(I)
            sage: Q.term_order()
            Degree reverse lexicographic term order
        """
        return self.__R.term_order()

    def random_element(self):
        r"""
        Return a random element of this quotient ring obtained by
        sampling a random element of the cover ring and reducing
        it modulo the defining ideal.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: S = R.quotient([x^3, y^2])
            sage: S.random_element()  # random
            -8/5*xbar^2 + 3/2*xbar*ybar + 2*xbar - 4/23

        TESTS:

        Make sure we are not just getting images of integers in this
        ring (which would be the case if the default implementation
        of this method was inherited from generic rings)::

            sage: any(S.random_element() not in ZZ for _ in range(999))
            True
        """
        return self.retract(self.cover_ring().random_element())


class QuotientRing_generic(QuotientRing_nc, ring.CommutativeRing):
    r"""
    Create a quotient ring of a *commutative* ring `R` by the ideal `I`.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(ZZ)
        sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
        sage: S = R.quotient_ring(I); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring
         by the ideal (x^2 + 3*x + 4, x^2 + 1)
    """

    def __init__(self, R, I, names, category=None):
        """
        Initialize ``self``.

        INPUT:

        - ``R`` -- a ring that is a :class:`~sage.rings.ring.CommutativeRing`

        - ``I`` -- an ideal of `R`

        - ``names`` -- list of generator names

        TESTS::

            sage: ZZ.quo(2) in Rings().Commutative()  # indirect doctest
            True
        """
        if R not in _CommRings:
            raise TypeError("This class is for quotients of commutative rings only.\n    For non-commutative rings, use <sage.rings.quotient_ring.QuotientRing_nc>")
        if not self._is_category_initialized():
            category = check_default_category(_CommutativeRingsQuotients,
                                              category)
        QuotientRing_nc.__init__(self, R, I, names, category=category)

    def _macaulay2_init_(self, macaulay2=None):
        r"""
        EXAMPLES:

        Quotients of multivariate polynomial rings over `\QQ`, `\ZZ` and
        a finite field::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: I = R.ideal([x^2 - y])
            sage: Q = R.quotient_ring(I); Q
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y)
            sage: Q._macaulay2_init_()                      # optional - macaulay2
            QQ[x...y]
            --------
              2
             x  - y

            sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
            sage: I = R.ideal([x*y - z^2, y^2 - w^2])
            sage: Q = R.quotient(I); Q
            Quotient of Multivariate Polynomial Ring in x, y, z, w over Integer Ring by the ideal (x*y - z^2, y^2 - w^2)
            sage: Q._macaulay2_init_()                      # optional - macaulay2
               ZZ[x...z, w]
            -------------------
                    2   2    2
            (x*y - z , y  - w )

            sage: R.<x,y> = PolynomialRing(GF(101), 2)
            sage: I = R.ideal([x^2 + x, y^2 + y])
            sage: Q = R.quotient_ring(I); Q
            Quotient of Multivariate Polynomial Ring in x, y over
             Finite Field of size 101 by the ideal (x^2 + x, y^2 + y)
            sage: Q._macaulay2_init_()                      # optional - macaulay2
                 ZZ
                ---[x...y]
                101
            ----------------
              2       2
            (x  + x, y  + y)

        Quotients of univariate polynomial rings::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: Q = R.quotient_ring(I); Q
            Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 3*x + 4, x^2 + 1)
            sage: Q._macaulay2_init_()                      # optional - macaulay2
                    ZZ[x]
            ---------------------
              2            2
            (x  + 3x + 4, x  + 1)
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default
        I = self.defining_ideal()._macaulay2_(macaulay2)
        return I.ring()._operator('/', I)

    def _ideal_class_(self, n=0):
        r"""
        Use a specialized class for quotient ring ideals.

        INPUT:

        - ``n`` -- integer (default: ``0``); the number of generators

        EXAMPLES::

            sage: type(Zmod(14).ideal(2))
            <class 'sage.rings.quotient_ring.QuotientRingIdeal_principal'>
            sage: type(Zmod(14).ideal([7]))
            <class 'sage.rings.quotient_ring.QuotientRingIdeal_principal'>
            sage: type(Zmod(14).ideal([2,7]))
            <class 'sage.rings.quotient_ring.QuotientRingIdeal_generic'>
        """
        if n == 1:
            return QuotientRingIdeal_principal
        return QuotientRingIdeal_generic


class QuotientRingIdeal_generic(ideal.Ideal_generic):
    r"""
    Specialized class for quotient-ring ideals.

    EXAMPLES::

        sage: Zmod(9).ideal([-6,9])
        Ideal (3, 0) of Ring of integers modulo 9
    """

    def _contains_(self, other):
        r"""
        Check whether this ideal contains the given element.

        EXAMPLES::

            sage: 1 in Zmod(15).ideal(2)
            True
            sage: 1 in Zmod(15).ideal(3)
            False
            sage: 1 in Zmod(15).ideal([5,10])
            False
            sage: 6 in Zmod(15).ideal([9])
            True
            sage: 1 in Zmod(15).ideal([3,5])
            True

        ::

            sage: # needs sage.libs.pari
            sage: R.<T> = QQ[]
            sage: S.<t> = R.quotient(T^3 - 1)
            sage: 1 in S.ideal(t^2 - 1)
            False
            sage: 7 in S.ideal(t^2 + 1)
            True
            sage: 5-5*t in S.ideal(t^2 - 1)
            True
        """
        R = self.ring()
        assert other in R
        if hasattr(R, 'defining_ideal'):
            Igens = list(R.defining_ideal().gens())
        else:
            Igens = [R.modulus()]
        Igens += [g.lift() for g in self.gens()]
        J = R.cover_ring().ideal(Igens)
        return other.lift() in J


class QuotientRingIdeal_principal(ideal.Ideal_principal, QuotientRingIdeal_generic):
    r"""
    Specialized class for principal quotient-ring ideals.

    EXAMPLES::

        sage: Zmod(9).ideal(-33)
        Principal ideal (3) of Ring of integers modulo 9
    """
