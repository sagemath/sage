# sage_setup: distribution = sagemath-singular
"""
Homogeneous ideals of free algebras

For twosided ideals and when the base ring is a field, this implementation
also provides Groebner bases and ideal containment tests.

EXAMPLES::

    sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
    sage: F
    Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
    sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
    sage: I
    Twosided Ideal (x*y + y*z, x*x + x*y - y*x - y*y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field

One can compute Groebner bases out to a finite degree, can compute normal
forms and can test containment in the ideal::

    sage: I.groebner_basis(degbound=3)
    Twosided Ideal (x*y + y*z,
        x*x - y*x - y*y - y*z,
        y*y*y - y*y*z + y*z*y - y*z*z,
        y*y*x + y*y*z + y*z*x + y*z*z) of Free Associative Unital Algebra
        on 3 generators (x, y, z) over Rational Field
    sage: (x*y*z*y*x).normal_form(I)
    y*z*z*y*z + y*z*z*z*x + y*z*z*z*z
    sage: x*y*z*y*x - (x*y*z*y*x).normal_form(I) in I
    True

AUTHOR:

- Simon King (2011-03-22):  See :issue:`7797`.
"""
# ****************************************************************************
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.noncommutative_ideals import Ideal_nc
from sage.libs.singular.function import lib, singular_function
from sage.algebras.letterplace.free_algebra_letterplace cimport FreeAlgebra_letterplace, FreeAlgebra_letterplace_libsingular
from sage.algebras.letterplace.free_algebra_element_letterplace cimport FreeAlgebraElement_letterplace
from sage.rings.infinity import Infinity

#####################
# Define some singular functions
lib("freegb.lib")
singular_twostd = singular_function("twostd")
poly_reduce = singular_function("NF")


class LetterplaceIdeal(Ideal_nc):
    """
    Graded homogeneous ideals in free algebras.

    In the two-sided case over a field, one can compute Groebner bases
    up to a degree bound, normal forms of graded homogeneous elements
    of the free algebra, and ideal containment.

    EXAMPLES::

        sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
        sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
        sage: I
        Twosided Ideal (x*y + y*z, x*x + x*y - y*x - y*y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
        sage: I.groebner_basis(2)
        Twosided Ideal (x*y + y*z, x*x - y*x - y*y - y*z) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
        sage: I.groebner_basis(4)
        Twosided Ideal (x*y + y*z,
            x*x - y*x - y*y - y*z,
            y*y*y - y*y*z + y*z*y - y*z*z,
            y*y*x + y*y*z + y*z*x + y*z*z,
            y*y*z*y - y*y*z*z + y*z*z*y - y*z*z*z,
            y*z*y*y - y*z*y*z + y*z*z*y - y*z*z*z,
            y*y*z*x + y*y*z*z + y*z*z*x + y*z*z*z,
            y*z*y*x + y*z*y*z + y*z*z*x + y*z*z*z) of Free Associative Unital
            Algebra on 3 generators (x, y, z) over Rational Field

    Groebner bases are cached. If one has computed a Groebner basis
    out to a high degree then it will also be returned if a Groebner
    basis with a lower degree bound is requested::

        sage: I.groebner_basis(2) is I.groebner_basis(4)
        True

    Of course, the normal form of any element has to satisfy the following::

        sage: x*y*z*y*x - (x*y*z*y*x).normal_form(I) in I
        True

    Left and right ideals can be constructed, but only twosided ideals provide
    Groebner bases::

        sage: JL = F*[x*y+y*z,x^2+x*y-y*x-y^2]; JL
        Left Ideal (x*y + y*z, x*x + x*y - y*x - y*y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
        sage: JR = [x*y+y*z,x^2+x*y-y*x-y^2]*F; JR
        Right Ideal (x*y + y*z, x*x + x*y - y*x - y*y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
        sage: JR.groebner_basis(2)
        Traceback (most recent call last):
        ...
        TypeError: This ideal is not two-sided. We can only compute two-sided Groebner bases
        sage: JL.groebner_basis(2)
        Traceback (most recent call last):
        ...
        TypeError: This ideal is not two-sided. We can only compute two-sided Groebner bases

    Also, it is currently not possible to compute a Groebner basis when the base
    ring is not a field::

        sage: FZ.<a,b,c> = FreeAlgebra(ZZ, implementation='letterplace')
        sage: J = FZ*[a^3-b^3]*FZ
        sage: J.groebner_basis(2)
        Traceback (most recent call last):
        ...
        TypeError: Currently, we can only compute Groebner bases if the ring of coefficients is a field

    The letterplace implementation of free algebras also provides integral degree weights
    for the generators, and we can compute Groebner bases for twosided graded homogeneous
    ideals::

        sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace',degrees=[1,2,3])
        sage: I = F*[x*y+z-y*x,x*y*z-x^6+y^3]*F
        sage: I.groebner_basis(Infinity)
        Twosided Ideal (x*y - y*x + z,
        x*x*x*x*x*x - y*x*z - y*y*y + z*z,
        x*z*z - y*x*x*z + y*x*z*x + y*y*z + y*z*y + z*x*z + z*y*y - z*z*x,
        x*x*x*x*x*z + x*x*x*x*z*x + x*x*x*z*x*x + x*x*z*x*x*x + x*z*x*x*x*x +
        y*x*z*y - y*y*x*z + y*z*z + z*x*x*x*x*x - z*z*y,
        x*x*x*x*z*y*y + x*x*x*z*y*y*x - x*x*x*z*y*z - x*x*z*y*x*z + x*x*z*y*y*x*x +
        x*x*z*y*y*y - x*x*z*y*z*x - x*z*y*x*x*z - x*z*y*x*z*x +
        x*z*y*y*x*x*x + 2*x*z*y*y*y*x - 2*x*z*y*y*z - x*z*y*z*x*x -
        x*z*y*z*y + y*x*z*x*x*x*x*x - 4*y*x*z*x*x*z - 4*y*x*z*x*z*x +
        4*y*x*z*y*x*x*x + 3*y*x*z*y*y*x - 4*y*x*z*y*z + y*y*x*x*x*x*z +
        y*y*x*x*x*z*x - 3*y*y*x*x*z*x*x - y*y*x*x*z*y +
        5*y*y*x*z*x*x*x + 4*y*y*x*z*y*x - 4*y*y*y*x*x*z +
        4*y*y*y*x*z*x + 3*y*y*y*y*z + 4*y*y*y*z*x*x + 6*y*y*y*z*y +
        y*y*z*x*x*x*x + y*y*z*x*z + 7*y*y*z*y*x*x + 7*y*y*z*y*y -
        7*y*y*z*z*x - y*z*x*x*x*z - y*z*x*x*z*x + 3*y*z*x*z*x*x +
        y*z*x*z*y + y*z*y*x*x*x*x - 3*y*z*y*x*z + 7*y*z*y*y*x*x +
        3*y*z*y*y*y - 3*y*z*y*z*x - 5*y*z*z*x*x*x - 4*y*z*z*y*x +
        4*y*z*z*z - z*y*x*x*x*z - z*y*x*x*z*x - z*y*x*z*x*x -
        z*y*x*z*y + z*y*y*x*x*x*x - 3*z*y*y*x*z + 3*z*y*y*y*x*x +
        z*y*y*y*y - 3*z*y*y*z*x - z*y*z*x*x*x - 2*z*y*z*y*x +
        2*z*y*z*z - z*z*x*x*x*x*x + 4*z*z*x*x*z + 4*z*z*x*z*x -
        4*z*z*y*x*x*x - 3*z*z*y*y*x + 4*z*z*y*z + 4*z*z*z*x*x +
        2*z*z*z*y)
        of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field

    Again, we can compute normal forms::

        sage: (z*I.0-I.1).normal_form(I)
        0
        sage: (z*I.0-x*y*z).normal_form(I)
        -y*x*z + z*z
    """
    def __init__(self, ring, gens, coerce=True, side='twosided'):
        """
        INPUT:

        - ``ring`` -- a free algebra in letterplace implementation
        - ``gens`` -- list, tuple or sequence of generators
        - ``coerce`` -- boolean (default: ``True``); whether ``gens`` shall be
          coerced first
        - ``side`` -- string; one of ``'twosided'`` (default),
          ``'left'`` or ``'right'``. Determines whether the ideal
          is a left, right or twosided ideal. Groebner bases or
          only supported in the twosided case.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: from sage.algebras.letterplace.letterplace_ideal import LetterplaceIdeal
            sage: LetterplaceIdeal(F,x)
            Twosided Ideal (x) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
            sage: LetterplaceIdeal(F,[x,y],side='left')
            Left Ideal (x, y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field

        It is not correctly detected that this class inherits from an
        extension class. Therefore, we have to skip one item of the
        test suite::

            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: TestSuite(I).run(skip=['_test_category'],verbose=True)
            running ._test_eq() . . . pass
            running ._test_new() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass
        """
        Ideal_nc.__init__(self, ring, gens, coerce=coerce, side=side)
        self.__GB = self
        self.__uptodeg = 0

    def groebner_basis(self, degbound=None):
        """
        Twosided Groebner basis with degree bound.

        INPUT:

        - ``degbound`` -- (optional) integer or Infinity; if it is provided,
          a Groebner basis at least out to that degree is returned. By
          default, the current degree bound of the underlying ring is used.

        ASSUMPTIONS:

        Currently, we can only compute Groebner bases for twosided
        ideals, and the ring of coefficients must be a field. A
        `TypeError` is raised if one of these conditions is violated.

        .. NOTE::

            - The result is cached. The same Groebner basis is returned
              if a smaller degree bound than the known one is requested.
            - If the degree bound ``Infinity`` is requested, it is attempted to
              compute a complete Groebner basis. But we cannot guarantee
              that the computation will terminate, since not all twosided
              homogeneous ideals of a free algebra have a finite Groebner
              basis.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F

        Since `F` was cached and since its degree bound cannot be
        decreased, it may happen that, as a side effect of other tests,
        it already has a degree bound bigger than 3. So, we cannot
        test against the output of ``I.groebner_basis()``::

            sage: F.set_degbound(3)
            sage: I.groebner_basis()   # not tested
            Twosided Ideal (y*y*y - y*y*z + y*z*y - y*z*z, y*y*x + y*y*z + y*z*x + y*z*z, x*y + y*z, x*x - y*x - y*y - y*z) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
            sage: I.groebner_basis(4)
            Twosided Ideal (x*y + y*z,
                x*x - y*x - y*y - y*z,
                y*y*y - y*y*z + y*z*y - y*z*z,
                y*y*x + y*y*z + y*z*x + y*z*z,
                y*y*z*y - y*y*z*z + y*z*z*y - y*z*z*z,
                y*z*y*y - y*z*y*z + y*z*z*y - y*z*z*z,
                y*y*z*x + y*y*z*z + y*z*z*x + y*z*z*z,
                y*z*y*x + y*z*y*z + y*z*z*x + y*z*z*z) of Free Associative
                Unital Algebra on 3 generators (x, y, z) over Rational Field
            sage: I.groebner_basis(2) is I.groebner_basis(4)
            True
            sage: G = I.groebner_basis(4)
            sage: G.groebner_basis(3) is G
            True

        If a finite complete Groebner basis exists, we can compute
        it as follows::

            sage: I = F*[x*y-y*x,x*z-z*x,y*z-z*y,x^2*y-z^3,x*y^2+z*x^2]*F
            sage: I.groebner_basis(Infinity)
            Twosided Ideal (-y*z + z*y,
                -x*z + z*x,
                -x*y + y*x,
                x*x*z + x*y*y,
                x*x*y - z*z*z,
                x*x*x*z + y*z*z*z,
                x*z*z*z*z + y*y*z*z*z) of Free Associative Unital Algebra
                on 3 generators (x, y, z) over Rational Field

        Since the commutators of the generators are contained in the ideal,
        we can verify the above result by a computation in a polynomial ring
        in negative lexicographic order::

            sage: P.<c,b,a> = PolynomialRing(QQ,order='neglex')
            sage: J = P*[a^2*b-c^3,a*b^2+c*a^2]
            sage: J.groebner_basis()
            [b*a^2 - c^3, b^2*a + c*a^2, c*a^3 + c^3*b, c^3*b^2 + c^4*a]

        Apparently, the results are compatible, by sending `a` to `x`, `b`
        to `y` and `c` to `z`.
        """
        cdef FreeAlgebra_letterplace A = self.ring()
        cdef FreeAlgebraElement_letterplace x
        if degbound is None:
            degbound = A.degbound()
        if self.__uptodeg >= degbound:
            return self.__GB
        if not A.base().is_field():
            raise TypeError("Currently, we can only compute Groebner bases if the ring of coefficients is a field")
        if self.side() != 'twosided':
            raise TypeError("This ideal is not two-sided. We can only compute two-sided Groebner bases")
        if degbound == Infinity:
            while self.__uptodeg < Infinity:
                test_bound = 2 * max([x._poly.degree()
                                      for x in self.__GB.gens()])
                self.groebner_basis(test_bound)
            return self.__GB
        # Set the options required by letterplace
        from sage.libs.singular.option import LibSingularOptions
        libsingular_options = LibSingularOptions()
        bck = (libsingular_options['redTail'], libsingular_options['redSB'])
        libsingular_options['redTail'] = True
        libsingular_options['redSB'] = True
        A.set_degbound(degbound)
        P = A._current_ring

        # note that degbound might be smaller than A._degbound due to caching,
        # but degbound must be large enough to map all generators to the
        # letterplace ring L
        if degbound < A._degbound:
            max_deg = max([x._poly.degree() for x in self.__GB.gens()])
            if degbound < max_deg:
                degbound = max_deg

        # The following is a workaround for calling Singular's new Letterplace
        # API (see :issue:`25993`). We construct a temporary polynomial ring L
        # with letterplace attributes set as required by the API. As L has
        # duplicate variable names, we need to handle this ring carefully; in
        # particular, we cannot coerce to and from L, so we use homomorphisms
        # for the conversion.

        cdef FreeAlgebra_letterplace_libsingular lp_ring = \
            FreeAlgebra_letterplace_libsingular(A._commutative_ring, degbound)
        L = lp_ring._lp_ring_internal
        to_L = P.hom(L.gens(), L, check=False)
        from_L = L.hom(P.gens(), P, check=False)
        I = L.ideal([to_L(x._poly) for x in self.__GB.gens()])
        gb = singular_twostd(I)
        out = [FreeAlgebraElement_letterplace(A, from_L(X), check=False)
               for X in gb]

        libsingular_options['redTail'] = bck[0]
        libsingular_options['redSB'] = bck[1]
        self.__GB = A.ideal(out, side='twosided', coerce=False)
        if degbound >= 2 * max([x._poly.degree() for x in out]):
            degbound = Infinity
        self.__uptodeg = degbound
        self.__GB.__uptodeg = degbound
        return self.__GB

    def __contains__(self, x):
        """
        The containment test is based on a normal form computation.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: x*I.0-I.1*y+I.0*y in I    # indirect doctest
            True
            sage: 1 in I
            False
        """
        R = self.ring()
        return (x in R) and R(x).normal_form(self).is_zero()

    def reduce(self, G):
        """
        Reduction of this ideal by another ideal,
        or normal form of an algebra element with respect to this ideal.

        INPUT:

        - ``G`` -- a list or tuple of elements, an ideal,
          the ambient algebra, or a single element

        OUTPUT:

        - The normal form of ``G`` with respect to this ideal, if
          ``G`` is an element of the algebra.
        - The reduction of this ideal by the elements resp. generators
          of ``G``, if ``G`` is a list, tuple or ideal.
        - The zero ideal, if ``G`` is the algebra containing this ideal.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: I.reduce(F)
            Twosided Ideal (0) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
            sage: I.reduce(x^3)
            -y*z*x - y*z*y - y*z*z
            sage: I.reduce([x*y])
            Twosided Ideal (y*z, x*x - y*x - y*y) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
            sage: I.reduce(F*[x^2+x*y,y^2+y*z]*F)
            Twosided Ideal (x*y + y*z, -y*x + y*z) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
        """
        P = self.ring()
        if not isinstance(G, (list, tuple)):
            if G == P:
                return P.ideal([P.zero()])
            if G in P:
                return G.normal_form(self)
            G = G.gens()
        C = P.current_ring()
        sI = C.ideal([C(X.letterplace_polynomial()) for X in self.gens()],
                     coerce=False)
        selfdeg = max([x.degree() for x in sI.gens()])
        gI = P._reductor_(G, selfdeg)
        from sage.libs.singular.option import LibSingularOptions
        libsingular_options = LibSingularOptions()
        bck = (libsingular_options['redTail'], libsingular_options['redSB'])
        libsingular_options['redTail'] = True
        libsingular_options['redSB'] = True
        sI = poly_reduce(sI, gI, ring=C, attributes={gI: {"isSB": 1}})
        libsingular_options['redTail'] = bck[0]
        libsingular_options['redSB'] = bck[1]
        return P.ideal([FreeAlgebraElement_letterplace(P, x, check=False)
                        for x in sI], coerce=False)
