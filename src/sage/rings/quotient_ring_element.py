"""
Quotient Ring Elements

AUTHORS:

- William Stein
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import RingElement
from sage.structure.richcmp import richcmp, rich_to_bool


try:
    from sage.interfaces.singular import singular as singular_default
except ImportError:
    singular_default = None


class QuotientRingElement(RingElement):
    """
    An element of a quotient ring `R/I`.

    INPUT:

    - ``parent`` - the ring `R/I`

    - ``rep`` - a representative of the element in `R`; this is used
      as the internal representation of the element

    - ``reduce`` - bool (optional, default: True) - if True, then the
      internal representation of the element is ``rep`` reduced modulo
      the ideal `I`

    EXAMPLES::

        sage: R.<x> = PolynomialRing(ZZ)
        sage: S.<xbar> = R.quo((4 + 3*x + x^2, 1 + x^2)); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring
         by the ideal (x^2 + 3*x + 4, x^2 + 1)
        sage: v = S.gens(); v
        (xbar,)

    ::

        sage: loads(v[0].dumps()) == v[0]
        True

    ::

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S = R.quo(x^2 + y^2); S
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field
         by the ideal (x^2 + y^2)
        sage: S.gens()                                                                  # optional - sage.libs.singular
        (xbar, ybar)

    We name each of the generators.

    ::

        sage: S.<a,b> = R.quotient(x^2 + y^2)                                           # optional - sage.libs.singular
        sage: a                                                                         # optional - sage.libs.singular
        a
        sage: b                                                                         # optional - sage.libs.singular
        b
        sage: a^2 + b^2 == 0                                                            # optional - sage.libs.singular
        True
        sage: b.lift()                                                                  # optional - sage.libs.singular
        y
        sage: (a^3 + b^2).lift()                                                        # optional - sage.libs.singular
        -x*y^2 + y^2
    """
    def __init__(self, parent, rep, reduce=True):
        """
        An element of a quotient ring `R/I`.  See
        ``QuotientRingElement`` for full documentation.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<xbar> = R.quo((4 + 3*x + x^2, 1 + x^2)); S
            Quotient of Univariate Polynomial Ring in x over Integer Ring
             by the ideal (x^2 + 3*x + 4, x^2 + 1)
            sage: v = S.gens(); v
            (xbar,)
        """
        RingElement.__init__(self, parent)
        self.__rep = rep
        if reduce:
            self._reduce_()

    def _reduce_(self):
        """
        Reduce the element modulo the defining ideal of the quotient
        ring.  This internal method replaces the cached representative
        by one in reduced form.

        (Note that this has nothing to do with pickling.)

        TESTS::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a._reduce_()                                                          # optional - sage.libs.singular
            sage: a._QuotientRingElement__rep                                           # optional - sage.libs.singular
            x
        """
        I = self.parent().defining_ideal()
        self.__rep = I.reduce(self.__rep)

    def lift(self):
        """
        If self is an element of `R/I`, then return self as an
        element of `R`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a.lift()                                                              # optional - sage.libs.singular
            x
            sage: (3/5*(a + a^2 + b^2)).lift()                                          # optional - sage.libs.singular
            3/5*x
        """
        return self.__rep

    def __bool__(self):
        """
        Return ``True`` if quotient ring element is non-zero in the
        quotient ring `R/I`, by determining whether the element
        is in `I`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: bool(a)     # indirect doctest                                        # optional - sage.libs.singular
            True
            sage: bool(S(0))                                                            # optional - sage.libs.singular
            False

        TESTS::

            sage: bool(a - a)                                                           # optional - sage.libs.singular
            False
        """
        return self.__rep not in self.parent().defining_ideal()

    def is_unit(self):
        """
        Return ``True`` if self is a unit in the quotient ring.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(1 - x*y); type(a)                     # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a*b                                                                   # optional - sage.libs.singular
            1
            sage: S(2).is_unit()                                                        # optional - sage.libs.singular
            True

        Check that :trac:`29469` is fixed::

            sage: a.is_unit()                                                           # optional - sage.libs.singular
            True
            sage: (a+b).is_unit()                                                       # optional - sage.libs.singular
            False
        """
        if self.__rep.is_unit():
            return True
        from sage.categories.fields import Fields
        if self.parent() in Fields():
            return not self.is_zero()
        try:
            self.__invert__()
            return True
        except ArithmeticError:
            return False
        raise NotImplementedError

    def _repr_(self):
        """
        String representation.

        TESTS::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a-2*a*b     # indirect doctest                                        # optional - sage.libs.singular
            -2*a*b + a

        In :trac:`11068`, the case of quotient rings without
        assigned names has been covered as well::

            sage: S = SteenrodAlgebra(2)                                                # optional - sage.libs.singular
            sage: I = S * [S.0 + S.1] * S                                               # optional - sage.libs.singular
            sage: Q = S.quo(I)                                                          # optional - sage.libs.singular
            sage: Q.0                                                                   # optional - sage.libs.singular
            Sq(1)

        """
        from sage.structure.parent_gens import localvars
        P = self.parent()
        R = P.cover_ring()
        # We print by temporarily (and safely!) changing the variable
        # names of the covering structure R to those of P.
        # These names get changed back, since we're using "with".
        # However, it may occur that no variable names are assigned.
        # That holds, in particular, if there are infinitely many
        # generators, as for Steenrod algebras.
        try:
            P.variable_names()
        except ValueError:
            return str(self.__rep)
        with localvars(R, P.variable_names(), normalize=False):
            return str(self.__rep)

    def _latex_(self):
        """
        Return the LaTeX representation as a string.

        EXAMPLES::

            sage: R = PolynomialRing(QQ, 'a, b, c')
            sage: a = R.gen(0)
            sage: I = R.ideal(a**2 + a + 1)
            sage: S = R.quotient(I, names=R.variable_names())
            sage: a = S.gen(0)                                                          # optional - sage.libs.singular
            sage: latex(a)                                                              # optional - sage.libs.singular
            a
        """
        from sage.structure.parent_gens import localvars
        P = self.parent()
        R = P.cover_ring()
        # see _repr_ above for the idea
        try:
            P.variable_names()
        except ValueError:
            return self.__rep._latex_()
        with localvars(R, P.variable_names(), normalize=False):
            return self.__rep._latex_()

    def __pari__(self):
        """
        The Pari representation of this quotient element.

        Since Pari does not support quotients by non-principal ideals,
        this function will raise an error in that case.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal(x^3, y^3)
            sage: S.<xb,yb> = R.quo(I)                                                  # optional - sage.libs.singular
            sage: pari(xb)                                                              # optional - sage.libs.pari sage.libs.singular
            Traceback (most recent call last):
            ...
            ValueError: Pari does not support quotients by non-principal ideals

        Note that the quotient does work in the case that the ideal is principal::

            sage: I = R.ideal(x^3 + y^3)
            sage: S.<xb,yb> = R.quo(I)                                                  # optional - sage.libs.singular
            sage: pari(xb)^4                                                            # optional - sage.libs.pari sage.libs.singular
            Mod(-y^3*x, x^3 + y^3)
            sage: pari(yb)^4                                                            # optional - sage.libs.pari sage.libs.singular
            Mod(y^4, x^3 + y^3)
        """
        gens = self.parent().defining_ideal().gens()
        if len(gens) != 1:
            raise ValueError("Pari does not support quotients by non-principal ideals")
        return self.__rep.__pari__().Mod(gens[0])

    def _add_(self, right):
        """
        Add quotient ring element ``self`` to another quotient ring
        element, ``right``. If the quotient is `R/I`, the addition is
        carried out in `R` and then reduced to `R/I`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a + b                                                                 # optional - sage.libs.singular
            a + b

        TESTS::

            sage: a._add_(b)                                                            # optional - sage.libs.singular
            a + b
        """
        return self.__class__(self.parent(), self.__rep + right.__rep)

    def _sub_(self, right):
        """
        Subtract quotient ring element ``right`` from quotient ring
        element ``self``. If the quotient is `R/I`, the subtraction is
        carried out in `R` and then reduced to `R/I`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a - b                                                                 # optional - sage.libs.singular
            a - b

        TESTS::

            sage: a._sub_(b)                                                            # optional - sage.libs.singular
            a - b
        """
        return self.__class__(self.parent(), self.__rep - right.__rep)

    def _mul_(self, right):
        """
        Multiply quotient ring element ``self`` by another quotient ring
        element, ``right``. If the quotient is `R/I`, the multiplication is
        carried out in `R` and then reduced to `R/I`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a * b                                                                 # optional - sage.libs.singular
            a*b

        TESTS::

            sage: a._mul_(b)                                                            # optional - sage.libs.singular
            a*b
            sage: a._mul_(a)                                                            # optional - sage.libs.singular
            -b^2
        """
        return self.__class__(self.parent(), self.__rep * right.__rep)

    def _div_(self, right):
        """
        Divide quotient ring element ``self`` by another quotient ring
        element, ``right``. If the quotient is `R/I`, the division is
        carried out in `R` and then reduced to `R/I`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([x^2 + 1, y^3 - 2])
            sage: S.<i,cuberoot> = R.quotient(I)                                        # optional - sage.libs.singular
            sage: 1/(1+i)                                                               # optional - sage.libs.singular
            -1/2*i + 1/2

        Confirm via symbolic computation::

            sage: 1/(1+sqrt(-1))                                                        # optional - sage.symbolic
            -1/2*I + 1/2

        Another more complicated quotient::

            sage: b = 1/(i+cuberoot); b                                                 # optional - sage.libs.singular
            1/5*i*cuberoot^2 - 2/5*i*cuberoot + 2/5*cuberoot^2 - 1/5*i + 1/5*cuberoot - 2/5
            sage: b*(i+cuberoot)                                                        # optional - sage.libs.singular
            1


        Another really easy example::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a / S(2)                                                              # optional - sage.libs.singular
            1/2*a
            sage: (a*b)._div_(b)                                                        # optional - sage.libs.singular
            a

        An example in which we try to divide in a ring that is not a
        field::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([x^2 - 1, y^3 - 2])
            sage: S.<a,cuberoot> = R.quotient(I)                                        # optional - sage.libs.singular
            sage: 1/cuberoot                                                            # optional - sage.libs.singular
            1/2*cuberoot^2
            sage: 1/a                                                                   # optional - sage.libs.singular
            a

        Check that :trac:`13670` is fixed (i.e. that the error message
        actually describes what happens when the result of division is not defined)::

            sage: R.<x1,x2> = QQ[]
            sage: S = R.quotient_ring( R.ideal(x2**2 + x1 - 2, x1**2 - 1) )             # optional - sage.libs.singular
            sage: 1 / S(x1 + x2)                                                        # optional - sage.libs.singular
            Traceback (most recent call last):
            ...
            ArithmeticError: Division failed. The numerator is not a multiple of the denominator.

        An example over a polynomial ring over a polynomial ring,
        which doesn't work (yet; obviously, this could be made to work
        by converting to a single polynomial quotient ring
        internally)::

            sage: R.<x> = QQ[]
            sage: S.<y,z> = R[]
            sage: Z.<ybar,zbar> = S.quotient([y^2 - 2, z^2 - 3])                        # optional - sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: Can only reduce polynomials over fields.
        """
        # Special case: if self==0 (and right is nonzero), just return self.
        if not self:
            if not right:
                raise ZeroDivisionError
            return self

        # We are computing L/R modulo the ideal.
        (L, R) = (self.__rep, right.__rep)
        P  = self.parent()
        I = P.defining_ideal()

        if not hasattr(I, 'groebner_basis'):
            # Try something very naive -- somebody will improve this
            # in the future.
            try:
                return L * R.inverse_mod(I)
            except NotImplementedError:
                if R.is_unit():
                    return L * ~R
                else:
                    raise

        # Now the parent is a polynomial ring, so we have an algorithm
        # at our disposal.

        # Our algorithm is to write L in terms of R and a Groebner
        # basis for the defining ideal.  We compute a Groebner basis
        # here explicitly purely for efficiency reasons, since it
        # makes the implicit Groebner basis computation of [R]+B
        # that is done in the lift command below faster.

        B = I.groebner_basis()
        try:
            XY = L.lift((R,) + tuple(B))
        except ValueError:
            raise ArithmeticError("Division failed. The numerator is not "
                                  "a multiple of the denominator.")
        return P(XY[0])

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` in ``codomain`` under the map
        that sends ``self.parent().gens()`` to ``im_gens``.

        INPUT:

        - ``codomain`` -- a ring

        - ``im_gens`` -- a tuple of elements `f(x)` in ``codomain``,
          one for each `x` in ``self.parent().gens()``, that define
          a homomorphism `f` from ``self.parent()`` to ``codomain``

        OUTPUT:

        The image of ``self`` in ``codomain`` under the above
        homomorphism `f`.

        EXAMPLES:

        Ring homomorphisms whose domain is the fraction field of a
        quotient ring work correctly (see :trac:`16135`)::

            sage: R.<x, y> = QQ[]
            sage: K = R.quotient(x^2 - y^3).fraction_field()                            # optional - sage.libs.singular
            sage: L.<t> = FunctionField(QQ)
            sage: f = K.hom((t^3, t^2))                                                 # optional - sage.libs.singular
            sage: list(map(f, K.gens()))                                                # optional - sage.libs.singular
            [t^3, t^2]
            sage: xbar, ybar = K.gens()                                                 # optional - sage.libs.singular
            sage: f(1/ybar)                                                             # optional - sage.libs.singular
            1/t^2
            sage: f(xbar/ybar)                                                          # optional - sage.libs.singular
            t
        """
        return self.lift()._im_gens_(codomain, im_gens, base_map=base_map)

    def __int__(self):
        """
        Try to convert self (an element of `R/I`) to an integer by
        converting its lift in `R` to an integer.  Return a TypeError
        if no such conversion can be found.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: int(S(-3))                # indirect doctest                          # optional - sage.libs.singular
            -3
            sage: type(int(S(-3)))                                                      # optional - sage.libs.singular
            <... 'int'>
            sage: int(a)                                                                # optional - sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to <class 'int'>
        """
        return int(self.lift())

    def _integer_(self, Z):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: ZZ(S(-3))                                                             # optional - sage.libs.singular
            -3

        TESTS::

            sage: type(ZZ(S(-3)))                                                       # optional - sage.libs.singular
            <class 'sage.rings.integer.Integer'>
        """
        return Z(self.lift())

    def _rational_(self):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: QQ(S(-2/3))                                                           # optional - sage.libs.singular
            -2/3

        TESTS::

            sage: type(S(-2/3)._rational_())                                            # optional - sage.libs.singular
            <class 'sage.rings.rational.Rational'>
        """
        from sage.rings.rational_field import QQ
        return QQ(self.lift())

    def __neg__(self):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: -a                     # indirect doctest                             # optional - sage.libs.singular
            -a
            sage: -(a+b)                                                                # optional - sage.libs.singular
            -a - b
        """
        return self.__class__(self.parent(), -self.__rep)

    def __pos__(self):
        """
        TESTS::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: (a+b).__pos__()                                                       # optional - sage.libs.singular
            a + b
            sage: c = a+b; c.__pos__() is c                                             # optional - sage.libs.singular
            True
        """
        return self

    def __invert__(self):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: ~S(2/3)                                                               # optional - sage.libs.singular
            3/2

        TESTS::

            sage: S(2/3).__invert__()                                                   # optional - sage.libs.singular
            3/2

        Note that a is not invertible as an element of R::

            sage: a.__invert__()                                                        # optional - sage.libs.singular
            Traceback (most recent call last):
            ...
            ArithmeticError: element is non-invertible
        """
        try:
            inv = self.__rep.inverse_mod(self.parent().defining_ideal())
        except NotImplementedError:
            return self.parent().one()/self
        return self.__class__(self.parent(), inv)

    def __float__(self):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: float(S(2/3))                                                         # optional - sage.libs.singular
            0.6666666666666666
            sage: float(a)                                                              # optional - sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: unable to convert non-constant polynomial x to <class 'float'>
        """
        return float(self.lift())

    def __hash__(self):
        r"""
        TESTS::

            sage: R.<x,y> = QQ[]                                                        # optional - sage.libs.singular
            sage: S.<a,b> = R.quo(x^2 + y^2)                                            # optional - sage.libs.singular
            sage: c = a*a + b                                                           # optional - sage.libs.singular
            sage: hash(a) != hash(b)                                                    # optional - sage.libs.singular
            True
        """
        return hash(self.__rep)

    def _richcmp_(self, other, op):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a > b    # indirect doctest                                           # optional - sage.libs.singular
            True
            sage: b > a                                                                 # optional - sage.libs.singular
            False
            sage: a == loads(dumps(a))                                                  # optional - sage.libs.singular
            True

        TESTS::

            sage: a == (a+1-1)                                                          # optional - sage.libs.singular
            True
            sage: a > b                                                                 # optional - sage.libs.singular
            True

        See :trac:`7797`::

            sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')             # optional - sage.combinat
            sage: I = F * [x*y + y*z, x^2 + x*y - y*x - y^2] * F                        # optional - sage.combinat
            sage: Q = F.quo(I)                                                          # optional - sage.combinat sage.libs.singular
            sage: Q.0^4    # indirect doctest                                           # optional - sage.combinat sage.libs.singular
            ybar*zbar*zbar*xbar + ybar*zbar*zbar*ybar + ybar*zbar*zbar*zbar

        The issue from :trac:`8005` was most likely fixed as part of
        :trac:`9138`::

            sage: F = GF(5)                                                             # optional - sage.rings.finite_rings
            sage: R.<x,y> = F[]                                                         # optional - sage.rings.finite_rings
            sage: I = Ideal(R, [x, y])                                                  # optional - sage.rings.finite_rings
            sage: S.<x1,y1> = QuotientRing(R, I)                                        # optional - sage.rings.finite_rings
            sage: x1^4                                                                  # optional - sage.rings.finite_rings
            0
        """
        # A containment test is not implemented for univariate polynomial
        # ideals. There are cases in which one would not like to add
        # elements of different degrees. The whole quotient stuff relies
        # in I.reduce(x) returning a normal form of x with respect to I.
        # Hence, we will not use more than that.

        # Since we have to compute normal forms anyway, it makes sense
        # to use it for comparison in the case of an inequality as well.
        if self.__rep == other.__rep:
            # Use a shortpath, so that we avoid expensive reductions
            return rich_to_bool(op, 0)
        I = self.parent().defining_ideal()
        return richcmp(I.reduce(self.__rep), I.reduce(other.__rep), op)

    def lt(self):
        """
        Return the leading term of this quotient ring element.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(GF(7), 3, order='lex')                     # optional - sage.rings.finite_rings
            sage: I = sage.rings.ideal.FieldIdeal(R)                                    # optional - sage.rings.finite_rings
            sage: Q = R.quo(I)                                                          # optional - sage.rings.finite_rings
            sage: f = Q(z*y + 2*x)                                                      # optional - sage.rings.finite_rings
            sage: f.lt()                                                                # optional - sage.rings.finite_rings
            2*xbar

        TESTS::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: (a + 3*a*b + b).lt()                                                  # optional - sage.libs.singular
            3*a*b
        """
        return self.__class__(self.parent(), self.__rep.lt())

    def lm(self):
        """
        Return the leading monomial of this quotient ring element.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(GF(7), 3, order='lex')                     # optional - sage.rings.finite_rings
            sage: I = sage.rings.ideal.FieldIdeal(R)                                    # optional - sage.rings.finite_rings
            sage: Q = R.quo(I)                                                          # optional - sage.rings.finite_rings
            sage: f = Q(z*y + 2*x)                                                      # optional - sage.rings.finite_rings
            sage: f.lm()                                                                # optional - sage.rings.finite_rings
            xbar

        TESTS::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: (a+3*a*b+b).lm()                                                      # optional - sage.libs.singular
            a*b

        """
        return self.__class__(self.parent(), self.__rep.lm())

    def lc(self):
        """
        Return the leading coefficient of this quotient ring element.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(GF(7), 3, order='lex')                     # optional - sage.rings.finite_rings
            sage: I = sage.rings.ideal.FieldIdeal(R)                                    # optional - sage.rings.finite_rings
            sage: Q = R.quo(I)                                                          # optional - sage.rings.finite_rings
            sage: f = Q(z*y + 2*x)                                                      # optional - sage.rings.finite_rings
            sage: f.lc()                                                                # optional - sage.rings.finite_rings
            2

        TESTS::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: (a + 3*a*b + b).lc()                                                  # optional - sage.libs.singular
            3
        """
        return self.__rep.lc()

    def variables(self):
        """
        Return all variables occurring in ``self``.

        OUTPUT:

        A tuple of linear monomials, one for each variable occurring
        in ``self``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a.variables()                                                         # optional - sage.libs.singular
            (a,)
            sage: b.variables()                                                         # optional - sage.libs.singular
            (b,)
            sage: s = a^2 + b^2 + 1; s                                                  # optional - sage.libs.singular
            1
            sage: s.variables()                                                         # optional - sage.libs.singular
            ()
            sage: (a + b).variables()                                                   # optional - sage.libs.singular
            (a, b)
        """
        return tuple(self.__class__(self.parent(), v) for v in self.__rep.variables())

    def monomials(self):
        """
        Return the monomials in ``self``.

        OUTPUT:

        A list of monomials.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: a.monomials()                                                         # optional - sage.libs.singular
            [a]
            sage: (a + a*b).monomials()                                                 # optional - sage.libs.singular
            [a*b, a]
            sage: R.zero().monomials()                                                  # optional - sage.libs.singular
            []
        """
        return [self.__class__(self.parent(), m) for m in self.__rep.monomials()]

    def _singular_(self, singular=singular_default):
        """
        Return Singular representation of self.

        INPUT:

        -  ``singular`` - a non-standard interpreter may be
           provided

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(GF(2), 2)                                    # optional - sage.rings.finite_rings
            sage: I = sage.rings.ideal.FieldIdeal(P)                                    # optional - sage.rings.finite_rings
            sage: Q = P.quo(I)                                                          # optional - sage.rings.finite_rings
            sage: Q._singular_()                                                        # optional - sage.rings.finite_rings
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/2
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            // quotient ring from ideal
            _[1]=x2+x
            _[2]=y2+y
            sage: xbar = Q(x); xbar                                                     # optional - sage.rings.finite_rings
            xbar
            sage: xbar._singular_()                                                     # optional - sage.rings.finite_rings
            x
            sage: Q(xbar._singular_()) # a round-trip                                   # optional - sage.rings.finite_rings
            xbar

        TESTS::

            sage: R.<x,y> = QQ[]; S.<a,b> = R.quo(x^2 + y^2); type(a)                   # optional - sage.libs.singular
            <class 'sage.rings.quotient_ring.QuotientRing_generic_with_category.element_class'>
            sage: (a - 2/3*b)._singular_()                                              # optional - sage.libs.singular
            x-2/3*y
            sage: S((a - 2/3*b)._singular_())                                           # optional - sage.libs.singular
            a - 2/3*b
        """
        if singular is None:
            raise ImportError("could not import singular")
        return self.__rep._singular_(singular)

    def _magma_init_(self, magma):
        """
        Returns the Magma representation of this quotient ring element.

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(GF(2))                                       # optional - sage.rings.finite_rings
            sage: Q = P.quotient(sage.rings.ideal.FieldIdeal(P))                        # optional - sage.rings.finite_rings
            sage: xbar, ybar = Q.gens()                                                 # optional - sage.rings.finite_rings
            sage: magma(xbar)                       # optional - magma                  # optional - sage.rings.finite_rings
            x
            sage: xbar._magma_init_(magma)          # optional - magma                  # optional - sage.rings.finite_rings
            '_sage_[...]!_sage_ref...'
        """
        g = magma(self.__rep)
        R = magma(self.parent())
        return '{}!{}'.format(R.name(), g._ref())

    def _macaulay2_(self, macaulay2=None):
        """
        The Macaulay2 element corresponding to this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(GF(7), 2)                                    # optional - sage.rings.finite_rings
            sage: Q = R.quotient([x^2 - y])                                             # optional - sage.rings.finite_rings
            sage: x, y = Q.gens()                                                       # optional - sage.rings.finite_rings
            sage: f = (x^3 + 2*y^2*x)^7; f                                              # optional - sage.rings.finite_rings
            2*xbar*ybar^17 + xbar*ybar^10
            sage: mf = macaulay2(f); mf             # optional - macaulay2              # optional - sage.rings.finite_rings
                17      10
            2x*y   + x*y
            sage: mf.sage()                         # optional - macaulay2              # optional - sage.rings.finite_rings
            2*x*y^17 + x*y^10
            sage: mf.sage() == f                    # optional - macaulay2              # optional - sage.rings.finite_rings
            True
            sage: Q(mf)                             # optional - macaulay2              # optional - sage.rings.finite_rings
            2*xbar*ybar^17 + xbar*ybar^10

        In Macaulay2, the variable names for a quotient ring are inherited from
        the variable names of the ambient ring. This is in contrast to Sage's
        default behaviour in which the variable names for the quotient ring are
        obtained by appending ``bar`` to the variable names of the ambient
        ring. This can be controlled using the ``names`` argument of the
        ``quotient`` method.

        ::

            sage: R.<x,y> = PolynomialRing(GF(7), 2)                                    # optional - sage.rings.finite_rings
            sage: Q = R.quotient([x^2 - y], names=R.gens())                             # optional - sage.rings.finite_rings
            sage: x, y = Q.gens()                                                       # optional - sage.rings.finite_rings
            sage: f = (x^3 + 2*y^2*x)^7; f                                              # optional - sage.rings.finite_rings
            2*x*y^17 + x*y^10
            sage: macaulay2(f)                      # optional - macaulay2              # optional - sage.rings.finite_rings
                17      10
            2x*y   + x*y
            sage: _.sage()                          # optional - macaulay2              # optional - sage.rings.finite_rings
            2*x*y^17 + x*y^10

        TESTS:

        Check that changing the currently defined global variables (`x`, `y`,
        ...) in Macaulay2 does not affect the result of this conversion::

            sage: R.<x,y> = PolynomialRing(GF(7), 2)                                    # optional - sage.rings.finite_rings
            sage: Q = R.quotient([x^2 - y], names=R.gens())                             # optional - sage.rings.finite_rings
            sage: x, y = Q.gens()                                                       # optional - sage.rings.finite_rings
            sage: f = (x^3 + 2*y^2*x)^7                                                 # optional - sage.rings.finite_rings
            sage: macaulay2(f)                      # optional - macaulay2              # optional - sage.rings.finite_rings
                17      10
            2x*y   + x*y
            sage: macaulay2.use(R.quotient([x, y])) # optional - macaulay2              # optional - sage.rings.finite_rings
            sage: macaulay2(f)                      # optional - macaulay2              # optional - sage.rings.finite_rings
                17      10
            2x*y   + x*y
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default
        m2_parent = self.parent()._macaulay2_(macaulay2)
        macaulay2.use(m2_parent)
        return macaulay2.substitute(repr(self.lift()), m2_parent)

    def reduce(self, G):
        r"""
        Reduce this quotient ring element by a set of quotient ring
        elements ``G``.

        INPUT:

        -  ``G`` - a list of quotient ring elements

        .. WARNING::

            This method is not guaranteed to return unique minimal results.
            For quotients of polynomial rings, use
            :meth:`~sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal.reduce`
            on the ideal generated by ``G``, instead.

        EXAMPLES::

            sage: P.<a,b,c,d,e> = PolynomialRing(GF(2), 5, order='lex')                 # optional - sage.rings.finite_rings
            sage: I1 = ideal([a*b + c*d + 1, a*c*e + d*e,                               # optional - sage.rings.finite_rings
            ....:             a*b*e + c*e, b*c + c*d*e + 1])
            sage: Q = P.quotient(sage.rings.ideal.FieldIdeal(P))                        # optional - sage.rings.finite_rings
            sage: I2 = ideal([Q(f) for f in I1.gens()])                                 # optional - sage.rings.finite_rings
            sage: f = Q((a*b + c*d + 1)^2  + e)                                         # optional - sage.rings.finite_rings
            sage: f.reduce(I2.gens())                                                   # optional - sage.rings.finite_rings
            ebar

        Notice that the result above is not minimal::

            sage: I2.reduce(f)                                                          # optional - sage.rings.finite_rings
            0
        """
        try:
            G = [f.lift() for f in G]
        except TypeError:
            pass
        # reduction w.r.t. the defining ideal is performed in the
        # constructor
        return self.__class__(self.parent(), self.__rep.reduce(G))
