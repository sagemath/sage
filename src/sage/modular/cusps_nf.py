# sage.doctest: needs sage.rings.number_field
r"""
The set `\mathbb{P}^1(K)` of cusps of a number field `K`

AUTHORS:

- Maite Aranes (2009): Initial version

EXAMPLES:

The space of cusps over a number field k::

    sage: x = polygen(ZZ, 'x')
    sage: k.<a> = NumberField(x^2 + 5)
    sage: kCusps = NFCusps(k); kCusps
    Set of all cusps of Number Field in a with defining polynomial x^2 + 5
    sage: kCusps is NFCusps(k)
    True

Define a cusp over a number field::

    sage: NFCusp(k, a, 2/(a+1))
    Cusp [a - 5: 2] of Number Field in a with defining polynomial x^2 + 5
    sage: kCusps((a,2))
    Cusp [a: 2] of Number Field in a with defining polynomial x^2 + 5
    sage: NFCusp(k,oo)
    Cusp Infinity of Number Field in a with defining polynomial x^2 + 5

Different operations with cusps over a number field::

    sage: alpha = NFCusp(k, 3, 1/a + 2); alpha
    Cusp [a + 10: 7] of Number Field in a with defining polynomial x^2 + 5
    sage: alpha.numerator()
    a + 10
    sage: alpha.denominator()
    7
    sage: alpha.ideal()
    Fractional ideal (7, a + 3)
    sage: M = alpha.ABmatrix(); M # random
    [a + 10, 2*a + 6, 7, a + 5]
    sage: NFCusp(k, oo).apply(M)
    Cusp [a + 10: 7] of Number Field in a with defining polynomial x^2 + 5

Check Gamma0(N)-equivalence of cusps::

    sage: N = k.ideal(3)
    sage: alpha = NFCusp(k, 3, a + 1)
    sage: beta = kCusps((2, a - 3))
    sage: alpha.is_Gamma0_equivalent(beta, N)
    True

Obtain transformation matrix for equivalent cusps::

    sage: t, M = alpha.is_Gamma0_equivalent(beta, N, Transformation=True)
    sage: M[2] in N
    True
    sage: M[0]*M[3] - M[1]*M[2] == 1
    True
    sage: alpha.apply(M) == beta
    True

List representatives for Gamma_0(N) - equivalence classes of cusps::

    sage: Gamma0_NFCusps(N)
    [Cusp [0: 1] of Number Field in a with defining polynomial x^2 + 5,
     Cusp [1: 3] of Number Field in a with defining polynomial x^2 + 5,
    ...]
"""
# ****************************************************************************
#       Copyright (C) 2009, Maite Aranes <M.T.Aranes@warwick.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import Element, InfinityElement
from sage.structure.richcmp import richcmp, rich_to_bool
from sage.structure.unique_representation import UniqueRepresentation

from sage.misc.cachefunc import cached_method, cached_function


@cached_function
def list_of_representatives(N):
    """
    Return a list of ideals, coprime to the ideal ``N``, representatives of
    the ideal classes of the corresponding number field.

    .. NOTE::

        This list, used every time we check `\\Gamma_0(N)` - equivalence of
        cusps, is cached.

    INPUT:

    - ``N`` -- an ideal of a number field

    OUTPUT:

    A list of ideals coprime to the ideal ``N``, such that they are
    representatives of all the ideal classes of the number field.

    EXAMPLES::

        sage: from sage.modular.cusps_nf import list_of_representatives
        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^4 + 13*x^3 - 11)
        sage: N = k.ideal(713, a + 208)
        sage: L = list_of_representatives(N); L
        (Fractional ideal (1),
         Fractional ideal (47, a - 9),
         Fractional ideal (53, a - 16))
    """
    return NFCusps_ideal_reps_for_levelN(N)[0]


@cached_function
def NFCusps(number_field):
    r"""
    The set of cusps of a number field `K`, i.e. `\mathbb{P}^1(K)`.

    INPUT:

    - ``number_field`` -- a number field

    OUTPUT: the set of cusps over the given number field

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^2 + 5)
        sage: kCusps = NFCusps(k); kCusps
        Set of all cusps of Number Field in a with defining polynomial x^2 + 5
        sage: kCusps is NFCusps(k)
        True

    Saving and loading works::

        sage: loads(kCusps.dumps()) == kCusps
        True
    """
    return NFCuspsSpace(number_field)


# *************************************************************************
#        NFCuspsSpace class                                               *
# *************************************************************************


class NFCuspsSpace(UniqueRepresentation, Parent):
    """
    The set of cusps of a number field. See ``NFCusps`` for full documentation.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^2 + 5)
        sage: kCusps = NFCusps(k); kCusps
        Set of all cusps of Number Field in a with defining polynomial x^2 + 5
    """
    def __init__(self, number_field):
        """
        See ``NFCusps`` for full documentation.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^3 + x^2 + 13)
            sage: kCusps = NFCusps(k); kCusps
            Set of all cusps of Number Field in a with defining polynomial x^3 + x^2 + 13
        """
        self.__number_field = number_field
        Parent.__init__(self, self)

    def __eq__(self, right):
        """
        Return equality only if right is the set of cusps for the same field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 5)
            sage: L.<a> = NumberField(x^2 + 23)
            sage: kCusps = NFCusps(k); kCusps
            Set of all cusps of Number Field in a with defining polynomial x^2 + 5
            sage: LCusps = NFCusps(L); LCusps
            Set of all cusps of Number Field in a with defining polynomial x^2 + 23
            sage: kCusps == NFCusps(k)
            True
            sage: LCusps == NFCusps(L)
            True
            sage: LCusps == kCusps
            False
        """
        if not isinstance(right, NFCuspsSpace):
            return False
        return self.number_field() == right.number_field()

    def __ne__(self, right):
        """
        Check that ``self`` is not equal to ``right``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 5)
            sage: L.<a> = NumberField(x^2 + 23)
            sage: kCusps = NFCusps(k); kCusps
            Set of all cusps of Number Field in a with defining polynomial x^2 + 5
            sage: LCusps = NFCusps(L); LCusps
            Set of all cusps of Number Field in a with defining polynomial x^2 + 23
            sage: kCusps != NFCusps(k)
            False
            sage: LCusps != NFCusps(L)
            False
            sage: LCusps != kCusps
            True
        """
        return not (self == right)

    def _repr_(self):
        """
        String representation of the set of cusps of a number field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 2)
            sage: kCusps = NFCusps(k)
            sage: kCusps
            Set of all cusps of Number Field in a with defining polynomial x^2 + 2
            sage: kCusps._repr_()
            'Set of all cusps of Number Field in a with defining polynomial x^2 + 2'
            sage: kCusps.rename('Number Field Cusps'); kCusps
            Number Field Cusps
            sage: kCusps.rename(); kCusps
            Set of all cusps of Number Field in a with defining polynomial x^2 + 2
        """
        return "Set of all cusps of %s" % self.number_field()

    def _latex_(self):
        r"""
        Return latex representation of ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 5)
            sage: kCusps = NFCusps(k)
            sage: latex(kCusps) # indirect doctest
            \mathbf{P}^1(\Bold{Q}[a]/(a^{2} + 5))
        """
        return r"\mathbf{P}^1(%s)" % self.number_field()._latex_()

    def __call__(self, x):
        """
        Convert x into the set of cusps of a number field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 5)
            sage: kCusps = NFCusps(k)
            sage: c = kCusps(a,2)
            Traceback (most recent call last):
            ...
            TypeError: ...__call__() takes 2 positional arguments but 3 were given

         ::

            sage: c = kCusps((a,2)); c
            Cusp [a: 2] of Number Field in a with defining polynomial x^2 + 5
            sage: kCusps(2/a)
            Cusp [-2*a: 5] of Number Field in a with defining polynomial x^2 + 5
            sage: kCusps(oo)
            Cusp Infinity of Number Field in a with defining polynomial x^2 + 5
        """
        return NFCusp(self.number_field(), x, parent=self)

    @cached_method
    def zero(self):
        """
        Return the zero cusp.

        .. NOTE::

            This method just exists to make some general algorithms work.
            It is not intended that the returned cusp is an additive
            neutral element.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 5)
            sage: kCusps = NFCusps(k)
            sage: kCusps.zero()
            Cusp [0: 1] of Number Field in a with defining polynomial x^2 + 5
        """
        return self(0)

    def number_field(self):
        """
        Return the number field that this set of cusps is attached to.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 1)
            sage: kCusps = NFCusps(k)
            sage: kCusps.number_field()
            Number Field in a with defining polynomial x^2 + 1
        """
        return self.__number_field

# *************************************************************************
#        NFCusp class                                                     *
# *************************************************************************


class NFCusp(Element):
    r"""
    Create a number field cusp, i.e., an element of `\mathbb{P}^1(k)`.

    A cusp on a number field is either an element of the field or infinity,
    i.e., an element of the projective line over the number field.  It is
    stored as a pair (a,b), where a, b are integral elements of the number
    field.

    INPUT:

    - ``number_field`` -- the number field over which the cusp is defined

    - ``a`` -- it can be a number field element (integral or not), or
      a number field cusp

    - ``b`` -- (optional) when present, it must be either Infinity or
      coercible to an element of the number field

    - ``lreps`` -- (optional) a list of chosen representatives for all the
      ideal classes of the field. When given, the representative of the cusp
      will be changed so its associated ideal is one of the ideals in the list.

    OUTPUT:

    ``[a: b]`` -- a number field cusp.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^2 + 5)
        sage: NFCusp(k, a, 2)
        Cusp [a: 2] of Number Field in a with defining polynomial x^2 + 5
        sage: NFCusp(k, (a,2))
        Cusp [a: 2] of Number Field in a with defining polynomial x^2 + 5
        sage: NFCusp(k, a, 2/(a+1))
        Cusp [a - 5: 2] of Number Field in a with defining polynomial x^2 + 5

    Cusp Infinity:

    ::

        sage: NFCusp(k, 0)
        Cusp [0: 1] of Number Field in a with defining polynomial x^2 + 5
        sage: NFCusp(k, oo)
        Cusp Infinity of Number Field in a with defining polynomial x^2 + 5
        sage: NFCusp(k, 3*a, oo)
        Cusp [0: 1] of Number Field in a with defining polynomial x^2 + 5
        sage: NFCusp(k, a + 5, 0)
        Cusp Infinity of Number Field in a with defining polynomial x^2 + 5

    Saving and loading works:

    ::

        sage: alpha = NFCusp(k, a, 2/(a+1))
        sage: loads(dumps(alpha))==alpha
        True

    Some tests:

    ::

        sage: I*I
        -1
        sage: NFCusp(k, I)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert I to a cusp of the number field

    ::

        sage: NFCusp(k, oo, oo)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert (+Infinity, +Infinity) to a cusp of the number field

    ::

        sage: NFCusp(k, 0, 0)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert (0, 0) to a cusp of the number field

    ::

        sage: NFCusp(k, "a + 2", a)
        Cusp [-2*a + 5: 5] of Number Field in a with defining polynomial x^2 + 5

    ::

        sage: NFCusp(k, NFCusp(k, oo))
        Cusp Infinity of Number Field in a with defining polynomial x^2 + 5
        sage: c = NFCusp(k, 3, 2*a)
        sage: NFCusp(k, c, a + 1)
        Cusp [-a - 5: 20] of Number Field in a with defining polynomial x^2 + 5
        sage: L.<b> = NumberField(x^2 + 2)
        sage: NFCusp(L, c)
        Traceback (most recent call last):
        ...
        ValueError: Cannot coerce cusps from one field to another
    """
    def __init__(self, number_field, a, b=None, parent=None, lreps=None):
        """
        Constructor of number field cusps. See ``NFCusp`` for full
        documentation.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 1)
            sage: c = NFCusp(k, 3, a+1); c
            Cusp [3: a + 1] of Number Field in a with defining polynomial x^2 + 1
            sage: c.parent()
            Set of all cusps of Number Field in a with defining polynomial x^2 + 1
            sage: kCusps = NFCusps(k)
            sage: c.parent() is kCusps
            True
        """
        if parent is None:
            parent = NFCusps(number_field)
        Element.__init__(self, parent)
        R = number_field.maximal_order()
        if b is None:
            if not a:  # that is cusp "0"
                self.__a = R.zero()
                self.__b = R.one()
                return
            if isinstance(a, NFCusp):
                if a.parent() == parent:
                    self.__a = R(a.__a)
                    self.__b = R(a.__b)
                else:
                    raise ValueError("Cannot coerce cusps from one field to another")
            elif a in R:
                self.__a = R(a)
                self.__b = R.one()
            elif a in number_field:
                self.__b = R(a.denominator())
                self.__a = R(a * self.__b)
            elif isinstance(a, InfinityElement):
                self.__a = R.one()
                self.__b = R.zero()
            elif isinstance(a, int):
                self.__a = R(a)
                self.__b = R.one()
            elif isinstance(a, (tuple, list)):
                if len(a) != 2:
                    raise TypeError("unable to convert %r to a cusp \
                                      of the number field" % a)
                if a[1].is_zero():
                    self.__a = R.one()
                    self.__b = R.zero()
                elif a[0] in R and a[1] in R:
                    self.__a = R(a[0])
                    self.__b = R(a[1])
                elif isinstance(a[0], NFCusp):  # we know that a[1] is not zero
                    if a[1] == 1:
                        self.__a = a[0].__a
                        self.__b = a[0].__b
                    else:
                        r = a[0].__a / (a[0].__b * a[1])
                        self.__b = R(r.denominator())
                        self.__a = R(r * self.__b)
                else:
                    try:
                        r = number_field(a[0] / a[1])
                        self.__b = R(r.denominator())
                        self.__a = R(r * self.__b)
                    except (ValueError, TypeError):
                        raise TypeError(f"unable to convert {a} to a cusp "
                                        "of the number field")
            else:
                try:
                    r = number_field(a)
                    self.__b = R(r.denominator())
                    self.__a = R(r * self.__b)
                except (ValueError, TypeError):
                    raise TypeError("unable to convert %r to a cusp "
                                    "of the number field" % a)
        else:  # 'b' is given
            if isinstance(b, InfinityElement):
                if isinstance(a, InfinityElement) or (isinstance(a, NFCusp) and a.is_infinity()):
                    raise TypeError("unable to convert (%r, %r) "
                                    "to a cusp of the number field" % (a, b))
                self.__a = R.zero()
                self.__b = R.one()
                return
            elif not b:
                if not a:
                    raise TypeError("unable to convert (%r, %r) "
                                    "to a cusp of the number field" % (a, b))
                self.__a = R.one()
                self.__b = R.zero()
                return
            if not a:
                self.__a = R.zero()
                self.__b = R.one()
                return
            if (b in R or isinstance(b, int)) and (a in R or isinstance(a, int)):
                self.__a = R(a)
                self.__b = R(b)
            else:
                if a in R or a in number_field:
                    r = a / b
                elif isinstance(a, InfinityElement):
                    self.__a = R.one()
                    self.__b = R.zero()
                    return
                elif isinstance(a, NFCusp):
                    if a.is_infinity():
                        self.__a = R.one()
                        self.__b = R.zero()
                        return
                    r = a.__a / (a.__b * b)
                elif isinstance(a, int):
                    r = R(a) / b
                elif isinstance(a, (tuple, list)):
                    if len(a) != 2:
                        raise TypeError("unable to convert (%r, %r) \
                                          to a cusp of the number field" % (a, b))
                    r = R(a[0]) / (R(a[1]) * b)
                else:
                    try:
                        r = number_field(a) / b
                    except (ValueError, TypeError):
                        raise TypeError("unable to convert (%r, %r) \
                                          to a cusp of the number field" % (a, b))
                self.__b = R(r.denominator())
                self.__a = R(r * self.__b)
        if lreps is not None:
            # Changes the representative of the cusp so the ideal associated
            # to the cusp is one of the ideals of the given list lreps.
            # Note: the trivial class is always represented by (1).
            I = self.ideal()
            for J in lreps:
                if (J / I).is_principal():
                    newI = J
            l = (newI / I).gens_reduced()[0]
            self.__a = R(l * self.__a)
            self.__b = R(l * self.__b)

    def _repr_(self):
        """
        String representation of this cusp.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 1)
            sage: c = NFCusp(k, a, 2); c
            Cusp [a: 2] of Number Field in a with defining polynomial x^2 + 1
            sage: c._repr_()
            'Cusp [a: 2] of Number Field in a with defining polynomial x^2 + 1'
            sage: c.rename('[a:2](cusp of a number field)');c
            [a:2](cusp of a number field)
            sage: c.rename();c
            Cusp [a: 2] of Number Field in a with defining polynomial x^2 + 1
        """
        if self.__b.is_zero():
            return "Cusp Infinity of %s" % self.parent().number_field()
        else:
            return "Cusp [%s: %s] of %s" % (self.__a, self.__b,
                                            self.parent().number_field())

    def number_field(self):
        """
        Return the number field of definition of the cusp ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 2)
            sage: alpha = NFCusp(k, 1, a + 1)
            sage: alpha.number_field()
            Number Field in a with defining polynomial x^2 + 2
        """
        return self.parent().number_field()

    def is_infinity(self):
        """
        Return ``True`` if this is the cusp infinity.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 1)
            sage: NFCusp(k, a, 2).is_infinity()
            False
            sage: NFCusp(k, 2, 0).is_infinity()
            True
            sage: NFCusp(k, oo).is_infinity()
            True
        """
        return self.__b == 0

    def numerator(self):
        """
        Return the numerator of the cusp ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 1)
            sage: c = NFCusp(k, a, 2)
            sage: c.numerator()
            a
            sage: d = NFCusp(k, 1, a)
            sage: d.numerator()
            1
            sage: NFCusp(k, oo).numerator()
            1
        """
        return self.__a

    def denominator(self):
        """
        Return the denominator of the cusp ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 1)
            sage: c = NFCusp(k, a, 2)
            sage: c.denominator()
            2
            sage: d = NFCusp(k, 1, a + 1);d
            Cusp [1: a + 1] of Number Field in a with defining polynomial x^2 + 1
            sage: d.denominator()
            a + 1
            sage: NFCusp(k, oo).denominator()
            0
        """
        return self.__b

    def _number_field_element_(self):
        """
        Coerce to an element of the number field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 2)
            sage: NFCusp(k, a, 2)._number_field_element_()
            1/2*a
            sage: NFCusp(k, 1, a + 1)._number_field_element_()
            -1/3*a + 1/3
        """
        if self.__b.is_zero():
            raise TypeError("%s is not an element of %s" % (self,
                                                            self.number_field()))
        k = self.number_field()
        return k(self.__a / self.__b)

    def _ring_of_integers_element_(self):
        """
        Coerce to an element of the ring of integers of the number field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 2)
            sage: NFCusp(k, a+1)._ring_of_integers_element_()
            a + 1
            sage: NFCusp(k, 1, a + 1)._ring_of_integers_element_()
            Traceback (most recent call last):
            ...
            TypeError: Cusp [1: a + 1] of Number Field in a with defining polynomial x^2 + 2 is not an integral element
        """
        if self.__b.is_one():
            return self.__a
        R = self.number_field().ring_of_integers()
        if self.__b.is_zero():
            raise TypeError("%s is not an element of %s" % (self, R))
        try:
            return R(self.__a / self.__b)
        except (ValueError, TypeError):
            raise TypeError("%s is not an integral element" % self)

    def _latex_(self):
        r"""
        Latex representation of this cusp.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 11)
            sage: latex(NFCusp(k, 3*a, a + 1)) # indirect doctest
            \[3 a: a + 1\]
            sage: latex(NFCusp(k, 3*a, a + 1)) == NFCusp(k, 3*a, a + 1)._latex_()
            True
            sage: latex(NFCusp(k, oo))
            \infty
        """
        if self.__b.is_zero():
            return "\\infty"
        else:
            return "\\[%s: %s\\]" % (self.__a._latex_(),
                                     self.__b._latex_())

    def _richcmp_(self, right, op):
        """
        Compare the cusps ``self`` and ``right``.

        Comparison is as for elements in the number field, except with
        the cusp oo which is greater than everything but itself.

        The ordering in comparison is only really meaningful for infinity.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^3 + x + 1)
            sage: kCusps = NFCusps(k)

        Comparing with infinity::

            sage: c = kCusps((a,2))
            sage: d = kCusps(oo)
            sage: c < d
            True
            sage: kCusps(oo) < d
            False

        Comparison as elements of the number field::

            sage: kCusps(2/3) < kCusps(5/2)
            False
            sage: k(2/3) < k(5/2)
            False
        """
        if self.__b.is_zero():
            if right.__b.is_zero():
                return rich_to_bool(op, 0)
            else:
                return rich_to_bool(op, 1)
        else:
            if right.__b.is_zero():
                return rich_to_bool(op, -1)
            else:
                return richcmp(self._number_field_element_(),
                               right._number_field_element_(), op)

    def __neg__(self):
        """
        The negative of this cusp.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: c = NFCusp(k, a, a+1); c
            Cusp [a: a + 1] of Number Field in a with defining polynomial x^2 + 23
            sage: -c
            Cusp [-a: a + 1] of Number Field in a with defining polynomial x^2 + 23
        """
        return NFCusp(self.parent().number_field(), -self.__a, self.__b)

    def apply(self, g):
        """
        Return g(``self``), where ``g`` is a 2x2 matrix, which we view as a
        linear fractional transformation.

        INPUT:

        - ``g`` -- list of integral elements [a, b, c, d] that are the
          entries of a 2x2 matrix

        OUTPUT:

        A number field cusp, obtained by the action of ``g`` on the cusp
        ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: beta = NFCusp(k, 0, 1)
            sage: beta.apply([0, -1, 1, 0])
            Cusp Infinity of Number Field in a with defining polynomial x^2 + 23
            sage: beta.apply([1, a, 0, 1])
            Cusp [a: 1] of Number Field in a with defining polynomial x^2 + 23
        """
        k = self.number_field()
        return NFCusp(k, g[0] * self.__a + g[1] * self.__b,
                      g[2] * self.__a + g[3] * self.__b)

    def ideal(self):
        """
        Return the ideal associated to the cusp ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: alpha = NFCusp(k, 3, a-1)
            sage: alpha.ideal()
            Fractional ideal (3, 1/2*a - 1/2)
            sage: NFCusp(k, oo).ideal()
            Fractional ideal (1)
        """
        k = self.number_field()
        return k.ideal(self.__a, self.__b)

    def ABmatrix(self):
        """
        Return AB-matrix associated to the cusp ``self``.

        Given R a Dedekind domain and A, B ideals of R in inverse classes, an
        AB-matrix is a matrix realizing the isomorphism between R+R and A+B.
        An AB-matrix associated to a cusp [a1: a2] is an AB-matrix with A the
        ideal associated to the cusp (A=<a1, a2>) and first column given by
        the coefficients of the cusp.

        EXAMPLES:

        ::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: alpha = NFCusp(k, oo)
            sage: alpha.ABmatrix()
            [1, 0, 0, 1]

        ::

            sage: alpha = NFCusp(k, 0)
            sage: alpha.ABmatrix()
            [0, -1, 1, 0]

        Note that the AB-matrix associated to a cusp is not unique, and the
        output of the ``ABmatrix`` function may change.

        ::

            sage: alpha = NFCusp(k, 3/2, a-1)
            sage: M = alpha.ABmatrix()
            sage: M # random
            [-a^2 - a - 1, -3*a - 7, 8, -2*a^2 - 3*a + 4]
            sage: M[0] == alpha.numerator() and M[2] == alpha.denominator()
            True

        An AB-matrix associated to a cusp alpha will send Infinity to alpha:

        ::

            sage: alpha = NFCusp(k, 3, a-1)
            sage: M = alpha.ABmatrix()
            sage: (k.ideal(M[1], M[3])*alpha.ideal()).is_principal()
            True
            sage: M[0] == alpha.numerator() and M[2] == alpha.denominator()
            True
            sage: NFCusp(k, oo).apply(M) == alpha
            True
        """
        k = self.number_field()
        A = self.ideal()

        if self.is_infinity():
            return [1, 0, 0, 1]
        if not self:
            return [0, -1, 1, 0]

        if A.is_principal():
            B = k.ideal(1)
        else:
            B = k.ideal(A.gens_reduced()[1]) / A
        assert (A * B).is_principal()

        a1 = self.__a
        a2 = self.__b

        g = (A * B).gens_reduced()[0]
        Ainv = A**(-1)
        A1 = a1 * Ainv
        A2 = a2 * Ainv
        r = A1.element_1_mod(A2)
        b1 = -(1 - r) / a2 * g
        b2 = (r / a1) * g
        ABM = [a1, b1, a2, b2]

        return ABM

    def is_Gamma0_equivalent(self, other, N, Transformation=False):
        r"""
        Check if cusps ``self`` and ``other`` are `\Gamma_0(N)`- equivalent.

        INPUT:

        - ``other`` -- a number field cusp or a list of two number field
          elements which define a cusp

        - ``N`` -- an ideal of the number field (level)

        OUTPUT: boolean; ``True`` if the cusps are equivalent

        - a transformation matrix -- (if ``Transformation=True``) a list of
          integral elements [a, b, c, d] which are the entries of a 2x2 matrix
          M in `\Gamma_0(N)` such that M * ``self`` = ``other`` if ``other``
          and ``self`` are `\Gamma_0(N)`- equivalent. If ``self`` and ``other``
          are not equivalent it returns zero.

        EXAMPLES:

        ::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 10)
            sage: N = K.ideal(a - 1)
            sage: alpha = NFCusp(K, 0)
            sage: beta = NFCusp(K, oo)
            sage: alpha.is_Gamma0_equivalent(beta, N)
            False
            sage: alpha.is_Gamma0_equivalent(beta, K.ideal(1))
            True
            sage: b, M = alpha.is_Gamma0_equivalent(beta, K.ideal(1),Transformation=True)
            sage: alpha.apply(M)
            Cusp Infinity of Number Field in a with defining polynomial x^3 - 10

        ::

            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3)
            sage: alpha1 = NFCusp(k, a+1, 4)
            sage: alpha2 = NFCusp(k, a-8, 29)
            sage: alpha1.is_Gamma0_equivalent(alpha2, N)
            True
            sage: b, M = alpha1.is_Gamma0_equivalent(alpha2, N, Transformation=True)
            sage: alpha1.apply(M) == alpha2
            True
            sage: M[2] in N
            True
        """
        k = self.number_field()
        other = NFCusp(k, other)
        if not (self.ideal() / other.ideal()).is_principal():
            if not Transformation:
                return False
            else:
                return False, 0

        reps = list_of_representatives(N)
        alpha1 = NFCusp(k, self, lreps=reps)
        alpha2 = NFCusp(k, other, lreps=reps)

        delta = k.ideal(alpha1.__b) + N
        if (k.ideal(alpha2.__b) + N) != delta:
            if not Transformation:
                return False
            else:
                return False, 0

        M1 = alpha1.ABmatrix()
        M2 = alpha2.ABmatrix()

        A = alpha1.ideal()
        B = k.ideal(M1[1], M1[3])

        ABdelta = A * B * delta * delta

        units = units_mod_ideal(ABdelta)
        for u in units:
            if (M2[2] * M1[3] - u * M1[2] * M2[3]) in ABdelta:
                if not Transformation:
                    return True
                else:
                    AuxCoeff = [1, 0, 0, 1]
                    Aux = M2[2] * M1[3] - u * M1[2] * M2[3]
                    if Aux in A * B * N:
                        if u != 1:
                            AuxCoeff[3] = u
                    else:
                        A1 = (A * B * N) / ABdelta
                        A2 = B * k.ideal(M1[2] * M2[2]) / (A * ABdelta)
                        f = A1.element_1_mod(A2)
                        w = ((1 - f) * Aux) / (M1[2] * M2[2])
                        AuxCoeff[3] = u
                        AuxCoeff[1] = w
                    from sage.matrix.constructor import Matrix
                    Maux = Matrix(k, 2, AuxCoeff)
                    M1inv = Matrix(k, 2, M1).inverse()
                    Mtrans = Matrix(k, 2, M2) * Maux * M1inv
                    assert Mtrans[1][0] in N
                    return True, Mtrans.list()
        if not Transformation:
            return False
        else:
            return False, 0

# *************************************************************************
#  Global functions:
#    - Gamma0_NFCusps --compute list of inequivalent cusps
#    Internal use only:
#    - number_of_Gamma0_NFCusps -- useful to test Gamma0_NFCusps
#    - NFCusps_ideal_reps_for_levelN -- lists of reps for ideal classes
#    - units_mod_ideal -- needed to check Gamma0(N)-equiv of cusps
# *************************************************************************


def Gamma0_NFCusps(N):
    r"""
    Return a list of inequivalent cusps for `\Gamma_0(N)`, i.e., a set of
    representatives for the orbits of ``self`` on `\mathbb{P}^1(k)`.

    INPUT:

    - ``N`` -- an integral ideal of the number field k (the level)

    OUTPUT: list of inequivalent number field cusps

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^2 + 5)
        sage: N = k.ideal(3)
        sage: L = Gamma0_NFCusps(N)

    The cusps in the list are inequivalent::

        sage: any(L[i].is_Gamma0_equivalent(L[j], N)
        ....:     for i in range(len(L)) for j in range(len(L)) if i < j)
        False

    We test that we obtain the right number of orbits::

        sage: from sage.modular.cusps_nf import number_of_Gamma0_NFCusps
        sage: len(L) == number_of_Gamma0_NFCusps(N)
        True

    Another example::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^4 - x^3 -21*x^2 + 17*x + 133)
        sage: N = k.ideal(5)
        sage: from sage.modular.cusps_nf import number_of_Gamma0_NFCusps
        sage: len(Gamma0_NFCusps(N)) == number_of_Gamma0_NFCusps(N) # long time (over 1 sec)
        True
    """
    # We create L a list of three lists, which are different and each a list of
    # prime ideals, coprime to N, representing the ideal classes of k
    L = NFCusps_ideal_reps_for_levelN(N, nlists=3)
    Laux = L[1] + L[2]
    Lreps = list_of_representatives(N)
    Lcusps = []

    k = N.number_field()

    for A in L[0]:
        # find B in inverse class:
        if A.is_trivial():
            B = k.ideal(1)
            # B = k.unit_ideal() produces an error because we need fract ideal
            g = 1
        else:
            Lbs = [P for P in Laux if (P * A).is_principal()]
            B = Lbs[0]
            g = (A * B).gens_reduced()[0]

        # for every divisor of N we have to find cusps
        from sage.arith.misc import divisors
        for d in divisors(N):
            # find delta prime coprime to B in inverse class of d*A
            # by searching in our list of auxiliary prime ideals
            Lds = [P for P in Laux
                   if (P * d * A).is_principal() and P.is_coprime(B)]
            deltap = Lds[0]
            a = (deltap * d * A).gens_reduced()[0]
            I = d + N / d
            # special case: A=B=d=<1>:
            if a.is_one() and I.is_trivial():
                Lcusps.append(NFCusp(k, 0, 1, lreps=Lreps))
            else:
                u = k.unit_group().gens()
                for b in I.invertible_residues_mod(u):
                    # Note: if I trivial, invertible_residues_mod returns [1]
                    # lift b to (R/a)star
                    # we need the part of d which is coprime to I, call it M
                    M = d.prime_to_idealM_part(I)
                    deltAM = deltap * A * M
                    u = (B * deltAM).element_1_mod(I)
                    v = (I * B).element_1_mod(deltAM)
                    newb = u * b + v
                    # build AB-matrix:
                    # ----> extended gcd for k.ideal(a), k.ideal(newb)
                    Y = k.ideal(newb).element_1_mod(k.ideal(a))
                    # if xa + yb = 1, cusp = y*g /a
                    Lcusps.append(NFCusp(k, Y * g, a, lreps=Lreps))
    return Lcusps


def number_of_Gamma0_NFCusps(N):
    """
    Return the total number of orbits of cusps under the action of the
    congruence subgroup `\\Gamma_0(N)`.

    INPUT:

    - ``N`` -- a number field ideal

    OUTPUT: integer; the number of orbits of cusps under Gamma0(N)-action

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^3 + 11)
        sage: N = k.ideal(2, a+1)
        sage: from sage.modular.cusps_nf import number_of_Gamma0_NFCusps
        sage: number_of_Gamma0_NFCusps(N)
        4
        sage: L = Gamma0_NFCusps(N)
        sage: len(L) == number_of_Gamma0_NFCusps(N)
        True
        sage: k.<a> = NumberField(x^2 + 7)
        sage: N = k.ideal(9)
        sage: number_of_Gamma0_NFCusps(N)
        6
        sage: N = k.ideal(a*9 + 7)
        sage: number_of_Gamma0_NFCusps(N)
        24
    """
    k = N.number_field()
    # The number of Gamma0(N)-sub-orbits for each Gamma-orbit:
    from sage.arith.misc import divisors
    Ugens = [k(u) for u in k.unit_group().gens()]
    s = sum([len((d + N / d).invertible_residues_mod(Ugens))
             for d in divisors(N)])
    # There are h Gamma-orbits, with h class number of underlying number field.
    return s * k.class_number()


def NFCusps_ideal_reps_for_levelN(N, nlists=1):
    """
    Return a list of lists (``nlists`` different lists) of prime ideals,
    coprime to ``N``, representing every ideal class of the number field.

    INPUT:

    - ``N`` -- number field ideal

    - ``nlists`` -- (default: 1) the number of lists of prime ideals
      we want

    OUTPUT:

    A list of lists of ideals representatives of the ideal classes, all coprime
    to ``N``, representing every ideal.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^3 + 11)
        sage: N = k.ideal(5, a + 1)
        sage: from sage.modular.cusps_nf import NFCusps_ideal_reps_for_levelN
        sage: NFCusps_ideal_reps_for_levelN(N)
        [(Fractional ideal (1), Fractional ideal (2, a + 1))]
        sage: L = NFCusps_ideal_reps_for_levelN(N, 3)
        sage: all(len(L[i]) == k.class_number() for i in range(len(L)))
        True

    ::

        sage: k.<a> = NumberField(x^4 - x^3 - 21*x^2 + 17*x + 133)
        sage: N = k.ideal(6)
        sage: from sage.modular.cusps_nf import NFCusps_ideal_reps_for_levelN
        sage: NFCusps_ideal_reps_for_levelN(N)
        [(Fractional ideal (1),
          Fractional ideal (67, -4/7*a^3 + 13/7*a^2 + 39/7*a - 43),
          Fractional ideal (127, -4/7*a^3 + 13/7*a^2 + 39/7*a - 42),
          Fractional ideal (157, -4/7*a^3 + 13/7*a^2 + 39/7*a + 48))]
        sage: L = NFCusps_ideal_reps_for_levelN(N, 5)
        sage: all(len(L[i]) == k.class_number() for i in range(len(L)))
        True
    """
    k = N.number_field()
    G = k.class_group()
    L = [[k.ideal(1)] for _ in range(nlists)]
    it = k.primes_of_degree_one_iter()
    for I in G.list():
        check = 0
        if not I.is_principal():
            Iinv = I.ideal()**(-1)
            while check < nlists:
                J = next(it)
                if (J * Iinv).is_principal() and J.is_coprime(N):
                    L[check].append(J)
                    check += 1
    return [tuple(l) for l in L]


def units_mod_ideal(I):
    """
    Return integral elements of the number field representing the images of
    the global units modulo the ideal ``I``.

    INPUT:

    - ``I`` -- number field ideal

    OUTPUT:

    A list of integral elements of the number field representing the images of
    the global units modulo the ideal ``I``. Elements of the list might be
    equivalent to each other mod ``I``.

    EXAMPLES::

        sage: from sage.modular.cusps_nf import units_mod_ideal
        sage: x = polygen(ZZ, 'x')
        sage: k.<a> = NumberField(x^2 + 1)
        sage: I = k.ideal(a + 1)
        sage: units_mod_ideal(I)
        [1]
        sage: I = k.ideal(3)
        sage: units_mod_ideal(I)
        [1, a, -1, -a]

    ::

        sage: from sage.modular.cusps_nf import units_mod_ideal
        sage: k.<a> = NumberField(x^3 + 11)
        sage: k.unit_group()
        Unit group with structure C2 x Z of
         Number Field in a with defining polynomial x^3 + 11
        sage: I = k.ideal(5, a + 1)
        sage: units_mod_ideal(I)
        [1,
        2*a^2 + 4*a - 1,
        ...]

    ::

        sage: from sage.modular.cusps_nf import units_mod_ideal
        sage: k.<a> = NumberField(x^4 - x^3 -21*x^2 + 17*x + 133)
        sage: k.unit_group()
        Unit group with structure C6 x Z of
         Number Field in a with defining polynomial x^4 - x^3 - 21*x^2 + 17*x + 133
        sage: I = k.ideal(3)
        sage: U = units_mod_ideal(I)
        sage: all(U[j].is_unit() and (U[j] not in I) for j in range(len(U)))
        True
    """
    k = I.number_field()
    Uk = k.unit_group()
    Istar = I.idealstar(2)
    ulist = Uk.gens_values()
    elist = [Istar(I.ideallog(u)).order() for u in ulist]

    from sage.misc.mrange import xmrange

    return [k.prod(u**e for u, e in zip(ulist, ei)) for ei in xmrange(elist)]
