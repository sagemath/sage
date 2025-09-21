"""
Divisor groups

AUTHORS:

- David Kohel (2006): Initial version

- Volker Braun (2010-07-16): Documentation, doctests, coercion fixes, bugfixes.
"""

#*******************************************************************************
#  Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.schemes.generic.divisor import Divisor_generic, Divisor_curve
from sage.structure.formal_sum import FormalSums
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

lazy_import('sage.schemes.curves.curve', 'Curve_generic')


def DivisorGroup(scheme, base_ring=None):
    r"""
    Return the group of divisors on the scheme.

    INPUT:

    - ``scheme`` -- a scheme

    - ``base_ring`` -- usually either `\ZZ` (default) or `\QQ`. The
      coefficient ring of the divisors. Not to be confused with the
      base ring of the scheme!

    OUTPUT: an instance of ``DivisorGroup_generic``

    EXAMPLES::

        sage: from sage.schemes.generic.divisor_group import DivisorGroup
        sage: DivisorGroup(Spec(ZZ))
        Group of ZZ-Divisors on Spectrum of Integer Ring
        sage: DivisorGroup(Spec(ZZ), base_ring=QQ)
        Group of QQ-Divisors on Spectrum of Integer Ring
    """
    if base_ring is None:
        base_ring = ZZ

    if isinstance(scheme, Curve_generic):
        DG = DivisorGroup_curve(scheme, base_ring)
    else:
        DG = DivisorGroup_generic(scheme, base_ring)

    return DG


def is_DivisorGroup(x):
    r"""
    Return whether ``x`` is a :class:`DivisorGroup_generic`.

    INPUT:

    - ``x`` -- anything

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.schemes.generic.divisor_group import is_DivisorGroup, DivisorGroup
        sage: Div = DivisorGroup(Spec(ZZ), base_ring=QQ)
        sage: is_DivisorGroup(Div)
        doctest:warning...
        DeprecationWarning: The function is_DivisorGroup is deprecated; use 'isinstance(..., DivisorGroup_generic)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        True
        sage: is_DivisorGroup('not a divisor')
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_DivisorGroup is deprecated; use 'isinstance(..., DivisorGroup_generic)' instead.")
    return isinstance(x, DivisorGroup_generic)


class DivisorGroup_generic(FormalSums):
    r"""
    The divisor group on a variety.
    """

    @staticmethod
    def __classcall__(cls, scheme, base_ring=ZZ):
        """
        Set the default value for the base ring.

        EXAMPLES::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup_generic
            sage: DivisorGroup_generic(Spec(ZZ),ZZ) == DivisorGroup_generic(Spec(ZZ))    # indirect test
            True
        """
        # Must not call super().__classcall__()!
        return UniqueRepresentation.__classcall__(cls, scheme, base_ring)

    def __init__(self, scheme, base_ring):
        r"""
        Construct a :class:`DivisorGroup_generic`.

        INPUT:

        - ``scheme`` -- a scheme

        - ``base_ring`` -- the coefficient ring of the divisor group

        Implementation note: :meth:`__classcall__` sets default value
        for ``base_ring``.

        EXAMPLES::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup_generic
            sage: DivisorGroup_generic(Spec(ZZ), QQ)
            Group of QQ-Divisors on Spectrum of Integer Ring

        TESTS::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: D1 = DivisorGroup(Spec(ZZ))
            sage: D2 = DivisorGroup(Spec(ZZ), base_ring=QQ)
            sage: D3 = DivisorGroup(Spec(QQ))
            sage: D1 == D1
            True
            sage: D1 == D2
            False
            sage: D1 != D3
            True
            sage: D2 == D2
            True
            sage: D2 == D3
            False
            sage: D3 != D3
            False
            sage: D1 == 'something'
            False
        """
        super().__init__(base_ring)
        self._scheme = scheme

    def _repr_(self):
        r"""
        Return a string representation of the divisor group.

        EXAMPLES::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: DivisorGroup(Spec(ZZ), base_ring=QQ)
            Group of QQ-Divisors on Spectrum of Integer Ring
        """
        ring = self.base_ring()
        if ring == ZZ:
            base_ring_str = 'ZZ'
        elif ring == QQ:
            base_ring_str = 'QQ'
        else:
            base_ring_str = '('+str(ring)+')'
        return 'Group of '+base_ring_str+'-Divisors on '+str(self._scheme)

    def _element_constructor_(self, x, check=True, reduce=True):
        r"""
        Construct an element of the divisor group.

        EXAMPLES::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: DivZZ = DivisorGroup(Spec(ZZ))
            sage: DivZZ([(2,5)])
            2*V(5)
        """
        if isinstance(x, Divisor_generic):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return Divisor_generic(x._data, check=False, reduce=False, parent=self)
            else:
                x = x._data
        if isinstance(x, list):
            return Divisor_generic(x, check=check, reduce=reduce, parent=self)
        if x == 0:
            return Divisor_generic([], check=False, reduce=False, parent=self)
        else:
            return Divisor_generic([(self.base_ring()(1), x)], check=False, reduce=False, parent=self)

    def _coerce_map_from_(self, other):
        r"""
        Return if there is a coercion map from ``other`` to ``self``.

        There is a coercion from another divisor group if there is
        a coercion map from the schemes and there is a coercion map from
        the base rings.

        TESTS::

            sage: C = EllipticCurve([2, 1])
            sage: E = EllipticCurve([1, 2])
            sage: C.divisor_group()._coerce_map_from_(E.divisor_group())
            False
            sage: E.divisor_group()._coerce_map_from_(C.divisor_group())
            False
            sage: E.divisor_group()._coerce_map_from_(E.divisor_group())
            True
            sage: C.divisor_group()._coerce_map_from_(C.divisor_group())
            True
            sage: D = 1/2 * E.divisor(E(1, 2))
            sage: D.parent()._coerce_map_from_(E.divisor_group())
            True
            sage: E.divisor_group()._coerce_map_from_(D.parent())
            False
        """
        return (isinstance(other, DivisorGroup_generic)
                and self.scheme().has_coerce_map_from(other.scheme())
                and super()._coerce_map_from_(other))

    def scheme(self):
        r"""
        Return the scheme supporting the divisors.

        EXAMPLES::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: Div = DivisorGroup(Spec(ZZ))   # indirect test
            sage: Div.scheme()
            Spectrum of Integer Ring
        """
        return self._scheme

    def _an_element_(self):
        r"""
        Return an element of the divisor group.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, CC)                                         # needs sage.rings.real_mpfr
            sage: C = Curve(y^2 - x^9 - x)                                              # needs sage.rings.real_mpfr sage.schemes
            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: DivisorGroup(C).an_element()    # indirect test                       # needs sage.rings.real_mpfr sage.schemes
            0
        """
        return self._scheme.divisor([], base_ring=self.base_ring(), check=False, reduce=False)

    def base_extend(self, R):
        """
        EXAMPLES::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: DivisorGroup(Spec(ZZ), ZZ).base_extend(QQ)
            Group of QQ-Divisors on Spectrum of Integer Ring
            sage: DivisorGroup(Spec(ZZ), ZZ).base_extend(GF(7))
            Group of (Finite Field of size 7)-Divisors on Spectrum of Integer Ring

        Divisor groups are unique::

            sage: A.<x, y> = AffineSpace(2, CC)                                         # needs sage.rings.real_mpfr
            sage: C = Curve(y^2 - x^9 - x)                                              # needs sage.rings.real_mpfr sage.schemes
            sage: DivisorGroup(C, ZZ).base_extend(QQ) is DivisorGroup(C, QQ)            # needs sage.rings.real_mpfr sage.schemes
            True
        """
        if self.base_ring().has_coerce_map_from(R):
            return self
        elif R.has_coerce_map_from(self.base_ring()):
            return DivisorGroup(self.scheme(), base_ring=R)


class DivisorGroup_curve(DivisorGroup_generic):
    r"""
    Special case of the group of divisors on a curve.
    """

    def _element_constructor_(self, x, check=True, reduce=True):
        r"""
        Construct an element of the divisor group.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr sage.schemes
            sage: A.<x, y> = AffineSpace(2, CC)
            sage: C = Curve(y^2 - x^9 - x)
            sage: DivZZ = C.divisor_group(ZZ)
            sage: DivQQ = C.divisor_group(QQ)
            sage: DivQQ(DivQQ.an_element())   # indirect test
            0
            sage: DivZZ(DivZZ.an_element())   # indirect test
            0
            sage: DivQQ(DivZZ.an_element())   # indirect test
            0
        """
        if isinstance(x, Divisor_curve):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return Divisor_curve(x._data, check=False, reduce=False, parent=self)
            else:
                x = x._data
        if isinstance(x, list):
            return Divisor_curve(x, check=check, reduce=reduce, parent=self)
        if x == 0:
            return Divisor_curve([], check=False, reduce=False, parent=self)
        else:
            return Divisor_curve([(self.base_ring()(1), x)], check=False, reduce=False, parent=self)
