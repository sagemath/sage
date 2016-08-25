"""
Monoid of ideals in a commutative ring
"""
from __future__ import absolute_import

from sage.structure.parent import Parent
import sage.rings.integer_ring
from . import ideal
from sage.categories.monoids import Monoids

def IdealMonoid(R):
    r"""
    Return the monoid of ideals in the ring ``R``.

    EXAMPLE::

        sage: R = QQ['x']
        sage: sage.rings.ideal_monoid.IdealMonoid(R)
        Monoid of ideals of Univariate Polynomial Ring in x over Rational Field
    """
    return IdealMonoid_c(R)

class IdealMonoid_c(Parent):
    r"""
    The monoid of ideals in a commutative ring.

    TESTS::

        sage: R = QQ['x']
        sage: M = sage.rings.ideal_monoid.IdealMonoid(R)
        sage: TestSuite(M).run()
          Failure in _test_category:
        ...
        The following tests failed: _test_elements

    (The "_test_category" test fails but I haven't the foggiest idea why.)
    """

    Element = ideal.Ideal_generic # this doesn't seem to do anything

    def __init__(self, R):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R = QuadraticField(-23, 'a')
            sage: M = sage.rings.ideal_monoid.IdealMonoid(R); M # indirect doctest
            Monoid of ideals of Number Field in a with defining polynomial x^2 + 23
        """
        self.__R = R
        Parent.__init__(self, base = sage.rings.integer_ring.ZZ, category = Monoids())
        self._populate_coercion_lists_()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: R = QuadraticField(-23, 'a')
            sage: M = sage.rings.ideal_monoid.IdealMonoid(R); M._repr_()
            'Monoid of ideals of Number Field in a with defining polynomial x^2 + 23'
        """
        return "Monoid of ideals of %s"%self.__R

    def ring(self):
        r"""
        Return the ring of which this is the ideal monoid.

        EXAMPLE::

            sage: R = QuadraticField(-23, 'a')
            sage: M = sage.rings.ideal_monoid.IdealMonoid(R); M.ring() is R
            True
        """
        return self.__R

    def _element_constructor_(self, x):
        r"""
        Create an ideal in this monoid from ``x``.

        EXAMPLES::

            sage: R.<a> = QuadraticField(-23)
            sage: M = sage.rings.ideal_monoid.IdealMonoid(R)
            sage: M(a)   # indirect doctest
            Fractional ideal (a)
            sage: M([a-4, 13])
            Fractional ideal (13, 1/2*a + 9/2)
        """
        try:
            side = x.side()
        except (AttributeError,TypeError):
            side = None
        try:
            x = x.gens()
        except AttributeError:
            pass
        if side is None:
            y = self.__R.ideal(x)
        else:
            y = self.__R.ideal(x,side=side)
        y._set_parent(self)
        return y

    def _coerce_map_from_(self, x):
        r"""
        Used by coercion framework.

        EXAMPLE::

            sage: R = QuadraticField(-23, 'a')
            sage: M = R.ideal_monoid()
            sage: M.has_coerce_map_from(R) # indirect doctest
            True
            sage: M.has_coerce_map_from(QQ.ideal_monoid())
            True
            sage: M.has_coerce_map_from(Zmod(6))
            False
            sage: M.has_coerce_map_from(loads(dumps(M)))
            True
        """
        if isinstance(x, IdealMonoid_c):
            return self.ring().has_coerce_map_from(x.ring())
        else:
            return self.ring().has_coerce_map_from(x)

    def __cmp__(self, other):
        r"""
        Comparison function.

        EXAMPLE::

            sage: R = QuadraticField(-23, 'a')
            sage: M = R.ideal_monoid()
            sage: M == QQ
            False
            sage: M == 17
            False
            sage: M == R.ideal_monoid()
            True
        """
        if not isinstance(other, IdealMonoid_c):
            return cmp(type(self), type(other))
        else:
            return cmp(self.ring(), other.ring())
