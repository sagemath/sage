"""
Monoid of ideals in a commutative ring

WARNING: This is used by some rings that are not commutative! ::

    sage: MS = MatrixSpace(QQ, 3, 3)                                                    # optional - sage.modules
    sage: type(MS.ideal(MS.one()).parent())                                             # optional - sage.modules
    <class 'sage.rings.ideal_monoid.IdealMonoid_c_with_category'>
"""

from sage.structure.parent import Parent
import sage.rings.integer_ring
from . import ideal
from sage.categories.monoids import Monoids


def IdealMonoid(R):
    r"""
    Return the monoid of ideals in the ring ``R``.

    EXAMPLES::

        sage: R = QQ['x']
        sage: from sage.rings.ideal_monoid import IdealMonoid
        sage: IdealMonoid(R)
        Monoid of ideals of Univariate Polynomial Ring in x over Rational Field
    """
    return IdealMonoid_c(R)


class IdealMonoid_c(Parent):
    r"""
    The monoid of ideals in a commutative ring.

    TESTS::

        sage: R = QQ['x']
        sage: from sage.rings.ideal_monoid import IdealMonoid
        sage: M = IdealMonoid(R)
        sage: TestSuite(M).run()
          Failure in _test_category:
        ...
        The following tests failed: _test_elements

    (The "_test_category" test fails but I haven't the foggiest idea why.)
    """

    Element = ideal.Ideal_generic  # this doesn't seem to do anything

    def __init__(self, R):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R = QuadraticField(-23, 'a')                                          # optional - sage.rings.number_field
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R); M # indirect doctest                              # optional - sage.rings.number_field
            Monoid of ideals of Number Field in a with defining polynomial x^2 + 23
             with a = 4.795831523312720?*I

            sage: id = QQ.ideal(6)
            sage: id.parent().category()
            Category of commutative monoids

            sage: MS = MatrixSpace(QQ, 3, 3)                                            # optional - sage.modules
            sage: MS.ideal(MS.one()).parent().category()                                # optional - sage.modules
            Category of monoids
        """
        self.__R = R
        cat = Monoids()
        if R.is_commutative():
            cat = cat.Commutative()
        Parent.__init__(self, base=sage.rings.integer_ring.ZZ,
                        category=cat)
        self._populate_coercion_lists_()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: R = QuadraticField(-23, 'a')                                          # optional - sage.rings.number_field
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R); M._repr_()                                        # optional - sage.rings.number_field
            'Monoid of ideals of Number Field in a with defining polynomial x^2 + 23
             with a = 4.795831523312720?*I'
        """
        return "Monoid of ideals of %s" % self.__R

    def ring(self):
        r"""
        Return the ring of which this is the ideal monoid.

        EXAMPLES::

            sage: R = QuadraticField(-23, 'a')                                          # optional - sage.rings.number_field
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R); M.ring() is R                                     # optional - sage.rings.number_field
            True
        """
        return self.__R

    def _element_constructor_(self, x):
        r"""
        Create an ideal in this monoid from ``x``.

        EXAMPLES::

            sage: R.<a> = QuadraticField(-23)                                           # optional - sage.rings.number_field
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R)                                                    # optional - sage.rings.number_field
            sage: M(a)   # indirect doctest                                             # optional - sage.rings.number_field
            Fractional ideal (a)
            sage: M([a-4, 13])                                                          # optional - sage.rings.number_field
            Fractional ideal (13, 1/2*a + 9/2)
        """
        try:
            side = x.side()
        except (AttributeError, TypeError):
            side = None
        try:
            x = x.gens()
        except AttributeError:
            pass
        if side is None:
            y = self.__R.ideal(x)
        else:
            y = self.__R.ideal(x, side=side)
        y._set_parent(self)
        return y

    def _coerce_map_from_(self, x):
        r"""
        Used by coercion framework.

        EXAMPLES::

            sage: R = QuadraticField(-23, 'a')                                          # optional - sage.rings.number_field
            sage: M = R.ideal_monoid()                                                  # optional - sage.rings.number_field
            sage: M.has_coerce_map_from(R) # indirect doctest                           # optional - sage.rings.number_field
            True
            sage: M.has_coerce_map_from(QQ.ideal_monoid())                              # optional - sage.rings.number_field
            True
            sage: M.has_coerce_map_from(Zmod(6))                                        # optional - sage.rings.number_field
            False
            sage: M.has_coerce_map_from(loads(dumps(M)))                                # optional - sage.rings.number_field
            True
        """
        if isinstance(x, IdealMonoid_c):
            return self.ring().has_coerce_map_from(x.ring())
        else:
            return self.ring().has_coerce_map_from(x)

    def __eq__(self, other):
        r"""
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: R = QuadraticField(-23, 'a')                                          # optional - sage.rings.number_field
            sage: M = R.ideal_monoid()                                                  # optional - sage.rings.number_field
            sage: M == QQ                                                               # optional - sage.rings.number_field
            False
            sage: M == 17                                                               # optional - sage.rings.number_field
            False
            sage: M == R.ideal_monoid()                                                 # optional - sage.rings.number_field
            True
        """
        if not isinstance(other, IdealMonoid_c):
            return False
        else:
            return self.ring() == other.ring()

    def __ne__(self, other):
        r"""
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: R = QuadraticField(-23, 'a')                                          # optional - sage.rings.number_field
            sage: M = R.ideal_monoid()                                                  # optional - sage.rings.number_field
            sage: M != QQ                                                               # optional - sage.rings.number_field
            True
            sage: M != 17                                                               # optional - sage.rings.number_field
            True
            sage: M != R.ideal_monoid()                                                 # optional - sage.rings.number_field
            False
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: R = QuadraticField(-23, 'a')                                          # optional - sage.rings.number_field
            sage: M = R.ideal_monoid()                                                  # optional - sage.rings.number_field
            sage: hash(M) == hash(QQ)                                                   # optional - sage.rings.number_field
            False
            sage: hash(M) == 17                                                         # optional - sage.rings.number_field
            False
            sage: hash(M) == hash(R.ideal_monoid())                                     # optional - sage.rings.number_field
            True
        """
        # uses a random number, to have a distinct hash
        return hash((1580963238588124931699, self.ring()))
