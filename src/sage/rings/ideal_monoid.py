# sage_setup: distribution = sagemath-categories
"""
Monoid of ideals in a commutative ring

WARNING: This is used by some rings that are not commutative! ::

    sage: MS = MatrixSpace(QQ, 3, 3)                                                    # needs sage.modules
    sage: type(MS.ideal(MS.one()).parent())                                             # needs sage.modules
    <class 'sage.rings.ideal_monoid.IdealMonoid_c_with_category'>
"""

import sage.rings.integer_ring
from sage.categories.monoids import Monoids
from sage.rings import ideal
from sage.structure.parent import Parent


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

            sage: R = QuadraticField(-23, 'a')                                          # needs sage.rings.number_field
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R); M  # indirect doctest                             # needs sage.rings.number_field
            Monoid of ideals of Number Field in a with defining polynomial x^2 + 23
             with a = 4.795831523312720?*I

            sage: id = QQ.ideal(6)
            sage: id.parent().category()
            Category of commutative monoids

            sage: MS = MatrixSpace(QQ, 3, 3)                                            # needs sage.modules
            sage: MS.ideal(MS.one()).parent().category()                                # needs sage.modules
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

            sage: R = QuadraticField(-23, 'a')                                          # needs sage.rings.number_field
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R); M._repr_()                                        # needs sage.rings.number_field
            'Monoid of ideals of Number Field in a with defining polynomial x^2 + 23
             with a = 4.795831523312720?*I'
        """
        return "Monoid of ideals of %s" % self.__R

    def ring(self):
        r"""
        Return the ring of which this is the ideal monoid.

        EXAMPLES::

            sage: R = QuadraticField(-23, 'a')                                          # needs sage.rings.number_field
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R); M.ring() is R                                     # needs sage.rings.number_field
            True
        """
        return self.__R

    def _element_constructor_(self, x):
        r"""
        Create an ideal in this monoid from ``x``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<a> = QuadraticField(-23)
            sage: from sage.rings.ideal_monoid import IdealMonoid
            sage: M = IdealMonoid(R)
            sage: M(a)   # indirect doctest
            Fractional ideal (a)
            sage: M([a-4, 13])
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

            sage: # needs sage.rings.number_field
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

    def __eq__(self, other):
        r"""
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
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
            return False
        else:
            return self.ring() == other.ring()

    def __ne__(self, other):
        r"""
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R = QuadraticField(-23, 'a')
            sage: M = R.ideal_monoid()
            sage: M != QQ
            True
            sage: M != 17
            True
            sage: M != R.ideal_monoid()
            False
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R = QuadraticField(-23, 'a')
            sage: M = R.ideal_monoid()
            sage: hash(M) == hash(QQ)
            False
            sage: hash(M) == 17
            False
            sage: hash(M) == hash(R.ideal_monoid())
            True
        """
        # uses a random number, to have a distinct hash
        return hash((1580963238588124931699, self.ring()))
