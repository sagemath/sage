# sage_setup: distribution = sagemath-combinat
r"""
Free abelian monoids

AUTHORS:

- David Kohel (2005-09)

Sage supports free abelian monoids on any prescribed finite number
`n\geq 0` of generators. Use the
``FreeAbelianMonoid`` function to create a free abelian
monoid, and the ``gen`` and ``gens``
functions to obtain the corresponding generators. You can print the
generators as arbitrary strings using the optional
``names`` argument to the
``FreeAbelianMonoid`` function.

EXAMPLE 1: It is possible to create an abelian monoid in zero or
more variables; the syntax T(1) creates the monoid identity
element even in the rank zero case.

::

    sage: T = FreeAbelianMonoid(0, '')
    sage: T
    Free abelian monoid on 0 generators ()
    sage: T.gens()
    ()
    sage: T(1)
    1

EXAMPLE 2: A free abelian monoid uses a multiplicative
representation of elements, but the underlying representation is
lists of integer exponents.

::

    sage: F = FreeAbelianMonoid(5,names='a,b,c,d,e')
    sage: (a,b,c,d,e) = F.gens()
    sage: a*b^2*e*d
    a*b^2*d*e
    sage: x = b^2*e*d*a^7
    sage: x
    a^7*b^2*d*e
    sage: x.list()
    [7, 2, 0, 1, 1]
"""
# ****************************************************************************
#       Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.category_object import normalize_names
from sage.structure.parent import Parent
from sage.categories.monoids import Monoids
from .free_abelian_monoid_element import FreeAbelianMonoidElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

from sage.structure.factory import UniqueFactory


class FreeAbelianMonoidFactory(UniqueFactory):
    """
    Create the free abelian monoid in `n` generators.

    INPUT:

    - ``n`` -- integer

    - ``names`` -- names of generators

    OUTPUT: free abelian monoid

    EXAMPLES::

        sage: FreeAbelianMonoid(0, '')
        Free abelian monoid on 0 generators ()
        sage: F = FreeAbelianMonoid(5,names = list("abcde"))
        sage: F
        Free abelian monoid on 5 generators (a, b, c, d, e)
        sage: F(1)
        1
        sage: (a, b, c, d, e) = F.gens()
        sage: mul([ a, b, a, c, b, d, c, d ], F(1))
        a^2*b^2*c^2*d^2
        sage: a**2 * b**3 * a**2 * b**4
        a^4*b^7

    ::

        sage: loads(dumps(F)) is F
        True
    """
    def create_key(self, n, names):
        n = int(n)
        names = normalize_names(n, names)
        return (n, names)

    def create_object(self, version, key):
        return FreeAbelianMonoid_class(*key)


FreeAbelianMonoid_factory = FreeAbelianMonoidFactory("sage.monoids.free_abelian_monoid.FreeAbelianMonoid_factory")


def FreeAbelianMonoid(index_set=None, names=None, **kwds):
    r"""
    Return a free abelian monoid on `n` generators or with the generators
    indexed by a set `I`.

    We construct free abelian monoids by specifying either:

    - the number of generators and/or the names of the generators
    - the indexing set for the generators (this ignores the other two inputs)

    INPUT:

    - ``index_set`` -- an indexing set for the generators; if an integer,
      then this becomes `\{0, 1, \ldots, n-1\}`

    - ``names`` -- names of generators

    OUTPUT: a free abelian monoid

    EXAMPLES::

        sage: F.<a,b,c,d,e> = FreeAbelianMonoid(); F
        Free abelian monoid on 5 generators (a, b, c, d, e)
        sage: FreeAbelianMonoid(index_set=ZZ)
        Free abelian monoid indexed by Integer Ring
        sage: FreeAbelianMonoid(names='x,y')
        Free abelian monoid on 2 generators (x, y)
    """
    if isinstance(index_set, str):  # Swap args (this works if names is None as well)
        names, index_set = index_set, names

    if index_set is None and names is not None:
        if isinstance(names, str):
            index_set = names.count(',') + 1
        else:
            index_set = len(names)

    if index_set not in ZZ:
        if names is not None:
            names = normalize_names(len(names), names)
        from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
        return IndexedFreeAbelianMonoid(index_set, names=names, **kwds)

    if names is None:
        raise ValueError("names must be specified")
    return FreeAbelianMonoid_factory(index_set, names)


def is_FreeAbelianMonoid(x):
    """
    Return ``True`` if `x` is a free abelian monoid.

    EXAMPLES::

        sage: from sage.monoids.free_abelian_monoid import is_FreeAbelianMonoid
        sage: is_FreeAbelianMonoid(5)
        doctest:warning...
        DeprecationWarning: the function is_FreeAbelianMonoid is deprecated;
        use 'isinstance(..., FreeAbelianMonoid_class)' instead
        See https://github.com/sagemath/sage/issues/37897 for details.
        False
        sage: is_FreeAbelianMonoid(FreeAbelianMonoid(7,'a'))
        True
        sage: is_FreeAbelianMonoid(FreeMonoid(7,'a'))
        False
        sage: is_FreeAbelianMonoid(FreeMonoid(0,''))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37897, "the function is_FreeAbelianMonoid is deprecated; use 'isinstance(..., FreeAbelianMonoid_class)' instead")
    return isinstance(x, FreeAbelianMonoid_class)


class FreeAbelianMonoid_class(Parent):
    """
    Free abelian monoid on `n` generators.
    """
    Element = FreeAbelianMonoidElement

    def __init__(self, n, names):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(6,'b')
            sage: TestSuite(F).run()
        """
        if not isinstance(n, (int, Integer)):
            raise TypeError("n (=%s) must be an integer" % n)
        if n < 0:
            raise ValueError("n (=%s) must be nonnegative" % n)
        self.__ngens = int(n)
        assert names is not None
        Parent.__init__(self, names=names, category=Monoids().Commutative())

    def __repr__(self):
        n = self.__ngens
        return f"Free abelian monoid on {n} generators {self.gens()}"

    def __call__(self, x):
        """
        Create an element of this abelian monoid from `x`.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(10,'x')
            sage: F(F.gen(2))
            x2
            sage: F(1)
            1
        """
        if isinstance(x, FreeAbelianMonoidElement) and x.parent() == self:
            return x
        return self.element_class(self, x)

    def __contains__(self, x):
        """
        Return ``True`` if `x` is an element of this abelian monoid.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(10,'b')
            sage: F.gen(2)*F.gen(3) in F
            True

        Note that a monoid on `9` generators is not considered a
        submonoid of one on `10` generators.

        ::

            sage: FreeAbelianMonoid(9,'c').gen(2) in F
            False

        However, multiple calls to the monoid constructor do *not* return
        multiple distinct monoids.

        ::

            sage: FreeAbelianMonoid(10,'b').gen(2) in F
            True
        """
        return isinstance(x, FreeAbelianMonoidElement) and x.parent() == self

    def gen(self, i=0):
        """
        The `i`-th generator of the abelian monoid.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5,'a')
            sage: F.gen(0)
            a0
            sage: F.gen(2)
            a2
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError(f"argument i (= {i}) must be between 0 and {n-1}")
        x = [0 for j in range(n)]
        x[int(i)] = 1
        return self.element_class(self, x)

    @cached_method
    def gens(self) -> tuple:
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5,'a')
            sage: F.gens()
            (a0, a1, a2, a3, a4)
        """
        return tuple(self.gen(i) for i in range(self.__ngens))

    def ngens(self):
        """
        The number of free generators of the abelian monoid.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(3000, 'a')
            sage: F.ngens()
            3000
        """
        return self.__ngens

    def cardinality(self):
        r"""
        Return the cardinality of ``self``, which is `\infty`.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(3000, 'a')
            sage: F.cardinality()
            +Infinity
        """
        if self.__ngens == 0:
            from sage.rings.integer_ring import ZZ
            return ZZ.one()
        from sage.rings.infinity import infinity
        return infinity
