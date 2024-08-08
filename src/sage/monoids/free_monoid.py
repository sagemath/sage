r"""
Free Monoids

AUTHORS:

- David Kohel (2005-09)
- Simon King (2011-04): Put free monoids into the category framework

Sage supports free monoids on any prescribed finite number
`n\geq 0` of generators. Use the ``FreeMonoid``
function to create a free monoid, and the ``gen`` and
``gens`` functions to obtain the corresponding
generators. You can print the generators as arbitrary strings using
the optional ``names`` argument to the
``FreeMonoid`` function.
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

from sage.rings.integer import Integer
from sage.structure.category_object import normalize_names
from .free_monoid_element import FreeMonoidElement

from .monoid import Monoid_class

from sage.combinat.words.finite_word import FiniteWord_class

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ


def is_FreeMonoid(x):
    """
    Return ``True`` if `x` is a free monoid.

    EXAMPLES::

        sage: from sage.monoids.free_monoid import is_FreeMonoid
        sage: is_FreeMonoid(5)
        doctest:warning...
        DeprecationWarning: the function is_FreeMonoid is deprecated;
        use 'isinstance(..., (FreeMonoid, IndexedFreeMonoid))' instead
        See https://github.com/sagemath/sage/issues/37897 for details.
        False
        sage: is_FreeMonoid(FreeMonoid(7,'a'))
        True
        sage: is_FreeMonoid(FreeAbelianMonoid(7,'a'))
        False
        sage: is_FreeMonoid(FreeAbelianMonoid(0,''))
        False
        sage: is_FreeMonoid(FreeMonoid(index_set=ZZ))
        True
        sage: is_FreeMonoid(FreeAbelianMonoid(index_set=ZZ))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37897, "the function is_FreeMonoid is deprecated; use 'isinstance(..., (FreeMonoid, IndexedFreeMonoid))' instead")
    if isinstance(x, FreeMonoid):
        return True
    from sage.monoids.indexed_free_monoid import IndexedFreeMonoid
    return isinstance(x, IndexedFreeMonoid)


class FreeMonoid(Monoid_class, UniqueRepresentation):
    r"""
    Return a free monoid on `n` generators or with the generators
    indexed by a set `I`.

    We construct free monoids by specifing either:

    - the number of generators and/or the names of the generators
    - the indexing set for the generators

    INPUT:

    - ``index_set`` -- an indexing set for the generators; if an
      integer `n`, than this becomes `\{0, 1, \ldots, n-1\}`

    - ``names`` -- names of generators

    - ``commutative`` -- boolean (default: ``False``); whether the free
      monoid is commutative or not

    OUTPUT: a free monoid

    EXAMPLES::

        sage: F = FreeMonoid(3,'x'); F
        Free monoid on 3 generators (x0, x1, x2)
        sage: x = F.gens()
        sage: x[0]*x[1]**5 * (x[0]*x[2])
        x0*x1^5*x0*x2
        sage: F = FreeMonoid(3, 'a')
        sage: F
        Free monoid on 3 generators (a0, a1, a2)

        sage: F.<a,b,c,d,e> = FreeMonoid(); F
        Free monoid on 5 generators (a, b, c, d, e)
        sage: FreeMonoid(index_set=ZZ)
        Free monoid indexed by Integer Ring

        sage: F.<x,y,z> = FreeMonoid(abelian=True); F
        Free abelian monoid on 3 generators (x, y, z)
        sage: FreeMonoid(index_set=ZZ, commutative=True)
        Free abelian monoid indexed by Integer Ring
    """
    @staticmethod
    def __classcall_private__(cls, index_set=None, names=None,
                              commutative=False, **kwds):
        r"""
        Construct a free monoid or a free abelian monoid, depending on the
        input. Also, normalize the input.

        EXAMPLES::

            sage: F.<a,b,c,d,e> = FreeMonoid(); F
            Free monoid on 5 generators (a, b, c, d, e)
            sage: FreeMonoid(index_set=ZZ)
            Free monoid indexed by Integer Ring
            sage: F.<x,y,z> = FreeMonoid(abelian=True); F
            Free abelian monoid on 3 generators (x, y, z)
            sage: FreeMonoid(index_set=ZZ, commutative=True)
            Free abelian monoid indexed by Integer Ring
            sage: F = FreeMonoid(index_set=ZZ, names='x,y,z')
            sage: G = FreeMonoid(index_set=ZZ, names=['x', 'y', 'z'])
            sage: F == G
            True
            sage: F is G
            True

            sage: FreeMonoid(2, names='a,b') is FreeMonoid(names=['a','b'])
            True

        Fix a bug when ``index_set`` is ``None`` and ``names`` is a
        string (:issue:`26221`)::

            sage: FreeMonoid(2, names=['a','b']) is FreeMonoid(names='a,b')
            True
        """
        if 'abelian' in kwds:
            commutative = kwds.pop('abelian')

        if commutative:
            from sage.monoids.free_abelian_monoid import FreeAbelianMonoid
            return FreeAbelianMonoid(index_set, names, **kwds)

        # Swap args (this works if names is None as well)
        if isinstance(index_set, str):
            names, index_set = index_set, names

        if index_set is None and names is not None:
            if isinstance(names, str):
                index_set = names.count(',') + 1
            else:
                index_set = len(names)

        if index_set not in ZZ:
            if names is not None:
                names = normalize_names(-1, names)
            from sage.monoids.indexed_free_monoid import IndexedFreeMonoid
            return IndexedFreeMonoid(index_set, names=names, **kwds)

        if names is None:
            raise ValueError("names must be specified")
        names = normalize_names(index_set, names)
        return super().__classcall__(cls, index_set, names)

    Element = FreeMonoidElement

    def __init__(self, n, names=None):
        """
        Return a free monoid on `n` generators or with the generators
        indexed by a set `I`.

        INPUT:

        - ``n`` -- number of generators

        - ``names`` -- names of generators

        TESTS::

            sage: M = FreeMonoid(3, names=['a','b','c'])
            sage: TestSuite(M).run()
            sage: F.<x,y> = FreeMonoid()
            sage: TestSuite(F).run()
        """
        if not isinstance(n, (int, Integer)):
            raise TypeError("n (=%s) must be an integer" % n)
        if n < 0:
            raise ValueError("n (=%s) must be nonnegative" % n)
        self.__ngens = int(n)
        Monoid_class.__init__(self, names)

    def _repr_(self):
        return f"Free monoid on {self.__ngens} generators {self.gens()}"

    def _element_constructor_(self, x, check=True):
        """
        Return ``x`` coerced into this free monoid.

        One can create a free monoid element from the integer 1, from
        a list of 2-tuples of integers `(i,j)`, where `(i,j)`
        corresponds to `x_i^j`, where `x_i` is the `i`-th generator,
        and from words in the same alphabet as the generators.

        EXAMPLES::

            sage: F = FreeMonoid(3, 'a')
            sage: F(1)
            1
            sage: F(F.gen(0))
            a0
            sage: F(0)
            Traceback (most recent call last):
            ...
            TypeError: argument x (= 0) is of the wrong type

        An example with a list::

            sage: F([(0,5),(1,2),(0,10),(0,2),(1,2)])
            a0^5*a1^2*a0^12*a1^2

        An example using words::

            sage: F = FreeMonoid(3, 'a,b,c')
            sage: w = Word('aabbcabac')
            sage: F(w)
            a^2*b^2*c*a*b*a*c
            sage: F(Word([]))
            1

        TESTS::

            sage: F(F(w), check=False)
            a^2*b^2*c*a*b*a*c
        """
        # There should really be some careful type checking here...
        if isinstance(x, FreeMonoidElement) and x.parent() is self:
            return x
        if isinstance(x, FreeMonoidElement) and x.parent() == self:
            return self.element_class(self, x._element_list, check)
        if isinstance(x, (int, Integer)) and x == 1:
            return self.element_class(self, x, check)
        if isinstance(x, FiniteWord_class):
            d = self.gens_dict()
            return self.prod([d[let] for let in x])
        if isinstance(x, list):
            return self.element_class(self, x, check)

        raise TypeError("argument x (= %s) is of the wrong type" % x)

    def __contains__(self, x):
        return isinstance(x, FreeMonoidElement) and x.parent() == self

    def gen(self, i=0):
        """
        The `i`-th generator of the monoid.

        INPUT:

        - ``i`` -- integer (default: 0)

        EXAMPLES::

            sage: F = FreeMonoid(3, 'a')
            sage: F.gen(1)
            a1
            sage: F.gen(2)
            a2
            sage: F.gen(5)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 5) must be between 0 and 2
        """
        n = self.__ngens
        if i < 0 or not i < n:
            msg = "argument i (= %s) must be between 0 and %s"
            raise IndexError(msg % (i, n - 1))
        return self.element_class(self, [(Integer(i), Integer(1))])

    def ngens(self):
        """
        The number of free generators of the monoid.

        EXAMPLES::

            sage: F = FreeMonoid(2005, 'a')
            sage: F.ngens()
            2005
        """
        return self.__ngens

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        This is `\infty` if there is at least one generator.

        EXAMPLES::

            sage: F = FreeMonoid(2005, 'a')
            sage: F.cardinality()
            +Infinity

            sage: F = FreeMonoid(0, [])
            sage: F.cardinality()
            1
        """
        if self.__ngens == 0:
            return Integer(1)
        from sage.rings.infinity import infinity
        return infinity
