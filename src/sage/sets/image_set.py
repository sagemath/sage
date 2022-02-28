"""
Image Sets
"""

# ****************************************************************************
#       Copyright (C) 2008      Mike Hansen <mhansen@gmail.com>
#                     2012      Christian Stump
#                     2020-2021 Frédéric Chapoton
#                     2021      Travis Scrimshaw
#                     2021      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from typing import Iterator

from sage.structure.parent import Parent, is_Parent
from sage.categories.map import is_Map
from sage.categories.morphism import IdentityMorphism
from sage.categories.poor_man_map import PoorManMap
from sage.categories.sets_cat import Sets
from sage.categories.enumerated_sets import EnumeratedSets
from sage.rings.integer import Integer
from sage.modules.free_module import FreeModule
from sage.symbolic.callable import is_CallableSymbolicExpression

from .set import Set_base, Set_add_sub_operators, Set_boolean_operators

class ImageSubobject(Parent):
    r"""
    The subset defined as the image of another set under a fixed map.

    Let `f: X \to Y` be a function. Then the image of `f` is defined as

    .. MATH::

        \{ f(x) | x \in X \} \subseteq Y.
    """
    def __init__(self, map, domain_subset, *, category=None, is_injective=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ, [0,1,2,3])
            sage: R.<x,y> = QQ[]
            sage: H = Hom(M, R, category=Sets())
            sage: f = H(lambda v: v[0]*x + v[1]*(x^2-y) + v[2]^2*(y+2) + v[3] - v[0]^2)
            sage: Im = f.image()
            sage: TestSuite(Im).run(skip=['_test_an_element', '_test_pickling',
            ....:                         '_test_some_elements'])
        """
        if not is_Parent(domain_subset):
            from sage.sets.set import Set
            domain_subset = Set(domain_subset)

        if not is_Map(map) and not isinstance(map, PoorManMap):
            map_name = f"The map {map}"
            if is_CallableSymbolicExpression(map):
                domain = map.parent().base()
                if len(map.arguments()) != 1:
                    domain = FreeModule(domain, len(map.arguments()))
                function = map
                def map(arg):
                    return function(*arg)
            else:
                domain = domain_subset
            map = PoorManMap(map, domain, name=map_name)

        if is_Map(map):
            map_category = map.category_for()
            if is_injective is None:
                try:
                    is_injective = map.is_injective()
                except NotImplementedError:
                    is_injective = False
        else:
            map_category = Sets()
            if is_injective is None:
                is_injective = False

        if category is None:
            category = map_category._meet_(domain_subset.category())

        category = category.Subobjects()
        if domain_subset in Sets().Finite() or map.codomain() in Sets().Finite():
            category = category.Finite()
        elif is_injective and domain_subset in Sets.Infinite():
            category = category.Infinite()

        if domain_subset in EnumeratedSets():
            category = category & EnumeratedSets()

        Parent.__init__(self, category=category)

        self._map = map
        self._domain_subset = domain_subset
        self._is_injective = is_injective

    def ambient(self):
        """
        Return the ambient set of ``self``, which is the codomain of
        the defining map.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(QQ, [0, 1, 2, 3])
            sage: R.<x,y> = ZZ[]
            sage: H = Hom(M, R, category=Sets())
            sage: f = H(lambda v: floor(v[0])*x + ceil(v[3] - v[0]^2))
            sage: Im = f.image()
            sage: Im.ambient() is R
            True

            sage: P = Partitions(3).map(attrcall('conjugate'))
            sage: P.ambient() is None
            True

            sage: R = Permutations(10).map(attrcall('reduced_word'))
            sage: R.ambient() is None
            True
        """
        return self._map.codomain()

    def lift(self, x):
        r"""
        Return the lift ``x`` to the ambient space, which is ``x``.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(QQ, [0, 1, 2, 3])
            sage: R.<x,y> = ZZ[]
            sage: H = Hom(M, R, category=Sets())
            sage: f = H(lambda v: floor(v[0])*x + ceil(v[3] - v[0]^2))
            sage: Im = f.image()
            sage: p = Im.lift(Im.an_element()); p
            2*x - 4
            sage: p.parent() is R
            True
        """
        return x

    def retract(self, x):
        """
        Return the retract of ``x`` from the ambient space, which is ``x``.

        .. WARNING::

            This does not check that ``x`` is actually in the image.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(QQ, [0, 1, 2, 3])
            sage: R.<x,y> = ZZ[]
            sage: H = Hom(M, R, category=Sets())
            sage: f = H(lambda v: floor(v[0])*x + ceil(v[3] - v[0]^2))
            sage: Im = f.image()
            sage: p = 2 * x - 4
            sage: Im.retract(p).parent()
            Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return x

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: Partitions(3).map(attrcall('conjugate'))
            Image of Partitions of the integer 3 by
             The map *.conjugate() from Partitions of the integer 3
        """
        return f"Image of {self._domain_subset} by {self._map}"

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: R = Permutations(10).map(attrcall('reduced_word'))
            sage: R.cardinality()
            3628800
        """
        if self._is_injective:
            return self._domain_subset.cardinality()
        return super().cardinality()

    def __iter__(self) -> Iterator:
        r"""
        Return an iterator over the elements of ``self``.

        EXAMPLES::

            sage: P = Partitions()
            sage: H = Hom(P, ZZ)
            sage: f = H(ZZ.sum)
            sage: X = f.image()
            sage: it = iter(X)
            sage: [next(it) for _ in range(5)]
            [0, 1, 2, 3, 4]
        """
        if self._is_injective:
            for x in self._domain_subset:
                yield self._map(x)
        else:
            visited = set()
            for x in self._domain_subset:
                y = self._map(x)
                if y in visited:
                    continue
                visited.add(y)
                yield y

    def _an_element_(self):
        r"""
        Return an element of this set.

        EXAMPLES::

            sage: R = SymmetricGroup(10).map(attrcall('reduced_word'))
            sage: R.an_element()
            [9, 8, 7, 6, 5, 4, 3, 2]
        """
        domain_element = self._domain_subset.an_element()
        return self._map(domain_element)

    def _sympy_(self):
        r"""
        Return an instance of a subclass of SymPy ``Set`` corresponding to ``self``.

        EXAMPLES::

            sage: from sage.sets.image_set import ImageSet
            sage: S = ImageSet(sin, RealSet.open(0, pi/4)); S
            Image of (0, 1/4*pi) by The map sin from (0, 1/4*pi)
            sage: S._sympy_()
            ImageSet(Lambda(x, sin(x)), Interval.open(0, pi/4))
        """
        from sympy import imageset
        try:
            sympy_map = self._map._sympy_()
        except AttributeError:
            sympy_map = self._map
        return imageset(sympy_map,
                        self._domain_subset._sympy_())


class ImageSet(ImageSubobject, Set_base, Set_add_sub_operators, Set_boolean_operators):
    r"""
    Image of a set by a map.

    EXAMPLES::

        sage: from sage.sets.image_set import ImageSet

    Symbolics::

        sage: ImageSet(sin, RealSet.open(0, pi/4))
        Image of (0, 1/4*pi) by The map sin from (0, 1/4*pi)
        sage: _.an_element()
        1/2*sqrt(-sqrt(2) + 2)

        sage: sos(x,y) = x^2 + y^2; sos
        (x, y) |--> x^2 + y^2
        sage: ImageSet(sos, ZZ^2)
        Image of
         Ambient free module of rank 2 over the principal ideal domain Integer Ring by
         The map (x, y) |--> x^2 + y^2 from Vector space of dimension 2 over Symbolic Ring
        sage: _.an_element()
        1
        sage: ImageSet(sos, Set([(3, 4), (3, -4)]))
        Image of {...(3, -4)...} by
         The map (x, y) |--> x^2 + y^2 from Vector space of dimension 2 over Symbolic Ring
        sage: _.an_element()
        25
    """
    pass

