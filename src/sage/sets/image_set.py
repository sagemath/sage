"""
Image Sets
"""

# ****************************************************************************
#       Copyright (C) 2008      Mike Hansen <mhansen@gmail.com>
#                     2012      Christian Stump
#                     2020-2021 Frédéric Chapoton
#                     2021      Travis Scrimshaw
#                     2021-2024 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Iterator

from sage.categories.map import Map
from sage.categories.poor_man_map import PoorManMap
from sage.categories.sets_cat import Sets
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.structure.element import Expression
from sage.structure.parent import Parent

from .set import Set_base, Set_add_sub_operators, Set_boolean_operators


class ImageSubobject(Parent):
    r"""
    The subset defined as the image of another set under a fixed map.

    Let `f: X \to Y` be a function. Then the image of `f` is defined as

    .. MATH::

        \{ f(x) | x \in X \} \subseteq Y.

    INPUT:

    - ``map`` -- a function

    - ``domain_subset`` -- the set `X`; optional if `f` has a domain

    - ``is_injective`` -- whether the ``map`` is injective:
      - ``None`` (default): infer from ``map`` or default to ``False``
      - ``False``: do not assume that ``map`` is injective
      - ``True``: ``map`` is known to be injective
      - ``'check'``: raise an error when ``map`` is not injective

    - ``inverse`` -- a function (optional); a map from `f(X)` to `X`

    EXAMPLES::

        sage: import itertools
        sage: from sage.sets.image_set import ImageSubobject
        sage: D = ZZ
        sage: I = ImageSubobject(abs, ZZ, is_injective='check')
        sage: list(itertools.islice(I, 10))
        Traceback (most recent call last):
        ...
        ValueError: The map <built-in function abs> from Integer Ring is not injective: 1
    """
    def __init__(self, map, domain_subset, *, category=None, is_injective=None, inverse=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: # needs sage.modules
            sage: M = CombinatorialFreeModule(ZZ, [0,1,2,3])
            sage: R.<x,y> = QQ[]
            sage: H = Hom(M, R, category=Sets())
            sage: f = H(lambda v: v[0]*x + v[1]*(x^2-y) + v[2]^2*(y+2) + v[3] - v[0]^2)
            sage: Im = f.image()
            sage: TestSuite(Im).run(skip=['_test_an_element', '_test_pickling',
            ....:                         '_test_some_elements', '_test_elements'])
        """
        if not isinstance(domain_subset, Parent):
            from sage.sets.set import Set
            domain_subset = Set(domain_subset)

        if not isinstance(map, Map) and not isinstance(map, PoorManMap):
            map_name = f"The map {map}"
            if isinstance(map, Expression) and map.is_callable():
                domain = map.parent().base()
                if len(map.arguments()) != 1:
                    from sage.modules.free_module import FreeModule
                    domain = FreeModule(domain, len(map.arguments()))
                function = map

                def map(arg):
                    return function(*arg)
            else:
                domain = domain_subset
            map = PoorManMap(map, domain, name=map_name)

        if isinstance(map, Map):
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
        self._inverse = inverse
        self._domain_subset = domain_subset
        self._is_injective = is_injective

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: from sage.sets.image_set import ImageSubobject
            sage: D = ZZ
            sage: def f(x):
            ....:     return 2 * x
            sage: I = ImageSubobject(f, ZZ)
            sage: I == ImageSubobject(f, ZZ)
            True

        This method does not take into account whether an inverse is provided,
        injectivity is declared, or the category::

            sage: def f_inv(y):
            ....:     return y // 2
            sage: I == ImageSubobject(f, ZZ, inverse=f_inv)
            True
            sage: I == ImageSubobject(f, ZZ, is_injective=True)
            True
            sage: I.category()
            Category of enumerated subobjects of sets
            sage: I == ImageSubobject(f, ZZ, category=EnumeratedSets().Infinite())
            True
        """
        if not isinstance(other, ImageSubobject):
            return False
        return (self._map == other._map
                and self._domain_subset == other._domain_subset)

    def __ne__(self, other):
        r"""
        EXAMPLES::

            sage: from sage.sets.image_set import ImageSubobject
            sage: D = ZZ
            sage: def f(x):
            ....:     return 2 * x
            sage: I = ImageSubobject(f, ZZ)
            sage: I != ImageSubobject(f, QQ)
            True
        """
        return not (self == other)

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.sets.image_set import ImageSubobject
            sage: def f(x):
            ....:     return 2 * x
            sage: I = ImageSubobject(f, ZZ)
            sage: hash(I) == hash(ImageSubobject(f, ZZ))
            True
        """
        return hash((self._map, self._domain_subset))

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: import itertools
            sage: from sage.sets.image_set import ImageSubobject
            sage: D = ZZ
            sage: I = ImageSubobject(lambda x: 2 * x, ZZ, inverse=lambda x: x/2)
            sage: I(8/2)
            4
            sage: _.parent()
            Integer Ring
            sage: I(10/2)
            Traceback (most recent call last):
            ...
            ValueError: 5 is not in Image of Integer Ring by The map <function <lambda> at ...> from Integer Ring
            sage: 6 in I
            True
            sage: 7 in I
            False
        """
        # Same as ImageManifoldSubset.__contains__
        codomain = self._map.codomain()
        if codomain is not None and x not in codomain:
            raise ValueError(f"{x} is not in {self}")
        if self._inverse is not None:
            preimage = self._inverse(x)
            if preimage not in self._domain_subset:
                raise ValueError(f"{x} is not in {self}")
            preimage = self._map.domain()(preimage)
            y = self._map(preimage)
            if y == x:
                return y
            raise ValueError(f"{x} is not in {self}")
        raise NotImplementedError

    def ambient(self):
        """
        Return the ambient set of ``self``, which is the codomain of
        the defining map.

        EXAMPLES::

            sage: # needs sage.modules
            sage: M = CombinatorialFreeModule(QQ, [0, 1, 2, 3])
            sage: R.<x,y> = ZZ[]
            sage: H = Hom(M, R, category=Sets())
            sage: f = H(lambda v: floor(v[0])*x + ceil(v[3] - v[0]^2))
            sage: Im = f.image()
            sage: Im.ambient() is R
            True

            sage: P = Partitions(3).map(attrcall('conjugate'))                          # needs sage.combinat
            sage: P.ambient() is None                                                   # needs sage.combinat
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

            sage: # needs sage.modules
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

            sage: # needs sage.modules
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

            sage: Partitions(3).map(attrcall('conjugate'))                              # needs sage.combinat
            Image of Partitions of the integer 3 by
             The map *.conjugate() from Partitions of the integer 3
        """
        return f"Image of {self._domain_subset} by {self._map}"

    @cached_method
    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        EXAMPLES:

        Injective case (note that
        :meth:`~sage.categories.enumerated_sets.EnumeratedSets.ParentMethods.map`
        defaults to ``is_injective=True``):

            sage: R = Permutations(10).map(attrcall('reduced_word'))
            sage: R.cardinality()
            3628800

            sage: Evens = ZZ.map(lambda x: 2 * x)
            sage: Evens.cardinality()
            +Infinity

        Non-injective case::

            sage: Z7 = Set(range(7))
            sage: from sage.sets.image_set import ImageSet
            sage: Z4711 = ImageSet(lambda x: x**4 % 11, Z7, is_injective=False)
            sage: Z4711.cardinality()
            6

            sage: Squares = ImageSet(lambda x: x^2, ZZ, is_injective=False,
            ....:                    category=Sets().Infinite())
            sage: Squares.cardinality()
            +Infinity

            sage: Mod2 = ZZ.map(lambda x: x % 2, is_injective=False)
            sage: Mod2.cardinality()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot determine cardinality of a
            non-injective image of an infinite set
        """
        domain_cardinality = self._domain_subset.cardinality()
        if self._is_injective and self._is_injective != 'check':
            return domain_cardinality
        if self in Sets().Infinite():
            return Infinity
        if domain_cardinality == Infinity:
            raise NotImplementedError('cannot determine cardinality of a non-injective image of an infinite set')
        # Fallback like EnumeratedSets.ParentMethods.__len__
        return Integer(len(list(iter(self))))

    def __iter__(self) -> Iterator:
        r"""
        Return an iterator over the elements of ``self``.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = Partitions()
            sage: H = Hom(P, ZZ)
            sage: f = H(ZZ.sum)
            sage: X = f.image()
            sage: it = iter(X)
            sage: [next(it) for _ in range(5)]
            [0, 1, 2, 3, 4]
        """
        if self._is_injective and self._is_injective != 'check':
            for x in self._domain_subset:
                yield self._map(x)
        else:
            visited = set()
            for x in self._domain_subset:
                y = self._map(x)
                if y in visited:
                    if self._is_injective == 'check':
                        raise ValueError(f'{self._map} is not injective: {y}')
                    continue
                visited.add(y)
                yield y

    def _an_element_(self):
        r"""
        Return an element of this set.

        EXAMPLES::

            sage: R = SymmetricGroup(10).map(attrcall('reduced_word'))                  # needs sage.groups
            sage: R.an_element()                                                        # needs sage.groups
            [9, 8, 7, 6, 5, 4, 3, 2]
        """
        domain_element = self._domain_subset.an_element()
        return self._map(domain_element)

    def _sympy_(self):
        r"""
        Return an instance of a subclass of SymPy ``Set`` corresponding to ``self``.

        EXAMPLES::

            sage: from sage.sets.image_set import ImageSet
            sage: S = ImageSet(sin, RealSet.open(0, pi/4)); S                           # needs sage.symbolic
            Image of (0, 1/4*pi) by The map sin from (0, 1/4*pi)
            sage: S._sympy_()                                                           # needs sage.symbolic
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

        sage: ImageSet(sin, RealSet.open(0, pi/4))                                      # needs sage.symbolic
        Image of (0, 1/4*pi) by The map sin from (0, 1/4*pi)
        sage: _.an_element()                                                            # needs sage.symbolic
        1/2*sqrt(-sqrt(2) + 2)

        sage: # needs sage.symbolic
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
