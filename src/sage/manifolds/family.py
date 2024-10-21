# sage_setup: distribution = sagemath-symbolics
r"""
Families of Manifold Objects

The class :class:`ManifoldObjectFiniteFamily` is a subclass of :class:`FiniteFamily`
that provides an associative container of manifold objects, indexed by their
``_name`` attributes.

:class:`ManifoldObjectFiniteFamily` instances are totally ordered according
to their lexicographically ordered element names.

The subclass :class:`ManifoldSubsetFiniteFamily` customizes the print
representation further.

AUTHORS:

- Matthias Koeppe (2021): initial version
"""
#*****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functools import total_ordering
from sage.sets.family import FiniteFamily

@total_ordering
class ManifoldObjectFiniteFamily(FiniteFamily):

    r"""
    Finite family of manifold objects, indexed by their names.

    The class :class:`ManifoldObjectFiniteFamily` inherits from
    :class:`FiniteFamily`.  Therefore it is an associative container.

    It provides specialized ``__repr__`` and ``_latex_`` methods.

    :class:`ManifoldObjectFiniteFamily` instances are totally ordered
    according to their lexicographically ordered element names.

    EXAMPLES::

        sage: from sage.manifolds.family import ManifoldObjectFiniteFamily
        sage: M = Manifold(2, 'M', structure='topological')
        sage: A = M.subset('A')
        sage: B = M.subset('B')
        sage: C = B.subset('C')
        sage: F = ManifoldObjectFiniteFamily([A, B, C]); F
        Set {A, B, C} of objects of the 2-dimensional topological manifold M
        sage: latex(F)
        \{A, B, C\}
        sage: F['B']
        Subset B of the 2-dimensional topological manifold M

    All objects must have the same base manifold::

        sage: N = Manifold(2, 'N', structure='topological')
        sage: ManifoldObjectFiniteFamily([M, N])
        Traceback (most recent call last):
        ...
        TypeError: all objects must have the same manifold
    """
    def __init__(self, objects=(), keys=None):
        r"""
        Initialize a new instance of :class:`ManifoldObjectFiniteFamily`.

        TESTS:

            sage: from sage.manifolds.family import ManifoldObjectFiniteFamily
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: C = B.subset('C')
            sage: F = ManifoldObjectFiniteFamily([A, B, C]); F
            Set {A, B, C} of objects of the 2-dimensional topological manifold M
            sage: TestSuite(F).run(skip='_test_elements')

        Like ``frozenset``, it can be created from any iterable::

            sage: from sage.manifolds.family import ManifoldSubsetFiniteFamily
            sage: M = Manifold(2, 'M', structure='topological')
            sage: I = M.subset('I')
            sage: gen = (subset for subset in (M, I, M, I, M, I)); gen
            <generator object ...>
            sage: ManifoldSubsetFiniteFamily(gen)
            Set {I, M} of subsets of the 2-dimensional topological manifold M
        """
        if isinstance(objects, dict):
            dictionary = objects
        else:
            dictionary = {object._name: object for object in objects}
        if keys is None:
            keys = sorted(dictionary.keys())
        FiniteFamily.__init__(self, dictionary, keys)
        names_and_latex_names = sorted((object._name, object._latex_name)
                                       for object in self)
        self._name = '{' + ', '.join(keys) + '}'
        latex_names = (latex_name for name, latex_name in names_and_latex_names)
        self._latex_name = r'\{' + ', '.join(latex_names) + r'\}'
        try:
            object_iter = iter(self)
            self._manifold = next(object_iter)._manifold
        except StopIteration:
            self._manifold = None
        else:
            if not all(object._manifold == self._manifold for object in object_iter):
                raise TypeError(f'all {self._repr_object_type()} must have the same manifold')

    def _repr_object_type(self):
        r"""
        String that describes the type of the elements (plural).

        TESTS::

            sage: from sage.manifolds.family import ManifoldObjectFiniteFamily
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: ManifoldObjectFiniteFamily([A, B]).__repr__()           # indirect doctest
            'Set {A, B} of objects of the 2-dimensional topological manifold M'
        """
        return "objects"

    def __lt__(self, other):
        r"""
        Implement the total order on instances of :class:`ManifoldObjectFiniteFamily`.

        TESTS::

            sage: from sage.manifolds.family import ManifoldSubsetFiniteFamily
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: sorted([ManifoldSubsetFiniteFamily([A, B]), ManifoldSubsetFiniteFamily([]),
            ....:         ManifoldSubsetFiniteFamily([B]), ManifoldSubsetFiniteFamily([A])])
            [{},
             Set {A} of subsets of the 2-dimensional topological manifold M,
             Set {A, B} of subsets of the 2-dimensional topological manifold M,
             Set {B} of subsets of the 2-dimensional topological manifold M]
        """
        if not isinstance(other, ManifoldSubsetFiniteFamily):
            return NotImplemented
        return self.keys() < other.keys()

    def __repr__(self):
        r"""
        String representation of the object.

        TESTS::

            sage: from sage.manifolds.family import ManifoldObjectFiniteFamily
            sage: ManifoldObjectFiniteFamily().__repr__()
            '{}'
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: ManifoldObjectFiniteFamily([A, B]).__repr__()
            'Set {A, B} of objects of the 2-dimensional topological manifold M'
        """
        if self:
            return "Set {} of {} of the {}".format(self._name, self._repr_object_type(), self._manifold)
        else:
            return "{}"

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: from sage.manifolds.family import ManifoldSubsetFiniteFamily
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: ManifoldSubsetFiniteFamily([B, A])._latex_()
            '\\{A, B\\}'
        """
        return self._latex_name

class ManifoldSubsetFiniteFamily(ManifoldObjectFiniteFamily):

    r"""
    Finite family of subsets of a topological manifold, indexed by their names.

    The class :class:`ManifoldSubsetFiniteFamily` inherits from
    :class:`ManifoldObjectFiniteFamily`.  It provides an associative
    container with specialized ``__repr__`` and ``_latex_`` methods.

    :class:`ManifoldSubsetFiniteFamily` instances are totally ordered according
    to their lexicographically ordered element (subset) names.

    EXAMPLES::

        sage: from sage.manifolds.family import ManifoldSubsetFiniteFamily
        sage: M = Manifold(2, 'M', structure='topological')
        sage: A = M.subset('A')
        sage: B = M.subset('B')
        sage: C = B.subset('C')
        sage: ManifoldSubsetFiniteFamily([A, B, C])
        Set {A, B, C} of subsets of the 2-dimensional topological manifold M
        sage: latex(_)
        \{A, B, C\}

    All subsets must have the same base manifold::

        sage: N = Manifold(2, 'N', structure='topological')
        sage: ManifoldSubsetFiniteFamily([M, N])
        Traceback (most recent call last):
        ...
        TypeError: all open subsets must have the same manifold
    """

    @classmethod
    def from_subsets_or_families(cls, *subsets_or_families):
        r"""
        Construct a ``ManifoldSubsetFiniteFamily`` from given subsets or
        iterables of subsets.

        EXAMPLES::

            sage: from sage.manifolds.family import ManifoldSubsetFiniteFamily
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: Bs = (M.subset(f'B{i}') for i in range(5))
            sage: Cs = ManifoldSubsetFiniteFamily([M.subset('C0'), M.subset('C1')])
            sage: ManifoldSubsetFiniteFamily.from_subsets_or_families(A, Bs, Cs)
            Set {A, B0, B1, B2, B3, B4, C0, C1} of subsets of the 2-dimensional topological manifold M
        """
        def generate_subsets():
            from sage.manifolds.subset import ManifoldSubset
            for arg in subsets_or_families:
                if isinstance(arg, ManifoldSubset):
                    yield arg
                else:
                    # arg must be an iterable of ManifoldSubset instances
                    yield from arg
        return cls(generate_subsets())

    def _repr_object_type(self):
        r"""
        String that describes the type of the elements (plural).

        TESTS::

            sage: from sage.manifolds.family import ManifoldSubsetFiniteFamily
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: ManifoldSubsetFiniteFamily([A, B]).__repr__()           # indirect doctest
            'Set {A, B} of subsets of the 2-dimensional topological manifold M'
        """
        if all(subset.is_open() for subset in self):
            return "open subsets"
        else:
            return "subsets"
