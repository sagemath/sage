# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Subcrystals

These are the crystals that are subsets of a larger ambient crystal.

AUTHORS:

- Travis Scrimshaw (2013-10-16): Initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

import collections.abc

from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.crystals import Crystals
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.supercrystals import SuperCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.structure.richcmp import richcmp


class Subcrystal(UniqueRepresentation, Parent):
    r"""
    A subcrystal `X` of an ambient crystal `Y` is a crystal formed by taking a
    subset of `Y` and whose crystal structure is induced by `Y`.

    INPUT:

    - ``ambient`` -- the ambient crystal
    - ``contained`` -- (optional) a set (or function) which specifies when an
      element is contained in the subcrystal; the default is everything
      possible is included
    - ``generators`` -- (optional) the generators for the subcrystal; the
      default is the generators for the ambient crystal
    - ``virtualization``, ``scaling_factors`` -- (optional)
      dictionaries whose key `i` corresponds to the sets `\sigma_i`
      and `\gamma_i` respectively used to define virtual crystals; see
      :class:`~sage.combinat.crystals.virtual_crystal.VirtualCrystal`
    - ``cartan_type`` -- (optional) the Cartan type for the subcrystal; the
      default is the Cartan type for the ambient crystal
    - ``index_set`` -- (optional) the index set for the subcrystal; the
      default is the index set for the Cartan type
    - ``category`` -- (optional) the category for the subcrystal; the
      default is the :class:`~sage.categories.crystals.Crystals` category

    .. SEEALSO::

        :meth:`~sage.categories.crystals.Crystals.ParentMethods.subcrystal`

    EXAMPLES:

    We build out a subcrystal starting from an element and only going
    to the lowest weight::

        sage: B = crystals.Tableaux(['A',3], shape=[2,1])
        sage: S = B.subcrystal(generators=[B(3,1,2)], direction='lower')
        sage: S.cardinality()
        11

    Here we build out in both directions starting from an element, but we
    also have restricted ourselves to type `A_2`::

        sage: T = B.subcrystal(index_set=[1,2], generators=[B(3,1,1)])
        sage: T.cardinality()
        8
        sage: list(T)
        [[[1, 1], [3]],
         [[1, 2], [3]],
         [[1, 1], [2]],
         [[2, 2], [3]],
         [[1, 2], [2]],
         [[2, 3], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]]]

    Now we take the crystal corresponding to the intersection of
    the previous two subcrystals::

        sage: U = B.subcrystal(contained=lambda x: x in S and x in T, generators=B)
        sage: list(U)
        [[[2, 3], [3]], [[1, 2], [3]], [[2, 2], [3]]]

    .. TODO::

        Include support for subcrystals which only contains certain arrows.

    TESTS:

    Check that the subcrystal respects being in the category
    of supercrystals (:issue:`27368`)::

        sage: T = crystals.Tableaux(['A',[1,1]], [2,1])
        sage: S = T.subcrystal(max_depth=3)
        sage: S.category()
        Category of finite super crystals
    """
    @staticmethod
    def __classcall_private__(cls, ambient, contained=None, generators=None,
                              virtualization=None, scaling_factors=None,
                              cartan_type=None, index_set=None, category=None):
        """
        Normalize arguments to ensure a (relatively) unique representation.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',4], shape=[2,1])
            sage: S1 = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            sage: S2 = B.subcrystal(generators=[B(2,1,1), B(5,2,4)], cartan_type=['A',4], index_set=(1,2))
            sage: S1 is S2
            True
        """
        if isinstance(contained, (collections.abc.Sequence, collections.abc.Set)):
            contained = frozenset(contained)
        #elif contained in Sets():

        if cartan_type is None:
            cartan_type = ambient.cartan_type()
        else:
            cartan_type = CartanType(cartan_type)
        if index_set is None:
            index_set = cartan_type.index_set
        if generators is None:
            generators = ambient.module_generators

        category = Crystals().or_subcategory(category)
        if ambient in SuperCrystals():
            category = category & SuperCrystals()
        if ambient in FiniteCrystals() or isinstance(contained, frozenset):
            category = category.Finite()

        if virtualization is not None:
            if scaling_factors is None:
                scaling_factors = {i:1 for i in index_set}
            from sage.combinat.crystals.virtual_crystal import VirtualCrystal
            return VirtualCrystal(ambient, virtualization, scaling_factors, contained,
                                  generators, cartan_type, index_set, category)
        if scaling_factors is not None:
            # virtualization must be None
            virtualization = {i:(i,) for i in index_set}
            from sage.combinat.crystals.virtual_crystal import VirtualCrystal
            return VirtualCrystal(ambient, virtualization, scaling_factors, contained,
                                  generators, cartan_type, index_set, category)

        # We need to give these as optional arguments so it unpickles correctly
        return super().__classcall__(cls, ambient, contained,
                                     tuple(generators),
                                     cartan_type=cartan_type,
                                     index_set=tuple(index_set),
                                     category=category)

    def __init__(self, ambient, contained, generators, cartan_type, index_set, category):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            sage: TestSuite(S).run()
        """
        self._ambient = ambient
        self._contained = contained
        self._cardinality = None # ``None`` means currently unknown
        self._cartan_type = cartan_type
        self._index_set = tuple(index_set)
        Parent.__init__(self, category=category)
        self.module_generators = tuple(self.element_class(self, g) for g in generators
                                       if self._containing(g))

        if isinstance(contained, frozenset):
            self._cardinality = Integer(len(contained))
            self._list = [self.element_class(self, x) for x in contained]

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',4], shape=[2,1])
            sage: B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            Subcrystal of The crystal of tableaux of type ['A', 4] and shape(s) [[2, 1]]
        """
        return "Subcrystal of {}".format(self._ambient)

    @lazy_attribute
    def _containing(self):
        """
        Check if ``x`` is contained in ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            sage: S._containing(B(5,2,4))
            True
            sage: S._containing(B(4,2,4))
            True
        """
        if self._contained is None:
            return lambda x: True
        if isinstance(self._contained, frozenset):
            return self._contained.__contains__
        return self._contained # Otherwise it should be a function

    def __contains__(self, x):
        """
        Check if ``x`` is in ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            sage: B(5,2,4) in S
            True
            sage: mg = B.module_generators[0]
            sage: mg in S
            True
            sage: mg.f(2).f(3) in S
            False
        """
        if isinstance(x, Subcrystal.Element) and x.parent() == self:
            return True

        if x in self._ambient:
            if not self._containing(x):
                return False
            x = self.element_class(self, x)

        if self in FiniteCrystals():
            return x in self.list()

        # TODO: make this work for infinite crystals
        import warnings
        warnings.warn("Testing containment in an infinite crystal"
                      " defaults to returning True")
        return True

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=[B(2,1,1)], index_set=[1,2])
            sage: S.cardinality()
            8
            sage: B = crystals.infinity.Tableaux(['A',2])
            sage: S = B.subcrystal(max_depth=4)
            sage: S.cardinality()
            22

        TESTS:

        Check that :issue:`19481` is fixed::

            sage: from sage.combinat.crystals.virtual_crystal import VirtualCrystal
            sage: A = crystals.infinity.Tableaux(['A',3])
            sage: V = VirtualCrystal(A, {1:(1,3), 2:(2,)}, {1:1, 2:2}, cartan_type=['C',2])
            sage: V.cardinality()
            Traceback (most recent call last):
            ...
            NotImplementedError: unknown cardinality
        """
        if self._cardinality is not None:
            return self._cardinality

        try:
            card = Integer(len(self._list))
            self._cardinality = card
            return self._cardinality
        except AttributeError:
            if self in FiniteCrystals():
                return Integer(len(self.list()))
            try:
                card = super().cardinality()
            except AttributeError:
                raise NotImplementedError("unknown cardinality")
            if card == infinity:
                self._cardinality = card
                return card
            self._cardinality = Integer(len(self.list()))
            return self._cardinality

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            sage: S.index_set()
            (1, 2)
        """
        return self._index_set

    class Element(ElementWrapper):
        """
        An element of a subcrystal. Wraps an element in the ambient crystal.
        """

        def _richcmp_(self, other, op):
            """
            EXAMPLES:

            For == operator::

                sage: A = crystals.KirillovReshetikhin(['C',2,1], 1,2).affinization()
                sage: S = A.subcrystal(max_depth=2)
                sage: sorted(S)
                [[[1, 1]](-1),
                 [[1, 2]](-1),
                 [](0),
                 [[1, 1]](0),
                 [[1, 2]](0),
                 [[1, -2]](0),
                 [[2, 2]](0),
                 [](1),
                 [[2, -1]](1),
                 [[-2, -1]](1),
                 [[-1, -1]](1),
                 [[-1, -1]](2)]

            For != operator::

                sage: ([(i,j) for i in range(len(S)) for j in range(len(S)) if S[i]!=S[j]]
                ....: == [(i,j) for i in range(len(S)) for j in range(len(S)) if
                ....: S[i].value!=S[j].value])
                True

            For < operator::

                sage: ([(i,j) for i in range(len(S)) for j in range(len(S)) if S[i]<S[j]]
                ....: == [(i,j) for i in range(len(S)) for j in range(len(S)) if
                ....: S[i].value<S[j].value])
                True

            For <= operator::

                sage: ([(i,j) for i in range(len(S)) for j in range(len(S)) if S[i]<=S[j]]
                ....: == [(i,j) for i in range(len(S)) for j in range(len(S)) if
                ....: S[i].value<=S[j].value])
                True

            For > operator::

                sage: ([(i,j) for i in range(len(S)) for j in range(len(S)) if S[i]>S[j]]
                ....: == [(i,j) for i in range(len(S)) for j in range(len(S)) if
                ....: S[i].value>S[j].value])
                True

            For >= operator::

                sage: ([(i,j) for i in range(len(S)) for j in range(len(S)) if S[i]>=S[j]]
                ....: == [(i,j) for i in range(len(S)) for j in range(len(S)) if
                ....: S[i].value>=S[j].value])
                True
            """
            return richcmp(self.value, other.value, op)

        def e(self, i):
            """
            Return `e_i` of ``self``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.e(2)
                sage: mg.e(1)
                [[1, 4], [5]]
            """
            ret = self.value.e(i)
            if ret is None or not self.parent()._containing(ret):
                return None
            return self.__class__(self.parent(), ret)

        def f(self, i):
            """
            Return `f_i` of ``self``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.f(1)
                sage: mg.f(2)
                [[3, 4], [5]]
            """
            ret = self.value.f(i)
            if ret is None or not self.parent()._containing(ret):
                return None
            return self.__class__(self.parent(), ret)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.epsilon(1)
                1
                sage: mg.epsilon(2)
                0
            """
            return self.value.epsilon(i)

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.phi(1)
                0
                sage: mg.phi(2)
                1
            """
            return self.value.phi(i)

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.weight()
                (0, 1, 0, 1, 1)
            """
            return self.value.weight()
