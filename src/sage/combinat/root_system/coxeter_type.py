"""
Coxeter types
"""
# ****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>,
#                     2015 Jean-Philippe Labbe <labbe at math.huji.ac.il>,
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc

from sage.combinat.root_system.cartan_type import CartanType
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.matrix.args import SparseEntry
from sage.matrix.constructor import Matrix
from sage.misc.lazy_import import lazy_import
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject

lazy_import('sage.rings.universal_cyclotomic_field', 'UniversalCyclotomicField')


class CoxeterType(SageObject, metaclass=ClasscallMetaclass):
    """
    Abstract class for Coxeter types.
    """
    @staticmethod
    def __classcall_private__(cls, *x):
        """
        Parse input ``x``.

        EXAMPLES::

            sage: CoxeterType(['A',3])
            Coxeter type of ['A', 3]
        """
        if len(x) == 1:
            x = x[0]

        if isinstance(x, CoxeterType):
            return x

        try:
            return CoxeterTypeFromCartanType(CartanType(x))
        except (ValueError, TypeError):
            pass

        if len(x) == 1:  # In case the input is similar to CoxeterType([['A',2]])
            return CoxeterType(x[0])

        raise NotImplementedError("Coxeter types not from Cartan types not yet implemented")

    @classmethod
    def samples(self, finite=None, affine=None, crystallographic=None):
        """
        Return a sample of the available Coxeter types.

        INPUT:

        - ``finite`` -- boolean or ``None`` (default: ``None``)

        - ``affine`` -- boolean or ``None`` (default: ``None``)

        - ``crystallographic`` -- boolean or ``None`` (default: ``None``)

        The sample contains all the exceptional finite and affine
        Coxeter types, as well as typical representatives of the
        infinite families.

        EXAMPLES::

            sage: CoxeterType.samples()
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
             Coxeter type of ['B', 1], Coxeter type of ['B', 5],
             Coxeter type of ['C', 1], Coxeter type of ['C', 5],
             Coxeter type of ['D', 4], Coxeter type of ['D', 5],
             Coxeter type of ['E', 6], Coxeter type of ['E', 7],
             Coxeter type of ['E', 8], Coxeter type of ['F', 4],
             Coxeter type of ['H', 3], Coxeter type of ['H', 4],
             Coxeter type of ['I', 10], Coxeter type of ['A', 2, 1],
             Coxeter type of ['B', 5, 1], Coxeter type of ['C', 5, 1],
             Coxeter type of ['D', 5, 1], Coxeter type of ['E', 6, 1],
             Coxeter type of ['E', 7, 1], Coxeter type of ['E', 8, 1],
             Coxeter type of ['F', 4, 1], Coxeter type of ['G', 2, 1],
             Coxeter type of ['A', 1, 1]]

        The finite, affine and crystallographic options allow
        respectively for restricting to (non) finite, (non) affine,
        and (non) crystallographic Cartan types::

            sage: CoxeterType.samples(finite=True)
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
             Coxeter type of ['B', 1], Coxeter type of ['B', 5],
             Coxeter type of ['C', 1], Coxeter type of ['C', 5],
             Coxeter type of ['D', 4], Coxeter type of ['D', 5],
             Coxeter type of ['E', 6], Coxeter type of ['E', 7],
             Coxeter type of ['E', 8], Coxeter type of ['F', 4],
             Coxeter type of ['H', 3], Coxeter type of ['H', 4],
             Coxeter type of ['I', 10]]

            sage: CoxeterType.samples(affine=True)
            [Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
             Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
             Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
             Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
             Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1]]

            sage: CoxeterType.samples(crystallographic=True)
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
             Coxeter type of ['B', 1], Coxeter type of ['B', 5],
             Coxeter type of ['C', 1], Coxeter type of ['C', 5],
             Coxeter type of ['D', 4], Coxeter type of ['D', 5],
             Coxeter type of ['E', 6], Coxeter type of ['E', 7],
             Coxeter type of ['E', 8], Coxeter type of ['F', 4],
             Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
             Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
             Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
             Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
             Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1]]

            sage: CoxeterType.samples(crystallographic=False)
            [Coxeter type of ['H', 3],
             Coxeter type of ['H', 4],
             Coxeter type of ['I', 10]]

        .. TODO:: add some reducible Coxeter types (suggestions?)

        TESTS::

            sage: for ct in CoxeterType.samples(): TestSuite(ct).run()
        """
        result = self._samples()
        if crystallographic is not None:
            result = [t for t in result if t.is_crystallographic() == crystallographic]
        if finite is not None:
            result = [t for t in result if t.is_finite() == finite]
        if affine is not None:
            result = [t for t in result if t.is_affine() == affine]
        return result

    @cached_method
    def _samples(self):
        """
        Return a sample of all implemented Coxeter types.

        .. NOTE::

            This is intended to be used through :meth:`samples`.

        EXAMPLES::

            sage: CoxeterType._samples()
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
             Coxeter type of ['B', 1], Coxeter type of ['B', 5],
             Coxeter type of ['C', 1], Coxeter type of ['C', 5],
             Coxeter type of ['D', 4], Coxeter type of ['D', 5],
             Coxeter type of ['E', 6], Coxeter type of ['E', 7],
             Coxeter type of ['E', 8], Coxeter type of ['F', 4],
             Coxeter type of ['H', 3], Coxeter type of ['H', 4],
             Coxeter type of ['I', 10], Coxeter type of ['A', 2, 1],
             Coxeter type of ['B', 5, 1], Coxeter type of ['C', 5, 1],
             Coxeter type of ['D', 5, 1], Coxeter type of ['E', 6, 1],
             Coxeter type of ['E', 7, 1], Coxeter type of ['E', 8, 1],
             Coxeter type of ['F', 4, 1], Coxeter type of ['G', 2, 1],
             Coxeter type of ['A', 1, 1]]
        """
        finite = [CoxeterType(t) for t in [['A', 1], ['A', 5], ['B', 1], ['B', 5],
                                           ['C', 1], ['C', 5], ['D', 4], ['D', 5],
                                           ['E', 6], ['E', 7], ['E', 8], ['F', 4],
                                           ['H', 3], ['H', 4], ['I', 10]]]

        affine = [CoxeterType(t) for t in [['A', 2, 1], ['B', 5, 1],
                                           ['C', 5, 1], ['D', 5, 1], ['E', 6, 1],
                                           ['E', 7, 1], ['E', 8, 1], ['F', 4, 1],
                                           ['G', 2, 1], ['A', 1, 1]]]

        return finite + affine

    @abstract_method
    def rank(self):
        """
        Return the rank of ``self``.

        This is the number of nodes of the associated Coxeter graph.

        EXAMPLES::

            sage: CoxeterType(['A', 4]).rank()
            4
            sage: CoxeterType(['A', 7, 2]).rank()
            5
            sage: CoxeterType(['I', 8]).rank()
            2
        """

    @abstract_method
    def index_set(self):
        """
        Return the index set for ``self``.

        This is the list of the nodes of the associated Coxeter graph.

        EXAMPLES::

            sage: CoxeterType(['A', 3, 1]).index_set()
            (0, 1, 2, 3)
            sage: CoxeterType(['D', 4]).index_set()
            (1, 2, 3, 4)
            sage: CoxeterType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CoxeterType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CoxeterType(['A', 6, 2]).index_set()
            (0, 1, 2, 3)
            sage: CoxeterType(['D', 6, 2]).index_set()
            (0, 1, 2, 3, 4, 5)
            sage: CoxeterType(['E', 6, 1]).index_set()
            (0, 1, 2, 3, 4, 5, 6)
            sage: CoxeterType(['E', 6, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CoxeterType(['A', 2, 2]).index_set()
            (0, 1)
            sage: CoxeterType(['G', 2, 1]).index_set()
            (0, 1, 2)
            sage: CoxeterType(['F', 4, 1]).index_set()
            (0, 1, 2, 3, 4)
        """

    @abstract_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::

            sage: CoxeterType(['A', 3]).coxeter_matrix()                                # needs sage.graphs
            [1 3 2]
            [3 1 3]
            [2 3 1]
            sage: CoxeterType(['A', 3, 1]).coxeter_matrix()                             # needs sage.graphs
            [1 3 2 3]
            [3 1 3 2]
            [2 3 1 3]
            [3 2 3 1]
        """

    @abstract_method
    def coxeter_graph(self):
        """
        Return the Coxeter graph associated to ``self``.

        EXAMPLES::

            sage: CoxeterType(['A', 3]).coxeter_graph()                                 # needs sage.graphs
            Graph on 3 vertices
            sage: CoxeterType(['A', 3, 1]).coxeter_graph()                              # needs sage.graphs
            Graph on 4 vertices
        """

    @abstract_method
    def is_finite(self):
        """
        Return whether ``self`` is finite.

        EXAMPLES::

            sage: CoxeterType(['A',4]).is_finite()
            True
            sage: CoxeterType(['A',4, 1]).is_finite()
            False
        """

    @abstract_method
    def is_affine(self):
        """
        Return whether ``self`` is affine.

        EXAMPLES::

            sage: CoxeterType(['A', 3]).is_affine()
            False
            sage: CoxeterType(['A', 3, 1]).is_affine()
            True
        """

    def is_crystallographic(self):
        """
        Return whether ``self`` is crystallographic.

        This returns ``False`` by default. Derived class should override this
        appropriately.

        EXAMPLES::

            sage: [ [t, t.is_crystallographic() ] for t in CartanType.samples(finite=True) ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], True],
             [['C', 1], True], [['C', 5], True],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], True], [['G', 2], True],
             [['I', 5], False], [['H', 3], False], [['H', 4], False]]
        """
        return False

    def is_simply_laced(self):
        """
        Return whether ``self`` is simply laced.

        This returns ``False`` by default. Derived class should override this
        appropriately.

        EXAMPLES::

            sage: [ [t, t.is_simply_laced() ] for t in CartanType.samples() ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], False],
             [['C', 1], True], [['C', 5], False],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], False], [['G', 2], False],
             [['I', 5], False], [['H', 3], False], [['H', 4], False],
             [['A', 1, 1], False], [['A', 5, 1], True],
             [['B', 1, 1], False], [['B', 5, 1], False],
             [['C', 1, 1], False], [['C', 5, 1], False],
             [['D', 3, 1], True], [['D', 5, 1], True],
             [['E', 6, 1], True], [['E', 7, 1], True], [['E', 8, 1], True],
             [['F', 4, 1], False], [['G', 2, 1], False],
             [['BC', 1, 2], False], [['BC', 5, 2], False],
             [['B', 5, 1]^*, False], [['C', 4, 1]^*, False],
             [['F', 4, 1]^*, False], [['G', 2, 1]^*, False],
             [['BC', 1, 2]^*, False], [['BC', 5, 2]^*, False]]
        """
        return False

    @cached_method
    def bilinear_form(self, R=None):
        """
        Return the bilinear form over ``R`` associated to ``self``.

        INPUT:

        - ``R`` -- (default: universal cyclotomic field) a ring used to
          compute the bilinear form

        EXAMPLES::

            sage: # needs sage.graphs sage.libs.gap
            sage: CoxeterType(['A', 2, 1]).bilinear_form()
            [   1 -1/2 -1/2]
            [-1/2    1 -1/2]
            [-1/2 -1/2    1]
            sage: CoxeterType(['H', 3]).bilinear_form()
            [                      1                    -1/2                       0]
            [                   -1/2                       1 1/2*E(5)^2 + 1/2*E(5)^3]
            [                      0 1/2*E(5)^2 + 1/2*E(5)^3                       1]
            sage: C = CoxeterMatrix([[1,-1,-1],[-1,1,-1],[-1,-1,1]])
            sage: C.bilinear_form()
            [ 1 -1 -1]
            [-1  1 -1]
            [-1 -1  1]
        """

        n = self.rank()
        mat = self.coxeter_matrix()._matrix

        if R is None:
            R = UniversalCyclotomicField()

        # Compute the matrix with entries `- \cos( \pi / m_{ij} )`.
        if isinstance(R, sage.rings.abc.UniversalCyclotomicField):
            E = R.gen

            def val(x):
                if x > -1:
                    return (E(2*x) + ~E(2*x)) / R(-2)
                else:
                    return R(x)
        elif isinstance(R, sage.rings.abc.NumberField_quadratic):
            E = UniversalCyclotomicField().gen

            def val(x):
                if x > -1:
                    return R((E(2*x) + ~E(2*x)).to_cyclotomic_field()) / R(-2)
                else:
                    return R(x)
        else:
            from sage.functions.trig import cos
            from sage.symbolic.constants import pi
            from sage.symbolic.ring import SR

            def val(x):
                if x > -1:
                    return -R(cos(pi / SR(x)))
                else:
                    return R(x)

        entries = [SparseEntry(i, j, val(mat[i, j]))
                   for i in range(n) for j in range(n)
                   if mat[i, j] != 2]
        bilinear = Matrix(R, n, entries)
        bilinear.set_immutable()
        return bilinear


class CoxeterTypeFromCartanType(UniqueRepresentation, CoxeterType):
    """
    A Coxeter type associated to a Cartan type.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.combinat.root_system.coxeter_type import CoxeterTypeFromCartanType
            sage: C1 = CoxeterTypeFromCartanType(['A',3])
            sage: C2 = CoxeterTypeFromCartanType(CartanType(['A',3]))
            sage: C1 is C2
            True
        """
        return super().__classcall__(cls, CartanType(cartan_type))

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['A',3])
            sage: TestSuite(C).run()
            sage: C = CoxeterType(['H',4])
            sage: TestSuite(C).run()
            sage: C = CoxeterType(['C',3,1])
            sage: TestSuite(C).run()
        """
        self._cartan_type = cartan_type

    def _repr_(self, compact=False):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CoxeterType(['A',3])
            Coxeter type of ['A', 3]
        """
        if compact:
            return self._cartan_type._repr_(compact)
        return "Coxeter type of {}".format(self._cartan_type)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CoxeterType(['A',3]))
            A_{3}
        """
        return self._cartan_type._latex_()

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['H',3])
            sage: C.coxeter_matrix()                                                    # needs sage.graphs
            [1 3 2]
            [3 1 5]
            [2 5 1]
        """
        return self._cartan_type.coxeter_matrix()

    def coxeter_graph(self):
        """
        Return the Coxeter graph of ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['H',3])
            sage: C.coxeter_graph().edges(sort=True)                                    # needs sage.graphs
            [(1, 2, 3), (2, 3, 5)]
        """
        return self._cartan_type.coxeter_diagram()

    def cartan_type(self):
        """
        Return the Cartan type used to construct ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['C',3])
            sage: C.cartan_type()
            ['C', 3]
        """
        return self._cartan_type

    def type(self):
        """
        Return the type of ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['A', 4])
            sage: C.type()
            'A'
        """
        return self._cartan_type.type()

    def rank(self):
        """
        Return the rank of ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['I', 16])
            sage: C.rank()
            2
        """
        return self._cartan_type.rank()

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['A', 4])
            sage: C.index_set()
            (1, 2, 3, 4)
        """
        return self._cartan_type.index_set()

    def is_finite(self):
        """
        Return if ``self`` is a finite type.

        EXAMPLES::

            sage: C = CoxeterType(['E', 6])
            sage: C.is_finite()
            True
        """
        return self._cartan_type.is_finite()

    def is_affine(self):
        """
        Return if ``self`` is an affine type.

        EXAMPLES::

            sage: C = CoxeterType(['F', 4, 1])
            sage: C.is_affine()
            True
        """
        return self._cartan_type.is_affine()

    def is_crystallographic(self):
        """
        Return if ``self`` is crystallographic.

        EXAMPLES::

            sage: C = CoxeterType(['C', 3])
            sage: C.is_crystallographic()
            True

            sage: C = CoxeterType(['H', 3])
            sage: C.is_crystallographic()
            False
        """
        return self._cartan_type.is_crystallographic()

    def is_simply_laced(self):
        """
        Return if ``self`` is simply-laced.

        EXAMPLES::

            sage: C = CoxeterType(['A', 5])
            sage: C.is_simply_laced()
            True

            sage: C = CoxeterType(['B', 3])
            sage: C.is_simply_laced()
            False
        """
        return self._cartan_type.is_simply_laced()

    def is_reducible(self):
        """
        Return if ``self`` is reducible.

        EXAMPLES::

            sage: C = CoxeterType(['A', 5])
            sage: C.is_reducible()
            False

            sage: C = CoxeterType('A2xA2')
            sage: C.is_reducible()
            True
        """
        return self._cartan_type.is_reducible()

    def is_irreducible(self):
        """
        Return if ``self`` is irreducible.

        EXAMPLES::

            sage: C = CoxeterType(['A', 5])
            sage: C.is_irreducible()
            True

            sage: C = CoxeterType('B3xB3')
            sage: C.is_irreducible()
            False
        """
        return self._cartan_type.is_irreducible()

    def component_types(self):
        """
        A list of Coxeter types making up the reducible type.

        EXAMPLES::

            sage: CoxeterType(['A',2],['B',2]).component_types()
            [Coxeter type of ['A', 2], Coxeter type of ['B', 2]]

            sage: CoxeterType('A4xB3').component_types()
            [Coxeter type of ['A', 4], Coxeter type of ['B', 3]]

            sage: CoxeterType(['A', 2]).component_types()
            Traceback (most recent call last):
            ...
            ValueError: component types only defined for reducible types
        """
        if self.is_irreducible():
            raise ValueError('component types only defined for reducible types')
        return [CoxeterType(t) for t in self._cartan_type.component_types()]

    def relabel(self, relabelling):
        """
        Return a relabelled copy of ``self``.

        EXAMPLES::

            sage: ct = CoxeterType(['A',2])
            sage: ct.relabel({1:-1, 2:-2})
            Coxeter type of ['A', 2] relabelled by {1: -1, 2: -2}
        """
        return CoxeterType(self._cartan_type.relabel(relabelling))
