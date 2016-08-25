"""
Tensor Products of Crystals

Main entry points:

- :class:`~sage.combinat.crystals.tensor_product.TensorProductOfCrystals`
- :class:`~sage.combinat.crystals.tensor_product.CrystalOfTableaux`

AUTHORS:

- Anne Schilling, Nicolas Thiery (2007): Initial version
- Ben Salisbury, Travis Scrimshaw (2013): Refactored tensor products to handle
  non-regular crystals and created new subclass to take advantage of
  the regularity
"""
#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
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
from __future__ import print_function
from __future__ import absolute_import

import operator
from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.structure.parent import Parent
from sage.structure.element import parent
from sage.structure.global_options import GlobalOptions
from sage.categories.category import Category
from sage.categories.cartesian_product import cartesian_product
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.sets_cat import Sets
from sage.categories.map import Map
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.combinat import CombinatorialElement
from sage.combinat.partition import Partition
from sage.combinat.tableau import Tableau
from .letters import CrystalOfLetters
from .spins import CrystalOfSpins, CrystalOfSpinsMinus, CrystalOfSpinsPlus
from sage.misc.flatten import flatten
from sage.structure.element import get_coercion_model
from sage.rings.all import ZZ

##############################################################################
# Until trunc gets implemented in sage.function.other

from sage.functions.other import floor, ceil
def trunc(i):
    """
    Truncates to the integer closer to zero

    EXAMPLES::

        sage: from sage.combinat.crystals.tensor_product import trunc
        sage: trunc(-3/2), trunc(-1), trunc(-1/2), trunc(0), trunc(1/2), trunc(1), trunc(3/2)
        (-1, -1, 0, 0, 0, 1, 1)
        sage: isinstance(trunc(3/2), Integer)
        True
    """
    if i>= 0:
        return floor(i)
    else:
        return ceil(i)

##############################################################################
# Support classes
##############################################################################

from sage.structure.unique_representation import UniqueRepresentation

class TestParent(UniqueRepresentation, Parent):
    """
    A parent for tests.
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import TestParent
            sage: TestParent()
            A parent for tests
        """
        return "A parent for tests"

class ImmutableListWithParent(CombinatorialElement):
    r"""
    A class for lists having a parent

    Specification: any subclass ``C`` should implement ``__init__`` which
    accepts the following form ``C(parent, list = list)``

    EXAMPLES:

    We create an immutable list whose parent is the class list::

        sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
        sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
        sage: l._list
        [1, 2, 3]
        sage: l.parent()
        A parent for tests
        sage: l.sibling([2,1]) == ImmutableListWithParent(TestParent(), [2,1])
        True
        sage: l.reversed()
        [3, 2, 1]
        sage: l.set_index(1,4)
        [1, 4, 3]

    TESTS::

        sage: TestSuite(l).run(skip = "_test_category")
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: l._repr_()
            '[1, 2, 3]'
        """
        return repr(self._list)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = ImmutableListWithParent(ZZ, [1,2,3])
            sage: n = ImmutableListWithParent(ZZ, [2,3,4])
            sage: l == l
            True
            sage: l == m
            False
            sage: m == n
            False
        """
        return self.__class__ is other.__class__ and \
               self.parent() == other.parent() and \
               self._list == other._list

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = ImmutableListWithParent(ZZ, [1,2,3])
            sage: n = ImmutableListWithParent(ZZ, [2,3,4])
            sage: l != l
            False
            sage: l != m
            True
            sage: m != n
            True
        """
        return not self == other

    def __lt__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = ImmutableListWithParent(ZZ, [1,2,3])
            sage: n = ImmutableListWithParent(ZZ, [2,3,4])
            sage: l < l
            False
            sage: l < m
            False
            sage: m < n
            True
        """
        if parent(self) is not parent(other):
            return NotImplemented
        return self._list.__lt__(other._list)

    def __le__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = ImmutableListWithParent(ZZ, [1,2,3])
            sage: n = ImmutableListWithParent(ZZ, [2,3,4])
            sage: l <= l
            True
            sage: l <= m
            True
            sage: m <= n
            True
        """
        return self == other or self.__lt__(other)

    def __gt__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = ImmutableListWithParent(ZZ, [1,2,3])
            sage: n = ImmutableListWithParent(ZZ, [2,3,4])
            sage: l > l
            False
            sage: l > m
            False
            sage: m > n
            False
        """
        if parent(self) is not parent(other):
            return NotImplemented
        return other.__lt__(self)

    def __ge__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = ImmutableListWithParent(ZZ, [1,2,3])
            sage: n = ImmutableListWithParent(ZZ, [2,3,4])
            sage: l >= l
            True
            sage: l >= m
            True
            sage: m >= n
            False
        """
        return self == other or self.__gt__(other)

    def sibling(self, l):
        """
        Returns an :class:`ImmutableListWithParent` object whose list is
        ``l`` and whose parent is the same as the parent of ``self``.

        Note that the implementation of this function makes an assumption
        about the constructor for subclasses.

        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = l.sibling([2,3,4]); m
            [2, 3, 4]
            sage: m.parent()
            A parent for tests
        """
        return self.__class__(self.parent(), list=l)

    def reversed(self):
        """
        Returns the sibling of ``self`` which is obtained by reversing the
        elements of`` self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: l.reversed()
            [3, 2, 1]
        """
        return self.sibling([ i for i in reversed(self._list)])

    def set_index(self, k, value):
        """
        Returns the sibling of ``self`` obtained by setting the
        `k^{th}` entry of self to value.

        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: l.set_index(0,2)
            [2, 2, 3]
            sage: l.set_index(1,4)
            [1, 4, 3]
            sage: _.parent()
            A parent for tests
        """
        l = [i for i in self._list]
        l[k] = value
        return self.sibling(l)

class CrystalOfWords(UniqueRepresentation, Parent):
    """
    Auxiliary class to provide a call method to create tensor product elements.
    This class is shared with several tensor product classes and is also used
    in :class:`~sage.combinat.crystals.tensor_product.CrystalOfTableaux`
    to allow tableaux of different tensor product structures in
    column-reading (and hence different shapes) to be considered elements
    in the same crystal.
    """
    def _element_constructor_(self, *crystalElements):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(1,1)
            [1, 1]
            sage: _.parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
            sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: T(C(2), C(1), C(1))
            [2, 1, 1]
        """
        return self.element_class(self, list(crystalElements))

    def one_dimensional_configuration_sum(self, q=None, group_components=True):
        r"""
        Computes the one-dimensional configuration sum.

        INPUT:

        - ``q`` -- (default: ``None``) a variable or ``None``; if ``None``,
          a variable `q` is set in the code
        - ``group_components`` -- (default: ``True``) boolean; if ``True``,
          then the terms are grouped by classical component

        The one-dimensional configuration sum is the sum of the weights of all
        elements in the crystal weighted by the energy function.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
            sage: T = crystals.TensorProduct(K,K)
            sage: T.one_dimensional_configuration_sum()
            B[-2*Lambda[1] + 2*Lambda[2]] + (q+1)*B[-Lambda[1]] + (q+1)*B[Lambda[1] - Lambda[2]]
            + B[2*Lambda[1]] + B[-2*Lambda[2]] + (q+1)*B[Lambda[2]]
            sage: R.<t> = ZZ[]
            sage: T.one_dimensional_configuration_sum(t, False)
            B[-2*Lambda[1] + 2*Lambda[2]] + (t+1)*B[-Lambda[1]] + (t+1)*B[Lambda[1] - Lambda[2]]
            + B[2*Lambda[1]] + B[-2*Lambda[2]] + (t+1)*B[Lambda[2]]

            sage: R = RootSystem(['A',2,1])
            sage: La = R.weight_space().basis()
            sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
            sage: LS.one_dimensional_configuration_sum() == T.one_dimensional_configuration_sum() # long time
            True

        TESTS::

            sage: K1 = crystals.KirillovReshetikhin(['A',2,1],1,1)
            sage: K2 = crystals.KirillovReshetikhin(['A',2,1],2,1)
            sage: T = crystals.TensorProduct(K1,K2)
            sage: T.one_dimensional_configuration_sum() == T.one_dimensional_configuration_sum(group_components=False)
            True

            sage: RC = RiggedConfigurations(['A',3,1],[[1,1],[1,2]])
            sage: B = crystals.KirillovReshetikhin(['A',3,1],1,1)
            sage: B1 = crystals.KirillovReshetikhin(['A',3,1],1,2)
            sage: T = crystals.TensorProduct(B,B1)
            sage: RC.fermionic_formula() == T.one_dimensional_configuration_sum()
            True
        """
        if q is None:
            from sage.rings.all import QQ
            q = QQ['q'].gens()[0]
        P0 = self.weight_lattice_realization().classical()
        B = P0.algebra(q.parent())
        if group_components:
            G = self.digraph(index_set = self.cartan_type().classical().index_set())
            C = G.connected_components()
            return sum(q**(c[0].energy_function())*B.sum(B(P0(b.weight())) for b in c) for c in C)
        return B.sum(q**(b.energy_function())*B(P0(b.weight())) for b in self)


class TensorProductOfCrystals(CrystalOfWords):
    r"""
    Tensor product of crystals.

    Given two crystals `B` and `B'` of the same Cartan type,
    one can form the tensor product `B \otimes B^{\prime}`. As a set
    `B \otimes B^{\prime}` is the Cartesian product
    `B \times B^{\prime}`. The crystal operators `f_i` and
    `e_i` act on `b \otimes b^{\prime} \in B \otimes B^{\prime}` as
    follows:

    .. MATH::

        f_i(b \otimes b^{\prime}) = \begin{cases}
        f_i(b) \otimes b^{\prime} & \text{if } \varepsilon_i(b) \geq
        \varphi_i(b^{\prime}) \\
        b \otimes f_i(b^{\prime}) & \text{otherwise}
        \end{cases}

    and

    .. MATH::

        e_i(b \otimes b') = \begin{cases}
        e_i(b) \otimes b' & \text{if } \varepsilon_i(b) >
        \varphi_i(b') \\ b \otimes e_i(b') & \text{otherwise.}
        \end{cases}

    We also define:

    .. MATH::

        \begin{aligned}
        \varphi_i(b \otimes b') & = \max\left( \varphi_i(b),
        \varphi_i(b') + \langle \alpha_i^{\vee}, \mathrm{wt}(b) \rangle
        \right),
        \\ \varepsilon_i(b \otimes b') & = \max\left( \varepsilon_i(b'),
        \varepsilon_i(b) - \langle \alpha_i^{\vee}, \mathrm{wt}(b') \rangle
        \right).
        \end{aligned}

    .. NOTE::

        This is the opposite of Kashiwara's convention for tensor
        products of crystals.

    Since tensor products are associative `(\mathcal{B} \otimes \mathcal{C})
    \otimes \mathcal{D} \cong \mathcal{B} \otimes (\mathcal{C} \otimes
    \mathcal{D})` via the natural isomorphism `(b \otimes c) \otimes d \mapsto
    b \otimes (c \otimes d)`, we can generalizing this to arbitrary tensor
    products. Thus consider `B_N \otimes \cdots \otimes B_1`, where each
    `B_k` is an abstract crystal. The underlying set of the tensor product is
    `B_N \times \cdots \times B_1`, while the crystal structure is given
    as follows. Let `I` be the index set, and fix some `i \in I` and `b_N
    \otimes \cdots \otimes b_1 \in B_N \otimes \cdots \otimes B_1`. Define

    .. MATH::

        a_i(k) := \varepsilon_i(b_k) - \sum_{j=1}^{k-1} \langle
        \alpha_i^{\vee}, \mathrm{wt}(b_j) \rangle.

    Then

    .. MATH::

        \begin{aligned}
        \mathrm{wt}(b_N \otimes \cdots \otimes b_1) &=
        \mathrm{wt}(b_N) + \cdots + \mathrm{wt}(b_1),
        \\ \varepsilon_i(b_N \otimes \cdots \otimes b_1) &= \max_{1 \leq k
        \leq n}\left( \sum_{j=1}^k \varepsilon_i(b_j) - \sum_{j=1}^{k-1}
        \varphi_i(b_j) \right)
        \\ & = \max_{1 \leq k \leq N}\bigl( a_i(k) \bigr),
        \\ \varphi_i(b_N \otimes \cdots \otimes b_1) &= \max_{1 \leq k \leq N}
        \left( \varphi_i(b_N) + \sum_{j=k}^{N-1} \big( \varphi_i(b_j) -
        \varepsilon_i(b_{j+1}) \big) \right)
        \\ & = \max_{1 \leq k \leq N}\bigl( \lambda_i + a_i(k) \bigr)
        \end{aligned}

    where `\lambda_i = \langle \alpha_i^{\vee}, \mathrm{wt}(b_N \otimes \cdots
    \otimes b_1) \rangle`. Then for `k = 1, \ldots, N` the action of the
    Kashiwara operators is determined as follows.

    - If `a_i(k) > a_i(j)` for `1 \leq j < k` and `a_i(k) \geq a_i(j)`
      for `k < j \leq N`:

      .. MATH::

          e_i(b_N \otimes \cdots \otimes b_1) = b_N \otimes \cdots \otimes
          e_i b_k \otimes \cdots \otimes b_1.

    - If `a_i(k) \geq a_i(j)` for `1 \leq j < k` and `a_i(k) > a_i(j)`
      for `k < j \leq N`:

      .. MATH::

          f_i(b_N \otimes \cdots \otimes b_1) = b_N \otimes \cdots \otimes
          f_i b_k \otimes \cdots \otimes b_1.

    Note that this is just recursively applying the definition of the tensor
    product on two crystals. Recall that `\langle \alpha_i^{\vee},
    \mathrm{wt}(b_j) \rangle = \varphi_i(b_j) - \varepsilon_i(b_j)` by the
    definition of a crystal.

    .. RUBRIC:: Regular crystals

    Now if all crystals `B_k` are regular crystals, all `\varepsilon_i` and
    `\varphi_i` are non-negative and we can
    define tensor product by the *signature rule*. We start by writing a word
    in `+` and `-` as follows:

    .. MATH::

        \underbrace{- \cdots -}_{\varphi_i(b_N) \text{ times}} \quad
        \underbrace{+ \cdots +}_{\varepsilon_i(b_N) \text{ times}}
        \quad \cdots \quad
        \underbrace{- \cdots -}_{\varphi_i(b_1) \text{ times}} \quad
        \underbrace{+ \cdots +}_{\varepsilon_i(b_1) \text{ times}},

    and then canceling ordered pairs of `+-` until the word is in the reduced
    form:

    .. MATH::

        \underbrace{- \cdots -}_{\varphi_i \text{ times}} \quad
        \underbrace{+ \cdots +}_{\varepsilon_i \text{ times}}.

    Here `e_i` acts on the factor corresponding to the leftmost `+` and `f_i`
    on the factor corresponding to the rightmost `-`. If there is no `+` or
    `-` respectively, then the result is `0` (``None``).

    EXAMPLES:

    We construct the type `A_2`-crystal generated by `2 \otimes 1 \otimes 1`::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])

    It has `8` elements::

        sage: T.list()
        [[2, 1, 1], [2, 1, 2], [2, 1, 3], [3, 1, 3], [3, 2, 3], [3, 1, 1], [3, 1, 2], [3, 2, 2]]

    One can also check the Cartan type of the crystal::

        sage: T.cartan_type()
        ['A', 2]

    Other examples include crystals of tableaux (which internally are
    represented as tensor products obtained by reading the tableaux
    columnwise)::

        sage: C = crystals.Tableaux(['A',3], shape=[1,1,0])
        sage: D = crystals.Tableaux(['A',3], shape=[1,0,0])
        sage: T = crystals.TensorProduct(C,D, generators=[[C(rows=[[1], [2]]), D(rows=[[1]])], [C(rows=[[2], [3]]), D(rows=[[1]])]])
        sage: T.cardinality()
        24
        sage: TestSuite(T).run()
        sage: T.module_generators
        [[[[1], [2]], [[1]]], [[[2], [3]], [[1]]]]
        sage: [x.weight() for x in T.module_generators]
        [(2, 1, 0, 0), (1, 1, 1, 0)]

    If no module generators are specified, we obtain the full tensor
    product::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C)
        sage: T.list()
        [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
        sage: T.cardinality()
        9

    For a tensor product of crystals without module generators, the
    default implementation of ``module_generators`` contains all elements
    in the tensor product of the crystals. If there is a subset of
    elements in the tensor product that still generates the crystal,
    this needs to be implemented for the specific crystal separately::

        sage: T.module_generators.list()
        [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    For classical highest weight crystals, it is also possible to list
    all highest weight elements::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)],[C(1),C(2),C(1)]])
        sage: T.highest_weight_vectors()
        ([2, 1, 1], [1, 2, 1])

    Examples with non-regular and infinite crystals (these did not work
    before :trac:`14402`)::

        sage: B = crystals.infinity.Tableaux(['D',10])
        sage: T = crystals.TensorProduct(B,B)
        sage: T
        Full tensor product of the crystals
        [The infinity crystal of tableaux of type ['D', 10],
         The infinity crystal of tableaux of type ['D', 10]]

        sage: B = crystals.infinity.GeneralizedYoungWalls(15)
        sage: T = crystals.TensorProduct(B,B,B)
        sage: T
        Full tensor product of the crystals
        [Crystal of generalized Young walls of type ['A', 15, 1],
         Crystal of generalized Young walls of type ['A', 15, 1],
         Crystal of generalized Young walls of type ['A', 15, 1]]

        sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
        sage: B = crystals.GeneralizedYoungWalls(2,La[0]+La[1])
        sage: C = crystals.GeneralizedYoungWalls(2,2*La[2])
        sage: D = crystals.GeneralizedYoungWalls(2,3*La[0]+La[2])
        sage: T = crystals.TensorProduct(B,C,D)
        sage: T
        Full tensor product of the crystals
        [Highest weight crystal of generalized Young walls of Cartan type ['A', 2, 1] and highest weight Lambda[0] + Lambda[1],
         Highest weight crystal of generalized Young walls of Cartan type ['A', 2, 1] and highest weight 2*Lambda[2],
         Highest weight crystal of generalized Young walls of Cartan type ['A', 2, 1] and highest weight 3*Lambda[0] + Lambda[2]]

    There is also a global option for setting the convention (by default Sage
    uses anti-Kashiwara)::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C)
        sage: elt = T(C(1), C(2)); elt
        [1, 2]
        sage: crystals.TensorProduct.options.convention = "Kashiwara"
        sage: elt
        [2, 1]
        sage: crystals.TensorProduct.options._reset()
    """
    @staticmethod
    def __classcall_private__(cls, *crystals, **options):
        """
        Create the correct parent object.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C, C)
            sage: T2 = crystals.TensorProduct(C, C, cartan_type=['A',2])
            sage: T is T2
            True
            sage: T.category()
            Category of tensor products of classical crystals

            sage: T3 = crystals.TensorProduct(C, C, C)
            sage: T3p = crystals.TensorProduct(T, C)
            sage: T3 is T3p
            True
            sage: B1 = crystals.TensorProduct(T, C)
            sage: B2 = crystals.TensorProduct(C, T)
            sage: B3 = crystals.TensorProduct(C, C, C)
            sage: B1 is B2 and B2 is B3
            True

            sage: B = crystals.infinity.Tableaux(['A',2])
            sage: T = crystals.TensorProduct(B, B)
            sage: T.category()
            Category of infinite tensor products of highest weight crystals

        TESTS:

        Check that mismatched Cartan types raise an error::

            sage: A2 = crystals.Letters(['A', 2])
            sage: A3 = crystals.Letters(['A', 3])
            sage: crystals.TensorProduct(A2, A3)
            Traceback (most recent call last):
            ...
            ValueError: all crystals must be of the same Cartan type
        """
        crystals = tuple(crystals)
        if "cartan_type" in options:
            cartan_type = CartanType(options.pop("cartan_type"))
        else:
            if not crystals:
                raise ValueError("you need to specify the Cartan type if the tensor product list is empty")
            else:
                cartan_type = crystals[0].cartan_type()

        if any(c.cartan_type() != cartan_type for c in crystals):
            raise ValueError("all crystals must be of the same Cartan type")

        if "generators" in options:
            generators = tuple(tuple(x) if isinstance(x, list) else x for x in options["generators"])

            if all(c in RegularCrystals() for c in crystals):
                return TensorProductOfRegularCrystalsWithGenerators(crystals, generators, cartan_type)
            return TensorProductOfCrystalsWithGenerators(crystals, generators, cartan_type)

        # Flatten out tensor products
        tp = sum([B.crystals if isinstance(B, FullTensorProductOfCrystals) else (B,)
                  for B in crystals], ())

        if all(c in RegularCrystals() for c in crystals):
            return FullTensorProductOfRegularCrystals(tp, cartan_type=cartan_type)
        return FullTensorProductOfCrystals(tp, cartan_type=cartan_type)

    # add options to class
    options=GlobalOptions('TensorProductOfCrystals', 
        module='sage.combinat.crystals',
        doc=r"""
        Sets the global options for tensor products of crystals. The default is to
        use the anti-Kashiwara convention.

        There are two conventions for how `e_i` and `f_i` act on tensor products,
        and the difference between the two is the order of the tensor factors
        are reversed. This affects both the input and output. See the example
        below.
        """,
        end_doc=r"""

        .. NOTE::

            Changing the ``convention`` also changes how the input is handled.

        .. WARNING::

            Internally, the crystals are always stored using the anti-Kashiwara
            convention.

        If no parameters are set, then the function returns a copy of the
        options dictionary.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: elt = T(C(1), C(2)); elt
            [1, 2]
            sage: crystals.TensorProduct.options.convention = "Kashiwara"
            sage: elt
            [2, 1]
            sage: T(C(1), C(2)) == elt
            False
            sage: T(C(2), C(1)) == elt
            True
            sage: crystals.TensorProduct.options._reset()
        """,
        convention=dict(default="antiKashiwara",
                        description='Sets the convention used for displaying/inputting tensor product of crystals',
                        values=dict(antiKashiwara='use the anti-Kashiwara convention',
                                    Kashiwara='use the Kashiwara convention'),
                            alias=dict(anti="antiKashiwara", opposite="antiKashiwara"),
                            case_sensitive=False)
    )

    def _element_constructor_(self, *crystalElements):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(1,1)
            [1, 1]
            sage: _.parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
            sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: T(C(2), C(1), C(1))
            [2, 1, 1]
        """
        if self.options.convention == "Kashiwara":
            crystalElements = reversed(crystalElements)
        return self.element_class(self, list(crystalElements))

class TensorProductOfCrystalsWithGenerators(TensorProductOfCrystals):
    """
    Tensor product of crystals with a generating set.

    .. TODO::

        Deprecate this class in favor of using
        :meth:`~sage.categories.crystals.Crystals.ParentMethods.subcrystal`.
    """
    def __init__(self, crystals, generators, cartan_type):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: TestSuite(T).run()
        """
        assert isinstance(crystals, tuple)
        assert isinstance(generators, tuple)
        category = Category.meet([crystal.category() for crystal in crystals])
        Parent.__init__(self, category = category)
        self.crystals = crystals
        self._cartan_type = cartan_type
        self.module_generators = [ self(*x) for x in generators ]

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: crystals.TensorProduct(C,C,generators=[[C(2),C(1)]])
            The tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
        """
        if self.options.convention == "Kashiwara":
            st = repr(list(reversed(self.crystals)))
        else:
            st = repr(list(self.crystals))
        return "The tensor product of the crystals {}".format(st)

class FullTensorProductOfCrystals(TensorProductOfCrystals):
    """
    Full tensor product of crystals.

    .. TODO::

        Merge this into :class:`TensorProductOfCrystals`.
    """
    def __init__(self, crystals, **options):
        """
        TESTS::

            sage: from sage.combinat.crystals.tensor_product import FullTensorProductOfCrystals
            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: isinstance(T, FullTensorProductOfCrystals)
            True
            sage: TestSuite(T).run()
        """
        category = Category.meet([crystal.category() for crystal in crystals])
        category = category.TensorProducts()
        if any(c in Sets().Infinite() for c in crystals):
            category = category.Infinite()
        Parent.__init__(self, category=category)
        self.crystals = crystals
        if 'cartan_type' in options:
            self._cartan_type = CartanType(options['cartan_type'])
        else:
            if not crystals:
                raise ValueError("you need to specify the Cartan type if the tensor product list is empty")
            else:
                self._cartan_type = crystals[0].cartan_type()
        self.cartesian_product = cartesian_product(self.crystals)
        self.module_generators = self

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: crystals.TensorProduct(C,C)
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
        """
        if self.options.convention == "Kashiwara":
            st = repr(list(reversed(self.crystals)))
        else:
            st = repr(list(self.crystals))
        return "Full tensor product of the crystals {}".format(st)

    # TODO: __iter__ and cardinality should be inherited from EnumeratedSets().CartesianProducts()
    def __iter__(self):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: list(T)
            [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
            sage: _[0].parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
        """
        for x in self.cartesian_product:
            yield self(*x)

#    list = CombinatorialClass._CombinatorialClass__list_from_iterator

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: T.cardinality()
            9
        """
        return self.cartesian_product.cardinality()

    @cached_method
    def weight_lattice_realization(self):
        r"""
        Return the weight lattice realization used to express weights.

        The weight lattice realization is the common parent which all
        weight lattice realizations of the crystals of ``self`` coerce
        into.

        EXAMPLES::

            sage: B = crystals.elementary.B(['A',4], 2)
            sage: B.weight_lattice_realization()
            Root lattice of the Root system of type ['A', 4]
            sage: T = crystals.infinity.Tableaux(['A',4])
            sage: T.weight_lattice_realization()
            Ambient space of the Root system of type ['A', 4]
            sage: TP = crystals.TensorProduct(B, T)
            sage: TP.weight_lattice_realization()
            Ambient space of the Root system of type ['A', 4]
        """
        cm = get_coercion_model()
        return cm.common_parent(*[crystal.weight_lattice_realization()
                                  for crystal in self.crystals])

class TensorProductOfCrystalsElement(ImmutableListWithParent):
    r"""
    A class for elements of tensor products of crystals.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',3])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(2))
            [1, 2]
        """
        if self.parent().options.convention == "Kashiwara":
            return repr(list(reversed(self._list)))
        return repr(self._list)

    def _latex_(self):
        r"""
        Return latex code for ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(["A",2])
            sage: D = crystals.Tableaux(["A",2], shape=[2])
            sage: E = crystals.TensorProduct(C,D)
            sage: latex(E.module_generators[0])
            1 \otimes {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{1}&\lr{1}\\\cline{1-2}
            \end{array}$}
            }
        """
        return ' \otimes '.join(latex(c) for c in self)

    def _ascii_art_(self):
        """
        Return an ASCII art representation of ``self``.

        EXAMPLES::

            sage: KT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D',4,1],[[3,3],[2,1],[1,2]])
            sage: ascii_art(KT.module_generators[0])
              1  1  1
              2  2  2 #   1 #   1  1
              3  3  3     2
             -4 -4 -4
        """
        from sage.typeset.ascii_art import ascii_art, AsciiArt
        s = ascii_art(self[0])
        s._baseline = s._h // 2
        ret = s
        for tableau in self[1:]:
            s = ascii_art(tableau)
            s._baseline = s._h // 2
            ret += AsciiArt([" # "]) + s
        return ret

    def __lt__(self, other):
        """
        Non elements of the crystal are incomparable with elements of the crystal
        (or should it return NotImplemented?).

        Comparison of two elements of this crystal:

        - different length: incomparable
        - otherwise lexicographicaly, considering ``self[i]`` and ``other[i]``
          as incomparable if ``self[i] < other[i]`` returns NotImplemented
        """
        if parent(self) is not parent(other):
            return False
        if len(self) != len(other):
            return False
        for i in range(len(self)):
            if (self[i] < other[i]):
                return True
            if (other[i] < self[i]):
                return False
        return False

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',3], shape=[3,1])
            sage: D = crystals.Tableaux(['A',3], shape=[1])
            sage: E = crystals.Tableaux(['A',3], shape=[2,2,2])
            sage: T = crystals.TensorProduct(C,D,E)
            sage: print(T.module_generators[0]._repr_diagram())
              1  1  1 (X)   1 (X)   1  1
              2                     2  2
                                    3  3
        """
        pplist = []
        max_widths = []
        num_cols = len(self)
        for c in self:
            try:
                pplist.append(c._repr_diagram().split('\n'))
            except AttributeError:
                pplist.append(c._repr_().split('\n'))
            max_widths.append(max(map(len, pplist[-1])))
        num_rows = max(map(len, pplist))
        ret = ""
        for i in range(num_rows):
            if i > 0:
                ret += '\n'
            for j in range(num_cols):
                if j > 0:
                    if i == 0:
                        ret += ' (X) '
                    else:
                        ret += '     '
                if i < len(pplist[j]):
                    ret += pplist[j][i]
                    ret += ' '*(max_widths[j] - len(pplist[j][i]))
                else:
                    ret += ' '*max_widths[j]
        return ret

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',3], shape=[3,1])
            sage: D = crystals.Tableaux(['A',3], shape=[1])
            sage: E = crystals.Tableaux(['A',3], shape=[2,2,2])
            sage: T = crystals.TensorProduct(C,D,E)
            sage: T.module_generators[0].pp()
              1  1  1 (X)   1 (X)   1  1
              2                     2  2
                                    3  3
        """
        print(self._repr_diagram())

    def weight(self):
        r"""
        Return the weight of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux("A3")
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([2,1,3])
            sage: b2 = B.highest_weight_vector().f(1)
            sage: t = T(b2, b1)
            sage: t
            [[[1, 1, 1, 2], [2, 2], [3]], [[1, 1, 1, 1, 2], [2, 2, 4], [3]]]
            sage: t.weight()
            (-2, 1, 0, 1)
        """
        return sum(self[i].weight() for i in range(len(self)))

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.

        INPUT:

        - ``i`` -- An element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux("G2")
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f(2)
            sage: b2 = B.highest_weight_vector().f_string([2,2,1])
            sage: t = T(b2, b1)
            sage: [t.epsilon(i) for i in B.index_set()]
            [0, 3]
        """
        return max(self._sig(i, k) for k in range(1, len(self)+1))

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.

        INPUT:

        - ``i`` -- An element of the index set

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: B = crystals.GeneralizedYoungWalls(2,La[0]+La[1])
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([1,0])
            sage: b2 = B.highest_weight_vector().f_string([0,1])
            sage: t = T(b2, b1)
            sage: [t.phi(i) for i in B.index_set()]
            [1, 1, 4]

        TESTS:

        Check that :trac:`15462` is fixed::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: La = RootSystem(['A',2]).ambient_space().fundamental_weights()
            sage: T = crystals.TensorProduct(crystals.elementary.T(['A',2], La[1]+La[2]), B)
            sage: t = T.an_element()
            sage: t.phi(1)
            2
            sage: t.phi(2)
            2
        """
        P = self[-1].parent().weight_lattice_realization()
        h = P.simple_coroots()
        omega = P(self.weight()).scalar(h[i])
        return max([omega + self._sig(i, k) for k in range(1, len(self)+1)])

    @cached_in_parent_method
    def _sig(self,i,k):
        r"""
        Return `a_i(k)` of ``self``.

        The value `a_i(k)` of a crystal `b = b_N \otimes \cdots \otimes b_1`
        is defined as:

        .. MATH::

            a_i(k) = \varepsilon_i(b_k) - \sum_{j=1}^{k-1} \langle h_i,
            \mathrm{wt}(b_j) \rangle

        where `\mathrm{wt}` is the :meth:`weight` of `b_j`.

        INPUT:

        - ``i`` -- An element of the index set

        - ``k`` -- The (1-based) index of the tensor factor of ``self``

        EXAMPLES::

            sage: B = crystals.infinity.GeneralizedYoungWalls(3)
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([0,3,1])
            sage: b2 = B.highest_weight_vector().f_string([3,2,1,0,2,3])
            sage: t = T(b1, b2)
            sage: [[t._sig(i,k) for k in range(1,len(t)+1)] for i in B.index_set()]
            [[0, -1], [0, 0], [0, 1], [1, 2]]

        TESTS:

        Check that :trac:`18469` is fixed::

            sage: E1 = crystals.elementary.B(['A',2], 1)
            sage: E2 = crystals.elementary.B(['A',2], 2)
            sage: T = crystals.TensorProduct(E1, E2)
            sage: x = T(E1.module_generators[0], E2.module_generators[0]); x
            [0, 0]
            sage: [[x._sig(i,k) for k in range(1,3)] for i in T.index_set()]
            [[-inf, 0], [0, -inf]]
            sage: x.f(1)
            [-1, 0]
            sage: x.e(1)
            [1, 0]
        """
        if k == 1:
            return self[-1].epsilon(i)
        ep = self[-k].epsilon(i)
        if ep == float("-inf"):
            return ep

        P = self[-1].parent().weight_lattice_realization()
        h = P.simple_coroots()
        wt = sum(P(self[-j].weight()) for j in range(1, k))
        return ep - P(wt).scalar(h[i])

    def e(self,i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- An element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux("D4")
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([1,4,3])
            sage: b2 = B.highest_weight_vector().f_string([2,2,3,1,4])
            sage: t = T(b2, b1)
            sage: t.e(1)
            [[[1, 1, 1, 1, 1], [2, 2, 3, -3], [3]], [[1, 1, 1, 1, 2], [2, 2, 2], [3, -3]]]
            sage: t.e(2)
            sage: t.e(3)
            [[[1, 1, 1, 1, 1, 2], [2, 2, 3, -4], [3]], [[1, 1, 1, 1, 2], [2, 2, 2], [3, -3]]]
            sage: t.e(4)
            [[[1, 1, 1, 1, 1, 2], [2, 2, 3, 4], [3]], [[1, 1, 1, 1, 2], [2, 2, 2], [3, -3]]]
        """
        N = len(self) + 1
        for k in range(1, N):
            if all(self._sig(i,k) > self._sig(i,j) for j in range(1, k)) and \
                    all(self._sig(i,k) >= self._sig(i,j) for j in range(k+1, N)):
                crystal = self[-k].e(i)
                if crystal is None:
                    return None
                return self.set_index(-k, crystal)
        return None

    def f(self,i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- An element of the index set

        EXAMPLES::

            sage: La = RootSystem(['A',3,1]).weight_lattice(extended=True).fundamental_weights()
            sage: B = crystals.GeneralizedYoungWalls(3,La[0])
            sage: T = crystals.TensorProduct(B,B,B)
            sage: b1 = B.highest_weight_vector().f_string([0,3])
            sage: b2 = B.highest_weight_vector().f_string([0])
            sage: b3 = B.highest_weight_vector()
            sage: t = T(b3, b2, b1)
            sage: t.f(0)
            [[[0]], [[0]], [[0, 3]]]
            sage: t.f(1)
            [[], [[0]], [[0, 3], [1]]]
            sage: t.f(2)
            [[], [[0]], [[0, 3, 2]]]
            sage: t.f(3)
            [[], [[0, 3]], [[0, 3]]]
        """
        N = len(self) + 1
        for k in range(1, N):
            if all(self._sig(i,k) >= self._sig(i,j) for j in range(1, k)) and \
                    all(self._sig(i,k) > self._sig(i,j) for j in range(k+1, N)):
                crystal = self[-k].f(i)
                if crystal is None:
                    return None
                return self.set_index(-k, crystal)
        return None

class TensorProductOfRegularCrystalsElement(TensorProductOfCrystalsElement):
    """
    Element class for a tensor product of regular crystals.

    TESTS::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C, C)
        sage: elt = T(C(1), C(2))
        sage: from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
        sage: isinstance(elt, TensorProductOfRegularCrystalsElement)
        True
    """
    def e(self, i):
        """
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(2)).e(1) == T(C(1),C(1))
            True
            sage: T(C(2),C(1)).e(1) is None
            True
            sage: T(C(2),C(2)).e(1) == T(C(1),C(2))
            True
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        position = self.positions_of_unmatched_plus(i)
        if position == []:
            return None
        k = position[0]
        return self.set_index(k, self[k].e(i))

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',3])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(2)).weight()
            (1, 1, 0, 0)
            sage: T = crystals.Tableaux(['D',4],shape=[])
            sage: T.list()[0].weight()
            (0, 0, 0, 0)
        """
        return sum((self[j].weight() for j in range(len(self))), self.parent().weight_lattice_realization().zero())

    def f(self, i):
        """
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(1)).f(1)
            [1, 2]
            sage: T(C(1),C(2)).f(1)
            [2, 2]
            sage: T(C(2),C(1)).f(1) is None
            True
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        position = self.positions_of_unmatched_minus(i)
        if position == []:
            return None
        k = position[len(position)-1]
        return self.set_index(k, self[k].f(i))

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(1)).phi(1)
            2
            sage: T(C(1),C(2)).phi(1)
            1
            sage: T(C(2),C(1)).phi(1)
            0
        """
        self = self.reversed()
        height = 0
        for j in range(len(self)):
            plus = self[j].epsilon(i)
            minus = self[j].phi(i)
            if height-plus < 0:
                height = minus
            else:
                height = height - plus + minus
        return height

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(1)).epsilon(1)
            0
            sage: T(C(1),C(2)).epsilon(1)
            1
            sage: T(C(2),C(1)).epsilon(1)
            0
        """
        height = 0
        for j in range(len(self)):
            minus = self[j].phi(i)
            plus = self[j].epsilon(i)
            if height-minus < 0:
                height = plus
            else:
                height = height - minus + plus
        return height

    def positions_of_unmatched_minus(self, i, dual=False, reverse=False):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(2),C(1)).positions_of_unmatched_minus(1)
            []
            sage: T(C(1),C(2)).positions_of_unmatched_minus(1)
            [0]
        """
        unmatched_plus = []
        height = 0
        if reverse:
            self = self.reversed()
        if not dual:
            for j in range(len(self)):
                minus = self[j].phi(i)
                plus = self[j].epsilon(i)
                if height-minus < 0:
                    unmatched_plus.append(j)
                    height = plus
                else:
                    height = height - minus + plus
        else:
            for j in range(len(self)):
                plus = self[j].epsilon(i)
                minus = self[j].phi(i)
                if height-plus < 0:
                    unmatched_plus.append(j)
                    height = minus
                else:
                    height = height - plus + minus
        return unmatched_plus

    def positions_of_unmatched_plus(self, i):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(2),C(1)).positions_of_unmatched_plus(1)
            []
            sage: T(C(1),C(2)).positions_of_unmatched_plus(1)
            [1]
        """
        l = self.positions_of_unmatched_minus(i, dual=True, reverse=True)
        l.reverse()
        return [len(self)-1-l[j] for j in range(len(l))]

    def energy_function(self, algorithm=None):
        r"""
        Return the energy function of ``self``.

        ALGORITHM:

        .. RUBRIC:: definition

        Let `T` be a tensor product of Kirillov-Reshetikhin
        crystals. Let `R_i` and `H_i` be the combinatorial
        `R`-matrix and local energy functions, respectively, acting
        on the `i` and `i+1` factors. Let `D_B` be the energy
        function of a single Kirillov-Reshetikhin crystal. The
        *energy function* is given by

        .. MATH::

            D = \sum_{j > i} H_i R_{i+1} R_{i+2} \cdots R_{j-1}
            + \sum_j D_B R_1 R_2 \cdots R_{j-1},

        where `D_B` acts on the rightmost factor.

        .. RUBRIC:: grading

        If  ``self`` is an element of `T`, a tensor product of
        perfect crystals of the same level, then use the affine
        grading to determine the energy. Specifically, let `g`
        denote the affine grading of ``self`` and `d` the affine
        grading of the maximal vector in `T`. Then the energy
        of ``self`` is given by `d - g`.

        For more details, see Theorem 7.5 in [SchillingTingley2011]_.

        INPUT:

        - ``algorithm`` -- (default: ``None``) use one of the
          following algorithms to determine the energy function:

          * ``'definition'`` - use the definition of the energy
            function;
          * ``'grading'`` - use the affine grading;

          if not specified, then this uses ``'grading'`` if all
          factors are perfect of the same level and otherwise
          this uses ``'definition'``

        OUTPUT: an integer

        REFERENCES:

        .. [SchillingTingley2011] \A. Schilling, P. Tingley.
           *Demazure crystals, Kirillov-Reshetikhin crystals, and
           the energy function*.
           Electronic Journal of Combinatorics. **19(2)**. 2012.
           :arXiv:`1104.2359`

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
            sage: T = crystals.TensorProduct(K,K,K)
            sage: hw = sorted([x for x in T if x.is_highest_weight([1,2])])
            sage: for b in hw:
            ....:     print("{} {}".format(b, b.energy_function()))
            [[[1]], [[1]], [[1]]] 0
            [[[1]], [[2]], [[1]]] 2
            [[[2]], [[1]], [[1]]] 1
            [[[3]], [[2]], [[1]]] 3

            sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
            sage: T = crystals.TensorProduct(K,K)
            sage: hw = [x for x in T if x.is_highest_weight([1,2])]
            sage: for b in hw:
            ....:     print("{} {}".format(b, b.energy_function()))
            [[], []] 4
            [[], [[1, 1]]] 1
            [[[1, 1]], []] 3
            [[[1, 1]], [[1, 1]]] 0
            [[[1, 2]], [[1, 1]]] 1
            [[[2, 2]], [[1, 1]]] 2
            [[[-1, -1]], [[1, 1]]] 2
            [[[1, -1]], [[1, 1]]] 2
            [[[2, -1]], [[1, 1]]] 2

            sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 1)
            sage: T = crystals.TensorProduct(K)
            sage: t = T.module_generators[0]
            sage: t.energy_function('grading')
            Traceback (most recent call last):
            ...
            NotImplementedError: all crystals in the tensor product need to be perfect of the same level

        TESTS::

            sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
            sage: K2 = crystals.KirillovReshetikhin(['C',2,1], 2, 2)
            sage: T = tensor([K, K2])
            sage: hw = [x for x in T if x.is_highest_weight([1,2])]
            sage: all(b.energy_function() == b.energy_function(algorithm='definition')
            ....:     for b in hw)
            True
        """
        C = self.parent().crystals[0]
        ell = ceil(C.s()/C.cartan_type().c()[C.r()])
        is_perfect = all(ell == K.s()/K.cartan_type().c()[K.r()]
                         for K in self.parent().crystals)
        if algorithm is None:
            if is_perfect:
                algorithm = 'grading'
            else:
                algorithm = 'definition'

        if algorithm == 'grading':
            if not is_perfect:
                raise NotImplementedError("all crystals in the tensor product need to be perfect of the same level")
            t = self.parent()(*[K.module_generator() for K in self.parent().crystals])
            d = t.affine_grading()
            return d - self.affine_grading()

        if algorithm == 'definition':
            # Setup
            energy = ZZ.zero()
            R_mats = [[K.R_matrix(Kp) for Kp in self.parent().crystals[i+1:]]
                      for i,K in enumerate(self.parent().crystals)]
            H_funcs = [[K.local_energy_function(Kp) for Kp in self.parent().crystals[i+1:]]
                       for i,K in enumerate(self.parent().crystals)]

            for i,b in enumerate(self):
                for j,R in enumerate(R_mats[i]):
                    H = H_funcs[i][j]
                    bp = self[i+j+1]
                    T = R.domain()
                    t = T(b, bp)
                    energy += H(t)
                    b = R(t)[1]
                energy += b.energy_function()  # D contribution
            return energy
        else:
            raise ValueError("invalid algorithm")

    def affine_grading(self):
        r"""
        Returns the affine grading of `self`.

        The affine grading is only defined when ``self`` is an element of a
        tensor product of affine Kirillov-Reshetikhin crystals. It is
        calculated by finding a path from ``self`` to a ground state path
        using the helper method :meth:`e_string_to_ground_state` and counting
        the number of affine Kashiwara operators `e_0` applied on the way.

        INPUT:

        - ``self`` -- an element of a tensor product of Kirillov-Reshetikhin
          crystals

        OUTPUT: an integer

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
            sage: T = crystals.TensorProduct(K,K)
            sage: t = T.module_generators[0]
            sage: t.affine_grading()
            1

            sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
            sage: T = crystals.TensorProduct(K,K,K)
            sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
            sage: for b in hw:
            ....:     print("{} {}".format(b, b.affine_grading()))
            [[[1]], [[1]], [[1]]] 3
            [[[1]], [[2]], [[1]]] 1
            [[[2]], [[1]], [[1]]] 2
            [[[3]], [[2]], [[1]]] 0

            sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
            sage: T = crystals.TensorProduct(K,K,K)
            sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
            sage: for b in hw:
            ....:     print("{} {}".format(b, b.affine_grading()))
            [[[1]], [[1]], [[1]]] 2
            [[[1]], [[2]], [[1]]] 1
            [[[1]], [[-1]], [[1]]] 0
            [[[2]], [[1]], [[1]]] 1
            [[[-2]], [[2]], [[1]]] 0
            [[[-1]], [[1]], [[1]]] 1
        """
        return self.e_string_to_ground_state().count(0)

    @cached_method
    def e_string_to_ground_state(self):
        r"""
        Returns a string of integers in the index set `(i_1,\ldots,i_k)` such
        that `e_{i_k} \cdots e_{i_1}` of ``self`` is the ground state.

        This method is only defined when ``self`` is an element of a tensor
        product of affine Kirillov-Reshetikhin crystals. It calculates a path
        from ``self`` to a ground state path using Demazure arrows as defined
        in Lemma 7.3 in [SchillingTingley2011]_.

        INPUT:

        - ``self`` -- an element of a tensor product of Kirillov-Reshetikhin
          crystals

        OUTPUT: a tuple of integers `(i_1,\ldots,i_k)`

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
            sage: T = crystals.TensorProduct(K,K)
            sage: t = T.module_generators[0]
            sage: t.e_string_to_ground_state()
            (0, 2)

            sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
            sage: T = crystals.TensorProduct(K,K)
            sage: t = T.module_generators[0]; t
            [[[1]], [[1]]]
            sage: t.e_string_to_ground_state()
            (0,)
            sage: x=t.e(0)
            sage: x.e_string_to_ground_state()
            ()
            sage: y=t.f_string([1,2,1,1,0]); y
            [[[2]], [[1]]]
            sage: y.e_string_to_ground_state()
            ()
        """
        from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux
        if self.parent().crystals[0].__module__ != 'sage.combinat.crystals.kirillov_reshetikhin' and \
                not isinstance(self.parent().crystals[0], KirillovReshetikhinTableaux):
            raise ValueError("All crystals in the tensor product need to be Kirillov-Reshetikhin crystals")
        I = self.cartan_type().classical().index_set()
        ell = max(ceil(K.s()/K.cartan_type().c()[K.r()]) for K in self.parent().crystals)
        for i in I:
            if self.epsilon(i) > 0:
                return (i,) + (self.e(i)).e_string_to_ground_state()
        if self.epsilon(0) > ell:
            return (0,) + (self.e(0)).e_string_to_ground_state()
        return ()

CrystalOfWords.Element = TensorProductOfCrystalsElement

class FullTensorProductOfRegularCrystals(FullTensorProductOfCrystals):
    """
    Full tensor product of regular crystals.
    """
    Element = TensorProductOfRegularCrystalsElement

class TensorProductOfRegularCrystalsWithGenerators(TensorProductOfCrystalsWithGenerators):
    """
    Tensor product of regular crystals with a generating set.
    """
    Element = TensorProductOfRegularCrystalsElement

#########################################################
## Crystal of tableaux

class CrystalOfTableaux(CrystalOfWords):
    r"""
    A class for crystals of tableaux with integer valued shapes

    INPUT:

    - ``cartan_type`` -- a Cartan type
    - ``shape`` -- a partition of length at most ``cartan_type.rank()``
    - ``shapes`` -- a list of such partitions

    This constructs a classical crystal with the given Cartan type and
    highest weight(s) corresponding to the given shape(s).

    If the type is `D_r`, the shape is permitted to have a negative
    value in the `r`-th position. Thus if the shape equals `[s_1,\ldots,s_r]`,
    then `s_r` may be negative but in any case `s_1 \geq \cdots \geq s_{r-1}
    \geq |s_r|`. This crystal is related to that of shape
    `[s_1,\ldots,|s_r|]` by the outer automorphism of `SO(2r)`.

    If the type is `D_r` or `B_r`, the shape is permitted to be of
    length `r` with all parts of half integer value. This corresponds
    to having one spin column at the beginning of the tableau. If
    several shapes are provided, they currently should all or none
    have this property.

    Crystals of tableaux are constructed using an embedding into
    tensor products following Kashiwara and Nakashima [KN94]_. Sage's tensor
    product rule for crystals differs from that of Kashiwara and Nakashima
    by reversing the order of the tensor factors. Sage produces the same
    crystals of tableaux as Kashiwara and Nakashima. With Sage's convention,
    the tensor product of crystals is the same as the monoid operation on
    tableaux and hence the plactic monoid.

    .. SEEALSO::

        :mod:`sage.combinat.crystals.crystals` for general help on
        crystals, and in particular plotting and `\LaTeX` output.

    EXAMPLES:

    We create the crystal of tableaux for type `A_2`, with
    highest weight given by the partition `[2,1,1]`::

        sage: T = crystals.Tableaux(['A',3], shape = [2,1,1])

    Here is the list of its elements::

        sage: T.list()
        [[[1, 1], [2], [3]], [[1, 2], [2], [3]], [[1, 3], [2], [3]],
         [[1, 4], [2], [3]], [[1, 4], [2], [4]], [[1, 4], [3], [4]],
         [[2, 4], [3], [4]], [[1, 1], [2], [4]], [[1, 2], [2], [4]],
         [[1, 3], [2], [4]], [[1, 3], [3], [4]], [[2, 3], [3], [4]],
         [[1, 1], [3], [4]], [[1, 2], [3], [4]], [[2, 2], [3], [4]]]

    Internally, a tableau of a given Cartan type is represented as a
    tensor product of letters of the same type. The order in which the
    tensor factors appear is by reading the columns of the tableaux
    left to right, top to bottom (in French notation). As an example::

        sage: T = crystals.Tableaux(['A',2], shape = [3,2])
        sage: T.module_generators[0]
        [[1, 1, 1], [2, 2]]
        sage: T.module_generators[0]._list
        [2, 1, 2, 1, 1]

    To create a tableau, one can use::

        sage: Tab = crystals.Tableaux(['A',3], shape = [2,2])
        sage: Tab(rows=[[1,2],[3,4]])
        [[1, 2], [3, 4]]
        sage: Tab(columns=[[3,1],[4,2]])
        [[1, 2], [3, 4]]

    .. TODO::

        FIXME:

        - Do we want to specify the columns increasingly or
          decreasingly? That is, should this be
          ``Tab(columns = [[1,3],[2,4]])``?
        - Make this fully consistent with
          :func:`~sage.combinat.tableau.Tableau`!

    We illustrate the use of a shape with a negative last entry in
    type `D`::

        sage: T = crystals.Tableaux(['D',4],shape=[1,1,1,-1])
        sage: T.cardinality()
        35
        sage: TestSuite(T).run()

    We illustrate the construction of crystals of spin tableaux when
    the partitions have half integer values in type `B` and `D`::

        sage: T = crystals.Tableaux(['B',3],shape=[3/2,1/2,1/2]); T
        The crystal of tableaux of type ['B', 3] and shape(s) [[3/2, 1/2, 1/2]]
        sage: T.cardinality()
        48
        sage: T.module_generators
        [[+++, [[1]]]]
        sage: TestSuite(T).run()

        sage: T = crystals.Tableaux(['D',3],shape=[3/2,1/2,-1/2]); T
        The crystal of tableaux of type ['D', 3] and shape(s) [[3/2, 1/2, -1/2]]
        sage: T.cardinality()
        20
        sage: T.module_generators
        [[++-, [[1]]]]
        sage: TestSuite(T).run()

    TESTS:

    Base cases::

        sage: T = crystals.Tableaux(['A',2], shape = [])
        sage: T.list()
        [[]]
        sage: TestSuite(T).run()

        sage: T = crystals.Tableaux(['C',2], shape = [1])
        sage: T.list()
        [[[1]], [[2]], [[-2]], [[-1]]]
        sage: TestSuite(T).run()

        sage: T = crystals.Tableaux(['A',2], shapes = [[],[1],[2]])
        sage: T.list()
        [[], [[1]], [[2]], [[3]], [[1, 1]], [[1, 2]], [[2, 2]], [[1, 3]], [[2, 3]], [[3, 3]]]
        sage: T.module_generators
        ([], [[1]], [[1, 1]])

        sage: T = crystals.Tableaux(['B',2], shape=[3])
        sage: T(rows=[[1,1,0]])
        [[1, 1, 0]]

    Input tests::

        sage: T = crystals.Tableaux(['A',3], shape = [2,2])
        sage: C = T.letters
        sage: Tab(rows    = [[1,2],[3,4]])._list == [C(3),C(1),C(4),C(2)]
        True
        sage: Tab(columns = [[3,1],[4,2]])._list == [C(3),C(1),C(4),C(2)]
        True

    For compatibility with
    :func:`~sage.combinat.crystals.tensor_product.TensorProductOfCrystals` we
    need to accept as input the internal list or sequence of elements::

        sage: Tab(list    = [3,1,4,2])._list     == [C(3),C(1),C(4),C(2)]
        True
        sage: Tab(3,1,4,2)._list                 == [C(3),C(1),C(4),C(2)]
        True

    The next example checks whether a given tableau is in fact a valid
    type `C` tableau or not::

        sage: T = crystals.Tableaux(['C',3], shape = [2,2,2])
        sage: Tab = T(rows=[[1,3],[2,-3],[3,-1]])
        sage: Tab in T.list()
        True
        sage: Tab = T(rows=[[2,3],[3,-3],[-3,-2]])
        sage: Tab in T.list()
        False
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, shapes = None, shape = None):
        """
        Normalizes the input arguments to ensure unique representation,
        and to delegate the construction of spin tableaux.

        EXAMPLES::

            sage: T1 = crystals.Tableaux(CartanType(['A',3]), shape  = [2,2])
            sage: T2 = crystals.Tableaux(['A',3],             shape  = (2,2))
            sage: T3 = crystals.Tableaux(['A',3],             shapes = ([2,2],))
            sage: T2 is T1, T3 is T1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)
        n = cartan_type.rank()
        # standardize shape/shapes input into a tuple of tuples
        assert operator.xor(shape is not None, shapes is not None)
        if shape is not None:
            shapes = (shape,)
        spin_shapes = tuple( tuple(shape) for shape in shapes )
        try:
            shapes = tuple( tuple(trunc(i) for i in shape) for shape in spin_shapes )
        except Exception:
            raise ValueError("shapes should all be partitions or half-integer partitions")
        if spin_shapes == shapes:
            return super(CrystalOfTableaux, cls).__classcall__(cls, cartan_type, shapes)

        # Handle the construction of a crystals of spin tableaux
        # Caveat: this currently only supports all shapes being half
        # integer partitions of length the rank for type B and D. In
        # particular, for type D, the spins all have to be plus or all
        # minus spins
        if any(len(sh) != n for sh in shapes):
            raise ValueError("the length of all half-integer partition shapes should be the rank")
        if any(2*i % 2 != 1 for shape in spin_shapes for i in shape):
            raise ValueError("shapes should be either all partitions or all half-integer partitions")
        if cartan_type.type() == 'D':
            if all( i >= 0 for shape in spin_shapes for i in shape):
                S = CrystalOfSpinsPlus(cartan_type)
            elif all(shape[-1]<0 for shape in spin_shapes):
                S = CrystalOfSpinsMinus(cartan_type)
            else:
                raise ValueError("In type D spins should all be positive or negative")
        else:
            if any( i < 0 for shape in spin_shapes for i in shape):
                raise ValueError("shapes should all be partitions")
            S = CrystalOfSpins(cartan_type)
        B = CrystalOfTableaux(cartan_type, shapes = shapes)
        T = TensorProductOfCrystals(S,B, generators=[[S.module_generators[0],x] for x in B.module_generators])
        T.rename("The crystal of tableaux of type %s and shape(s) %s"%(cartan_type, list(list(shape) for shape in spin_shapes)))
        T.shapes = spin_shapes
        return T


    def __init__(self, cartan_type, shapes):
        """
        Construct the crystal of all tableaux of the given shapes.

        INPUT:

        - ``cartan_type`` -- (data coercible into) a Cartan type
        - ``shapes``      -- a list (or iterable) of shapes
        - ``shape``       -- a shape

        Shapes themselves are lists (or iterable) of integers.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: TestSuite(T).run()
        """
#        super(CrystalOfTableaux, self).__init__(category = FiniteEnumeratedSets())
        Parent.__init__(self, category = ClassicalCrystals())
        self.letters = CrystalOfLetters(cartan_type)
        self.shapes = shapes
        self.module_generators = tuple(self.module_generator(la) for la in shapes)
        self.rename("The crystal of tableaux of type %s and shape(s) %s"%(cartan_type, list(list(shape) for shape in shapes)))

    def cartan_type(self):
        """
        Returns the Cartan type of the associated crystal

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: T.cartan_type()
            ['A', 3]
        """
        return self.letters.cartan_type()

    def module_generator(self, shape):
        """
        This yields the module generator (or highest weight element) of a classical
        crystal of given shape. The module generator is the unique tableau with equal
        shape and content.

        EXAMPLE::

            sage: T = crystals.Tableaux(['D',3], shape = [1,1])
            sage: T.module_generator([1,1])
            [[1], [2]]

            sage: T = crystals.Tableaux(['D',4],shape=[2,2,2,-2])
            sage: T.module_generator(tuple([2,2,2,-2]))
            [[1, 1], [2, 2], [3, 3], [-4, -4]]
            sage: T.cardinality()
            294
            sage: T = crystals.Tableaux(['D',4],shape=[2,2,2,2])
            sage: T.module_generator(tuple([2,2,2,2]))
            [[1, 1], [2, 2], [3, 3], [4, 4]]
            sage: T.cardinality()
            294
        """
        type = self.cartan_type()
        if type[0] == 'D' and len(shape) == type[1] and shape[type[1]-1] < 0:
            invert = True
            shape = shape[:-1]+(-shape[type[1]-1],)
        else:
            invert = False
        p = Partition(shape).conjugate()
        # The column canonical tableau, read by columns
        module_generator = flatten([[p[j]-i for i in range(p[j])] for j in range(len(p))])
        if invert:
            module_generator = [(-x if x == type[1] else x) for x in module_generator]
        return self(list=[self.letters(x) for x in module_generator])

    def _element_constructor_(self, *args, **options):
        """
        Return a
        :class:`~sage.combinat.crystals.tensor_product.CrystalOfTableauxElement`.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: T(rows=[[1,2],[3,4]])
            [[1, 2], [3, 4]]
            sage: T(columns=[[3,1],[4,2]])
            [[1, 2], [3, 4]]
        """
        return self.element_class(self, *args, **options)



class CrystalOfTableauxElement(TensorProductOfRegularCrystalsElement):
    """
    Element in a crystal of tableaux.
    """
    def __init__(self, parent, *args, **options):
        """
        There are several ways to input tableaux, by rows,
        by columns, as the list of column elements, or as a sequence of numbers
        in column reading.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t
            [[1, 2], [3, 4]]
            sage: TestSuite(t).run()

            sage: t = T(columns=[[3,1],[4,2]])
            sage: t
            [[1, 2], [3, 4]]
            sage: TestSuite(t).run()

            sage: t = T(list=[3,1,4,2])
            sage: t
            [[1, 2], [3, 4]]

            sage: t = T(3,1,4,2)
            sage: t
            [[1, 2], [3, 4]]

        Currently inputting the empty tableau as an empty sequence is
        broken due to a bug in the generic __call__ method (see :trac:`8648`).

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape=[])
            sage: t = T()
            sage: t._list
            [0]

        TESTS:

        Integer types that are not a Sage ``Integer`` (such as a Python ``int``
        and typically arise from compiled code) were not converted into a
        letter. This caused certain functions to fail. This is fixed in
        :trac:`13204`::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(list=[int(3),1,4,2])
            sage: type(t[0])
            <type 'sage.combinat.crystals.letters.Crystal_of_letters_type_A_element'>
            sage: t = T(list=[3,int(1),4,2])
            sage: type(t[1])
            <type 'sage.combinat.crystals.letters.Crystal_of_letters_type_A_element'>
            sage: C = crystals.KirillovReshetikhin(['A',int(3),1], 1,1)
            sage: C[0].e(0)
            [[4]]
        """
        if len(args) == 1:
            if isinstance(args[0], Tableau):
                options['rows'] = args[0]
        if 'list' in options:
            the_list = options['list']
        elif 'rows' in options:
            rows = options['rows']
#            the_list=Tableau(rows).to_word_by_column()
            rows = Tableau(rows).conjugate()
            the_list = []
            for col in rows:
                the_list += reversed(col)
        elif 'columns' in options:
            columns = options['columns']
            the_list = []
            for col in columns:
                the_list += col
        else:
            the_list = [i for i in args]
        TensorProductOfRegularCrystalsElement.__init__(self, parent, [parent.letters(_) for _ in the_list])

    def _repr_(self):
        """
        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t._repr_()
            '[[1, 2], [3, 4]]'
        """
        return repr(self.to_tableau())

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A', 4], shape=[4,2,1])
            sage: elt = C(rows=[[1,1,1,2], [2,3], [4]])
            sage: print(elt._repr_diagram())
              1  1  1  2
              2  3
              4
        """
        return self.to_tableau()._repr_diagram()

    def pp(self):
        """
        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t.pp()
            1  2
            3  4
        """
        return self.to_tableau().pp()

    def _ascii_art_(self):
        """
        Return an ascii art version of ``self``.

        EXAMPLES:

        We check that :trac:`16486` is fixed::

            sage: T = crystals.Tableaux(['B',6], shape=[1]*5)
            sage: ascii_art(T.module_generators[0])
              1
              2
              3
              4
              5
            sage: T = crystals.Tableaux(['D',4], shape=[2,1])
            sage: t = T.module_generators[0].f_string([1,2,3,4,2,2,3,4])
            sage: ascii_art(t)
              1 -2
             -3
        """
        return self.to_tableau()._ascii_art_()

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [4,2])
            sage: t = T(rows=[[1,1,2,3],[2,3]])
            sage: latex(t) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{1}&\lr{1}&\lr{2}&\lr{3}\\\cline{1-4}
            \lr{2}&\lr{3}\\\cline{1-2}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        # Modified version of to_tableau() to have the entrys be letters
        #   rather than their values
        if self._list == []:
            return "{\\emptyset}"

        tab = [ [self[0]] ]
        for i in range(1,len(self)):
            if self[i-1] < self[i] or (self[i-1].value != 0 and self[i-1] == self[i]):
                tab.append([self[i]])
            else:
                l = len(tab)-1
                tab[l].append(self[i])
        for x in tab:
            x.reverse()
        T = Tableau(tab).conjugate()
        return tex_from_array([[letter._latex_() for letter in row] for row in T])

    @cached_method
    def to_tableau(self):
        """
        Returns the Tableau object corresponding to self.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]]).to_tableau(); t
            [[1, 2], [3, 4]]
            sage: type(t)
            <class 'sage.combinat.tableau.Tableaux_all_with_category.element_class'>
            sage: type(t[0][0])
            <type 'int'>
            sage: T = crystals.Tableaux(['D',3], shape = [1,1])
            sage: t=T(rows=[[-3],[3]]).to_tableau(); t
            [[-3], [3]]
            sage: t=T(rows=[[3],[-3]]).to_tableau(); t
            [[3], [-3]]
            sage: T = crystals.Tableaux(['B',2], shape = [1,1])
            sage: t = T(rows=[[0],[0]]).to_tableau(); t
            [[0], [0]]
        """
        if self._list == []:
            return Tableau([])
        tab = [ [self[0].value] ]
        for i in range(1,len(self)):
            if self[i-1] < self[i] or (self[i-1].value != 0 and self[i-1] == self[i]):
                tab.append([self[i].value])
            else:
                l = len(tab)-1
                tab[l].append(self[i].value)
        for x in tab:
            x.reverse()
        return Tableau(tab).conjugate()

    def promotion(self):
        """
        Promotion for type A crystals of tableaux of rectangular shape

        Returns the result of applying promotion on this tableau.

        This method only makes sense in type A with rectangular shapes.

        EXAMPLES::

            sage: C = crystals.Tableaux(["A",3], shape = [3,3,3])
            sage: t = C(Tableau([[1,1,1],[2,2,3],[3,4,4]]))
            sage: t
            [[1, 1, 1], [2, 2, 3], [3, 4, 4]]
            sage: t.promotion()
            [[1, 1, 2], [2, 2, 3], [3, 4, 4]]
            sage: t.promotion().parent()
            The crystal of tableaux of type ['A', 3] and shape(s) [[3, 3, 3]]
        """
        crystal = self.parent()
        cartan_type = crystal.cartan_type()
        assert cartan_type.type() == 'A'
        return crystal(self.to_tableau().promotion(cartan_type.rank()))

    def promotion_inverse(self):
        """
        Inverse promotion for type A crystals of tableaux of rectangular shape

        Returns the result of applying inverse promotion on this tableau.

        This method only makes sense in type A with rectangular shapes.

        EXAMPLES::

            sage: C = crystals.Tableaux(["A",3], shape = [3,3,3])
            sage: t = C(Tableau([[1,1,1],[2,2,3],[3,4,4]]))
            sage: t
            [[1, 1, 1], [2, 2, 3], [3, 4, 4]]
            sage: t.promotion_inverse()
            [[1, 1, 2], [2, 3, 3], [4, 4, 4]]
            sage: t.promotion_inverse().parent()
            The crystal of tableaux of type ['A', 3] and shape(s) [[3, 3, 3]]
        """
        crystal = self.parent()
        cartan_type = crystal.cartan_type()
        assert cartan_type.type() == 'A'
        return crystal(self.to_tableau().promotion_inverse(cartan_type.rank()))

CrystalOfTableaux.Element = CrystalOfTableauxElement

#####################################################################
## Local energy function

class LocalEnergyFunction(Map):
    r"""
    The local energy function.

    Let `B` and `B'` be Kirillov-Reshetikhin crystals with maximal
    vectors `u_B` and `u_{B'}` respectively. The *local energy function*
    `H : B \otimes B' \to \ZZ` is the function which satisfies

    .. MATH::

        H(e_0(b \otimes b')) = H(b \otimes b') + \begin{cases}
        1 & \text{if } i = 0 \text{ and LL}, \\
        -1 & \text{if } i = 0 \text{ and RR}, \\
        0 & \text{otherwise,}
        \end{cases}

    where LL (resp. RR) denote `e_0` acts on the left (resp. right)
    on both `b \otimes b'` and `R(b \otimes b')`, and
    normalized by `H(u_B \otimes u_{B'}) = 0`.

    INPUT:

    - ``B`` -- a Kirillov-Reshetikhin crystal
    - ``Bp`` -- a Kirillov-Reshetikhin crystal
    - ``normalization`` -- (default: 0) the normalization value

    EXAMPLES::

        sage: K = crystals.KirillovReshetikhin(['C',2,1], 1,2)
        sage: K2 = crystals.KirillovReshetikhin(['C',2,1], 2,1)
        sage: H = K.local_energy_function(K2)
        sage: T = tensor([K, K2])
        sage: hw = [x for x in T if x.is_highest_weight([1,2])]
        sage: for b in hw:
        ....:     b, H(b)
        ([[], [[1], [2]]], 1)
        ([[[1, 1]], [[1], [2]]], 0)
        ([[[2, -2]], [[1], [2]]], 1)
        ([[[1, -2]], [[1], [2]]], 1)

    REFERENCES:

    .. [KKMMNN92] S-J. Kang, M. Kashiwara, K. C. Misra, T. Miwa,
       T. Nakashima, and A. Nakayashiki.
       *Affine crystals and vertex models*.
       Int. J. Mod. Phys. A, **7** (suppl. 1A), (1992) pp. 449-484.
    """
    def __init__(self, B, Bp, normalization=0):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A',7,2], 1,2)
            sage: K2 = crystals.KirillovReshetikhin(['A',7,2], 2,1)
            sage: H = K.local_energy_function(K2)
            sage: TestSuite(H).run(skip=['_test_category', '_test_pickling'])
        """
        self._B = B
        self._Bp = Bp
        self._R_matrix = self._B.R_matrix(self._Bp)
        T = B.tensor(Bp)
        self._known_values = {T(*[K.module_generator() for K in T.crystals]):
                              ZZ(normalization)}
        self._I0 = T.cartan_type().classical().index_set()
        from sage.categories.homset import Hom
        Map.__init__(self, Hom(T, ZZ))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A', 6, 2], 2, 1)
            sage: Kp = crystals.KirillovReshetikhin(['A', 6, 2], 1, 1)
            sage: H = K.local_energy_function(Kp); H
            Local energy function of
             Kirillov-Reshetikhin crystal of type ['BC', 3, 2] with (r,s)=(2,1)
            tensor
             Kirillov-Reshetikhin crystal of type ['BC', 3, 2] with (r,s)=(1,1)
        """
        return "Local energy function of {} tensor {}".format(self._B, self._Bp)

    def _call_(self, x):
        """
        Return the local energy of ``x``.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['B',4,1], 1,2)
            sage: K2 = crystals.KirillovReshetikhin(['B',4,1], 2,1)
            sage: H = K.local_energy_function(K2)
            sage: T = tensor([K, K2])
            sage: hw = [x for x in T if x.is_highest_weight([1,2])]
            sage: H(hw[0])
            1
        """
        # Setup variables
        visited = {x: 0}
        check0 = [x]

        # Helper function
        def to_classical_hw(cur):
            for i in self._I0:
                b = cur.e(i)
                if b is not None and b not in visited:
                    visited[b] = visited[cur] # No change
                    return b
            return None # is classically HW or all have been visited

        cur = x
        # Get the affine node (it might not be 0 if the type
        #   has been relabeled)
        i0 = x.parent().cartan_type().special_node()
        while cur not in self._known_values:
            # We first go towards the classically highest weight since
            #   the maximal vector is classically highest weight
            b = to_classical_hw(cur)

            # If classically HW, then try 0 arrows
            while b is None:
                b = check0.pop()
                c = b.e(i0)
                # If there is no 0 arrow or we have already seen c, move along
                if c is None or c in visited:
                    b = None
                    continue

                bp = self._R_matrix(b)
                cp = bp.e(i0)
                if b[1] == c[1] and bp[1] == cp[1]: # LL case
                    visited[c] = visited[b] + 1
                elif b[0] == c[0] and bp[0] == cp[0]: # RR case
                    visited[c] = visited[b] - 1
                else:
                    visited[c] = visited[b] # Otherwise no change
                b = c

            cur = b
            check0.append(b)

        baseline = self._known_values[cur] - visited[cur]
        for y in visited:
            self._known_values[y] = baseline + visited[y]

        return self._known_values[x]


# deprecations from trac:18555
from sage.misc.superseded import deprecated_function_alias
TensorProductOfCrystals.global_options=deprecated_function_alias(18555, TensorProductOfCrystals.options)
TensorProductOfCrystalsOptions=deprecated_function_alias(18555, TensorProductOfCrystals.options)
