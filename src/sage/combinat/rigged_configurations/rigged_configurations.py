# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Rigged Configurations

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2010-2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

import itertools

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.global_options import GlobalOptions
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.combinat.misc import IterableFunctionCall
import sage.combinat.tableau as tableau
from sage.rings.rational_field import QQ
from sage.categories.loop_crystals import KirillovReshetikhinCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.kleber_tree import KleberTree, VirtualKleberTree
from sage.combinat.rigged_configurations.rigged_configuration_element import (
     RiggedConfigurationElement, KRRCSimplyLacedElement, KRRCNonSimplyLacedElement,
     KRRCTypeA2DualElement)
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition


# Used in the KR crystals catalog so that there is a common interface
def KirillovReshetikhinCrystal(cartan_type, r, s):
    r"""
    Return the KR crystal `B^{r,s}` using
    :class:`rigged configurations <RiggedConfigurations>`.

    This is the rigged configuration `RC(B^{r,s})` or `RC(L)` with
    `L = (L_i^{(a)})` and `L_i^{(a)} = \delta_{a,r} \delta_{i,s}`.

    EXAMPLES::

        sage: K1 = crystals.kirillov_reshetikhin.RiggedConfigurations(['A',6,2], 2, 1)
        sage: K2 = crystals.kirillov_reshetikhin.LSPaths(['A',6,2], 2, 1)
        sage: K1.digraph().is_isomorphic(K2.digraph(), edge_labels=True)
        True

    TESTS:

    We explicitly import and check we get the same crystal::

        sage: from sage.combinat.rigged_configurations.rigged_configurations import KirillovReshetikhinCrystal
        sage: K1 = crystals.kirillov_reshetikhin.RiggedConfigurations(['A',6,2], 2, 1)
        sage: K1 is KirillovReshetikhinCrystal(['A',6,2], 2, 1)
        True
    """
    return RiggedConfigurations(cartan_type, [[r, s]])


# Note on implementation, this class is used for simply-laced types only
class RiggedConfigurations(UniqueRepresentation, Parent):
    r"""
    Rigged configurations as `U_q^{\prime}(\mathfrak{g})`-crystals.

    Let `\overline{I}` denote the classical index set associated to the Cartan
    type of the rigged configurations. A rigged configuration of multiplicity
    array `L_i^{(a)}` and dominant weight `\Lambda` is a sequence of partitions
    `\{ \nu^{(a)} \mid a \in \overline{I} \}` such that

    .. MATH::

        \sum_{\overline{I} \times \ZZ_{>0}} i m_i^{(a)} \alpha_a
        = \sum_{\overline{I} \times \ZZ_{>0}} i L_i^{(a)} \Lambda_a
        - \Lambda

    where `\alpha_a` is a simple root, `\Lambda_a` is a fundamental weight,
    and `m_i^{(a)}` is the number of rows of length `i` in the partition
    `\nu^{(a)}`.

    Each partition `\nu^{(a)}`, in the sequence also comes with a sequence of
    statistics `p_i^{(a)}` called *vacancy numbers* and a weakly decreasing
    sequence `J_i^{(a)}` of length `m_i^{(a)}` called *riggings*.
    Vacancy numbers are computed based upon the partitions and `L_i^{(a)}`,
    and the riggings must satisfy `\max J_i^{(a)} \leq p_i^{(a)}`. We call
    such a partition a *rigged partition*. For more, see
    [RigConBijection]_ [CrysStructSchilling06]_ [BijectionLRT]_.

    Rigged configurations form combinatorial objects first introduced by
    Kerov, Kirillov and Reshetikhin that arose from studies of statistical
    mechanical models using the Bethe Ansatz. They are sequences of rigged
    partitions. A rigged partition is a partition together with a label
    associated to each part that satisfy certain constraints. The labels
    are also called riggings.

    Rigged configurations exist for all affine Kac-Moody Lie algebras. See
    for example [HKOTT2002]_. In Sage they are specified by providing a Cartan
    type and a list of rectangular shapes `B`. The list of all (highest
    weight) rigged configurations for given `B` is computed via the (virtual)
    Kleber algorithm (see also
    :class:`~sage.combinat.rigged_configurations.kleber_tree.KleberTree` and
    :class:`~sage.combinat.rigged_configurations.kleber_tree.VirtualKleberTree`).

    Rigged configurations in simply-laced types all admit a classical crystal
    structure [CrysStructSchilling06]_. For non-simply-laced types, the
    crystal is given by using virtual rigged configurations [OSS03]_. The
    highest weight rigged configurations are those where all riggings are
    nonnegative. The list of all rigged configurations is computed from the
    highest weight ones using the crystal operators.

    Rigged configurations are conjecturally in bijection with
    :class:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux.TensorProductOfKirillovReshetikhinTableaux`
    of non-exceptional affine types where the list `B` corresponds to the
    tensor factors `B^{r,s}`. The bijection has been proven in types `A_n^{(1)}`
    and `D_n^{(1)}` and when the only nonzero entries of `L_i^{(a)}` are either
    only `L_1^{(a)}` or only `L_i^{(1)}` (corresponding to single columns or
    rows respectively) [RigConBijection]_, [BijectionLRT]_, [BijectionDn]_.

    KR crystals are implemented in Sage, see
    :func:`~sage.combinat.crystals.kirillov_reshetikhin.KirillovReshetikhinCrystal`,
    however, in the bijection with rigged configurations a different
    realization of the elements in the crystal are obtained, which are
    coined KR tableaux, see
    :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableaux`.
    For more details see [OSS2011]_.

    .. NOTE::

        All non-simply-laced rigged configurations have not been proven to
        give rise to aligned virtual crystals (i.e. have the correct crystal
        structure or isomorphic as affine crystals to the tensor product of
        KR tableaux).

    INPUT:

    - ``cartan_type`` -- a Cartan type

    - ``B`` -- list of positive integer tuples `(r,s)` corresponding to the
      tensor factors in the bijection with tensor product of
      Kirillov-Reshetikhin tableaux or equivalently the sequence of width `s`
      and height `r` rectangles

    REFERENCES:

    .. [HKOTT2002] \G. Hatayama, A. Kuniba, M. Okado, T. Takagi, Z. Tsuboi.
       Paths, Crystals and Fermionic Formulae.
       Prog. Math. Phys. **23** (2002) Pages 205-272.

    .. [CrysStructSchilling06] Anne Schilling.
       Crystal structure on rigged configurations.
       International Mathematics Research Notices.
       Volume 2006. (2006) Article ID 97376. Pages 1-27.

    .. [RigConBijection] Masato Okado, Anne Schilling, Mark Shimozono.
       A crystal to rigged configuration bijection for non-exceptional affine
       algebras.
       Algebraic Combinatorics and Quantum Groups.
       Edited by N. Jing. World Scientific. (2003) Pages 85-124.

    .. [BijectionDn] Anne Schilling.
       A bijection between type `D_n^{(1)}` crystals and rigged configurations.
       J. Algebra. **285** (2005) 292-334

    .. [BijectionLRT] Anatol N. Kirillov, Anne Schilling, Mark Shimozono.
       A bijection between Littlewood-Richardson tableaux and rigged
       configurations.
       Selecta Mathematica (N.S.). **8** (2002) Pages 67-135.
       (:mathscinet:`MR1890195`).

    EXAMPLES::

        sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]])
        sage: RC
        Rigged configurations of type ['A', 3, 1] and factor(s) ((3, 2), (1, 2), (1, 1))

        sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1]]); RC
        Rigged configurations of type ['A', 3, 1] and factor(s) ((2, 1),)
        sage: RC.cardinality()
        6
        sage: len(RC.list()) == RC.cardinality()
        True
        sage: RC.list() # random
        [
        <BLANKLINE>
                                                 0[ ]0
        (/)  (/)      (/)      -1[ ]-1  -1[ ]-1
                                                 -1[ ]-1
        (/)  -1[ ]-1  0[ ]0    0[ ]0    1[ ]1    -1[ ]-1
        <BLANKLINE>
        (/)  (/)      -1[ ]-1  (/)      -1[ ]-1  0[ ]0
           ,        ,        ,        ,        ,
        ]

    A rigged configuration element with all riggings equal to the vacancy
    numbers can be created as follows::

        sage: RC = RiggedConfigurations(['A', 3, 1], [[3,2], [2,1], [1,1], [1,1]]); RC
        Rigged configurations of type ['A', 3, 1] and factor(s) ((3, 2), (2, 1), (1, 1), (1, 1))
        sage: elt = RC(partition_list=[[1],[],[]]); elt
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>
        (/)
        <BLANKLINE>
        (/)
        <BLANKLINE>

    If on the other hand we also want to specify the riggings, this can be
    achieved as follows::

        sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]])
        sage: RC(partition_list=[[2],[2],[2]])
        <BLANKLINE>
        1[ ][ ]1
        <BLANKLINE>
        0[ ][ ]0
        <BLANKLINE>
        0[ ][ ]0
        sage: RC(partition_list=[[2],[2],[2]], rigging_list=[[0],[0],[0]])
        <BLANKLINE>
        1[ ][ ]0
        <BLANKLINE>
        0[ ][ ]0
        <BLANKLINE>
        0[ ][ ]0

    A larger example::

        sage: RC = RiggedConfigurations(['D', 7, 1], [[3,3],[5,2],[4,3],[2,3],[4,4],[3,1],[1,4],[2,2]])
        sage: elt = RC(partition_list=[[2],[3,2,1],[2,2,1,1],[2,2,1,1,1,1],[3,2,1,1,1,1],[2,1,1],[2,2]],
        ....:          rigging_list=[[2],[1,0,0],[4,1,2,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0],[0,0]])
        sage: elt
        <BLANKLINE>
        3[ ][ ]2
        <BLANKLINE>
        1[ ][ ][ ]1
        2[ ][ ]0
        1[ ]0
        <BLANKLINE>
        4[ ][ ]4
        4[ ][ ]1
        3[ ]2
        3[ ]1
        <BLANKLINE>
        2[ ][ ]1
        2[ ][ ]0
        0[ ]0
        0[ ]0
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ][ ][ ]0
        2[ ][ ]1
        0[ ]0
        0[ ]0
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ][ ]0
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ][ ]0
        0[ ][ ]0
        <BLANKLINE>

    To obtain the KR tableau under the bijection between rigged configurations
    and KR tableaux, we can type the following. This example was checked
    against Reiho Sakamoto's Mathematica program on rigged configurations::

        sage: output = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); output
        [[1, 1, 1], [2, 3, 3], [3, 4, -5]] (X) [[1, 1], [2, 2], [3, 3], [5, -6], [6, -5]] (X)
        [[1, 1, 2], [2, 2, 3], [3, 3, 7], [4, 4, -7]] (X) [[1, 1, 1], [2, 2, 2]] (X)
        [[1, 1, 1, 3], [2, 2, 3, 4], [3, 3, 4, 5], [4, 4, 5, 6]] (X) [[1], [2], [3]] (X) [[1, 1, 1, 1]] (X) [[1, 1], [2, 2]]
        sage: elt.to_tensor_product_of_kirillov_reshetikhin_tableaux().to_rigged_configuration() == elt
        True
        sage: output.to_rigged_configuration().to_tensor_product_of_kirillov_reshetikhin_tableaux() == output
        True

    We can also convert between rigged configurations and tensor products of
    KR crystals::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
        sage: elt = RC(partition_list=[[1],[1,1],[1],[1]])
        sage: tp_krc = elt.to_tensor_product_of_kirillov_reshetikhin_crystals(); tp_krc
        [[]]
        sage: ret = RC(tp_krc)
        sage: ret == elt
        True

    ::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[4,1], [3,3]])
        sage: KR1 = crystals.KirillovReshetikhin(['D', 4, 1], 4, 1)
        sage: KR2 = crystals.KirillovReshetikhin(['D', 4, 1], 3, 3)
        sage: T = crystals.TensorProduct(KR1, KR2)
        sage: t = T[1]; t
        [[++++, []], [+++-, [[1], [2], [4], [-4]]]]
        sage: ret = RC(t)
        sage: ret.to_tensor_product_of_kirillov_reshetikhin_crystals()
        [[++++, []], [+++-, [[1], [2], [4], [-4]]]]

    TESTS::

        sage: RC = RiggedConfigurations(['A', 3, 1], [[3,2], [2,1], [1,1], [1,1]])
        sage: len(RC.module_generators)
        17
        sage: RC = RiggedConfigurations(['D', 4, 1], [[1, 1]])
        sage: RC.cardinality()
        8

        sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
        sage: c = RC.cardinality(); c
        29
        sage: K = crystals.KirillovReshetikhin(['D',4,1],2,1)
        sage: K.cardinality() == c
        True
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, B):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = RiggedConfigurations(CartanType(['A',3,1]), [[2,2]])
            sage: RC2 = RiggedConfigurations(['A',3,1], [(2,2)])
            sage: RC3 = RiggedConfigurations(['A',3,1], ((2,2),))
            sage: RC2 is RC1, RC3 is RC1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_affine():
            raise ValueError("The Cartan type must be affine")

        # Standardize B input into a tuple of tuples
        B = tuple(tuple(factor) for factor in B)
        if not B:
            raise ValueError("must contain at least one factor")

        if cartan_type.type() == 'BC': # Type `A_{2n}^{(2)}`
            return RCTypeA2Even(cartan_type, B)
        if cartan_type.dual().type() == 'BC': # Type 'A_{2n}^{(2)\dagger`
            return RCTypeA2Dual(cartan_type, B)
        # We check the classical type to account for A^{(1)}_1 which is not
        #    a virtual rigged configuration.
        if not cartan_type.classical().is_simply_laced():
            return RCNonSimplyLaced(cartan_type, B)

        return super().__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, B):
        r"""
        Initialize the RiggedConfigurations class.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3,1], [1,2]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['A',1,1], [[1,1], [1,1]])
            sage: TestSuite(RC).run()
            sage: RC = RiggedConfigurations(['A',2,1], [[1,1], [2,1]])
            sage: TestSuite(RC).run()
            sage: RC = RiggedConfigurations(['D', 4, 1], [[2,1], [1,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D', 4, 1], [[3,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D', 4, 1], [[4,2]])
            sage: TestSuite(RC).run() # long time
        """
        self._cartan_type = cartan_type
        self.dims = B
        cl = cartan_type.classical()
        self._rc_index = cl.index_set()
        self._rc_index_inverse = {i: ii for ii,i in enumerate(self._rc_index)}
        # We store the Cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = cl.cartan_matrix()
        Parent.__init__(self, category=KirillovReshetikhinCrystals().TensorProducts())

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the options for rigged configurations.
        If no parameters are set, then the function returns a copy of
        the options dictionary.

        The ``options`` to partitions can be accessed as the method
        :obj:`RiggedConfigurations.options` of
        :class:`RiggedConfigurations`.

        @OPTIONS@

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',3,1], [[2,2],[1,1],[1,1]])
            sage: elt = RC(partition_list=[[3,1], [3], [1]])
            sage: elt
            <BLANKLINE>
            -3[ ][ ][ ]-3
            -1[ ]-1
            <BLANKLINE>
            1[ ][ ][ ]1
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: RiggedConfigurations.options(display='horizontal', convention='french')
            sage: elt
            -1[ ]-1         1[ ][ ][ ]1   -1[ ]-1
            -3[ ][ ][ ]-3

        Changing the ``convention`` for rigged configurations also changes the
        ``convention`` option for tableaux and vice versa::

            sage: T = Tableau([[1,2,3],[4,5]])
            sage: T.pp()
              4  5
              1  2  3
            sage: Tableaux.options.convention="english"
            sage: elt
            -3[ ][ ][ ]-3   1[ ][ ][ ]1   -1[ ]-1
            -1[ ]-1
            sage: T.pp()
              1  2  3
              4  5
            sage: RiggedConfigurations.options._reset()
        """
        NAME = 'RiggedConfigurations'
        module = 'sage.combinat.rigged_configurations.rigged_configurations'
        display = dict(default='vertical',
                     description='Specifies how rigged configurations should be printed',
                     values=dict(vertical='displayed vertically',
                                 horizontal='displayed horizontally'),
                     case_sensitive=False)
        element_ascii_art = dict(default=True,
                         description='display using the repr option ``element_ascii_art``',
                         checker=lambda x: isinstance(x, bool))
        half_width_boxes_type_B = dict(default=True,
                description='display the last rigged partition in affine type B as half width boxes',
                checker=lambda x: isinstance(x, bool))
        convention = dict(link_to=(tableau.Tableaux.options,'convention'))
        notation = dict(alt_name='convention')

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]])
            Rigged configurations of type ['A', 3, 1] and factor(s) ((3, 2), (1, 2), (1, 1))
        """
        return "Rigged configurations of type {} and factor(s) {}".format(self._cartan_type, self.dims)

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1]])
            sage: RC._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return self.options.element_ascii_art
        return super()._repr_option(key)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1], [1,1]])
            sage: L = [x for x in RC]
            sage: len(L)
            24
        """
        index_set = self._rc_index
        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
        return RecursivelyEnumeratedSet(self.module_generators,
                    lambda x: [x.f(i) for i in index_set],
                    structure='graded').breadth_first_search_iterator()

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for this set of rigged configurations.

        Iterate over the highest weight rigged configurations by moving
        through the
        :class:`~sage.combinat.rigged_configurations.kleber_tree.KleberTree`
        and then setting appropriate values of the partitions.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2,1]])
            sage: for x in RC.module_generators: x
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>

        TESTS:

        We check that this works with relabelled Cartan types (:issue:`16876`)::

            sage: ct = CartanType(['A',3,1]).relabel(lambda x: x+2)
            sage: RC = RiggedConfigurations(ct, [[4,1],[5,1]])
            sage: len(RC.module_generators)
            2
            sage: ct = CartanType(['A',3,1]).relabel(lambda x: (x+2) % 4)
            sage: RC = RiggedConfigurations(ct, [[0,1],[1,1]])
            sage: len(RC.module_generators)
            2
        """
        module_gens = []
        n = len(self._rc_index)

        for tree_node in self.kleber_tree():
            shapes = []
            cur = tree_node
            path_lambda = [cur.up_root.to_vector()] # Build the lambda values
            # Note that these are not same lambda as in the paper,
            #   but a less computational version.
            while cur.parent_node is not None:
                path_lambda.insert(0, (cur.parent_node.up_root - cur.up_root).to_vector())
                cur = cur.parent_node

            for a in range(n):
                shapes.append([])
                for i, cur_lambda in enumerate(path_lambda):
                    for j in range(cur_lambda[a]):
                        shapes[-1].insert(0, i)

            # Start with a base to calculate the vacancy numbers
            # Make a copy just to be safe
            base = self.element_class(self, partition_list=shapes[:])

            # Build out the blocks for the partition values
            vac_nums = []
            blocks = []
            L = []

            for partition in base:
                vac_nums.append(partition.vacancy_numbers)
                blocks.append([[]])

                # If the partition is empty, there's nothing to do
                if not partition._list:
                    L.append([[]])
                    continue

                # Setup the first block
                block_len = partition[0]
                for i, rowLen in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != rowLen:
                        blocks[-1].append([])
                        block_len = rowLen

                    blocks[-1][-1].append(partition.vacancy_numbers[i])

                L2 = []
                for block in blocks[-1]:
                    L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(itertools.product(*L2))

            C = itertools.product(*L)
            for curBlocks in C:
                module_gens.append( self.element_class(self, KT_constructor=[shapes[:],
                                         self._blocks_to_values(curBlocks[:]),
                                         vac_nums[:]]) )

        return tuple(module_gens)

    def _block_iterator(self, container):
        r"""
        Iterate over all possible riggings for a particular block.

        Helper iterator which iterates over all possible partitions contained
        within the container.

        INPUT:

        - ``container`` -- list of widths of the rows of the container

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: for x in RC._block_iterator([]): x
            []
            sage: for x in RC._block_iterator([2,3]): x
            [0, 0]
            [1, 0]
            [1, 1]
            [2, 0]
            [2, 1]
            [2, 2]
        """
        if not container:
            yield []
            return

        pos = 0
        length = len(container)
        ret_part = [-1] * length
        while pos >= 0:
            ret_part[pos] += 1

            if ret_part[pos] > container[pos] or (pos != 0 and ret_part[pos] > ret_part[pos - 1]):
                ret_part[pos] = -1
                pos -= 1
            else:
                pos += 1

            if pos == length:
                yield ret_part[:]
                pos -= 1

    def _blocks_to_values(self, blocks):
        r"""
        Convert an array of blocks into a list of partition values.

        INPUT:

        - ``blocks`` -- the (2-dim) array blocks of the partition values

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC._blocks_to_values([[[2, 1]]])
            [[2, 1]]
        """
        values = []
        for part_block in blocks:
            if not part_block:
                values.append([])
            else:
                values.append(part_block[0][:]) # Need to make a copy
                for block in part_block[1:]:
                    values[-1].extend(block)
        return values

    def classically_highest_weight_vectors(self):
        """
        Return the classically highest weight elements of ``self``.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: ascii_art(RC.classically_highest_weight_vectors())
            ( (/)  (/)  (/)  (/) )
        """
        return self.module_generators

    def _element_constructor_(self, *lst, **options):
        """
        Construct a ``RiggedConfigurationElement``.

        Typically the user should not call this method since it does not check
        if it is a valid configuration. Instead the user should use the
        iterator methods.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: RC(partition_list=[[1], [1], [], []], rigging_list=[[-1], [0], [], []])
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>

        TESTS::

            sage: KT = crystals.TensorProductOfKirillovReshetikhinTableaux(['C',2,1], [[2,4],[1,2]])
            sage: t = KT(pathlist=[[2,1,2,1,-2,2,-1,-2],[2,-2]])
            sage: rc = t.to_rigged_configuration(); rc
            <BLANKLINE>
            -1[ ][ ][ ]-1
             0[ ][ ]0
            <BLANKLINE>
            -1[ ][ ]-1
            -1[ ]-1
            -1[ ]-1
            <BLANKLINE>
            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[2,4]])
            sage: RC(rc)
            <BLANKLINE>
            -1[ ][ ][ ]-1
             0[ ][ ]0
            <BLANKLINE>
            -1[ ][ ]-1
            -1[ ]-1
            -1[ ]-1
            <BLANKLINE>

        TESTS:

        Check that :issue:`17054` is fixed::

            sage: B = crystals.infinity.RiggedConfigurations(['A',2])
            sage: RC = RiggedConfigurations(['A',2,1], [[1,1]]*4 + [[2,1]]*4)
            sage: x = B.an_element().f_string([2,2,1,1,2,1,2,1])
            sage: ascii_art(x)
            -4[ ][ ][ ][ ]-4  -4[ ][ ][ ][ ]0
            sage: ascii_art(RC(x))
            0[ ][ ][ ][ ]-4  0[ ][ ][ ][ ]0
            sage: x == B.an_element().f_string([2,2,1,1,2,1,2,1])
            True
        """
        if not lst:
            return self.element_class(self, [], **options)

        from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element import TensorProductOfKirillovReshetikhinTableauxElement
        if isinstance(lst[0], TensorProductOfKirillovReshetikhinTableauxElement):
            if self != lst[0].parent().rigged_configurations():
                raise ValueError("incorrect bijection image")
            return lst[0].to_rigged_configuration()

        from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
        if isinstance(lst[0], TensorProductOfRegularCrystalsElement):
            lst = lst[0]
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinGenericCrystalElement
        if isinstance(lst[0], KirillovReshetikhinGenericCrystalElement):
            KRT = self.tensor_product_of_kirillov_reshetikhin_tableaux()
            krt_elt = KRT(*[x.to_kirillov_reshetikhin_tableau() for x in lst])
            return krt_elt.to_rigged_configuration()

        if isinstance(lst[0], (list, tuple)):
            lst = lst[0]

        if isinstance(lst[0], RiggedPartition):
            lst = [p._clone() for p in lst] # Make a deep copy
        elif isinstance(lst[0], RiggedConfigurationElement):
            lst = [p._clone() for p in lst[0]] # Make a deep copy

        return self.element_class(self, list(lst), **options)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number `p_i^{(a)}` in ``self``.

        This assumes that `\gamma_a = 1` for all `a` and `(\alpha_a \mid
        \alpha_b ) = A_{ab}`.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row length

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 1)
            0
        """
        vac_num = 0
        if "B" in options:
            for tab in options["B"]:
                if len(tab) == self._rc_index[a]:
                    vac_num += min(i, len(tab[0]))
        elif "L" in options:
            L = options["L"]
            if a in L:
                for kvp in L[a].items():
                    vac_num += min(kvp[0], i) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)
        else:
            for dim in self.dims:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)

        if i == float('inf'):
            vac_num -= sum(self._cartan_matrix[a,b] * sum(nu)
                           for b,nu in enumerate(partitions))
        else:
            vac_num -= sum(self._cartan_matrix[a,b] * nu.get_num_cells_to_column(i)
                           for b,nu in enumerate(partitions))

        return vac_num

    def kleber_tree(self):
        r"""
        Return the underlying Kleber tree used to generate all highest
        weight rigged configurations.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',3,1], [[1,1], [2,1]])
            sage: RC.kleber_tree()
            Kleber tree of Cartan type ['A', 3, 1] and B = ((1, 1), (2, 1))
        """
        return KleberTree(self._cartan_type, self.dims)

    @cached_method
    def tensor_product_of_kirillov_reshetikhin_tableaux(self):
        """
        Return the corresponding tensor product of Kirillov-Reshetikhin
        tableaux.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2]])
            sage: RC.tensor_product_of_kirillov_reshetikhin_tableaux()
            Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and factor(s) ((3, 2), (1, 2))
        """
        from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
        return TensorProductOfKirillovReshetikhinTableaux(self._cartan_type, self.dims)

    @cached_method
    def tensor_product_of_kirillov_reshetikhin_crystals(self):
        """
        Return the corresponding tensor product of Kirillov-Reshetikhin
        crystals.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3,1],[2,2]])
            sage: RC.tensor_product_of_kirillov_reshetikhin_crystals()
            Full tensor product of the crystals
            [Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(3,1),
             Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(2,2)]
        """
        return self.tensor_product_of_kirillov_reshetikhin_tableaux().tensor_product_of_kirillov_reshetikhin_crystals()

    def fermionic_formula(self, q=None, only_highest_weight=False, weight=None):
        r"""
        Return the fermionic formula associated to ``self``.

        Given a set of rigged configurations `RC(\lambda, L)`, the fermionic
        formula is defined as:

        .. MATH::

            M(\lambda, L; q) = \sum_{(\nu,J)} q^{cc(\nu, J)}

        where we sum over all (classically highest weight) rigged
        configurations of weight `\lambda` where `cc` is the
        :meth:`cocharge statistic
        <sage.combinat.rigged_configurations.rigged_configuration_element.RiggedConfigurationElement.cc>`.
        This is known to reduce to

        .. MATH::

            M(\lambda, L; q) = \sum_{\nu} q^{cc(\nu)} \prod_{(a,i) \in
            I \times \ZZ} \begin{bmatrix} p_i^{(a)} + m_i^{(a)} \\ m_i^{(a)}
            \end{bmatrix}_q.

        The generating function of `M(\lambda, L; q)` in the weight algebra
        subsumes all fermionic formulas:

        .. MATH::

            M(L; q) = \sum_{\lambda \in P} M(\lambda, L; q) \lambda.

        This is conjecturally equal to the
        :meth:`one dimensional configuration sum
        <sage.combinat.crystals.tensor_product.CrystalOfWords.one_dimensional_configuration_sum>`
        of the corresponding tensor product of Kirillov-Reshetikhin crystals, see [HKOTT2002]_.
        This has been proven in general for type `A_n^{(1)}` [BijectionLRT]_,
        single factors `B^{r,s}` in type `D_n^{(1)}` [OSS2011]_ with the result
        from [Sakamoto13]_, as well as for a tensor product of single columns
        [OSS2003]_, [BijectionDn]_ or a tensor product of single rows [OSS03]_
        for all non-exceptional types.

        INPUT:

        - ``q`` -- the variable `q`
        - ``only_highest_weight`` -- use only the classically highest weight
          rigged configurations
        - ``weight`` -- return the fermionic formula `M(\lambda, L; q)` where
          `\lambda` is the classical weight ``weight``

        REFERENCES:

        .. [OSS2003] Masato Okado, Anne Schilling, and Mark Shimozono.
           Virtual crystals and fermionic formulas of type `D_{n+1}^{(2)}`,
           `A_{2n}^{(2)}`, and `C_n^{(1)}`. Representation Theory. **7** (2003)
           :arxiv:`math.QA/0105017`.

        .. [Sakamoto13] Reiho Sakamoto.
           Rigged configurations and Kashiwara operators.
           (2013) :arxiv:`1302.4562v1`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 2, 1], [[1,1], [1,1]])
            sage: RC.fermionic_formula()
            B[-2*Lambda[1] + 2*Lambda[2]] + (q+1)*B[-Lambda[1]]
             + (q+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
             + B[-2*Lambda[2]] + (q+1)*B[Lambda[2]]
            sage: t = QQ['t'].gen(0)
            sage: RC.fermionic_formula(t)
            B[-2*Lambda[1] + 2*Lambda[2]] + (t+1)*B[-Lambda[1]]
             + (t+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
             + B[-2*Lambda[2]] + (t+1)*B[Lambda[2]]
            sage: La = RC.weight_lattice_realization().classical().fundamental_weights()
            sage: RC.fermionic_formula(weight=La[2])
            q + 1
            sage: RC.fermionic_formula(only_highest_weight=True, weight=La[2])
            q

        Only using the highest weight elements on other types::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3,1], [2,2]])
            sage: RC.fermionic_formula(only_highest_weight=True)
            q*B[Lambda[1] + Lambda[2]] + B[2*Lambda[2] + Lambda[3]]
            sage: RC = RiggedConfigurations(['D', 4, 1], [[3,1], [4,1], [2,1]])
            sage: RC.fermionic_formula(only_highest_weight=True)
            (q^4+q^3+q^2)*B[Lambda[1]] + (q^2+q)*B[Lambda[1] + Lambda[2]]
             + q*B[Lambda[1] + 2*Lambda[3]] + q*B[Lambda[1] + 2*Lambda[4]]
             + B[Lambda[2] + Lambda[3] + Lambda[4]] + (q^3+2*q^2+q)*B[Lambda[3] + Lambda[4]]
            sage: RC = RiggedConfigurations(['E', 6, 1], [[2,2]])
            sage: RC.fermionic_formula(only_highest_weight=True)
            q^2*B[0] + q*B[Lambda[2]] + B[2*Lambda[2]]
            sage: RC = RiggedConfigurations(['B', 3, 1], [[3,1], [2,2]])
            sage: RC.fermionic_formula(only_highest_weight=True) # long time
            q*B[Lambda[1] + Lambda[2] + Lambda[3]] + q^2*B[Lambda[1]
             + Lambda[3]] + (q^2+q)*B[Lambda[2] + Lambda[3]] + B[2*Lambda[2]
             + Lambda[3]] + (q^3+q^2)*B[Lambda[3]]
            sage: RC = RiggedConfigurations(['C', 3, 1], [[3,1], [2,2]])
            sage: RC.fermionic_formula(only_highest_weight=True) # long time
            (q^3+q^2)*B[Lambda[1] + Lambda[2]] + q*B[Lambda[1] + 2*Lambda[2]]
             + (q^2+q)*B[2*Lambda[1] + Lambda[3]] + B[2*Lambda[2] + Lambda[3]]
             + (q^4+q^3+q^2)*B[Lambda[3]]
            sage: RC = RiggedConfigurations(['D', 4, 2], [[3,1], [2,2]])
            sage: RC.fermionic_formula(only_highest_weight=True) # long time
            (q^2+q)*B[Lambda[1] + Lambda[2] + Lambda[3]] + (q^5+2*q^4+q^3)*B[Lambda[1]
             + Lambda[3]] + (q^3+q^2)*B[2*Lambda[1] + Lambda[3]] + (q^4+q^3+q^2)*B[Lambda[2]
             + Lambda[3]] + B[2*Lambda[2] + Lambda[3]] + (q^6+q^5+q^4)*B[Lambda[3]]
            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[1,1],[2,2]])
            sage: RC.fermionic_formula(only_highest_weight=True)
            (q^3+q^2)*B[Lambda[1]] + (q^2+q)*B[Lambda[1] + 2*Lambda[2]]
             + B[Lambda[1] + 4*Lambda[2]] + q*B[3*Lambda[1]] + q*B[4*Lambda[2]]

        TESTS::

            sage: RC = RiggedConfigurations(['A', 2, 1], [[1,1], [1,1]])
            sage: KR = RC.tensor_product_of_kirillov_reshetikhin_crystals()
            sage: RC.fermionic_formula() == KR.one_dimensional_configuration_sum()
            True
            sage: KT = RC.tensor_product_of_kirillov_reshetikhin_tableaux()
            sage: RC.fermionic_formula() == KT.one_dimensional_configuration_sum()
            True
            sage: RC = RiggedConfigurations(['C', 2, 1], [[2,1], [2,1]])
            sage: KR = RC.tensor_product_of_kirillov_reshetikhin_crystals()
            sage: RC.fermionic_formula() == KR.one_dimensional_configuration_sum() # long time
            True
            sage: t = QQ['t'].gen(0)
            sage: RC = RiggedConfigurations(['D', 4, 1], [[1,1], [2,1]])
            sage: KR = RC.tensor_product_of_kirillov_reshetikhin_crystals()
            sage: RC.fermionic_formula(t) == KR.one_dimensional_configuration_sum(t) # long time
            True
        """
        if q is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            q = PolynomialRing(QQ, 'q').gen(0)

        if only_highest_weight:
            L = self.module_generators
        else:
            L = self

        P = q.parent()
        WLR = self.weight_lattice_realization().classical()

        if weight is not None:
            weight = WLR(weight)
            return P.sum(q**x.cc() for x in L if WLR(x.weight()) == weight)

        B = WLR.algebra(P)
        return B.sum(q**x.cc() * B(WLR(x.weight())) for x in L)

    def _test_bijection(self, **options):
        r"""
        Test function to make sure that the bijection between rigged
        configurations and Kirillov-Reshetikhin tableaux is correct.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1],[1,1]])
            sage: RC._test_bijection()
        """
        tester = self._tester(**options)
        rejects = []
        for x in self:
            y = x.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            z = y.to_rigged_configuration()
            if z != x:
                rejects.append((x, z))

        tester.assertEqual(len(rejects), 0, "Bijection is not correct: {}".format(rejects))
        if rejects:
            return rejects

    def tensor(self, *crystals, **options):
        """
        Return the tensor product of ``self`` with ``crystals``.

        If ``crystals`` is a list of rigged configurations of the same
        Cartan type, then this returns a new :class:`RiggedConfigurations`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1],[1,3]])
            sage: RC2 = RiggedConfigurations(['A', 3, 1], [[1,1], [3,3]])
            sage: RC.tensor(RC2, RC2)
            Rigged configurations of type ['A', 3, 1]
             and factor(s) ((2, 1), (1, 3), (1, 1), (3, 3), (1, 1), (3, 3))

            sage: K = crystals.KirillovReshetikhin(['A', 3, 1], 2, 2, model='KR')
            sage: RC.tensor(K)
            Full tensor product of the crystals
             [Rigged configurations of type ['A', 3, 1] and factor(s) ((2, 1), (1, 3)),
              Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and shape (2, 2)]
        """
        ct = self._cartan_type
        if all(isinstance(B, RiggedConfigurations) and B.cartan_type() == ct for B in crystals):
            dims = self.dims
            for B in crystals:
                dims += B.dims
            return RiggedConfigurations(ct, dims)
        return super().tensor(*crystals, **options)

    Element = KRRCSimplyLacedElement


class RCNonSimplyLaced(RiggedConfigurations):
    r"""
    Rigged configurations in non-simply-laced types.

    These are rigged configurations which lift to virtual rigged configurations
    in a simply-laced type.

    For more on rigged configurations, see :class:`RiggedConfigurations`.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, B):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = RiggedConfigurations(CartanType(['A',4,2]), [[2,2]])
            sage: RC2 = RiggedConfigurations(['A',4,2], [(2,2)])
            sage: RC3 = RiggedConfigurations(['BC',2,2], ((2,2),))
            sage: RC2 is RC1, RC3 is RC1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)

        # Standardize B input into a tuple of tuples
        B = tuple(map(tuple, B))
        return super().__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, dims):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,1]])
            sage: TestSuite(RC).run()
            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[2,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['B',3,1], [[3,1],[1,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D',4,2], [[2,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['A',5,2], [[2,1]])
            sage: TestSuite(RC).run() # long time
        """
        self._folded_ct = cartan_type.as_folding()
        RiggedConfigurations.__init__(self, cartan_type, dims)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number `p_i^{(a)}` in ``self``.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row length

        TESTS::

            sage: RC = RiggedConfigurations(['C', 4, 1], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [2], [2], [1]])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 2)
            0
        """
        vac_num = 0
        if "B" in options:
            for tab in options["B"]:
                if len(tab) == self._rc_index[a]:
                    vac_num += min(i, len(tab[0]))
        elif "L" in options:
            L = options["L"]
            if a in L:
                for kvp in L[a].items():
                    vac_num += min(kvp[0], i) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)
        else:
            for dim in self.dims:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)

        if i == float('inf'):
            vac_num -= sum(self._cartan_matrix[a,b] * sum(nu)
                           for b,nu in enumerate(partitions))
        else:
            gamma = self._folded_ct.scaling_factors()
            vac_num -= sum(self._cartan_matrix[a,b]
                           * nu.get_num_cells_to_column(gamma[a+1]*i, gamma[b+1])
                           // gamma[b+1]
                           for b,nu in enumerate(partitions))

        return vac_num

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for this set of rigged configurations.

        Iterate over the highest weight rigged configurations by moving
        through the
        :class:`~sage.combinat.rigged_configurations.kleber_tree.KleberTree`
        and then setting appropriate values of the partitions.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C', 3, 1], [[1,2]])
            sage: for x in RC.module_generators: x
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>

            sage: RC = RiggedConfigurations(['D',4,3], [[1,1]])
            sage: RC.module_generators
            (
            <BLANKLINE>
                 0[ ]0
            (/)  0[ ]0
            <BLANKLINE>
            (/)  0[ ]0
               ,
            )
        """
        module_gens = []
        vec_len = len(self.kleber_tree().root.up_root.to_vector())

        for tree_node in self.kleber_tree():
            shapes = []
            cur = tree_node
            path_lambda = [cur.up_root.to_vector()] # Build the lambda values
            # Note that these are not same lambda as in the paper,
            #   but a less computational version.
            while cur.parent_node is not None:
                path_lambda.insert(0, (cur.parent_node.up_root - cur.up_root).to_vector())
                cur = cur.parent_node

            for a in range(vec_len):
                shapes.append([])
                for i, cur_lambda in enumerate(path_lambda):
                    for j in range(cur_lambda[a]):
                        shapes[-1].insert(0, i)

            # Convert from the virtual rigged configuration
            # As a special case, we do not need to do anything for type `A_{2n}^{(2)}`
            sigma = self._folded_ct.folding_orbit()
            vindex = self.virtual._rc_index
            shapes = [shapes[vindex.index(sigma[a][0])] for a in self._rc_index]
            if self._cartan_type.type() != 'BC':
                gamma = self._folded_ct.scaling_factors()
                for a,shape in enumerate(shapes):
                    for i in range(len(shape)):
                        shape[i] = shape[i] // gamma[self._rc_index[a]]

            # Start with a base to calculate the vacancy numbers
            # Make a copy just to be safe
            base = self.element_class(self, partition_list=shapes[:])

            # Build out the blocks for the partition values
            vac_nums = []
            blocks = []
            L = []

            for partition in base:
                vac_nums.append(partition.vacancy_numbers)
                blocks.append([[]])

                # If the partition is empty, there's nothing to do
                if not partition._list:
                    L.append([[]])
                    continue

                # Setup the first block
                block_len = partition[0]
                for i, rowLen in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != rowLen:
                        blocks[-1].append([])
                        block_len = rowLen

                    blocks[-1][-1].append(partition.vacancy_numbers[i])

                L2 = []
                for block in blocks[-1]:
                    L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(itertools.product(*L2))

            C = itertools.product(*L)
            for cur_blocks in C:
                module_gens.append( self.element_class(self, KT_constructor=[shapes[:],
                                         self._blocks_to_values(cur_blocks[:]), vac_nums[:]]) )

        return tuple(module_gens)

    def kleber_tree(self):
        r"""
        Return the underlying (virtual) Kleber tree used to generate all
        highest weight rigged configurations.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',3,1], [[1,1], [2,1]])
            sage: RC.kleber_tree()
            Virtual Kleber tree of Cartan type ['C', 3, 1] and B = ((1, 1), (2, 1))
        """
        return VirtualKleberTree(self._cartan_type, self.dims)

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
            sage: RC
            Rigged configurations of type ['C', 2, 1] and factor(s) ((1, 2), (1, 1), (2, 1))
            sage: RC.virtual
            Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 2))
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        virtual_dims = []
        for r,s in self.dims:
            for a in sigma[r]:
                virtual_dims.append([a, s*gamma[r]])
        return RiggedConfigurations(self._folded_ct._folding, virtual_dims)

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
            sage: elt = RC(partition_list=[[3],[2]]); elt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ]0
            sage: velt.parent()
            Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 2))
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = len(self.virtual._rc_index)
        # +/- 1 for indexing
        partitions = [None] * n
        for a,rp in enumerate(rc):
            g = gamma[a+1]
            for i in sigma[a+1]:
                partitions[i-1] = RiggedPartition([row_len*g for row_len in rp._list],
                                                  [rig_val*g for rig_val in rp.rigging],
                                                  [vac_num*g for vac_num in rp.vacancy_numbers])
        return self.virtual.element_class(self.virtual, partitions, use_vacancy_numbers=True)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configuration of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
            sage: elt = RC(partition_list=[[3],[2]])
            sage: vrc_elt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(vrc_elt); ret
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: ret == elt
            True
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = len(self._rc_index)
        partitions = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            rp = vrc[sigma[a+1][0] - 1]
            g = gamma[a+1]
            partitions[a] = RiggedPartition([row_len//g for row_len in rp._list],
                                            [rig_val//g for rig_val in rp.rigging],
                                            [vac_val//g for vac_val in rp.vacancy_numbers])
        return self.element_class(self, partitions, use_vacancy_numbers=True)

    def _test_virtual_vacancy_numbers(self, **options):
        """
        Test to make sure that the vacancy numbers obtained from the virtual
        rigged configuration agree with the explicit computation of the
        vacancy numbers done here.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['B', 3, 1], [[2,1]])
            sage: RC._test_virtual_vacancy_numbers()
        """
        tester = self._tester(**options)
        for x in self:
            parts_list = [p._list[:] for p in x]
            elt = self.element_class(self, partition_list=parts_list)
            for i, p in enumerate(elt):
                for j, vac_num in enumerate(p.vacancy_numbers):
                    tester.assertEqual(vac_num, x[i].vacancy_numbers[j],
                      "Incorrect vacancy number: {}\nComputed: {}\nFor: {}".format(
                       x[i].vacancy_numbers[j], vac_num, x))

    Element = KRRCNonSimplyLacedElement


class RCTypeA2Even(RCNonSimplyLaced):
    """
    Rigged configurations for type `A_{2n}^{(2)}`.

    For more on rigged configurations, see :class:`RiggedConfigurations`.

    EXAMPLES::

        sage: RC = RiggedConfigurations(['A',4,2], [[2,1], [1,2]])
        sage: RC.cardinality()
        150
        sage: RC = RiggedConfigurations(['A',2,2], [[1,1]])
        sage: RC.cardinality()
        3
        sage: RC = RiggedConfigurations(['A',2,2], [[1,2],[1,1]])
        sage: TestSuite(RC).run() # long time
        sage: RC = RiggedConfigurations(['A',4,2], [[2,1]])
        sage: TestSuite(RC).run() # long time
    """

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',4,2], [[1,1], [2,2]])
            sage: RC.cardinality()
            250
        """
        return self.tensor_product_of_kirillov_reshetikhin_tableaux().cardinality()

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',4,2], [[1,2],[1,1],[2,1]])
            sage: RC
            Rigged configurations of type ['BC', 2, 2] and factor(s) ((1, 2), (1, 1), (2, 1))
            sage: RC.virtual
            Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 1), (2, 1))
        """
        sigma = self._folded_ct.folding_orbit()
        n = len(sigma) - 1
        virtual_dims = []
        for r,s in self.dims:
            if r == n:
                virtual_dims.extend([[n, s], [n, s]])
            else:
                for a in sigma[r]:
                    virtual_dims.append([a, s])
        return RiggedConfigurations(self._folded_ct._folding, virtual_dims)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number `p_i^{(a)}` in ``self``.

        This is a special implementation for type `A_{2n}^{(2)}`.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row length

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 2], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [2]])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 2)
            0
        """
        vac_num = 0
        if "B" in options:
            for tab in options["B"]:
                if len(tab) == self._rc_index[a]:
                    vac_num += min(i, len(tab[0]))
        elif "L" in options:
            L = options["L"]
            if a in L:
                for kvp in L[a].items():
                    vac_num += min(kvp[0], i) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)
        else:
            for dim in self.dims:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)

        gamma = self._folded_ct.scaling_factors()
        if i == float('inf'):
            vac_num -= sum(self._cartan_matrix[a,b] * sum(nu) // gamma[b+1]
                           for b, nu in enumerate(partitions))
        else:
            vac_num -= sum(self._cartan_matrix[a,b] * nu.get_num_cells_to_column(i) // gamma[b+1]
                           for b, nu in enumerate(partitions))

        return vac_num

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',4,2], [[2,2]])
            sage: elt = RC(partition_list=[[1],[1]]); elt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            2[ ]2
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: velt.parent()
            Rigged configurations of type ['A', 3, 1] and factor(s) ((2, 2), (2, 2))
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = len(self.virtual._rc_index)
        partitions = [None] * n
        for a,rp in enumerate(rc):
            g = gamma[a+1]
            for i in sigma[a+1]:
                partitions[i-1] = RiggedPartition(list(rp._list),
                                                  [rig_val*g for rig_val in rp.rigging],
                                                  [vac_num*g for vac_num in rp.vacancy_numbers])
        return self.virtual.element_class(self.virtual, partitions, use_vacancy_numbers=True)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configuration of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',4,2], [[2,2]])
            sage: elt = RC(partition_list=[[1],[1]])
            sage: velt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(velt); ret
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: ret == elt
            True
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = len(self._rc_index)
        partitions = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            rp = vrc[sigma[a+1][0] - 1]
            g = gamma[a+1]
            partitions[a] = RiggedPartition(list(rp._list),
                                            [rig_val//g for rig_val in rp.rigging],
                                            [vac_val//g for vac_val in rp.vacancy_numbers])
        return self.element_class(self, partitions, use_vacancy_numbers=True)


class RCTypeA2Dual(RCTypeA2Even):
    r"""
    Rigged configurations of type `A_{2n}^{(2)\dagger}`.

    For more on rigged configurations, see :class:`RiggedConfigurations`.

    EXAMPLES::

        sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[1,2],[1,1],[2,1]])
        sage: RC
        Rigged configurations of type ['BC', 2, 2]^* and factor(s) ((1, 2), (1, 1), (2, 1))
        sage: RC.cardinality()
        750
        sage: RC.virtual
        Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 1), (2, 1))
        sage: RC = RiggedConfigurations(CartanType(['A',2,2]).dual(), [[1,1]])
        sage: RC.cardinality()
        3
        sage: RC = RiggedConfigurations(CartanType(['A',2,2]).dual(), [[1,2],[1,1]])
        sage: TestSuite(RC).run() # long time
        sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[2,1]])
        sage: TestSuite(RC).run() # long time
    """

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number `p_i^{(a)}` in ``self``. A special case
        is needed for the `n`-th partition for type `A_{2n}^{(2)\dagger}`.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row length

        TESTS::

            sage: RC = RiggedConfigurations(CartanType(['A', 6, 2]).dual(), [[2,1]])
            sage: elt = RC(partition_list=[[1], [2], [2]])
            sage: RC._calc_vacancy_number(elt.nu(), 0, 1)
            -1
        """
        if a != len(self._rc_index) - 1:
            return RCTypeA2Even._calc_vacancy_number(self, partitions, a, i, **options)

        vac_num = 0
        if "B" in options:
            for tab in options["B"]:
                if len(tab) == self._rc_index[a]:
                    vac_num += min(i, len(tab[0]))
        elif "L" in options:
            L = options["L"]
            if a in L:
                for kvp in L[a].items():
                    vac_num += min(kvp[0], i) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)
        else:
            for dim in self.dims:
                if dim[0] == self._rc_index[a]:
                    vac_num += min(dim[1], i)

        if i == float('inf'):
            vac_num -= sum(self._cartan_matrix[a,b] * sum(nu) / 2
                           for b,nu in enumerate(partitions))
        else:
            vac_num -= sum(self._cartan_matrix[a,b] * nu.get_num_cells_to_column(i) / 2
                           for b,nu in enumerate(partitions))

        return vac_num

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for rigged configurations of type
        `A_{2n}^{(2)\dagger}`.

        Iterate over the highest weight rigged configurations by moving
        through the
        :class:`~sage.combinat.rigged_configurations.kleber_tree.KleberTree`
        and then setting appropriate values of the partitions. This also
        skips rigged configurations where `P_i^{(n)} < 1` when `i` is odd.

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A', 4, 2]).dual(), [[1,1]])
            sage: for x in RC.module_generators: x
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        module_gens = []
        # This is for the non-simply-laced types
        vec_len = len(self.kleber_tree().root.up_root.to_vector())

        for tree_node in self.kleber_tree():
            shapes = []
            cur = tree_node
            path_lambda = [cur.up_root.to_vector()] # Build the lambda values
            # Note that these are not same lambda as in the paper,
            #   but a less computational version.
            while cur.parent_node is not None:
                path_lambda.insert(0, (cur.parent_node.up_root - cur.up_root).to_vector())
                cur = cur.parent_node

            for a in range(vec_len):
                shapes.append([])
                for i, cur_lambda in enumerate(path_lambda):
                    for j in range(cur_lambda[a]):
                        shapes[-1].insert(0, i)

            # We are not simply-laced, so convert from the virtual rigged configuration
            sigma = self._folded_ct.folding_orbit()
            vindex = self.virtual._rc_index
            shapes = [shapes[vindex.index(sigma[a][0])] for a in self._rc_index]
            # Nothing more to do since gamma[i] == 1 for all i >= 1

            # Start with a base to calculate the vacancy numbers
            # Make a copy just to be safe
            base = self.element_class(self, partition_list=shapes[:])

            # Check the special condition of odd rows in the n-th partition
            invalid_RC = False
            for i in range(len(base[-1]._list)):
                if base[-1]._list[i] % 2 == 1 and base[-1].vacancy_numbers[i] < 1:
                    invalid_RC = True
                    break
            # If it is invalid, skip it
            if invalid_RC:
                continue

            # Build out the blocks for the partition values
            vac_nums = []
            blocks = []
            L = []

            for partition in base[:-1]:
                vac_nums.append(partition.vacancy_numbers)
                blocks.append([[]])

                # If the partition is empty, there's nothing to do
                if len(partition) <= 0:
                    L.append([[]])
                    continue

                # Setup the first block
                block_len = partition[0]
                for i, row_len in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != row_len:
                        blocks[-1].append([])
                        block_len = row_len

                    blocks[-1][-1].append(partition.vacancy_numbers[i])

                L2 = []
                for block in blocks[-1]:
                    L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(itertools.product(*L2))

            # Special case for the final tableau
            partition = base[-1]
            vac_nums.append(partition.vacancy_numbers)

            # If the partition is empty, there's nothing to do
            if len(partition) <= 0:
                L.append([[]])
            else:
                # Setup the first block
                block_len = partition[0]
                blocks = [[]]
                odd_block = []
                for i, row_len in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != row_len:
                        blocks.append([])
                        odd_block.append(block_len % 2 == 1)
                        block_len = row_len

                    blocks[-1].append(partition.vacancy_numbers[i])
                odd_block.append(block_len % 2 == 1)

                L2 = []
                for i, block in enumerate(blocks):
                    if odd_block[i]:
                        L2.append(IterableFunctionCall(self._block_iterator_n_odd, block))
                    else:
                        L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(itertools.product(*L2))

            C = itertools.product(*L)
            for curBlocks in C:
                module_gens.append( self.element_class(self, KT_constructor=[shapes[:],
                                         self._blocks_to_values(curBlocks[:]), vac_nums[:]]) )

        return tuple(module_gens)

    def _block_iterator_n_odd(self, container):
        r"""
        Iterate over all possible riggings for a block of odd length in the
        `n`-th rigged partition for type `A_{2n}^{(2)\dagger}`.

        Helper iterator which iterates over all possible partitions of
        `\frac{2k+1}{2}` sizes contained within the container.

        INPUT:

        - ``container`` -- list the widths of the rows of the container

        TESTS::

            sage: RC = RiggedConfigurations(CartanType(['A', 4, 2]).dual(), [[2, 2]])
            sage: for x in RC._block_iterator_n_odd([]): x
            []
            sage: for x in RC._block_iterator_n_odd([2,2]): x
            [1/2, 1/2]
            [3/2, 1/2]
            [3/2, 3/2]
        """
        if len(container) == 0:
            yield []
            return

        pos = 0
        length = len(container)
        ret_part = [-1] * length
        while pos >= 0:
            ret_part[pos] += 2

            if ret_part[pos] > container[pos]*2 or (pos != 0 and ret_part[pos] > ret_part[pos - 1]):
                ret_part[pos] = -1
                pos -= 1
            else:
                pos += 1

            if pos == length:
                yield [QQ(n) / QQ(2) for n in ret_part]
                pos -= 1

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[2,2]])
            sage: elt = RC(partition_list=[[1],[1]]); elt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            2[ ]2
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: velt.parent()
            Rigged configurations of type ['A', 3, 1] and factor(s) ((2, 2), (2, 2))
        """
        gammatilde = list(self._folded_ct.scaling_factors())
        gammatilde[-1] = 2
        sigma = self._folded_ct.folding_orbit()
        n = len(self.virtual._rc_index)
        partitions = [None] * n
        for a,rp in enumerate(rc):
            g = gammatilde[a+1]
            for i in sigma[a+1]:
                partitions[i-1] = RiggedPartition(list(rp._list),
                                                  [rig_val*g for rig_val in rp.rigging])
        return self.virtual.element_class(self.virtual, partitions)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configuration of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[2,2]])
            sage: elt = RC(partition_list=[[1],[1]])
            sage: velt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(velt); ret
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: ret == elt
            True
        """
        gammatilde = list(self._folded_ct.scaling_factors())
        gammatilde[-1] = QQ(2)
        sigma = self._folded_ct.folding_orbit()
        n = len(self._rc_index)
        partitions = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            rp = vrc[sigma[a+1][0] - 1]
            g = gammatilde[a+1]
            partitions[a] = RiggedPartition(list(rp._list),
                                            [rig_val/g for rig_val in rp.rigging])
        return self.element_class(self, partitions)

    Element = KRRCTypeA2DualElement
