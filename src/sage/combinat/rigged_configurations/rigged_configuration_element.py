# sage.doctest: needs sage.combinat sage.modules
r"""
Rigged configuration elements

A rigged configuration element is a sequence of
:class:`~sage.combinat.rigged_configurations.rigged_partition.RiggedPartition`
objects.

AUTHORS:

- Travis Scrimshaw (2010-09-26): initial version
- Travis Scrimshaw (2012-10-25): added virtual rigged configurations
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
from sage.misc.cachefunc import cached_method
from sage.structure.list_clone import ClonableArray
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition, RiggedPartitionTypeB


####################################################
#  Base classes for rigged configuration elements  #
####################################################

class RiggedConfigurationElement(ClonableArray):
    """
    A rigged configuration for simply-laced types.

    For more information on rigged configurations, see
    :class:`RiggedConfigurations`. For rigged configurations for
    non-simply-laced types, use :class:`RCNonSimplyLacedElement`.

    Typically to create a specific rigged configuration, the user will pass in
    the optional argument ``partition_list`` and if the user wants to specify
    the rigging values, give the optional argument ``rigging_list`` as well.
    If ``rigging_list`` is not passed, the rigging values are set to the
    corresponding vacancy numbers.

    INPUT:

    - ``parent`` -- the parent of this element

    - ``rigged_partitions`` -- list of rigged partitions

    There are two optional arguments to explicitly construct a rigged
    configuration. The first is ``partition_list`` which gives a list of
    partitions, and the second is ``rigging_list`` which is a list of
    corresponding lists of riggings. If only partition_list is specified,
    then it sets the rigging equal to the calculated vacancy numbers.

    If we are constructing a rigged configuration from a rigged configuration
    (say of another type) and we don't want to recompute the vacancy numbers,
    we can use the ``use_vacancy_numbers`` to avoid the recomputation.

    EXAMPLES:

    Type `A_n^{(1)}` examples::

        sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
        sage: RC(partition_list=[[2], [2, 2], [2], [2]])
        <BLANKLINE>
        0[ ][ ]0
        <BLANKLINE>
        -2[ ][ ]-2
        -2[ ][ ]-2
        <BLANKLINE>
        2[ ][ ]2
        <BLANKLINE>
        -2[ ][ ]-2
        <BLANKLINE>

        sage: RC = RiggedConfigurations(['A', 4, 1], [[1, 1], [1, 1]])
        sage: RC(partition_list=[[], [], [], []])
        <BLANKLINE>
        (/)
        <BLANKLINE>
        (/)
        <BLANKLINE>
        (/)
        <BLANKLINE>
        (/)
        <BLANKLINE>

    Type `D_n^{(1)}` examples::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
        sage: RC(partition_list=[[3], [3,2], [4], [3]])
        <BLANKLINE>
        -1[ ][ ][ ]-1
        <BLANKLINE>
        1[ ][ ][ ]1
        0[ ][ ]0
        <BLANKLINE>
        -3[ ][ ][ ][ ]-3
        <BLANKLINE>
        -1[ ][ ][ ]-1
        <BLANKLINE>

        sage: RC = RiggedConfigurations(['D', 4, 1], [[1, 1], [2, 1]])
        sage: RC(partition_list=[[1], [1,1], [1], [1]])
        <BLANKLINE>
        1[ ]1
        <BLANKLINE>
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>
        sage: elt = RC(partition_list=[[1], [1,1], [1], [1]], rigging_list=[[0], [0,0], [0], [0]]); elt
        <BLANKLINE>
        1[ ]0
        <BLANKLINE>
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>

        sage: from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition
        sage: RC2 = RiggedConfigurations(['D', 5, 1], [[2, 1], [3, 1]])
        sage: l = [RiggedPartition()] + list(elt)
        sage: ascii_art(RC2(*l))
        (/)  1[ ]0  0[ ]0  0[ ]0  0[ ]0
                    0[ ]0
        sage: ascii_art(RC2(*l, use_vacancy_numbers=True))
        (/)  1[ ]0  0[ ]0  0[ ]0  0[ ]0
                    0[ ]0
    """

    def __init__(self, parent, rigged_partitions=[], **options):
        r"""
        Construct a rigged configuration element.

        .. WARNING::

            This changes the vacancy numbers of the rigged partitions, so
            if the rigged partitions comes from another rigged configuration,
            a deep copy should be made before being passed here. We do not
            make a deep copy here because the crystal operators generate
            their own rigged partitions. See :issue:`17054`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: RC(partition_list=[[], [], [], []])
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: RC(partition_list=[[1], [1], [], []])
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: elt = RC(partition_list=[[1], [1], [], []], rigging_list=[[-1], [0], [], []]); elt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: TestSuite(elt).run()
        """
        n = options.get('n', parent._cartan_type.rank())
        if "partition_list" in options:
            data = options["partition_list"]
            if len(data) == 0:
                # Create a size n array of empty rigged tableau since no tableau
                #   were given
                nu = [RiggedPartition() for _ in range(n)]
            else:
                if len(data) != n:  # otherwise n should be equal to the number of tableaux
                    raise ValueError("incorrect number of partitions")

                nu = []
                if "rigging_list" in options:
                    rigging_data = options["rigging_list"]

                    if len(rigging_data) != n:
                        raise ValueError("incorrect number of riggings")

                    for i in range(n):
                        nu.append(RiggedPartition(tuple(data[i]),
                                                  list(rigging_data[i])))
                else:
                    for partition_data in data:
                        nu.append(RiggedPartition(tuple(partition_data)))
        elif n == len(rigged_partitions) and isinstance(rigged_partitions[0], RiggedPartition):
            # The isinstance check is to make sure we are not in the n == 1 special case because
            #   Parent's __call__ always passes at least 1 argument to the element constructor

            if options.get('use_vacancy_numbers', False):
                ClonableArray.__init__(self, parent, rigged_partitions)
                return
            nu = rigged_partitions
        else:
            # Otherwise we did not receive any info, create a size n array of
            #   empty rigged partitions
            nu = []
            for i in range(n):
                nu.append(RiggedPartition())
            # raise ValueError("Invalid input")
            # raise ValueError("Incorrect number of rigged partitions")

        # Set the vacancy numbers
        for a, partition in enumerate(nu):
            # If the partition is empty, there's nothing to do
            if len(partition) <= 0:
                continue

            # Setup the first block
            block_len = partition[0]
            vac_num = parent._calc_vacancy_number(nu, a, block_len)

            for i, row_len in enumerate(partition):
                # If we've gone to a different sized block, then update the
                #   values which change when moving to a new block size
                if block_len != row_len:
                    vac_num = parent._calc_vacancy_number(nu, a, row_len)
                    block_len = row_len

                partition.vacancy_numbers[i] = vac_num
                if partition.rigging[i] is None:
                    partition.rigging[i] = partition.vacancy_numbers[i]

        ClonableArray.__init__(self, parent, nu)

    def check(self):
        """
        Check the rigged configuration is properly defined.

        There is nothing to check here.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A', 4])
            sage: b = RC.module_generators[0].f_string([1,2,1,1,2,4,2,3,3,2])
            sage: b.check()
        """
        pass

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[2], [3,1], [3], [3]]); elt
            <BLANKLINE>
            -1[ ][ ]-1
            <BLANKLINE>
            2[ ][ ][ ]2
            0[ ]0
            <BLANKLINE>
            -2[ ][ ][ ]-2
            <BLANKLINE>
            -2[ ][ ][ ]-2
            <BLANKLINE>
            sage: RC.options(display='horizontal')
            sage: elt
            -1[ ][ ]-1   2[ ][ ][ ]2   -2[ ][ ][ ]-2   -2[ ][ ][ ]-2
                         0[ ]0
            sage: RC.options._reset()
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_vertical(self):
        """
        Return the string representation of ``self`` vertically.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: print(RC(partition_list=[[2], [3,1], [3], [3]])._repr_vertical())
            <BLANKLINE>
            -1[ ][ ]-1
            <BLANKLINE>
            2[ ][ ][ ]2
            0[ ]0
            <BLANKLINE>
            -2[ ][ ][ ]-2
            <BLANKLINE>
            -2[ ][ ][ ]-2
            <BLANKLINE>
            sage: print(RC(partition_list=[[],[],[],[]])._repr_vertical())
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        ret_str = ""
        for tableau in self:
            ret_str += "\n" + repr(tableau)
        return ret_str

    def _repr_horizontal(self):
        """
        Return the string representation of ``self`` horizontally.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: print(RC(partition_list=[[2], [3,1], [3], [3]])._repr_horizontal())
            -1[ ][ ]-1   2[ ][ ][ ]2   -2[ ][ ][ ]-2   -2[ ][ ][ ]-2
                         0[ ]0
            sage: print(RC(partition_list=[[],[],[],[]])._repr_horizontal())
            (/)   (/)   (/)   (/)
        """
        tab_str = [repr(x).splitlines() for x in self]
        height = max(len(t) for t in tab_str)
        widths = [max(len(x) for x in t) for t in tab_str]
        ret_str = ''
        for i in range(height):
            if i != 0:
                ret_str += '\n'
            for j,t in enumerate(tab_str):
                if j != 0:
                    ret_str += '   '
                if i < len(t):
                    ret_str += t[i] + ' ' * (widths[j]-len(t[i]))
                else:
                    ret_str += ' ' * widths[j]
        return ret_str

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: latex(RC(partition_list=[[2], [3,1], [3], [3]]))
            {
            \begin{array}[t]{r|c|c|l}
            \cline{2-3} -1 &\phantom{|}&\phantom{|}& -1 \\
             \cline{2-3}
            \end{array}
            }
            \quad
             {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} 2 &\phantom{|}&\phantom{|}&\phantom{|}& 2 \\
             \cline{2-4} 0 &\phantom{|}& \multicolumn{3 }{l}{ 0 } \\
             \cline{2-2}
            \end{array}
            }
            \quad
             {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} -2 &\phantom{|}&\phantom{|}&\phantom{|}& -2 \\
             \cline{2-4}
            \end{array}
            }
            \quad
             {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} -2 &\phantom{|}&\phantom{|}&\phantom{|}& -2 \\
             \cline{2-4}
            \end{array}
            }
            sage: latex(RC(partition_list=[[],[],[],[]]))
            {\emptyset}
            \quad
            {\emptyset}
            \quad
            {\emptyset}
            \quad
            {\emptyset}
        """
        ret_string = self[0]._latex_()

        for partition in self[1:]:
            ret_string += "\n\\quad\n" + partition._latex_()

        return ret_string

    def _ascii_art_(self):
        """
        Return an ASCII art representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: ascii_art(RC(partition_list=[[2], [3,1], [3], [3]]))
            -1[ ][ ]-1  2[ ][ ][ ]2  -2[ ][ ][ ]-2  -2[ ][ ][ ]-2
                        0[ ]0
            sage: ascii_art(RC(partition_list=[[],[],[],[]]))
            (/)  (/)  (/)  (/)
            sage: RC = RiggedConfigurations(['D', 7, 1], [[3,3],[5,2],[4,3],[2,3],[4,4],[3,1],[1,4],[2,2]])
            sage: elt = RC(partition_list=[[2],[3,2,1],[2,2,1,1],[2,2,1,1,1,1],[3,2,1,1,1,1],[2,1,1],[2,2]],
            ....:          rigging_list=[[2],[1,0,0],[4,1,2,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0],[0,0]])
            sage: ascii_art(elt)
            3[ ][ ]2  1[ ][ ][ ]1  4[ ][ ]4  2[ ][ ]1  0[ ][ ][ ]0  0[ ][ ]0  0[ ][ ]0
                      2[ ][ ]0     4[ ][ ]1  2[ ][ ]0  2[ ][ ]1     0[ ]0     0[ ][ ]0
                      1[ ]0        3[ ]2     0[ ]0     0[ ]0        0[ ]0
                                   3[ ]1     0[ ]0     0[ ]0
                                             0[ ]0     0[ ]0
                                             0[ ]0     0[ ]0
            sage: Partitions.options(convention='French')
            sage: ascii_art(elt)
                                             0[ ]0     0[ ]0
                                             0[ ]0     0[ ]0
                                   3[ ]1     0[ ]0     0[ ]0
                      1[ ]0        3[ ]2     0[ ]0     0[ ]0        0[ ]0
                      2[ ][ ]0     4[ ][ ]1  2[ ][ ]0  2[ ][ ]1     0[ ]0     0[ ][ ]0
            3[ ][ ]2  1[ ][ ][ ]1  4[ ][ ]4  2[ ][ ]1  0[ ][ ][ ]0  0[ ][ ]0  0[ ][ ]0
            sage: Partitions.options._reset()
        """
        from sage.combinat.partition import Partitions
        if Partitions.options.convention == "French":
            baseline = lambda s: 0
        else:
            baseline = len
        from sage.typeset.ascii_art import AsciiArt
        s = repr(self[0]).splitlines()
        ret = AsciiArt(s, baseline=baseline(s))
        for tableau in self[1:]:
            s = repr(tableau).splitlines()
            ret += AsciiArt(["  "], baseline=baseline(s)) + AsciiArt(s, baseline=baseline(s))
        return ret

    def nu(self):
        r"""
        Return the list `\nu` of rigged partitions of this rigged
        configuration element.

        OUTPUT: the `\nu` array as a list

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [2,2], [2], [2]]).nu()
            [0[ ][ ]0
             , -2[ ][ ]-2
             -2[ ][ ]-2
             , 2[ ][ ]2
             , -2[ ][ ]-2
             ]
        """
        return list(self)

    # TODO: Change e/f to work for all types
    def e(self, a):
        r"""
        Return the action of the crystal operator `e_a` on ``self``.

        This implements the method defined in [CrysStructSchilling06]_ which
        finds the value `k` which is  the length of the string with the
        smallest negative rigging of smallest length. Then it removes a box
        from a string of length `k` in the `a`-th rigged partition, keeping all
        colabels fixed and increasing the new label by one. If no such string
        exists, then `e_a` is undefined.

        This method can also be used when the underlying Cartan matrix is a
        Borcherds-Cartan matrix.  In this case, then method of [SS2018]_ is
        used, where the new label is increased by half of the `a`-th diagonal
        entry of the underlying Borcherds-Cartan matrix.  This method will also
        return ``None`` if `a` is imaginary and the smallest rigging in the
        `a`-th rigged partition is not exactly half of the `a`-th diagonal entry
        of the Borcherds-Cartan matrix.

        INPUT:

        - ``a`` -- the index of the partition to remove a box

        OUTPUT: the resulting rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: elt = RC(partition_list=[[1], [1], [1], [1]])
            sage: elt.e(3)
            sage: elt.e(1)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>

            sage: A = CartanMatrix([[-2,-1],[-1,-2]], borcherds=True)
            sage: RC = crystals.infinity.RiggedConfigurations(A)
            sage: nu0 = RC(partition_list=[[],[]])
            sage: nu = nu0.f_string([1,0,0,0])
            sage: ascii_art(nu.e(0))
            5[ ]3  4[ ]3
            5[ ]1
        """
        if a not in self.parent()._rc_index_inverse:
            raise ValueError("{} is not in the index set".format(a))
        a = self.parent()._rc_index_inverse[a]
        M = self.parent()._cartan_matrix

        new_list = self[a][:]
        new_vac_nums = self[a].vacancy_numbers[:]
        new_rigging = self[a].rigging[:]

        # Separate out one of the Borcherds cases
        if M[a,a] != 2:
            k = None
            set_vac_num = True
            if new_rigging[-1] != -M[a,a] // 2:
                return None
            new_list.pop()
            new_vac_nums.pop()
            new_rigging.pop()
        else:
            # Find k and perform e_a
            k = None
            num_rows = len(new_list)
            cur_rigging = -1
            rigging_index = None
            for i in range(num_rows):
                if new_rigging[i] <= cur_rigging:
                    cur_rigging = new_rigging[i]
                    rigging_index = i

            # If we've not found a valid k
            if rigging_index is None:
                return None

            # Note that because the riggings are weakly decreasing, we will always
            #   remove the last box on of a block
            k = new_list[rigging_index]
            set_vac_num = False
            if k == 1:
                new_list.pop()
                new_vac_nums.pop()
                new_rigging.pop()
            else:
                new_list[rigging_index] -= 1
                cur_rigging += M[a,a] // 2
                # Properly sort the riggings
                j = rigging_index + 1
                # Update the vacancy number if the row lengths are the same
                if j < num_rows and new_list[j] == new_list[rigging_index]:
                    new_vac_nums[rigging_index] = new_vac_nums[j]
                    set_vac_num = True
                while j < num_rows and new_list[j] == new_list[rigging_index] \
                  and new_rigging[j] > cur_rigging:
                    new_rigging[j-1] = new_rigging[j] # Shuffle it along
                    j += 1
                new_rigging[j-1] = cur_rigging

        new_partitions = []
        for b in range(len(self)):
            if b != a:
                new_partitions.append(self._generate_partition_e(a, b, k))
            else:
                # Update the vacancy numbers and the rigging
                for i in range(len(new_vac_nums)):
                    if k is not None and new_list[i] < k:
                        break

                    new_vac_nums[i] += M[a,b]
                    new_rigging[i] += M[a,b]

                if k != 1 and not set_vac_num: # If we did not remove a row nor found another row of length k-1
                    new_vac_nums[rigging_index] += 2

                new_partitions.append(RiggedPartition(new_list, new_rigging, new_vac_nums))

        ret_RC = self.__class__(self.parent(), new_partitions, use_vacancy_numbers=True)
        nu = ret_RC.nu()
        if k != 1 and not set_vac_num: # If we did not remove a row nor found another row of length k-1
            # Update that row's vacancy number
            ret_RC[a].vacancy_numbers[rigging_index] = \
              self.parent()._calc_vacancy_number(nu, a, nu[a][rigging_index])
        return ret_RC

    def _generate_partition_e(self, a, b, k):
        r"""
        Generate a new partition for a given value of `a` by updating the
        vacancy numbers and preserving co-labels for the map `e_a`.

        INPUT:

        - ``a`` -- the index of the partition we operated on
        - ``b`` -- the index of the partition to generate
        - ``k`` -- the length of the string with the smallest negative
          rigging of smallest length

        OUTPUT: the constructed rigged partition

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: RC(partition_list=[[1], [1], [1], [1]])._generate_partition_e(1, 2, 1)
            -1[ ]-1
            <BLANKLINE>
        """
        # Check to make sure we will do something
        if not self.parent()._cartan_matrix[a,b]:
            return self[b]

        new_list = self[b]._list
        new_vac_nums = self[b].vacancy_numbers[:]
        new_rigging = self[b].rigging[:]

        # Update the vacancy numbers and the rigging
        value = self.parent()._cartan_matrix[b,a]
        for i in range(len(new_vac_nums)):
            if k is not None and new_list[i] < k:
                break

            new_vac_nums[i] += value
            new_rigging[i] += value

        return RiggedPartition(new_list, new_rigging, new_vac_nums)

    def f(self, a):
        r"""
        Return the action of the crystal operator `f_a` on ``self``.

        This implements the method defined in [CrysStructSchilling06]_ which
        finds the value `k` which is the length of the string with the
        smallest nonpositive rigging of largest length. Then it adds a box from
        a string of length `k` in the `a`-th rigged partition, keeping all
        colabels fixed and decreasing the new label by one. If no such string
        exists, then it adds a new string of length 1 with label `-1`. However
        we need to modify the definition to work for `B(\infty)` by removing
        the condition that the resulting rigged configuration is valid.

        This method can also be used when the underlying Cartan matrix is a
        Borcherds-Cartan matrix.  In this case, then method of [SS2018]_ is
        used, where the new label is decreased by half of the `a`-th diagonal
        entry of the underlying Borcherds-Cartan matrix.

        INPUT:

        - ``a`` -- the index of the partition to add a box

        OUTPUT: the resulting rigged configuration element

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A', 3])
            sage: nu0 = RC.module_generators[0]
            sage: nu0.f(2)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            -2[ ]-1
            <BLANKLINE>
            (/)
            <BLANKLINE>

            sage: A = CartanMatrix([[-2,-1],[-1,-2]], borcherds=True)
            sage: RC = crystals.infinity.RiggedConfigurations(A)
            sage: nu0 = RC(partition_list=[[],[]])
            sage: nu = nu0.f_string([1,0,0,0])
            sage: ascii_art(nu.f(0))
            9[ ]7  6[ ]5
            9[ ]5
            9[ ]3
            9[ ]1
        """
        if a not in self.parent()._rc_index_inverse:
            raise ValueError("{} is not in the index set".format(a))
        a = self.parent()._rc_index_inverse[a]
        M = self.parent()._cartan_matrix

        new_list = self[a][:]
        new_vac_nums = self[a].vacancy_numbers[:]
        new_rigging = self[a].rigging[:]

        # Find k and perform f_a
        k = None
        add_index = -1 # Index where we will add our row too
        rigging_index = None # Index which we will pull the rigging from
        cur_rigging = ZZ.zero()
        num_rows = len(new_list)
        for i in reversed(range(num_rows)):
            # If we need to increment a row, look for when we change rows for
            #   the correct index.
            if add_index is None and new_list[i] != new_list[rigging_index]:
                add_index = i+1

            if new_rigging[i] <= cur_rigging:
                cur_rigging = new_rigging[i]
                k = new_list[i]
                rigging_index = i
                add_index = None

        # If we've not found a valid k
        if k is None:
            new_list.append(1)
            new_rigging.append(-M[a,a] // 2)
            new_vac_nums.append(None)
            k = 0
            add_index = num_rows
            num_rows += 1 # We've added a row
        else:
            if add_index is None: # We are adding to the first row in the list
                add_index = 0
            new_list[add_index] += 1
            new_rigging.insert(add_index, new_rigging[rigging_index] - M[a,a] // 2)
            new_vac_nums.insert(add_index, None)
            new_rigging.pop(rigging_index + 1) # add 1 for the insertion
            new_vac_nums.pop(rigging_index + 1)

        new_partitions = []
        for b in range(len(self)):
            if b != a:
                new_partitions.append(self._generate_partition_f(a, b, k))
            else:
                # Update the vacancy numbers and the rigging
                for i in range(num_rows):
                    if new_list[i] <= k:
                        break

                    if i != add_index:
                        new_vac_nums[i] -= M[a,b]
                        new_rigging[i] -= M[a,b]

                new_partitions.append(RiggedPartition(new_list, new_rigging, new_vac_nums))

        new_partitions[a].vacancy_numbers[add_index] = \
          self.parent()._calc_vacancy_number(new_partitions, a,
                                             new_partitions[a][add_index])

        # Note that we do not need to sort the rigging since if there was a
        #   smaller rigging in a larger row, then `k` would be larger.
        return self.__class__(self.parent(), new_partitions, use_vacancy_numbers=True)

    def _generate_partition_f(self, a, b, k):
        r"""
        Generate a new partition for a given value of `a` by updating the
        vacancy numbers and preserving co-labels for the map `f_a`.

        INPUT:

        - ``a`` -- the index of the partition we operated on
        - ``b`` -- the index of the partition to generate
        - ``k`` -- the length of the string with smallest nonpositive rigging
          of largest length

        OUTPUT: the constructed rigged partition

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: RC(partition_list=[[1], [1], [1], [1]])._generate_partition_f(1, 2, 1)
            0[ ]0
            <BLANKLINE>
        """
        # Check to make sure we will do something
        if not self.parent()._cartan_matrix[a,b]:
            return self[b]

        new_list = self[b]._list
        new_vac_nums = self[b].vacancy_numbers[:]
        new_rigging = self[b].rigging[:]

        # Update the vacancy numbers and the rigging
        value = self.parent()._cartan_matrix[b,a]
        for i in range(len(new_vac_nums)):
            if new_list[i] <= k:
                break

            new_vac_nums[i] -= value
            new_rigging[i] -= value

        return RiggedPartition(new_list, new_rigging, new_vac_nums)

    def epsilon(self, a):
        r"""
        Return `\varepsilon_a` of ``self``.

        Let `x_{\ell}` be the smallest string of `\nu^{(a)}` or `0` if
        `\nu^{(a)} = \emptyset`, then we have
        `\varepsilon_a = -\min(0, x_{\ell})`.

        EXAMPLES::

            sage: La = RootSystem(['B',2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1]+La[2])
            sage: I = RC.index_set()
            sage: matrix([[rc.epsilon(i) for i in I] for rc in RC[:4]])
            [0 0]
            [1 0]
            [0 1]
            [0 2]
        """
        a = self.parent()._rc_index_inverse[a]
        if not self[a]:
            return ZZ.zero()
        return Integer(-min(0, min(self[a].rigging)))

    def phi(self, a):
        r"""
        Return `\varphi_a` of ``self``.

        Let `x_{\ell}` be the smallest string of `\nu^{(a)}` or `0` if
        `\nu^{(a)} = \emptyset`, then we have
        `\varepsilon_a = p_{\infty}^{(a)} - \min(0, x_{\ell})`.

        EXAMPLES::

            sage: La = RootSystem(['B',2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1]+La[2])
            sage: I = RC.index_set()
            sage: matrix([[rc.phi(i) for i in I] for rc in RC[:4]])
            [1 1]
            [0 3]
            [0 2]
            [1 1]
        """
        a = self.parent()._rc_index_inverse[a]
        p_inf = self.parent()._calc_vacancy_number(self, a, float("inf"))
        if not self[a]:
            return Integer(p_inf)
        return Integer(p_inf - min(0, min(self[a].rigging)))

    def vacancy_number(self, a, i):
        r"""
        Return the vacancy number `p_i^{(a)}`.

        INPUT:

        - ``a`` -- the index of the rigged partition

        - ``i`` -- the row of the rigged partition

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[1], [2,1], [1], []])
            sage: elt.vacancy_number(2, 3)
            -2
            sage: elt.vacancy_number(2, 2)
            -2
            sage: elt.vacancy_number(2, 1)
            -1

            sage: RC = RiggedConfigurations(['D',4,1], [[2,1], [2,1]])
            sage: x = RC(partition_list=[[3], [3,1,1], [2], [3,1]]); ascii_art(x)
            -1[ ][ ][ ]-1  1[ ][ ][ ]1  0[ ][ ]0  -3[ ][ ][ ]-3
                           0[ ]0                  -1[ ]-1
                           0[ ]0
            sage: x.vacancy_number(2,2)
            1
        """
        a = self.parent()._rc_index_inverse[a]
        return self.parent()._calc_vacancy_number(self, a, i)

    def partition_rigging_lists(self):
        """
        Return the list of partitions and the associated list of riggings
        of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',3,1], [[1,2],[2,2]])
            sage: rc = RC(partition_list=[[2],[1],[1]], rigging_list=[[-1],[0],[-1]]); rc
            <BLANKLINE>
            -1[ ][ ]-1
            <BLANKLINE>
            1[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: rc.partition_rigging_lists()
            [[[2], [1], [1]], [[-1], [0], [-1]]]
        """
        partitions = []
        riggings = []
        for p in self:
            partitions.append(list(p))
            riggings.append(list(p.rigging))
        return [partitions, riggings]


class RCNonSimplyLacedElement(RiggedConfigurationElement):
    """
    Rigged configuration elements for non-simply-laced types.

    TESTS::

        sage: vct = CartanType(['C',2,1]).as_folding()
        sage: RC = crystals.infinity.RiggedConfigurations(vct)
        sage: elt = RC.module_generators[0].f_string([1,0,2,2,0,1]); elt
        <BLANKLINE>
        -2[ ][ ]-1
        <BLANKLINE>
        -2[ ]-1
        -2[ ]-1
        <BLANKLINE>
        -2[ ][ ]-1
        <BLANKLINE>
        sage: TestSuite(elt).run()
    """

    def to_virtual_configuration(self):
        """
        Return the corresponding rigged configuration in the virtual crystal.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
            sage: elt = RC(partition_list=[[3],[2]]); elt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: elt.to_virtual_configuration()
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ]0
        """
        return self.parent().to_virtual(self)

    def e(self, a):
        r"""
        Return the action of `e_a` on ``self``.

        This works by lifting into the virtual configuration, then applying

        .. MATH::

            e^v_a = \prod_{j \in \iota(a)} \hat{e}_j^{\gamma_j}

        and pulling back.

        EXAMPLES::

            sage: vct = CartanType(['C',2,1]).as_folding()
            sage: RC = crystals.infinity.RiggedConfigurations(vct)
            sage: elt = RC(partition_list=[[2],[1,1],[2]], rigging_list=[[-1],[-1,-1],[-1]])
            sage: ascii_art(elt.e(0))
            0[ ]0  -2[ ]-1  -2[ ][ ]-1
                   -2[ ]-1
            sage: ascii_art(elt.e(1))
            -3[ ][ ]-2  0[ ]1  -3[ ][ ]-2
            sage: ascii_art(elt.e(2))
            -2[ ][ ]-1  -2[ ]-1  0[ ]0
                        -2[ ]-1
        """
        vct = self.parent()._folded_ct
        L = []
        gamma = vct.scaling_factors()
        for i in vct.folding_orbit()[a]:
            L.extend([i]*gamma[a])
        virtual_rc = self.parent().to_virtual(self).e_string(L)
        if virtual_rc is None:
            return None
        return self.parent().from_virtual(virtual_rc)

    def f(self, a):
        r"""
        Return the action of `f_a` on ``self``.

        This works by lifting into the virtual configuration, then applying

        .. MATH::

            f^v_a = \prod_{j \in \iota(a)} \hat{f}_j^{\gamma_j}

        and pulling back.

        EXAMPLES::

            sage: vct = CartanType(['C',2,1]).as_folding()
            sage: RC = crystals.infinity.RiggedConfigurations(vct)
            sage: elt = RC(partition_list=[[2],[1,1],[2]], rigging_list=[[-1],[-1,-1],[-1]])
            sage: ascii_art(elt.f(0))
            -4[ ][ ][ ]-2  -2[ ]-1  -2[ ][ ]-1
                           -2[ ]-1
            sage: ascii_art(elt.f(1))
            -1[ ][ ]0  -2[ ][ ]-2  -1[ ][ ]0
                       -2[ ]-1
            sage: ascii_art(elt.f(2))
            -2[ ][ ]-1  -2[ ]-1  -4[ ][ ][ ]-2
                        -2[ ]-1
        """
        vct = self.parent()._folded_ct
        L = []
        gamma = vct.scaling_factors()
        for i in vct.folding_orbit()[a]:
            L.extend([i]*gamma[a])
        virtual_rc = self.parent().to_virtual(self).f_string(L)
        if virtual_rc is None:
            return None
        return self.parent().from_virtual(virtual_rc)

##########################################################
#  Highest weight crystal rigged configuration elements  #
##########################################################


class RCHighestWeightElement(RiggedConfigurationElement):
    """
    Rigged configurations in highest weight crystals.

    TESTS::

        sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
        sage: RC = crystals.RiggedConfigurations(['A',2,1], La[0])
        sage: elt = RC(partition_list=[[1,1],[1],[2]]); elt
        <BLANKLINE>
        -1[ ]-1
        -1[ ]-1
        <BLANKLINE>
        1[ ]1
        <BLANKLINE>
        -1[ ][ ]-1
        <BLANKLINE>
        sage: TestSuite(elt).run()
    """

    def check(self):
        """
        Make sure all of the riggings are less than or equal to the
        vacancy number.

        TESTS::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(['A',2,1], La[0])
            sage: elt = RC(partition_list=[[1,1],[1],[2]])
            sage: elt.check()
        """
        for a, partition in enumerate(self):
            for i, vac_num in enumerate(partition.vacancy_numbers):
                if vac_num < partition.rigging[i]:
                    raise ValueError("rigging can be at most the vacancy number")

    def f(self, a):
        r"""
        Return the action of the crystal operator `f_a` on ``self``.

        This implements the method defined in [CrysStructSchilling06]_ which
        finds the value `k` which is  the length of the string with the
        smallest nonpositive rigging of largest length. Then it adds a box
        from a string of length `k` in the `a`-th rigged partition, keeping
        all colabels fixed and decreasing the new label by one. If no such
        string exists, then it adds a new string of length 1 with label `-1`.
        If any of the resulting vacancy numbers are larger than the labels
        (i.e. it is an invalid rigged configuration), then `f_a` is
        undefined.

        INPUT:

        - ``a`` -- the index of the partition to add a box

        OUTPUT: the resulting rigged configuration element

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(['A',2,1], La[0])
            sage: elt = RC(partition_list=[[1,1],[1],[2]])
            sage: elt.f(0)
            <BLANKLINE>
            -2[ ][ ]-2
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            sage: elt.f(1)
            <BLANKLINE>
            0[ ]0
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            -1[ ]-1
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            sage: elt.f(2)
        """
        if not self.phi(a):
            return None
        return RiggedConfigurationElement.f(self, a)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: B = crystals.RiggedConfigurations(['A',2,1], La[0])
            sage: mg = B.module_generators[0]
            sage: mg.f_string([0,1,2,0]).weight()
            -Lambda[0] + Lambda[1] + Lambda[2] - 2*delta
        """
        P = self.parent().weight_lattice_realization()
        alpha = list(P.simple_roots())
        return self.parent()._wt - sum(sum(x) * alpha[i] for i,x in enumerate(self))


class RCHWNonSimplyLacedElement(RCNonSimplyLacedElement):
    """
    Rigged configurations in highest weight crystals.

    TESTS::

        sage: La = RootSystem(['C',2,1]).weight_lattice(extended=True).fundamental_weights()
        sage: vct = CartanType(['C',2,1]).as_folding()
        sage: RC = crystals.RiggedConfigurations(vct, La[0])
        sage: elt = RC(partition_list=[[1,1],[2],[2]]); ascii_art(elt)
        -1[ ]-1  2[ ][ ]2  -2[ ][ ]-2
        -1[ ]-1
        sage: TestSuite(elt).run()
    """

    def check(self):
        """
        Make sure all of the riggings are less than or equal to the
        vacancy number.

        TESTS::

            sage: La = RootSystem(['C',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: vct = CartanType(['C',2,1]).as_folding()
            sage: RC = crystals.RiggedConfigurations(vct, La[0])
            sage: elt = RC(partition_list=[[1,1],[2],[2]])
            sage: elt.check()
        """
        for partition in self:
            for i, vac_num in enumerate(partition.vacancy_numbers):
                if vac_num < partition.rigging[i]:
                    raise ValueError("rigging can be at most the vacancy number")

    def f(self, a):
        r"""
        Return the action of `f_a` on ``self``.

        This works by lifting into the virtual configuration, then applying

        .. MATH::

            f^v_a = \prod_{j \in \iota(a)} \hat{f}_j^{\gamma_j}

        and pulling back.

        EXAMPLES::

            sage: La = RootSystem(['C',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: vct = CartanType(['C',2,1]).as_folding()
            sage: RC = crystals.RiggedConfigurations(vct, La[0])
            sage: elt = RC(partition_list=[[1,1],[2],[2]])
            sage: elt.f(0)
            sage: ascii_art(elt.f(1))
            0[ ]0   0[ ][ ]0  -1[ ][ ]-1
            0[ ]0  -1[ ]-1
            sage: elt.f(2)
        """
        if not self.phi(a):
            return None
        return RCNonSimplyLacedElement.f(self, a)

    # FIXME: Do not duplicate with the simply-laced HW RC element class
    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['C',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: vct = CartanType(['C',2,1]).as_folding()
            sage: B = crystals.RiggedConfigurations(vct, La[0])
            sage: mg = B.module_generators[0]
            sage: mg.f_string([0,1,2]).weight()
            2*Lambda[1] - Lambda[2] - delta
        """
        P = self.parent().weight_lattice_realization()
        alpha = list(P.simple_roots())
        return self.parent()._wt - sum(sum(x) * alpha[i] for i,x in enumerate(self))

##############################################
#  KR crystal rigged configuration elements  #
##############################################


class KRRiggedConfigurationElement(RiggedConfigurationElement):
    r"""
    `U_q^{\prime}(\mathfrak{g})` rigged configurations.

    EXAMPLES:

    We can go between :class:`rigged configurations <RiggedConfigurations>`
    and tensor products of :class:`tensor products of KR tableaux
    <sage.combinat.rigged_configurations.tensor_product_kr_tableaux.TensorProductOfKirillovReshetikhinTableaux>`::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[1,1], [2,1]])
        sage: rc_elt = RC(partition_list=[[1], [1,1], [1], [1]])
        sage: tp_krtab = rc_elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); tp_krtab
        [[-2]] (X) [[1], [2]]
        sage: tp_krcrys = rc_elt.to_tensor_product_of_kirillov_reshetikhin_crystals(); tp_krcrys
        [[[-2]], [[1], [2]]]
        sage: tp_krcrys == tp_krtab.to_tensor_product_of_kirillov_reshetikhin_crystals()
        True
        sage: RC(tp_krcrys) == rc_elt
        True
        sage: RC(tp_krtab) == rc_elt
        True
        sage: tp_krtab.to_rigged_configuration() == rc_elt
        True
    """

    def __init__(self, parent, rigged_partitions=[], **options):
        r"""
        Construct a rigged configuration element.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: RC(partition_list=[[], [], [], []])
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: RC(partition_list=[[1], [1], [], []])
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: elt = RC(partition_list=[[1], [1], [], []], rigging_list=[[-1], [0], [], []]); elt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: TestSuite(elt).run()
        """
        n = len(parent._rc_index)
        if "KT_constructor" in options:
            # Used only by the Kleber tree
            # Not recommended to be called by the user since it avoids safety
            #   checks for speed
            data = options["KT_constructor"]
            shape_data = data[0]
            rigging_data = data[1]
            vac_data = data[2]
            nu = [RiggedPartition(a, b, c)
                  for a, b, c in zip(shape_data, rigging_data, vac_data)]
            # Special display case
            if parent.cartan_type().type() == 'B':
                nu[-1] = RiggedPartitionTypeB(nu[-1])
            ClonableArray.__init__(self, parent, nu)
            return
        RiggedConfigurationElement.__init__(self, parent, rigged_partitions, n=n, **options)
        # Special display case
        if parent.cartan_type().type() == 'B':
            self._set_mutable()
            self[-1] = RiggedPartitionTypeB(self[-1])
            self.set_immutable()

    def check(self):
        """
        Make sure all of the riggings are less than or equal to the
        vacancy number.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: elt.check()
        """
        for partition in self:
            for i, vac_num in enumerate(partition.vacancy_numbers):
                if vac_num < partition.rigging[i]:
                    raise ValueError("rigging can be at most the vacancy number")

    def e(self, a):
        r"""
        Return the action of the crystal operator `e_a` on ``self``.

        For the classical operators, this implements the method defined
        in [CrysStructSchilling06]_. For `e_0`, this converts the class to
        a tensor product of KR tableaux and does the corresponding `e_0`
        and pulls back.

        .. TODO::

            Implement `e_0` without appealing to tensor product of
            KR tableaux.

        INPUT:

        - ``a`` -- the index of the partition to remove a box

        OUTPUT: the resulting rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: elt = RC(partition_list=[[1], [1], [1], [1]])
            sage: elt.e(3)
            sage: elt.e(1)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
        """
        if a not in self.parent()._cartan_type.index_set():
            raise ValueError("{} is not in the index set".format(a))
        if a == self.parent()._cartan_type.special_node():
            try:
                ret = self.to_tensor_product_of_kirillov_reshetikhin_tableaux().e(a)
                if ret is None:
                    return None
                return ret.to_rigged_configuration()
            except NotImplementedError:
                # We haven't implemented the bijection yet, so return None
                # This is to make sure we can at least view it as a classical
                #   crystal if there is no bijection.
                return None

        return RiggedConfigurationElement.e(self, a)

    def f(self, a):
        r"""
        Return the action of the crystal operator `f_a` on ``self``.

        For the classical operators, this implements the method defined
        in [CrysStructSchilling06]_. For `f_0`, this converts the class to
        a tensor product of KR tableaux and does the corresponding `f_0`
        and pulls back.

        .. TODO::

            Implement `f_0` without appealing to tensor product of
            KR tableaux.

        INPUT:

        - ``a`` -- the index of the partition to add a box

        OUTPUT: the resulting rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: elt = RC(partition_list=[[1], [1], [1], [1]])
            sage: elt.f(1)
            sage: elt.f(2)
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
        """
        ct = self.parent()._cartan_type
        if a not in ct.index_set():
            raise ValueError("{} is not in the index set".format(a))
        if a == ct.special_node():
            try:
                ret = self.to_tensor_product_of_kirillov_reshetikhin_tableaux().f(a)
                if ret is None:
                    return None
                return ret.to_rigged_configuration()
            except NotImplementedError:
                # We haven't implemented the bijection yet, so return None
                # This is to make sure we can at least view it as a classical
                #   crystal if there is no bijection.
                return None

        if not self.phi(a):
            return None

        return RiggedConfigurationElement.f(self, a)

    def epsilon(self, a):
        r"""
        Return `\varepsilon_a` of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: I = RC.index_set()
            sage: matrix([[mg.epsilon(i) for i in I] for mg in RC.module_generators])
            [4 0 0 0 0]
            [3 0 0 0 0]
            [2 0 0 0 0]
        """
        if a == self.parent()._cartan_type.special_node():
            return self.to_tensor_product_of_kirillov_reshetikhin_tableaux().epsilon(a)
        return RiggedConfigurationElement.epsilon(self, a)

    def phi(self, a):
        r"""
        Return `\varphi_a` of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: I = RC.index_set()
            sage: matrix([[mg.phi(i) for i in I] for mg in RC.module_generators])
            [0 0 2 0 0]
            [1 0 1 0 0]
            [2 0 0 0 0]
        """
        if a == self.parent()._cartan_type.special_node():
            return self.to_tensor_product_of_kirillov_reshetikhin_tableaux().phi(a)
        return RiggedConfigurationElement.phi(self, a)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['E', 6, 1], [[2,2]])
            sage: [x.weight() for x in RC.module_generators]
            [-4*Lambda[0] + 2*Lambda[2], -2*Lambda[0] + Lambda[2], 0]
            sage: KR = crystals.KirillovReshetikhin(['E',6,1], 2,2)
            sage: [x.weight() for x in KR.module_generators]  # long time
            [0, -2*Lambda[0] + Lambda[2], -4*Lambda[0] + 2*Lambda[2]]

            sage: RC = RiggedConfigurations(['D', 6, 1], [[4,2]])
            sage: [x.weight() for x in RC.module_generators]
            [-4*Lambda[0] + 2*Lambda[4], -4*Lambda[0] + Lambda[2] + Lambda[4],
             -2*Lambda[0] + Lambda[4], -4*Lambda[0] + 2*Lambda[2],
             -2*Lambda[0] + Lambda[2], 0]
        """
        WLR = self.parent().weight_lattice_realization()
        La = WLR.fundamental_weights()
        cl_index = self.parent()._rc_index
        wt = WLR.sum((self.phi(i) - self.epsilon(i)) * La[i] for i in cl_index)
        return -wt.level() * La[0] + wt

    @cached_method
    def classical_weight(self):
        r"""
        Return the classical weight of ``self``.

        The classical weight `\Lambda` of a rigged configuration is

        .. MATH::

            \Lambda = \sum_{a \in \overline{I}} \sum_{i > 0}
            i L_i^{(a)} \Lambda_a - \sum_{a \in \overline{I}} \sum_{i > 0}
            i m_i^{(a)} \alpha_a.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D',4,1], [[2,2]])
            sage: elt = RC(partition_list=[[2],[2,1],[1],[1]])
            sage: elt.classical_weight()
            (0, 1, 1, 0)

        This agrees with the corresponding classical weight as KR tableaux::

            sage: krt = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); krt
            [[2, 1], [3, -1]]
            sage: krt.classical_weight() == elt.classical_weight()
            True

        TESTS:

        We check the classical weights agree in an entire crystal::

            sage: RC = RiggedConfigurations(['A',2,1], [[2,1], [1,1]])
            sage: for x in RC:
            ....:     y = x.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            ....:     assert x.classical_weight() == y.classical_weight()
        """
        F = self.cartan_type().classical().root_system()
        if F.ambient_space() is None:
            WLR = F.weight_lattice()
        else:
            WLR = F.ambient_space()
        La = WLR.fundamental_weights()
        wt = WLR.sum(La[r] * s for r,s in self.parent().dims)

        alpha = WLR.simple_roots()
        rc_index = self.parent()._rc_index
        for a, nu in enumerate(self):
            wt -= sum(nu) * alpha[rc_index[a]]
        return wt

    def to_tensor_product_of_kirillov_reshetikhin_tableaux(self, display_steps=False, build_graph=False):
        r"""
        Perform the bijection from this rigged configuration to a tensor
        product of Kirillov-Reshetikhin tableaux given in [RigConBijection]_
        for single boxes and with [BijectionLRT]_ and [BijectionDn]_ for
        multiple columns and rows.

        .. NOTE::

            This is only proven to be a bijection in types `A_n^{(1)}`
            and `D_n^{(1)}`, as well as `\bigotimes_i B^{r_i,1}` and
            `\bigotimes_i B^{1,s_i}` for general affine types.

        INPUT:

        - ``display_steps`` -- boolean (default: ``False``); indicates whether
          to print each step in the algorithm
        - ``build_graph`` -- boolean (default: ``False``); indicates whether
          to construct and return a graph of the bijection whose
          vertices are rigged configurations obtained at each step and edges
          are labeled by either the return value of `\delta` or the
          doubling/halving map

        OUTPUT:

        - The tensor product of KR tableaux element corresponding to this
          rigged configuration.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [2,2], [2], [2]]).to_tensor_product_of_kirillov_reshetikhin_tableaux()
            [[3, 3], [5, 5]]
            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[2], [2,2], [1], [1]])
            sage: tp_krt = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); tp_krt
            [[2, 3], [3, -2]]

        This is invertible by calling
        :meth:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element.TensorProductOfKirillovReshetikhinTableauxElement.to_rigged_configuration()`::

            sage: ret = tp_krt.to_rigged_configuration(); ret
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            -2[ ][ ]-2
            -2[ ][ ]-2
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            sage: elt == ret
            True

        To view the steps of the bijection in the output, run with
        the ``display_steps=True`` option::

            sage: elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(True)
            ====================
            ...
            ====================
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -2[ ][ ]-2
             0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            --------------------
            [[3, 2]]
            --------------------
            ...
            [[2, 3], [3, -2]]

        We can also construct and display a graph of the bijection
        as follows::

            sage: ret, G = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(build_graph=True)
            sage: view(G) # not tested
        """
        from sage.combinat.rigged_configurations.bijection import RCToKRTBijection
        bij = RCToKRTBijection(self)
        ret = bij.run(display_steps, build_graph)
        if build_graph:
            return (ret, bij._graph)
        return ret

    def to_tensor_product_of_kirillov_reshetikhin_crystals(self, display_steps=False, build_graph=False):
        r"""
        Return the corresponding tensor product of Kirillov-Reshetikhin
        crystals.

        This is a composition of the map to a tensor product of KR tableaux,
        and then to a tensor product of KR crystals.

        INPUT:

        - ``display_steps`` -- boolean (default: ``False``); indicates whether
          to print each step in the algorithm
        - ``build_graph`` -- boolean (default: ``False``); indicates whether
          to construct and return a graph of the bijection whose
          vertices are rigged configurations obtained at each step and edges
          are labeled by either the return value of `\delta` or the
          doubling/halving map

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[2], [2,2], [1], [1]])
            sage: krc = elt.to_tensor_product_of_kirillov_reshetikhin_crystals(); krc
            [[[2, 3], [3, -2]]]

        We can recover the rigged configuration::

            sage: ret = RC(krc); ret
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            -2[ ][ ]-2
            -2[ ][ ]-2
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            sage: elt == ret
            True

        We can also construct and display a graph of the bijection
        as follows::

            sage: ret, G = elt.to_tensor_product_of_kirillov_reshetikhin_crystals(build_graph=True)
            sage: view(G) # not tested
        """
        if build_graph:
            kr_tab, G = self.to_tensor_product_of_kirillov_reshetikhin_tableaux(display_steps, build_graph)
            return (kr_tab.to_tensor_product_of_kirillov_reshetikhin_crystals(), G)
        kr_tab = self.to_tensor_product_of_kirillov_reshetikhin_tableaux(display_steps)
        return kr_tab.to_tensor_product_of_kirillov_reshetikhin_crystals()

    # TODO: Move the morphisms to a lazy attribute of RiggedConfigurations
    #   once #15463 is done
    def left_split(self):
        r"""
        Return the image of ``self`` under the left column splitting
        map `\beta`.

        Consider the map `\beta : RC(B^{r,s} \otimes B) \to RC(B^{r,1}
        \otimes B^{r,s-1} \otimes B)` for `s > 1` which is a natural classical
        crystal injection. On rigged configurations, the map `\beta` does
        nothing (except possibly changing the vacancy numbers).

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',4,1], [[3,3]])
            sage: mg = RC.module_generators[-1]
            sage: ascii_art(mg)
            0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ]0
                      0[ ][ ]0  0[ ][ ]0  0[ ]0
                                0[ ][ ]0  0[ ]0
            sage: ascii_art(mg.left_split())
            0[ ][ ]0  0[ ][ ]0  1[ ][ ]0  0[ ]0
                      0[ ][ ]0  1[ ][ ]0  0[ ]0
                                1[ ][ ]0  0[ ]0
        """
        P = self.parent()
        if P.dims[0][1] == 1:
            raise ValueError("cannot split a single column")
        r,s = P.dims[0]
        B = [[r,1], [r,s-1]]
        B.extend(P.dims[1:])
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        RC = RiggedConfigurations(P._cartan_type, B)
        return RC(*[x._clone() for x in self]) # Make a deep copy

    def right_split(self):
        r"""
        Return the image of ``self`` under the right column splitting
        map `\beta^*`.

        Let `\theta` denote the
        :meth:`complement rigging map<complement_rigging>` which reverses
        the tensor factors and `\beta` denote the
        :meth:`left splitting map<left_split>`, we define the right
        splitting map by `\beta^* := \theta \circ \beta \circ \theta`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',4,1], [[3,3]])
            sage: mg = RC.module_generators[-1]
            sage: ascii_art(mg)
            0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ]0
                      0[ ][ ]0  0[ ][ ]0  0[ ]0
                                0[ ][ ]0  0[ ]0
            sage: ascii_art(mg.right_split())
            0[ ][ ]0  0[ ][ ]0  1[ ][ ]1  0[ ]0
                      0[ ][ ]0  1[ ][ ]1  0[ ]0
                                1[ ][ ]1  0[ ]0

            sage: RC = RiggedConfigurations(['D',4,1], [[2,2],[1,2]])
            sage: elt = RC(partition_list=[[3,1], [2,2,1], [2,1], [2]])
            sage: ascii_art(elt)
            -1[ ][ ][ ]-1  0[ ][ ]0  -1[ ][ ]-1  1[ ][ ]1
             0[ ]0         0[ ][ ]0  -1[ ]-1
                           0[ ]0
            sage: ascii_art(elt.right_split())
            -1[ ][ ][ ]-1  0[ ][ ]0  -1[ ][ ]-1  1[ ][ ]1
             1[ ]0         0[ ][ ]0  -1[ ]-1
                           0[ ]0

        We check that the bijection commutes with the right splitting map::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[1,1], [2,2]])
            sage: all(rc.right_split().to_tensor_product_of_kirillov_reshetikhin_tableaux()
            ....:     == rc.to_tensor_product_of_kirillov_reshetikhin_tableaux().right_split() for rc in RC)
            True
        """
        return self.complement_rigging(True).left_split().complement_rigging(True)

    def left_box(self, return_b=False):
        r"""
        Return the image of ``self`` under the left box removal map `\delta`.

        The map `\delta : RC(B^{r,1} \otimes B) \to RC(B^{r-1,1}
        \otimes B)` (if `r = 1`, then we remove the left-most factor) is the
        basic map in the bijection `\Phi` between rigged configurations and
        tensor products of Kirillov-Reshetikhin tableaux. For more
        information, see
        :meth:`to_tensor_product_of_kirillov_reshetikhin_tableaux()`.
        We can extend `\delta` when the left-most factor is not a single
        column by precomposing with a :meth:`left_split()`.

        .. NOTE::

            Due to the special nature of the bijection for the spinor cases in
            types `D_n^{(1)}`, `B_n^{(1)}`, and `A_{2n-1}^{(2)}`, this map is
            not defined in these cases.

        INPUT:

        - ``return_b`` -- boolean (default: ``False``); whether to return the
          resulting letter from `\delta`

        OUTPUT:

        The resulting rigged configuration or if ``return_b`` is ``True``,
        then a tuple of the resulting rigged configuration and the letter.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',4,1], [[3,2]])
            sage: mg = RC.module_generators[-1]
            sage: ascii_art(mg)
            0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ]0
                      0[ ][ ]0  0[ ][ ]0  0[ ]0
                                0[ ][ ]0  0[ ]0
            sage: ascii_art(mg.left_box())
            0[ ]0  0[ ][ ]0  0[ ][ ]0  0[ ]0
                   0[ ]0     0[ ][ ]0  0[ ]0
            sage: x,b = mg.left_box(True)
            sage: b
            -1
            sage: b.parent()
            The crystal of letters for type ['C', 4]
        """
        # Don't do spinor cases
        P = self.parent()
        ct = P.cartan_type()
        if ct.type() == 'D':
            if P.dims[0][0] >= ct.rank() - 2:
                raise ValueError("only for non-spinor cases")
        elif ct.type() == 'B' or ct.dual().type() == 'B':
            if P.dims[0][0] == ct.rank() - 1:
                raise ValueError("only for non-spinor cases")

        from sage.combinat.rigged_configurations.bijection import RCToKRTBijection
        rc = self
        if P.dims[0][1] != 1:
            rc = self.left_split()
        bij = RCToKRTBijection(rc)
        ht = bij.cur_dims[0][0]
        bij.cur_dims[0][0] = bij._next_index(ht)
        b = bij.next_state(ht)
        if bij.cur_dims[0][0] == 0:
            bij.cur_dims.pop(0)
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        RC = RiggedConfigurations(ct, bij.cur_dims)
        rc = RC(*bij.cur_partitions)
        if return_b:
            from sage.combinat.crystals.letters import CrystalOfLetters
            L = CrystalOfLetters(self.parent()._cartan_type.classical())
            return (rc, L(b))
        return rc

    delta = left_box

    def left_column_box(self):
        r"""
        Return the image of ``self`` under the left column box splitting
        map `\gamma`.

        Consider the map `\gamma : RC(B^{r,1} \otimes B) \to RC(B^{1,1}
        \otimes B^{r-1,1} \otimes B)` for `r > 1`, which is a natural strict
        classical crystal injection. On rigged configurations, the map
        `\gamma` adds a singular string of length `1` to `\nu^{(a)}`.

        We can extend `\gamma` when the left-most factor is not a single
        column by precomposing with a :meth:`left_split()`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',3,1], [[3,1], [2,1]])
            sage: mg = RC.module_generators[-1]
            sage: ascii_art(mg)
            0[ ]0  0[ ][ ]0  0[ ]0
                   0[ ]0     0[ ]0
            sage: ascii_art(mg.left_column_box())
            0[ ]0  0[ ][ ]0  0[ ]0
            0[ ]0  0[ ]0     0[ ]0
                   0[ ]0

            sage: RC = RiggedConfigurations(['C',3,1], [[2,1], [1,1], [3,1]])
            sage: mg = RC.module_generators[7]
            sage: ascii_art(mg)
            1[ ]0  0[ ][ ]0  0[ ]0
                   0[ ]0     0[ ]0
            sage: ascii_art(mg.left_column_box())
            1[ ]1  0[ ][ ]0  0[ ]0
            1[ ]0  0[ ]0     0[ ]0
        """
        P = self.parent()
        r = P.dims[0][0]
        if r == 1:
            raise ValueError("cannot split a single box")
        ct = P.cartan_type()
        if ct.type() == 'D':
            if P.dims[0][0] >= ct.rank() - 2:
                raise ValueError("only for non-spinor cases")
        elif ct.type() == 'B' or ct.dual().type() == 'B':
            if P.dims[0][0] == ct.rank() - 1:
                raise ValueError("only for non-spinor cases")

        if P.dims[0][1] > 1:
            return self.left_split().left_column_box()

        B = [[1,1], [r-1,1]]
        B.extend(P.dims[1:])
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        RC = RiggedConfigurations(P._cartan_type, B)
        parts = [x._clone() for x in self] # Make a deep copy
        for nu in parts[:r-1]:
            nu._list.append(1)
        for a, nu in enumerate(parts[:r-1]):
            vac_num = RC._calc_vacancy_number(parts, a, 1)
            i = nu._list.index(1)
            nu.vacancy_numbers.insert(i, vac_num)
            nu.rigging.insert(i, vac_num)
        return RC(*parts)

    def right_column_box(self):
        r"""
        Return the image of ``self`` under the right column box splitting
        map `\gamma^*`.

        Consider the map `\gamma^* : RC(B \otimes B^{r,1}) \to RC(B \otimes
        B^{r-1,1} \otimes B^{1,1})` for `r > 1`, which is a natural strict
        classical crystal injection. On rigged configurations, the map
        `\gamma` adds a string of length `1` with rigging 0 to `\nu^{(a)}`
        for all `a < r` to a classically highest weight element and extended
        as a classical crystal morphism.

        We can extend `\gamma^*` when the right-most factor is not a single
        column by precomposing with a :meth:`right_split()`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',3,1], [[2,1], [1,1], [3,1]])
            sage: mg = RC.module_generators[7]
            sage: ascii_art(mg)
            1[ ]0  0[ ][ ]0  0[ ]0
                   0[ ]0     0[ ]0
            sage: ascii_art(mg.right_column_box())
            1[ ]0  0[ ][ ]0  0[ ]0
            1[ ]0  0[ ]0     0[ ]0
                   0[ ]0
        """
        P = self.parent()
        r = P.dims[-1][0]
        if r == 1:
            raise ValueError("cannot split a single box")
        ct = P.cartan_type()
        if ct.type() == 'D':
            if P.dims[-1][0] >= ct.rank() - 2:
                raise ValueError("only for non-spinor cases")
        elif ct.type() == 'B' or ct.dual().type() == 'B':
            if P.dims[-1][0] == ct.rank() - 1:
                raise ValueError("only for non-spinor cases")

        if P.dims[-1][1] > 1:
            return self.right_split().right_column_box()

        rc, e_string = self.to_highest_weight(P._rc_index)

        B = P.dims[:-1] + ([r-1,1], [1,1])
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        RC = RiggedConfigurations(P._cartan_type, B)
        parts = [x._clone() for x in rc] # Make a deep copy
        for nu in parts[:r-1]:
            nu._list.append(1)
        for a, nu in enumerate(parts[:r-1]):
            vac_num = RC._calc_vacancy_number(parts, a, -1)
            nu.vacancy_numbers.append(vac_num)
            nu.rigging.append(0)
        return RC(*parts).f_string(reversed(e_string))

    def complement_rigging(self, reverse_factors=False):
        r"""
        Apply the complement rigging morphism `\theta` to ``self``.

        Consider a highest weight rigged configuration `(\nu, J)`, the
        complement rigging morphism `\theta : RC(L) \to RC(L)` is given by
        sending `(\nu, J) \mapsto (\nu, J')`, where `J'` is obtained by
        taking the coriggings `x' = p_i^{(a)} - x`, and then extending as
        a crystal morphism. (The name comes from taking the complement
        partition for the riggings in a `m_i^{(a)} \times p_i^{(a)}` box.)

        INPUT:

        - ``reverse_factors`` -- boolean (default: ``False``); if ``True``, then this
          returns an element in `RC(B')` where `B'` is the tensor factors
          of ``self`` in reverse order

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D',4,1], [[1,1],[2,2]])
            sage: mg = RC.module_generators[-1]
            sage: ascii_art(mg)
            1[ ][ ]1  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0
                      0[ ][ ]0
            sage: ascii_art(mg.complement_rigging())
            1[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0
                      0[ ][ ]0

            sage: lw = mg.to_lowest_weight([1,2,3,4])[0]
            sage: ascii_art(lw)
            -1[ ][ ]-1  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0
            -1[ ]-1     0[ ][ ]0  0[ ]0     0[ ]0
            -1[ ]-1     0[ ]0
                        0[ ]0
            sage: ascii_art(lw.complement_rigging())
            -1[ ][ ][ ]-1  0[ ][ ][ ]0  0[ ][ ][ ]0  0[ ][ ][ ]0
            -1[ ]-1        0[ ][ ][ ]0
            sage: lw.complement_rigging() == mg.complement_rigging().to_lowest_weight([1,2,3,4])[0]
            True

            sage: mg.complement_rigging(True).parent()
            Rigged configurations of type ['D', 4, 1] and factor(s) ((2, 2), (1, 1))

        We check that the Lusztig involution (under the modification of also
        mapping to the highest weight element) intertwines with the
        complement map `\theta` (that reverses the tensor factors)
        under the bijection `\Phi`::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2], [2, 1], [1, 2]])
            sage: for mg in RC.module_generators: # long time
            ....:     y = mg.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            ....:     hw = y.lusztig_involution().to_highest_weight([1,2,3,4])[0]
            ....:     c = mg.complement_rigging(True)
            ....:     hwc = c.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            ....:     assert hw == hwc

        TESTS:

        We check that :issue:`18898` is fixed::

            sage: RC = RiggedConfigurations(['D',4,1], [[2,1], [2,1], [2,3]])
            sage: x = RC(partition_list=[[1], [1,1], [1], [1]], rigging_list=[[0], [2,1], [0], [0]])
            sage: ascii_art(x)
            0[ ]0  2[ ]2  0[ ]0  0[ ]0
                   2[ ]1
            sage: ascii_art(x.complement_rigging())
            0[ ]0  2[ ]1  0[ ]0  0[ ]0
                   2[ ]0
        """
        P = self.parent()
        if reverse_factors:
            from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
            P = RiggedConfigurations(P._cartan_type, reversed(P.dims))

        mg, e_str = self.to_highest_weight(P._rc_index)
        nu = []
        rig = []
        for a,p in enumerate(mg):
            nu.append(list(p))
            vac_nums = p.vacancy_numbers
            riggings = [vac - p.rigging[i] for i,vac in enumerate(vac_nums)]
            block = 0
            for j,i in enumerate(p):
                if p[block] != i:
                    riggings[block:j] = sorted(riggings[block:j], reverse=True)
                    block = j
            riggings[block:] = sorted(riggings[block:], reverse=True)
            rig.append(riggings)

        rc = P(partition_list=nu, rigging_list=rig)
        return rc.f_string(reversed(e_str))


class KRRCSimplyLacedElement(KRRiggedConfigurationElement):
    r"""
    `U_q^{\prime}(\mathfrak{g})` rigged configurations in simply-laced types.

    TESTS::

        sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [2,1], [1,1]])
        sage: elt = RC(partition_list=[[1], [1], []]); elt
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>
        (/)
        <BLANKLINE>
        sage: TestSuite(elt).run()
    """
    @cached_method
    def cocharge(self):
        r"""
        Compute the cocharge statistic of ``self``.

        Computes the cocharge statistic [CrysStructSchilling06]_ on this
        rigged configuration `(\nu, J)`. The cocharge statistic is defined as:

        .. MATH::

            cc(\nu, J) = \frac{1}{2} \sum_{a, b \in I_0}
            \sum_{j,k > 0} \left( \alpha_a \mid \alpha_b \right)
            \min(j, k) m_j^{(a)} m_k^{(b)}
            + \sum_{a \in I} \sum_{i > 0} \left\lvert J^{(a, i)} \right\rvert.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [2,1], [1,1]])
            sage: RC(partition_list=[[1], [1], []]).cocharge()
            1
        """
        cc = 0
        rigging_sum = 0
        for a, p in enumerate(self):
            for pos, i in enumerate(p._list):
                # Add the rigging
                rigging_sum += p.rigging[pos]
                # Add the L matrix contribution
                for dim in self.parent().dims:
                    if dim[0] == a + 1:
                        cc += min(dim[1], i)
                # Subtract the vacancy number
                cc -= p.vacancy_numbers[pos]
        return cc // 2 + rigging_sum

    cc = cocharge

    @cached_method
    def charge(self):
        r"""
        Compute the charge statistic of ``self``.

        Let `B` denote a set of rigged configurations. The *charge* `c` of
        a rigged configuration `b` is computed as

        .. MATH::

            c(b) = \max(cc(b) \mid b \in B) - cc(b).

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [2,1], [1,1]])
            sage: RC(partition_list=[[],[],[]]).charge()
            2
            sage: RC(partition_list=[[1], [1], []]).charge()
            1
        """
        B = self.parent()
        if not hasattr(B, "_max_charge"):
            B._max_charge = max(b.cocharge() for b in B.module_generators)
        return B._max_charge - self.cocharge()


class KRRCNonSimplyLacedElement(KRRiggedConfigurationElement, RCNonSimplyLacedElement):
    r"""
    `U_q^{\prime}(\mathfrak{g})` rigged configurations in non-simply-laced
    types.

    TESTS::

        sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
        sage: elt = RC(partition_list=[[3],[2]]); elt
        <BLANKLINE>
        0[ ][ ][ ]0
        <BLANKLINE>
        0[ ][ ]0
        sage: TestSuite(elt).run()
    """

    def e(self, a):
        r"""
        Return the action of `e_a` on ``self``.

        This works by lifting into the virtual configuration, then applying

        .. MATH::

            e^v_a = \prod_{j \in \iota(a)} \hat{e}_j^{\gamma_j}

        and pulling back.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',6,2], [[1,1]]*7)
            sage: elt = RC(partition_list=[[1]*5,[2,1,1],[3,2]])
            sage: elt.e(3)
            <BLANKLINE>
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            <BLANKLINE>
            0[ ][ ]0
            1[ ]1
            1[ ]1
            <BLANKLINE>
            1[ ][ ]1
            1[ ]0
            <BLANKLINE>
        """
        if a == self.parent()._cartan_type.special_node():
            try:
                ret = self.to_tensor_product_of_kirillov_reshetikhin_tableaux().e(a)
                if ret is None:
                    return None
                return ret.to_rigged_configuration()
            except (NotImplementedError, TypeError):
                # We haven't implemented the bijection yet, so try by lifting
                #   to the simply-laced case
                return RCNonSimplyLacedElement.e(self, a)

        if not self.epsilon(a):
            return None
        return RCNonSimplyLacedElement.e(self, a)

    def f(self, a):
        r"""
        Return the action of `f_a` on ``self``.

        This works by lifting into the virtual configuration, then applying

        .. MATH::

            f^v_a = \prod_{j \in \iota(a)} \hat{f}_j^{\gamma_j}

        and pulling back.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',6,2], [[1,1]]*7)
            sage: elt = RC(partition_list=[[1]*5,[2,1,1],[2,1]], rigging_list=[[0]*5,[0,1,1],[1,0]])
            sage: elt.f(3)
            <BLANKLINE>
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            <BLANKLINE>
            1[ ][ ]1
            1[ ]1
            1[ ]1
            <BLANKLINE>
            -1[ ][ ][ ]-1
             0[ ][ ]0
            <BLANKLINE>
        """
        if a == self.parent()._cartan_type.special_node():
            try:
                ret = self.to_tensor_product_of_kirillov_reshetikhin_tableaux().f(a)
                if ret is None:
                    return None
                return ret.to_rigged_configuration()
            except (NotImplementedError, TypeError):
                # We haven't implemented the bijection yet, so try by lifting
                #   to the simply-laced case
                return RCNonSimplyLacedElement.f(self, a)

        if not self.phi(a):
            return None
        return RCNonSimplyLacedElement.f(self, a)

    @cached_method
    def cocharge(self):
        r"""
        Compute the cocharge statistic.

        Computes the cocharge statistic [OSS03]_ on this
        rigged configuration `(\nu, J)` by computing the cocharge as a virtual
        rigged configuration `(\hat{\nu}, \hat{J})` and then using the
        identity `cc(\hat{\nu}, \hat{J}) = \gamma_0 cc(\nu, J)`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C', 3, 1], [[2,1], [1,1]])
            sage: RC(partition_list=[[1,1],[2,1],[1,1]]).cocharge()
            1
        """
        # return self.to_virtual_configuration().cocharge() / self.parent()._folded_ct.gamma[0]
        vct = self.parent()._folded_ct
        cc = ZZ.zero()
        rigging_sum = ZZ.zero()
        sigma = vct.folding_orbit()
        gamma = vct.scaling_factors()
        for a, p in enumerate(self):
            t_check = len(sigma[a + 1]) * gamma[a+1] // gamma[0]
            for pos, i in enumerate(p._list):
                # Add the rigging
                rigging_sum += t_check * p.rigging[pos]
                # Add the L matrix contribution
                for dim in self.parent().dims:
                    if dim[0] == a + 1:
                        cc += t_check * min(dim[1], i)
                # Subtract the vacancy number
                cc -= t_check * p.vacancy_numbers[pos]
        return cc // 2 + rigging_sum

    cc = cocharge


class KRRCTypeA2DualElement(KRRCNonSimplyLacedElement):
    r"""
    `U_q^{\prime}(\mathfrak{g})` rigged configurations in type
    `A_{2n}^{(2)\dagger}`.
    """

    def epsilon(self, a):
        r"""
        Return the value of `\varepsilon_a` of ``self``.

        Here we need to modify the usual definition by
        `\varepsilon_n^{\prime} := 2 \varepsilon_n`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[1,1], [2,2]])
            sage: def epsilon(x, i):
            ....:     x = x.e(i)
            ....:     eps = 0
            ....:     while x is not None:
            ....:         x = x.e(i)
            ....:         eps = eps + 1
            ....:     return eps
            sage: all(epsilon(rc, 2) == rc.epsilon(2) for rc in RC)
            True
        """
        if a == self.parent()._cartan_type.special_node():
            return self.to_tensor_product_of_kirillov_reshetikhin_tableaux().epsilon(a)

        a = self.parent()._rc_index_inverse[a]
        if not self[a]:
            epsilon = 0
        else:
            epsilon = -min(0, min(self[a].rigging))
        n = len(self.parent()._rc_index)
        if a == n-1: # -1 for indexing
            epsilon *= 2
        return Integer(epsilon)

    def phi(self, a):
        r"""
        Return the value of `\varphi_a` of ``self``.

        Here we need to modify the usual definition by
        `\varphi_n^{\prime} := 2 \varphi_n`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[1,1], [2,2]])
            sage: def phi(x, i):
            ....:     x = x.f(i)
            ....:     ph = 0
            ....:     while x is not None:
            ....:         x = x.f(i)
            ....:         ph = ph + 1
            ....:     return ph
            sage: all(phi(rc, 2) == rc.phi(2) for rc in RC)
            True
        """
        if a == self.parent()._cartan_type.special_node():
            return self.to_tensor_product_of_kirillov_reshetikhin_tableaux().phi(a)

        a = self.parent()._rc_index_inverse[a]
        p_inf = self.parent()._calc_vacancy_number(self, a, float("inf"))
        if not self[a]:
            phi = p_inf
        else:
            phi = p_inf - min(0, min(self[a].rigging))
        n = len(self.parent()._rc_index)
        if a == n-1: # -1 for indexing
            phi *= 2
        return Integer(phi)

    @cached_method
    def cocharge(self):
        r"""
        Compute the cocharge statistic.

        Computes the cocharge statistic [RigConBijection]_ on this
        rigged configuration `(\nu, J)`. The cocharge statistic is
        computed as:

        .. MATH::

            cc(\nu, J) = \frac{1}{2} \sum_{a \in I_0} \sum_{i > 0}
            t_a^{\vee} m_i^{(a)} \left( \sum_{j > 0} \min(i, j) L_j^{(a)}
            - p_i^{(a)} \right) + \sum_{a \in I} t_a^{\vee} \sum_{i > 0}
            \left\lvert J^{(a, i)} \right\rvert.

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[1,1],[2,2]])
            sage: sc = RC.cartan_type().as_folding().scaling_factors()
            sage: all(mg.cocharge() * sc[0] == mg.to_virtual_configuration().cocharge()
            ....:     for mg in RC.module_generators)
            True
        """
        # return self.to_virtual_configuration().cocharge() / self.parent()._folded_ct.gamma[0]
        cc = ZZ.zero()
        rigging_sum = ZZ.zero()
        # vct = self.parent()._folded_ct
        # sigma = vct.folding_orbit()
        # gammatilde = list(vct.scaling_factors())
        # gammatilde[-1] = 2
        for a, p in enumerate(self):
            t_check = 1  # == len(sigma[a+1]) * gammatilde[a+1] / gammatilde[0]
            for pos, i in enumerate(p._list):
                # Add the rigging
                rigging_sum += t_check * p.rigging[pos]
                # Add the L matrix contribution
                for dim in self.parent().dims:
                    if dim[0] == a + 1:
                        cc += t_check * min(dim[1], i)
                # Subtract the vacancy number
                cc -= t_check * p.vacancy_numbers[pos]
        return cc / ZZ(2) + rigging_sum

    cc = cocharge
