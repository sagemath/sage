# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection classes for type `C_n^{(1)}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `C_n^{(1)}`.

AUTHORS:

- Travis Scrimshaw (2012-12-21): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['C', 3, 1], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_C import KRTToRCBijectionTypeC
    sage: bijection = KRTToRCBijectionTypeC(KRT(pathlist=[[-1,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['C', 3, 1], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_C import RCToKRTBijectionTypeC
    sage: bijection = RCToKRTBijectionTypeC(RC(partition_list=[[],[],[]]))
    sage: TestSuite(bijection).run()
"""

# ****************************************************************************
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA


class KRTToRCBijectionTypeC(KRTToRCBijectionTypeA):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `C_n^{(1)}`.

    This inherits from type `A_n^{(1)}` because we use the same methods in
    some places.
    """

    def next_state(self, val):
        r"""
        Build the next state for type `C_n^{(1)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['C', 3, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_C import KRTToRCBijectionTypeC
            sage: bijection = KRTToRCBijectionTypeC(KRT(pathlist=[[-1,2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [2])
            sage: bijection.next_state(2)
        """
        # Note that we must subtract 1 from n to match the indices.
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        # If it is a regular value, we follow the A_n rules
        if val > 0:
            KRTToRCBijectionTypeA.next_state(self, val)
            return

        pos_val = -val
        case_S = [None] * n

        # Always add a cell to the first singular value in the first
        #   tableau we are updating.
        if len(self.ret_rig_con[pos_val - 1]) > 0:
            max_width = self.ret_rig_con[pos_val - 1][0]
        else:
            max_width = 1

        # Special case for inserting -n
        if pos_val == n:
            max_width *= 2

        # Add cells similar to type A_n but we move to the right until we
        #   reach the value of n
        for a in range(pos_val - 1, n - 1):
            max_width = self.ret_rig_con[a].insert_cell(max_width)
            case_S[a] = max_width

        # Special case for n
        max_width = self.ret_rig_con[n-1].insert_cell(max_width // 2) * 2

        # Now go back following the special C_n rules
        for a in reversed(range(tableau_height, n - 1)):
            if case_S[a] == max_width:
                self._insert_cell_case_S(self.ret_rig_con[a])
            else:
                max_width = self.ret_rig_con[a].insert_cell(max_width)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged partitions
        if tableau_height < n:
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

        if pos_val <= tableau_height:
            for a in range(pos_val-1, tableau_height):
                self._update_vacancy_nums(a)
                self._update_partition_values(a)
            if pos_val > 1:
                self._update_vacancy_nums(pos_val - 2)
                self._update_partition_values(pos_val - 2)
        elif tableau_height > 0:
            self._update_vacancy_nums(tableau_height - 1)
            self._update_partition_values(tableau_height - 1)

    def _insert_cell_case_S(self, partition):
        """
        Insert a cell when case `(S)` holds.

        TESTS::

            sage: RC = RiggedConfigurations(['C', 2, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2]])[1]
            sage: RP
            -4[ ][ ]-4
            -4[ ][ ]-4
            <BLANKLINE>
            sage: RP.rigging[0] = None
            sage: from sage.combinat.rigged_configurations.bij_type_C import KRTToRCBijectionTypeC
            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['C', 3, 1], [[2,1]])
            sage: bijection = KRTToRCBijectionTypeC(KRT(pathlist=[[-1,2]]))
            sage: bijection._insert_cell_case_S(RP)
            sage: RP
            -4[ ][ ][ ]None
            -4[ ][ ]-4
            <BLANKLINE>
        """
        # Special case when adding twice to the first row
        if partition.rigging[0] is None:
            partition._list[0] += 1
            return

        num_rows = len(partition)
        for i in reversed(range(1, num_rows)):
            if partition.rigging[i] is None:
                j = i - 1
                while j >= 0 and partition._list[j] == partition._list[i]:
                    partition.rigging[j+1] = partition.rigging[j]  # Shuffle it along
                    j -= 1
                partition._list[j+1] += 1
                partition.rigging[j+1] = None
                return


class RCToKRTBijectionTypeC(RCToKRTBijectionTypeA):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `C_n^{(1)}`.
    """

    def next_state(self, height):
        r"""
        Build the next state for type `C_n^{(1)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['C', 3, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_C import RCToKRTBijectionTypeC
            sage: bijection = RCToKRTBijectionTypeC(RC(partition_list=[[2],[2],[1]]))
            sage: bijection.next_state(1)
            -1
        """
        height -= 1  # indexing
        n = self.n
        ell = [None] * (2*n)
        case_S = [False] * n
        b = None

        # Calculate the rank and ell values

        last_size = 0
        for a in range(height, n-1):
            ell[a] = self._find_singular_string(self.cur_partitions[a], last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = self.cur_partitions[a][ell[a]]

        # Special case for n
        if b is None:
            # Since we are dividing by 2, we can use the identity of
            #   ceiling = floor + remainder
            ell[n-1] = self._find_singular_string(self.cur_partitions[n-1],
                                                 (last_size // 2) + (last_size % 2))

            if ell[n-1] is None:
                b = n
            else:
                last_size = self.cur_partitions[n-1][ell[n-1]] * 2

        if b is None:
            # Now go back
            ell[2*n-1] = ell[n-1]
            case_S[n-1] = True
            for a in reversed(range(n-1)):
                if a >= height and self.cur_partitions[a][ell[a]] == last_size:
                    ell[n+a] = ell[a]
                    case_S[a] = True
                else:  # note last_size > 1
                    ell[n+a] = self._find_singular_string(self.cur_partitions[a], last_size)

                    if ell[n + a] is None:
                        b = -(a + 2)
                        break
                    else:
                        last_size = self.cur_partitions[a][ell[n + a]]

        if b is None:
            b = -1

        # Determine the new rigged configuration by removing boxes from the
        #   selected string and then making the new string singular
        if case_S[0]:
            row_num = self.cur_partitions[0].remove_cell(ell[0], 2)
            row_num_bar = None
        else:
            row_num = self.cur_partitions[0].remove_cell(ell[0])
            row_num_bar = self.cur_partitions[0].remove_cell(ell[n])
        for a in range(1, n-1):
            if case_S[a]:
                row_num_next = self.cur_partitions[a].remove_cell(ell[a], 2)
                row_num_bar_next = None
            else:
                row_num_next = self.cur_partitions[a].remove_cell(ell[a])
                row_num_bar_next = self.cur_partitions[a].remove_cell(ell[n+a])

            self._update_vacancy_numbers(a - 1)
            if row_num is not None:
                self.cur_partitions[a-1].rigging[row_num] = self.cur_partitions[a-1].vacancy_numbers[row_num]
            if row_num_bar is not None:
                self.cur_partitions[a-1].rigging[row_num_bar] = self.cur_partitions[a-1].vacancy_numbers[row_num_bar]
            row_num = row_num_next
            row_num_bar = row_num_bar_next

        row_num_next = self.cur_partitions[n-1].remove_cell(ell[n-1])

        self._update_vacancy_numbers(n - 2)
        if row_num is not None:
            self.cur_partitions[n-2].rigging[row_num] = self.cur_partitions[n-2].vacancy_numbers[row_num]
        if row_num_bar is not None:
            self.cur_partitions[n-2].rigging[row_num_bar] = self.cur_partitions[n-2].vacancy_numbers[row_num_bar]

        self._update_vacancy_numbers(n - 1)
        if row_num_next is not None:
            self.cur_partitions[n-1].rigging[row_num_next] = self.cur_partitions[n-1].vacancy_numbers[row_num_next]

        return b
