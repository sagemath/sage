# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection classes for type `A_{2n}^{(2)}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `A_{2n}^{(2)}`.

AUTHORS:

- Travis Scrimshaw (2012-12-21): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 2], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A2_even import KRTToRCBijectionTypeA2Even
    sage: bijection = KRTToRCBijectionTypeA2Even(KRT(pathlist=[[-1,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['A', 4, 2], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A2_even import RCToKRTBijectionTypeA2Even
    sage: bijection = RCToKRTBijectionTypeA2Even(RC(partition_list=[[],[]]))
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
from sage.combinat.rigged_configurations.bij_type_C import KRTToRCBijectionTypeC
from sage.combinat.rigged_configurations.bij_type_C import RCToKRTBijectionTypeC


class KRTToRCBijectionTypeA2Even(KRTToRCBijectionTypeC):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `A_{2n}^{(2)}`.

    This inherits from type `C_n^{(1)}` because we use the same methods in
    some places.
    """

    def next_state(self, val):
        r"""
        Build the next state for type `A_{2n}^{(2)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 2], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A2_even import KRTToRCBijectionTypeA2Even
            sage: bijection = KRTToRCBijectionTypeA2Even(KRT(pathlist=[[-1,-2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [-2])
            sage: bijection.next_state(-2)
        """
        # Note that we must subtract 1 from n to match the indices.
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        # If it is a regular value, we follow the A_n rules
        if val != 'E' and val > 0:
            KRTToRCBijectionTypeA.next_state(self, val)
            return

        case_S = [None] * n
        if val == 'E':
            pos_val = n
            max_width = self.ret_rig_con[n-1].insert_cell(0)
        else:
            pos_val = -val

            # Always add a cell to the first singular value in the first
            #   tableau we are updating.
            if len(self.ret_rig_con[pos_val - 1]) > 0:
                max_width = self.ret_rig_con[pos_val - 1][0]
            else:
                max_width = 1

            # Add cells similar to type A_n but we move to the right until we
            #   reach the value of n
            for a in range(pos_val - 1, n):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                case_S[a] = max_width

            # Special case for n
            self._insert_cell_case_S(self.ret_rig_con[n-1])

        # Now go back following the regular C_n (ish) rules
        for a in reversed(range(tableau_height, n-1)):
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


class RCToKRTBijectionTypeA2Even(RCToKRTBijectionTypeC):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `A_{2n}^{(2)}`.
    """

    def next_state(self, height):
        r"""
        Build the next state for type `A_{2n}^{(2)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 2], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A2_even import RCToKRTBijectionTypeA2Even
            sage: bijection = RCToKRTBijectionTypeA2Even(RC(partition_list=[[2],[2,2]]))
            sage: bijection.next_state(2)
            -1
        """
        height -= 1  # indexing
        n = self.n
        ell = [None] * (2*n)
        case_S = [False] * n
        b = None

        # Calculate the rank and ell values

        last_size = 0
        for a in range(height, n):
            ell[a] = self._find_singular_string(self.cur_partitions[a], last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = self.cur_partitions[a][ell[a]]

        if b is None and last_size == 1:
            b = 'E'
            # This is a slight hack since remove_cell() will just delete the
            #   appropriate row
            case_S[n-1] = True

        if b is None:
            # Now go back
            ell[2*n-1] = ell[n-1]
            case_S[n-1] = True
            for a in reversed(range(n-1)):
                if a >= height and self.cur_partitions[a][ell[a]] == last_size:
                    ell[n+a] = ell[a]
                    case_S[a] = True
                else:
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
        for a in range(1, n):
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

        self._update_vacancy_numbers(n - 1)
        if row_num is not None:
            self.cur_partitions[n-1].rigging[row_num] = self.cur_partitions[n-1].vacancy_numbers[row_num]
        if row_num_bar is not None:
            self.cur_partitions[n-1].rigging[row_num_bar] = self.cur_partitions[n-1].vacancy_numbers[row_num_bar]

        return b
