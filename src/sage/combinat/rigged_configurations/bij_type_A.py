r"""
Bijection classes for type `A_n^{(1)}`

Part of the (internal) classes which run the bijection between rigged
configurations and tensor products of Kirillov-Reshetikhin tableaux of
type `A_n^{(1)}`.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
    sage: bijection = KRTToRCBijectionTypeA(KRT(pathlist=[[5,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
    sage: bijection = RCToKRTBijectionTypeA(RC(partition_list=[[],[],[],[]]))
    sage: TestSuite(bijection).run()
"""

#*****************************************************************************
#       Copyright (C) 2011, 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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
#*****************************************************************************

from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract

class KRTToRCBijectionTypeA(KRTToRCBijectionAbstract):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `A_n^{(1)}`.
    """

    def next_state(self, val):
        r"""
        Build the next state for type `A_n^{(1)}`.

        EXAMPLES::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
            sage: bijection = KRTToRCBijectionTypeA(KRT(pathlist=[[4,3]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [3])
            sage: bijection.next_state(3)
        """
        tableau_height = len(self.cur_path[0]) - 1
        n = self.n

        # Note first we subtract off for the n = max value (in the path) - 1,
        #   then we remove 1 to match the indices between math and programming.
        if val - 1 > tableau_height:
            # Always add a cell to the first singular value in the first
            #   tableau we are updating.
            if len(self.ret_rig_con[val - 2]) > 0:
                max_width = self.ret_rig_con[val - 2][0]
            else:
                max_width = 1

            # Insert a cell into the rightmost rigged partition
            max_width = self.ret_rig_con[val - 2].insert_cell(max_width)

            # Move to the left and update values as we have finished modifying
            #   everything which affects its vacancy/partition values
            for a in reversed(range(tableau_height, val - 2)):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                self._update_vacancy_nums(a + 1)
                self._update_partition_values(a + 1)

            # Update the final rigged tableau
            # Note if tabelauHeight = n+1, then we must have val = n+1 (in order
            #   to be column strict increasing), but then tableau_height is never
            #   greater than val, so we don't enter into this statement.
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

            if val - 1 < n:
                self._update_vacancy_nums(val - 1)

            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height - 1)
        elif tableau_height - 1 < n:
            # Otherwise we just need to update the vacancy numbers that are affected
            if tableau_height < n:
                self._update_vacancy_nums(tableau_height)

            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height - 1)

class RCToKRTBijectionTypeA(RCToKRTBijectionAbstract):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `A_n^{(1)}`.
    """

    def next_state(self, height):
        r"""
        Build the next state for type `A_n^{(1)}`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
            sage: bijection = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection.next_state(1)
            5
        """
        height -= 1 # indexing
        n = self.n
        ell = [None] * n
        b = None

        # Calculate the rank and ell values
        last_size = 0
        a = height
        for partition in self.cur_partitions[height:]:
            ell[a] = self._find_singular_string(partition, last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = partition[ell[a]]
            a += 1

        if b is None:
            b = n + 1

        # Determine the new rigged configuration by removing a box from the selected
        #   string and then making the new string singular
        row_num = self.cur_partitions[0].remove_cell(ell[0])
        for a in range(1, n):
            row_num_next = self.cur_partitions[a].remove_cell(ell[a])

            self._update_vacancy_numbers(a - 1)
            if row_num is not None:
                self.cur_partitions[a - 1].rigging[row_num] = self.cur_partitions[a - 1].vacancy_numbers[row_num]
            row_num = row_num_next

        self._update_vacancy_numbers(n - 1)
        if row_num is not None:
            self.cur_partitions[n - 1].rigging[row_num] = self.cur_partitions[n - 1].vacancy_numbers[row_num]

        return(b)

