r"""
Rigged Partitions

Class and methods of the rigged partition which are used by the rigged
configuration class. This is an internal class used by the rigged
configurations and KR tableaux during the bijection, and is not to be used by
the end-user.

We hold the partitions as an 1-dim array of positive integers where each
value corresponds to the length of the row. This is the shape of the
partition which can be accessed by the regular index.

The data for the vacancy number is also stored in a 1-dim array which each
entry corresponds to the row of the tableau, and similarly for the
partition values.

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

.. TODO::

    Convert this to using multiplicities `m_i` (perhaps with a dictionary?)?
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.combinat import CombinatorialObject
from sage.misc.latex import latex

class RiggedPartition(CombinatorialObject):
    r"""
    The RiggedPartition class which is the data structure of a rigged (i.e.
    marked or decorated) Young diagram of a partition.

    Note that this class as a stand-alone object does not make sense since the
    vacancy numbers are calculated using the entire rigged configuration. For
    more, see :class:`RiggedConfigurations`.

    EXAMPLES::

        sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
        sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
        sage: RP
         0[ ][ ]0
        -1[ ]-1
        <BLANKLINE>
        sage: type(RP)
        <class 'sage.combinat.rigged_configurations.rigged_partition.RiggedPartition'>
    """

    def __init__(self, shape=None, rigging_list=None, vacancy_nums=None):
        r"""
        Initialize by the rigged partition.

        Note that this only performs checks to see that the sizes match up.

        INPUT:

        - ``shape`` -- (Default: ``None``) The shape

        - ``rigging_list`` -- (Default: ``None``) The riggings

        - ``vacancy_nums`` -- (Default: ``None``) The vacancy numbers

        TESTS::

            sage: sage.combinat.rigged_configurations.rigged_partition.RiggedPartition()
            (/)
            <BLANKLINE>
            sage: RP = sage.combinat.rigged_configurations.rigged_partition.RiggedPartition([2,1], [0,0], [1, 0])
            sage: RP
            1[ ][ ]0
            0[ ]0
            <BLANKLINE>
            sage: TestSuite(RP).run()
        """
        if shape is None:
            CombinatorialObject.__init__(self, [])
            self.vacancy_numbers = []
            self.rigging = []
        else:
            CombinatorialObject.__init__(self, shape)

            if vacancy_nums is not None:
                if len(shape) != len(vacancy_nums):
                    raise ValueError("Mismatch between shape and vacancy numbers")

                self.vacancy_numbers = list(vacancy_nums)
            else:
                self.vacancy_numbers = [None] * len(shape)

            if rigging_list is not None:

                if len(shape) != len(rigging_list):
                    raise ValueError("Mismatch between shape and rigging list")

                self.rigging = list(rigging_list)
            else:
                self.rigging = [None] * len(shape)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: elt =  RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
            sage: elt
             0[ ][ ]0
            -1[ ]-1
            <BLANKLINE>
            sage: Partitions.global_options(convention="french")
            sage: elt
            -1[ ]-1
             0[ ][ ]0
            <BLANKLINE>
            sage: Partitions.global_options.reset()
        """
        # If it is empty, return saying so
        if len(self._list) == 0:
            return("(/)\n")

        from sage.combinat.partition import Partitions
        if Partitions.global_options("convention") == "French":
            itr = reversed(list(enumerate(self._list)))
        else:
            itr = enumerate(self._list)
        ret_str = ""
        vac_num_width = max(len(str(vac_num)) for vac_num in self.vacancy_numbers)
        for i, val in itr:
            ret_str += ("{:>" + str(vac_num_width) + "}").format(self.vacancy_numbers[i])
            ret_str += "[ ]"*val
            ret_str += str(self.rigging[i])
            ret_str += "\n"
        return(ret_str)

    def _latex_(self):
        r"""
        Returns LaTeX representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: latex(RC(partition_list=[[2],[2,2],[2,1],[2]])[2])
            {
            \begin{array}[t]{r|c|c|l}
             \cline{2-3} 0 &\phantom{|}&\phantom{|}& 0 \\
             \cline{2-3} -1 &\phantom{|}& \multicolumn{2 }{l}{ -1 } \\
             \cline{2-2} 
            \end{array}
            }

        TESTS:

        Check that this prints using the French convention::

            sage: RC = RiggedConfigurations(['D',5,1], [[2,1], [1,2]])
            sage: RiggedConfigurations.global_options(convention='French')
            sage: latex(RC(partition_list=[[3],[3,1],[1,1],[1],[1]])[1])
            {
            \begin{array}[t]{r|c|c|c|l}
             \cline{2-2} 0 &\phantom{|}& \multicolumn{3 }{l}{ 0 } \\
             \cline{2-4} -2 &\phantom{|}&\phantom{|}&\phantom{|}& -2 \\
             \cline{2-4}
            \end{array}
            }
            sage: RiggedConfigurations.global_options.reset()
        """
        num_rows = len(self._list)
        if num_rows == 0:
            return "{\\emptyset}"

        num_cols = self._list[0]
        ret_string = "{\n\\begin{array}[t]{r|" + "c|"*num_cols + "l}\n"

        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurationOptions
        if RiggedConfigurationOptions['convention'] == 'English':
            ret_string += "\\cline{2-%s} "%(1+num_cols) + latex(self.vacancy_numbers[0])
            for i, row_len in enumerate(self._list):

                ret_string += " &" + "\\phantom{|}&"*row_len

                if num_cols == row_len:
                    ret_string += " " + latex(self.rigging[i])
                else:
                    ret_string += " \\multicolumn{" + repr(num_cols - row_len + 1)
                    ret_string += "}{l}{" + latex(self.rigging[i]) + "}"

                ret_string += " \\\\\n"

                ret_string += "\\cline{2-" + repr(1 + row_len) + "} "
                if i != num_rows - 1 and row_len != self._list[i + 1]:
                    ret_string += latex(self.vacancy_numbers[i + 1])
            ret_string += "\n\\end{array}\n}"
        else:
            for i, row_len in enumerate(reversed(self._list)):
                ret_string += "\\cline{2-%s} "%(1 + row_len) + latex(self.vacancy_numbers[-i-1])
                ret_string += " &" + "\\phantom{|}&"*row_len

                if num_cols == row_len:
                    ret_string += " " + latex(self.rigging[-i-1])
                else:
                    ret_string += " \\multicolumn{" + repr(num_cols - row_len + 1)
                    ret_string += "}{l}{" + latex(self.rigging[-i-1]) + "}"

                ret_string += " \\\\\n"
            ret_string += "\\cline{2-%s}\n\\end{array}\n}"%(1 + num_cols)

        return ret_string

    def _clone(self):
        r"""
        Makes a (deep) copy of this rigged partition.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]; RP
             0[ ][ ]0
            -1[ ]-1
            <BLANKLINE>
            sage: RP2 = RP._clone(); RP2
             0[ ][ ]0
            -1[ ]-1
            <BLANKLINE>
            sage: RP == RP2
            True
            sage: RP is RP2
            False
        """
        return self.__class__(self._list[:], self.rigging[:], self.vacancy_numbers[:])

    def __eq__(self, rhs):
        r"""
        Return true if ``self`` equals ``rhs``.

        TESTS::

            sage: RC = RiggedConfigurations(['A',2,1], [[1,1],[1,1],[1,1]])
            sage: x = RC(partition_list=[[1], []], rigging_list=[[0], []])
            sage: y = RC(partition_list=[[1], []], rigging_list=[[1], []])
            sage: x == y
            False
        """
        if isinstance(rhs, RiggedPartition):
            return(self._list.__eq__(rhs._list) and self.rigging.__eq__(rhs.rigging))

        return False

    # Should we move these functions to the CP -> RC bijections?

    def get_num_cells_to_column(self, end_column, t=1):
        r"""
        Get the number of cells in all columns before the ``end_column``.

        INPUT:

        - ``end_column`` -- The index of the column to end at

        - ``t`` -- The scaling factor

        OUTPUT:

        - The number of cells

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
            sage: RP.get_num_cells_to_column(1)
            2
            sage: RP.get_num_cells_to_column(2)
            3
            sage: RP.get_num_cells_to_column(3)
            3
            sage: RP.get_num_cells_to_column(3, 2)
            5
        """
        sum_cells = 0
        # Sum up from the reverse (the smallest row sizes)
        i = len(self._list) - 1
        while i >= 0 and self._list[i]*t < end_column:
            sum_cells += self._list[i]*t
            i -= 1

        # Add the remaining cells
        if i > -1:
            sum_cells += end_column * (i + 1)

        return sum_cells

    def insert_cell(self, max_width):
        r"""
        Insert a cell given at a singular value as long as its less than the
        specified width.

        Note that :meth:`insert_cell` does not update riggings or vacancy
        numbers, but it does prepare the space for them. Returns the width of
        the row we inserted at.

        INPUT:

        - ``max_width`` -- The maximum width (i.e. row length) that we can
          insert the cell at

        OUTPUT:

        - The width of the row we inserted at.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
            sage: RP.insert_cell(2)
            2
            sage: RP
             0[ ][ ][ ]None
            -1[ ]-1
            <BLANKLINE>
        """
        max_pos = -1
        if max_width > 0:
            for i, vac_num in enumerate(self.vacancy_numbers):
                if self._list[i] <= max_width and vac_num == self.rigging[i]:
                    max_pos = i
                    break

        if max_pos == -1: # No singular values, then add a new row
            self._list.append(1)
            self.vacancy_numbers.append(None)
            # Go through our partition until we find a length of greater than 1
            i = len(self._list) - 1
            while i >= 0 and self._list[i] == 1:
                i -= 1
            self.rigging.insert(i + 1, None)
            return 0

        self._list[max_pos] += 1
        self.rigging[max_pos] = None # State that we've changed this row
        return self._list[max_pos] - 1

    def remove_cell(self, row, num_cells=1):
        r"""
        Removes a cell at the specified ``row``.

        Note that :meth:`remove_cell` does not set/update the vacancy numbers
        or the riggings, but guarantees that the location has been allocated
        in the returned index.

        INPUT:

        - ``row`` -- The row to remove the cell from

        - ``num_cells`` -- (Default: 1) The number of cells to remove

        OUTPUT:

        - The location of the newly constructed row or ``None`` if unable to
          remove row or if deleted a row.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
            sage: RP.remove_cell(0)
            0
            sage: RP
             0[ ]0
            -1[ ]-1
            <BLANKLINE>
        """
        if row is None:
            return None

        if self._list[row] <= num_cells:
            self._list.pop(row)
            self.vacancy_numbers.pop(row)
            self.rigging.pop(row)
            return None

        # Find the beginning of the next block we want
        block_len = self._list[row] - num_cells # The length of the desired block
        if row + 1 == len(self._list):
            # If we are at the end, just do a simple remove
            self._list[row] = block_len
            return row

        for i in range(row + 1, len(self._list)):
            if self._list[i] <= block_len:
                if i == row + 1:
                    # If the next row is a block change, just reduce by num_cells
                    self._list[row] = block_len
                    return row

                # Otherwise we need to "move" the row
                self._list.insert(i, block_len)
                # These should be updated (so there should be no need to carry them over)
                self.vacancy_numbers.insert(i, None)
                self.rigging.insert(i, None)

                self._list.pop(row)
                self.vacancy_numbers.pop(row)
                self.rigging.pop(row)
                return i - 1

        # We need to "move" the row to the end of the partition
        self._list.pop(row)
        self.vacancy_numbers.pop(row)
        self.rigging.pop(row)

        self._list.append(block_len)
        # Placeholders as above
        self.vacancy_numbers.append(None)
        self.rigging.append(None)
        return len(self._list) - 1

class RiggedPartitionTypeB(RiggedPartition):
    r"""
    Rigged partitions for type `B_n^{(1)}` which has special printing rules
    which comes from the fact that the `n`-th partition can have columns of
    width `\frac{1}{2}`.
    """
    def __init__(self, arg0, arg1=None, arg2=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RP = sage.combinat.rigged_configurations.rigged_partition.RiggedPartition([2,1], [0,0], [1, 0])
            sage: B = sage.combinat.rigged_configurations.rigged_partition.RiggedPartitionTypeB(RP); B
            1[][]0
            0[]0
            <BLANKLINE>
            sage: TestSuite(B).run()
        """
        if arg1 is not None:
            RiggedPartition.__init__(self, arg0, arg1, arg2)
            return

        RiggedPartition.__init__(self,
                                 arg0._list,
                                 arg0.rigging,
                                 arg0.vacancy_numbers)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        INPUT:

        - ``half_width_boxes`` -- (Default: ``True``) Display the partition
          using half width boxes

        EXAMPLES::

            sage: RC = RiggedConfigurations(['B', 2, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[2],[2,1]])[1]
            sage: elt
            -2[][]-2
            -2[]-2
            <BLANKLINE>
            sage: RiggedConfigurations.global_options(half_width_boxes_type_B=False)
            sage: elt
            -2[ ][ ]-2
            -2[ ]-2
            <BLANKLINE>
            sage: RiggedConfigurations.global_options.reset()
        """
        # If it is empty, return saying so
        if len(self._list) == 0:
            return("(/)\n")

        from sage.combinat.partition import Partitions
        if Partitions.global_options("convention") == "french":
            itr = reversed(list(enumerate(self._list)))
        else:
            itr = enumerate(self._list)
        ret_str = ""

        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurationOptions
        if RiggedConfigurationOptions['half_width_boxes_type_B']:
            box_str = "[]"
        else:
            box_str = "[ ]"

        vac_num_width = max(len(str(vac_num)) for vac_num in self.vacancy_numbers)
        for i, val in itr:
            ret_str += ("{:>" + str(vac_num_width) + "}").format(self.vacancy_numbers[i])
            ret_str += box_str*val
            ret_str += str(self.rigging[i])
            ret_str += "\n"
        return(ret_str)

    def _latex_(self):
        r"""
        Returns LaTeX representation of ``self``.

        INPUT:

        - ``half_width_boxes`` -- (Default: ``True``) Display the partition
          using half width boxes

        EXAMPLES::

            sage: RC = RiggedConfigurations(['B', 2, 1], [[1, 1]])
            sage: RP = RC(partition_list=[[],[2]])[1]
            sage: latex(RP)
            {
            \begin{array}[t]{r|c|c|l}
             \cline{2-3} -4 &\phantom{a}&\phantom{a}& -4 \\
             \cline{2-3} 
            \end{array}
            }
            sage: RiggedConfigurations.global_options(half_width_boxes_type_B=False)
            sage: latex(RP)
            {
            \begin{array}[t]{r|c|c|l}
             \cline{2-3} -4 &\phantom{X|}&\phantom{X|}& -4 \\
             \cline{2-3} 
            \end{array}
            }
            sage: RiggedConfigurations.global_options.reset()
        """
        num_rows = len(self._list)
        if num_rows == 0:
            return "{\\emptyset}"
        
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurationOptions
        if RiggedConfigurationOptions['half_width_boxes_type_B']:
            box_str = "\\phantom{a}&"
        else:
            box_str = "\\phantom{X|}&"

        num_cols = self._list[0]
        ret_string = "{\n\\begin{array}[t]{r|" + "c|"*num_cols + "l}\n"

        if RiggedConfigurationOptions['convention'] == 'English':
            ret_string += "\\cline{2-%s} "%(1+num_cols) + latex(self.vacancy_numbers[0])
            for i, row_len in enumerate(self._list):
                ret_string += " &" + box_str*row_len

                if num_cols == row_len:
                    ret_string += " " + latex(self.rigging[i])
                else:
                    ret_string += " \\multicolumn{" + repr(num_cols - row_len + 1)
                    ret_string += "}{l}{" + latex(self.rigging[i]) + "}"

                ret_string += " \\\\\n"

                ret_string += "\\cline{2-" + repr(1 + row_len) + "} "
                if i != num_rows - 1 and row_len != self._list[i + 1]:
                    ret_string += latex(self.vacancy_numbers[i + 1])
            ret_string += "\n\\end{array}\n}"
        else:
            for i, row_len in enumerate(reversed(self._list)):
                ret_string += "\\cline{2-%s} "%(1 + row_len)
                ret_string += latex(self.vacancy_numbers[-i-1])
                ret_string += " &" + box_str*row_len

                if num_cols == row_len:
                    ret_string += " " + latex(self.rigging[-i-1])
                else:
                    ret_string += " \\multicolumn{" + repr(num_cols - row_len + 1)
                    ret_string += "}{l}{" + latex(self.rigging[-i-1]) + "}"

                ret_string += " \\\\\n"
            ret_string += "\\cline{2-%s}\n\\end{array}\n}"%(1 + num_cols)

        return ret_string

