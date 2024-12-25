# sage_setup: distribution = sagemath-repl
"""
Utility functions for pretty-printing

These utility functions are used in the implementations of ``_repr_``
methods elsewhere.
"""

# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


class TallListFormatter:
    """
    Special representation for lists with tall entries (e.g. matrices).

    .. automethod:: __call__
    """

    # This is used to wrap lines when printing "tall" lists.
    MAX_COLUMN = 70

    def _tall_list_row(self, running_lines, last_row=False):
        """
        Helper for :meth:`_check_tall_list_and_format`.

        This helper function processes and outputs the contents of the
        running_lines array.

        TESTS::

            sage: from sage.repl.display.util import format_list
            sage: format_list._tall_list_row(['a   b', 'b  c', 'c'])
            ['a           b', 'b        c', 'c,', '']
        """
        s = []
        for i, line in enumerate(running_lines):
            if i + 1 != len(running_lines):
                sep, tail = '  ', ''
            else:
                # The commas go on the bottom line of this row.
                sep, tail = ', ', '' if last_row else ','
            s.append(sep.join(line) + tail)
        # Separate rows with a newline to make them stand out.
        if not last_row:
            s.append("")
        return s

    def try_format(self, the_list):
        """
        First check whether a list is "tall" -- whether the reprs of the
        elements of the list will span multiple lines and cause the list
        to be printed awkwardly.  If not, this function returns ``None`` and
        does nothing; you should revert back to the normal method for
        printing an object (its repr). If so, return the string in the
        special format. Note that the special format isn't just for
        matrices. Any object with a multiline repr will be formatted.

        INPUT:

        - ``the_list`` -- the list (or a tuple)

        OUTPUT:

        String or ``None``. The latter is returned if the list is not
        deemed to be tall enough and another formatter should be used.

        TESTS::

            sage: from sage.repl.display.util import format_list
            sage: print(format_list.try_format(                                         # needs sage.modules
            ....:        [matrix([[1, 2, 3, 4], [5, 6, 7, 8]]) for i in range(7)]))
            [
            [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]
            [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8],
            <BLANKLINE>
            [1 2 3 4]
            [5 6 7 8]
            ]

            sage: format_list.try_format(['not', 'tall']) is None
            True
        """
        # For every object to be printed, split its repr on newlines and store the
        # result in this list.
        split_reprs = []
        tall = False
        for elem in the_list:
            split_reprs.append(repr(elem).split('\n'))
            if len(split_reprs[-1]) > 1:
                # Meanwhile, check to make sure the list is actually "tall".
                tall = True
        if not tall:
            return None
        # Figure out which type of parenthesis to use, based on the type of the_list.
        if isinstance(the_list, tuple):
            parens = '()'
        elif isinstance(the_list, list):
            parens = '[]'
        else:
            raise TypeError('expected list or tuple')

        # running_lines is a list of lines, which are stored as lists of strings
        # to be joined later. For each split repr, we add its lines to the
        # running_lines array. When current_column exceeds MAX_COLUMN, process
        # and output running_lines using _print_tall_list_row.
        running_lines = [[]]
        current_column = 0
        s = [parens[0]]
        for split_repr in split_reprs:
            width = max(len(x) for x in split_repr)
            if current_column + width > self.MAX_COLUMN and not (width > self.MAX_COLUMN):
                s.extend(self._tall_list_row(running_lines))
                running_lines = [[]]
                current_column = 0
            current_column += width + 2
            # Add the lines from split_repr to the running_lines array. It may
            # be necessary to add or remove lines from either one so that the
            # number of lines matches up.
            for i in range(len(running_lines), len(split_repr)):
                running_lines.insert(0, [' ' * len(x) for x in running_lines[-1]])
            line_diff = len(running_lines) - len(split_repr)
            for i, x in enumerate(split_repr):
                running_lines[i + line_diff].append(x.ljust(width))
            for i in range(line_diff):
                running_lines[i].append(' ' * width)
        # Output any remaining entries.
        if len(running_lines[0]) > 0:
            s.extend(self._tall_list_row(running_lines, True))
        s.append(parens[1])
        return "\n".join(s)

    def __call__(self, the_list):
        """
        Return "tall list" string representation.

        See also :meth:`try_format`.

        INPUT:

        - ``the_list`` -- list or tuple

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.display.util import format_list
            sage: format_list(['not', 'tall'])
            "['not', 'tall']"
        """
        output = self.try_format(the_list)
        if output is None:
            output = repr(the_list)
        return output


format_list = TallListFormatter()
