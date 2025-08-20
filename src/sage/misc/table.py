r"""
Tables

Display a rectangular array as a table, either in plain text, LaTeX,
or html.  See the documentation for :class:`table` for details and
examples.

AUTHORS:

- John H. Palmieri (2012-11)
"""

from io import StringIO

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method


class table(SageObject):
    r"""
    Display a rectangular array as a table, either in plain text, LaTeX,
    or html.

    INPUT:

    - ``rows`` -- (default: ``None``) list of lists (or list of tuples,
      etc.), containing the data to be displayed
    - ``columns`` -- (default: ``None``) list of lists (etc.), containing
      the data to be displayed, but stored as columns. Set either ``rows``
      or ``columns``, but not both.
    - ``header_row`` -- (default: ``False``) if ``True``, first row is
      highlighted
    - ``header_column`` -- (default: ``False``) if ``True``, first column is
      highlighted
    - ``frame`` -- (default: ``False``) if ``True``, put a box around each
      cell
    - ``align`` -- (default: ``'left'``) the alignment of each entry: either
      ``'left'``, ``'center'``, or ``'right'``

    EXAMPLES::

        sage: rows = [['a', 'b', 'c'], [100,2,3], [4,5,60]]
        sage: table(rows)
          a     b   c
          100   2   3
          4     5   60
        sage: latex(table(rows))
        \begin{tabular}{lll}
        a & b & c \\
        $100$ & $2$ & $3$ \\
        $4$ & $5$ & $60$ \\
        \end{tabular}

    If ``header_row`` is ``True``, then the first row is highlighted. If
    ``header_column`` is ``True``, then the first column is
    highlighted. If ``frame`` is ``True``, then print a box around every
    "cell". ::

        sage: table(rows, header_row=True)
            a     b   c
          ├─────┼───┼────┤
            100   2   3
            4     5   60
        sage: latex(table(rows, header_row=True))
        \begin{tabular}{lll}
        a & b & c \\ \hline
        $100$ & $2$ & $3$ \\
        $4$ & $5$ & $60$ \\
        \end{tabular}
        sage: table(rows=rows, frame=True)
        ┌─────┬───┬────┐
        │ a   │ b │ c  │
        ├─────┼───┼────┤
        │ 100 │ 2 │ 3  │
        ├─────┼───┼────┤
        │ 4   │ 5 │ 60 │
        └─────┴───┴────┘
        sage: latex(table(rows=rows, frame=True))
        \begin{tabular}{|l|l|l|} \hline
        a & b & c \\ \hline
        $100$ & $2$ & $3$ \\ \hline
        $4$ & $5$ & $60$ \\ \hline
        \end{tabular}
        sage: table(rows, header_column=True, frame=True)
        ┌─────╥───┬────┐
        │ a   ║ b │ c  │
        ├─────╫───┼────┤
        │ 100 ║ 2 │ 3  │
        ├─────╫───┼────┤
        │ 4   ║ 5 │ 60 │
        └─────╨───┴────┘
        sage: latex(table(rows, header_row=True, frame=True))
        \begin{tabular}{|l|l|l|} \hline
        a & b & c \\ \hline \hline
        $100$ & $2$ & $3$ \\ \hline
        $4$ & $5$ & $60$ \\ \hline
        \end{tabular}
        sage: table(rows, header_column=True)
          a   │ b   c
          100 │ 2   3
          4   │ 5   60

    The argument ``header_row`` can, instead of being ``True`` or
    ``False``, be the contents of the header row, so that ``rows``
    consists of the data, while ``header_row`` is the header
    information.  The same goes for ``header_column``. Passing lists
    for both arguments simultaneously is not supported. ::

        sage: table([(x,n(sin(x), digits=2)) for x in [0..3]],                          # needs sage.symbolic
        ....:       header_row=["$x$", r"$\sin(x)$"], frame=True)
        ┌─────┬───────────┐
        │ $x$ │ $\sin(x)$ │
        ╞═════╪═══════════╡
        │ 0   │ 0.00      │
        ├─────┼───────────┤
        │ 1   │ 0.84      │
        ├─────┼───────────┤
        │ 2   │ 0.91      │
        ├─────┼───────────┤
        │ 3   │ 0.14      │
        └─────┴───────────┘

    You can create the transpose of this table in several ways, for
    example, "by hand," that is, changing the data defining the table::

        sage: table(rows=[[x for x in [0..3]],                                          # needs sage.symbolic
        ....:             [n(sin(x), digits=2) for x in [0..3]]],
        ....:       header_column=['$x$', r'$\sin(x)$'], frame=True)
        ┌───────────╥──────┬──────┬──────┬──────┐
        │ $x$       ║ 0    │ 1    │ 2    │ 3    │
        ├───────────╫──────┼──────┼──────┼──────┤
        │ $\sin(x)$ ║ 0.00 │ 0.84 │ 0.91 │ 0.14 │
        └───────────╨──────┴──────┴──────┴──────┘

    or by passing the original data as the ``columns`` of the table
    and using ``header_column`` instead of ``header_row``::

        sage: table(columns=[(x, n(sin(x), digits=2)) for x in [0..3]],                 # needs sage.symbolic
        ....:       header_column=['$x$', r'$\sin(x)$'], frame=True)
        ┌───────────╥──────┬──────┬──────┬──────┐
        │ $x$       ║ 0    │ 1    │ 2    │ 3    │
        ├───────────╫──────┼──────┼──────┼──────┤
        │ $\sin(x)$ ║ 0.00 │ 0.84 │ 0.91 │ 0.14 │
        └───────────╨──────┴──────┴──────┴──────┘

    or by taking the :meth:`transpose` of the original table::

        sage: table(rows=[(x, n(sin(x), digits=2)) for x in [0..3]],                    # needs sage.symbolic
        ....:       header_row=['$x$', r'$\sin(x)$'], frame=True).transpose()
        ┌───────────╥──────┬──────┬──────┬──────┐
        │ $x$       ║ 0    │ 1    │ 2    │ 3    │
        ├───────────╫──────┼──────┼──────┼──────┤
        │ $\sin(x)$ ║ 0.00 │ 0.84 │ 0.91 │ 0.14 │
        └───────────╨──────┴──────┴──────┴──────┘

    In either plain text or LaTeX, entries in tables can be aligned to the
    left (default), center, or right::

        sage: table(rows, align='left')
          a     b   c
          100   2   3
          4     5   60
        sage: table(rows, align='center')
          a    b   c
         100   2   3
          4    5   60
        sage: table(rows, align='right', frame=True)
        ┌─────┬───┬────┐
        │   a │ b │  c │
        ├─────┼───┼────┤
        │ 100 │ 2 │  3 │
        ├─────┼───┼────┤
        │   4 │ 5 │ 60 │
        └─────┴───┴────┘

    To generate HTML you should use ``html(table(...))``::

        sage: # needs sage.symbolic
        sage: data = [["$x$", r"$\sin(x)$"]] + [(x, n(sin(x), digits=2))
        ....:                                   for x in [0..3]]
        sage: output = html(table(data, header_row=True, frame=True))
        sage: type(output)
        <class 'sage.misc.html.HtmlFragment'>
        sage: print(output)
        <div class="notruncate">
        <table border="1" class="table_form">
        <tbody>
        <tr>
        <th style="text-align:left">\(x\)</th>
        <th style="text-align:left">\(\sin(x)\)</th>
        </tr>
        <tr class ="row-a">
        <td style="text-align:left">\(0\)</td>
        <td style="text-align:left">\(0.00\)</td>
        </tr>
        <tr class ="row-b">
        <td style="text-align:left">\(1\)</td>
        <td style="text-align:left">\(0.84\)</td>
        </tr>
        <tr class ="row-a">
        <td style="text-align:left">\(2\)</td>
        <td style="text-align:left">\(0.91\)</td>
        </tr>
        <tr class ="row-b">
        <td style="text-align:left">\(3\)</td>
        <td style="text-align:left">\(0.14\)</td>
        </tr>
        </tbody>
        </table>
        </div>

    It is an error to specify both ``rows`` and ``columns``::

        sage: table(rows=[[1,2,3], [4,5,6]], columns=[[0,0,0], [0,0,1024]])
        Traceback (most recent call last):
        ...
        ValueError: do not set both 'rows' and 'columns' when defining a table

        sage: table(columns=[[0,0,0], [0,0,1024]])
        0 0
        0 0
        0 1024

    Note that if ``rows`` is just a list or tuple, not nested, then
    it is treated as a single row::

        sage: table([1,2,3])
        1   2   3

    Also, if you pass a non-rectangular array, the longer rows or
    columns get truncated::

        sage: table([[1,2,3,7,12], [4,5]])
        1   2
        4   5
        sage: table(columns=[[1,2,3], [4,5,6,7]])
        1   4
        2   5
        3   6

    TESTS::

        sage: TestSuite(table([["$x$", r"$\sin(x)$"]] +                                 # needs sage.symbolic
        ....:                  [(x, n(sin(x), digits=2)) for x in [0..3]],
        ....:                 header_row=True, frame=True)).run()

    .. automethod:: _rich_repr_
    """
    def __init__(self, rows=None, columns=None, header_row=False,
                 header_column=False, frame=False, align='left'):
        r"""
        EXAMPLES::

            sage: table([1,2,3], frame=True)
            ┌───┬───┬───┐
            │ 1 │ 2 │ 3 │
            └───┴───┴───┘
        """
        # If both rows and columns are set, raise an error.
        if rows and columns:
            raise ValueError("do not set both 'rows' and 'columns' when defining a table")
        # If columns is set, use its transpose for rows.
        if columns:
            rows = list(zip(*columns))
        # Set the rest of the options.
        self._options = {}
        if header_row is True:
            self._options['header_row'] = True
        elif header_row:
            self._options['header_row'] = True
            rows = [header_row] + rows
        else:
            self._options['header_row'] = False
        if header_column is True:
            self._options['header_column'] = True
        elif header_column:
            self._options['header_column'] = True
            rows = [(a,) + tuple(x) for (a, x) in zip(header_column, rows)]
        else:
            self._options['header_column'] = False

        self._options['frame'] = frame
        self._options['align'] = align
        # Store rows as a tuple.
        if not isinstance(rows[0], (list, tuple)):
            rows = (rows,)
        self._rows = tuple(rows)

    def __eq__(self, other):
        r"""
        Two tables are equal if and only if their data rowss and
        their options are the same.

        EXAMPLES::

            sage: # needs sage.modules sage.plot sage.symbolic
            sage: rows = [['a', 'b', 'c'], [1,plot(sin(x)),3], [4,5,identity_matrix(2)]]
            sage: T = table(rows, header_row=True)
            sage: T2 = table(rows, header_row=True)
            sage: T is T2
            False
            sage: T == T2
            True
            sage: T2.options(frame=True)
            sage: T == T2
            False
        """
        return (self._rows == other._rows and self.options() == other.options())

    def options(self, **kwds):
        r"""
        With no arguments, return the dictionary of options for this
        table. With arguments, modify options.

        INPUT:

        - ``header_row`` -- if ``True``, first row is highlighted
        - ``header_column`` -- if ``True``, first column is highlighted
        - ``frame`` -- if ``True``, put a box around each cell
        - ``align`` -- the alignment of each entry: either ``'left'``,
          ``'center'``, or ``'right'``

        EXAMPLES::

            sage: T = table([['a', 'b', 'c'], [1,2,3]])
            sage: T.options()['align'], T.options()['frame']
            ('left', False)
            sage: T.options(align='right', frame=True)
            sage: T.options()['align'], T.options()['frame']
            ('right', True)

        Note that when first initializing a table, ``header_row`` or
        ``header_column`` can be a list. In this case, during the
        initialization process, the header is merged with the rest of
        the data, so changing the header option later using
        ``table.options(...)`` doesn't affect the contents of the
        table, just whether the row or column is highlighted. When
        using this :meth:`options` method, no merging of data occurs,
        so here ``header_row`` and ``header_column`` should just be
        ``True`` or ``False``, not a list. ::

            sage: T = table([[1,2,3], [4,5,6]], header_row=['a', 'b', 'c'], frame=True)
            sage: T
            ┌───┬───┬───┐
            │ a │ b │ c │
            ╞═══╪═══╪═══╡
            │ 1 │ 2 │ 3 │
            ├───┼───┼───┤
            │ 4 │ 5 │ 6 │
            └───┴───┴───┘
            sage: T.options(header_row=False)
            sage: T
            ┌───┬───┬───┐
            │ a │ b │ c │
            ├───┼───┼───┤
            │ 1 │ 2 │ 3 │
            ├───┼───┼───┤
            │ 4 │ 5 │ 6 │
            └───┴───┴───┘

        If you do specify a list for ``header_row``, an error is raised::

            sage: T.options(header_row=['x', 'y', 'z'])
            Traceback (most recent call last):
            ...
            TypeError: header_row should be either True or False.
        """
        if kwds:
            for option in ['align', 'frame']:
                if option in kwds:
                    self._options[option] = kwds[option]
            for option in ['header_row', 'header_column']:
                if option in kwds:
                    if not kwds[option]:
                        self._options[option] = kwds[option]
                    elif kwds[option] is True:
                        self._options[option] = kwds[option]
                    else:
                        raise TypeError("%s should be either True or False." % option)
        else:
            return self._options

    def transpose(self):
        r"""
        Return a table which is the transpose of this one:
        rows and columns have been interchanged. Several of the
        properties of the original table are preserved: whether a
        frame is present and any alignment setting. On the other hand,
        header rows are converted to header columns, and vice versa.

        EXAMPLES::

            sage: T = table([[1,2,3], [4,5,6]])
            sage: T.transpose()
              1   4
              2   5
              3   6
            sage: T = table([[1,2,3], [4,5,6]], header_row=['x', 'y', 'z'], frame=True)
            sage: T.transpose()
            ┌───╥───┬───┐
            │ x ║ 1 │ 4 │
            ├───╫───┼───┤
            │ y ║ 2 │ 5 │
            ├───╫───┼───┤
            │ z ║ 3 │ 6 │
            └───╨───┴───┘
        """
        return table(list(zip(*self._rows)),
                     header_row=self._options['header_column'],
                     header_column=self._options['header_row'],
                     frame=self._options['frame'],
                     align=self._options['align'])

    @cached_method
    def _widths(self):
        r"""
        The maximum widths for (the string representation of) each
        column. Used by the :meth:`_repr_` method.

        EXAMPLES::

            sage: table([['a', 'bb', 'ccccc'], [10, -12, 0], [1, 2, 3]])._widths()
            (2, 3, 5)
        """
        nc = len(self._rows[0])

        widths = [0] * nc
        for row in self._rows:
            w = []
            for (idx, x) in zip(range(nc), row):
                w.append(max(widths[idx], len(str(x))))
            widths = w
        return tuple(widths)

    def _repr_(self) -> str:
        r"""
        String representation of a table.

        The class docstring has many examples; here is one more.

        EXAMPLES::

            sage: table([['a', 'bb', 'ccccc'], [10, -12, 0], [1, 2, 3]], align='right') # indirect doctest
               a    bb   ccccc
              10   -12       0
               1     2       3
        """
        rows = self._rows
        nc = len(rows[0])
        if len(rows) == 0 or nc == 0:
            return ""

        top_line = "┌" + "┬".join("─" * (x + 2) for x in self._widths()) + "┐\n"
        frame_line = "├" + "┼".join("─" * (x + 2) for x in self._widths()) + "┤\n"
        bottom_line = "└" + "┴".join("─" * (x + 2) for x in self._widths()) + "┘\n"

        if self._options['header_column'] and self._options['frame']:
            top_line = top_line.replace('┬', '╥', 1)
            frame_line = frame_line.replace('┼', '╫', 1)
            bottom_line = bottom_line.replace('┴', '╨', 1)

        if self._options['frame']:
            s = top_line
        else:
            s = ""

        if self._options['header_row']:
            s += self._str_table_row(rows[0], header_row=True)
            rows = rows[1:]

        for row in rows[:-1]:
            s += self._str_table_row(row, header_row=False)
        if rows:
            s += self._str_table_row(rows[-1], header_row=False, last_row=True)
        return s.strip("\n")

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: t = table([1, 2, 3])
            sage: t._rich_repr_(dm)    # the doctest backend does not support html
        """
        OutputHtml = display_manager.types.OutputHtml
        if OutputHtml in display_manager.supported_output():
            return OutputHtml(self._html_())

    def _str_table_row(self, row, header_row=False, last_row=False):
        r"""
        String representation of a row of a table. Used by the
        :meth:`_repr_` method.

        EXAMPLES::

            sage: T = table([['a', 'bb', 'ccccc'], [10, -12, 0], [1, 2, 3]], align='right')
            sage: T._str_table_row([1,2,3])
            '   1     2       3\n'
            sage: T._str_table_row([1,2,3], True)
            '   1     2       3\n├────┼─────┼───────┤\n'
            sage: T.options(header_column=True)
            sage: T._str_table_row([1,2,3], True)
            '   1 │   2       3\n├────┼─────┼───────┤\n'
            sage: T.options(frame=True)
            sage: T._str_table_row([1,2,3], False)
            '│  1 ║   2 │     3 │\n├────╫─────┼───────┤\n'

        Check that :issue:`14601` has been fixed::

            sage: table([['111111', '222222', '333333']])._str_table_row([False,True,None], False)
            '  False    True     None\n'
        """
        frame = self._options['frame']
        widths = self._widths()

        top_line = "┌" + "┬".join("─" * (x + 2) for x in self._widths()) + "┐\n"
        frame_line = "├" + "┼".join("─" * (x + 2) for x in self._widths()) + "┤\n"
        bottom_line = "└" + "┴".join("─" * (x + 2) for x in self._widths()) + "┘\n"

        align = self._options['align']
        if align == 'right':
            align_char = '>'
        elif align == 'center':
            align_char = '^'
        else:
            align_char = '<'

        if frame:
            s = "│ "
        else:
            s = "  "

        if self._options['header_column']:
            if frame:
                top_line = top_line.replace('┬', '╥', 1)
                frame_line = frame_line.replace('┼', '╫', 1)
                bottom_line = bottom_line.replace('┴', '╨', 1)

            s += ("{!s:" + align_char + str(widths[0]) + "}").format(row[0])
            if frame:
                s += " ║ "
            else:
                s += " │ "
            row = row[1:]
            widths = widths[1:]

        for entry, width in zip(row, widths):
            s += ("{!s:" + align_char + str(width) + "}").format(entry)
            if frame:
                s += " │ "
            else:
                s += "   "
        s = s.rstrip(' ')
        s += "\n"
        if frame and header_row:
            double_line = frame_line.replace('─', '═').replace('├', '╞')
            double_line = double_line.replace('┤', '╡').replace('┼', '╪')
            s += double_line
        elif frame or header_row:
            s += bottom_line if last_row else frame_line
        return s

    def _latex_(self):
        r"""
        LaTeX representation of a table.

        If an entry is a Sage object, it is replaced by its LaTeX
        representation, delimited by dollar signs (i.e., ``x`` is
        replaced by ``$latex(x)$``). If an entry is a string, the
        dollar signs are not automatically added, so tables can
        include both plain text and mathematics.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.misc.table import table
            sage: a = [[r'$\sin(x)$', '$x$', 'text'],                                   # needs sage.modules
            ....:      [1, 34342, 3], [identity_matrix(2), 5, 6]]
            sage: latex(table(a))  # indirect doctest                                   # needs sage.modules
            \begin{tabular}{lll}
            $\sin(x)$ & $x$ & text \\
            $1$ & $34342$ & $3$ \\
            $\left(\begin{array}{rr}
            1 & 0 \\
            0 & 1
            \end{array}\right)$ & $5$ & $6$ \\
            \end{tabular}
            sage: latex(table(a, frame=True, align='center'))                           # needs sage.modules
            \begin{tabular}{|c|c|c|} \hline
            $\sin(x)$ & $x$ & text \\ \hline
            $1$ & $34342$ & $3$ \\ \hline
            $\left(\begin{array}{rr}
            1 & 0 \\
            0 & 1
            \end{array}\right)$ & $5$ & $6$ \\ \hline
            \end{tabular}
        """
        from .latex import latex, LatexExpr

        rows = self._rows
        nc = len(rows[0])
        if len(rows) == 0 or nc == 0:
            return ""

        align_char = self._options['align'][0]   # 'l', 'c', 'r'
        if self._options['frame']:
            frame_char = '|'
            frame_str = ' \\hline'
        else:
            frame_char = ''
            frame_str = ''
        if self._options['header_column']:
            head_col_char = '|'
        else:
            head_col_char = ''
        if self._options['header_row']:
            head_row_str = ' \\hline'
        else:
            head_row_str = ''

        # table header
        s = "\\begin{tabular}{"
        s += frame_char + align_char + frame_char + head_col_char
        s += frame_char.join([align_char] * (nc - 1))
        s += frame_char + "}" + frame_str + "\n"
        # first row
        s += " & ".join(LatexExpr(x) if isinstance(x, (str, LatexExpr))
                        else '$' + latex(x).strip() + '$' for x in rows[0])
        s += " \\\\" + frame_str + head_row_str + "\n"
        # other rows
        for row in rows[1:]:
            s += " & ".join(LatexExpr(x) if isinstance(x, (str, LatexExpr))
                            else '$' + latex(x).strip() + '$' for x in row)
            s += " \\\\" + frame_str + "\n"
        s += "\\end{tabular}"
        return s

    def _html_(self):
        r"""
        HTML representation of a table.

        Strings of html will be parsed for math inside dollar and
        double-dollar signs.  2D graphics will be displayed in the
        cells.  Expressions will be latexed.

        The ``align`` option for tables is ignored in HTML
        output. Specifying ``header_column=True`` may not have any
        visible effect in the Sage notebook, depending on the version
        of the notebook.

        OUTPUT: :class:`~sage.misc.html.HtmlFragment`

        EXAMPLES::

            sage: T = table([[r'$\sin(x)$', '$x$', 'text'],                             # needs sage.modules
            ....:            [1, 34342, 3], [identity_matrix(2), 5, 6]])
            sage: T._html_()                                                            # needs sage.modules
            '<div.../div>'
            sage: print(T._html_())                                                     # needs sage.modules
            <div class="notruncate">
            <table  class="table_form">
            <tbody>
            <tr class ="row-a">
            <td style="text-align:left">\(\sin(x)\)</td>
            <td style="text-align:left">\(x\)</td>
            <td style="text-align:left">text</td>
            </tr>
            <tr class ="row-b">
            <td style="text-align:left">\(1\)</td>
            <td style="text-align:left">\(34342\)</td>
            <td style="text-align:left">\(3\)</td>
            </tr>
            <tr class ="row-a">
            <td style="text-align:left">\(\left(\begin{array}{rr}
            1 & 0 \\
            0 & 1
            \end{array}\right)\)</td>
            <td style="text-align:left">\(5\)</td>
            <td style="text-align:left">\(6\)</td>
            </tr>
            </tbody>
            </table>
            </div>

        Note that calling ``html(table(...))`` has the same effect as
        calling ``table(...)._html_()``::

            sage: T = table([["$x$", r"$\sin(x)$"]]                                     # needs sage.symbolic
            ....:            + [(x,n(sin(x), digits=2)) for x in [0..3]],
            ....:           header_row=True, frame=True)
            sage: T                                                                     # needs sage.symbolic
            ┌─────┬───────────┐
            │ $x$ │ $\sin(x)$ │
            ╞═════╪═══════════╡
            │ 0   │ 0.00      │
            ├─────┼───────────┤
            │ 1   │ 0.84      │
            ├─────┼───────────┤
            │ 2   │ 0.91      │
            ├─────┼───────────┤
            │ 3   │ 0.14      │
            └─────┴───────────┘
            sage: print(html(T))                                                        # needs sage.symbolic
            <div class="notruncate">
            <table border="1" class="table_form">
            <tbody>
            <tr>
            <th style="text-align:left">\(x\)</th>
            <th style="text-align:left">\(\sin(x)\)</th>
            </tr>
            <tr class ="row-a">
            <td style="text-align:left">\(0\)</td>
            <td style="text-align:left">\(0.00\)</td>
            </tr>
            <tr class ="row-b">
            <td style="text-align:left">\(1\)</td>
            <td style="text-align:left">\(0.84\)</td>
            </tr>
            <tr class ="row-a">
            <td style="text-align:left">\(2\)</td>
            <td style="text-align:left">\(0.91\)</td>
            </tr>
            <tr class ="row-b">
            <td style="text-align:left">\(3\)</td>
            <td style="text-align:left">\(0.14\)</td>
            </tr>
            </tbody>
            </table>
            </div>
        """
        from itertools import cycle
        rows = self._rows
        header_row = self._options['header_row']
        if self._options['frame']:
            frame = 'border="1"'
        else:
            frame = ''
        s = StringIO()
        if rows:
            s.writelines([
                # If the table has < 100 rows, don't truncate the output in the notebook
                '<div class="notruncate">\n' if len(rows) <= 100 else '<div class="truncate">',
                '<table {} class="table_form">\n'.format(frame),
                '<tbody>\n',
            ])
            # First row:
            if header_row:
                s.write('<tr>\n')
                self._html_table_row(s, rows[0], header=header_row)
                s.write('</tr>\n')
                rows = rows[1:]

            # Other rows:
            for row_class, row in zip(cycle(["row-a", "row-b"]), rows):
                s.write('<tr class ="{}">\n'.format(row_class))
                self._html_table_row(s, row, header=False)
                s.write('</tr>\n')
            s.write('</tbody>\n</table>\n</div>')
        return s.getvalue()

    def _html_table_row(self, file, row, header=False):
        r"""
        Write table row.

        Helper method used by the :meth:`_html_` method.

        INPUT:

        - ``file`` -- file-like object; the table row data will be
          written to it

        - ``row`` -- list with the same number of entries as each row
          of the table

        - ``header`` -- boolean (default: ``False``); if ``True``, treat this
          as a header row, using ``<th>`` instead of ``<td>``

        OUTPUT: this method returns nothing; all output is written to ``file``

        Strings are written verbatim unless they seem to be LaTeX
        code, in which case they are enclosed in a ``script`` tag
        appropriate for MathJax. Sage objects are printed using their
        LaTeX representations.

        EXAMPLES::

            sage: T = table([['a', 'bb', 'ccccc'], [10, -12, 0], [1, 2, 3]])
            sage: from io import StringIO
            sage: s = StringIO()
            sage: T._html_table_row(s, ['a', 2, '$x$'])
            sage: print(s.getvalue())
            <td style="text-align:left">a</td>
            <td style="text-align:left">\(2\)</td>
            <td style="text-align:left">\(x\)</td>
        """
        from .latex import latex
        from .html import math_parse
        import types

        try:
            from sage.plot.all import Graphics
        except ImportError:
            Graphics = ()

        if isinstance(row, types.GeneratorType):
            row = list(row)
        elif not isinstance(row, (list, tuple)):
            row = [row]

        align_char = self._options['align'][0]   # 'l', 'c', 'r'

        if align_char == 'l':
            style = 'text-align:left'
        elif align_char == 'c':
            style = 'text-align:center'
        elif align_char == 'r':
            style = 'text-align:right'
        else:
            style = ''

        style_attr = f' style="{style}"' if style else ''

        column_tag = f'<th{style_attr}>%s</th>\n' if header else f'<td{style_attr}>%s</td>\n'

        if self._options['header_column']:
            first_column_tag = '<th class="ch"{style_attr}>%s</th>\n' if header else '<td class="ch"{style_attr}>%s</td>\n'
        else:
            first_column_tag = column_tag

        # first entry of row
        entry = row[0]
        if isinstance(entry, Graphics):
            file.write(first_column_tag % entry.show(linkmode=True))
        elif isinstance(entry, str):
            file.write(first_column_tag % math_parse(entry))
        else:
            file.write(first_column_tag % (r'\(%s\)' % latex(entry)))

        # other entries
        for column in range(1, len(row)):
            if isinstance(row[column], Graphics):
                file.write(column_tag % row[column].show(linkmode=True))
            elif isinstance(row[column], str):
                file.write(column_tag % math_parse(row[column]))
            else:
                file.write(column_tag % (r'\(%s\)' % latex(row[column])))
