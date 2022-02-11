# -*- coding: utf-8 -*-
r"""
Output functions

These are the output functions for latexing and ascii/unicode art versions of
partitions and tableaux.

AUTHORS:

- Mike Hansen (?): initial version
- Andrew Mathas (2013-02-14): Added support for displaying conventions and
  lines, and tableaux of skew partition, composition, and
  skew/composition/partition/tableaux tuple shape.
- Travis Scrimshaw (2020-08): Added support for ascii/unicode art
"""


from string import Template
from sage.combinat.tableau import Tableaux

# The tex macro used to latex individual cells in an array (as a template).
# When using bar should be replaced by '|' or ''.
lr_macro = Template(r'\def\lr#1{\multicolumn{1}{$bar@{\hspace{.6ex}}c@{\hspace{.6ex}}$bar}{\raisebox{-.3ex}{$$#1$$}}}')

def tex_from_array(array, with_lines=True):
    r"""
    Return a latex string for a two dimensional array of partition, composition or skew composition shape

    INPUT:

    - ``array`` -- a list of list
    - ``with_lines`` -- a boolean (default: ``True``)
       Whether to draw a line to separate the entries in the array.

    Empty rows are allowed; however, such rows should be given as
    ``[None]`` rather than ``[]``.

    The array is drawn using either the English or French convention
    following :meth:`Tableaux.options`.

    .. SEEALSO:: :meth:`tex_from_array_tuple`

    EXAMPLES::

        sage: from sage.combinat.output import tex_from_array
        sage: print(tex_from_array([[1,2,3],[4,5]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
        \lr{1}&\lr{2}&\lr{3}\\\cline{1-3}
        \lr{4}&\lr{5}\\\cline{1-2}
        \end{array}$}
        }
        sage: print(tex_from_array([[1,2,3],[4,5]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\\
        \lr{1}&\lr{2}&\lr{3}\\
        \lr{4}&\lr{5}\\
        \end{array}$}
        }
        sage: print(tex_from_array([[1,2,3],[4,5,6,7],[8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-3}
        \lr{1}&\lr{2}&\lr{3}\\\cline{1-4}
        \lr{4}&\lr{5}&\lr{6}&\lr{7}\\\cline{1-4}
        \lr{8}\\\cline{1-1}
        \end{array}$}
        }
        sage: print(tex_from_array([[1,2,3],[4,5,6,7],[8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\\
        \lr{1}&\lr{2}&\lr{3}\\
        \lr{4}&\lr{5}&\lr{6}&\lr{7}\\
        \lr{8}\\
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{3-3}
        &&\lr{3}\\\cline{2-4}
        &\lr{5}&\lr{6}&\lr{7}\\\cline{1-4}
        \lr{8}\\\cline{1-1}
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[None,8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{3-3}
        &&\lr{3}\\\cline{2-4}
        &\lr{5}&\lr{6}&\lr{7}\\\cline{2-4}
        &\lr{8}\\\cline{2-2}
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\\
        &&\lr{3}\\
        &\lr{5}&\lr{6}&\lr{7}\\
        \lr{8}\\
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[None,8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\\
        &&\lr{3}\\
        &\lr{5}&\lr{6}&\lr{7}\\
        &\lr{8}\\
        \end{array}$}
        }
        sage: Tableaux.options.convention="french"
        sage: print(tex_from_array([[1,2,3],[4,5]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\cline{1-2}
        \lr{4}&\lr{5}\\\cline{1-3}
        \lr{1}&\lr{2}&\lr{3}\\\cline{1-3}
        \end{array}$}
        }
        sage: print(tex_from_array([[1,2,3],[4,5]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\\
        \lr{4}&\lr{5}\\
        \lr{1}&\lr{2}&\lr{3}\\
        \end{array}$}
        }
        sage: print(tex_from_array([[1,2,3],[4,5,6,7],[8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{4}c}\cline{1-1}
        \lr{8}\\\cline{1-4}
        \lr{4}&\lr{5}&\lr{6}&\lr{7}\\\cline{1-4}
        \lr{1}&\lr{2}&\lr{3}\\\cline{1-3}
        \end{array}$}
        }
        sage: print(tex_from_array([[1,2,3],[4,5,6,7],[8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{4}c}\\
        \lr{8}\\
        \lr{4}&\lr{5}&\lr{6}&\lr{7}\\
        \lr{1}&\lr{2}&\lr{3}\\
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{4}c}\cline{1-1}
        \lr{8}\\\cline{1-4}
        &\lr{5}&\lr{6}&\lr{7}\\\cline{2-4}
        &&\lr{3}\\\cline{3-3}
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[None,8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{4}c}\cline{2-2}
        &\lr{8}\\\cline{2-4}
        &\lr{5}&\lr{6}&\lr{7}\\\cline{2-4}
        &&\lr{3}\\\cline{3-3}
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{4}c}\\
        \lr{8}\\
        &\lr{5}&\lr{6}&\lr{7}\\
        &&\lr{3}\\
        \end{array}$}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[None,8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{4}c}\\
        &\lr{8}\\
        &\lr{5}&\lr{6}&\lr{7}\\
        &&\lr{3}\\
        \end{array}$}
        }
        sage: Tableaux.options._reset()
    """
    lr=lr_macro.substitute(bar='|' if with_lines else '')
    if Tableaux.options.convention == "English":
        return '{%s\n%s\n}' % (lr, tex_from_skew_array(array, with_lines))
    else:
        return '{%s\n%s\n}' % (lr, tex_from_skew_array(array[::-1], with_lines, align='t'))


def tex_from_array_tuple(a_tuple, with_lines=True):
    r"""
    Return a latex string for a tuple of two dimensional array of partition,
    composition or skew composition shape.

    INPUT:

    - ``a_tuple`` -- a tuple of lists of lists
    - ``with_lines`` -- a boolean (default: ``True``)
      Whether to draw lines to separate the entries in the components of ``a_tuple``.

    .. SEEALSO:: :meth:`tex_from_array` for the description of each array

    EXAMPLES::

        sage: from sage.combinat.output import tex_from_array_tuple
        sage: print(tex_from_array_tuple([[[1,2,3],[4,5]],[],[[None,6,7],[None,8],[9]]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
        \lr{1}&\lr{2}&\lr{3}\\\cline{1-3}
        \lr{4}&\lr{5}\\\cline{1-2}
        \end{array}$},\emptyset,\raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{2-3}
        &\lr{6}&\lr{7}\\\cline{2-3}
        &\lr{8}\\\cline{1-2}
        \lr{9}\\\cline{1-1}
        \end{array}$}
        }
        sage: print(tex_from_array_tuple([[[1,2,3],[4,5]],[],[[None,6,7],[None,8],[9]]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\\
        \lr{1}&\lr{2}&\lr{3}\\
        \lr{4}&\lr{5}\\
        \end{array}$},\emptyset,\raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\\
        &\lr{6}&\lr{7}\\
        &\lr{8}\\
        \lr{9}\\
        \end{array}$}
        }
        sage: Tableaux.options.convention="french"
        sage: print(tex_from_array_tuple([[[1,2,3],[4,5]],[],[[None,6,7],[None,8],[9]]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\cline{1-2}
        \lr{4}&\lr{5}\\\cline{1-3}
        \lr{1}&\lr{2}&\lr{3}\\\cline{1-3}
        \end{array}$},\emptyset,\raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\cline{1-1}
        \lr{9}\\\cline{1-2}
        &\lr{8}\\\cline{2-3}
        &\lr{6}&\lr{7}\\\cline{2-3}
        \end{array}$}
        }
        sage: print(tex_from_array_tuple([[[1,2,3],[4,5]],[],[[None,6,7],[None,8],[9]]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\\
        \lr{4}&\lr{5}\\
        \lr{1}&\lr{2}&\lr{3}\\
        \end{array}$},\emptyset,\raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\\
        \lr{9}\\
        &\lr{8}\\
        &\lr{6}&\lr{7}\\
        \end{array}$}
        }
        sage: Tableaux.options._reset()
    """
    lr=lr_macro.substitute(bar='|' if with_lines else '')
    if Tableaux.options.convention == "English":
        return '{%s\n%s\n}' % (lr, ','.join(
            r'\emptyset' if comp==[] else tex_from_skew_array(comp, with_lines) for comp in a_tuple))
    else:
        return '{%s\n%s\n}' % (lr, ','.join(
            r'\emptyset' if comp==[] else tex_from_skew_array(comp[::-1], with_lines, align='t') for comp in a_tuple))


def tex_from_skew_array(array, with_lines=False, align='b'):
    r"""
    This function creates latex code for a "skew composition" ``array``.
    That is, for a two dimensional array in which each row can begin with
    an arbitrary number ``None``'s and the remaining entries could, in
    principle, be anything but probably should be strings or integers of similar
    width. A row consisting completely of ``None``'s is allowed.

    INPUT:

    - ``array`` -- The array

    - ``with_lines`` -- (Default: ``False``) If ``True`` lines are drawn, if
      ``False`` they are not

    - ``align`` -- (Default: ``'b'``) Determines the alignment on the latex
      array environments

    EXAMPLES::

        sage: array=[[None, 2,3,4],[None,None],[5,6,7,8]]
        sage: print(sage.combinat.output.tex_from_skew_array(array))
        \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\\
        &\lr{2}&\lr{3}&\lr{4}\\
        &\\
        \lr{5}&\lr{6}&\lr{7}&\lr{8}\\
        \end{array}$}

    TESTS::

        sage: sage.combinat.output.tex_from_skew_array([(1,2,3), (2,3,4)])
        '\\raisebox{-.6ex}{$\\begin{array}[b]{*{3}c}\\\\\n\\lr{1}&\\lr{2}&\\lr{3}\\\\\n\\lr{2}&\\lr{3}&\\lr{4}\\\\\n\\end{array}$}'
        sage: sage.combinat.output.tex_from_skew_array([((1,2,),)])
        '\\raisebox{-.6ex}{$\\begin{array}[b]{*{1}c}\\\\\n\\lr{(1, 2)}\\\\\n\\end{array}$}'
    """
    # first identify where the None's appear in ``array`` and define a
    # function end_line which puts in the required \cline's.
    if with_lines:
        # last position of None in each row
        nones = [1 if None not in row else 1 + len(row) - row[::-1].index(None)
                 for row in array]

        def end_line(r):
            # in a slightly unpythonic way, we label the lines as 0, 1, ..., len(array)
            if r==0:
                return r'\cline{%s-%s}'%(nones[0],len(array[0]))
            elif r==len(array):
                start=nones[r-1]
                finish=len(array[r-1])
            else:
                start=min(nones[r], nones[r-1])
                finish=max(len(array[r]), len(array[r-1]))
            return r'\\' if start>finish else r'\\\cline{%s-%s}'%(start, finish)
    else:
        end_line=lambda r: r'\\'

    # now we draw the array
    tex=r'\raisebox{-.6ex}{$\begin{array}[%s]{*{%s}c}'%(align,max(map(len,array)))
    tex+=end_line(0)+'\n'
    for r in range(len(array)):
        tex+='&'.join('' if c is None else r'\lr{%s}'%(c,) for c in array[r])
        tex+=end_line(r+1)+'\n'
    return tex+r'\end{array}$}'


def ascii_art_table(data, use_unicode=False, convention="English"):
    r"""
    Return an ascii art table of ``data``.

    EXAMPLES::

        sage: from sage.combinat.output import ascii_art_table

        sage: data = [[None, None, 1], [2, 2], [3,4,5], [None, None, 10], [], [6]]
        sage: print(ascii_art_table(data))
                +----+
                | 1  |
        +---+---+----+
        | 2 | 2 |
        +---+---+----+
        | 3 | 4 | 5  |
        +---+---+----+
                | 10 |
                +----+
        <BLANKLINE>
        +---+
        | 6 |
        +---+
        sage: print(ascii_art_table(data, use_unicode=True))
                ┌────┐
                │ 1  │
        ┌───┬───┼────┘
        │ 2 │ 2 │
        ├───┼───┼────┐
        │ 3 │ 4 │ 5  │
        └───┴───┼────┤
                │ 10 │
                └────┘
        <BLANKLINE>
        ┌───┐
        │ 6 │
        └───┘

        sage: data = [[1, None, 2], [None, 2]]
        sage: print(ascii_art_table(data))
        +---+   +---+
        | 1 |   | 2 |
        +---+---+---+
            | 2 |
            +---+
        sage: print(ascii_art_table(data, use_unicode=True))
        ┌───┐   ┌───┐
        │ 1 │   │ 2 │
        └───┼───┼───┘
            │ 2 │
            └───┘
    """
    if use_unicode:
        import unicodedata
        v  = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL')
        h  = unicodedata.lookup('BOX DRAWINGS LIGHT HORIZONTAL')
        dl = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND LEFT')
        dr = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND RIGHT')
        ul = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND LEFT')
        ur = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND RIGHT')
        vr = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND RIGHT')
        vl = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND LEFT')
        uh = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND HORIZONTAL')
        dh = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND HORIZONTAL')
        vh = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND HORIZONTAL')
        from sage.typeset.unicode_art import unicode_art as art
    else:
        v = '|'
        h = '-'
        dl = dr = ul = ur = vr = vl = uh = dh = vh = '+'
        from sage.typeset.ascii_art import ascii_art as art

    if not data:
        return dr + dl + '\n' + ur + ul

    # Convert the input into a rectangular array with the top and bottom row
    #   being all None's for ease later on.
    ncols = max(len(row) for row in data)
    str_tab = [[None]*ncols] + [[art(val) if val is not None else None for val in row] + [None]*(ncols-len(row))
                                for row in data]
    str_tab.append([None]*ncols)
    # Get the widths of the columns
    col_widths = [1]*len(str_tab[0])
    if use_unicode:
        # Special handling of overline not adding to printed length
        def get_len(e):
            if e is None:
                return 0
            return len(e) - list(str(e)).count(u"\u0304")
    else:
        def get_len(e):
            if e is None:
                return 0
            return len(e)
    for row in str_tab:
        for i,e in enumerate(row):
            col_widths[i] = max(col_widths[i], get_len(e))

    matr = []  # just the list of lines
    for nrow,row in enumerate(str_tab):
        if nrow == 0: # skip the first row
            continue

        l1 = ""
        l2 = ""
        for i,(e,w) in enumerate(zip(row, col_widths)):
            prev_row = str_tab[nrow-1]
            if i == 0:
                if e is None:
                    if prev_row[i] is None:
                        l1 += " "*(3+w)
                    else:
                        l1 += ur + h*(2+w)
                    l2 += " "*(3+w)
                else:
                    if prev_row[i] is None:
                        l1 += dr + h*(2+w)
                    else:
                        l1 += vr + h*(2+w)
                    l2 += "{} {:^{width}} ".format(v, e, width=w)
            else:
                if e is None:
                    if row[i-1] is None:
                        if prev_row[i-1] is None:
                            if prev_row[i] is None:
                                l1 += " "*(3+w)
                            else:
                                l1 += ur + h*(2+w)
                        else:
                            if prev_row[i] is None:
                                l1 += ul + " "*(2+w)
                            else:
                                l1 += uh + h*(2+w)
                        l2 += " "*(3+w)
                    else:
                        if prev_row[i-1] is None:
                            if prev_row[i] is None:
                                l1 += dl + " "*(2+w)
                            else:
                                l1 += vh + h*(2+w)
                        else:
                            if prev_row[i] is None:
                                l1 += vl + " "*(2+w)
                            else:
                                l1 += vh + h*(2+w)
                        l2 += v + " "*(2+w)
                else:
                    if row[i-1] is None:
                        if prev_row[i-1] is None:
                            if prev_row[i] is None:
                                l1 += dr + h*(2+w)
                            else:
                                l1 += vr + h*(2+w)
                        else:
                            l1 += vh + h*(2+w)
                    else:
                        if prev_row[i-1] is None and prev_row[i] is None:
                                l1 += dh + h*(2+w)
                        else:
                            l1 += vh + h*(2+w)
                    l2 += "{} {:^{width}} ".format(v, e, width=w)

        if row[-1] is None:
            if prev_row[-1] is None:
                l1 += " "
            else:
                l1 += ul
            l2 += " "
        else:
            if prev_row[-1] is None:
                l1 += dl
            else:
                l1 += vl
            l2 += v

        matr.append(l1)
        matr.append(l2)

    matr.pop()  # Remove the last row (which is blank)

    if convention == "English":
        return "\n".join(matr)
    else:
        output = "\n".join(reversed(matr))
        if use_unicode:
            tr = {
                ord(dl): ul, ord(dr): ur,
                ord(ul): dl, ord(ur): dr,
                ord(dh): uh, ord(uh): dh}
            return output.translate(tr)
        else:
            return output

