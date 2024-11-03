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
    Return a latex string for a two dimensional array of partition, composition
    or skew composition shape.

    INPUT:

    - ``array`` -- list of list
    - ``with_lines`` -- boolean (default: ``True``); whether to draw a line to
      separate the entries in the array

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
        sage: Tableaux.options.convention="russian"
        sage: print(tex_from_array([[1,2,3],[4,5]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{3}c}\cline{1-2}
        \lr{\rotatebox{-45}{4}}&\lr{\rotatebox{-45}{5}}\\\cline{1-3}
        \lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{2}}&\lr{\rotatebox{-45}{3}}\\\cline{1-3}
        \end{array}$}}
        }
        sage: print(tex_from_array([[1,2,3],[4,5]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{3}c}\\
        \lr{\rotatebox{-45}{4}}&\lr{\rotatebox{-45}{5}}\\
        \lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{2}}&\lr{\rotatebox{-45}{3}}\\
        \end{array}$}}
        }
        sage: print(tex_from_array([[1,2,3],[4,5,6,7],[8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{4}c}\cline{1-1}
        \lr{\rotatebox{-45}{8}}\\\cline{1-4}
        \lr{\rotatebox{-45}{4}}&\lr{\rotatebox{-45}{5}}&\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\\cline{1-4}
        \lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{2}}&\lr{\rotatebox{-45}{3}}\\\cline{1-3}
        \end{array}$}}
        }
        sage: print(tex_from_array([[1,2,3],[4,5,6,7],[8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{4}c}\\
        \lr{\rotatebox{-45}{8}}\\
        \lr{\rotatebox{-45}{4}}&\lr{\rotatebox{-45}{5}}&\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\
        \lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{2}}&\lr{\rotatebox{-45}{3}}\\
        \end{array}$}}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{4}c}\cline{1-1}
        \lr{\rotatebox{-45}{8}}\\\cline{1-4}
        &\lr{\rotatebox{-45}{5}}&\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\\cline{2-4}
        &&\lr{\rotatebox{-45}{3}}\\\cline{3-3}
        \end{array}$}}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[None,8]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{4}c}\cline{2-2}
        &\lr{\rotatebox{-45}{8}}\\\cline{2-4}
        &\lr{\rotatebox{-45}{5}}&\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\\cline{2-4}
        &&\lr{\rotatebox{-45}{3}}\\\cline{3-3}
        \end{array}$}}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{4}c}\\
        \lr{\rotatebox{-45}{8}}\\
        &\lr{\rotatebox{-45}{5}}&\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\
        &&\lr{\rotatebox{-45}{3}}\\
        \end{array}$}}
        }
        sage: print(tex_from_array([[None,None,3],[None,5,6,7],[None,8]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{4}c}\\
        &\lr{\rotatebox{-45}{8}}\\
        &\lr{\rotatebox{-45}{5}}&\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\
        &&\lr{\rotatebox{-45}{3}}\\
        \end{array}$}}
        }

        sage: Tableaux.options._reset()
    """
    lr = lr_macro.substitute(bar='|' if with_lines else '')
    if Tableaux.options.convention == "English":
        return '{%s\n%s\n}' % (lr, tex_from_skew_array(array, with_lines))
    else:
        return '{%s\n%s\n}' % (lr, tex_from_skew_array(array[::-1], with_lines, align='t'))


def tex_from_array_tuple(a_tuple, with_lines=True):
    r"""
    Return a latex string for a tuple of two dimensional array of partition,
    composition or skew composition shape.

    INPUT:

    - ``a_tuple`` -- tuple of lists of lists
    - ``with_lines`` -- boolean (default: ``True``); whether to draw lines to
      separate the entries in the components of ``a_tuple``

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
        sage: Tableaux.options.convention="russian"
        sage: print(tex_from_array_tuple([[[1,2,3],[4,5]],[],[[None,6,7],[None,8],[9]]]))
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{3}c}\cline{1-2}
        \lr{\rotatebox{-45}{4}}&\lr{\rotatebox{-45}{5}}\\\cline{1-3}
        \lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{2}}&\lr{\rotatebox{-45}{3}}\\\cline{1-3}
        \end{array}$}},\emptyset,\raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{3}c}\cline{1-1}
        \lr{\rotatebox{-45}{9}}\\\cline{1-2}
        &\lr{\rotatebox{-45}{8}}\\\cline{2-3}
        &\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\\cline{2-3}
        \end{array}$}}
        }
        sage: print(tex_from_array_tuple([[[1,2,3],[4,5]],[],[[None,6,7],[None,8],[9]]], with_lines=False))
        {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{3}c}\\
        \lr{\rotatebox{-45}{4}}&\lr{\rotatebox{-45}{5}}\\
        \lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{2}}&\lr{\rotatebox{-45}{3}}\\
        \end{array}$}},\emptyset,\raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{3}c}\\
        \lr{\rotatebox{-45}{9}}\\
        &\lr{\rotatebox{-45}{8}}\\
        &\lr{\rotatebox{-45}{6}}&\lr{\rotatebox{-45}{7}}\\
        \end{array}$}}
        }

        sage: Tableaux.options._reset()
    """
    lr = lr_macro.substitute(bar='|' if with_lines else '')
    if Tableaux.options.convention == "English":
        return '{%s\n%s\n}' % (lr, ','.join(
            r'\emptyset' if comp == [] else tex_from_skew_array(comp, with_lines) for comp in a_tuple))
    else:
        return '{%s\n%s\n}' % (lr, ','.join(
            r'\emptyset' if comp == [] else tex_from_skew_array(comp[::-1], with_lines, align='t') for comp in a_tuple))


def tex_from_skew_array(array, with_lines=False, align='b'):
    r"""
    This function creates latex code for a "skew composition" ``array``.
    That is, for a two dimensional array in which each row can begin with
    an arbitrary number ``None``'s and the remaining entries could, in
    principle, be anything but probably should be strings or integers of similar
    width. A row consisting completely of ``None``'s is allowed.

    INPUT:

    - ``array`` -- the array

    - ``with_lines`` -- (default: ``False``) if ``True`` lines are drawn, if
      ``False`` they are not

    - ``align`` -- (default: ``'b'``) determine the alignment on the latex
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
            if r == 0:
                return r'\cline{%s-%s}' % (nones[0], len(array[0]))
            elif r == len(array):
                start = nones[r-1]
                finish = len(array[r-1])
            else:
                start = min(nones[r], nones[r-1])
                finish = max(len(array[r]), len(array[r-1]))
            return r'\\' if start > finish else r'\\\cline{%s-%s}' % (start, finish)
    else:
        end_line = lambda r: r'\\'

    # now we draw the array
    raisebox_start = r'\raisebox{-.6ex}{'
    raisebox_end = r'}'
    lr_start = r'\lr{'
    lr_end = r'}'
    if Tableaux.options.convention == "Russian":
        raisebox_start += r'\rotatebox{45}{'
        raisebox_end += r'}'
        lr_start += r'\rotatebox{-45}{'
        lr_end += r'}'

    tex = r'%s$\begin{array}[%s]{*{%s}c}' % (raisebox_start, align, max(map(len, array)))
    tex += end_line(0)+'\n'
    for r in range(len(array)):
        tex += '&'.join('' if c is None else r'%s%s%s' % (lr_start, c, lr_end) for c in array[r])
        tex += end_line(r+1)+'\n'
    return tex+r'\end{array}$'+raisebox_end


def ascii_art_table(data, use_unicode=False, convention='English'):
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
    if convention == "Russian":
        return ascii_art_table_russian(data, use_unicode)

    if use_unicode:
        import unicodedata
        v = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL')
        h = unicodedata.lookup('BOX DRAWINGS LIGHT HORIZONTAL')
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
        # return dr + dl + '\n' + ur + ul
        return ''

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
            return len(e) - list(str(e)).count("\u0304")
    else:
        def get_len(e):
            if e is None:
                return 0
            return len(e)
    for row in str_tab:
        for i, e in enumerate(row):
            col_widths[i] = max(col_widths[i], get_len(e))

    matr = []  # just the list of lines
    for nrow, row in enumerate(str_tab):
        if nrow == 0:  # skip the first row
            continue

        l1 = ""
        l2 = ""
        for i, (e, w) in enumerate(zip(row, col_widths)):
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


def ascii_art_table_russian(data, use_unicode=False, compact=False):
    r"""
    Return an ascii art table of ``data`` for the russian convention.

    EXAMPLES::

        sage: from sage.combinat.output import ascii_art_table_russian
        sage: data = [[None, None, 1], [2, 2], [3,4,5], [None, None, 10], [], [6]]
        sage: print(ascii_art_table_russian(data))
           / \         / \
          /   \       /   \
         \  6  /     \ 10  \
          \   /       \   / \
           \ /         \ /   \
                        X  5  /
                       / \   /
                      /   \ /
                     /  4  X
                    / \   / \   / \
                   /   \ /   \ /   \
                  \  3  X  2  X  1  /
                   \   / \   / \   /
                    \ /   \ /   \ /
                     \  2  /
                      \   /
                       \ /
        sage: print(ascii_art_table_russian(data, use_unicode=True))
           ╱ ╲         ╱ ╲
          ╱   ╲       ╱   ╲
         ╲  6  ╱     ╲ 10  ╲
          ╲   ╱       ╲   ╱ ╲
           ╲ ╱         ╲ ╱   ╲
                        ╳  5  ╱
                       ╱ ╲   ╱
                      ╱   ╲ ╱
                     ╱  4  ╳
                    ╱ ╲   ╱ ╲   ╱ ╲
                   ╱   ╲ ╱   ╲ ╱   ╲
                  ╲  3  ╳  2  ╳  1  ╱
                   ╲   ╱ ╲   ╱ ╲   ╱
                    ╲ ╱   ╲ ╱   ╲ ╱
                     ╲  2  ╱
                      ╲   ╱
                       ╲ ╱
        sage: data = [[1, None, 2], [None, 2]]
        sage: print(ascii_art_table_russian(data))
          / \ / \
         \ 2 X 2 /
          \ / \ /
           X
          / \
         \ 1 /
          \ /
        sage: print(ascii_art_table_russian(data, use_unicode=True))
          ╱ ╲ ╱ ╲
         ╲ 2 ╳ 2 ╱
          ╲ ╱ ╲ ╱
           ╳
          ╱ ╲
         ╲ 1 ╱
          ╲ ╱
    """
    if use_unicode:
        import unicodedata
        urdl = unicodedata.lookup('BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT')
        uldr = unicodedata.lookup('BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT')
        x = unicodedata.lookup('BOX DRAWINGS LIGHT DIAGONAL CROSS')
    else:
        urdl = '/'
        uldr = '\\'
        x = 'X'

    if not data:
        # return urdl + uldr + '\n' + uldr + urdl
        return ''

    if use_unicode:
        # Special handling of overline not adding to printed length
        def get_len(e):
            if e is None:
                return 0
            return len(e) - list(str(e)).count("\u0304")
    else:
        def get_len(e):
            if e is None:
                return 0
            return len(e)

    # Length of max string (ensure it's odd)
    str_tab = [[str(val) if val is not None else None for val in row] for row in data]
    max_str = max([max([1] + [get_len(e) for e in row]) for row in str_tab])
    max_str = max_str + 1 - (max_str % 2)

    if compact:
        diag_length = max_str
    else:
        diag_length = max_str + 2  # space on both sides

    row_height = int((diag_length + 1) // 2)
    max_height = max(a + len(val) for a, val in enumerate(str_tab))
    str_list = []
    for k in range(max_height, -1, -1):
        for i in range(row_height):
            if k == max_height and i == 0:
                continue
            st = ' ' * ((max_height - k) * row_height)
            for j in range(k + 1):
                N_box = box_exists(str_tab, k-j+1, j)
                S_box = box_exists(str_tab, k-j, j-1)
                SE_box = box_exists(str_tab, k-j-1, j)
                E_box = box_exists(str_tab, k-j, j)
                W_box = box_exists(str_tab, k-j+1, j-1)
                if i == 0:
                    if (N_box and S_box) or (W_box and E_box):
                        st += x
                    elif (E_box and S_box) or (W_box and N_box):
                        st += urdl
                    elif (E_box and N_box) or (W_box and S_box):
                        st += uldr
                    elif E_box:
                        st += uldr
                    elif W_box:
                        st += urdl
                    else:
                        st += ' '
                    if E_box:
                        st_num = str_tab[k-j][j]
                        ln_left = len(st_num) // 2
                        st += st_num.rjust(row_height - 1 - ln_left + len(st_num), ' ').ljust(diag_length, ' ')
                    else:
                        st += ' ' * diag_length
                    if j == k and E_box:
                        st += urdl
                else:
                    lstr = ' '
                    rstr = ' '
                    if E_box or S_box:
                        lstr = uldr
                    if E_box or SE_box:
                        rstr = urdl
                    st += ' ' * i
                    st += lstr
                    st += ' ' * (2 * (row_height - i) - 1)
                    st += rstr
                    st += ' ' * (i-1)
            str_list.append(st)

    import re
    mm = min(len(re.search('^ +', l)[0]) for l in str_list) - 1
    str_list = [l[mm:].rstrip() for l in str_list]
    while not str_list[-1]:
        str_list.pop()
    return "\n".join(str_list)


def box_exists(tab, i, j) -> bool:
    r"""
    Return ``True`` if ``tab[i][j]`` exists and is not ``None``; in particular this
    allows for `tab[i][j]` to be ``''`` or ``0``.

    INPUT:

    - ``tab`` -- list of lists
    - ``i`` -- first coordinate
    - ``j`` -- second coordinate

    TESTS::

        sage: from sage.combinat.output import box_exists
        sage: tab = [[1,None,'', 0],[None]]
        sage: box_exists(tab, 0, 0)
        True
        sage: box_exists(tab, 0, 1)
        False
        sage: box_exists(tab, 0, 2)
        True
        sage: box_exists(tab, 0, 3)
        True
        sage: box_exists(tab, 0, 4)
        False
        sage: box_exists(tab, 1, 0)
        False
        sage: box_exists(tab, 1, 1)
        False
        sage: box_exists(tab, 0, -1)
        False
    """
    if j < 0 or i < 0:
        return False
    return len(tab) > i and len(tab[i]) > j and tab[i][j] is not None
