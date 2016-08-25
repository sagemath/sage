"""
Output functions

These are the output functions for latexing partitions and tableaux.

AUTHORS:

- Mike Hansen (?): initial version
- Andrew Mathas (2013-02-14): Added support for displaying conventions and
  lines, and tableaux of skew partition, composition, and
  skew/composition/partition/tableaux tuple shape.
"""
from __future__ import print_function

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
    following :meth:`Tableaux.options``.

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
    principe, be anything but probably should be strings or integers of similar
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
    """
    # first identify where the None's appear in ``array`` and define a
    # function end_line which puts in the required \cline's.
    if with_lines:
        # last position of None in each row
        nones=[1 if not None in row else 1+len(row)-row[::-1].index(None) for row in array]
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
    for r in xrange(len(array)):
        tex+='&'.join('' if c is None else r'\lr{%s}'%c for c in array[r])
        tex+=end_line(r+1)+'\n'
    return tex+r'\end{array}$}'

