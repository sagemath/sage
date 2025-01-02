r"""
Skew Partitions

A skew partition ``skp`` of size `n` is a pair of
partitions `[p_1, p_2]` where `p_1` is a
partition of the integer `n_1`, `p_2` is a
partition of the integer `n_2`, `p_2` is an inner
partition of `p_1`, and `n = n_1 - n_2`. We say
that `p_1` and `p_2` are respectively the *inner*
and *outer* partitions of ``skp``.

A skew partition can be depicted by a diagram made of rows of
cells, in the same way as a partition. Only the cells of the outer
partition `p_1` which are not in the inner partition
`p_2` appear in the picture. For example, this is the
diagram of the skew partition [[5,4,3,1],[3,3,1]].

::

    sage: print(SkewPartition([[5,4,3,1],[3,3,1]]).diagram())
       **
       *
     **
    *

A skew partition can be *connected*, which can easily be described
in graphic terms: for each pair of consecutive rows, there are at
least two cells (one in each row) which have a common edge. This is
the diagram of the connected skew partition ``[[5,4,3,1],[3,1]]``::

    sage: print(SkewPartition([[5,4,3,1],[3,1]]).diagram())
       **
     ***
    ***
    *
    sage: SkewPartition([[5,4,3,1],[3,1]]).is_connected()
    True

The first example of a skew partition is not a connected one.

Applying a reflection with respect to the main diagonal yields the
diagram of the *conjugate skew partition*, here
``[[4,3,3,2,1],[3,3,2]]``::

    sage: SkewPartition([[5,4,3,1],[3,3,1]]).conjugate()
    [4, 3, 3, 2, 1] / [3, 2, 2]
    sage: print(SkewPartition([[5,4,3,1],[3,3,1]]).conjugate().diagram())
       *
      *
      *
    **
    *

The *outer corners* of a skew partition are the corners of its
outer partition. The *inner corners* are the internal corners of
the outer partition when the inner partition is taken off. Shown
below are the coordinates of the inner and outer corners.

::

    sage: SkewPartition([[5,4,3,1],[3,3,1]]).outer_corners()
    [(0, 4), (1, 3), (2, 2), (3, 0)]
    sage: SkewPartition([[5,4,3,1],[3,3,1]]).inner_corners()
    [(0, 3), (2, 1), (3, 0)]

EXAMPLES:

There are 9 skew partitions of size 3, with no empty row nor empty
column::

    sage: SkewPartitions(3).cardinality()
    9
    sage: SkewPartitions(3).list()
    [[3] / [],
     [2, 1] / [],
     [3, 1] / [1],
     [2, 2] / [1],
     [3, 2] / [2],
     [1, 1, 1] / [],
     [2, 2, 1] / [1, 1],
     [2, 1, 1] / [1],
     [3, 2, 1] / [2, 1]]

There are 4 connected skew partitions of size 3::

    sage: SkewPartitions(3, overlap=1).cardinality()
    4
    sage: SkewPartitions(3, overlap=1).list()
    [[3] / [], [2, 1] / [], [2, 2] / [1], [1, 1, 1] / []]

This is the conjugate of the skew partition ``[[4,3,1], [2]]``

::

    sage: SkewPartition([[4,3,1], [2]]).conjugate()
    [3, 2, 2, 1] / [1, 1]

Geometrically, we just applied a reflection with respect to the main
diagonal on the diagram of the partition. Of course, this operation
is an involution::

    sage: SkewPartition([[4,3,1],[2]]).conjugate().conjugate()
    [4, 3, 1] / [2]

The :meth:`jacobi_trudi()` method computes the Jacobi-Trudi matrix. See
[Mac1995]_ for a definition and discussion.

::

    sage: SkewPartition([[4,3,1],[2]]).jacobi_trudi()                                   # needs sage.modules
    [h[2]  h[]    0]
    [h[5] h[3]  h[]]
    [h[6] h[4] h[1]]

This example shows how to compute the corners of a skew partition.

::

    sage: SkewPartition([[4,3,1],[2]]).inner_corners()
    [(0, 2), (1, 0)]
    sage: SkewPartition([[4,3,1],[2]]).outer_corners()
    [(0, 3), (1, 2), (2, 0)]

AUTHORS:

- Mike Hansen: Initial version
- Travis Scrimshaw (2013-02-11): Factored out ``CombinatorialClass``
- Trevor K. Karn (2022-08-03): Add ``outside_corners``
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.structure.global_options import GlobalOptions
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.set import Set
from sage.misc.lazy_import import lazy_import

from sage.combinat.combinat import CombinatorialElement
from sage.combinat.partition import Partitions, _Partitions
from sage.combinat.tableau import Tableaux
from sage.combinat.composition import Compositions

lazy_import('sage.graphs.digraph', 'DiGraph')
lazy_import('sage.matrix.matrix_space', 'MatrixSpace')


class SkewPartition(CombinatorialElement):
    r"""
    A skew partition.

    A skew partition of shape `\lambda / \mu` is the Young diagram from the
    partition `\lambda` and removing the partition `\mu` from the upper-left
    corner in English convention.
    """
    @staticmethod
    def __classcall_private__(cls, skp):
        """
        Return the skew partition object corresponding to ``skp``.

        EXAMPLES::

            sage: skp = SkewPartition([[3,2,1],[2,1]]); skp
            [3, 2, 1] / [2, 1]
            sage: skp.inner()
            [2, 1]
            sage: skp.outer()
            [3, 2, 1]
        """
        skp = [_Partitions(p) for p in skp]
        if skp not in SkewPartitions():
            raise ValueError("invalid skew partition: %s" % skp)
        return SkewPartitions()(skp)

    def __init__(self, parent, skp):
        """
        TESTS::

            sage: skp = SkewPartition([[3,2,1],[2,1]])
            sage: TestSuite(skp).run()
        """
        CombinatorialElement.__init__(self, parent,
                [_Partitions(skp[0]), _Partitions(skp[1])])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        For more on the display options, see
        :obj:`SkewPartitions.options`.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2,1]])
            [3, 2, 1] / [2, 1]
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_quotient(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: print(SkewPartition([[3,2,1],[2,1]])._repr_quotient())
            [3, 2, 1] / [2, 1]
        """
        return "%s / %s" % (self[0], self[1])

    def _repr_lists(self):
        """
        Return a string representation of ``self`` as a pair of lists.

        EXAMPLES::

            sage: print(SkewPartition([[3,2,1],[2,1]])._repr_lists())
            [[3, 2, 1], [2, 1]]
        """
        return repr([list(r) for r in self])

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        For more on the latex options, see
        :obj:`SkewPartitions.options`.

        EXAMPLES::

            sage: s = SkewPartition([[5,4,3],[3,1,1]])
            sage: latex(s)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}c}\cline{4-5}
            &&&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{2-5}
            &\lr{\phantom{x}}&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{2-4}
            &\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{2-3}
            \end{array}$}
            }
        """
        return self.parent().options._dispatch(self, '_latex_', 'latex')

    def _latex_diagram(self):
        r"""
        Return a `\LaTeX` representation as a young diagram.

        EXAMPLES::

            sage: print(SkewPartition([[5,4,3],[3,1,1]])._latex_diagram())
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}c}\cline{4-5}
            &&&\lr{\ast}&\lr{\ast}\\\cline{2-5}
            &\lr{\ast}&\lr{\ast}&\lr{\ast}\\\cline{2-4}
            &\lr{\ast}&\lr{\ast}\\\cline{2-3}
            \end{array}$}
            }

        TESTS:

        Check that :issue:`34760` is fixed::

            sage: print(SkewPartition([[],[]])._latex_diagram())
            {\emptyset}
        """
        if not any(self._list):
            return "{\\emptyset}"

        char = self.parent().options.latex_diagram_str

        from sage.combinat.output import tex_from_array
        arr = [[char]*row_size for row_size in self[0]]
        for i, skew_size in enumerate(self[1]): # This is always smaller by containment
            for j in range(skew_size):
                arr[i][j] = None
        return tex_from_array(arr)

    def _latex_young_diagram(self):
        r"""
        Return a `\LaTeX` representation of ``self`` as a young diagram.

        EXAMPLES::

            sage: print(SkewPartition([[5,4,3],[3,1,1]])._latex_young_diagram())
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}c}\cline{4-5}
            &&&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{2-5}
            &\lr{\phantom{x}}&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{2-4}
            &\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{2-3}
            \end{array}$}
            }

        TESTS:

        Check that :issue:`34760` is fixed::

            sage: print(SkewPartition([[],[]])._latex_young_diagram())
            {\emptyset}
        """
        if not any(self._list):
            return "{\\emptyset}"

        from sage.combinat.output import tex_from_array
        arr = [["\\phantom{x}"]*row_size for row_size in self[0]]
        for i, skew_size in enumerate(self[1]): # This is always smaller by containment
            for j in range(skew_size):
                arr[i][j] = None
        return tex_from_array(arr)

    def _latex_marked(self):
        r"""
        Return a `\LaTeX` representation as a marked partition.

        EXAMPLES::

            sage: print(SkewPartition([[5,4,3],[3,1,1]])._latex_marked())
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}c}\cline{1-5}
            \lr{X}&\lr{X}&\lr{X}&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-5}
            \lr{X}&\lr{\phantom{x}}&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-4}
            \lr{X}&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-3}
            \end{array}$}
            }

        TESTS:

        Check that :issue:`34760` is fixed::

            sage: print(SkewPartition([[],[]])._latex_marked())
            {\emptyset}
        """
        if not any(self._list):
            return "{\\emptyset}"

        from sage.combinat.output import tex_from_array
        char = self.parent().options.latex_marking_str
        arr = [["\\phantom{x}"]*row_size for row_size in self[0]]
        for i, skew_size in enumerate(self[1]): # This is always smaller by containment
            for j in range(skew_size):
                arr[i][j] = char
        return tex_from_array(arr)

    def __setstate__(self, state):
        r"""
        In order to maintain backwards compatibility and be able to unpickle
        a old pickle from ``SkewPartition_class`` we have to override the
        default ``__setstate__``.

        EXAMPLES::

            sage: loads(b'x\x9c\x85P\xcbN\xc2@\x14\r\x08>\x06\xf1\xfd~\xbb+.\x9a\xa8\xdf\xe0\xc2McJ\xba4\x93i\xb9v&\xb4\x03w\x1e!,Ht!\xfe\xb6Sh1\xb0qw\xce}\x9c{\xee\xf9\xac\'\x9a\xa5\xe0\'\x83<\x16\x92\x19_\xf7aD\x87L\x19a\xc4@\x92\xae\xa3o\x15\xa3I\xc6\xb4&X\xeb|a}\x82k^\xd4\xa4\x9ci\x8e\x8d\xc0\xa1Lh\x83\xcdw\\\xf7\xe6\x92\xda(\x9b\x18\xab\xc0\xef\x8d%\xcbER\xae/3\xdc\xf0\xa2\x87\xc5\x05MY\x96\xd1\x910\x9c&\xcc@:Pc\x1f2\xc8A\x9a\xf9<n\xae\xf8\xfd\xb3\xba\x10!\xb8\x95P\x1a[\x91\x19!)%)\x18f\x8c"HV\x8dY)\xd0\x02U0T\xa0\xdd\r6[\xb7RA\xcf&@\xb0U\x1e\x9b[\x11\xa0}!?\x84\x14\x06(H\x9b\x83r\x8d\x1e\xd5`4y-\x1b/\x8bz\xb7(\xe3vg\xf2\x83\xed\x10w\xa2\xf6\xf2#\xbb\xd3\x10\xf7\xa6\xb8\x1f\x04\x81\t\xf1\xc0Ez\xc8[\xff?7K\x88\xe0Q!{\x1c\xe2\xc9\x04O=\xde\x08\xb8\x0b\xfe\xac\x0c^\t\x99\x16N\x9diP$g}\xa0\x15\xc1\xf3\xa8\xf6\xfc\x1d\xe2\x05w\xe0\xc9\x81\xcb\x02<:p\x05v\x1a\xf3\xc2\xc65w\xa27\x95\xe8\xadWM\xdcU\xe0\xbe\x18\x05\x1b\xfb\xbf\x8e\x7f\xcc\xbb')
            [3, 2, 1] / [1, 1]
            sage: loads(dumps( SkewPartition([[3,2,1], [1,1]]) ))
            [3, 2, 1] / [1, 1]
        """
        if isinstance(state, dict):   # for old pickles from SkewPartition_class
            self._set_parent(SkewPartitions())
            self.__dict__ = state
        else:
            self._set_parent(state[0])
            self.__dict__ = state[1]

    def ferrers_diagram(self):
        """
        Return the Ferrers diagram of ``self``.

        EXAMPLES::

            sage: print(SkewPartition([[5,4,3,1],[3,3,1]]).ferrers_diagram())
               **
               *
             **
            *
            sage: print(SkewPartition([[5,4,3,1],[3,1]]).diagram())
               **
             ***
            ***
            *
            sage: SkewPartitions.options(diagram_str='#', convention='French')
            sage: print(SkewPartition([[5,4,3,1],[3,1]]).diagram())
            #
            ###
             ###
               ##
            sage: SkewPartitions.options._reset()
        """
        char, convention = self.parent().options('diagram_str', 'convention')

        if convention == "English":
            L = range(len(self[0]))
        else:
            L = reversed(range(len(self[0])))
        s = ""
        for i in L:
            if len(self[1]) > i:
                s += " "*self[1][i]
                s += char*(self[0][i]-self[1][i])
            else:
                s += char*self[0][i]
            s += "\n"
        return s[:-1]

    diagram = ferrers_diagram
    _repr_diagram = ferrers_diagram

    def pp(self):
        """
        Pretty-print ``self``.

        EXAMPLES::

            sage: SkewPartition([[5,4,3,1],[3,3,1]]).pp()
               **
               *
             **
            *
        """
        print(self.ferrers_diagram())

    def _ascii_art_(self):
        """
        TESTS::

            sage: ascii_art(SkewPartitions(3).list())
            [                        *   *   *    * ]
            [      **   **   *    *  *   *  *    *  ]
            [ ***, * , *  , **, ** , *, * , * , *   ]
            sage: SkewPartitions.options(diagram_str='#', convention='French')
            sage: ascii_art(SkewPartitions(3).list())
            [                        #  #   #   #   ]
            [      #   #    ##  ##   #   #  #    #  ]
            [ ###, ##,  ##,  #,   #, #,  #,  #,   # ]
            sage: SkewPartitions.options._reset()
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self.diagram().splitlines())

    def _unicode_art_(self):
        """
        .. WARNING::

            not working in presence of empty lines

        TESTS::

            sage: unicode_art(SkewPartitions(3).list())
            ⎡                             ┌┐   ┌┐   ┌┐    ┌┐ ⎤
            ⎢       ┌┬┐   ┌┬┐   ┌┐    ┌┐  ├┤   ├┤  ┌┼┘   ┌┼┘ ⎥
            ⎢ ┌┬┬┐  ├┼┘  ┌┼┴┘  ┌┼┤  ┌┬┼┘  ├┤  ┌┼┘  ├┤   ┌┼┘  ⎥
            ⎣ └┴┴┘, └┘ , └┘  , └┴┘, └┴┘ , └┘, └┘ , └┘ , └┘   ⎦
            sage: SkewPartitions.options.convention = "French"
            sage: unicode_art(SkewPartitions(3).list())
            ⎡                             ┌┐  ┌┐   ┌┐   ┌┐   ⎤
            ⎢       ┌┐   ┌┐    ┌┬┐  ┌┬┐   ├┤  └┼┐  ├┤   └┼┐  ⎥
            ⎢ ┌┬┬┐  ├┼┐  └┼┬┐  └┼┤  └┴┼┐  ├┤   ├┤  └┼┐   └┼┐ ⎥
            ⎣ └┴┴┘, └┴┘,  └┴┘,  └┘,   └┘, └┘,  └┘,  └┘,   └┘ ⎦
            sage: SkewPartitions.options._reset()

            sage: unicode_art(SkewPartition([[3,1],[2]]))
              ┌┐
            ┌┬┴┘
            └┘

        TESTS::

            sage: unicode_art(SkewPartition([[],[]]))
            ∅
        """
        from sage.typeset.unicode_art import UnicodeArt
        out, inn = self
        inn = inn + [0] * (len(out) - len(inn))
        if not any(self._list):
            return UnicodeArt('∅')
        if self.parent().options.convention == "French":
            s, t, b, l, r, tr, tl, br, bl, x, h = list(' ┴┬├┤┘└┐┌┼─')
        else:
            s, t, b, l, r, tr, tl, br, bl, x, h = list(' ┬┴├┤┐┌┘└┼─')

        # working with English conventions
        txt = [s * inn[0] + tl + t * (out[0] - inn[0] - 1) + tr]
        for i in range(len(out) - 1):
            o0 = out[i]
            o1 = out[i + 1]
            i0 = inn[i]
            i1 = inn[i + 1]

            if i0 == i1:
                start = ' ' * i1 + l
                d0 = 1
            else:
                start = ' ' * i1 + tl
                d0 = 0

            if o0 == o1:
                end = r
                d1 = 1
            else:
                end = br
                d1 = 0

            if i0 <= o1:
                middle = t * (i0 - i1 - 1 + d0)
                middle += x * (o1 - i0 + 1 - d0 - d1)
                middle += b * (o0 - o1 - 1 + d1)
            else:
                middle = t * (i0 - i1 - 1)
                middle += h * (i0 - o1 - 1)
                middle += b * (o0 - o1 - 1)

            txt += [start + middle + end]
        txt += [s * inn[-1] + bl + b * (out[-1] - inn[-1] - 1) + br]

        if self.parent().options.convention == "French":
            txt = list(reversed(txt))
        return UnicodeArt(txt, baseline=0)

    def inner(self):
        """
        Return the inner partition of ``self``.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).inner()
            [1, 1]
        """
        return self[1]

    def outer(self):
        """
        Return the outer partition of ``self``.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).outer()
            [3, 2, 1]
        """
        return self[0]

    def column_lengths(self):
        """
        Return the column lengths of ``self``.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).column_lengths()
            [1, 2, 1]
            sage: SkewPartition([[5,2,2,2],[2,1]]).column_lengths()
            [2, 3, 1, 1, 1]
        """
        return self.conjugate().row_lengths()

    def row_lengths(self):
        """
        Return the row lengths of ``self``.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).row_lengths()
            [2, 1, 1]
        """
        skp = self
        o = skp[0]
        i = skp[1]+[0]*(len(skp[0])-len(skp[1]))
        return [x[0]-x[1] for x in zip(o,i)]

    def size(self):
        """
        Return the size of ``self``.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).size()
            4
        """
        return sum(self.row_lengths())

    def is_connected(self):
        """
        Return ``True`` if ``self`` is a connected skew partition.

        A skew partition is said to be *connected* if for each pair of
        consecutive rows, there are at least two cells (one in each row)
        which have a common edge.

        EXAMPLES::

            sage: SkewPartition([[5,4,3,1],[3,3,1]]).is_connected()
            False
            sage: SkewPartition([[5,4,3,1],[3,1]]).is_connected()
            True
        """
        return self.is_overlap(1)

    def overlap(self):
        """
        Return the overlap of ``self``.

        The overlap of two consecutive rows in a skew partition is the
        number of pairs of cells (one in each row) that share a common
        edge. This number can be positive, zero, or negative.

        The overlap of a skew partition is the minimum of the overlap of
        the consecutive rows, or infinity in the case of at most one row.
        If the overlap is positive, then the skew partition is called
        *connected*.

        EXAMPLES::

            sage: SkewPartition([[],[]]).overlap()
            +Infinity
            sage: SkewPartition([[1],[]]).overlap()
            +Infinity
            sage: SkewPartition([[10],[]]).overlap()
            +Infinity
            sage: SkewPartition([[10],[2]]).overlap()
            +Infinity
            sage: SkewPartition([[10,1],[2]]).overlap()
            -1
            sage: SkewPartition([[10,10],[1]]).overlap()
            9
        """
        p,q = self
        if len(p) <= 1:
            from sage.rings.infinity import PlusInfinity
            return PlusInfinity()
        if len(q) == 0:
            return min(p)
        q = [q[0]] + list(q)
        return min(row_lengths_aux([p,q]))

    def is_overlap(self, n):
        r"""
        Return ``True`` if the overlap of ``self`` is at most ``n``.

        .. SEEALSO::

           :meth:`overlap`

        EXAMPLES::

            sage: SkewPartition([[5,4,3,1],[3,1]]).is_overlap(1)
            True
        """
        return n <= self.overlap()

    def is_ribbon(self):
        r"""
        Return ``True`` if and only if ``self`` is a ribbon.

        This means that if it has exactly one cell in each of `q`
        consecutive diagonals for some nonnegative integer `q`.

        EXAMPLES::

            sage: P = SkewPartition([[4,4,3,3],[3,2,2]])
            sage: P.pp()
               *
              **
              *
            ***
            sage: P.is_ribbon()
            True

            sage: P = SkewPartition([[4,3,3],[1,1]])
            sage: P.pp()
             ***
             **
            ***
            sage: P.is_ribbon()
            False

            sage: P = SkewPartition([[4,4,3,2],[3,2,2]])
            sage: P.pp()
               *
              **
              *
            **
            sage: P.is_ribbon()
            False

            sage: P = SkewPartition([[4,4,3,3],[4,2,2,1]])
            sage: P.pp()
            <BLANKLINE>
              **
              *
             **
            sage: P.is_ribbon()
            True

            sage: P = SkewPartition([[4,4,3,3],[4,2,2]])
            sage: P.pp()
            <BLANKLINE>
              **
              *
            ***
            sage: P.is_ribbon()
            True

            sage: SkewPartition([[2,2,1],[2,2,1]]).is_ribbon()
            True
        """
        lam = self[0]
        mu = self[1]
        l_out = len(lam)
        l_in = len(mu)
        mu += [0]*(l_out-l_in)

        if l_out == 0:
            return True
        else:
            # Find the least u for which lam[u]>mu[u], if it exists
            # If it does not exist then u will equal l_out
            u = 0
            u_test = True
            while u_test:
                if u >= l_out or lam[u] > mu[u]:
                    u_test = False
                else:
                    u += 1

            # Find the least v strictly greater than u for which
            # lam[v] != mu[v-1]+1
            v = u + 1
            v_test = True
            while v_test:
                if v >= l_out or lam[v] != mu[v-1] + 1:
                    v_test = False
                else:
                    v += 1

            # Check if lam[i]==mu[i] for all i >= v
            for i in range(v, l_out):
                if lam[i] != mu[i]:
                    return False

        return True

    def conjugate(self):
        """
        Return the conjugate of the skew partition skp.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2]]).conjugate()
            [3, 2, 1] / [1, 1]
        """
        return SkewPartition([x.conjugate() for x in self])

    def outer_corners(self):
        """
        Return a list of the outer corners of ``self``.

        These are corners that are contained inside of the shape.
        For the corners which are outside of the shape,
        use :meth:`outside_corners`.

        .. WARNING::

            In the case that ``self`` is an honest (rather than skew) partition,
            these are the :meth:`~sage.combinat.partition.Partition.corners`
            of the outer partition. In the language of [Sag2001]_ these would
            be the "inner corners" of the outer partition.

        .. SEEALSO::

            - :meth:`sage.combinat.skew_partition.SkewPartition.outside_corners`
            - :meth:`sage.combinat.partition.Partition.outside_corners`

        EXAMPLES::

            sage: SkewPartition([[4, 3, 1], [2]]).outer_corners()
            [(0, 3), (1, 2), (2, 0)]
        """
        return self.outer().corners()

    def inner_corners(self):
        """
        Return a list of the inner corners of ``self``.

        EXAMPLES::

            sage: SkewPartition([[4, 3, 1], [2]]).inner_corners()
            [(0, 2), (1, 0)]
            sage: SkewPartition([[4, 3, 1], []]).inner_corners()
            [(0, 0)]
        """
        inner = self.inner()
        outer = self.outer()
        if inner == []:
            if outer == []:
                return []
            else:
                return [(0,0)]
        icorners = [(0, inner[0])]
        nn = len(inner)
        for i in range(1,nn):
            if inner[i] != inner[i-1]:
                icorners += [ (i, inner[i]) ]

        icorners += [(nn, 0)]
        return icorners

    def cell_poset(self, orientation='SE'):
        """
        Return the Young diagram of ``self`` as a poset. The optional
        keyword variable ``orientation`` determines the order relation
        of the poset.

        The poset always uses the set of cells of the Young diagram
        of ``self`` as its ground set. The order relation of the poset
        depends on the ``orientation`` variable (which defaults to
        ``'SE'``). Concretely, ``orientation`` has to be specified to
        one of the strings ``'NW'``, ``'NE'``, ``'SW'``, and ``'SE'``,
        standing for "northwest", "northeast", "southwest" and
        "southeast", respectively. If ``orientation`` is ``'SE'``, then
        the order relation of the poset is such that a cell `u` is
        greater or equal to a cell `v` in the poset if and only if `u`
        lies weakly southeast of `v` (this means that `u` can be
        reached from `v` by a sequence of south and east steps; the
        sequence is allowed to consist of south steps only, or of east
        steps only, or even be empty). Similarly the order relation is
        defined for the other three orientations. The Young diagram is
        supposed to be drawn in English notation.

        The elements of the poset are the cells of the Young diagram
        of ``self``, written as tuples of zero-based coordinates (so
        that `(3, 7)` stands for the `8`-th cell of the `4`-th row,
        etc.).

        EXAMPLES::

            sage: # needs sage.graphs
            sage: p = SkewPartition([[3,3,1], [2,1]])
            sage: Q = p.cell_poset(); Q
            Finite poset containing 4 elements
            sage: sorted(Q)
            [(0, 2), (1, 1), (1, 2), (2, 0)]
            sage: sorted(Q.maximal_elements())
            [(1, 2), (2, 0)]
            sage: sorted(Q.minimal_elements())
            [(0, 2), (1, 1), (2, 0)]
            sage: sorted(Q.upper_covers((1, 1)))
            [(1, 2)]
            sage: sorted(Q.upper_covers((0, 2)))
            [(1, 2)]

            sage: # needs sage.graphs
            sage: P = p.cell_poset(orientation="NW"); P
            Finite poset containing 4 elements
            sage: sorted(P)
            [(0, 2), (1, 1), (1, 2), (2, 0)]
            sage: sorted(P.minimal_elements())
            [(1, 2), (2, 0)]
            sage: sorted(P.maximal_elements())
            [(0, 2), (1, 1), (2, 0)]
            sage: sorted(P.upper_covers((1, 2)))
            [(0, 2), (1, 1)]

            sage: # needs sage.graphs
            sage: R = p.cell_poset(orientation="NE"); R
            Finite poset containing 4 elements
            sage: sorted(R)
            [(0, 2), (1, 1), (1, 2), (2, 0)]
            sage: R.maximal_elements()
            [(0, 2)]
            sage: R.minimal_elements()
            [(2, 0)]
            sage: R.upper_covers((2, 0))
            [(1, 1)]
            sage: sorted([len(R.upper_covers(v)) for v in R])
            [0, 1, 1, 1]

        TESTS:

        We check that the posets are really what they should be for size
        up to `6`::

            sage: def check_NW(n):
            ....:     for p in SkewPartitions(n):
            ....:         P = p.cell_poset(orientation="NW")
            ....:         for c in p.cells():
            ....:             for d in p.cells():
            ....:                 if P.le(c, d) != (c[0] >= d[0]
            ....:                                   and c[1] >= d[1]):
            ....:                     return False
            ....:     return True
            sage: all( check_NW(n) for n in range(7) )                                  # needs sage.graphs
            True

            sage: def check_NE(n):
            ....:     for p in SkewPartitions(n):
            ....:         P = p.cell_poset(orientation="NE")
            ....:         for c in p.cells():
            ....:             for d in p.cells():
            ....:                 if P.le(c, d) != (c[0] >= d[0]
            ....:                                   and c[1] <= d[1]):
            ....:                     return False
            ....:     return True
            sage: all( check_NE(n) for n in range(7) )                                  # needs sage.graphs
            True

            sage: def test_duality(n, ori1, ori2):
            ....:     for p in SkewPartitions(n):
            ....:         P = p.cell_poset(orientation=ori1)
            ....:         Q = p.cell_poset(orientation=ori2)
            ....:         for c in p.cells():
            ....:             for d in p.cells():
            ....:                 if P.lt(c, d) != Q.lt(d, c):
            ....:                     return False
            ....:     return True
            sage: all( test_duality(n, "NW", "SE") for n in range(7) )                  # needs sage.graphs
            True
            sage: all( test_duality(n, "NE", "SW") for n in range(7) )                  # needs sage.graphs
            True
            sage: all( test_duality(n, "NE", "SE") for n in range(4) )                  # needs sage.graphs
            False
        """
        from sage.combinat.posets.posets import Poset
        # Getting the cover relations seems hard, so let's just compute
        # the comparison function.
        if orientation == "NW":
            def poset_le(u, v):
                return u[0] >= v[0] and u[1] >= v[1]
        elif orientation == "NE":
            def poset_le(u, v):
                return u[0] >= v[0] and u[1] <= v[1]
        elif orientation == "SE":
            def poset_le(u, v):
                return u[0] <= v[0] and u[1] <= v[1]
        elif orientation == "SW":
            def poset_le(u, v):
                return u[0] <= v[0] and u[1] >= v[1]
        return Poset((self.cells(), poset_le))

    def frobenius_rank(self):
        r"""
        Return the Frobenius rank of the skew partition ``self``.

        The Frobenius rank of a skew partition `\lambda / \mu` can be
        defined in various ways. The quickest one is probably the
        following: Writing `\lambda` as
        `(\lambda_1, \lambda_2, \cdots , \lambda_N)`, and writing `\mu`
        as `(\mu_1, \mu_2, \cdots , \mu_N)`, we define the Frobenius
        rank of `\lambda / \mu` to be the number of all
        `1 \leq i \leq N` such that

        .. MATH::

            \lambda_i - i
            \not\in \{ \mu_1 - 1, \mu_2 - 2, \cdots , \mu_N - N \}.

        In other words, the Frobenius rank of `\lambda / \mu` is the
        number of rows in the Jacobi-Trudi matrix of `\lambda / \mu`
        which don't contain `h_0`. Further definitions have been
        considered in [Sta2002]_ (where Frobenius rank is just being
        called rank).

        If `\mu` is the empty shape, then the Frobenius rank of
        `\lambda / \mu` is just the usual Frobenius rank of the
        partition `\lambda` (see
        :meth:`~sage.combinat.partition.Partition.frobenius_rank()`).

        EXAMPLES::

            sage: SkewPartition([[8,8,7,4], [4,1,1]]).frobenius_rank()
            4
            sage: SkewPartition([[2,1], [1]]).frobenius_rank()
            2
            sage: SkewPartition([[2,1,1], [1]]).frobenius_rank()
            2
            sage: SkewPartition([[2,1,1], [1,1]]).frobenius_rank()
            2
            sage: SkewPartition([[5,4,3,2], [2,1,1]]).frobenius_rank()
            3
            sage: SkewPartition([[4,2,1], [3,1,1]]).frobenius_rank()
            2
            sage: SkewPartition([[4,2,1], [3,2,1]]).frobenius_rank()
            1

        If the inner shape is empty, then the Frobenius rank of the skew
        partition is just the standard Frobenius rank of the partition::

            sage: all( SkewPartition([lam, Partition([])]).frobenius_rank()
            ....:      == lam.frobenius_rank() for i in range(6)
            ....:      for lam in Partitions(i) )
            True

        If the inner and outer shapes are equal, then the Frobenius rank
        is zero::

            sage: all( SkewPartition([lam, lam]).frobenius_rank() == 0
            ....:      for i in range(6) for lam in Partitions(i) )
            True
        """
        N = len(self[0])
        mu_betas = [x - j for j, x in enumerate(self[1])]
        mu_betas.extend([- j for j in range(len(self[1]), N)])
        res = 0
        for i, x in enumerate(self[0]):
            if (x - i) not in mu_betas:
                res += 1
        return res

    def cells(self):
        """
        Return the coordinates of the cells of ``self``. Coordinates are
        given as ``(row-index, column-index)`` and are 0 based.

        EXAMPLES::

            sage: SkewPartition([[4, 3, 1], [2]]).cells()
            [(0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (2, 0)]
            sage: SkewPartition([[4, 3, 1], []]).cells()
            [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (2, 0)]
            sage: SkewPartition([[2], []]).cells()
            [(0, 0), (0, 1)]
        """
        outer = self.outer()
        inner = self.inner()[:]
        inner += [0] * (len(outer) - len(inner))
        return [(i, j) for i, outi in enumerate(outer)
                for j in range(inner[i], outi)]

    def to_list(self):
        """
        Return ``self`` as a list of lists.

        EXAMPLES::

            sage: s = SkewPartition([[4,3,1],[2]])
            sage: s.to_list()
            [[4, 3, 1], [2]]
            sage: type(s.to_list())
            <class 'list'>
        """
        return [list(r) for r in list(self)]

    def to_dag(self, format='string'):
        """
        Return a directed acyclic graph corresponding to the skew
        partition ``self``.

        The directed acyclic graph corresponding to a skew partition
        `p` is the digraph whose vertices are the cells of `p`, and
        whose edges go from each cell to its lower and right
        neighbors (in English notation).

        INPUT:

        - ``format`` -- either ``'string'`` or ``'tuple'`` (default:
          ``'string'``); determines whether the vertices of the
          resulting dag will be strings or 2-tuples of coordinates

        EXAMPLES::

            sage: # needs sage.graphs
            sage: dag = SkewPartition([[3, 3, 1], [1, 1]]).to_dag()
            sage: dag.edges(sort=True)
            [('0,1', '0,2', None),
            ('0,1', '1,1', None),
            ('0,2', '1,2', None),
            ('1,1', '1,2', None)]
            sage: dag.vertices(sort=True)
            ['0,1', '0,2', '1,1', '1,2', '2,0']
            sage: dag = SkewPartition([[3, 2, 1], [1, 1]]).to_dag(format='tuple')
            sage: dag.edges(sort=True)
            [((0, 1), (0, 2), None), ((0, 1), (1, 1), None)]
            sage: dag.vertices(sort=True)
            [(0, 1), (0, 2), (1, 1), (2, 0)]
        """
        outer = list(self.outer())
        inner = list(self.inner())
        inner += [0] * (len(outer) - len(inner))

        G = DiGraph()
        for i, outer_i in enumerate(outer):
            for j in range(inner[i], outer_i):
                if format == "string":
                    string = "%d,%d" % (i, j)
                else:
                    string = (i, j)
                G.add_vertex(string)
                #Check to see if there is a node to the right
                if j != outer_i - 1:
                    if format == "string":
                        newstring = "%d,%d" % (i, j + 1)
                    else:
                        newstring = (i, j + 1)
                    G.add_edge(string, newstring)

                #Check to see if there is anything below
                if i != len(outer) - 1:
                    if outer[i+1] > j:
                        if format == "string":
                            newstring = "%d,%d" % (i + 1, j)
                        else:
                            newstring = (i + 1, j)
                        G.add_edge(string, newstring)
        return G

    def quotient(self, k):
        """
        The quotient map extended to skew partitions.

        EXAMPLES::

            sage: SkewPartition([[3, 3, 2, 1], [2, 1]]).quotient(2)
            [[3] / [], [] / []]
        """
        ## k-th element is the skew partition built using the k-th partition of the
        ## k-quotient of the outer and the inner partition.
        ## This bijection is only defined if the inner and the outer partition
        ## have the same core
        if self.inner().core(k) == self.outer().core(k):
            rqinner = self.inner().quotient(k)
            rqouter = self.outer().quotient(k)
            return [ SkewPartitions()([rqouter[i],rqinner[i]]) for i in range(k) ]
        else:
            raise ValueError("quotient map is only defined for skew partitions with inner and outer partitions having the same core")

    def rows_intersection_set(self):
        r"""
        Return the set of cells in the rows of the outer shape of
        ``self`` which rows intersect the skew diagram of ``self``.

        EXAMPLES::

            sage: skp = SkewPartition([[3,2,1],[2,1]])
            sage: cells = Set([ (0,0), (0, 1), (0,2), (1, 0), (1, 1), (2, 0)])
            sage: skp.rows_intersection_set() == cells
            True
        """
        outer = self.outer()
        inner = self.inner()
        inner += [0] * (len(outer) - len(inner))

        res = [(i, j) for i, outi in enumerate(outer)
               for j in range(outi)
               if outi != inner[i]]
        return Set(res)

    def columns_intersection_set(self):
        """
        Return the set of cells in the columns of the outer shape of
        ``self`` which columns intersect the skew diagram of ``self``.

        EXAMPLES::

            sage: skp = SkewPartition([[3,2,1],[2,1]])
            sage: cells = Set([ (0,0), (0, 1), (0,2), (1, 0), (1, 1), (2, 0)])
            sage: skp.columns_intersection_set() == cells
            True
        """
        res = [ (x[1], x[0]) for x in self.conjugate().rows_intersection_set()]
        return Set(res)

    def pieri_macdonald_coeffs(self):
        """
        Computation of the coefficients which appear in the Pieri formula
        for Macdonald polynomials given in his book ( Chapter 6.6 formula
        6.24(ii) )

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2,1]]).pieri_macdonald_coeffs()
            1
            sage: SkewPartition([[3,2,1],[2,2]]).pieri_macdonald_coeffs()
            (q^2*t^3 - q^2*t - t^2 + 1)/(q^2*t^3 - q*t^2 - q*t + 1)
            sage: SkewPartition([[3,2,1],[2,2,1]]).pieri_macdonald_coeffs()
            (q^6*t^8 - q^6*t^6 - q^4*t^7 - q^5*t^5 + q^4*t^5 - q^3*t^6 + q^5*t^3 + 2*q^3*t^4 + q*t^5 - q^3*t^2 + q^2*t^3 - q*t^3 - q^2*t - t^2 + 1)/(q^6*t^8 - q^5*t^7 - q^5*t^6 - q^4*t^6 + q^3*t^5 + 2*q^3*t^4 + q^3*t^3 - q^2*t^2 - q*t^2 - q*t + 1)
            sage: SkewPartition([[3,3,2,2],[3,2,2,1]]).pieri_macdonald_coeffs()
            (q^5*t^6 - q^5*t^5 + q^4*t^6 - q^4*t^5 - q^4*t^3 + q^4*t^2 - q^3*t^3 - q^2*t^4 + q^3*t^2 + q^2*t^3 - q*t^4 + q*t^3 + q*t - q + t - 1)/(q^5*t^6 - q^4*t^5 - q^3*t^4 - q^3*t^3 + q^2*t^3 + q^2*t^2 + q*t - 1)
        """

        set_prod = self.rows_intersection_set() - self.columns_intersection_set()
        res = 1
        for s in set_prod:
            res *= self.inner().arms_legs_coeff(s[0],s[1])
            res /= self.outer().arms_legs_coeff(s[0],s[1])
        return res

    def k_conjugate(self, k):
        """
        Return the `k`-conjugate of the skew partition.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2,1]]).k_conjugate(3)
            [2, 1, 1, 1, 1] / [2, 1]
            sage: SkewPartition([[3,2,1],[2,1]]).k_conjugate(4)
            [2, 2, 1, 1] / [2, 1]
            sage: SkewPartition([[3,2,1],[2,1]]).k_conjugate(5)
            [3, 2, 1] / [2, 1]
        """
        return SkewPartition([ self.outer().k_conjugate(k), self.inner().k_conjugate(k) ])

    def jacobi_trudi(self):
        """
        Return the Jacobi-Trudi matrix of ``self``.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2,1]]).jacobi_trudi()                         # needs sage.modules
            [h[1]    0    0]
            [h[3] h[1]    0]
            [h[5] h[3] h[1]]
            sage: SkewPartition([[4,3,2],[2,1]]).jacobi_trudi()                         # needs sage.modules
            [h[2]  h[]    0]
            [h[4] h[2]  h[]]
            [h[6] h[4] h[2]]
        """
        p = self.outer()
        q = self.inner()
        from sage.combinat.sf.sf import SymmetricFunctions
        nn = len(p)
        if nn == 0:
            return MatrixSpace(SymmetricFunctions(QQ).homogeneous(), 0)(0)
        h = SymmetricFunctions(QQ).homogeneous()
        H = MatrixSpace(h, nn)

        q = q + [0] * (nn - len(q))
        m = []
        for i in range(1,nn+1):
            row = []
            for j in range(1,nn+1):
                v = p[j-1]-q[i-1]-j+i
                if v < 0:
                    row.append(h.zero())
                elif v == 0:
                    row.append(h([]))
                else:
                    row.append(h([v]))
            m.append(row)
        return H(m)

    def outside_corners(self):
        r"""
        Return the outside corners of ``self``.

        The outside corners are corners which are outside of the shape. This
        should not be confused with :meth:`outer_corners` which consists of
        corners inside the shape. It returns a result analogous to the
        ``.outside_corners()`` method on (non-skew) ``Partitions``.

        .. SEEALSO::

            - :meth:`sage.combinat.skew_partition.SkewPartition.outer_corners`
            - :meth:`sage.combinat.partition.Partition.outside_corners`

        EXAMPLES::

            sage: mu = SkewPartition([[3,2,1],[2,1]])
            sage: mu.pp()
              *
             *
            *
            sage: mu.outside_corners()
            [(0, 3), (1, 2), (2, 1), (3, 0)]
        """
        return self.outer().outside_corners()

    def specht_module(self, base_ring=None):
        r"""
        Return the Specht module corresponding to ``self``.

        EXAMPLES::

            sage: mu = SkewPartition([[3,2,1], [2]])
            sage: SM = mu.specht_module(QQ)                                             # needs sage.modules
            sage: s = SymmetricFunctions(QQ).s()                                        # needs sage.modules
            sage: s(SM.frobenius_image())                                               # needs sage.modules
            s[2, 1, 1] + s[2, 2] + s[3, 1]

        We verify that the Frobenius image is the corresponding
        skew Schur function::

            sage: s[3,2,1].skew_by(s[2])                                                # needs sage.modules
            s[2, 1, 1] + s[2, 2] + s[3, 1]

        ::

            sage: mu = SkewPartition([[4,2,1], [2,1]])
            sage: SM = mu.specht_module(QQ)                                             # needs sage.modules
            sage: s(SM.frobenius_image())                                               # needs sage.modules
            s[2, 1, 1] + s[2, 2] + 2*s[3, 1] + s[4]
            sage: s(mu)                                                                 # needs sage.modules
            s[2, 1, 1] + s[2, 2] + 2*s[3, 1] + s[4]
        """
        from sage.combinat.specht_module import SpechtModule
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ
        R = SymmetricGroupAlgebra(base_ring, self.size())
        return SpechtModule(R, self.cells())

    def specht_module_dimension(self, base_ring=None):
        r"""
        Return the dimension of the Specht module corresponding to ``self``.

        This is equal to the number of standard (skew) tableaux of
        shape ``self``.

        EXAMPLES::

            sage: mu = SkewPartition([[3,2,1], [2]])
            sage: mu.specht_module_dimension()                                          # needs sage.modules
            8
            sage: mu.specht_module_dimension(GF(2))                                     # needs sage.modules sage.rings.finite_rings
            8
        """
        from sage.categories.fields import Fields
        if base_ring is None or (base_ring in Fields() and base_ring.characteristic() == 0):
            from sage.combinat.skew_tableau import StandardSkewTableaux
            return StandardSkewTableaux(self).cardinality()
        from sage.combinat.specht_module import specht_module_rank
        return specht_module_rank(self, base_ring)


def row_lengths_aux(skp):
    """
    EXAMPLES::

        sage: from sage.combinat.skew_partition import row_lengths_aux
        sage: row_lengths_aux([[5,4,3,1],[3,3,1]])
        [2, 1, 2]
        sage: row_lengths_aux([[5,4,3,1],[3,1]])
        [2, 3]
    """
    if skp[0] == []:
        return []
    else:
        return [x[0] - x[1] for x in zip(skp[0], skp[1])]


class SkewPartitions(UniqueRepresentation, Parent):
    """
    Skew partitions.

    .. WARNING::

        The iterator of this class only yields skew partitions which
        are reduced, in the sense that there are no empty rows
        before the last nonempty row, and there are no empty columns
        before the last nonempty column.

    EXAMPLES::

        sage: SkewPartitions(4)
        Skew partitions of 4
        sage: SkewPartitions(4).cardinality()
        28
        sage: SkewPartitions(row_lengths=[2,1,2])
        Skew partitions with row lengths [2, 1, 2]
        sage: SkewPartitions(4, overlap=2)
        Skew partitions of 4 with a minimum overlap of 2
        sage: SkewPartitions(4, overlap=2).list()
        [[4] / [], [2, 2] / []]
    """
    @staticmethod
    def __classcall_private__(self, n=None, row_lengths=None, overlap=0):
        """
        Return the correct parent based upon the input.

        EXAMPLES::

            sage: SP1 = SkewPartitions(row_lengths=(2,1,2))
            sage: SP2 = SkewPartitions(row_lengths=[2,1,2])
            sage: SP1 is SP2
            True
        """
        if n is not None:
            if row_lengths is not None:
                raise ValueError("you can only specify one of n or row_lengths")
            return SkewPartitions_n(n, overlap)
        elif row_lengths is not None:
            return SkewPartitions_rowlengths(row_lengths, overlap)
        else:
            return SkewPartitions_all()

    def __init__(self, is_infinite=False):
        """
        TESTS::

            sage: S = SkewPartitions()
            sage: TestSuite(S).run()
        """
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

    # add options to class
    class options(GlobalOptions):
        """
        Set and display the options for elements of the skew partition
        classes.  If no parameters are set, then the function returns a copy of
        the options dictionary.

        The ``options`` to skew partitions can be accessed as the method
        :obj:`SkewPartitions.options` of :class:`SkewPartitions` and
        related parent classes.

        @OPTIONS@

        EXAMPLES::

            sage: SP = SkewPartition([[4,2,2,1], [3, 1, 1]])
            sage: SP
            [4, 2, 2, 1] / [3, 1, 1]
            sage: SkewPartitions.options.display="lists"
            sage: SP
            [[4, 2, 2, 1], [3, 1, 1]]

        Changing the ``convention`` for skew partitions also changes the
        ``convention`` option for partitions and tableaux and vice versa::

            sage: SkewPartitions.options(display='diagram', convention='French')
            sage: SP
            *
             *
             *
               *
            sage: T = Tableau([[1,2,3],[4,5]])
            sage: T.pp()
              4  5
              1  2  3
            sage: P = Partition([4, 2, 2, 1])
            sage: P.pp()
            *
            **
            **
            ****
            sage: Tableaux.options.convention="english"
            sage: SP
               *
             *
             *
            *
            sage: T.pp()
              1  2  3
              4  5
            sage: SkewPartitions.options._reset()
        """
        NAME = 'SkewPartitions'
        module = 'sage.combinat.skew_partition'
        display = dict(default='quotient',
                     description='Specifies how skew partitions should be printed',
                     values=dict(lists='displayed as a pair of lists',
                                 quotient='displayed as a quotient of partitions',
                                 diagram='as a skew Ferrers diagram'),
                     alias=dict(array='diagram', ferrers_diagram='diagram',
                                young_diagram='diagram', pair='lists'),
                     case_sensitive=False)
        latex = dict(default='young_diagram',
                   description='Specifies how skew partitions should be latexed',
                   values=dict(diagram='latex as a skew Ferrers diagram',
                               young_diagram='latex as a skew Young diagram',
                               marked='latex as a partition where the skew shape is marked'),
                   alias=dict(array='diagram', ferrers_diagram='diagram'),
                   case_sensitive=False)
        diagram_str = dict(link_to=(Partitions.options,'diagram_str'))
        latex_diagram_str = dict(link_to=(Partitions.options,'latex_diagram_str'))
        latex_marking_str = dict(default='X',
                         description='The character used to marked the deleted cells when latexing marked partitions',
                         checker=lambda char: isinstance(char, str))
        convention = dict(link_to=(Tableaux.options,'convention'))
        notation = dict(alt_name='convention')

    Element = SkewPartition

    def _element_constructor_(self, skp):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: S = SkewPartitions()
            sage: S([[3,1], [1]])
            [3, 1] / [1]
        """
        return self.element_class(self, skp)

    def __contains__(self, x):
        """
        TESTS::

            sage: [[], []] in SkewPartitions()
            True
            sage: [[], [1]] in SkewPartitions()
            False
            sage: [[], [-1]] in SkewPartitions()
            False
            sage: [[], [0]] in SkewPartitions()
            True
            sage: [[3,2,1],[]] in SkewPartitions()
            True
            sage: [[3,2,1],[1]] in SkewPartitions()
            True
            sage: [[3,2,1],[2]] in SkewPartitions()
            True
            sage: [[3,2,1],[3]] in SkewPartitions()
            True
            sage: [[3,2,1],[4]] in SkewPartitions()
            False
            sage: [[3,2,1],[1,1]] in SkewPartitions()
            True
            sage: [[3,2,1],[1,2]] in SkewPartitions()
            False
            sage: [[3,2,1],[2,1]] in SkewPartitions()
            True
            sage: [[3,2,1],[2,2]] in SkewPartitions()
            True
            sage: [[3,2,1],[3,2]] in SkewPartitions()
            True
            sage: [[3,2,1],[1,1,1]] in SkewPartitions()
            True
            sage: [[7, 4, 3, 2], [8, 2, 1]] in SkewPartitions()
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions()
            True
            sage: [[4,2,1],[1,1,1,1]] in SkewPartitions()
            False
            sage: [[1,1,1,0],[1,1,0,0]] in SkewPartitions()
            True
        """
        if isinstance(x, SkewPartition):
            return True

        try:
            if len(x) != 2:
                return False
        except TypeError:
            return False

        p = _Partitions
        if x[0] not in p:
            return False
        if x[1] not in p:
            return False

        if not p(x[0]).contains(p(x[1])):
            return False

        return True

    def from_row_and_column_length(self, rowL, colL):
        """
        Construct a partition from its row lengths and column lengths.

        INPUT:

        - ``rowL`` -- a composition or a list of positive integers

        - ``colL`` -- a composition or a list of positive integers

        OUTPUT:

        - If it exists the unique skew-partitions with row lengths ``rowL``
          and column lengths ``colL``.
        - Raise a :exc:`ValueError` if ``rowL`` and ``colL`` are not compatible.

        EXAMPLES::

            sage: S = SkewPartitions()
            sage: print(S.from_row_and_column_length([3,1,2,2],[2,3,1,1,1]).diagram())
              ***
             *
            **
            **
            sage: S.from_row_and_column_length([],[])
            [] / []
            sage: S.from_row_and_column_length([1],[1])
            [1] / []
            sage: S.from_row_and_column_length([2,1],[2,1])
            [2, 1] / []
            sage: S.from_row_and_column_length([1,2],[1,2])
            [2, 2] / [1]
            sage: S.from_row_and_column_length([1,2],[1,3])
            Traceback (most recent call last):
            ...
            ValueError: sum mismatch: [1, 2] and [1, 3]
            sage: S.from_row_and_column_length([3,2,1,2],[2,3,1,1,1])
            Traceback (most recent call last):
            ...
            ValueError: incompatible row and column length : [3, 2, 1, 2] and [2, 3, 1, 1, 1]

        .. WARNING::

            If some rows and columns have length zero, there is no way to retrieve
            unambiguously the skew partition. We therefore raise
            a :exc:`ValueError`.
            For examples here are two skew partitions with the same row and column
            lengths::

                sage: skp1 = SkewPartition([[2,2],[2,2]])
                sage: skp2 = SkewPartition([[2,1],[2,1]])
                sage: skp1.row_lengths(), skp1.column_lengths()
                ([0, 0], [0, 0])
                sage: skp2.row_lengths(), skp2.column_lengths()
                ([0, 0], [0, 0])
                sage: SkewPartitions().from_row_and_column_length([0,0], [0,0])
                Traceback (most recent call last):
                ...
                ValueError: row and column length must be positive

        TESTS::

            sage: all(SkewPartitions().from_row_and_column_length(p.row_lengths(), p.column_lengths()) == p
            ....:       for i in range(7) for p in SkewPartitions(i))
            True
        """
        if sum(rowL) != sum(colL):
            raise ValueError(f"sum mismatch: {rowL} and {colL}")
        if not all(i > 0 for i in rowL) or not all(i > 0 for i in colL):
            raise ValueError("row and column length must be positive")
        if not rowL:
            return self.element_class(self, [[], []])
        colL_new = colL[:]
        resIn = []
        resOut = []
        inPOld = len(colL)
        for row in rowL:
            inP = len(colL_new) - row
            if inP < 0 or inP > inPOld:
                raise ValueError("incompatible row and column length : %s and %s" % (rowL, colL))
            inPOld = inP
            resIn.append(inP)
            resOut.append(len(colL_new))
            for iCol in range(inP, len(colL_new)):
                colL_new[iCol] -= 1
                if colL_new[iCol] < 0:
                    raise ValueError("incompatible row and column length : %s and %s" % (rowL, colL))
            while colL_new and colL_new[-1] == 0:
                colL_new.pop()
        return self.element_class(self, [resOut, [x for x in resIn if x]])


class SkewPartitions_all(SkewPartitions):
    """
    Class of all skew partitions.
    """

    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SkewPartitions()
            sage: TestSuite(S).run()
        """
        SkewPartitions.__init__(self, True)

    def _repr_(self):
        """
        TESTS::

            sage: SkewPartitions()
            Skew partitions
        """
        return "Skew partitions"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SP = SkewPartitions()
            sage: it = SP.__iter__()
            sage: [next(it) for x in range(10)]
            [[] / [],
             [1] / [],
             [2] / [],
             [1, 1] / [],
             [2, 1] / [1],
             [3] / [],
             [2, 1] / [],
             [3, 1] / [1],
             [2, 2] / [1],
             [3, 2] / [2]]
        """
        n = 0
        while True:
            for p in SkewPartitions_n(n):
                yield self.element_class(self, p)
            n += 1


class SkewPartitions_n(SkewPartitions):
    """
    The set of skew partitions of ``n`` with overlap at least
    ``overlap`` and no empty row.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``overlap`` -- integer (default: `0`)

    Caveat: this set is stable under conjugation only for ``overlap`` equal
    to 0 or 1. What exactly happens for negative overlaps is not yet
    well specified and subject to change (we may want to
    introduce vertical overlap constraints as well).

    .. TODO::

        As is, this set is essentially the composition of
        ``Compositions(n)`` (which give the row lengths) and
        ``SkewPartition(n, row_lengths=...)``, and one would want to
        "inherit" list and cardinality from this composition.
    """
    @staticmethod
    def __classcall_private__(cls, n, overlap=0):
        """
        Normalize input so we have a unique representation.

        EXAMPLES::

            sage: S = SkewPartitions(3, overlap=1)
            sage: S2 = SkewPartitions(int(3), overlap='connected')
            sage: S is S2
            True
        """
        if overlap == 'connected':
            overlap = 1
        return super().__classcall__(cls, n, overlap)

    def __init__(self, n, overlap):
        """
        Return the set of the skew partitions of ``n`` with overlap
        at least ``overlap``, and no empty row.

        The iteration order is not specified yet.

        Caveat: this set is stable under conjugation only for overlap=
        0 or 1. What exactly happens for negative overlaps is not yet
        well specified, and subject to change (we may want to
        introduce vertical overlap constraints as well). ``overlap`` would
        also better be named ``min_overlap``.

        Todo: as is, this set is essentially the composition of
        ``Compositions(n)`` (which give the row lengths) and
        ``SkewPartition(n, row_lengths=...)``, and one would want to
        "inherit" list and cardinality from this composition.

        INPUT:

        - ``n`` -- nonnegative integer
        - ``overlap`` -- integer

        TESTS::

            sage: S = SkewPartitions(3)
            sage: TestSuite(S).run()
            sage: S = SkewPartitions(3, overlap=1)
            sage: TestSuite(S).run()
        """
        self.n = n
        self.overlap = overlap
        SkewPartitions.__init__(self, False)

    def __contains__(self, x):
        """
        TESTS::

            sage: [[],[]] in SkewPartitions(0)
            True
            sage: [[3,2,1], []] in SkewPartitions(6)
            True
            sage: [[3,2,1], []] in SkewPartitions(7)
            False
            sage: [[3,2,1], []] in SkewPartitions(5)
            False
            sage: [[7, 4, 3, 2], [8, 2, 1]] in SkewPartitions(8)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(5)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(5, overlap=-1)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap=-1)
            True
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap=0)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap='connected')
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap=-2)
            True
        """
        return x in SkewPartitions() \
            and sum(x[0])-sum(x[1]) == self.n \
            and self.overlap <= SkewPartition(x).overlap()

    def _repr_(self):
        """
        TESTS::

            sage: SkewPartitions(3)
            Skew partitions of 3
            sage: SkewPartitions(3, overlap=1)
            Skew partitions of 3 with a minimum overlap of 1
        """
        string = "Skew partitions of %s" % self.n
        if self.overlap:
            string += " with a minimum overlap of %s" % self.overlap
        return string

    def _count_slide(self, co, overlap=0):
        """
        Return the number of skew partitions related to the composition
        ``co`` by 'sliding'. The composition ``co`` is the list of row
        lengths of the skew partition.

        EXAMPLES::

            sage: s = SkewPartitions(3)
            sage: s._count_slide([2,1])
            2
            sage: [ sp for sp in s if sp.row_lengths() == [2,1] ]
            [[2, 1] / [], [3, 1] / [1]]
            sage: s = SkewPartitions(3, overlap=1)
            sage: s._count_slide([2,1], overlap=1)
            1
            sage: [ sp for sp in s if sp.row_lengths() == [2,1] ]
            [[2, 1] / []]
        """
        result = 1
        shift = 1 - overlap
        for i in range(len(co) - 1):
            result *= min(co[i], co[i + 1]) + shift
        return result

    def cardinality(self):
        """
        Return the number of skew partitions of the integer `n`
        (with given overlap, if specified; and with no empty rows before
        the last row).

        EXAMPLES::

            sage: SkewPartitions(0).cardinality()
            1
            sage: SkewPartitions(4).cardinality()
            28
            sage: SkewPartitions(5).cardinality()
            87
            sage: SkewPartitions(4, overlap=1).cardinality()
            9
            sage: SkewPartitions(5, overlap=1).cardinality()
            20
            sage: s = SkewPartitions(5, overlap=-1)
            sage: s.cardinality() == len(s.list())
            True
        """
        if self.n == 0:
            return ZZ.one()

        if self.overlap > 0:
            gg = Compositions(self.n, min_part=max(1, self.overlap))
        else:
            gg = Compositions(self.n)

        sum_a = 0
        for co in gg:
            sum_a += self._count_slide(co, overlap=self.overlap)

        return ZZ(sum_a)

    def __iter__(self):
        """
        Iterate through the skew partitions of `n`
        (with given overlap, if specified; and with no empty rows before
        the last row).

        EXAMPLES::

            sage: SkewPartitions(3).list()
            [[3] / [],
             [2, 1] / [],
             [3, 1] / [1],
             [2, 2] / [1],
             [3, 2] / [2],
             [1, 1, 1] / [],
             [2, 2, 1] / [1, 1],
             [2, 1, 1] / [1],
             [3, 2, 1] / [2, 1]]

            sage: SkewPartitions(3, overlap=0).list()
            [[3] / [],
             [2, 1] / [],
             [3, 1] / [1],
             [2, 2] / [1],
             [3, 2] / [2],
             [1, 1, 1] / [],
             [2, 2, 1] / [1, 1],
             [2, 1, 1] / [1],
             [3, 2, 1] / [2, 1]]
            sage: SkewPartitions(3, overlap=1).list()
            [[3] / [],
             [2, 1] / [],
             [2, 2] / [1],
             [1, 1, 1] / []]
            sage: SkewPartitions(3, overlap=2).list()
            [[3] / []]
            sage: SkewPartitions(3, overlap=3).list()
            [[3] / []]
            sage: SkewPartitions(3, overlap=4).list()
            []
        """
        for co in Compositions(self.n, min_part=max(1, self.overlap)):
            for sp in SkewPartitions(row_lengths=co, overlap=self.overlap):
                yield self.element_class(self, sp)

######################################
# Skew Partitions (from row lengths) #
######################################


class SkewPartitions_rowlengths(SkewPartitions):
    """
    All skew partitions with given row lengths.
    """
    @staticmethod
    def __classcall_private__(cls, co, overlap=0):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = SkewPartitions(row_lengths=[2,1], overlap=1)
            sage: S2 = SkewPartitions(row_lengths=(2,1), overlap='connected')
            sage: S is S2
            True
        """
        co = Compositions()(co)
        if overlap == 'connected':
            overlap = 1
        return super().__classcall__(cls, co, overlap)

    def __init__(self, co, overlap):
        """
        TESTS::

            sage: S = SkewPartitions(row_lengths=[2,1])
            sage: TestSuite(S).run()
        """
        self.co = co
        self.overlap = overlap
        SkewPartitions.__init__(self, False)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[4,3,1],[2]] in SkewPartitions(row_lengths=[2,3,1])
            True
            sage: [[4,3,1],[2]] in SkewPartitions(row_lengths=[2,1,3])
            False
            sage: [[5,4,3,1],[3,3,1]] in SkewPartitions(row_lengths=[2,1,1,2])
            False
            sage: [[5,4,3,1],[3,3,1]] in SkewPartitions(row_lengths=[2,1,2,1])
            True
        """
        if x in SkewPartitions():
            o = x[0]
            i = x[1]+[0]*(len(x[0])-len(x[1]))
            return [u[0]-u[1] for u in zip(o,i)] == self.co
        return False

    def _repr_(self):
        """
        TESTS::

            sage: SkewPartitions(row_lengths=[2,1])
            Skew partitions with row lengths [2, 1]
        """
        return "Skew partitions with row lengths %s" % self.co

    def _from_row_lengths_aux(self, sskp, ck_1, ck, overlap=0):
        """
        EXAMPLES::

            sage: s = SkewPartitions(row_lengths=[2,1])
            sage: list(s._from_row_lengths_aux([[1], []], 1, 1, overlap=0))
            [[1, 1] / [], [2, 1] / [1]]
            sage: list(s._from_row_lengths_aux([[1, 1], []], 1, 1, overlap=0))
            [[1, 1, 1] / [], [2, 2, 1] / [1, 1]]
            sage: list(s._from_row_lengths_aux([[2, 1], [1]], 1, 1, overlap=0))
            [[2, 1, 1] / [1], [3, 2, 1] / [2, 1]]
            sage: list(s._from_row_lengths_aux([[1], []], 1, 2, overlap=0))
            [[2, 2] / [1], [3, 2] / [2]]
            sage: list(s._from_row_lengths_aux([[2], []], 2, 1, overlap=0))
            [[2, 1] / [], [3, 1] / [1]]
        """
        nn = min(ck_1, ck)
        mm = max(0, ck-ck_1)
        # nn should be >= 0.  In the case of the positive overlap,
        # the min_part condition insures ck>=overlap for all k

        nn -= overlap
        for i in range(nn+1):
            skp1, skp2 = sskp
            skp2 += [0]*(len(skp1)-len(skp2))
            skp1 = [x + i + mm for x in skp1]
            skp1 += [ck]
            skp2 = [x + i + mm for x in skp2]
            skp2 = [x for x in skp2 if x != 0]
            yield SkewPartition([skp1, skp2])

    def __iter__(self):
        """
        Iterate through all the skew partitions that have row lengths
        given by the composition ``self.co``.

        EXAMPLES::

            sage: SkewPartitions(row_lengths=[2,2]).list()
            [[2, 2] / [], [3, 2] / [1], [4, 2] / [2]]
            sage: SkewPartitions(row_lengths=[2,2], overlap=1).list()
            [[2, 2] / [], [3, 2] / [1]]
        """
        if self.co == []:
            yield self.element_class(self, [[],[]])
            return

        nn = len(self.co)
        if nn == 1:
            yield self.element_class(self, [[self.co[0]],[]])
            return

        for sskp in SkewPartitions(row_lengths=self.co[:-1], overlap=self.overlap):
            for sp in self._from_row_lengths_aux(sskp, self.co[-2], self.co[-1], self.overlap):
                yield self.element_class(self, sp)


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.skew_partition', 'SkewPartition_class', SkewPartition)
