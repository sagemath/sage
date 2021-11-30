r"""
Fully packed loops

AUTHORS:

- Vincent Knight, James Campbell, Kevin Dilks, Emily Gunawan (2015): Initial version
- Vincent Delecroix (2017): cleaning and enhanced plotting function
"""
# ****************************************************************************
#       Copyright (C) 2015 Vincent Knight <vincent.knight@gmail.com>
#                          James Campbell <james.campbell@tanti.org.uk>
#                          Kevin Dilks <kdilks@gmail.com>
#                          Emily Gunawan <egunawan@umn.edu>
#                     2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import parent, Element

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.six_vertex_model import (SquareIceModel,
                                            SixVertexConfiguration,
                                            SixVertexModel)
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix

from sage.misc.decorators import options
from sage.matrix.constructor import matrix
from sage.arith.all import factorial
from sage.rings.integer import Integer
from sage.misc.misc_c import prod

# edges of a fpl in terms of the six vertex possible configurations
R = (1, 0)
L = (-1, 0)
U = (0, 1)
D = (0, -1)

FPL_edges = (
#   0 UD   1 RD,  2 UR,  3 LR,  4 LD   5 LU
   ((D,U), (L,D), (D,R), (R,L), (L,U), (R,U)),  # even
   ((R,L), (R,U), (L,U), (D,U), (D,R), (L,D))   # odd
   )

FPL_turns = (
# 0 UD          1 RD          2 UR          3 LR          4 LD          5 LU
 ({U: U, D: D}, {R: D, U: L}, {U: R, L: D}, {L: L, R: R}, {R: U, D: L}, {L: U, D: R}), # even
 ({L: L, R: R}, {L: U, D: R}, {R: U, D: L}, {U: U, D: D}, {U: R, L: D}, {R: D, U: L})  # odd
 )

def _make_color_list(n, colors=None,  color_map=None, randomize=False):
    r"""
    TESTS::

        sage: from sage.combinat.fully_packed_loop import _make_color_list
        sage: _make_color_list(5)
        sage: _make_color_list(5, ['blue', 'red'])
        ['blue', 'red', 'blue', 'red', 'blue']
        sage: _make_color_list(5, color_map='summer')
        [(0.0, 0.5, 0.4),
         (0.25098039215686274, 0.6254901960784314, 0.4),
         (0.5019607843137255, 0.7509803921568627, 0.4),
         (0.7529411764705882, 0.8764705882352941, 0.4),
         (1.0, 1.0, 0.4)]
        sage: l = _make_color_list(8, ['blue', 'red'], randomize=True)
        sage: len(l)
        8
        sage: l.count('blue')
        4
        sage: l.count('red')
        4
    """
    if colors:
        dim = len(colors)
        colors = [colors[i % dim] for i in range(n)]

    elif color_map:
        from matplotlib import cm
        if color_map not in cm.datad:
            raise ValueError('unknown color map %s' % color_map)
        cmap = cm.__dict__[color_map]
        colors = [cmap(i / float(n - 1))[:3] for i in range(n)]

    if colors and randomize:
        from sage.misc.prandom import shuffle
        shuffle(colors)

    return colors


class FullyPackedLoop(Element, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A class for fully packed loops.

    A fully packed loop is a collection of non-intersecting lattice paths on a square
    grid such that every vertex is part of some path, and the paths are either closed
    internal loops or have endpoints corresponding to alternate points on the
    boundary [Pro2001]_. They are known to be in bijection with alternating sign
    matrices.

    .. SEEALSO::

        :class:`AlternatingSignMatrix`

    To each fully packed loop, we assign a link pattern, which is the non-crossing
    matching attained by seeing which points on the boundary are connected
    by open paths in the fully packed loop.

    We can create a fully packed loop using the corresponding alternating sign
    matrix and also extract the link pattern::

        sage: A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        sage: fpl = FullyPackedLoop(A)
        sage: fpl.link_pattern()
        [(1, 4), (2, 3), (5, 6)]
        sage: fpl
                |         |
                |         |
                + -- +    +
                     |    |
                     |    |
             -- +    +    + --
                |    |
                |    |
                +    + -- +
                |         |
                |         |
        sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sage: fplb = FullyPackedLoop(B)
        sage: fplb.link_pattern()
        [(1, 6), (2, 5), (3, 4)]
        sage: fplb
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |

    The class also has a plot method::

        sage: fpl.plot()
        Graphics object consisting of 3 graphics primitives

    which gives:

    .. PLOT::
        :width: 200 px

        A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        fpl = FullyPackedLoop(A)
        p = fpl.plot()
        sphinx_plot(p)

    Note that we can also create a fully packed loop from a six vertex model configuration::

        sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
        sage: S
            ^    ^    ^
            |    |    |
        --> # -> # -> # <--
            ^    ^    |
            |    |    V
        --> # -> # <- # <--
            ^    |    |
            |    V    V
        --> # <- # <- # <--
            |    |    |
            V    V    V
        sage: fpl = FullyPackedLoop(S)
        sage: fpl
            |         |
            |         |
            + -- +    +
                 |    |
                 |    |
         -- +    +    + --
            |    |
            |    |
            +    + -- +
            |         |
            |         |

    Once we have a fully packed loop we can obtain the corresponding alternating sign matrix::

        sage: fpl.to_alternating_sign_matrix()
        [0 0 1]
        [0 1 0]
        [1 0 0]

    Here are some more examples using bigger ASMs::

        sage: A = AlternatingSignMatrix([[0,1,0,0],[0,0,1,0],[1,-1,0,1],[0,1,0,0]])
        sage: S = SixVertexModel(4, boundary_conditions='ice').from_alternating_sign_matrix(A)
        sage: fpl = FullyPackedLoop(S)
        sage: fpl.link_pattern()
        [(1, 2), (3, 6), (4, 5), (7, 8)]
        sage: fpl
            |         |
            |         |
            + -- + -- +    + --
                           |
                           |
         -- +    + -- + -- +
            |    |
            |    |
            +    +    + -- + --
            |    |    |
            |    |    |
         -- +    +    + -- +
                 |         |
                 |         |

        sage: m = AlternatingSignMatrix([[0,0,1,0,0,0],
        ....:                            [1,0,-1,0,1,0],
        ....:                            [0,0,0,1,0,0],
        ....:                            [0,1,0,0,-1,1],
        ....:                            [0,0,0,0,1,0],
        ....:                            [0,0,1,0,0,0]])
        sage: fpl = FullyPackedLoop(m)
        sage: fpl.link_pattern()
        [(1, 12), (2, 7), (3, 4), (5, 6), (8, 9), (10, 11)]
        sage: fpl
            |         |         |
            |         |         |
            + -- +    +    + -- +    + --
                 |    |    |         |
                 |    |    |         |
         -- + -- +    +    + -- + -- +
                      |
                      |
            + -- +    + -- + -- +    + --
            |    |              |    |
            |    |              |    |
         -- +    +    + -- +    +    +
                 |    |    |    |    |
                 |    |    |    |    |
            + -- +    + -- +    +    + --
            |                   |
            |                   |
         -- +    + -- + -- +    + -- +
                 |         |         |
                 |         |         |

        sage: m = AlternatingSignMatrix([[0,1,0,0,0,0,0],
        ....:                            [1,-1,0,0,1,0,0],
        ....:                            [0,0,0,1,0,0,0],
        ....:                            [0,1,0,0,-1,1,0],
        ....:                            [0,0,0,0,1,0,0],
        ....:                            [0,0,1,0,-1,0,1],
        ....:                            [0,0,0,0,1,0,0]])
        sage: fpl = FullyPackedLoop(m)
        sage: fpl.link_pattern()
        [(1, 2), (3, 4), (5, 6), (7, 8), (9, 14), (10, 11), (12, 13)]
        sage: fpl
            |         |         |         |
            |         |         |         |
            + -- + -- +    + -- +    + -- +
                           |         |
                           |         |
         -- + -- + -- +    + -- + -- +    + --
                      |                   |
                      |                   |
            + -- +    + -- + -- +    + -- +
            |    |              |    |
            |    |              |    |
         -- +    +    + -- +    +    +    + --
                 |    |    |    |    |    |
                 |    |    |    |    |    |
            + -- +    + -- +    +    + -- +
            |                   |
            |                   |
         -- +    + -- + -- +    +    + -- + --
                 |         |    |    |
                 |         |    |    |
            + -- +    + -- +    +    + -- +
            |         |         |         |
            |         |         |         |

    Gyration on an alternating sign matrix/fully packed loop ``fpl``
    of the link pattern corresponding to ``fpl``::

        sage: ASMs = AlternatingSignMatrices(3).list()
        sage: ncp = FullyPackedLoop(ASMs[1]).link_pattern() # fpl's gyration orbit size is 2
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(5):
        ....:         a,b=a%6+1,b%6+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(ASMs[1].gyration().to_fully_packed_loop().link_pattern()) ==\
        ....:     PerfectMatching(rotated_ncp)
        True

        sage: fpl = FullyPackedLoop(ASMs[0])
        sage: ncp = fpl.link_pattern() # fpl's gyration size is 3
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(5):
        ....:         a,b=a%6+1,b%6+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(ASMs[0].gyration().to_fully_packed_loop().link_pattern()) ==\
        ....:     PerfectMatching(rotated_ncp)
        True

        sage: mat = AlternatingSignMatrix([[0,0,1,0,0,0,0],[1,0,-1,0,1,0,0],
        ....:     [0,0,1,0,0,0,0],[0,1,-1,0,0,1,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,0,1]])
        sage: fpl = FullyPackedLoop(mat) # n=7
        sage: ncp = fpl.link_pattern()
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(13):
        ....:         a,b=a%14+1,b%14+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
        ....:     PerfectMatching(rotated_ncp)
        True

        sage: mat = AlternatingSignMatrix([[0,0,0,1,0,0], [0,0,1,-1,1,0],
        ....:     [0,1,0,0,-1,1], [1,0,-1,1,0,0], [0,0,1,0,0,0], [0,0,0,0,1,0]])
        sage: fpl = FullyPackedLoop(mat) # n =6
        sage: ncp = fpl.link_pattern()
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(11):
        ....:         a,b=a%12+1,b%12+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
        ....:     PerfectMatching(rotated_ncp)
        True

    More examples:

    We can initiate a fully packed loop using an alternating sign matrix::

        sage: A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        sage: fpl = FullyPackedLoop(A)
        sage: fpl
            |         |
            |         |
            + -- +    +
                 |    |
                 |    |
         -- +    +    + --
            |    |
            |    |
            +    + -- +
            |         |
            |         |
        sage: FullyPackedLoops(3)(A) == fpl
        True

    We can also input a matrix::

        sage: FullyPackedLoop([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
            |         |
            |         |
            + -- +    +
                 |    |
                 |    |
         -- +    +    + --
            |    |
            |    |
            +    + -- +
            |         |
            |         |
        sage: FullyPackedLoop([[0, 0, 1], [0, 1, 0], [1, 0, 0]]) ==\
        ....: FullyPackedLoops(3)([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        True

    Otherwise we initiate a fully packed loop using a six vertex model::

        sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
        sage: fpl = FullyPackedLoop(S)
        sage: fpl
            |         |
            |         |
            + -- +    +
                 |    |
                 |    |
         -- +    +    + --
            |    |
            |    |
            +    + -- +
            |         |
            |         |

        sage: FullyPackedLoops(3)(S) == FullyPackedLoop(S)
        True

        sage: fpl.six_vertex_model().to_alternating_sign_matrix()
        [0 0 1]
        [0 1 0]
        [1 0 0]

    We can also input the matrix associated to a six vertex model::

        sage: SixVertexModel(2)([[3,1],[5,3]])
            ^    ^
            |    |
        --> # <- # <--
            |    ^
            V    |
        --> # -> # <--
            |    |
            V    V

        sage: FullyPackedLoop([[3,1],[5,3]])
            |
            |
            +    + --
            |    |
            |    |
         -- +    +
                 |
                 |

        sage: FullyPackedLoops(2)([[3,1],[5,3]]) == FullyPackedLoop([[3,1],[5,3]])
        True

    Note that the matrix corresponding to a six vertex model without
    the ice boundary condition is not allowed::

        sage: SixVertexModel(2)([[3,1],[5,5]])
            ^    ^
            |    |
        --> # <- # <--
            |    ^
            V    V
        --> # -> # -->
            |    |
            V    V

        sage: FullyPackedLoop([[3,1],[5,5]])
        Traceback (most recent call last):
        ...
        ValueError: invalid alternating sign matrix

        sage: FullyPackedLoops(2)([[3,1],[5,5]])
        Traceback (most recent call last):
        ...
        ValueError: invalid alternating sign matrix

    Note that if anything else is used to generate the fully packed loop an error will occur::

        sage: fpl = FullyPackedLoop(5)
        Traceback (most recent call last):
        ...
        ValueError: invalid alternating sign matrix

        sage: fpl = FullyPackedLoop((1, 2, 3))
        Traceback (most recent call last):
        ...
        ValueError: The alternating sign matrices must be square

        sage: SVM = SixVertexModel(3)[0]
        sage: FullyPackedLoop(SVM)
        Traceback (most recent call last):
        ...
        ValueError: invalid alternating sign matrix

    REFERENCES:

    - [Pro2001]_
    - [Str2015]_
    """
    @staticmethod
    def __classcall_private__(cls, generator):
        """
        Create a FPL.

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: FullyPackedLoop(A)
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |

            sage: SVM = SixVertexModel(4, boundary_conditions='ice')[0]
            sage: FullyPackedLoop(SVM)
                |         |
                |         |
                +    + -- +    + --
                |    |         |
                |    |         |
             -- +    +    + -- +
                     |    |
                     |    |
                + -- +    +    + --
                |         |    |
                |         |    |
             -- +    + -- +    +
                     |         |
                     |         |
        """
        if isinstance(generator, AlternatingSignMatrix):
            SVM = generator.to_six_vertex_model()
        elif isinstance(generator, SquareIceModel.Element):
            SVM = generator
        elif isinstance(generator, SixVertexConfiguration):
            # Check that this is an ice square model
            generator = SixVertexModel(generator.parent()._nrows,
                                       boundary_conditions='ice')(generator)
            M = generator.to_alternating_sign_matrix().to_matrix()
            AlternatingSignMatrix(M)
            SVM = generator
        else: # Not ASM nor SVM
            try:
                SVM = AlternatingSignMatrix(generator).to_six_vertex_model()
            except (TypeError, ValueError):
                generator = matrix(generator)
                generator = SixVertexModel(generator.nrows(), boundary_conditions='ice')(generator)
                # Check that this is an ice square model
                generator.to_alternating_sign_matrix()
                SVM = generator

        if not SVM:
            raise TypeError('generator for FullyPackedLoop must either be an '
                            'AlternatingSignMatrix or a SquareIceModel.Element')
        FPLs = FullyPackedLoops(len(SVM))
        return FPLs(generator)

    def __init__(self, parent, generator):
        """
        Initialise object, can take ASM of FPL as generator.

        TESTS::

            sage: A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: TestSuite(fpl).run()

        """
        if isinstance(generator, AlternatingSignMatrix):
            self._six_vertex_model = generator.to_six_vertex_model()
        elif isinstance(generator, SquareIceModel.Element):
            self._six_vertex_model = generator

        Element.__init__(self, parent)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |

            sage: A = AlternatingSignMatrix([[0,1,0,0],[0,0,1,0],[1,-1,0,1],[0,1,0,0]])
            sage: S = SixVertexModel(4, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl
                |         |
                |         |
                + -- + -- +    + --
                               |
                               |
             -- +    + -- + -- +
                |    |
                |    |
                +    +    + -- + --
                |    |    |
                |    |    |
             -- +    +    + -- +
                     |         |
                     |         |

        """
        # List are in the order of URDL
        # One set of rules for how to draw around even vertex, one set of rules for odd vertex
        n = len(self._six_vertex_model) - 1
        ascii1 = [[r'     ', ' -', r'     ', '- '], # LR
                 [r'  |  ', '  ', r'     ', '- '], # LU
                 [r'     ', '  ', r'  |  ', '- '], # LD
                 [r'  |  ', '  ', r'  |  ', '  '], # UD
                 [r'  |  ', ' -', r'     ', '  '], # UR
                 [r'     ', ' -', r'  |  ', '  ']] # RD

        ascii2 = [[r'  |  ', '  ', r'  |  ', '  '], # LR
                 [r'     ', ' -', r'  |  ', '  '], # LU
                 [r'  |  ', ' -', r'     ', '  '], # LD
                 [r'     ', ' -', r'     ', '- '], # UD
                 [r'     ', '  ', r'  |  ', '- '], # UR
                 [r'  |  ', '  ', r'     ', '- ']] # RD
        ret = '  '
        # Do the top line
        for i,entry in enumerate(self._six_vertex_model[0]):
            if i % 2 == 0:
                ret += '  |  '
            else:
                ret += '     '

        plus_sign = '+'

        # Do the meat of the ascii art
        for j,row in enumerate(self._six_vertex_model):
            ret += '\n  '
            # Do the top row
            for i,entry in enumerate(row):
                if (i + j) % 2 == 0:
                    ret += ascii1[entry][0]
                else:
                    ret += ascii2[entry][0]
            ret += '\n'

            # Do the left-most entry
            if (j) % 2 == 0:
                ret += '  '
            else:
                ret += ' -'

            # Do the middle row
            for i,entry in enumerate(row):
                if (i + j) % 2 == 0:
                    ret += ascii1[entry][3] + plus_sign + ascii1[entry][1]
                else:
                    ret += ascii2[entry][3] + plus_sign + ascii2[entry][1]

            # Do the right-most entry
            if (j+n) % 2 ==0:
                ret += '  '
            else:
                ret += '- '

            # Do the bottom row
            ret += '\n  '
            for i,entry in enumerate(row):
                if (i + j) % 2 ==0:
                    ret += ascii1[entry][2]
                else:
                    ret += ascii2[entry][2]

        # Do the bottom line
        ret += '\n  '
        for i,entry in enumerate(self._six_vertex_model[-1]):
            if (i+n+1) % 2 ==0:
                ret += '     '
            else:
                ret += '  |  '

        return ret

    def _richcmp_(self, other, op):
        """
        Check equality or inequality.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A.random_element()
            sage: FullyPackedLoop(M) == M.to_fully_packed_loop()
            True

            sage: FullyPackedLoop(A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])) ==\
            ....:    FullyPackedLoop(A([[1, 0, 0],[0, 0, 1],[0, 1, 0]]))
            False

            sage: FullyPackedLoop(M) == M
            False

            sage: M = A.random_element()
            sage: FullyPackedLoop(M) != M.to_fully_packed_loop()
            False

            sage: f0 = FullyPackedLoop(A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]))
            sage: f1 = FullyPackedLoop(A([[1, 0, 0],[0, 0, 1],[0, 1, 0]]))
            sage: f0 != f1
            True
        """
        return self._six_vertex_model._richcmp_(other._six_vertex_model, op)

    def to_alternating_sign_matrix(self):
        """
        Return the alternating sign matrix corresponding to this class.

        .. SEEALSO::

            :class:`AlternatingSignMatrix`

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl.to_alternating_sign_matrix()
            [ 0  1  0]
            [ 1 -1  1]
            [ 0  1  0]
            sage: A = AlternatingSignMatrix([[0,1,0,0],[0,0,1,0],[1,-1,0,1],[0,1,0,0]])
            sage: S = SixVertexModel(4, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl.to_alternating_sign_matrix()
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 1 -1  0  1]
            [ 0  1  0  0]
        """
        return self._six_vertex_model.to_alternating_sign_matrix()


    @options(link=True, loop=True, loop_fill=False)
    def plot(self, **options):
        r"""
        Return a graphical object of the Fully Packed Loop.

        Each option can be specified separately for links (the curves that join
        boundary points) and the loops. In order to do so, you need to prefix
        its name with either ``'link_'`` or ``'loop_'``. As an example, setting
        ``color='red'`` will color both links and loops in red while setting
        ``link_color='red'`` will only apply the color option for the links.

        INPUT:

        - ``link``, ``loop`` - (boolean, default ``True``) whether to plot the links
          or the loops

        - ``color``, ``link_color``, ``loop_color`` - (optional, a string or a
          RGB triple)

        - ``colors``, ``link_colors``, ``loop_colors`` - (optional, list) a list of
          colors

        - ``color_map``, ``link_color_map``, ``loop_color_map`` - (string,
          optional) a name of a matplotlib color map for the link or the loop

        - ``link_color_randomize`` - (boolean, default ``False``) when
          ``link_colors`` or ``link_color_map`` is specified it randomizes
          its order. Setting this option to ``True`` makes it unlikely to
          have two neighboring links with the same color.

        - ``loop_fill`` - (boolean, optional) whether to fill the interior of the loops

        EXAMPLES:

        To plot the fully packed loop associated to the following alternating sign
        matrix

        .. MATH::

            \begin{pmatrix} 0&1&1 \\ 1&-1&1 \\ 0&1&0 \end{pmatrix}

        simply do::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl.plot()
            Graphics object consisting of 3 graphics primitives

        The resulting graphics is as follows

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        You can also have the three links in different colors with::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl.plot(link_color_map='rainbow')
            Graphics object consisting of 3 graphics primitives

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot(link_color_map='rainbow')
            sphinx_plot(p)

        You can plot the 42 fully packed loops of size `4 \times 4` using::

            sage: G = [fpl.plot(link_color_map='winter', loop_color='black') for fpl in FullyPackedLoops(4)]
            sage: graphics_array(G, 7, 6)
            Graphics Array of size 7 x 6

        .. PLOT::
            :width: 600 px

            G = [fpl.plot(link_color_map='winter', loop_color='black') for fpl in FullyPackedLoops(4)]
            p = graphics_array(G, 7, 6)
            sphinx_plot(p)

        Here is an example of a `20 \times 20` fully packed loop::

            sage: s = "00000000000+0000000000000000+00-0+00000000000+00-00+0-+00000\
            ....: 0000+-00+00-+00000000+00-0000+0000-+00000000+000-0+0-+0-+000\
            ....: 000+-000+-00+0000000+-+-000+00-+0-000+000+-000+-0+0000000-0+\
            ....: 0000+0-+0-+00000-+00000+-+0-0+-00+0000000000+-0000+0-00+0000\
            ....: 000000+0-000+000000000000000+0000-00+00000000+0000-000+00000\
            ....: 00+0-00+0000000000000000+-0000+000000-+000000+00-0000+-00+00\
            ....: 00000000+-0000+00000000000000+0000000000"
            sage: a = matrix(20, [{'0':0, '+':1, '-': -1}[i] for i in s])
            sage: fpl = FullyPackedLoop(a)
            sage: fpl.plot(loop_fill=True, loop_color_map='rainbow')
            Graphics object consisting of 27 graphics primitives

        .. PLOT::
            :width: 400 px

            s = "00000000000+0000000000000000+00-0+00000000000+00-00+0-+00000\
            0000+-00+00-+00000000+00-0000+0000-+00000000+000-0+0-+0-+000\
            000+-000+-00+0000000+-+-000+00-+0-000+000+-000+-0+0000000-0+\
            0000+0-+0-+00000-+00000+-+0-0+-00+0000000000+-0000+0-00+0000\
            000000+0-000+000000000000000+0000-00+00000000+0000-000+00000\
            00+0-00+0000000000000000+-0000+000000-+000000+00-0000+-00+00\
            00000000+-0000+00000000000000+0000000000"
            a = matrix(20, [{'0':0, '+':1, '-': -1}[i] for i in s])
            p = FullyPackedLoop(a).plot(loop_fill=True, loop_color_map='rainbow')
            sphinx_plot(p)
        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line2d
        from sage.plot.polygon import polygon2d

        link_options = {}
        loop_options = {}
        for k, v in options.items():
            if k == 'link':
                link = v
            elif k == 'loop':
                loop = v
            elif k.startswith('link_'):
                link_options[k[5:]] = v
            elif k.startswith('loop_'):
                loop_options[k[5:]] = v
            else:
                link_options[k] = v
                loop_options[k] = v

        sv = self._six_vertex_model
        n = len(sv)

        # LR boundaries => odd sum
        # UD boundaries => even sum
        rank = self.parent()._boundary_index
        unrank = self.parent()._boundary
        seen = [False] * (2*n)

        squares = set((i,j) for i in range(n) for j in range(n))

        colors = _make_color_list(2*n,
                colors = link_options.pop('colors', None),
                color_map = link_options.pop('color_map', None),
                randomize = link_options.pop('color_randomize', False))

        G = Graphics()
        for i in range(2*n):
            if seen[i]:
                continue
            orbit = self._link_or_loop_from(unrank(i))
            j = rank(orbit[-1])
            seen[i] = seen[j] = True
            squares.difference_update(orbit)

            if link:
                if colors:
                    link_options['color'] = colors.pop()

                # make it upside down
                orbit = [(j, n - i - 1) for i, j in orbit]
                G += line2d(orbit, **link_options)

        loops = []
        while squares:
            orbit = self._link_or_loop_from(squares.pop())
            loops.append(orbit)
            squares.difference_update(orbit)

        if loop:
            colors = _make_color_list(len(loops),
                    colors = loop_options.pop('colors', None),
                    color_map = loop_options.pop('color_map', None),
                    randomize = loop_options.pop('color_randomize', False))

            fill = loop_options.pop('fill')

            for orbit in loops:
                if colors:
                    loop_options['color'] = colors.pop()

                # make it upside down
                orbit = [(j, n - i - 1) for i,j in orbit]

                if fill:
                    G += polygon2d(orbit, **loop_options)
                else:
                    G += line2d(orbit, **loop_options)

        G.axes(False)
        G.set_aspect_ratio(1)
        return G

    def gyration(self):
        r"""
        Return the fully packed loop obtained by applying gyration
        to the alternating sign matrix in bijection with ``self``.

        Gyration was first defined in [Wie2000]_ as an action on
        fully-packed loops.

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl.gyration().to_alternating_sign_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: asm = AlternatingSignMatrix([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: f = FullyPackedLoop(asm)
            sage: f.gyration().to_alternating_sign_matrix()
            [0 1 0]
            [0 0 1]
            [1 0 0]
        """
        return FullyPackedLoop(self.to_alternating_sign_matrix().gyration())

    def _link_or_loop_from(self, pos, d0=None):
        r"""
        Return the coordinates of the line passing through ``pos``.

        EXAMPLES:

        A link::

            sage: fpl = FullyPackedLoops(4).first()
            sage: fpl._link_or_loop_from((2,2))
            [(0, 4), (0, 3), (1, 3), (1, 2), (2, 2), (3, 2), (3, 1), (4, 1)]
            sage: fpl._link_or_loop_from((-1, 0))
            [(-1, 0), (0, 0), (1, 0), (1, -1)]

        A loop::

            sage: a = AlternatingSignMatrix([[0,1,0,0], [0,0,0,1], [1,0,0,0], [0,0,1,0]])
            sage: fpl = FullyPackedLoop(a)
            sage: fpl._link_or_loop_from((1,1))
            [(1, 1), (2, 1), (2, 2), (1, 2), (1, 1)]
        """
        global R, L, U, D, FPL_turns, FPL_edges

        orbit = [pos]
        sv = self._six_vertex_model
        n = len(sv)
        i, j = pos

        # deal with boundary cases
        if i < -1 or i > n or j < -1 or j > n:
            raise ValueError('indices out of range')
        if (i == -1 or i == n) and not (i + j) % 2:
            raise ValueError('left and right boundary values must have odd sum')
        if (j == -1 or j == n) and (i + j) % 2:
            raise ValueError('up and down boundary values must have even sum')

        if i == -1:
            d = R
        elif i == n:
            d = L
        elif j == -1:
            d = U
        elif j == n:
            d = D
        elif d0 is None:
            d = FPL_edges[(i + j) % 2][sv[i][j]][0]
        elif d0 in FPL_edges[(i + j) % 2][sv[i][j]]:
            d = d0
        else:
            raise ValueError('invalid direction')

        # compute the link or loop
        while True:
            i += d[0]
            j += d[1]
            orbit.append((i, j))
            if (i, j) == orbit[0] or i == -1 or j == -1 or i == n or j == n:
                break

            conf = sv[i][j]
            parity = (i + j) % 2

            d = FPL_turns[parity][conf][d]
            if d is None:
                raise RuntimeError

        if i == -1 or j == -1 or i == n or j == n:
            i0,j0 = orbit[0]
            if d0 is None and i0 != -1 and i0 != n and j0 != -1 and j0 != n:
                # only half of a link -> compute the other half
                i1,j1 = orbit[1]
                d = (i0-i1, j0-j1)
                orbit2 = self._link_or_loop_from(orbit[1], d)
                assert orbit2[0] == (i1,j1) and orbit2[1] == (i0,j0)
                return orbit2[:1:-1] + orbit
            return orbit
        else:
            return orbit

    def link_pattern(self):
        r"""
        Return a link pattern corresponding to a fully packed loop.

        Here we define a link pattern `LP` to be a partition of the list
        `[1, ..., 2k]` into 2-element sets (such a partition is also known as
        a perfect matching) such that the following non-crossing condition holds:
        Let the numbers `1, ..., 2k` be written on the perimeter of a circle.
        For every 2-element set `(a,b)` of the partition `LP`, draw an arc
        linking the two numbers `a` and `b`. We say that `LP` is non-crossing
        if every arc can be drawn so that no two arcs intersect.

        Since every endpoint of a fully packed loop `fpl` is connected to a different
        endpoint, there is a natural surjection from the fully packed loops on an
        nxn grid onto the link patterns on the list `[1, \dots, 2n]`.
        The pairs of connected endpoints of a fully packed loop `fpl` correspond to
        the 2-element tuples of the corresponding link pattern.

        .. SEEALSO::

            :class:`PerfectMatching`

        .. NOTE::

            by convention, we choose the top left vertex to be even.
            See [Pro2001]_ and [Str2015]_.

        EXAMPLES:

        We can extract the underlying link pattern (a non-crossing
        partition) from a fully packed loop::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl.link_pattern()
            [(1, 2), (3, 6), (4, 5)]

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl.link_pattern()
            [(1, 6), (2, 5), (3, 4)]

        Gyration on an alternating sign matrix/fully packed loop ``fpl``
        corresponds to a rotation (i.e. a becomes a-1 mod 2n)
        of the link pattern corresponding to ``fpl``::

            sage: ASMs = AlternatingSignMatrices(3).list()
            sage: ncp = FullyPackedLoop(ASMs[1]).link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(5):
            ....:         a,b=a%6+1,b%6+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(ASMs[1].gyration().to_fully_packed_loop().link_pattern()) ==\
            ....:     PerfectMatching(rotated_ncp)
            True

            sage: fpl = FullyPackedLoop(ASMs[0])
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(5):
            ....:         a,b=a%6+1,b%6+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(ASMs[0].gyration().to_fully_packed_loop().link_pattern()) ==\
            ....:     PerfectMatching(rotated_ncp)
            True

            sage: mat = AlternatingSignMatrix([[0,0,1,0,0,0,0],[1,0,-1,0,1,0,0],[0,0,1,0,0,0,0],
            ....:     [0,1,-1,0,0,1,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,0,1]])
            sage: fpl = FullyPackedLoop(mat) # n=7
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(13):
            ....:         a,b=a%14+1,b%14+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
            ....:     PerfectMatching(rotated_ncp)
            True

            sage: mat = AlternatingSignMatrix([[0,0,0,1,0,0], [0,0,1,-1,1,0], [0,1,0,0,-1,1], [1,0,-1,1,0,0],
            ....:     [0,0,1,0,0,0], [0,0,0,0,1,0]])
            sage: fpl = FullyPackedLoop(mat)
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(11):
            ....:         a,b=a%12+1,b%12+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
            ....:     PerfectMatching(rotated_ncp)
            True

        TESTS:

        We test previous two bugs which showed up when this method is called twice::

            sage: A = AlternatingSignMatrices(6)
            sage: B = A.random_element()
            sage: C = FullyPackedLoop(B)
            sage: D = C.link_pattern()
            sage: E = C.link_pattern()
            sage: D == E
            True
        """
        global L, R, U, D, FPL_turns

        link_pattern = []
        n = len(self._six_vertex_model)
        seen = [False] * (2*n)
        unrank = self.parent()._boundary
        rank = self.parent()._boundary_index
        sv = self._six_vertex_model

        for k in range(2*n):
            if seen[k]:
                continue

            i,j = unrank(k)

            # initial direction
            if i == -1:
                d = R
            elif i == n:
                d = L
            elif j == -1:
                d = U
            elif j == n:
                d = D

            # go through the link
            while True:
                i += d[0]
                j += d[1]
                if i == -1 or j == -1 or i == n or j == n:
                    break

                conf = sv[i][j]
                parity = (i + j) % 2
                d = FPL_turns[parity][conf][d]

            # update seen and link_pattern
            l = rank((i,j))
            seen[k] = seen[l] = True
            link_pattern.append((k+1, l+1))

        return link_pattern

    def six_vertex_model(self):
        """
        Return the underlying six vertex model configuration.

        EXAMPLES::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |
            sage: fpl.six_vertex_model()
                ^    ^    ^
                |    |    |
            --> # <- # <- # <--
                |    ^    ^
                V    |    |
            --> # -> # <- # <--
                |    |    ^
                V    V    |
            --> # -> # -> # <--
                |    |    |
                V    V    V
        """
        return self._six_vertex_model

class FullyPackedLoops(Parent, UniqueRepresentation):
    r"""
    Class of all fully packed loops on an  `n \times n` grid.

    They are known to be in bijection with alternating sign matrices.

    .. SEEALSO::

        :class:`AlternatingSignMatrices`

    INPUT:

    - ``n`` -- the number of row (and column) or grid

    EXAMPLES:

    This will create an instance to manipulate the fully packed loops of size 3::

        sage: FPLs = FullyPackedLoops(3)
        sage: FPLs
        Fully packed loops on a 3x3 grid
        sage: FPLs.cardinality()
        7

    When using the square ice model, it is known that the number of
    configurations is equal to the number of alternating sign matrices::

        sage: M = FullyPackedLoops(1)
        sage: len(M)
        1
        sage: M = FullyPackedLoops(4)
        sage: len(M)
        42
        sage: all(len(SixVertexModel(n, boundary_conditions='ice'))
        ....:     == FullyPackedLoops(n).cardinality() for n in range(1, 7))
        True
    """
    def __init__(self, n):
        r"""
        Initialize ``self``.

        TESTS::

            sage: FPLs = FullyPackedLoops(3)
            sage: TestSuite(FPLs).run()
        """
        self._n = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: FPLs = FullyPackedLoops(2)
            sage: len(FPLs)
            2
        """
        for X in SixVertexModel(self._n, boundary_conditions='ice'):
            yield self.element_class(self, X)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: FPLs = FullyPackedLoops(4); FPLs
            Fully packed loops on a 4x4 grid
        """
        return "Fully packed loops on a %sx%s grid" % (self._n,self._n)

    def __contains__(self, fpl):
        """
        Check if ``fpl`` is in ``self``.

        TESTS::

            sage: FPLs = FullyPackedLoops(3)
            sage: FullyPackedLoop(AlternatingSignMatrix([[0,1,0],[1,0,0],[0,0,1]])) in FPLs
            True
            sage: FullyPackedLoop(AlternatingSignMatrix([[0,1,0],[1,-1,1],[0,1,0]])) in FPLs
            True
            sage: FullyPackedLoop(AlternatingSignMatrix([[0, 1],[1,0]])) in FPLs
            False
            sage: FullyPackedLoop(AlternatingSignMatrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])) in FPLs
            False
            sage: [1,2,3] in FPLs
            False
        """
        return parent(fpl) is self

    def _element_constructor_(self, generator):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: FPLs = FullyPackedLoops(4)
            sage: M = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
            sage: A = AlternatingSignMatrix(M)
            sage: elt = FullyPackedLoop(A)
            sage: FPL = FPLs(elt); FPL
                |         |
                |         |
                +    + -- +    + --
                |    |         |
                |    |         |
             -- +    +    + -- +
                     |    |
                     |    |
                + -- +    +    + --
                |         |    |
                |         |    |
             -- +    + -- +    +
                     |         |
                     |         |

            sage: FPLs(A) == FPL
            True

            sage: FPLs(M) == FPL
            True

            sage: FPLs(FPL._six_vertex_model) == FPL
            True

            sage: FPL.parent() is FPLs
            True

            sage: FPL = FullyPackedLoops(2)
            sage: FPL([[3,1],[5,3]])
                |
                |
                +    + --
                |    |
                |    |
             -- +    +
                     |
                     |
        """
        if isinstance(generator, AlternatingSignMatrix):
            SVM = generator.to_six_vertex_model()
        elif isinstance(generator, SquareIceModel.Element) or \
        isinstance(generator, SixVertexConfiguration):
            SVM = generator
        else:  # Not ASM nor SVM
            try:
                SVM = AlternatingSignMatrix(generator).to_six_vertex_model()
            except (TypeError, ValueError):
                SVM = SixVertexModel(self._n, boundary_conditions='ice')(generator)
                SVM.to_alternating_sign_matrix()
        if len(SVM) != self._n:
            raise ValueError("invalid size")
        return self.element_class(self, SVM)

    Element = FullyPackedLoop

    def size(self):
        r"""
        Return the size of the matrices in ``self``.

        TESTS::

            sage: FPLs = FullyPackedLoops(4)
            sage: FPLs.size()
            4
        """
        return self._n

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of fully packed loops on  `n \times n` grid

        .. MATH::

            \prod_{k=0}^{n-1} \frac{(3k+1)!}{(n+k)!} = \frac{1! 4! 7! 10!
            \cdots (3n-2)!}{n! (n+1)! (n+2)! (n+3)! \cdots (2n-1)!}.

        EXAMPLES::

            sage: [AlternatingSignMatrices(n).cardinality() for n in range(11)]
            [1, 1, 2, 7, 42, 429, 7436, 218348, 10850216, 911835460, 129534272700]
        """
        return Integer(prod( [ factorial(3*k+1)/factorial(self._n+k)
                       for k in range(self._n)] ))

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: FPLs = FullyPackedLoops(3)
            sage: FPLs.an_element()
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |
        """
        #ASM = AlternatingSignMatrix(matrix.identity(self._n))
        #SVM = ASM.to_six_vertex_model()
        SVM = SixVertexModel(self._n,boundary_conditions='ice').an_element()
        return self.element_class(self, SVM)

    def _boundary(self, k):
        r"""
        Return the coordinates of the ``k``-th boundary.

        TESTS::

            sage: F = FullyPackedLoops(5)
            sage: [F._boundary(k) for k in range(10)] == F._boundaries()
            True
            sage: all(F._boundary_index(F._boundary(k)) == k for k in range(10))
            True

            sage: F = FullyPackedLoops(6)
            sage: [F._boundary(k) for k in range(12)] == F._boundaries()
            True
            sage: all(F._boundary_index(F._boundary(k)) == k for k in range(12))
            True
        """
        n = self._n
        n_LR = n//2 if n%2 == 0 else (n+1) // 2
        n_TB = n//2 if n%2 == 0 else (n-1) // 2
        if k < n_LR:
            return (-1, 2*k)
        k -= n_LR
        if k < n_TB:
            return (n%2 + 2*k, n)
        k -= n_TB
        if k < n_LR:
            return (n, n - 1 - 2*k)
        k -= n_LR
        if k < n_TB:
            return (n - 1 - n%2 - 2*k, -1)

    def _boundary_index(self, pos):
        r"""
        Return the index of the boundary at position ``pos``.

        TESTS::

            sage: F = FullyPackedLoops(5)
            sage: [F._boundary_index(b) for b in F._boundaries()]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: all(F._boundary(F._boundary_index(b)) == b for b in F._boundaries())
            True

            sage: F = FullyPackedLoops(6)
            sage: [F._boundary_index(b) for b in F._boundaries()]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            sage: all(F._boundary(F._boundary_index(b)) == b for b in F._boundaries())
            True
        """
        n = self._n
        i, j = pos
        if i == -1:
            return j//2
        elif j == n:
            return (n + 1) // 2 + i // 2
        elif i == n:
            return n + (n - j) // 2
        elif j == -1:
            return 3 * n // 2 + (n - i) // 2

    def _boundaries(self):
        r"""
        Return the list of coordinates for the link in the boundaries.

        TESTS::

            sage: FullyPackedLoops(5)._boundaries()
            [(-1, 0), (-1, 2), (-1, 4), (1, 5), (3, 5),
             (5, 4),  (5, 2), (5, 0), (3, -1), (1, -1)]

            sage: FullyPackedLoops(6)._boundaries()
            [(-1, 0), (-1, 2), (-1, 4), (0, 6), (2, 6), (4, 6),
             (6, 5), (6, 3), (6, 1), (5, -1), (3, -1), (1, -1)]
        """
        n = self._n
        boundaries = []
        # left side: j = 0 mod 2
        boundaries.extend((-1, j) for j in range(0, n, 2))
        # top side: i = n mod 2
        boundaries.extend((i, n) for i in range(n % 2, n, 2))
        # right side: j = n+1 mod 2
        boundaries.extend((n, j) for j in range(n - 1, -1, -2))
        # bottom side: i = 1 mod 2
        boundaries.extend((i, -1) for i in range(n - 1 - n % 2, -1, -2))
        return boundaries
