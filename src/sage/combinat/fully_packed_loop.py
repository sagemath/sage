r"""
Fully packed loops
"""
# python3
from __future__ import division, print_function

from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.six_vertex_model import SquareIceModel, \
SixVertexConfiguration, SixVertexModel
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix
from sage.plot.graphics import Graphics
from sage.matrix.constructor import matrix
from sage.plot.line import line
from sage.combinat.perfect_matching import PerfectMatching
from sage.arith.all import factorial
from sage.rings.integer import Integer
from sage.misc.all import prod
from sage.misc.lazy_attribute import lazy_attribute


class FullyPackedLoop(Element):
    r"""
    A class for fully packed loops.

    A fully packed loop is a collection of non-intersecting lattice paths on a square
    grid such that every vertex is part of some path, and the paths are either closed
    internal loops or have endpoints corresponding to alternate points on the
    boundary [Propp2001]_. They are known to be in bijection with alternating sign
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
        Graphics object consisting of 15 graphics primitives

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
        ....:     for i in range(0,5):
        ....:         a,b=a%6+1,b%6+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(ASMs[1].gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

        sage: fpl = FullyPackedLoop(ASMs[0])
        sage: ncp = fpl.link_pattern() # fpl's gyration size is 3
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(0,5):
        ....:         a,b=a%6+1,b%6+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(ASMs[0].gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

        sage: mat = AlternatingSignMatrix([[0,0,1,0,0,0,0],[1,0,-1,0,1,0,0],\
        [0,0,1,0,0,0,0],[0,1,-1,0,0,1,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,0,1]])
        sage: fpl = FullyPackedLoop(mat) # n=7
        sage: ncp = fpl.link_pattern()
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(0,13):
        ....:         a,b=a%14+1,b%14+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

        sage: mat = AlternatingSignMatrix([[0,0,0,1,0,0], [0,0,1,-1,1,0],\
        [0,1,0,0,-1,1], [1,0,-1,1,0,0], [0,0,1,0,0,0], [0,0,0,0,1,0]])
        sage: fpl = FullyPackedLoop(mat) # n =6
        sage: ncp = fpl.link_pattern()
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(0,11):
        ....:         a,b=a%12+1,b%12+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

    More examples::

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
            ValueError: Invalid alternating sign matrix

            sage: FullyPackedLoops(2)([[3,1],[5,5]])
            Traceback (most recent call last):
            ...
            ValueError: Invalid alternating sign matrix

        Note that if anything else is used to generate the fully packed loop an error will occur::

            sage: fpl = FullyPackedLoop(5)
            Traceback (most recent call last):
            ...
            ValueError: Invalid alternating sign matrix

            sage: fpl = FullyPackedLoop((1, 2, 3))
            Traceback (most recent call last):
            ...
            ValueError: The alternating sign matrices must be square

            sage: SVM = SixVertexModel(3)[0]
            sage: FullyPackedLoop(SVM)
            Traceback (most recent call last):
            ...
            ValueError: Invalid alternating sign matrix

    REFERENCES:

    .. [Propp2001] James Propp.
       *The Many Faces of Alternating Sign Matrices*,
       Discrete Mathematics and Theoretical Computer Science 43 (2001): 58
       :arxiv:`math/0208125`

    .. [Striker2015] Jessica Striker.
       *The toggle group, homomesy, and the Razumov-Stroganov correspondence*,
       Electron. J. Combin. 22 (2015) no. 2
       :arxiv:`1503.08898`
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

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
            generator = SixVertexModel(generator.parent()._nrows, \
            boundary_conditions='ice')(generator)
            M = generator.to_alternating_sign_matrix().to_matrix()
            M = AlternatingSignMatrix(M)
            SVM = generator
        else: # Not ASM nor SVM
            try:
                SVM = AlternatingSignMatrix(generator).to_six_vertex_model()
            except (TypeError, ValueError):
                generator = matrix(generator)
                generator = SixVertexModel(generator.nrows(), boundary_conditions='ice')(generator)
                # Check that this is an ice square model
                M = generator.to_alternating_sign_matrix()
                SVM = generator

        if not SVM:
            raise TypeError('generator for FullyPackedLoop must either be an \
            AlternatingSignMatrix or a SquareIceModel.Element')
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

        self.configuration = matrix(list(self._six_vertex_model))
        self._n = len(self._end_point_dictionary) // 2
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
        n=len(self._six_vertex_model)-1
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
                if (i+j) % 2 == 0:
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
                if (i+j) % 2 == 0:
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
                if (i+j) % 2 ==0:
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

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A.random_element()
            sage: FullyPackedLoop(M) == M.to_fully_packed_loop()
            True

            sage: FullyPackedLoop(A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])) ==\
            FullyPackedLoop(A([[1, 0, 0],[0, 0, 1],[0, 1, 0]]))
            False

            sage: FullyPackedLoop(M) == M
            False
        """
        return repr(self) == repr(other) and \
        self._end_point_dictionary == other._end_point_dictionary\
        and self._six_vertex_model == other._six_vertex_model

    def to_alternating_sign_matrix(self):
        """
        Returns the alternating sign matrix corresponding to this class.

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

    def plot(self):
        r"""
        Return a graphical object of the Fully Packed Loop

        EXAMPLES:

        Here is the fully packed loop for

        .. MATH::

            \begin{pmatrix} 0&1&1 \\ 1&-1&1 \\ 0&1&0 \end{pmatrix}:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        Here is how Sage represents this::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: print(fpl.plot().description())
            Line defined by 2 points:       [(-1.0, 1.0), (0.0, 1.0)]
            Line defined by 2 points:       [(0.0, 0.0), (0.0, -1.0)]
            Line defined by 2 points:       [(0.0, 0.0), (1.0, 0.0)]
            Line defined by 2 points:       [(0.0, 2.0), (0.0, 3.0)]
            Line defined by 2 points:       [(0.0, 2.0), (0.0, 3.0)]
            Line defined by 2 points:       [(0.0, 2.0), (1.0, 2.0)]
            Line defined by 2 points:       [(1.0, 1.0), (0.0, 1.0)]
            Line defined by 2 points:       [(1.0, 1.0), (2.0, 1.0)]
            Line defined by 2 points:       [(2.0, 0.0), (1.0, 0.0)]
            Line defined by 2 points:       [(2.0, 0.0), (2.0, -1.0)]
            Line defined by 2 points:       [(2.0, 2.0), (1.0, 2.0)]
            Line defined by 2 points:       [(2.0, 2.0), (2.0, 3.0)]
            Line defined by 2 points:       [(2.0, 2.0), (2.0, 3.0)]
            Line defined by 2 points:       [(3.0, 1.0), (2.0, 1.0)]
            Line defined by 2 points:       [(3.0, 1.0), (2.0, 1.0)]

        Here are the other 3 by 3 Alternating Sign Matrices and their corresponding
        fully packed loops:

        .. MATH::

            A = \begin{pmatrix} 1&0&0 \\ 0&1&0 \\ 0&0&1 \\ \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. MATH::

            A = \begin{pmatrix} 1&0&0 \\ 0&0&1 \\ 0&1&0 \\ \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. MATH::

            A = \begin{pmatrix} 0&1&0\\ 1&0&0\\ 0&0&1\\ \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. MATH::

            A = \begin{pmatrix} 0&1&0\\ 0&0&1\\ 1&0&0\\ \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. MATH::

            A = \begin{pmatrix} 0&0&1\\ 1&0&0\\ 0&1&0\\ \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. MATH::

            A = \begin{pmatrix} 0&0&1\\ 0&1&0\\ 1&0&0\\ \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[0, 1, 0, 0], [1, -1, 0, 1], \
            [0, 1, 0, 0],[0, 0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: print(fpl.plot().description())
            Line defined by 2 points:       [(-1.0, 0.0), (0.0, 0.0)]
            Line defined by 2 points:       [(-1.0, 2.0), (0.0, 2.0)]
            Line defined by 2 points:       [(0.0, 1.0), (0.0, 0.0)]
            Line defined by 2 points:       [(0.0, 1.0), (1.0, 1.0)]
            Line defined by 2 points:       [(0.0, 3.0), (0.0, 4.0)]
            Line defined by 2 points:       [(0.0, 3.0), (0.0, 4.0)]
            Line defined by 2 points:       [(0.0, 3.0), (1.0, 3.0)]
            Line defined by 2 points:       [(1.0, 0.0), (1.0, -1.0)]
            Line defined by 2 points:       [(1.0, 0.0), (2.0, 0.0)]
            Line defined by 2 points:       [(1.0, 2.0), (0.0, 2.0)]
            Line defined by 2 points:       [(1.0, 2.0), (2.0, 2.0)]
            Line defined by 2 points:       [(2.0, 1.0), (1.0, 1.0)]
            Line defined by 2 points:       [(2.0, 1.0), (2.0, 2.0)]
            Line defined by 2 points:       [(2.0, 3.0), (1.0, 3.0)]
            Line defined by 2 points:       [(2.0, 3.0), (2.0, 4.0)]
            Line defined by 2 points:       [(2.0, 3.0), (2.0, 4.0)]
            Line defined by 2 points:       [(3.0, 0.0), (2.0, 0.0)]
            Line defined by 2 points:       [(3.0, 0.0), (3.0, -1.0)]
            Line defined by 2 points:       [(3.0, 2.0), (3.0, 1.0)]
            Line defined by 2 points:       [(3.0, 2.0), (3.0, 3.0)]
            Line defined by 2 points:       [(4.0, 1.0), (3.0, 1.0)]
            Line defined by 2 points:       [(4.0, 1.0), (3.0, 1.0)]
            Line defined by 2 points:       [(4.0, 3.0), (3.0, 3.0)]
            Line defined by 2 points:       [(4.0, 3.0), (3.0, 3.0)]

        Here is the plot:

        .. PLOT::
            :width: 300 px

            A = AlternatingSignMatrix([[0, 1, 0, 0], [1, -1, 0, 1], [0, 1, 0, 0],[0, 0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        """
        G = Graphics()
        n=len(self._six_vertex_model)-1
        for j,row in enumerate(reversed(self._six_vertex_model)):
            for i,entry in enumerate(row):
                if i == 0 and (i+j+n+1) % 2 ==0:
                    G+= line([(i-1,j),(i,j)])
                if i == n and (i+j+n+1) % 2 ==0:
                    G+= line([(i+1,j),(i,j)])
                if j == 0 and (i+j+n) % 2 ==0:
                    G+= line([(i,j),(i,j-1)])
                if j == n and (i+j+n) % 2 ==0:
                    G+= line([(i,j),(i,j+1)])
                if entry == 0: # LR
                    if (i+j+n) % 2==0:
                        G += line([(i,j), (i+1,j)])
                    else:
                        G += line([(i,j),(i,j+1)])
                elif entry == 1: # LU
                    if (i+j+n) % 2 ==0:
                        G += line([(i,j), (i,j+1)])
                    else:
                        G += line([(i+1,j), (i,j)])
                elif entry == 2: # LD
                    if (i+j+n) % 2 == 0:
                        pass
                    else:
                        G += line([(i,j+1), (i,j)])
                        G += line([(i+1,j), (i,j)])
                elif entry == 3: # UD
                    if (i+j+n) % 2 == 0:
                        G += line([(i,j), (i,j+1)])
                    else:
                        G += line([(i+1,j), (i,j)])
                elif entry == 4: # UR
                    if (i+j+n) % 2 ==0:
                        G += line([(i,j), (i,j+1)])
                        G += line([(i,j), (i+1,j)])
                    else:
                        pass
                elif entry == 5: # RD
                    if (i+j+n) % 2 ==0:
                        G += line([(i,j), (i+1,j)])
                    else:
                        G += line([(i,j+1), (i,j)])
        G.axes(False)
        return G

    def gyration(self):
        r"""
        Return the fully packed loop obtained by applying gyration
        to the alternating sign matrix in bijection with ``self``.

        Gyration was first defined in [Wieland00]_ as an action on 
        fully-packed loops.

        REFERENCES:

        .. [Wieland00] \B. Wieland. *A large dihedral symmetry of the set of
           alternating sign matrices*. Electron. J. Combin. 7 (2000).

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


    def link_pattern(self):
        """
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
            See [Propp2001]_ and [Striker2015]_.

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
            ....:     for i in range(0,5):
            ....:         a,b=a%6+1,b%6+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(ASMs[1].gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
            True

            sage: fpl = FullyPackedLoop(ASMs[0])
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(0,5):
            ....:         a,b=a%6+1,b%6+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(ASMs[0].gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
            True

            sage: mat = AlternatingSignMatrix([[0,0,1,0,0,0,0],[1,0,-1,0,1,0,0],[0,0,1,0,0,0,0],\
            [0,1,-1,0,0,1,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,0,1]])
            sage: fpl = FullyPackedLoop(mat) # n=7
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(0,13):
            ....:         a,b=a%14+1,b%14+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
            True

            sage: mat = AlternatingSignMatrix([[0,0,0,1,0,0], [0,0,1,-1,1,0], [0,1,0,0,-1,1], [1,0,-1,1,0,0], \
            [0,0,1,0,0,0], [0,0,0,0,1,0]])
            sage: fpl = FullyPackedLoop(mat)
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(0,11):
            ....:         a,b=a%12+1,b%12+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
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
        link_pattern = []
        boundary_d = self._end_point_dictionary.copy()
        vertices_d = self._vertex_dictionary.copy()

        while len(boundary_d) > 2:
            startpoint = boundary_d.keys()[0]
            position = boundary_d[startpoint]

            boundary_d.pop(startpoint)
            vertices_d[position] = False # allows us to start

            while not vertices_d[position]:
                vertices_d.pop(position)
                choices = self._get_coordinates(position)
                if choices[0] in vertices_d:
                    position = choices[0]
                elif choices[1] in vertices_d:
                    position = choices[1]
                else:
                    raise ValueError('No valid choices')

            endpoint = vertices_d[position]
            vertices_d.pop(position)
            link_pattern.append((startpoint, endpoint))
            boundary_d.pop(endpoint)

        link_pattern.append(tuple(boundary_d.keys()))

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

    @lazy_attribute
    def _vertex_dictionary(self):
        """
        A helper function for :meth:`link_pattern`.
        Return a dictionary where the keys are the coordinates of each vertex and
        the values are either 0 or the values 1 ,..., 2n. The vertices connected to
        endpoints of paths have values 1, ..., 2n. Other vertices have values 0.

        TESTS::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl._vertex_dictionary
            {(0, 0): 1,
             (0, 1): 0,
             (0, 2): 2,
             (1, 0): 6,
             (1, 1): 0,
             (1, 2): 3,
             (2, 0): 5,
             (2, 1): 0,
             (2, 2): 4}
             sage: fpl._vertex_dictionary[(2,2)]
             4
        """
        n = len(self._six_vertex_model)
        vertices = {}
        for i in range(n):
            for j in range(n):
                vertices[(i, j)] = 0

        for end, vertex in self._end_point_dictionary.iteritems():
            vertices[vertex] = end

        return vertices

    def _get_coordinates(self, current_pos):
        """
        Return a list of 2 coordinates that refer to the moves that could
        potentialy be made.

        TESTS::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: matrix(list(fpl._six_vertex_model))
            [3 1 1]
            [5 3 1]
            [5 5 3]
            sage: fpl._get_coordinates((0, 1))
            [(1, 1), (0, 2)]
            sage: fpl._get_coordinates((1, 1))
            [(2, 1), (0, 1)]
            sage: fpl._get_coordinates((0, 0))
            [(1, 0), (-1, 0)]
            sage: fpl._get_coordinates((0, 2))
            [(-1, 2), (0, 1)]
            sage: fpl._get_coordinates((2, 1))
            [(1, 1), (2, 0)]

            sage: B = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(B)
            sage: matrix(list(fpl._six_vertex_model))
            [4 3 1]
            [3 0 3]
            [5 3 2]
            sage: fpl._get_coordinates((1, 1))
            [(1, 0), (1, 2)]
            sage: fpl._get_coordinates((0, 0))
            [(-1, 0), (0, 1)]
            sage: fpl._get_coordinates((0, 2))
            [(-1, 2), (0, 1)]
            sage: fpl._get_coordinates((2, 1))
            [(2, 0), (2, 2)]
        """
        # 0 UD, 1 RD, 2 UR, 3 LR, 4 LD, 5 LU
        odd = {0: [(1, 0), (-1, 0)],
                1: [(1, 0), (0, 1)],
                2: [(-1, 0), (0, 1)],
                3: [(0, -1), (0, 1)],
                4: [(0, -1), (1, 0)],
                5: [(-1, 0), (0, -1)]
                }

        # 0 LR, 1 LU, 2 LD, 3 UD, 4 UR, 5 RD
        even = {0: [(0, -1), (0, 1)],
               1: [(-1, 0), (0, -1)],
               2: [(0, -1), (1, 0)],
               3: [(1, 0), (-1, 0)],
               4: [(-1, 0), (0, 1)],
               5: [(1, 0), (0, 1)]
               }

        parity = sum(current_pos)
        conf = self.configuration[current_pos]

        if parity % 2 == 0:
            potential_directions = even[conf]
        else:
            potential_directions = odd[conf]

        return [(current_pos[0] + d[0], current_pos[1] + d[1]) for d in potential_directions]

    @lazy_attribute
    def _end_point_dictionary(self):
        r"""
        A helper function for :meth:`link_pattern`.
        Return a dictionary where the keys are the labels 1, ..., 2n and the values
        are the coordinates of the endpoints of the paths of ``self``.

        TESTS::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl._end_point_dictionary
            {1: (0, 0), 2: (0, 2), 3: (1, 2), 4: (2, 2), 5: (2, 0), 6: (1, 0)}
            sage: fpl._end_point_dictionary[2]
            (0, 2)

        """
        n = len(self._six_vertex_model)
        end_points = {}

        for k in range(n):
            if k % 2 == 0:
                # top row
                end_points[1 + k // 2] = (0, k)

                # bottom row
                end_points[n + 1 + k // 2] = (n-1, n-1-k)

        # sides for even case
        if n % 2 == 0:
            for k in range(n):
                if k % 2 == 0:
                    # left side
                    end_points[((3*n + 2 + k) // 2)] = (n-1-k, 0)

                    # right side
                    end_points[(n + 2 + k) // 2] = (k, n-1)

        # side for odd case
        if n % 2 == 1:
            for k in range(n):
                if k % 2 == 1:
                    # left side
                    end_points[(3*n + 2 + k) // 2] = (n-1-k, 0)

                    # right side
                    end_points[(n + 2 + k) // 2] = (k, n-1)

        return end_points

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
        if isinstance(fpl, FullyPackedLoop):
            return fpl._n == self._n
        return False

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
        else: # Not ASM nor SVM
            try:
                SVM = AlternatingSignMatrix(generator).to_six_vertex_model()
            except (TypeError, ValueError):
                SVM = SixVertexModel(self._n, boundary_conditions='ice')(generator)
                M = SVM.to_alternating_sign_matrix()

        if len(SVM) != self._n:
            raise ValueError("invalid size")
        return self.element_class(self,SVM)

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

            sage: [AlternatingSignMatrices(n).cardinality() for n in range(0, 11)]
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
