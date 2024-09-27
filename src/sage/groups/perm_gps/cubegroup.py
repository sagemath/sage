r"""
Rubik's cube group functions

.. _sec-rubik:

.. NOTE::

    "Rubiks cube" is trademarked. We shall omit the trademark
    symbol below for simplicity.

NOTATION:

`B` denotes a clockwise quarter turn of the back face, `D`
denotes a clockwise quarter turn of the down face, and similarly for
`F` (front), `L` (left), `R` (right), and `U` (up). Products of moves are read
right to left, so for example, `R \cdot U` means move `U` first and then `R`.

See ``CubeGroup.parse()`` for all possible input notations.

The "Singmaster notation":

- moves: `U, D, R, L, F, B` as in the
  diagram below,

- corners: `xyz` means the facet is on face `x` (in `R,F,L,U,D,B`) and the
  clockwise rotation of the corner sends `x-y-z`

- edges: `xy` means the facet is on face `x` and a flip of the edge sends
  `x-y`.

::

    sage: rubik = CubeGroup()
    sage: rubik.display2d("")
                 ┌──────────────┐
                 │  1    2    3 │
                 │  4   top   5 │
                 │  6    7    8 │
    ┌────────────┼──────────────┼─────────────┬────────────┐
    │  9  10  11 │ 17   18   19 │ 25   26  27 │ 33  34  35 │
    │ 12 left 13 │ 20  front 21 │ 28 right 29 │ 36 rear 37 │
    │ 14  15  16 │ 22   23   24 │ 30   31  32 │ 38  39  40 │
    └────────────┼──────────────┼─────────────┴────────────┘
                 │ 41   42   43 │
                 │ 44 bottom 45 │
                 │ 46   47   48 │
                 └──────────────┘

AUTHORS:

- David Joyner (2006-10-21): first version

- David Joyner (2007-05): changed faces, added legal and solve

- David Joyner(2007-06): added plotting functions

- David Joyner (2007, 2008): colors corrected, "solve" rewritten
  (again),typos fixed.

- Robert Miller (2007, 2008): editing, cleaned up display2d

- Robert Bradshaw (2007, 2008): RubiksCube object, 3d
  plotting.

- David Joyner (2007-09): rewrote docstring for CubeGroup's
  "solve".

- Robert Bradshaw (2007-09): Versatile parse function for
  all input types.

- Robert Bradshaw (2007-11): Cleanup.

REFERENCES:

- Cameron, P., Permutation Groups. New York: Cambridge
  University Press, 1999.

- Wielandt, H., Finite Permutation Groups.
  New York: Academic Press, 1964.

- Dixon, J. and Mortimer, B.,
  Permutation Groups, Springer-Verlag, Berlin/New York, 1996.

- Joyner, D., Adventures in Group Theory, Johns Hopkins Univ Press,
  2002.
"""

# *****************************************************************************
#       Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from random import randint

from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp, richcmp_method
from sage.libs.gap.libgap import libgap
from sage.rings.real_double import RDF
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.polygon", "polygon")
lazy_import("sage.plot.text", "text")
pi = RDF.pi()
from sage.rings.rational_field import QQ

lazy_import("sage.plot.plot3d.shapes", "Box")
lazy_import("sage.plot.plot3d.texture", "Texture")

####################### predefined colors ##################

named_colors = {
    'red': (1, 0, 0),              ## F face
    'green': (0, 1, 0),            ## R face
    'blue': (0, 0, 1),             ## D face
    'yellow': (1, 1, 0),           ## L face
    'white': (1, 1, 1),            ## none
    'orange': (1, 0.6, 0.3),       ## B face
    'purple': (1, 0, 1),           ## none
    'lpurple': (1, 0.63, 1),       ## U face
    'lightblue': (0, 1, 1),        ## none
    'lgrey': (0.75, 0.75, 0.75),   ## sagemath.org color
}
globals().update(named_colors)

#########################################################
# written by Tom Boothby, placed in the public domain

def xproj(x, y, z, r):
    r"""
    Return the `x`-projection of `(x,y,z)` rotated by `r`.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import rotation_list, xproj
        sage: rot = rotation_list(30, 45)
        sage: xproj(1,2,3,rot)
        0.6123724356957945
    """
    return (y*r[1] - x*r[3])*r[2]


def yproj(x, y, z, r):
    r"""
    Return the `y`-projection of `(x,y,z)` rotated by `r`.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import rotation_list, yproj
        sage: rot = rotation_list(30, 45)
        sage: yproj(1,2,3,rot)
        1.378497416975604
    """
    return z*r[2] - (x*r[1] + y*r[2])*r[0]


def rotation_list(tilt, turn):
    r"""
    Return a list `[\sin(\theta), \sin(\phi), \cos(\theta), \cos(\phi)]` of
    rotations where `\theta` is ``tilt`` and `\phi` is ``turn``.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import rotation_list
        sage: rotation_list(30, 45)
        [0.49999999999999994, 0.7071067811865475, 0.8660254037844387, 0.7071067811865476]
    """
    from sage.functions.all import sin, cos
    return [sin(tilt*pi/180.0), sin(turn*pi/180.0),
            cos(tilt*pi/180.0), cos(turn*pi/180.0)]


def polygon_plot3d(points, tilt=30, turn=30, **kwargs):
    r"""
    Plot a polygon viewed from an angle determined by ``tilt``, ``turn``, and
    vertices ``points``.

    .. WARNING::

       The ordering of the points is important to get "correct"
       and if you add several of these plots together, the one added first
       is also drawn first (ie, addition of Graphics objects is not
       commutative).

    The following example produced a green-colored square with vertices
    at the points indicated.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import polygon_plot3d,green
        sage: P = polygon_plot3d([[1,3,1],[2,3,1],[2,3,2],[1,3,2],[1,3,1]],             # needs sage.plot
        ....:                    rgbcolor=green)
    """
    rot = rotation_list(tilt, turn)
    points2 = [(xproj(x, y, z, rot), yproj(x, y, z, rot))
               for (x, y, z) in points]
    return polygon(points2, **kwargs)

###########################################################

#############  lots of "internal" utility plot functions #########


def inv_list(lst):
    r"""
    Input a list of ints `1, \ldots, m` (in any order), outputs inverse
    perm.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import inv_list
        sage: L = [2,3,1]
        sage: inv_list(L)
        [3, 1, 2]
    """
    return [lst.index(i) + 1 for i in range(1, 1 + len(lst))]


face_polys = {
### bottom layer L, F, R, B
    'ldb': [[-3,0],[-2,0], [-2,1], [-3,1]],      #square labeled 14
    'ld': [[-2,0],[-1,0], [-1,1], [-2,1]],      #square labeled 15
    'lfd': [[-1,0],[0,0], [0,1], [-1,1]],      #square labeled 16
    'fdl': [[0,0],[1,0], [1,1], [0,1]],      #square labeled 22
    'fd': [[1,0],[2,0], [2,1], [1,1]],      #square labeled 23
    'frd': [[2,0],[3,0], [3,1], [2,1]],      #square labeled 24
    'rdf': [[3,0],[4,0], [4,1], [3,1]],      #square labeled 30
    'rd': [[4,0],[5,0], [5,1], [4,1]],      #square labeled 31
    'rbd': [[5,0],[6,0], [6,1], [5,1]],      #square labeled 32
    'bdr': [[6,0],[7,0], [7,1], [6,1]],      #square labeled 38
    'bd': [[7,0],[8,0], [8,1], [7,1]],      #square labeled 39
    'bld': [[8,0],[9,0], [9,1], [8,1]],      #square labeled 40
### middle layer L,F,R, B
    'lb': [[-3,1],[-2,1], [-2,2], [-3,2]],      #square labeled 12
    'l_center': [[-2,1],[-1,1], [-1,2], [-2,2]],        #center square
    'lf': [[-1,1],[0,1], [0,2], [-1,2]],      #square labeled 13
    'fl': [[0,1],[1,1], [1,2], [0,2]],      #square labeled 20
    'f_center': [[1,1],[2,1], [2,2], [1,2]],        #center square
    'fr': [[2,1],[3,1], [3,2], [2,2]],      #square labeled 21
    'rf': [[3,1],[4,1], [4,2], [3,2]],      #square labeled 28
    'r_center': [[4,1],[5,1], [5,2], [4,2]],        #center square
    'rb': [[5,1],[6,1], [6,2], [5,2]],      #square labeled 29
    'br': [[6,1],[7,1], [7,2], [6,2]],      #square labeled 36
    'b_center': [[7,1],[8,1], [8,2], [7,2]],        #center square
    'bl': [[8,1],[9,1], [9,2], [8,2]],      #square labeled 37
## top layer L, F, R, B
    'lbu': [[-3,2],[-2,2], [-2,3], [-3,3]],      #square labeled 9
    'lu': [[-2,2],[-1,2], [-1,3], [-2,3]],      #square labeled 10
    'luf': [[-1,2],[0,2], [0,3], [-1,3]],      #square labeled 11
    'flu': [[0,2],[1,2], [1,3], [0,3]],      #square labeled 17
    'fu': [[1,2],[2,2], [2,3], [1,3]],      #square labeled 18
    'fur': [[2,2],[3,2], [3,3], [2,3]],      #square labeled 19
    'ruf': [[3,2],[4,2], [4,3], [3,3]],      #square labeled 25
    'ru': [[4,2],[5,2], [5,3], [4,3]],      #square labeled 26
    'rub': [[5,2],[6,2], [6,3], [5,3]],      #square labeled 27
    'bur': [[6,2],[7,2], [7,3], [6,3]],      #square labeled 33
    'bu': [[7,2],[8,2], [8,3], [7,3]],      #square labeled 34
    'bul': [[8,2],[9,2], [9,3], [8,3]],      #square labeled 35
# down face
    'dlf': [[0,-1],[1,-1], [1,0], [0,0]],      #square labeled 41
    'df': [[1,-1],[2,-1], [2,0], [1,0]],      #square labeled 42
    'dfr': [[2,-1],[3,-1], [3,0], [2,0]],      #square labeled 43
    'dl': [[0,-2],[1,-2], [1,-1], [0,-1]],      #square labeled 44
    'd_center': [[1,-2],[2,-2], [2,-1], [1,-1]],        #center square
    'dr': [[2,-2],[3,-2], [3,-1], [2,-1]],      #square labeled 45
    'dlb': [[0,-3],[1,-3], [1,-2], [0,-2]],      #square labeled 46
    'db': [[1,-3],[2,-3], [2,-2], [1,-2]],      #square labeled 47
    'drb': [[2,-3],[3,-3], [3,-2], [2,-2]],      #square labeled 48
# up face
    'ufl': [[0,3],[1,3], [1,4], [0,4]],      #square labeled 6
    'uf': [[1,3],[2,3], [2,4], [1,4]],      #square labeled 7
    'urf': [[2,3],[3,3], [3,4], [2,4]],      #square labeled 8
    'ul': [[0,4],[1,4], [1,5], [0,5]],      #square labeled 4
    'u_center': [[1,4],[2,4], [2,5], [1,5]],        #center square
    'ur': [[2,4],[3,4], [3,5], [2,5]],      #square labeled 5
    'ulb': [[0,6],[1,6], [1,5], [0,5]],      #square labeled 1
    'ub': [[1,6],[2,6], [2,5], [1,5]],      #square labeled 2
    'ubr': [[2,6],[3,6], [3,5], [2,5]],      #square labeled 3
}


def create_poly(face, color):
    """
    Create the polygon given by ``face`` with color ``color``.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import create_poly, red
        sage: create_poly('ur', red)                                                    # needs sage.plot
        Graphics object consisting of 1 graphics primitive
    """
    return polygon(face_polys[face], rgbcolor=color)

####################################################


singmaster_indices = {
    1: "ulb",
    2: "ub",
    3: "ubr",
    4: "ul",
    5: "ur",
    6: "ufl",
    7: "uf",
    8: "urf",
    14: "ldb",
    15: "ld",
    16: "lfd",
    12: "lb",
    13: "lf",
    9: "lbu",
    10: "lu",
    11: "luf",
    17: "flu",
    18: "fu",
    19: "fur",
    20: "fl",
    21: "fr",
    22: "fdl",
    23: "fd",
    24: "frd",
    41: "dlf",
    42: "df",
    43: "dfr",
    44: "dl",
    45: "dr",
    46: "dlb",
    47: "db",
    48: "drb",
    33: "bur",
    34: "bu",
    35: "bul",
    36: "br",
    37: "bl",
    38: "bdr",
    39: "bd",
    40: "bld",
    25: "ruf",
    26: "ru",
    27: "rub",
    28: "rf",
    29: "rb",
    30: "rdf",
    31: "rd",
    32: "rbd",
}

def index2singmaster(facet):
    """
    Translate index used (eg, 43) to Singmaster facet notation (eg,
    fdr).

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import index2singmaster
        sage: index2singmaster(41)
        'dlf'
    """
    return singmaster_indices[facet]


def color_of_square(facet, colors=['lpurple', 'yellow', 'red', 'green', 'orange', 'blue']):
    """
    Return the color the facet has in the solved state.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import color_of_square
        sage: color_of_square(41)
        'blue'
    """
    return colors[(facet - 1) // 8]


cubie_center_list = {
    #  centers of the cubies on the F,U, R faces
    1:  [1//2, 1//2, 5//2], # ulb
    2:  [1//2, 3//2, 5//2], #  ub
    3:  [1//2, 5//2, 5//2], #  ubr
    4:  [3//2, 1//2, 5//2], #  ul
    5:  [3//2, 5//2, 5//2], #  ur
    6:  [5//2, 1//2, 5//2], #  ufl
    7:  [5//2, 3//2, 5//2], #  uf
    8:  [5//2, 5//2, 5//2], #  urf
    17: [5//2, 1//2, 5//2], #  flu
    18: [5//2, 3//2, 5//2], #  fu
    19: [5//2, 5//2, 5//2], #  fur
    20: [5//2, 1//2, 3//2], #  fl
    21: [5//2, 5//2, 3//2], #  fr
    22: [5//2, 1//2, 1//2], #  fdl
    23: [5//2, 3//2, 1//2], #  fd
    24: [5//2, 5//2, 1//2], #  frd
    25: [5//2, 5//2, 5//2], # rfu
    26: [3//2, 5//2, 5//2], #  ru
    27: [1//2, 5//2, 5//2], #  rub
    28: [5//2, 5//2, 3//2], #  rf
    29: [1//2, 5//2, 3//2], #  rb
    30: [5//2, 5//2, 1//2], #  rdf
    31: [3//2, 5//2, 1//2], #  rd
    32: [1//2, 5//2, 1//2], # rbd
}


def cubie_centers(label):
    r"""
    Return the cubie center list element given by ``label``.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import cubie_centers
        sage: cubie_centers(3)
        [0, 2, 2]
    """
    return cubie_center_list[label]


def cubie_colors(label, state0):
    r"""
    Return the color of the cubie given by ``label`` at ``state0``.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import cubie_colors
        sage: G = CubeGroup()
        sage: g = G.parse("R*U")
        sage: cubie_colors(3, G.facets(g))
        [(1, 1, 1), (1, 0.63, 1), (1, 0.6, 0.3)]
    """
    #  colors of the cubies on the F,U, R faces
    clr_any = named_colors['white']
    state = inv_list(state0)
    if label == 1:
        return [clr_any, named_colors[color_of_square(state[1-1])], clr_any] #ulb,
    if label == 2:
        return [clr_any,named_colors[color_of_square(state[2-1])],clr_any] # ub,
    if label == 3:
        return [clr_any, named_colors[color_of_square(state[3-1])], named_colors[color_of_square(state[27-1])]] # ubr,
    if label == 4:
        return [clr_any, named_colors[color_of_square(state[4-1])], clr_any] # ul,
    if label == 5:
        return [clr_any, named_colors[color_of_square(state[5-1])], named_colors[color_of_square(state[26-1])]] # ur,
    if label == 6:
        return [named_colors[color_of_square(state[17-1])], named_colors[color_of_square(state[6-1])], clr_any] # ufl,
    if label == 7:
        return [named_colors[color_of_square(state[18-1])], named_colors[color_of_square(state[7-1])], clr_any] # uf,
    if label == 8:
        return [named_colors[color_of_square(state[19-1])], named_colors[color_of_square(state[8-1])], named_colors[color_of_square(state[25-1])]] # urf,
    if label == 17:
        return [named_colors[color_of_square(state[17-1])], named_colors[color_of_square(state[6-1])], clr_any] # flu
    if label == 18:
        return [named_colors[color_of_square(state[18-1])], named_colors[color_of_square(state[7-1])], clr_any] # fu
    if label == 19:
        return [named_colors[color_of_square(state[19-1])], named_colors[color_of_square(state[8-1])], named_colors[color_of_square(state[25-1])]] # fur
    if label == 20:
        return [named_colors[color_of_square(state[20-1])], clr_any, clr_any] # fl
    if label == 21:
        return [named_colors[color_of_square(state[21-1])], clr_any, named_colors[color_of_square(state[28-1])]] # fr
    if label == 22:
        return [named_colors[color_of_square(state[22-1])], clr_any, clr_any] # fdl
    if label == 23:
        return [named_colors[color_of_square(state[23-1])], clr_any, clr_any] # fd
    if label == 24:
        return [named_colors[color_of_square(state[24-1])], clr_any, named_colors[color_of_square(state[30-1])]] # frd
    if label == 25:
        return [named_colors[color_of_square(state[19-1])],named_colors[color_of_square(state[8-1])],named_colors[color_of_square(state[25-1])]]  #rfu,
    if label == 26:
        return [clr_any,named_colors[color_of_square(state[5-1])],named_colors[color_of_square(state[26-1])]] # ru,
    if label == 27:
        return [clr_any,named_colors[color_of_square(state[3-1])],named_colors[color_of_square(state[27-1])]] # rub,
    if label == 28:
        return [named_colors[color_of_square(state[21-1])],clr_any,named_colors[color_of_square(state[28-1])]] # rf,
    if label == 29:
        return [clr_any,clr_any,named_colors[color_of_square(state[29-1])]] # rb,
    if label == 30:
        return [named_colors[color_of_square(state[24-1])],clr_any,named_colors[color_of_square(state[30-1])]] # rdf,
    if label == 31:
        return [clr_any,clr_any,named_colors[color_of_square(state[31-1])]] # rd,
    if label == 32:
        return [clr_any,clr_any,named_colors[color_of_square(state[32-1])]] #rbd,


def plot3d_cubie(cnt, clrs):
    r"""
    Plot the front, up and right face of a cubie centered at cnt and
    rgbcolors given by clrs (in the order FUR).

    Type ``P.show()`` to view.

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import plot3d_cubie, blue, red, green
        sage: clrF = blue; clrU = red; clrR = green
        sage: P = plot3d_cubie([1/2,1/2,1/2],[clrF,clrU,clrR])                          # needs sage.plot
    """
    half = QQ((1, 2))
    x = cnt[0] - half
    y = cnt[1] - half
    z = cnt[2] - half
    #ptsD = [[x+0,y+0,0+z],[x+1,y+0,0+z],[x+1,y+1,0+z],[x+0,y+1,0+z],[x+0,y+0,0+z]]
    ptsF = [[x+1,y+0,0+z],[x+1,y+1,0+z],[x+1,y+1,1+z],[x+1,y+0,1+z],[x+1,y+0,0+z]]
    #ptsB = [[x+0,y+0,0+z],[x+0,y+1,0+z],[x+0,y+1,1+z],[x+0,y+0,1+z],[x+0,y+0,0+z]]
    ptsU = [[x+0,y+0,1+z],[x+1,y+0,1+z],[x+1,y+1,1+z],[x+0,y+1,1+z],[x+0,y+0,1+z]]
    #ptsL = [[x+0,y+0,0+z],[x+1,y+0,0+z],[x+1,y+0,1+z],[x+0,y+0,1+z],[x+0,y+0,0+z]]
    ptsR = [[x+0,y+1,0+z],[x+1,y+1,0+z],[x+1,y+1,1+z],[x+0,y+1,1+z],[x+0,y+1,0+z]]
    P = polygon_plot3d(ptsR, rgbcolor=clrs[2])
    P += polygon_plot3d(ptsU, rgbcolor=clrs[1])
    P += polygon_plot3d(ptsF, rgbcolor=clrs[0])
    P.axes(show=False)
    return P


####################### end of "internal" utility plot functions  #################


class CubeGroup(PermutationGroup_generic):
    r"""
    A python class to help compute Rubik's cube group actions.

    .. NOTE::

        This group is also available via ``groups.permutation.RubiksCube()``.

    EXAMPLES:

    If G denotes the cube group then it may be regarded as a
    subgroup of ``SymmetricGroup(48)``, where the 48 facets are labeled as
    follows.

    ::

        sage: rubik = CubeGroup()
        sage: rubik.display2d("")
                     ┌──────────────┐
                     │  1    2    3 │
                     │  4   top   5 │
                     │  6    7    8 │
        ┌────────────┼──────────────┼─────────────┬────────────┐
        │  9  10  11 │ 17   18   19 │ 25   26  27 │ 33  34  35 │
        │ 12 left 13 │ 20  front 21 │ 28 right 29 │ 36 rear 37 │
        │ 14  15  16 │ 22   23   24 │ 30   31  32 │ 38  39  40 │
        └────────────┼──────────────┼─────────────┴────────────┘
                     │ 41   42   43 │
                     │ 44 bottom 45 │
                     │ 46   47   48 │
                     └──────────────┘

    ::

        sage: rubik
        The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48).

        TESTS::

            sage: groups.permutation.RubiksCube()
            The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48).
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: TestSuite(rubik).run(skip='_test_enumerated_set_contains') # because the group is very large

        TESTS:

        Check that :issue:`11360` is fixed::

            sage: rubik = CubeGroup()
            sage: rubik.order()
            43252003274489856000
        """
        U = "( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19)"  # U = top
        L = "( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35)"  # L = left
        F = "(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11)"  # F = front
        R = "(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24)"  # R = right
        B = "(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27)"  # B = back or rear
        D = "(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40)"  # D = down or bottom
        PermutationGroup_generic.__init__(self, gens=[B, D, F, L, R, U], canonicalize=False)

    def gen_names(self):
        """
        Return the names of the generators.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.gen_names()
            ['B', 'D', 'F', 'L', 'R', 'U']
        """
        return ['B','D','F','L','R','U']

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik
            The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48).
        """
        return "The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48)."

    def B(self):
        """
        Return the generator `B` in Singmaster notation.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.B()
            (1,14,48,27)(2,12,47,29)(3,9,46,32)(33,35,40,38)(34,37,39,36)
        """
        return self.gens()[0]

    def D(self):
        """
        Return the generator `D` in Singmaster notation.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.D()
            (14,22,30,38)(15,23,31,39)(16,24,32,40)(41,43,48,46)(42,45,47,44)
        """
        return self.gens()[1]

    def F(self):
        """
        Return the generator `F` in Singmaster notation.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.F()
            (6,25,43,16)(7,28,42,13)(8,30,41,11)(17,19,24,22)(18,21,23,20)
        """
        return self.gens()[2]

    def L(self):
        """
        Return the generator `L` in Singmaster notation.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.L()
            (1,17,41,40)(4,20,44,37)(6,22,46,35)(9,11,16,14)(10,13,15,12)
        """
        return self.gens()[3]

    def R(self):
        """
        Return the generator `R` in Singmaster notation.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.R()
            (3,38,43,19)(5,36,45,21)(8,33,48,24)(25,27,32,30)(26,29,31,28)
        """
        return self.gens()[4]

    def U(self):
        """
        Return the generator `U` in Singmaster notation.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.U()
            (1,3,8,6)(2,5,7,4)(9,33,25,17)(10,34,26,18)(11,35,27,19)
        """
        return self.gens()[5]

    def parse(self, mv, check=True):
        r"""
        This function allows one to create the permutation group element
        from a variety of formats.

        INPUT:

        - ``mv`` -- can one of the following:

          - ``list`` -- list of facets (as returned by
            self.facets())

          - ``dict`` -- list of faces (as returned by
            ``self.faces()``)

          - ``str`` -- either cycle notation (passed to GAP) or
            a product of generators or Singmaster notation

          - ``perm_group element`` -- returned as an element of ``self``

        - ``check`` -- check if the input is valid

        EXAMPLES::

            sage: C = CubeGroup()
            sage: C.parse(list(range(1,49)))
            ()
            sage: g = C.parse("L"); g
            (1,17,41,40)(4,20,44,37)(6,22,46,35)(9,11,16,14)(10,13,15,12)
            sage: C.parse(str(g)) == g
            True
            sage: facets = C.facets(g); facets
            [17, 2, 3, 20, 5, 22, 7, 8, 11, 13, 16, 10, 15, 9, 12, 14, 41, 18, 19, 44, 21, 46, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 6, 36, 4, 38, 39, 1, 40, 42, 43, 37, 45, 35, 47, 48]
            sage: C.parse(facets)
            (1,17,41,40)(4,20,44,37)(6,22,46,35)(9,11,16,14)(10,13,15,12)
            sage: C.parse(facets) == g
            True
            sage: faces = C.faces("L"); faces
            {'back': [[33, 34, 6], [36, 0, 4], [38, 39, 1]],
             'down': [[40, 42, 43], [37, 0, 45], [35, 47, 48]],
             'front': [[41, 18, 19], [44, 0, 21], [46, 23, 24]],
             'left': [[11, 13, 16], [10, 0, 15], [9, 12, 14]],
             'right': [[25, 26, 27], [28, 0, 29], [30, 31, 32]],
             'up': [[17, 2, 3], [20, 0, 5], [22, 7, 8]]}
            sage: C.parse(faces) == C.parse("L")
            True
            sage: C.parse("L' R2") == C.parse("L^(-1)*R^2")
            True
            sage: C.parse("L' R2")
            (1,40,41,17)(3,43)(4,37,44,20)(5,45)(6,35,46,22)(8,48)(9,14,16,11)(10,12,15,13)(19,38)(21,36)(24,33)(25,32)(26,31)(27,30)(28,29)
            sage: C.parse("L^4")
            ()
            sage: C.parse("L^(-1)*R")
            (1,40,41,17)(3,38,43,19)(4,37,44,20)(5,36,45,21)(6,35,46,22)(8,33,48,24)(9,14,16,11)(10,12,15,13)(25,27,32,30)(26,29,31,28)
        """
        if isinstance(mv, PermutationGroupElement):
            # mv is a perm_group element, return mv
            return mv if mv.parent() is self else PermutationGroup_generic.__call__(self, mv, check)
        elif isinstance(mv, str):
            # It is a string: may be in cycle notation or Rubik's notation
            if '(' in mv and '^' not in mv:
                return PermutationGroup_generic.__call__(self, mv, check)
            else:
                gens = self.gens()
                names = self.gen_names()
                map = {}
                for i in range(6):
                    map[names[i]] = gens[i]
                g = self.identity()
                mv = mv.strip().replace(" ","*").replace("**", "*").replace("'", "-1").replace('^','').replace('(','').replace(')','')
                M = mv.split("*")
                for m in M:
                    if not m:
                        pass
                    elif len(m) == 1:
                        g *= map[m[0]]
                    else:
                        g *= map[m[0]]**int(m[1:])
                return g
        elif isinstance(mv, dict):
            state = mv
            state_facets = []
            keyss = sorted(state.keys())
            for k in keyss:
                r = state[k][0]+state[k][1]+state[k][2]
                r.remove(0)
                state_facets = state_facets + r
            state0 = self.faces("")
            state0_facets = []
            keyss = sorted(state0.keys())
            for k in keyss:
                r = state0[k][0]+state0[k][1]+state0[k][2]
                r.remove(0)
                state0_facets = state0_facets + r
            p1 = [state0_facets.index(x) for x in range(1,49)]
            p2 = [state_facets[j] for j in p1]
            return PermutationGroup_generic.__call__(self, p2, check)
        else:
            return PermutationGroup_generic.__call__(self, mv, check)

    __call__ = parse

    def facets(self, g=None):
        r"""
        Return the set of facets on which the group acts. This function is
        a "constant".

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.facets() == list(range(1,49))
            True
        """
        fcts = range(1, 49)
        if g is not None:
            return [g(i) for i in fcts]
        else:
            return list(fcts)

    def faces(self, mv):
        r"""
        Return the dictionary of faces created by the effect of the move
        mv, which is a string of the form `X^a*Y^b*...`, where
        `X, Y, \ldots` are in `\{R,L,F,B,U,D\}` and
        `a, b, \ldots` are integers. We call this ordering of the faces
        the "BDFLRU, L2R, T2B ordering".

        EXAMPLES::

            sage: rubik = CubeGroup()

        Here is the dictionary of the solved state::

            sage: sorted(rubik.faces("").items())
            [('back', [[33, 34, 35], [36, 0, 37], [38, 39, 40]]),
             ('down', [[41, 42, 43], [44, 0, 45], [46, 47, 48]]),
             ('front', [[17, 18, 19], [20, 0, 21], [22, 23, 24]]),
             ('left', [[9, 10, 11], [12, 0, 13], [14, 15, 16]]),
             ('right', [[25, 26, 27], [28, 0, 29], [30, 31, 32]]),
             ('up', [[1, 2, 3], [4, 0, 5], [6, 7, 8]])]

        Now the dictionary of the state obtained after making the move `R`
        followed by `L`::

            sage: sorted(rubik.faces("R*U").items())
            [('back', [[48, 26, 27], [45, 0, 37], [43, 39, 40]]),
             ('down', [[41, 42, 11], [44, 0, 21], [46, 47, 24]]),
             ('front', [[9, 10, 8], [20, 0, 7], [22, 23, 6]]),
             ('left', [[33, 34, 35], [12, 0, 13], [14, 15, 16]]),
             ('right', [[19, 29, 32], [18, 0, 31], [17, 28, 30]]),
             ('up', [[3, 5, 38], [2, 0, 36], [1, 4, 25]])]
        """
        fcts = self.facets(self.parse(mv))
        faceR = [[fcts[24],fcts[25],fcts[26]],[fcts[27],0,fcts[28]],[fcts[29],fcts[30],fcts[31]]]
        faceL = [[fcts[8],fcts[9],fcts[10]],[fcts[11],0,fcts[12]],[fcts[13],fcts[14],fcts[15]]]
        faceU = [[fcts[0],fcts[1],fcts[2]],[fcts[3],0,fcts[4]],[fcts[5],fcts[6],fcts[7]]]
        faceD = [[fcts[40],fcts[41],fcts[42]],[fcts[43],0,fcts[44]],[fcts[45],fcts[46],fcts[47]]]
        faceF = [[fcts[16],fcts[17],fcts[18]],[fcts[19],0,fcts[20]],[fcts[21],fcts[22],fcts[23]]]
        faceB = [[fcts[32],fcts[33],fcts[34]],[fcts[35],0,fcts[36]],[fcts[37],fcts[38],fcts[39]]]
        return {'right': faceR, 'left': faceL, 'up': faceU, 'down': faceD, 'front': faceF, 'back': faceB}

    def move(self, mv):
        r"""
        Return the group element and the reordered list of facets, as
        moved by the list ``mv`` (read left-to-right)

        INPUT:

        - ``mv`` -- string of the form ``Xa*Yb*...``,
          where ``X``, ``Y``, ... are in ``R``, ``L``, ``F``, ``B``, ``U``,
          ``D`` and ``a``, ``b``, ... are integers

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.move("")[0]
            ()
            sage: rubik.move("R")[0]
            (3,38,43,19)(5,36,45,21)(8,33,48,24)(25,27,32,30)(26,29,31,28)
            sage: rubik.R()
            (3,38,43,19)(5,36,45,21)(8,33,48,24)(25,27,32,30)(26,29,31,28)
        """
        g = self.parse(mv)
        return g, self.facets(g)

    def display2d(self, mv):
        r"""
        Print the 2d representation of ``self``.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: rubik.display2d("R")
                         ┌──────────────┐
                         │  1    2   38 │
                         │  4   top  36 │
                         │  6    7   33 │
            ┌────────────┼──────────────┼─────────────┬────────────┐
            │  9  10  11 │ 17   18    3 │ 27   29  32 │ 48  34  35 │
            │ 12 left 13 │ 20  front  5 │ 26 right 31 │ 45 rear 37 │
            │ 14  15  16 │ 22   23    8 │ 25   28  30 │ 43  39  40 │
            └────────────┼──────────────┼─────────────┴────────────┘
                         │ 41   42   19 │
                         │ 44 bottom 21 │
                         │ 46   47   24 │
                         └──────────────┘
        """
        print(self.repr2d(mv))

    def repr2d(self, mv):
        r"""
        Displays a 2D map of the Rubik's cube after the move mv has been
        made. Nothing is returned.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: print(rubik.repr2d(""))
                         ┌──────────────┐
                         │  1    2    3 │
                         │  4   top   5 │
                         │  6    7    8 │
            ┌────────────┼──────────────┼─────────────┬────────────┐
            │  9  10  11 │ 17   18   19 │ 25   26  27 │ 33  34  35 │
            │ 12 left 13 │ 20  front 21 │ 28 right 29 │ 36 rear 37 │
            │ 14  15  16 │ 22   23   24 │ 30   31  32 │ 38  39  40 │
            └────────────┼──────────────┼─────────────┴────────────┘
                         │ 41   42   43 │
                         │ 44 bottom 45 │
                         │ 46   47   48 │
                         └──────────────┘

        ::

            sage: print(rubik.repr2d("R"))
                         ┌──────────────┐
                         │  1    2   38 │
                         │  4   top  36 │
                         │  6    7   33 │
            ┌────────────┼──────────────┼─────────────┬────────────┐
            │  9  10  11 │ 17   18    3 │ 27   29  32 │ 48  34  35 │
            │ 12 left 13 │ 20  front  5 │ 26 right 31 │ 45 rear 37 │
            │ 14  15  16 │ 22   23    8 │ 25   28  30 │ 43  39  40 │
            └────────────┼──────────────┼─────────────┴────────────┘
                         │ 41   42   19 │
                         │ 44 bottom 21 │
                         │ 46   47   24 │
                         └──────────────┘

        You can see the right face has been rotated but not the left face.
        """
        g = self.parse(mv)
        lst = self.facets(g)
        line1 = "             ┌──────────────┐\n"
        line2 = "             │%3d  %3d  %3d │\n" % (lst[0], lst[1], lst[2])
        line3 = "             │%3d   top %3d │\n" % (lst[3], lst[4])
        line4 = "             │%3d  %3d  %3d │\n" % (lst[5], lst[6], lst[7])
        line5 = "┌────────────┼──────────────┼─────────────┬────────────┐\n"
        line6 = "│%3d %3d %3d │%3d  %3d  %3d │%3d  %3d %3d │%3d %3d %3d │\n" % (lst[8], lst[9], lst[10], lst[16], lst[17], lst[18], lst[24], lst[25], lst[26], lst[32], lst[33], lst[34])
        line7 = "│%3d left%3d │%3d  front%3d │%3d right%3d │%3d rear%3d │\n" % (lst[11], lst[12], lst[19], lst[20], lst[27], lst[28], lst[35], lst[36])
        line8 = "│%3d %3d %3d │%3d  %3d  %3d │%3d  %3d %3d │%3d %3d %3d │\n" % (lst[13], lst[14], lst[15], lst[21], lst[22], lst[23], lst[29], lst[30], lst[31], lst[37], lst[38], lst[39])
        line9 = "└────────────┼──────────────┼─────────────┴────────────┘\n"
        line10 = "             │%3d  %3d  %3d │\n" % (lst[40], lst[41], lst[42])
        line11 = "             │%3d bottom%3d │\n" % (lst[43], lst[44])
        line12 = "             │%3d  %3d  %3d │\n" % (lst[45], lst[46], lst[47])
        line13 = "             └──────────────┘\n"
        return line1+line2+line3+line4+line5+line6+line7+line8+line9+line10+line11+line12+line13

    def plot_cube(self, mv, title=True, colors=[lpurple, yellow, red, green, orange, blue]):
        r"""
        Input the move mv, as a string in the Singmaster notation, and
        output the 2D plot of the cube in that state.

        Type ``P.show()`` to display any of the plots below.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: P = rubik.plot_cube("R^2*U^2*R^2*U^2*R^2*U^2", title=False)           # needs sage.plot
            sage: # (R^2U^2)^3  permutes 2 pairs of edges (uf,ub)(fr,br)
            sage: P = rubik.plot_cube("R*L*D^2*B^3*L^2*F^2*R^2*U^3*D*R^3*D^2*F^3*B^3*D^3*F^2*D^3*R^2*U^3*F^2*D^3")      # needs sage.plot
            sage: # the superflip (in 20f* moves)
            sage: P = rubik.plot_cube("U^2*F*U^2*L*R^(-1)*F^2*U*F^3*B^3*R*L*U^2*R*D^3*U*L^3*R*D*R^3*L^3*D^2")           # needs sage.plot
            sage: # "superflip+4 spot" (in 26q* moves)
        """
        g = self.parse(mv)
        state = self.facets(g)
        cubies = [create_poly(index2singmaster(state[x]), color_of_square(x+1, colors)) for x in range(48)]
        centers = [create_poly('%s_center' % "ulfrbd"[i], colors[i]) for i in range(6)]
        clrs = sum(cubies) + sum(centers)
        clrs.axes(show=False)
        if title:
            t = text('sagemath.org', (7.8, -3.5),rgbcolor=lgrey)
            P = clrs + t
            P.axes(show=False)
            return P
        return clrs

    def plot3d_cube(self, mv, title=True):
        r"""
        Displays `F,U,R` faces of the cube after the given move ``mv``. Mostly
        included for the purpose of drawing pictures and checking moves.

        INPUT:

        - ``mv`` -- string in the Singmaster notation
        - ``title`` -- boolean (default: ``True``); display the title information

        The first one below is "superflip+4 spot" (in 26q\* moves) and the
        second one is the superflip (in 20f\* moves). Type show(P) to view
        them.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: P = rubik.plot3d_cube("U^2*F*U^2*L*R^(-1)*F^2*U*F^3*B^3*R*L*U^2*R*D^3*U*L^3*R*D*R^3*L^3*D^2")         # needs sage.plot
            sage: P = rubik.plot3d_cube("R*L*D^2*B^3*L^2*F^2*R^2*U^3*D*R^3*D^2*F^3*B^3*D^3*F^2*D^3*R^2*U^3*F^2*D^3")    # needs sage.plot
        """
        g = self.parse(mv)
        state = self.facets(g)
        cubiesR = [plot3d_cubie(cubie_centers(c),cubie_colors(c,state)) for c in [32,31,30,29,28,27,26,25]]
        cubeR = sum(cubiesR)
        cubiesU = [plot3d_cubie(cubie_centers(c),cubie_colors(c,state)) for c in range(1,9)]
        cubeU = sum(cubiesU)
        cubiesF = [plot3d_cubie(cubie_centers(c),cubie_colors(c,state)) for c in [22,23,24,20,21]]
        cubeF = sum(cubiesF)
        centerR = polygon_plot3d([[1,3,1],[2,3,1],[2,3,2],[1,3,2],[1,3,1]],rgbcolor=green)
        centerF = polygon_plot3d([[3,1,1],[3,2,1],[3,2,2],[3,1,2],[3,1,1]],rgbcolor=red)
        centerU = polygon_plot3d([[1,1,3],[1,2,3],[2,2,3],[2,1,3],[1,1,3]],rgbcolor=lpurple)
        centers = centerF+centerR+centerU
        P = cubeR+cubeF+cubeU+centers
        P.axes(show=False)
        if title:
            t1 = text('Up, Front, and Right faces. '   , (-0.2, -2.5))
            t2 = text('      sagemath.org', (0.8, -3.1),rgbcolor=lgrey)
            t3 = text("     ", (3.5, 0), rgbcolor=white)
            P = P + t1 + t2 + t3
            P.axes(show=False)
            return P
        return P

    def legal(self, state, mode='quiet'):
        r"""
        Return 1 (true) if the dictionary ``state`` (in the
        same format as returned by the faces method) represents a legal
        position (or state) of the Rubik's cube or 0 (false)
        otherwise.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: r0 = rubik.faces("")
            sage: r1 = {'back': [[33, 34, 35], [36, 0, 37], [38, 39, 40]], 'down': [[41, 42, 43], [44, 0, 45], [46, 47, 48]],'front': [[17, 18, 19], [20, 0, 21], [22, 23, 24]],'left': [[9, 10, 11], [12, 0, 13], [14, 15, 16]],'right': [[25, 26, 27], [28, 0, 29], [30, 31, 32]],'up': [[1, 2, 3], [4, 0, 5], [6, 8, 7]]}
            sage: rubik.legal(r0)
            1
            sage: rubik.legal(r0,"verbose")
            (1, ())
            sage: rubik.legal(r1)
            0
        """
        try:
            g = self.parse(state)
            res = 1
        except ValueError:
            res = 0
            g = self.identity()
        if mode != 'quiet':
            return res, g
        else:
            return res

    def solve(self, state, algorithm='default'):
        r"""
        Solve the cube in the ``state``, given as a dictionary
        as in ``legal``. See the ``solve`` method
        of the RubiksCube class for more details.

        This may use GAP's ``EpimorphismFromFreeGroup`` and
        ``PreImagesRepresentative`` as explained below, if
        'gap' is passed in as the algorithm.

        This algorithm


        #. constructs the free group on 6 generators then computes a
           reasonable set of relations which they satisfy

        #. computes a homomorphism from the cube group to this free group
           quotient

        #. takes the cube position, regarded as a group element, and maps
           it over to the free group quotient

        #. using those relations and tricks from combinatorial group theory
           (stabilizer chains), solves the "word problem" for that element.

        #. uses python string parsing to rewrite that in cube notation.


        The Rubik's cube group has about `4.3 \times 10^{19}`
        elements, so this process is time-consuming. See
        https://www.gap-system.org/Doc/Examples/rubik.html for an
        interesting discussion of some GAP code analyzing the Rubik's
        cube.

        EXAMPLES::

            sage: rubik = CubeGroup()
            sage: state = rubik.faces("R")
            sage: rubik.solve(state)
            'R'
            sage: state = rubik.faces("R*U")
            sage: rubik.solve(state, algorithm='gap')       # long time
            'R*U'

        You can also check this another (but similar) way using the
        ``word_problem`` method (eg, G = rubik.group(); g =
        G("(3,38,43,19)(5,36,45,21)(8,33,48,24)(25,27,32,30)(26,29,31,28)");
        g.word_problem([b,d,f,l,r,u]), though the output will be less
        intuitive).
        """
        try:
            g = self.parse(state)
        except TypeError:
            return "Illegal or syntactically incorrect state. No solution."
        if algorithm != 'gap':
            C = RubiksCube(g)
            return C.solve(algorithm)

        group = libgap(self)
        hom = group.EpimorphismFromFreeGroup()
        sol = str(hom.PreImagesRepresentative(libgap(g)))
        for i, name in enumerate(self.gen_names()):
            sol = sol.replace("x%s" % (i + 1), name)
        return sol


##########################################################
#              3d object generation
##########################################################

def cubie_faces():
    """
    This provides a map from the 6 faces of the 27 cubies to the 48
    facets of the larger cube.

    -1,-1,-1 is left, top, front

    EXAMPLES::

        sage: from sage.groups.perm_gps.cubegroup import cubie_faces
        sage: sorted(cubie_faces().items())
        [((-1, -1, -1), [6, 17, 11, 0, 0, 0]),
         ((-1, -1, 0), [4, 0, 10, 0, 0, 0]),
         ((-1, -1, 1), [1, 0, 9, 0, 35, 0]),
         ((-1, 0, -1), [0, 20, 13, 0, 0, 0]),
         ((-1, 0, 0), [0, 0, -5, 0, 0, 0]),
         ((-1, 0, 1), [0, 0, 12, 0, 37, 0]),
         ((-1, 1, -1), [0, 22, 16, 41, 0, 0]),
         ((-1, 1, 0), [0, 0, 15, 44, 0, 0]),
         ((-1, 1, 1), [0, 0, 14, 46, 40, 0]),
         ((0, -1, -1), [7, 18, 0, 0, 0, 0]),
         ((0, -1, 0), [-6, 0, 0, 0, 0, 0]),
         ((0, -1, 1), [2, 0, 0, 0, 34, 0]),
         ((0, 0, -1), [0, -4, 0, 0, 0, 0]),
         ((0, 0, 0), [0, 0, 0, 0, 0, 0]),
         ((0, 0, 1), [0, 0, 0, 0, -2, 0]),
         ((0, 1, -1), [0, 23, 0, 42, 0, 0]),
         ((0, 1, 0), [0, 0, 0, -1, 0, 0]),
         ((0, 1, 1), [0, 0, 0, 47, 39, 0]),
         ((1, -1, -1), [8, 19, 0, 0, 0, 25]),
         ((1, -1, 0), [5, 0, 0, 0, 0, 26]),
         ((1, -1, 1), [3, 0, 0, 0, 33, 27]),
         ((1, 0, -1), [0, 21, 0, 0, 0, 28]),
         ((1, 0, 0), [0, 0, 0, 0, 0, -3]),
         ((1, 0, 1), [0, 0, 0, 0, 36, 29]),
         ((1, 1, -1), [0, 24, 0, 43, 0, 30]),
         ((1, 1, 0), [0, 0, 0, 45, 0, 31]),
         ((1, 1, 1), [0, 0, 0, 48, 38, 32])]
    """
    faceR = [[25,26,27], [28,-3,29], [30,31,32]] # green
    faceL = [[ 9,10,11], [12,-5,13], [14,15,16]] # orange
    faceU = [[ 1, 2, 3], [ 4,-6, 5], [ 6, 7, 8]] # red
    faceD = [[41,42,43], [44,-1,45], [46,47,48]] # purple
    faceF = [[17,18,19], [20,-4,21], [22,23,24]] # yellow
    faceB = [[33,34,35], [36,-2,37], [38,39,40]] # blue
    cubies = {}
    for x in [-1,0,1]:
        for y in [-1,0,1]:
            for z in [-1,0,1]:
                cubies[x,y,z] = [0,0,0,0,0,0]

    for i in [-1,0,1]:
        for j in [-1,0,1]:
            cubies[  i,  j, -1][1] = faceF[1+j][1+i]
            cubies[  i,  j,  1][4] = faceB[1+j][1-i]
            cubies[  i, -1,  j][0] = faceU[1-j][1+i]
            cubies[  i,  1,  j][3] = faceD[1+j][1+i]
            cubies[ -1,  i,  j][2] = faceL[1+i][1-j]
            cubies[  1,  i,  j][5] = faceR[1+i][1+j]

    return cubies


cubie_face_list = cubie_faces()


rand_colors = [(RDF.random_element(), RDF.random_element(), RDF.random_element()) for _ in range(56)]


@richcmp_method
class RubiksCube(SageObject):
    r"""
    The Rubik's cube (in a given state).

    EXAMPLES::

        sage: C = RubiksCube().move("R U R'")
        sage: C.show3d()                                                                # needs sage.plot

    ::

        sage: C = RubiksCube("R*L"); C
                     ┌──────────────┐
                     │ 17    2   38 │
                     │ 20   top  36 │
                     │ 22    7   33 │
        ┌────────────┼──────────────┼─────────────┬────────────┐
        │ 11  13  16 │ 41   18    3 │ 27   29  32 │ 48  34   6 │
        │ 10 left 15 │ 44  front  5 │ 26 right 31 │ 45 rear  4 │
        │  9  12  14 │ 46   23    8 │ 25   28  30 │ 43  39   1 │
        └────────────┼──────────────┼─────────────┴────────────┘
                     │ 40   42   19 │
                     │ 37 bottom 21 │
                     │ 35   47   24 │
                     └──────────────┘
        sage: C.show()                                                                  # needs sage.plot
        sage: C.solve(algorithm='gap')  # long time
        'L*R'
        sage: C == RubiksCube("L*R")
        True
    """
    def __init__(self, state=None, history=[], colors=[lpurple,yellow,red,green,orange,blue]):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = RubiksCube().move("R*U")
            sage: TestSuite(C).run()
        """
        self.colors = colors
        self._history = history
        self._group = CubeGroup()
        if state is None:
            self._state = self._group.identity()
        else:
            if isinstance(state, str):
                state = self._group.faces(state)
            if not isinstance(state, PermutationGroupElement):
                legal, state = self._group.legal(state, mode='gimme_group_element')
                if not legal:
                    raise ValueError("not a legal cube")
            self._state = state

    def move(self, g):
        r"""
        Move the Rubik's cube by ``g``.

        EXAMPLES::

            sage: RubiksCube().move("R*U") == RubiksCube("R*U")
            True
        """
        if not isinstance(g, self._group.element_class):
            g = self._group.move(g)[0]
        return RubiksCube(self._state * g, self._history + [g], self.colors)

    def undo(self):
        r"""
        Undo the last move of the Rubik's cube.

        EXAMPLES::

            sage: C = RubiksCube()
            sage: D = C.move("R*U")
            sage: D.undo() == C
            True
        """
        if not self._history:
            raise ValueError("no moves to undo")
        g = self._history[-1]
        return RubiksCube(self._state * ~g, self._history[:-1], self.colors)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RubiksCube().move("R*U")
                         ┌──────────────┐
                         │  3    5   38 │
                         │  2   top  36 │
                         │  1    4   25 │
            ┌────────────┼──────────────┼─────────────┬────────────┐
            │ 33  34  35 │  9   10    8 │ 19   29  32 │ 48  26  27 │
            │ 12 left 13 │ 20  front  7 │ 18 right 31 │ 45 rear 37 │
            │ 14  15  16 │ 22   23    6 │ 17   28  30 │ 43  39  40 │
            └────────────┼──────────────┼─────────────┴────────────┘
                         │ 41   42   11 │
                         │ 44 bottom 21 │
                         │ 46   47   24 │
                         └──────────────┘
        """
        return self._group.repr2d(self._state)

    def facets(self):
        r"""
        Return the facets of ``self``.

        EXAMPLES::

            sage: C = RubiksCube("R*U")
            sage: C.facets()
            [3, 5, 38, 2, 36, 1, 4, 25, 33, 34, 35, 12, 13, 14, 15, 16, 9, 10,
             8, 20, 7, 22, 23, 6, 19, 29, 32, 18, 31, 17, 28, 30, 48, 26, 27,
             45, 37, 43, 39, 40, 41, 42, 11, 44, 21, 46, 47, 24]
        """
        return self._group.facets(self._state)

    def plot(self):
        r"""
        Return a plot of ``self``.

        EXAMPLES::

            sage: C = RubiksCube("R*U")
            sage: C.plot()                                                              # needs sage.plot
            Graphics object consisting of 55 graphics primitives
        """
        return self._group.plot_cube(self._state)

    def show(self):
        r"""
        Show a plot of ``self``.

        EXAMPLES::

            sage: C = RubiksCube("R*U")
            sage: C.show()                                                              # needs sage.plot
        """
        self.plot().show()

    def cubie(self, size, gap, x, y, z, colors, stickers=True):
        """
        Return the cubie at `(x,y,z)`.

        INPUT:

        - ``size`` -- the size of the cubie
        - ``gap`` -- the gap between cubies
        - ``x``, ``y``, ``z`` -- the position of the cubie
        - ``colors`` -- the list of colors
        - ``stickers`` -- boolean (default: ``True``); whether to display stickers

        EXAMPLES::

            sage: C = RubiksCube("R*U")
            sage: C.cubie(0.15, 0.025, 0,0,0, C.colors*3)                               # needs sage.plot
            Graphics3d Object
        """
        sides = cubie_face_list[x, y, z]
        t = 2 * size + gap
        my_colors = [colors[sides[i]+6] for i in range(6)]
        if stickers:
            B = Box(size, size, size, color=(.1, .1, .1))
            S = B + B.stickers(my_colors, size*.1, size*.01)
            return S.translate(-t*x, -t*z, -t*y)
        else:
            return ColorCube(size, [colors[sides[i]+6] for i in range(6)]).translate(-t*x, -t*z, -t*y)

    def plot3d(self, stickers=True):
        r"""
        Return a 3D plot of ``self``.

        EXAMPLES::

            sage: C = RubiksCube("R*U")
            sage: C.plot3d()                                                            # needs sage.plot
            Graphics3d Object
        """
        while len(self.colors) < 7:
            self.colors.append((.1, .1, .1))
        side_colors = [Texture(color=c, ambient=.75) for c in self.colors]
        start_colors = sum([[c]*8 for c in side_colors], [])
        facets = self._group.facets(self._state)
        facet_colors = [0] * 48
        for i in range(48):
            facet_colors[facets[i]-1] = start_colors[i]
        all_colors = side_colors + facet_colors
        pm = [-1,0,1]
        C = sum([self.cubie(.15, .025, x, y, z, all_colors, stickers) for x in pm for y in pm for z in pm], Box(.35, .35, .35, color=self.colors[-1]))
        return C.rotateZ(1.5) #.scale([1,-1,1]).rotateZ(1.5)

    def show3d(self):
        r"""
        Show a 3D plot of ``self``.

        EXAMPLES::

            sage: C = RubiksCube("R*U")
            sage: C.show3d()                                                            # needs sage.plot
        """
        return self.plot3d().show()

    def __richcmp__(self, other, op):
        """
        Comparison.

        INPUT:

        - ``other`` -- anything

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: C = RubiksCube()
            sage: D = RubiksCube("R*U")
            sage: C < D
            True
            sage: C > D
            False
            sage: C == C
            True
            sage: C != D
            True
        """
        if not isinstance(other, RubiksCube):
            return NotImplemented
        return richcmp(self._state, other._state, op)

    def solve(self, algorithm='hybrid', timeout=15):
        r"""
        Solve the Rubik's cube.

        INPUT:

        - ``algorithm`` -- must be one of the following:

          - ``hybrid`` -- try ``kociemba`` for timeout seconds, then ``dietz``
          - ``kociemba`` -- use Dik T. Winter's program
            (reasonable speed, few moves)
          - ``dietz`` -- use Eric Dietz's cubex program
            (fast but lots of moves)
          - ``optimal`` -- use Michael Reid's optimal program
            (may take a long time)
          - ``gap`` -- use GAP word solution (can be slow)

        Any choice other than ``gap`` requires the optional package
        ``rubiks``. Otherwise, the ``gap`` algorithm is used.

        EXAMPLES::

            sage: C = RubiksCube("R U F L B D")
            sage: C.solve()           # optional - rubiks
            'R U F L B D'

        Dietz's program is much faster, but may give highly non-optimal
        solutions::

            sage: s = C.solve('dietz'); s   # optional - rubiks
            "U' L' L' U L U' L U D L L D' L' D L' D' L D L' U' L D' L' U L' B' U' L' U B L D L D' U' L' U L B L B' L' U L U' L' F' L' F L' F L F' L' D' L' D D L D' B L B' L B' L B F' L F F B' L F' B D' D' L D B' B' L' D' B U' U' L' B' D' F' F' L D F'"
            sage: C2 = RubiksCube(s)  # optional - rubiks
            sage: C == C2             # optional - rubiks
            True
        """
        from sage.features.rubiks import Rubiks
        if Rubiks().is_present():
            import sage.interfaces.rubik  # here to avoid circular referencing
        else:
            algorithm = 'gap'

        if algorithm == "default":
            algorithm = "hybrid"

        if algorithm == "hybrid":
            try:
                solver = sage.interfaces.rubik.DikSolver()
                return solver.solve(self.facets(), timeout=timeout)
            except RuntimeError:
                solver = sage.interfaces.rubik.CubexSolver()
                return solver.solve(self.facets())

        elif algorithm == "kociemba":
            solver = sage.interfaces.rubik.DikSolver()
            return solver.solve(self.facets(), timeout=timeout)

        elif algorithm == "dietz":
            solver = sage.interfaces.rubik.CubexSolver()
            return solver.solve(self.facets())

        elif algorithm == "optimal":
            # TODO: cache this, startup is expensive
            solver = sage.interfaces.rubik.OptimalSolver()
            return solver.solve(self.facets())

        elif algorithm == "gap":
            solver = CubeGroup()
            return solver.solve(self._state, algorithm='gap')

        else:
            raise ValueError(f"Unrecognized algorithm: {algorithm}")

    def scramble(self, moves=30):
        """
        Scramble the Rubik's cube.

        EXAMPLES::

            sage: C = RubiksCube()
            sage: C.scramble() # random
                         ┌──────────────┐
                         │  3    5   38 │
                         │  2   top  36 │
                         │  1    4   25 │
            ┌────────────┼──────────────┼─────────────┬────────────┐
            │ 33  34  35 │  9   10    8 │ 19   29  32 │ 48  26  27 │
            │ 12 left 13 │ 20  front  7 │ 18 right 31 │ 45 rear 37 │
            │ 14  15  16 │ 22   23    6 │ 17   28  30 │ 43  39  40 │
            └────────────┼──────────────┼─────────────┴────────────┘
                         │ 41   42   11 │
                         │ 44 bottom 21 │
                         │ 46   47   24 │
                         └──────────────┘
        """
        last_move = move = "  "
        all = []
        for i in range(moves):
            while move[0] == last_move[0]:
                move = "RLUDBF"[randint(0, 5)] + " '2"[randint(0, 2)]
            last_move = move
            all.append(move)
        return self.move(' '.join(all))
