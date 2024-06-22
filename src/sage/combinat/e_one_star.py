# sage.doctest: needs sage.combinat sage.modules
r"""
Substitutions over unit cube faces (Rauzy fractals)

This module implements the `E_1^*(\sigma)` substitution
associated with a one-dimensional substitution `\sigma`,
that acts on unit faces of dimension `(d-1)` in `\RR^d`.

This module defines the following classes and functions:

- ``Face`` -- a class to model a face

- ``Patch`` -- a class to model a finite set of faces

- ``E1Star`` -- a class to model the `E_1^*(\sigma)` application
  defined by the substitution sigma

See the documentation of these objects for more information.

The convention for the choice of the unit faces and the
definition of `E_1^*(\sigma)` varies from article to article.
Here, unit faces are defined by

.. MATH::

    \begin{array}{ccc}
    \,[x, 1]^* & = & \{x + \lambda e_2 + \mu e_3 : \lambda, \mu \in [0,1]\} \\
    \,[x, 2]^* & = & \{x + \lambda e_1 + \mu e_3 : \lambda, \mu \in [0,1]\} \\
    \,[x, 3]^* & = & \{x + \lambda e_1 + \mu e_2 : \lambda, \mu \in [0,1]\}
    \end{array}

and the dual substitution `E_1^*(\sigma)` is defined by

.. MATH::

    E_1^*(\sigma)([x,i]^*) =
    \bigcup_{k = 1,2,3} \; \bigcup_{s | \sigma(k) = pis}
    [M^{-1}(x + \ell(s)), k]^*,

where `\ell(s)` is the abelianized of `s`, and `M` is the matrix of `\sigma`.

AUTHORS:

- Franco Saliola (2009): initial version
- Vincent Delecroix, Timo Jolivet, Stepan Starosta, Sebastien Labbe (2010-05): redesign
- Timo Jolivet (2010-08, 2010-09, 2011): redesign

REFERENCES:

.. [AI] \P. Arnoux, S. Ito,
   Pisot substitutions and Rauzy fractals,
   Bull. Belg. Math. Soc. 8 (2), 2001, pp. 181--207

.. [SAI] \Y. Sano, P. Arnoux, S. Ito,
   Higher dimensional extensions of substitutions and their dual maps,
   J. Anal. Math. 83, 2001, pp. 183--206

EXAMPLES:

We start by drawing a simple three-face patch::

    sage: from sage.combinat.e_one_star import E1Star, Face, Patch
    sage: x = [Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)]
    sage: P = Patch(x)
    sage: P
    Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
    sage: P.plot()                   #not tested

We apply a substitution to this patch, and draw the result::

    sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
    sage: E = E1Star(sigma)
    sage: E(P)
    Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(0, 1, -1), 2]*, [(1, 0, -1), 1]*]
    sage: E(P).plot()                #not tested

.. NOTE::

    - The type of a face is given by an integer in ``[1, ..., d]``
      where ``d`` is the length of the vector of the face.

    - The alphabet of the domain and the codomain of `\sigma` must be
      equal, and they must be of the form ``[1, ..., d]``, where ``d``
      is a positive integer corresponding to the length of the vectors
      of the faces on which `E_1^*(\sigma)` will act.

::

    sage: P = Patch([Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)])
    sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
    sage: E = E1Star(sigma)
    sage: E(P)
    Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(0, 1, -1), 2]*, [(1, 0, -1), 1]*]

The application of an ``E1Star`` substitution assigns to each new face the color of its preimage.
The ``repaint`` method allows us to repaint the faces of a patch.
A single color can also be assigned to every face, by specifying a list of a single color::

    sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: P = E(P, 5)
    sage: P.repaint(['green'])
    sage: P.plot()                   #not tested

A list of colors allows us to color the faces sequentially::

    sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: P = E(P)
    sage: P.repaint(['red', 'yellow', 'green', 'blue', 'black'])
    sage: P = E(P, 3)
    sage: P.plot()                   #not tested

All the color schemes from ``list(matplotlib.cm.datad)`` can be used::

    sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: P.repaint(cmap='summer')
    sage: P = E(P, 3)
    sage: P.plot()                   #not tested
    sage: P.repaint(cmap='hsv')
    sage: P = E(P, 2)
    sage: P.plot()                   #not tested

It is also possible to specify a dictionary to color the faces according to their type::

    sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: P = E(P, 5)
    sage: P.repaint({1:(0.7, 0.7, 0.7), 2:(0.5,0.5,0.5), 3:(0.3,0.3,0.3)})
    sage: P.plot()                   #not tested
    sage: P.repaint({1:'red', 2:'yellow', 3:'green'})
    sage: P.plot()                   #not tested

Let us look at a nice big patch in 3D::

    sage: sigma = WordMorphism({1:[1,2], 2:[3], 3:[1]})
    sage: E = E1Star(sigma)
    sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: P = P + P.translate([-1,1,0])
    sage: P = E(P, 11)
    sage: P.plot3d()                 #not tested

Plotting with TikZ pictures is possible::

    sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: s = P.plot_tikz()
    sage: print(s)                    #not tested
    \begin{tikzpicture}
    [x={(-0.216506cm,-0.125000cm)}, y={(0.216506cm,-0.125000cm)}, z={(0.000000cm,0.250000cm)}]
    \definecolor{facecolor}{rgb}{0.000,1.000,0.000}
    \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
    (0, 0, 0) -- (0, 0, 1) -- (1, 0, 1) -- (1, 0, 0) -- cycle;
    \definecolor{facecolor}{rgb}{1.000,0.000,0.000}
    \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
    (0, 0, 0) -- (0, 1, 0) -- (0, 1, 1) -- (0, 0, 1) -- cycle;
    \definecolor{facecolor}{rgb}{0.000,0.000,1.000}
    \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
    (0, 0, 0) -- (1, 0, 0) -- (1, 1, 0) -- (0, 1, 0) -- cycle;
    \end{tikzpicture}

Plotting patches made of unit segments instead of unit faces::

    sage: P = Patch([Face([0,0], 1), Face([0,0], 2)])
    sage: E = E1Star(WordMorphism({1:[1,2],2:[1]}))
    sage: F = E1Star(WordMorphism({1:[1,1,2],2:[2,1]}))
    sage: E(P,5).plot()                                                                 # needs sage.plot
    Graphics object consisting of 21 graphics primitives
    sage: F(P,3).plot()                                                                 # needs sage.plot
    Graphics object consisting of 34 graphics primitives

Everything works in any dimension (except for the plotting features
which only work in dimension two or three)::

    sage: P = Patch([Face((0,0,0,0),1), Face((0,0,0,0),4)])
    sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1,4], 4:[1]})
    sage: E = E1Star(sigma)
    sage: E(P)
    Patch: [[(0, 0, 0, 0), 3]*, [(0, 0, 0, 0), 4]*, [(0, 0, 1, -1), 3]*, [(0, 1, 0, -1), 2]*, [(1, 0, 0, -1), 1]*]

::

    sage: sigma = WordMorphism({1:[1,2],2:[1,3],3:[1,4],4:[1,5],5:[1,6],6:[1,7],7:[1,8],8:[1,9],9:[1,10],10:[1,11],11:[1,12],12:[1]})
    sage: E = E1Star(sigma)
    sage: E
    E_1^*(1->12, 10->1,11, 11->1,12, 12->1, 2->13, 3->14, 4->15, 5->16, 6->17, 7->18, 8->19, 9->1,10)
    sage: P = Patch([Face((0,0,0,0,0,0,0,0,0,0,0,0),t) for t in [1,2,3]])
    sage: for x in sorted(E(P), key=lambda x : (x.vector(),x.type())): print(x)
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 1]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 2]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 12]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1), 11]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1), 10]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1), 9]*
    [(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1), 8]*
    [(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1), 7]*
    [(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1), 6]*
    [(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1), 5]*
    [(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1), 4]*
    [(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1), 3]*
    [(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1), 2]*
    [(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1), 1]*
"""
# ****************************************************************************
#       Copyright (C) 2010 Franco Saliola <saliola@gmail.com>
#                          Vincent Delecroix <20100.delecroix@gmail.com>
#                          Timo Jolivet <timo.jolivet@gmail.com>
#                          Stepan Starosta <stepan.starosta@gmail.com>
#                          Sébastien Labbé <slabqc at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations

from sage.misc.functional import det
from sage.structure.sage_object import SageObject
from sage.combinat.words.morphism import WordMorphism
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.all", "Graphics")
lazy_import("sage.plot.colors", "Color")
lazy_import("sage.plot.polygon", "polygon")
lazy_import("sage.plot.line", "line")
from sage.rings.integer_ring import ZZ
from sage.misc.latex import LatexExpr
from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import richcmp_by_eq_and_lt, richcmp_method

lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.modules.free_module_element', 'vector')


# matplotlib color maps, loaded on-demand
cm = None


@richcmp_method
class Face(SageObject):
    r"""
    A class to model a unit face of arbitrary dimension.

    A unit face in dimension `d` is represented by
    a `d`-dimensional vector ``v`` and a type ``t`` in `\{1, \ldots, d\}`.
    The type of the face corresponds to the canonical unit vector
    to which the face is orthogonal.
    The optional ``color`` argument is used in plotting functions.

    INPUT:

    - ``v`` -- tuple of integers
    - ``t`` -- integer in ``[1, ..., len(v)]``, type of the face. The face of type `i`
      is orthogonal to the canonical vector `e_i`.
    - ``color`` -- color (default: ``None``); color of the face,
      used for plotting only. If ``None``, its value is guessed from the
      face type.

    EXAMPLES::

        sage: from sage.combinat.e_one_star import Face
        sage: f = Face((0,2,0), 3)
        sage: f.vector()
        (0, 2, 0)
        sage: f.type()
        3

    ::

        sage: f = Face((0,2,0), 3, color=(0.5, 0.5, 0.5))
        sage: f.color()
        RGB color (0.5, 0.5, 0.5)
    """

    def __init__(self, v, t, color=None):
        r"""
        Face constructor. See class doc for more information.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,2,0), 3)
            sage: f.vector()
            (0, 2, 0)
            sage: f.type()
            3

        TESTS:

        We test that types can be given by an int (see :issue:`10699`)::

            sage: f = Face((0,2,0), int(1))
        """
        self._vector = (ZZ**len(v))(v)
        self._vector.set_immutable()

        if not ((t in ZZ) and 1 <= t <= len(v)):
            raise ValueError('the type must be an integer between 1 and len(v)')
        self._type = t

        if color is None:
            if self._type == 1:
                color = Color((1, 0, 0))
            elif self._type == 2:
                color = Color((0, 1, 0))
            elif self._type == 3:
                color = Color((0, 0, 1))
            else:
                color = Color()
        self._color = Color(color)

    def __repr__(self) -> str:
        r"""
        String representation of a face.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,0,0,3), 3)
            sage: f
            [(0, 0, 0, 3), 3]*

        ::

            sage: f = Face((0,0,0,3), 3)
            sage: f
            [(0, 0, 0, 3), 3]*
        """
        return "[%s, %s]*" % (self.vector(), self.type())

    __richcmp__ = richcmp_by_eq_and_lt('_eq', '_lt')

    def _eq(self, other) -> bool:
        r"""
        Equality of faces.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,0,0,3), 3)
            sage: g = Face((0,0,0,3), 3)
            sage: f == g
            True
        """
        return (isinstance(other, Face) and
                self.vector() == other.vector() and
                self.type() == other.type())

    def _lt(self, other) -> bool:
        r"""
        Compare ``self`` and ``other``.

        The vectors of the faces are first compared,
        and the types of the faces are compared if the vectors are equal.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: Face([-2,1,0], 2) < Face([-1,2,2],3)
            True
            sage: Face([-2,1,0], 2) < Face([-2,1,0],3)
            True
            sage: Face([-2,1,0], 2) < Face([-2,1,0],2)
            False
        """
        if self.vector() < other.vector():
            return True
        if self.vector() == other.vector():
            return self.type() < other.type()
        return False

    def __hash__(self) -> int:
        r"""
        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,0,0,3), 3)
            sage: g = Face((0,0,0,3), 3)
            sage: hash(f) == hash(g)
            True
        """
        return hash((self.vector(), self.type()))

    def __add__(self, other):
        r"""
        Addition of ``self`` with a Face, a Patch or a finite iterable of faces.

        INPUT:

        - ``other`` -- a Patch or a Face or a finite iterable of faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: f = Face([0,0,0], 3)
            sage: g = Face([0,1,-1], 2)
            sage: f + g
            Patch: [[(0, 0, 0), 3]*, [(0, 1, -1), 2]*]
            sage: P = Patch([Face([0,0,0], 1), Face([0,0,0], 2)])
            sage: f + P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]

        Adding a finite iterable of faces::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face([0,0,0], 3)
            sage: f + [f,f]
            Patch: [[(0, 0, 0), 3]*]
        """
        if isinstance(other, Face):
            return Patch([self, other])
        else:
            return Patch(other).union(self)

    def vector(self):
        r"""
        Return the vector of the face.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,2,0), 3)
            sage: f.vector()
            (0, 2, 0)
        """
        return self._vector

    def type(self):
        r"""
        Return the type of the face.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,2,0), 3)
            sage: f.type()
            3

        ::

            sage: f = Face((0,2,0), 3)
            sage: f.type()
            3
        """
        return self._type

    def color(self, color=None):
        r"""
        Return or change the color of the face.

        INPUT:

        - ``color`` -- string, rgb tuple, color (default: ``None``)
          the new color to assign to the face. If ``None``, it returns the
          color of the face.

        OUTPUT: color or None

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,2,0), 3)
            sage: f.color()
            RGB color (0.0, 0.0, 1.0)
            sage: f.color('red')
            sage: f.color()
            RGB color (1.0, 0.0, 0.0)
        """
        if color is not None:
            self._color = Color(color)
        else:
            return self._color

    def _plot(self, projmat, face_contour, opacity) -> Graphics:
        r"""
        Return a 2D graphic object representing the face.

        INPUT:

        - ``projmat`` -- 2*3 projection matrix (used only for faces in three dimensions)
        - ``face_contour`` -- dict, maps the face type to vectors describing
          the contour of unit faces (used only for faces in three dimensions)
        - ``opacity`` -- the alpha value for the color of the face

        OUTPUT: 2D graphic object

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,0,3), 3)
            sage: projmat = matrix(2, [-1.7320508075688772*0.5, 1.7320508075688772*0.5, 0, -0.5, -0.5, 1])
            sage: face_contour = {}
            sage: face_contour[1] = map(vector, [(0,0,0),(0,1,0),(0,1,1),(0,0,1)])
            sage: face_contour[2] = map(vector, [(0,0,0),(0,0,1),(1,0,1),(1,0,0)])
            sage: face_contour[3] = map(vector, [(0,0,0),(1,0,0),(1,1,0),(0,1,0)])
            sage: G = f._plot(projmat, face_contour, 0.75)                              # needs sage.plot

        ::

            sage: f = Face((0,0), 2)
            sage: f._plot(None, None, 1)                                                # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        v = self.vector()
        t = self.type()
        G = Graphics()

        if len(v) == 2:
            if t == 1:
                G += line([v, v + vector([0, 1])], rgbcolor=self.color(), thickness=1.5, alpha=opacity)
            elif t == 2:
                G += line([v, v + vector([1, 0])], rgbcolor=self.color(), thickness=1.5, alpha=opacity)

        elif len(v) == 3:
            G += polygon([projmat * (u + v)
                          for u in face_contour[t]], alpha=opacity,
                         thickness=1, rgbcolor=self.color())

        else:
            raise NotImplementedError("plotting is implemented only for patches in two or three dimensions.")

        return G

    def _plot3d(self, face_contour):
        r"""
        3D representation of a unit face (Jmol).

        INPUT:

        - ``face_contour`` -- dict, maps the face type to vectors describing
          the contour of unit faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face
            sage: f = Face((0,0,3), 3)
            sage: face_contour = {1: map(vector, [(0,0,0),(0,1,0),(0,1,1),(0,0,1)]), 2: map(vector, [(0,0,0),(0,0,1),(1,0,1),(1,0,0)]), 3: map(vector, [(0,0,0),(1,0,0),(1,1,0),(0,1,0)])}
            sage: G = f._plot3d(face_contour)      #not tested
        """
        v = self.vector()
        t = self.type()
        c = self.color()
        G = polygon([u + v for u in face_contour[t]], rgbcolor=c)
        return G


class Patch(SageObject):
    r"""
    A class to model a collection of faces. A patch is represented by an immutable set of Faces.

    .. NOTE::

        The dimension of a patch is the length of the vectors of the faces in the patch,
        which is assumed to be the same for every face in the patch.

    .. NOTE::

        Since version 4.7.1, Patches are immutable, except for the colors of the faces,
        which are not taken into account for equality tests and hash functions.

    INPUT:

    - ``faces`` -- finite iterable of faces
    - ``face_contour`` -- dictionary (default: ``None``); maps the face
      type to vectors describing the contour of unit faces. If ``None``,
      defaults contour are assumed for faces of type 1, 2, 3 or 1, 2, 3.
      Used in plotting methods only.

    EXAMPLES::

        sage: from sage.combinat.e_one_star import Face, Patch
        sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
        sage: P
        Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]

    ::

        sage: face_contour = {}
        sage: face_contour[1] = map(vector, [(0,0,0),(0,1,0),(0,1,1),(0,0,1)])
        sage: face_contour[2] = map(vector, [(0,0,0),(0,0,1),(1,0,1),(1,0,0)])
        sage: face_contour[3] = map(vector, [(0,0,0),(1,0,0),(1,1,0),(0,1,0)])
        sage: Patch([Face((0,0,0),t) for t in [1,2,3]], face_contour=face_contour)
        Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
    """

    def __init__(self, faces, face_contour=None):
        r"""
        Constructor of a patch (set of faces).

        See class doc for more information.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]

        TESTS:

        We test that colors are not anymore mixed up between
        Patches (see :issue:`11255`)::

            sage: P = Patch([Face([0,0,0],2)])
            sage: Q = Patch(P)
            sage: next(iter(P)).color()
            RGB color (0.0, 1.0, 0.0)
            sage: next(iter(Q)).color('yellow')
            sage: next(iter(P)).color()
            RGB color (0.0, 1.0, 0.0)
        """
        self._faces = frozenset(Face(f.vector(), f.type(), f.color()) for f in faces)

        try:
            f0 = next(iter(self._faces))
        except StopIteration:
            self._dimension = None
        else:
            self._dimension = len(f0.vector())

        if face_contour is not None:
            self._face_contour = face_contour

        else:
            self._face_contour = {
                1: [vector(t) for t in [(0, 0, 0), (0, 1, 0),
                                        (0, 1, 1), (0, 0, 1)]],
                2: [vector(t) for t in [(0, 0, 0), (0, 0, 1),
                                        (1, 0, 1), (1, 0, 0)]],
                3: [vector(t) for t in [(0, 0, 0), (1, 0, 0),
                                        (1, 1, 0), (0, 1, 0)]]
            }

    def __eq__(self, other) -> bool:
        r"""
        Equality test for Patch.

        INPUT:

        - ``other`` -- an object

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)])
            sage: Q = Patch([Face((0,1,0),1), Face((0,0,0),3)])
            sage: P == P
            True
            sage: P == Q
            False
            sage: P == 4
            False

        ::

            sage: s = WordMorphism({1:[1,3], 2:[1,2,3], 3:[3]})
            sage: t = WordMorphism({1:[1,2,3], 2:[2,3], 3:[3]})
            sage: P = Patch([Face((0,0,0), 1), Face((0,0,0), 2), Face((0,0,0), 3)])
            sage: E1Star(s)(P) == E1Star(t)(P)
            False
            sage: E1Star(s*t)(P) == E1Star(t)(E1Star(s)(P))
            True
        """
        return (isinstance(other, Patch) and self._faces == other._faces)

    def __hash__(self) -> int:
        r"""
        Hash function of Patch.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: hash(P)      #random
            -4839605361791007520

        TESTS:

        We test that two equal patches have the same hash (see :issue:`11255`)::

            sage: P = Patch([Face([0,0,0],1), Face([0,0,0],2)])
            sage: Q = Patch([Face([0,0,0],2), Face([0,0,0],1)])
            sage: P == Q
            True
            sage: hash(P) == hash(Q)
            True

        Changing the color does not affect the hash value::

            sage: p = Patch([Face((0,0,0), t) for t in [1,2,3]])
            sage: H1 = hash(p)
            sage: p.repaint(['blue'])
            sage: H2 = hash(p)
            sage: H1 == H2
            True
        """
        return hash(self._faces)

    def __len__(self) -> int:
        r"""
        Return the number of faces contained in the patch.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: len(P)       #indirect doctest
            3
        """
        return len(self._faces)

    def __iter__(self):
        r"""
        Return an iterator over the faces of the patch.

        OUTPUT: iterator

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: it = iter(P)
            sage: type(next(it))
            <class 'sage.combinat.e_one_star.Face'>
            sage: type(next(it))
            <class 'sage.combinat.e_one_star.Face'>
            sage: type(next(it))
            <class 'sage.combinat.e_one_star.Face'>
            sage: type(next(it))
            Traceback (most recent call last):
            ...
            StopIteration
        """
        return iter(self._faces)

    def __add__(self, other):
        r"""
        Addition of patches (union).

        INPUT:

        - ``other`` -- a Patch or a Face or a finite iterable of faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face([0,0,0], 1), Face([0,0,0], 2)])
            sage: Q = P.translate([1,-1,0])
            sage: P + Q
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(1, -1, 0), 1]*, [(1, -1, 0), 2]*]
            sage: P + Face([0,0,0],3)
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
            sage: P + [Face([0,0,0],3), Face([1,1,1],2)]
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(1, 1, 1), 2]*]
        """
        return self.union(other)

    def __sub__(self, other):
        r"""
        Subtraction of patches (difference).

        INPUT:

        - ``other`` -- a Patch or a Face or a finite iterable of faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P - Face([0,0,0],2)
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 3]*]
            sage: P - P
            Patch: []
        """
        return self.difference(other)

    def __repr__(self) -> str:
        r"""
        String representation of a patch.

        Displays all the faces if there less than 20,
        otherwise displays only the number of faces.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]

        ::

            sage: x = [Face((0,0,a),1) for a in range(25)]
            sage: P = Patch(x)
            sage: P
            Patch of 25 faces
        """
        if len(self) <= 20:
            L = list(self)
            L.sort(key=lambda x: (x.vector(), x.type()))
            return "Patch: %s" % L
        else:
            return "Patch of %s faces" % len(self)

    def union(self, other) -> Patch:
        r"""
        Return a Patch consisting of the union of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a Patch or a Face or a finite iterable of faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),1), Face((0,0,0),2)])
            sage: P.union(Face((1,2,3), 3))
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(1, 2, 3), 3]*]
            sage: P.union([Face((1,2,3), 3), Face((2,3,3), 2)])
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(1, 2, 3), 3]*, [(2, 3, 3), 2]*]
        """
        if isinstance(other, Face):
            return Patch(self._faces.union([other]))
        else:
            return Patch(self._faces.union(other))

    def difference(self, other) -> Patch:
        r"""
        Return the difference of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a finite iterable of faces or a single face

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.difference(Face([0,0,0],2))
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 3]*]
            sage: P.difference(P)
            Patch: []
        """
        if isinstance(other, Face):
            return Patch(self._faces.difference([other]))
        else:
            return Patch(self._faces.difference(other))

    def dimension(self) -> None | int:
        r"""
        Return the dimension of the vectors of the faces of ``self``.

        It returns ``None`` if ``self`` is the empty patch.

        The dimension of a patch is the length of the vectors of the faces in the patch,
        which is assumed to be the same for every face in the patch.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.dimension()
            3

        TESTS::

            sage: from sage.combinat.e_one_star import Patch
            sage: p = Patch([])
            sage: p.dimension() is None
            True

        It works when the patch is created from an iterator::

            sage: p = Patch(Face((0,0,0),t) for t in [1,2,3])
            sage: p.dimension()
            3
        """
        return self._dimension

    def faces_of_vector(self, v) -> list[Face]:
        r"""
        Return a list of the faces whose vector is ``v``.

        INPUT:

        - ``v`` -- a vector

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),1), Face((1,2,0),3), Face((1,2,0),1)])
            sage: sorted(P.faces_of_vector([1,2,0]))
            [[(1, 2, 0), 1]*, [(1, 2, 0), 3]*]
        """
        v = vector(v)
        return [f for f in self if f.vector() == v]

    def faces_of_type(self, t) -> list[Face]:
        r"""
        Return a list of the faces that have type ``t``.

        INPUT:

        - ``t`` -- integer or any other type

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),1), Face((1,2,0),3), Face((1,2,0),1)])
            sage: sorted(P.faces_of_type(1))
            [[(0, 0, 0), 1]*, [(1, 2, 0), 1]*]
        """
        return [f for f in self if f.type() == t]

    def faces_of_color(self, color) -> list[Face]:
        r"""
        Return a list of the faces that have the given color.

        INPUT:

        - ``color`` -- color

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),1, 'red'), Face((1,2,0),3, 'blue'), Face((1,2,0),1, 'red')])
            sage: sorted(P.faces_of_color('red'))
            [[(0, 0, 0), 1]*, [(1, 2, 0), 1]*]
        """
        color = tuple(Color(color))
        return [f for f in self if tuple(f.color()) == color]

    def translate(self, v) -> Patch:
        r"""
        Return a translated copy of ``self`` by vector `v`.

        INPUT:

        - ``v`` -- vector or tuple

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: P = Patch([Face((0,0,0),1), Face((1,2,0),3), Face((1,2,0),1)])
            sage: P.translate([-1,-2,0])
            Patch: [[(-1, -2, 0), 1]*, [(0, 0, 0), 1]*, [(0, 0, 0), 3]*]
        """
        v = vector(v)
        return Patch(Face(f.vector() + v, f.type(), f.color()) for f in self)

    def occurrences_of(self, other) -> list:
        r"""
        Return all positions at which other appears in self, that is,
        all vectors v such that ``set(other.translate(v)) <= set(self)``.

        INPUT:

        - ``other`` -- a Patch

        OUTPUT: list of vectors

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch, E1Star
            sage: P = Patch([Face([0,0,0], 1), Face([0,0,0], 2), Face([0,0,0], 3)])
            sage: Q = Patch([Face([0,0,0], 1), Face([0,0,0], 2)])
            sage: P.occurrences_of(Q)
            [(0, 0, 0)]
            sage: Q = Q.translate([1,2,3])
            sage: P.occurrences_of(Q)
            [(-1, -2, -3)]

        ::

            sage: E = E1Star(WordMorphism({1:[1,2], 2:[1,3], 3:[1]}))
            sage: P = Patch([Face([0,0,0], 1), Face([0,0,0], 2), Face([0,0,0], 3)])
            sage: P = E(P,4)
            sage: Q = Patch([Face([0,0,0], 1), Face([0,0,0], 2)])
            sage: L = P.occurrences_of(Q)
            sage: sorted(L)
            [(0, 0, 0), (0, 0, 1), (0, 1, -1), (1, 0, -1), (1, 1, -3), (1, 1, -2)]
        """
        f0 = next(iter(other))
        x = f0.vector()
        t = f0.type()
        L = self.faces_of_type(t)
        positions = []
        for f in L:
            y = f.vector()
            if other.translate(y - x)._faces.issubset(self._faces):
                positions.append(y - x)
        return positions

    def repaint(self, cmap='Set1') -> None:
        r"""
        Repaint all the faces of ``self`` from the given color map.

        This only changes the colors of the faces of ``self``.

        INPUT:

        - ``cmap`` -- color map (default: ``'Set1'``). It can be one of the
          following:

           - ``string`` -- a coloring map; for available coloring map names type:
             ``sorted(colormaps)``
           - ``list`` -- list of colors to assign cyclically to the faces
             A list of a single color colors all the faces with the same color
           - ``dict`` -- dictionary of face types mapped to colors, to color the
             faces according to their type
           - ``{}``, the empty dict -- shortcut for
             ``{1:'red', 2:'green', 3:'blue'}``

        EXAMPLES:

        Using a color map::

            sage: from sage.combinat.e_one_star import Face, Patch
            sage: color = (0, 0, 0)
            sage: P = Patch([Face((0,0,0),t,color) for t in [1,2,3]])
            sage: for f in P: f.color()
            RGB color (0.0, 0.0, 0.0)
            RGB color (0.0, 0.0, 0.0)
            RGB color (0.0, 0.0, 0.0)
            sage: P.repaint()
            sage: next(iter(P)).color()    #random
            RGB color (0.498..., 0.432..., 0.522...)

        Using a list of colors::

            sage: P = Patch([Face((0,0,0),t,color) for t in [1,2,3]])
            sage: P.repaint([(0.9, 0.9, 0.9), (0.65,0.65,0.65), (0.4,0.4,0.4)])
            sage: for f in P: f.color()
            RGB color (0.9, 0.9, 0.9)
            RGB color (0.65, 0.65, 0.65)
            RGB color (0.4, 0.4, 0.4)

        Using a dictionary to color faces according to their type::

            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.repaint({1:'black', 2:'yellow', 3:'green'})
            sage: P.plot()                   #not tested
            sage: P.repaint({})
            sage: P.plot()                   #not tested
        """
        if cmap == {}:
            cmap = {1: 'red', 2: 'green', 3: 'blue'}

        if isinstance(cmap, dict):
            for f in self:
                f.color(cmap[f.type()])

        elif isinstance(cmap, list):
            L = len(cmap)
            for i, f in enumerate(self):
                f.color(cmap[i % L])

        elif isinstance(cmap, str):
            # matplotlib color maps
            global cm
            if cm is None:
                from matplotlib import cm
            assert cm is not None
            if cmap not in cm.datad:
                raise RuntimeError("color map %s not known (type sorted(colors) for valid names)" % cmap)
            cmap = cm.__dict__[cmap]
            dim = float(len(self))
            for i, f in enumerate(self):
                f.color(cmap(i / dim)[:3])

        else:
            raise TypeError("type of cmap (=%s) must be dict, list or str" % cmap)

    def plot(self, projmat=None, opacity=0.75) -> Graphics:
        r"""
        Return a 2D graphic object depicting the patch.

        INPUT:

        - ``projmat`` -- matrix (default: ``None``); the projection
          matrix. Its number of lines must be two. Its number of columns
          must equal the dimension of the ambient space of the faces. If
          ``None``, the isometric projection is used by default.

        - ``opacity`` -- float between ``0`` and ``1`` (default: ``0.75``)
          opacity of the face

        .. WARNING::

            Plotting is implemented only for patches in two or three dimensions.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.plot()                                                              # needs sage.plot
            Graphics object consisting of 3 graphics primitives

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P = E(P, 5)
            sage: P.plot()                                                              # needs sage.plot
            Graphics object consisting of 57 graphics primitives

        Plot with a different projection matrix::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: M = matrix(2, 3, [1,0,-1,0.3,1,-3])
            sage: P = E(P, 3)
            sage: P.plot(projmat=M)                                                     # needs sage.plot
            Graphics object consisting of 17 graphics primitives

        Plot patches made of unit segments::

            sage: P = Patch([Face([0,0], 1), Face([0,0], 2)])
            sage: E = E1Star(WordMorphism({1:[1,2],2:[1]}))
            sage: F = E1Star(WordMorphism({1:[1,1,2],2:[2,1]}))
            sage: E(P,5).plot()                                                         # needs sage.plot
            Graphics object consisting of 21 graphics primitives
            sage: F(P,3).plot()                                                         # needs sage.plot
            Graphics object consisting of 34 graphics primitives
        """
        if self.dimension() == 2:
            G = Graphics()
            for face in self:
                G += face._plot(None, None, 1)
            G.set_aspect_ratio(1)
            return G

        if self.dimension() == 3:
            if projmat is None:
                projmat = matrix(2, [-1.7320508075688772 * 0.5,
                                     1.7320508075688772 * 0.5,
                                     0, -0.5, -0.5, 1])

            G = Graphics()
            for face in self:
                G += face._plot(projmat, self._face_contour, opacity)
            G.set_aspect_ratio(1)
            return G

        else:
            raise NotImplementedError("plotting is implemented only for patches in two or three dimensions.")

    def plot3d(self):
        r"""
        Return a 3D graphics object depicting the patch.

        .. WARNING::

            3D plotting is implemented only for patches in three dimensions.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.plot3d()                #not tested

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P = E(P, 5)
            sage: P.repaint()
            sage: P.plot3d()                #not tested
        """
        if self.dimension() != 3:
            raise NotImplementedError("3D plotting is implemented only for patches in three dimensions")

        face_list = [face._plot3d(self._face_contour) for face in self]
        G = sum(face_list)
        return G

    def plot_tikz(self, projmat=None, print_tikz_env=True, edgecolor='black',
            scale=0.25, drawzero=False, extra_code_before='', extra_code_after='') -> str:
        r"""
        Return a string containing some TikZ code to be included into
        a LaTeX document, depicting the patch.

        .. WARNING::

            Tikz Plotting is implemented only for patches in three dimensions.

        INPUT:

        - ``projmat`` -- matrix (default: ``None``); the projection
          matrix. Its number of lines must be two. Its number of columns
          must equal the dimension of the ambient space of the faces. If
          ``None``, the isometric projection is used by default.
        - ``print_tikz_env`` -- boolean (default: ``True``); if ``True``,
          the tikzpicture environment are printed
        - ``edgecolor`` -- string (default: ``'black'``); either
          ``'black'`` or ``'facecolor'`` (color of unit face edges)
        - ``scale`` -- real number (default: ``0.25``) scaling
          constant for the whole figure
        - ``drawzero`` -- boolean (default: ``False``); if ``True``,
          mark the origin by a black dot
        - ``extra_code_before`` -- string (default: ``''``); extra code to
          include in the tikz picture
        - ``extra_code_after`` -- string (default: ``''``); extra code to
          include in the tikz picture

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: s = P.plot_tikz()
            sage: len(s)
            602
            sage: print(s)       #not tested
            \begin{tikzpicture}
            [x={(-0.216506cm,-0.125000cm)}, y={(0.216506cm,-0.125000cm)}, z={(0.000000cm,0.250000cm)}]
            \definecolor{facecolor}{rgb}{0.000,1.000,0.000}
            \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
            (0, 0, 0) -- (0, 0, 1) -- (1, 0, 1) -- (1, 0, 0) -- cycle;
            \definecolor{facecolor}{rgb}{1.000,0.000,0.000}
            \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
            (0, 0, 0) -- (0, 1, 0) -- (0, 1, 1) -- (0, 0, 1) -- cycle;
            \definecolor{facecolor}{rgb}{0.000,0.000,1.000}
            \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
            (0, 0, 0) -- (1, 0, 0) -- (1, 1, 0) -- (0, 1, 0) -- cycle;
            \end{tikzpicture}

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P = E(P, 4)
            sage: from sage.misc.latex import latex             #not tested
            sage: latex.add_to_preamble('\\usepackage{tikz}')   #not tested
            sage: view(P)                       #not tested

        Plot using shades of gray (useful for article figures)::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.repaint([(0.9, 0.9, 0.9), (0.65,0.65,0.65), (0.4,0.4,0.4)])
            sage: P = E(P, 4)
            sage: s = P.plot_tikz()

        Plotting with various options::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: M = matrix(2,3,[float(u) for u in [1,0,-0.7071,0,1,-0.7071]])
            sage: P = E(P, 3)
            sage: s = P.plot_tikz(projmat=M, edgecolor='facecolor', scale=0.6, drawzero=True)

        Adding X, Y, Z axes using the extra code feature::

            sage: length = 1.5
            sage: space = 0.3
            sage: axes = ''
            sage: axes += "\\draw[->, thick, black] (0,0,0) -- (%s, 0, 0);\n" % length
            sage: axes += "\\draw[->, thick, black] (0,0,0) -- (0, %s, 0);\n" % length
            sage: axes += "\\node at (%s,0,0) {$x$};\n" % (length + space)
            sage: axes += "\\node at (0,%s,0) {$y$};\n" % (length + space)
            sage: axes += "\\node at (0,0,%s) {$z$};\n" % (length + space)
            sage: axes += "\\draw[->, thick, black] (0,0,0) -- (0, 0, %s);\n" % length
            sage: cube = Patch([Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)])
            sage: options = dict(scale=0.5,drawzero=True,extra_code_before=axes)
            sage: s = cube.plot_tikz(**options)
            sage: len(s)
            986
            sage: print(s)   #not tested
            \begin{tikzpicture}
            [x={(-0.433013cm,-0.250000cm)}, y={(0.433013cm,-0.250000cm)}, z={(0.000000cm,0.500000cm)}]
            \draw[->, thick, black] (0,0,0) -- (1.50000000000000, 0, 0);
            \draw[->, thick, black] (0,0,0) -- (0, 1.50000000000000, 0);
            \node at (1.80000000000000,0,0) {$x$};
            \node at (0,1.80000000000000,0) {$y$};
            \node at (0,0,1.80000000000000) {$z$};
            \draw[->, thick, black] (0,0,0) -- (0, 0, 1.50000000000000);
            \definecolor{facecolor}{rgb}{0.000,1.000,0.000}
            \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
            (0, 0, 0) -- (0, 0, 1) -- (1, 0, 1) -- (1, 0, 0) -- cycle;
            \definecolor{facecolor}{rgb}{1.000,0.000,0.000}
            \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
            (0, 0, 0) -- (0, 1, 0) -- (0, 1, 1) -- (0, 0, 1) -- cycle;
            \definecolor{facecolor}{rgb}{0.000,0.000,1.000}
            \fill[fill=facecolor, draw=black, shift={(0,0,0)}]
            (0, 0, 0) -- (1, 0, 0) -- (1, 1, 0) -- (0, 1, 0) -- cycle;
            \node[circle,fill=black,draw=black,minimum size=1.5mm,inner sep=0pt] at (0,0,0) {};
            \end{tikzpicture}
        """
        if self.dimension() != 3:
            raise NotImplementedError("Tikz plotting is implemented only for patches in three dimensions")

        if projmat is None:
            projmat = matrix(2, [-1.7320508075688772 * 0.5,
                                 1.7320508075688772 * 0.5,
                                 0, -0.5, -0.5, 1]) * scale

        e1 = projmat * vector([1, 0, 0])
        e2 = projmat * vector([0, 1, 0])
        e3 = projmat * vector([0, 0, 1])
        face_contour = self._face_contour
        color = None

        # string s contains the TiKZ code of the patch
        s = ''

        if print_tikz_env:
            s += '\\begin{tikzpicture}\n'
            s += '[x={(%fcm,%fcm)}, y={(%fcm,%fcm)}, z={(%fcm,%fcm)}]\n' % (e1[0], e1[1], e2[0], e2[1], e3[0], e3[1])

        s += extra_code_before

        for f in self:
            t = f.type()
            x, y, z = f.vector()

            if color is None or color != f.color():
                color = f.color()
                s += '\\definecolor{facecolor}{rgb}{%.3f,%.3f,%.3f}\n' % (color[0], color[1], color[2])

            s += '\\fill[fill=facecolor, draw=%s, shift={(%d,%d,%d)}]\n' % (edgecolor, x, y, z)
            s += ' -- '.join(map(str, face_contour[t])) + ' -- cycle;\n'

        s += extra_code_after

        if drawzero:
            s += '\\node[circle,fill=black,draw=black,minimum size=1.5mm,inner sep=0pt] at (0,0,0) {};\n'

        if print_tikz_env:
            s += '\\end{tikzpicture}'

        return LatexExpr(s)

    _latex_ = plot_tikz


class E1Star(SageObject):
    r"""
    A class to model the `E_1^*(\sigma)` map associated with
    a unimodular substitution `\sigma`.

    INPUT:

    - ``sigma`` -- unimodular ``WordMorphism``, i.e. such that its incidence
      matrix has determinant `\pm 1`

    - ``method`` -- 'prefix' or 'suffix' (default: ``'suffix'``);
      enables to use an alternative definition `E_1^*(\sigma)` substitutions,
      where the abelianized of the prefix` is used instead of the suffix

    .. NOTE::

        The alphabet of the domain and the codomain of `\sigma` must be
        equal, and they must be of the form ``[1, ..., d]``, where ``d``
        is a positive integer corresponding to the length of the vectors
        of the faces on which `E_1^*(\sigma)` will act.

    EXAMPLES::

        sage: from sage.combinat.e_one_star import E1Star, Face, Patch
        sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
        sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
        sage: E = E1Star(sigma)
        sage: E(P)
        Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(0, 1, -1), 2]*, [(1, 0, -1), 1]*]

    ::

        sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
        sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
        sage: E = E1Star(sigma, method='prefix')
        sage: E(P)
        Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(0, 0, 1), 1]*, [(0, 0, 1), 2]*]

    ::

        sage: x = [Face((0,0,0,0),1), Face((0,0,0,0),4)]
        sage: P = Patch(x)
        sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1,4], 4:[1]})
        sage: E = E1Star(sigma)
        sage: E(P)
        Patch: [[(0, 0, 0, 0), 3]*, [(0, 0, 0, 0), 4]*, [(0, 0, 1, -1), 3]*, [(0, 1, 0, -1), 2]*, [(1, 0, 0, -1), 1]*]
    """

    def __init__(self, sigma, method='suffix'):
        r"""
        E1Star constructor. See class doc for more information.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E
            E_1^*(1->12, 2->13, 3->1)
        """
        if not isinstance(sigma, WordMorphism):
            raise TypeError("sigma (=%s) must be an instance of WordMorphism" % sigma)

        if sigma.domain().alphabet() != sigma.codomain().alphabet():
            raise ValueError("the domain and codomain of (%s) must be the same" % sigma)

        if abs(det(matrix(sigma))) != 1:
            raise ValueError("the substitution (%s) must be unimodular" % sigma)

        first_letter = sigma.codomain().alphabet()[0]
        if first_letter not in ZZ or first_letter < 1:
            raise ValueError(f"the substitution ({sigma}) must be defined on positive integers")

        self._sigma = WordMorphism(sigma)
        self._d = self._sigma.domain().alphabet().cardinality()

        # self._base_iter is a base for the iteration of the application of self on set
        # of faces. (Exploits the linearity of `E_1^*(\sigma)` to optimize computation.)
        alphabet = self._sigma.domain().alphabet()
        X = {}
        for k in alphabet:
            subst_im = self._sigma.image(k)
            for n, letter in enumerate(subst_im):
                if method == 'suffix':
                    image_word = subst_im[n + 1:]
                elif method == 'prefix':
                    image_word = subst_im[:n]
                else:
                    raise ValueError("option 'method' can only be 'prefix' or 'suffix'")
                if letter not in X:
                    X[letter] = []
                v = self.inverse_matrix() * vector(image_word.abelian_vector())
                X[letter].append((v, k))
        self._base_iter = X

    def __eq__(self, other) -> bool:
        r"""
        Equality test for E1Star morphisms.

        INPUT:

        - ``other`` -- an object

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: s = WordMorphism({1:[1,3], 2:[1,2,3], 3:[3]})
            sage: t = WordMorphism({1:[1,2,3], 2:[2,3], 3:[3]})
            sage: S = E1Star(s)
            sage: T = E1Star(t)
            sage: S == T
            False
            sage: S2 = E1Star(s, method='prefix')
            sage: S == S2
            False
        """
        return (isinstance(other, E1Star) and self._base_iter == other._base_iter)

    def __call__(self, patch, iterations=1) -> Patch:
        r"""
        Applies a generalized substitution to a Patch; this returns a new object.

        The color of every new face in the image is given the same color as its preimage.

        INPUT:

        - ``patch`` -- a patch
        - ``iterations`` -- integer (default: 1); number of iterations

        OUTPUT: a patch

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E(P)
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(0, 1, -1), 2]*, [(1, 0, -1), 1]*]
            sage: E(P, iterations=4)
            Patch of 31 faces

        TESTS:

        We test that iterations=0 works (see :issue:`10699`)::

            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E(P, iterations=0)
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
        """
        if iterations == 0:
            return Patch(patch)
        elif iterations < 0:
            raise ValueError("iterations (=%s) must be >= 0" % iterations)
        else:
            old_faces = patch
            for _ in range(iterations):
                new_faces = []
                for f in old_faces:
                    new_faces.extend(self._call_on_face(f, color=f.color()))
                old_faces = new_faces
            return Patch(new_faces)

    def __mul__(self, other) -> E1Star:
        r"""
        Return the product of ``self`` and ``other``.

        The product satisfies the following rule:
        `E_1^*(\sigma\circ\sigma') = E_1^*(\sigma')` \circ  E_1^*(\sigma)`

        INPUT:

        - ``other`` -- an instance of E1Star

        OUTPUT: an instance of E1Star

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: s = WordMorphism({1:[2],2:[3],3:[1,2]})
            sage: t = WordMorphism({1:[1,3,1],2:[1],3:[1,1,3,2]})
            sage: E1Star(s) * E1Star(t)
            E_1^*(1->1, 2->1132, 3->1311)
            sage: E1Star(t * s)
            E_1^*(1->1, 2->1132, 3->1311)
        """
        if not isinstance(other, E1Star):
            raise TypeError("other (=%s) must be an instance of E1Star" % other)
        return E1Star(other.sigma() * self.sigma())

    def __repr__(self) -> str:
        r"""
        String representation of a patch.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E
            E_1^*(1->12, 2->13, 3->1)
        """
        return "E_1^*(%s)" % str(self._sigma)

    def _call_on_face(self, face, color=None):
        r"""
        Return an iterator of faces obtained by applying ``self`` on the face.

        INPUT:

        - ``face`` -- a face
        - ``color`` -- string (default: ``None``); RGB tuple or color

        OUTPUT: iterator of faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,2,0), 1)
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: list(E._call_on_face(f))
            [[(3, 0, -3), 1]*, [(2, 1, -3), 2]*, [(2, 0, -2), 3]*]
        """
        if len(face.vector()) != self._d:
            raise ValueError("the dimension of the faces must be equal to the size of the alphabet of the substitution")
        x_new = self.inverse_matrix() * face.vector()
        t = face.type()
        return (Face(x_new + v, k, color=color) for v, k in self._base_iter[t])

    @cached_method
    def matrix(self):
        r"""
        Return the matrix associated with ``self``.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E.matrix()
            [1 1 1]
            [1 0 0]
            [0 1 0]
        """
        return self._sigma.incidence_matrix()

    @cached_method
    def inverse_matrix(self):
        r"""
        Return the inverse of the matrix associated with ``self``.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E.inverse_matrix()
            [ 0  1  0]
            [ 0  0  1]
            [ 1 -1 -1]
        """
        return self.matrix().inverse()

    def sigma(self) -> WordMorphism:
        r"""
        Return the ``WordMorphism`` associated with ``self``.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E.sigma()
            WordMorphism: 1->12, 2->13, 3->1
        """
        return self._sigma
