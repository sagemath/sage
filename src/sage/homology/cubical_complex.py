r"""
Finite cubical complexes

AUTHORS:

- John H. Palmieri (2009-08)

This module implements the basic structure of finite cubical
complexes.  For full mathematical details, see Kaczynski, Mischaikow,
and Mrozek [KMM]_, for example.

Cubical complexes are topological spaces built from gluing together
cubes of various dimensions; the collection of cubes must be closed
under taking faces, just as with a simplicial complex.  In this
context, a "cube" means a product of intervals of length 1 or length 0
(degenerate intervals), with integer endpoints, and its faces are
obtained by using the nondegenerate intervals: if `C` is a cube -- a
product of degenerate and nondegenerate intervals -- and if `[i,i+1]`
is the `k`-th nondegenerate factor, then `C` has two faces indexed by
`k`: the cubes obtained by replacing `[i, i+1]` with `[i, i]` or
`[i+1, i+1]`.

So to construct a space homeomorphic to a circle as a cubical complex,
we could take for example the four line segments in the plane from
`(0,2)` to `(0,3)` to `(1,3)` to `(1,2)` to `(0,2)`.  In Sage, this is
done with the following command::

    sage: S1 = CubicalComplex([([0,0], [2,3]), ([0,1], [3,3]), ([0,1], [2,2]), ([1,1], [2,3])]); S1
    Cubical complex with 4 vertices and 8 cubes

The argument to ``CubicalComplex`` is a list of the maximal "cubes" in
the complex.  Each "cube" can be an instance of the class ``Cube`` or
a list (or tuple) of "intervals", and an "interval" is a pair of
integers, of one of the two forms `[i, i]` or `[i, i+1]`.  So the
cubical complex ``S1`` above has four maximal cubes::

    sage: S1.maximal_cells()
    {[0,0] x [2,3], [1,1] x [2,3], [0,1] x [3,3], [0,1] x [2,2]}

The first of these, for instance, is the product of the degenerate
interval `[0,0]` with the unit interval `[2,3]`: this is the line
segment in the plane from `(0,2)` to `(0,3)`.  We could form a
topologically equivalent space by inserting some degenerate simplices::

    sage: S1.homology()
    {0: 0, 1: Z}
    sage: X = CubicalComplex([([0,0], [2,3], [2]), ([0,1], [3,3], [2]), ([0,1], [2,2], [2]), ([1,1], [2,3], [2])])
    sage: X.homology()
    {0: 0, 1: Z}

Topologically, the cubical complex ``X`` consists of four edges of a
square in `\RR^3`: the same unit square as ``S1``, but embedded in
`\RR^3` with `z`-coordinate equal to 2.  Thus ``X`` is homeomorphic to
``S1`` (in fact, they're "cubically equivalent"), and this is
reflected in the fact that they have isomorphic homology groups.

REFERENCES:

.. [KMM] Tomasz Kaczynski, Konstantin Mischaikow, and Marian Mrozek,
         "Computational Homology", Springer-Verlag (2004).

.. note::

   This class derives from
   :class:`~sage.homology.cell_complex.GenericCellComplex`, and so
   inherits its methods.  Some of those methods are not listed here;
   see the :mod:`Generic Cell Complex <sage.homology.cell_complex>`
   page instead.
"""

from copy import copy
from sage.homology.cell_complex import GenericCellComplex
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.sets.set import Set
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.homology.chain_complex import ChainComplex
from sage.graphs.graph import Graph

class Cube(SageObject):
    r"""
    Define a cube for use in constructing a cubical complex.

    "Elementary cubes" are products of intervals with integer
    endpoints, each of which is either a unit interval or a degenerate
    (length 0) interval; for example,

    .. math::

       [0,1] \times [3,4] \times [2,2] \times [1,2]

    is a 3-dimensional cube (since one of the intervals is degenerate)
    embedded in `\RR^4`.

    :param data: list or tuple of terms of the form ``(i,i+1)`` or
      ``(i,i)`` or ``(i,)`` -- the last two are degenerate intervals.
    :return: an elementary cube

    Each cube is stored in a standard form: a tuple of tuples, with a
    nondegenerate interval ``[j,j]`` represented by ``(j,j)``, not
    ``(j,)``.  (This is so that for any interval ``I``, ``I[1]`` will
    produce a value, not an ``IndexError``.)

    EXAMPLES::

        sage: from sage.homology.cubical_complex import Cube
        sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]]); C
        [1,2] x [5,5] x [6,7] x [-1,0]
        sage: C.dimension() # number of nondegenerate intervals
        3
        sage: C.nondegenerate_intervals()  # indices of these intervals
        [0, 2, 3]
        sage: C.face(1, upper=False)
        [1,2] x [5,5] x [6,6] x [-1,0]
        sage: C.face(1, upper=True)
        [1,2] x [5,5] x [7,7] x [-1,0]
        sage: Cube(()).dimension()  # empty cube has dimension -1
        -1
    """
    def __init__(self, data):
        """
        Define a cube for use in constructing a cubical complex.

        See ``Cube`` for more information.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]]); C # indirect doctest
            [1,2] x [5,5] x [6,7] x [-1,0]
            sage: C == loads(dumps(C))
            True
        """
        if isinstance(data, Cube):
            data = tuple(data)
        new_data = []
        nondegenerate = []
        i = 0
        for x in data:
            if len(x) == 2:
                try:
                    Integer(x[0])
                except TypeError:
                    raise ValueError, "The interval %s is not of the correct form" % x
                if x[0] + 1 == x[1]:
                    nondegenerate.append(i)
                elif x[0] != x[1]:
                    raise ValueError, "The interval %s is not of the correct form" % x
                new_data.append(tuple(x))
            elif len(x) == 1:
                y = tuple(x)
                new_data.append(y+y)
            elif len(x) != 1:
                raise ValueError, "The interval %s is not of the correct form" % x
            i += 1
        self.__tuple = tuple(new_data)
        self.__nondegenerate = nondegenerate

    def tuple(self):
        """
        The tuple attached to this cube.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]])
            sage: C.tuple()
            ((1, 2), (5, 5), (6, 7), (-1, 0))
        """
        return self.__tuple

    def is_face(self, other):
        """
        Return True iff this cube is a face of other.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C1 = Cube([[1,2], [5,], [6,7], [-1, 0]])
            sage: C2 = Cube([[1,2], [5,], [6,], [-1, 0]])
            sage: C1.is_face(C2)
            False
            sage: C1.is_face(C1)
            True
            sage: C2.is_face(C1)
            True
        """
        def is_subinterval(i1, i2):
            return ((i1[0] == i2[0] and i1[1] == i2[1]) or
                    (i1[0] == i2[0] and i1[1] == i2[0]) or
                    (i1[0] == i2[1] and i1[1] == i2[1]))

        t = self.tuple()
        u = other.tuple()
        embed = len(u)
        if len(t) == embed: # these must be equal for self to be a face of other
            return all([is_subinterval(t[i], u[i]) for i in range(embed)])
        else:
            return False

    def _translate(self, vec):
        """
        Translate ``self`` by ``vec``.

        :param vec: anything which can be converted to a tuple of integers
        :return: the translation of ``self`` by ``vec``
        :rtype: Cube

        If ``vec`` is shorter than the list of intervals forming the
        cube, pad with zeroes, and similarly if the cube's defining
        tuple is too short.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]])
            sage: C._translate((-12,))
            [-11,-10] x [5,5] x [6,7] x [-1,0]
            sage: C._translate((0, 0, 0, 0, 0, 5))
            [1,2] x [5,5] x [6,7] x [-1,0] x [0,0] x [5,5]
        """
        t = self.__tuple
        embed = max(len(t), len(vec))
        t = t + ((0,0),) * (embed-len(t))
        vec = tuple(vec) + (0,) * (embed-len(vec))
        new = []
        for (a,b) in zip(t, vec):
            new.append([a[0]+b, a[1]+b])
        return Cube(new)

    def __getitem__(self, n):
        """
        Return the nth interval in this cube.

        :param n: an integer
        :return: tuple representing the `n`-th interval in the cube.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]])
            sage: C[2]
            (6, 7)
            sage: C[1]
            (5, 5)
        """
        return self.__tuple.__getitem__(n)

    def __iter__(self):
        """
        Iterator for the intervals of this cube.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]])
            sage: [x[0] for x in C]
            [1, 5, 6, -1]
        """
        return self.__tuple.__iter__()

    def __add__(self, other):
        """
        Cube obtained by concatenating the underlying tuples of the
        two arguments.

        :param other: another cube
        :return: the product of ``self`` and ``other``, as a Cube

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [3,]])
            sage: D = Cube([[4], [0,1]])
            sage: C.product(D)
            [1,2] x [3,3] x [4,4] x [0,1]

        You can also use ``__add__`` or ``+`` or ``__mul__`` or ``*``::

            sage: D * C
            [4,4] x [0,1] x [1,2] x [3,3]
            sage: D + C * C
            [4,4] x [0,1] x [1,2] x [3,3] x [1,2] x [3,3]
        """
        return Cube(self.__tuple.__add__(other.__tuple))

    # the __add__ operation actually produces the product of the two cubes
    __mul__ = __add__
    product = __add__

    def nondegenerate_intervals(self):
        """
        The list of indices of nondegenerate intervals of this cube.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]])
            sage: C.nondegenerate_intervals()
            [0, 2, 3]
            sage: C = Cube([[1,], [5,], [6,], [-1,]])
            sage: C.nondegenerate_intervals()
            []
        """
        return self.__nondegenerate

    def dimension(self):
        """
        The dimension of this cube: the number of its nondegenerate intervals.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]])
            sage: C.dimension()
            3
            sage: C = Cube([[1,], [5,], [6,], [-1,]])
            sage: C.dimension()
            0
            sage: Cube([]).dimension()  # empty cube has dimension -1
            -1
        """
        if len(self.__tuple) == 0:  # empty cube
            return -1
        return len(self.nondegenerate_intervals())

    def face(self, n, upper=True):
        """
        The nth primary face of this cube.

        :param n: an integer between 0 and one less than the dimension
          of this cube
        :param upper: if True, return the "upper" nth primary face;
          otherwise, return the "lower" nth primary face.
        :type upper: boolean; optional, default=True
        :return: the cube obtained by replacing the nth non-degenrate
          interval with either its upper or lower endpoint.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [5,], [6,7], [-1, 0]]); C
            [1,2] x [5,5] x [6,7] x [-1,0]
            sage: C.face(0)
            [2,2] x [5,5] x [6,7] x [-1,0]
            sage: C.face(0, upper=False)
            [1,1] x [5,5] x [6,7] x [-1,0]
            sage: C.face(1)
            [1,2] x [5,5] x [7,7] x [-1,0]
            sage: C.face(2, upper=False)
            [1,2] x [5,5] x [6,7] x [-1,-1]
            sage: C.face(3)
            Traceback (most recent call last):
            ...
            ValueError: Can only compute the nth face if 0 <= n < dim.
        """
        if n < 0 or n >= self.dimension():
            raise ValueError, "Can only compute the nth face if 0 <= n < dim."
        idx = self.nondegenerate_intervals()[n]
        t = self.__tuple
        if upper:
            new = t[idx][1]
        else:
            new = t[idx][0]
        return Cube(t[0:idx] + ((new,new),) + t[idx+1:])

    def faces(self):
        """
        The list of faces (of codimension 1) of this cube.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [3,4]])
            sage: C.faces()
            [[2,2] x [3,4], [1,2] x [4,4], [1,1] x [3,4], [1,2] x [3,3]]
        """
        upper = [self.face(i,True) for i in range(self.dimension())]
        lower = [self.face(i,False) for i in range(self.dimension())]
        return upper + lower

    def faces_as_pairs(self):
        """
        The list of faces (of codimension 1) of this cube, as pairs
        (upper, lower).

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [3,4]])
            sage: C.faces_as_pairs()
            [([2,2] x [3,4], [1,1] x [3,4]), ([1,2] x [4,4], [1,2] x [3,3])]
        """
        upper = [self.face(i,True) for i in range(self.dimension())]
        lower = [self.face(i,False) for i in range(self.dimension())]
        return zip(upper,lower)

    def _compare_for_gluing(self, other):
        r"""
        Given two cubes ``self`` and ``other``, describe how to
        transform them so that they become equal.

        :param other: a cube of the same dimension as ``self``
        :return: a triple ``(insert_self, insert_other, translate)``.
          ``insert_self`` is a tuple with entries ``(index, (list of
          degenerate intervals))``.  ``insert_other`` is similar.
          ``translate`` is a tuple of integers, suitable as a second
          argument for the ``_translate`` method.

        To do this, ``self`` and ``other`` must have the same
        dimension; degenerate intervals from ``other`` are added to
        ``self``, and vice versa.  Intervals in ``other`` are
        translated so that they coincide with the intervals in
        ``self``.  The output is a triple, as noted above: in the
        tuple ``insert_self``, an entry like ``(3, (3, 4, 0))`` means
        that in position 3 in ``self``, insert the degenerate
        intervals ``[3,3]``, ``[4,4]``, and ``[0,0]``.  The same goes
        for ``insert_other``.  After applying the translations to the
        cube ``other``, call ``_translate`` with argument the tuple
        ``translate``.

        This is used in forming connected sums of cubical complexes:
        the two complexes are modified, via this method, so that they
        have a cube which matches up, then those matching cubes are
        removed.

        In the example below, this method is called with arguments
        ``C1`` and ``C2``, where

        .. math::

            C1 = [0,1] \times [3] \times [4] \times [6,7] \\
            C2 = [2] \times [7,8] \times [9] \times [1,2] \times [0] \times [5]

        To C1, we need to add [2] in position 0 and [0] and [5] in
        position 5.  To C2, we need to add [4] in position 3.  Once
        this has been done, we need to translate the new C2 by the
        vector ``(0, -7, -6, 0, 5, 0, 0)``.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C1 = Cube([[0,1], [3], [4], [6,7]])
            sage: C2 = Cube([[2], [7,8], [9], [1,2], [0], [5]])
            sage: C1._compare_for_gluing(C2)
            ([(0, ((2, 2),)), (5, ((0, 0), (5, 5)))], [(3, ((4, 4),))], [0, -7, -6, 0, 5, 0, 0])

            sage: C1 = Cube([[1,1], [0,1]])
            sage: C2 = Cube([[2,3], [4,4], [5,5]])
            sage: C1._compare_for_gluing(C2)
            ([(2, ((4, 4), (5, 5)))], [(0, ((1, 1),))], [0, -2, 0, 0])
        """
        d = self.dimension()
        if d != other.dimension():
            raise ValueError, "Cubes must be of the same dimension."
        insert_self = []
        insert_other = []
        translate = []
        self_tuple = self.tuple()
        other_tuple = other.tuple()
        nondegen = (zip(self.nondegenerate_intervals(),
                        other.nondegenerate_intervals())
                    + [(len(self_tuple), len(other_tuple))])
        old = (-1, -1)
        self_added = 0
        other_added = 0

        for current in nondegen:
            # number of positions between nondegenerate intervals:
            self_diff = current[0] - old[0]
            other_diff = current[1] - old[1]
            diff = self_diff - other_diff

            if diff < 0:
                insert_self.append((old[0] + self_diff + self_added,
                                    other.tuple()[current[1]+diff:current[1]]))
                common_terms = self_diff
                diff = -diff
                self_added += diff
            elif diff > 0:
                insert_other.append((old[1] + other_diff + other_added,
                                     self.tuple()[current[0]-diff:current[0]]))
                common_terms = other_diff
                other_added += diff
            else:
                common_terms = other_diff

            if old[0] > -1:
                translate.extend([self_tuple[old[0]+idx][0] -
                                  other_tuple[old[1]+idx][0] for idx in
                                  range(common_terms)])
            translate.extend(diff*[0])
            old = current

        return (insert_self, insert_other, translate)

    def _triangulation_(self):
        r"""
        Triangulate this cube by "pulling vertices," as described by
        Hetyei.  Return a list of simplices which triangulate
        ``self``.

        ALGORITHM:

        If the cube is given by

        .. math::

           C = [i_1, j_1] \times [i_2, j_2] \times ... \times [i_k, j_k]

        let `v_1` be the "upper" corner of `C`: `v` is the point
        `(j_1, ..., j_k)`.  Choose a coordinate `n` where the interval
        `[i_n, j_n]` is non-degenerate and form `v_2` by replacing
        `j_n` by `i_n`; repeat to define `v_3`, etc.  The last vertex
        so defined will be `(i_1, ..., i_k)`.  These vertices define a
        simplex, as do the vertices obtained by making different
        choices at each stage.  Return the list of such simplices;
        thus if `C` is `n`-dimensional, then it is subdivided into
        `n!` simplices.

        REFERENCES:

        - G. Hetyei, "On the Stanley ring of a cubical complex",
          Discrete Comput. Geom. 14 (1995), 305-330.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C = Cube([[1,2], [3,4]]); C
            [1,2] x [3,4]
            sage: C._triangulation_()
            [((1, 3), (1, 4), (2, 4)), ((1, 3), (2, 3), (2, 4))]
            sage: C = Cube([[1,2], [3,4], [8,9]])
            sage: len(C._triangulation_())
            6
        """
        from sage.homology.simplicial_complex import Simplex
        if self.dimension() < 0: # the empty cube
            return [Simplex(())] # the empty simplex
        v = tuple([max(j) for j in self.tuple()])
        if self.dimension() == 0: # just v
            return [Simplex((v,))]
        simplices = []
        for i in range(self.dimension()):
            for S in self.face(i, upper=False)._triangulation_():
                simplices.append(S.join(Simplex((v,)), rename_vertices=False))
        return simplices

    def __cmp__(self, other):
        """
        Return True iff this cube is the same as ``other``: that is,
        if they are the product of the same intervals in the same
        order.

        :param other: another cube

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C1 = Cube([[1,1], [2,3], [4,5]])
            sage: C2 = Cube([[1], [2,3], [4,5]])
            sage: C3 = Cube([[0], [2,3], [4,5]])
            sage: C1 == C2  # indirect doctest
            True
            sage: C1 == C3  # indirect doctest
            False
        """
        if self.__tuple == other.__tuple:
            return 0
        else:
            return -1

    def __hash__(self):
        """
        Hash value for this cube.  This computes the hash value of the
        underlying tuple, since this is what's important when testing
        equality.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C1 = Cube([[1,1], [2,3], [4,5]])
            sage: C1.__hash__()
            837272820736660832  # 64-bit
            -1004989088  # 32-bit
        """
        return hash(self.__tuple)

    def _repr_(self):
        """
        Print representation of a cube.

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C1 = Cube([[1,1], [2,3], [4,5]])
            sage: C1
            [1,1] x [2,3] x [4,5]
            sage: C1._repr_()
            '[1,1] x [2,3] x [4,5]'
        """
        s = ["[%s,%s]"%(str(x), str(y)) for (x,y) in self.__tuple]
        return " x ".join(s)

    def _latex_(self):
        r"""
        LaTeX representation of a cube..

        EXAMPLES::

            sage: from sage.homology.cubical_complex import Cube
            sage: C1 = Cube([[1,1], [2,3], [4,5]])
            sage: latex(C1)
            [1,1] \times [2,3] \times [4,5]
            sage: C1._latex_()
            '[1,1] \\times [2,3] \\times [4,5]'
        """
        return self._repr_().replace('x', r'\times')


class CubicalComplex(GenericCellComplex):
    r"""
    Define a cubical complex.

    :param maximal_faces: set of maximal faces
    :param maximality_check: see below
    :type maximality_check: boolean; optional, default True
    :return: a cubical complex

    ``maximal_faces`` should be a list or tuple or set (or anything
    which may be converted to a set) of "cubes": instances of the
    class :class:`Cube`, or lists or tuples suitable for conversion to
    cubes.  These cubes are the maximal cubes in the complex.

    In addition, ``maximal_faces`` may be a cubical complex, in which
    case that complex is returned.  Also, ``maximal_faces`` may
    instead be any object which has a ``_cubical_`` method (e.g., a
    simplicial complex); then that method is used to convert the
    object to a cubical complex.

    If ``maximality_check`` is True, check that each maximal face is,
    in fact, maximal. In this case, when producing the internal
    representation of the cubical complex, omit those that are not.
    It is highly recommended that this be True; various methods for
    this class may fail if faces which are claimed to be maximal are
    in fact not.

    EXAMPLES:

    The empty complex, consisting of one cube, the empty cube::

        sage: CubicalComplex()
        Cubical complex with 0 vertices and 1 cube

    A "circle" (four edges connecting the vertices (0,2), (0,3),
    (1,2), and (1,3))::

        sage: S1 = CubicalComplex([([0,0], [2,3]), ([0,1], [3,3]), ([0,1], [2,2]), ([1,1], [2,3])])
        sage: S1
        Cubical complex with 4 vertices and 8 cubes
        sage: S1.homology()
        {0: 0, 1: Z}

    A set of five points and its product with ``S1``::

        sage: pts = CubicalComplex([([0],), ([3],), ([6],), ([-12],), ([5],)])
        sage: pts
        Cubical complex with 5 vertices and 5 cubes
        sage: pts.homology()
        {0: Z x Z x Z x Z}
        sage: X = S1.product(pts); X
        Cubical complex with 20 vertices and 40 cubes
        sage: X.homology()
        {0: Z x Z x Z x Z, 1: Z^5}

    Converting a simplicial complex to a cubical complex::

        sage: S2 = simplicial_complexes.Sphere(2)
        sage: C2 = CubicalComplex(S2)
        sage: all([C2.homology(n) == S2.homology(n) for n in range(3)])
        True

    You can get the set of maximal cells or a dictionary of all cells::

        sage: X.maximal_cells()
        {[0,0] x [2,3] x [-12,-12], [0,1] x [3,3] x [5,5], [0,1] x [2,2] x [3,3], [0,1] x [2,2] x [0,0], [0,1] x [3,3] x [6,6], [1,1] x [2,3] x [0,0], [0,1] x [2,2] x [-12,-12], [0,0] x [2,3] x [6,6], [1,1] x [2,3] x [-12,-12], [1,1] x [2,3] x [5,5], [0,1] x [2,2] x [5,5], [0,1] x [3,3] x [3,3], [1,1] x [2,3] x [3,3], [0,0] x [2,3] x [5,5], [0,1] x [3,3] x [0,0], [1,1] x [2,3] x [6,6], [0,1] x [2,2] x [6,6], [0,0] x [2,3] x [0,0], [0,0] x [2,3] x [3,3], [0,1] x [3,3] x [-12,-12]}
        sage: S1.cells()
        {0: set([[1,1] x [2,2], [0,0] x [2,2], [1,1] x [3,3], [0,0] x [3,3]]), 1: set([[0,1] x [3,3], [1,1] x [2,3], [0,0] x [2,3], [0,1] x [2,2]]), -1: set([])}

    Chain complexes, homology, and cohomology::

        sage: T = S1.product(S1); T
        Cubical complex with 16 vertices and 64 cubes
        sage: T.chain_complex()
        Chain complex with at most 3 nonzero terms over Integer Ring
        sage: T.homology(base_ring=QQ)
        {0: Vector space of dimension 0 over Rational Field,
         1: Vector space of dimension 2 over Rational Field,
         2: Vector space of dimension 1 over Rational Field}
        sage: RP2 = cubical_complexes.RealProjectivePlane()
        sage: RP2.cohomology(dim=[1, 2], base_ring=GF(2))
        {1: Vector space of dimension 1 over Finite Field of size 2,
         2: Vector space of dimension 1 over Finite Field of size 2}

    Joins are not implemented::

        sage: S1.join(S1)
        Traceback (most recent call last):
        ...
        NotImplementedError: Joins are not implemented for cubical complexes.

    Therefore, neither are cones or suspensions.
    """
    def __init__(self, maximal_faces=[], **kwds):
        r"""
        Define a cubical complex.  See ``CubicalComplex`` for more
        documentation.

        EXAMPLES::

            sage: X = CubicalComplex([([0,0], [2,3]), ([0,1], [3,3]), ([0,1], [2,2]), ([1,1], [2,3])]); X
            Cubical complex with 4 vertices and 8 cubes
            sage: X == loads(dumps(X))
            True
        """
        maximality_check = kwds.get('maximality_check', True)

        C = None
        if isinstance(maximal_faces, CubicalComplex):
            C = maximal_faces
        try:
            C = maximal_faces._cubical_()
        except AttributeError:
            pass
        if C is not None:
            self._facets = copy(C._facets)
            self._cells = copy(C._cells)
            self._complex = copy(C._complex)
            return

        good_faces = []
        maximal_cubes = [Cube(f) for f in maximal_faces]
        for face in maximal_cubes:
            # check whether each given face is actually maximal
            face_is_maximal = True
            if maximality_check:
                faces_to_be_removed = []
                for other in good_faces:
                    if other.is_face(face):
                        faces_to_be_removed.append(other)
                    elif face_is_maximal:
                        face_is_maximal = not face.is_face(other)
                for x in faces_to_be_removed:
                    good_faces.remove(x)
            if face_is_maximal:
                good_faces += [face]
        # if no maximal faces, add the empty face as a facet
        if len(maximal_cubes) == 0:
            good_faces.append(Cube(()))
        # self._facets: tuple of facets
        self._facets = tuple(good_faces)
        # self._cells: dictionary of dictionaries of faces.  The main
        # dictionary is keyed by subcomplexes, and each value is a
        # dictionary keyed by dimension.  This should be empty until
        # needed -- that is, until the faces method is called
        self._cells = {}
        # self._complex: dictionary indexed by dimension d, base_ring,
        # etc.: differential from dim d to dim d-1 in the associated
        # chain complex.  thus to get the differential in the cochain
        # complex from dim d-1 to dim d, take the transpose of this
        # one.
        self._complex = {}

    def maximal_cells(self):
        """
        The set of maximal cells (with respect to inclusion) of this
        cubical complex.

        :return: Set of maximal cells

        This just returns the set of cubes used in defining the
        cubical complex, so if the complex was defined with no
        maximality checking, none is done here, either.

        EXAMPLES::

            sage: interval = cubical_complexes.Cube(1)
            sage: interval
            Cubical complex with 2 vertices and 3 cubes
            sage: interval.maximal_cells()
            {[0,1]}
            sage: interval.product(interval).maximal_cells()
            {[0,1] x [0,1]}
        """
        return Set(self._facets)

    def __cmp__(self,other):
        r"""
        Return True if the set of maximal cells is the same for
        ``self`` and ``other``.

        :param other: another cubical complex
        :return: True if the set of maximal cells is the same for ``self`` and ``other``
        :rtype: bool

        EXAMPLES::

            sage: I1 = cubical_complexes.Cube(1)
            sage: I2 = cubical_complexes.Cube(1)
            sage: I1.product(I2) == I2.product(I1)
            True
            sage: I1.product(I2.product(I2)) == I2.product(I1.product(I1))
            True
            sage: S1 = cubical_complexes.Sphere(1)
            sage: I1.product(S1) == S1.product(I1)
            False
        """
        if self.maximal_cells() == other.maximal_cells():
            return 0
        else:
            return -1

    def is_subcomplex(self, other):
        r"""
        Return True if ``self`` is a subcomplex of ``other``.

        :param other: a cubical complex

        Each maximal cube of ``self`` must be a face of a maximal cube
        of ``other`` for this to be True.

        EXAMPLES::

            sage: S1 = cubical_complexes.Sphere(1)
            sage: C0 = cubical_complexes.Cube(0)
            sage: C1 = cubical_complexes.Cube(1)
            sage: cyl = S1.product(C1)
            sage: end = S1.product(C0)
            sage: end.is_subcomplex(cyl)
            True
            sage: cyl.is_subcomplex(end)
            False

        The embedding of the cubical complex is important here::

            sage: C2 = cubical_complexes.Cube(2)
            sage: C1.is_subcomplex(C2)
            False
            sage: C1.product(C0).is_subcomplex(C2)
            True

        ``C1`` is not a subcomplex of ``C2`` because it's not embedded
        in `\RR^2`.  On the other hand, ``C1 x C0`` is a face of
        ``C2``.  Look at their maximal cells::

            sage: C1.maximal_cells()
            {[0,1]}
            sage: C2.maximal_cells()
            {[0,1] x [0,1]}
            sage: C1.product(C0).maximal_cells()
            {[0,1] x [0,0]}
        """
        other_facets = other.maximal_cells()
        answer = True
        for cube in self.maximal_cells():
            answer = answer and any([cube.is_face(other_cube)
                                     for other_cube in other_facets])
        return answer

    def cells(self, subcomplex=None):
        """
        The cells of this cubical complex, in the form of a dictionary:
        the keys are integers, representing dimension, and the value
        associated to an integer d is the list of d-cells.

        If the optional argument ``subcomplex`` is present, then
        return only the faces which are *not* in the subcomplex.

        :param subcomplex: a subcomplex of this cubical complex
        :type subcomplex: a cubical complex; optional, default None
        :return: cells of this complex not contained in ``subcomplex``
        :rtype: dictionary

        EXAMPLES::

            sage: S2 = cubical_complexes.Sphere(2)
            sage: S2.cells()[2]
            set([[1,1] x [0,1] x [0,1], [0,1] x [0,0] x [0,1], [0,1] x [1,1] x [0,1], [0,0] x [0,1] x [0,1], [0,1] x [0,1] x [1,1], [0,1] x [0,1] x [0,0]])
        """
        if subcomplex not in self._cells:
            if subcomplex is not None and subcomplex.dimension() > -1:
                if not subcomplex.is_subcomplex(self):
                    raise ValueError, "The 'subcomplex' is not actually a subcomplex."
            # Cells is the dictionary of cells in self but not in
            # subcomplex, indexed by dimension
            Cells = {}
            # sub_facets is the dictionary of facets in the subcomplex
            sub_facets = {}
            dimension = max([cube.dimension() for cube in self._facets])
            # initialize the lists: add each maximal cube to Cells and sub_facets
            for i in range(-1,dimension+1):
                Cells[i] = set([])
                sub_facets[i] = set([])
            for f in self._facets:
                Cells[f.dimension()].add(f)
            if subcomplex is not None:
                for g in subcomplex._facets:
                    dim = g.dimension()
                    Cells[dim].discard(g)
                    sub_facets[dim].add(g)
            # bad_faces is the set of faces in the subcomplex in the
            # current dimension
            bad_faces = sub_facets[dimension]
            for dim in range(dimension, -1, -1):
                # bad_bdries = boundaries of bad_faces: things to be
                # discarded in dim-1
                bad_bdries = sub_facets[dim-1]
                for f in bad_faces:
                    bad_bdries.update(f.faces())
                for f in Cells[dim]:
                    Cells[dim-1].update(set(f.faces()).difference(bad_bdries))
                bad_faces = bad_bdries
            self._cells[subcomplex] = Cells
        return self._cells[subcomplex]

    def n_cubes(self, n, subcomplex=None):
        """
        The set of cubes of dimension n of this cubical complex.
        If the optional argument ``subcomplex`` is present, then
        return the ``n``-dimensional cubes which are *not* in the
        subcomplex.

        :param n: dimension
        :type n: integer
        :param subcomplex: a subcomplex of this cubical complex
        :type subcomplex: a cubical complex; optional, default None
        :return: cells in dimension ``n``
        :rtype: set

        EXAMPLES::

            sage: C = cubical_complexes.Cube(3)
            sage: C.n_cubes(3)
            set([[0,1] x [0,1] x [0,1]])
            sage: C.n_cubes(2)
            set([[1,1] x [0,1] x [0,1], [0,1] x [0,0] x [0,1], [0,1] x [1,1] x [0,1], [0,0] x [0,1] x [0,1], [0,1] x [0,1] x [1,1], [0,1] x [0,1] x [0,0]])
        """
        return set(self.n_cells(n, subcomplex))

    def chain_complex(self, **kwds):
        r"""
        The chain complex associated to this cubical complex.

        :param dimensions: if None, compute the chain complex in all
           dimensions.  If a list or tuple of integers, compute the
           chain complex in those dimensions, setting the chain groups
           in all other dimensions to zero.  NOT IMPLEMENTED YET: this
           function always returns the entire chain complex
        :param base_ring: commutative ring
        :type base_ring: optional, default ZZ
        :param subcomplex: a subcomplex of this cubical complex.
           Compute the chain complex relative to this subcomplex.
        :type subcomplex: optional, default empty
        :param augmented: If True, return the augmented chain complex
           (that is, include a class in dimension `-1` corresponding
           to the empty cell).  This is ignored if ``dimensions`` is
           specified.
        :type augmented: boolean; optional, default False
        :param cochain: If True, return the cochain complex (that is,
           the dual of the chain complex).
        :type cochain: boolean; optional, default False
        :param verbose: If True, print some messages as the chain
           complex is computed.
        :type verbose: boolean; optional, default False
        :param check_diffs: If True, make sure that the chain complex
           is actually a chain complex: the differentials are
           composable and their product is zero.
        :type check_diffs: boolean; optional, default False

        .. note::

           If subcomplex is nonempty, then the argument ``augmented``
           has no effect: the chain complex relative to a nonempty
           subcomplex is zero in dimension `-1`.

        EXAMPLES::

            sage: S2 = cubical_complexes.Sphere(2)
            sage: S2.chain_complex()
            Chain complex with at most 3 nonzero terms over Integer Ring
            sage: Prod = S2.product(S2); Prod
            Cubical complex with 64 vertices and 676 cubes
            sage: Prod.chain_complex()
            Chain complex with at most 5 nonzero terms over Integer Ring
            sage: Prod.chain_complex(base_ring=QQ)
            Chain complex with at most 5 nonzero terms over Rational Field
            sage: C1 = cubical_complexes.Cube(1)
            sage: S0 = cubical_complexes.Sphere(0)
            sage: C1.chain_complex(subcomplex=S0)
            Chain complex with at most 1 nonzero terms over Integer Ring
            sage: C1.homology(subcomplex=S0)
            {0: 0, 1: Z}
        """
        augmented = kwds.get('augmented', False)
        cochain = kwds.get('cochain', False)
        verbose = kwds.get('verbose', False)
        check_diffs = kwds.get('check_diffs', False)
        base_ring = kwds.get('base_ring', ZZ)
        dimensions = kwds.get('dimensions', None)
        subcomplex = kwds.get('subcomplex', None)

        # initialize subcomplex
        if subcomplex is None:
            subcomplex = CubicalComplex()
        else:
            # subcomplex is not empty, so don't augment the chain complex
            augmented = False
        differentials = {}
        if augmented:
            empty_cell = 1  # number of (-1)-dimensional cubes
        else:
            empty_cell = 0
        vertices = self.n_cubes(0, subcomplex=subcomplex)
        n = len(vertices)
        mat = matrix(base_ring, empty_cell, n, n*empty_cell*[1])
        if cochain:
            differentials[-1] = mat.transpose()
        else:
            differentials[0] = mat
        current = vertices
        # now loop from 1 to dimension of the complex
        for dim in range(1,self.dimension()+1):
            if verbose:
                print "  starting dimension %s" % dim
            if (dim, subcomplex) in self._complex:
                if cochain:
                    differentials[dim-1] = self._complex[(dim, subcomplex)].transpose().change_ring(base_ring)
                    mat = differentials[dim-1]
                else:
                    differentials[dim] = self._complex[(dim, subcomplex)].change_ring(base_ring)
                    mat = differentials[dim]
                if verbose:
                    print "    boundary matrix (cached): it's %s by %s." % (mat.nrows(), mat.ncols())
            else:
                # 'current' is the list of cells in dimension n
                #
                # 'old' is a dictionary, with keys the cells in the
                # previous dimension, values the integers 0, 1, 2,
                # ... (the index of the face).  finding an entry in a
                # dictionary seems to be faster than finding the index
                # of an entry in a list.
                old = dict(zip(current, range(len(current))))
                current = list(self.n_cubes(dim, subcomplex=subcomplex))
                # construct matrix.  it is easiest to construct it as
                # a sparse matrix, specifying which entries are
                # nonzero via a dictionary.
                matrix_data = {}
                col = 0
                if len(old) and len(current):
                    for cube in current:
                        faces = cube.faces_as_pairs()
                        sign = 1
                        for (upper, lower) in faces:
                            try:
                                matrix_data[(old[upper], col)] = sign
                                sign *= -1
                                matrix_data[(old[lower], col)] = sign
                            except KeyError:
                                pass
                        col += 1
                mat = matrix(ZZ, len(old), len(current), matrix_data)
                self._complex[(dim, subcomplex)] = mat
                if cochain:
                    differentials[dim-1] = mat.transpose().change_ring(base_ring)
                else:
                    differentials[dim] = mat.change_ring(base_ring)
                if verbose:
                    print "    boundary matrix computed: it's %s by %s." % (mat.nrows(), mat.ncols())
        # finally, return the chain complex
        if cochain:
            return ChainComplex(data=differentials, base_ring=base_ring,
                                degree=1, check=check_diffs)
        else:
            return ChainComplex(data=differentials, base_ring=base_ring,
                                degree=-1, check=check_diffs)

    def n_skeleton(self, n):
        r"""
        The n-skeleton of this cubical complex.

        :param n: dimension
        :type n: non-negative integer
        :return: cubical complex

        EXAMPLES::

            sage: S2 = cubical_complexes.Sphere(2)
            sage: C3 = cubical_complexes.Cube(3)
            sage: S2 == C3.n_skeleton(2)
            True
        """
        if n >= self.dimension():
            return self
        else:
            data = []
            for d in range(n+1):
                data.extend(list(self.cells()[d]))
            return CubicalComplex(data)

    def graph(self):
        """
        The 1-skeleton of this cubical complex, as a graph.

        EXAMPLES::

            sage: cubical_complexes.Sphere(2).graph()
            Graph on 8 vertices
        """
        data = {}
        vertex_dict = {}
        i = 0
        for vertex in self.n_cells(0):
            vertex_dict[vertex] = i
            data[i] = []
            i += 1
        for edge in self.n_cells(1):
            start = edge.face(0, False)
            end = edge.face(0, True)
            data[vertex_dict[start]].append(vertex_dict[end])
        return Graph(data)

    def is_pure(self):
        """
        True iff this cubical complex is pure: that is,
        all of its maximal faces have the same dimension.

        .. warning::

           This may give the wrong answer if the cubical complex
           was constructed with ``maximality_check`` set to False.

        EXAMPLES::

            sage: S4 = cubical_complexes.Sphere(4)
            sage: S4.is_pure()
            True
            sage: C = CubicalComplex([([0,0], [3,3]), ([1,2], [4,5])])
            sage: C.is_pure()
            False
        """
        dims = [face.dimension() for face in self._facets]
        return max(dims) == min(dims)

    def join(self, other):
        r"""
        The join of this cubical complex with another one.

        NOT IMPLEMENTED.

        :param other: another cubical complex

        EXAMPLES::

            sage: C1 = cubical_complexes.Cube(1)
            sage: C1.join(C1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Joins are not implemented for cubical complexes.
        """
        raise NotImplementedError, "Joins are not implemented for cubical complexes."

    # Use * to mean 'join':
    # __mul__ = join

    def cone(self):
        r"""
        The cone on this cubical complex.

        NOT IMPLEMENTED

        The cone is the complex formed by taking the join of the
        original complex with a one-point complex (that is, a
        0-dimensional cube).  Since joins are not implemented for
        cubical complexes, neither are cones.

        EXAMPLES::

            sage: C1 = cubical_complexes.Cube(1)
            sage: C1.cone()
            Traceback (most recent call last):
            ...
            NotImplementedError: Cones are not implemented for cubical complexes.
        """
        #return self.join(cubical_complexes.Cube(0))
        raise NotImplementedError, "Cones are not implemented for cubical complexes."

    def suspension(self, n=1):
        r"""
        The suspension of this cubical complex.

        NOT IMPLEMENTED

        :param n: suspend this many times
        :type n: positive integer; optional, default 1

        The suspension is the complex formed by taking the join of the
        original complex with a two-point complex (the 0-sphere).
        Since joins are not implemented for cubical complexes, neither
        are suspensions.

        EXAMPLES::

            sage: C1 = cubical_complexes.Cube(1)
            sage: C1.suspension()
            Traceback (most recent call last):
            ...
            NotImplementedError: Suspensions are not implemented for cubical complexes.
        """
#         if n<0:
#             raise ValueError, "n must be non-negative."
#         if n==0:
#             return self
#         if n==1:
#             return self.join(cubical_complexes.Sphere(0))
#         return self.suspension().suspension(int(n-1))
        raise NotImplementedError, "Suspensions are not implemented for cubical complexes."

    def product(self, other):
        r"""
        The product of this cubical complex with another one.

        :param other: another cubical complex

        EXAMPLES::

            sage: RP2 = cubical_complexes.RealProjectivePlane()
            sage: S1 = cubical_complexes.Sphere(1)
            sage: RP2.product(S1).homology()[1]
            Z x C2
        """
        facets = []
        for f in self._facets:
            for g in other._facets:
                facets.append(f.product(g))
        return CubicalComplex(facets)

    def disjoint_union(self, other):
        """
        The disjoint union of this cubical complex with another one.

        :param right: the other cubical complex (the right-hand factor)

        Algorithm: first embed both complexes in d-dimensional
        Euclidean space.  Then embed in (1+d)-dimensional space,
        calling the new axis `x`, and putting the first complex at
        `x=0`, the second at `x=1`.

        EXAMPLES::

            sage: S1 = cubical_complexes.Sphere(1)
            sage: S2 = cubical_complexes.Sphere(2)
            sage: S1.disjoint_union(S2).homology()
            {0: Z, 1: Z, 2: Z}
        """
        embedded_left = len(tuple(self.maximal_cells()[0]))
        embedded_right = len(tuple(other.maximal_cells()[0]))
        zero = [0] * max(embedded_left, embedded_right)
        facets = []
        for f in self.maximal_cells():
            facets.append(Cube([[0,0]]).product(f._translate(zero)))
        for f in other.maximal_cells():
            facets.append(Cube([[1,1]]).product(f._translate(zero)))
        return CubicalComplex(facets)

    def wedge(self, other):
        """
        The wedge (one-point union) of this cubical complex with
        another one.

        :param right: the other cubical complex (the right-hand factor)

        Algorithm: if ``self`` is embedded in `d` dimensions and
        ``other`` in `n` dimensions, embed them in `d+n` dimensions:
        ``self`` using the first `d` coordinates, ``other`` using the
        last `n`, translating them so that they have the origin as a
        common vertex.

        .. note::

            This operation is not well-defined if ``self`` or
            ``other`` is not path-connected.

        EXAMPLES::

            sage: S1 = cubical_complexes.Sphere(1)
            sage: S2 = cubical_complexes.Sphere(2)
            sage: S1.wedge(S2).homology()
            {0: 0, 1: Z, 2: Z}
        """
        embedded_left = len(tuple(self.maximal_cells()[0]))
        embedded_right = len(tuple(other.maximal_cells()[0]))
        translate_left = [-a[0] for a in self.maximal_cells()[0]] + [0] * embedded_right
        translate_right = [-a[0] for a in other.maximal_cells()[0]]
        point_right = Cube([[0,0]] * embedded_left)

        facets = []
        for f in self.maximal_cells():
            facets.append(f._translate(translate_left))
        for f in other.maximal_cells():
            facets.append(point_right.product(f._translate(translate_right)))
        return CubicalComplex(facets)

    def connected_sum(self, other):
        """
        Return the connected sum of self with other.

        :param other: another cubical complex
        :return: the connected sum ``self # other``

        .. warning::

           This does not check that self and other are manifolds, only
           that their facets all have the same dimension.  Since a
           (more or less) random facet is chosen from each complex and
           then glued together, this method may return random
           results if applied to non-manifolds, depending on which
           facet is chosen.

        EXAMPLES::

            sage: T = cubical_complexes.Torus()
            sage: S2 = cubical_complexes.Sphere(2)
            sage: T.connected_sum(S2).cohomology() == T.cohomology()
            True
            sage: RP2 = cubical_complexes.RealProjectivePlane()
            sage: T.connected_sum(RP2).homology(1)
            Z x Z x C2
            sage: RP2.connected_sum(RP2).connected_sum(RP2).homology(1)
            Z x Z x C2
        """
        # connected_sum: first check whether the complexes are pure
        # and of the same dimension.  Then insert degenerate intervals
        # and translate them so that they have a common cube C.  Add one
        # more dimension, embedding the first complex as (..., 0) and
        # the second as (..., 1).  Keep all of the other facets, but remove
        # C x 0 and C x 1, putting in its place (its boundary) x (0,1).
        if not (self.is_pure() and other.is_pure() and
                self.dimension() == other.dimension()):
            raise ValueError, "Complexes are not pure of the same dimension."

        self_facets = list(self.maximal_cells())
        other_facets = list(other.maximal_cells())

        C1 = self_facets.pop()
        C2 = other_facets.pop()
        (insert_self, insert_other, translate) = C1._compare_for_gluing(C2)

        CL = list(C1.tuple())
        for (idx, L) in insert_self:
            CL[idx:idx] = L
        removed = Cube(CL)

        # start assembling the facets in the connected sum: first, the
        # cylinder on the removed face.
        new_facets = []
        cylinder = removed.product(Cube([[0,1]]))
        # don't want to include the ends of the cylinder, so don't
        # include the last pair of faces.  therefore, choose faces up
        # to removed.dimension(), not cylinder.dimension().
        for n in range(removed.dimension()):
            new_facets.append(cylinder.face(n, upper=False))
            new_facets.append(cylinder.face(n, upper=True))

        for cube in self_facets:
            CL = list(cube.tuple())
            for (idx, L) in insert_self:
                CL[idx:idx] = L
            CL.append((0,0))
            new_facets.append(Cube(CL))
        for cube in other_facets:
            CL = list(cube.tuple())
            for (idx, L) in insert_other:
                CL[idx:idx] = L
            CL.append((1,1))
            new_facets.append(Cube(CL)._translate(translate))
        return CubicalComplex(new_facets)

    def _translate(self, vec):
        """
        Translate ``self`` by ``vec``.

        :param vec: anything which can be converted to a tuple of integers
        :return: the translation of ``self`` by ``vec``
        :rtype: cubical complex

        If ``vec`` is shorter than the list of intervals forming the
        complex, pad with zeroes, and similarly if the complexes
        defining tuples are too short.

        EXAMPLES::

            sage: C1 = cubical_complexes.Cube(1)
            sage: C1.maximal_cells()
            {[0,1]}
            sage: C1._translate([2,6]).maximal_cells()
            {[2,3] x [6,6]}
        """
        return CubicalComplex([f._translate(vec) for f in self.maximal_cells()])

    def _chomp_repr_(self):
        r"""
        String representation of self suitable for use by the CHomP
        program.  This lists each maximal cube on its own line.

        EXAMPLES::

            sage: C = cubical_complexes.Cube(0).product(cubical_complexes.Cube(2))
            sage: C.maximal_cells()
            {[0,0] x [0,1] x [0,1]}
            sage: C._chomp_repr_()
            '[0,0] x [0,1] x [0,1]\n'
        """
        s = ""
        for c in self.maximal_cells():
            s += str(c)
            s += "\n"
        return s

    def _simplicial_(self):
        r"""
        Simplicial complex constructed from self.

        ALGORITHM:

        This is constructed as described by Hetyei: choose a total
        ordering of the vertices of the cubical complex.  Then for
        each maximal face

        .. math::

           C = [i_1, j_1] \times [i_2, j_2] \times ... \times [i_k, j_k]

        let `v_1` be the "upper" corner of `C`: `v` is the point
        `(j_1, ..., j_k)`.  Choose a coordinate `n` where the interval
        `[i_n, j_n]` is non-degenerate and form `v_2` by replacing
        `j_n` by `i_n`; repeat to define `v_3`, etc.  The last vertex
        so defined will be `(i_1, ..., i_k)`.  These vertices define a
        simplex, and do the vertices obtained by making different
        choices at each stage.  Thus each `n`-cube is subdivided into
        `n!` simplices.

        REFERENCES:

        - G. Hetyei, "On the Stanley ring of a cubical complex",
          Discrete Comput. Geom. 14 (1995), 305-330.

        EXAMPLES::

            sage: T = cubical_complexes.Torus(); T
            Cubical complex with 16 vertices and 64 cubes
            sage: len(T.maximal_cells())
            16

        When this is triangulated, each maximal 2-dimensional cube
        gets turned into a pair of triangles.  Since there were 16
        maximal cubes, this results in 32 facets in the simplicial
        complex::

            sage: Ts = T._simplicial_(); Ts
            Simplicial complex with 16 vertices and 32 facets
            sage: T.homology() == Ts.homology()
            True

        Each `n`-dimensional cube produces `n!` `n`-simplices::

            sage: S4 = cubical_complexes.Sphere(4)
            sage: len(S4.maximal_cells())
            10
            sage: SimplicialComplex(S4) # calls S4._simplicial_()
            Simplicial complex with 32 vertices and 240 facets
        """
        from sage.homology.simplicial_complex import SimplicialComplex
        simplices = []
        for C in self.maximal_cells():
            simplices.extend(C._triangulation_())
        return SimplicialComplex(simplices)

    def _string_constants(self):
        """
        Tuple containing the name of the type of complex, and the
        singular and plural of the name of the cells from which it is
        built.  This is used in constructing the string representation.

        EXAMPLES::

            sage: S3 = cubical_complexes.Sphere(3)
            sage: S3._string_constants()
            ('Cubical', 'cube', 'cubes')
            sage: S3._repr_()  # indirect doctest
            'Cubical complex with 16 vertices and 80 cubes'
        """
        return ('Cubical', 'cube', 'cubes')


class CubicalComplexExamples():
    r"""
    Some examples of cubical complexes.

    Here are the available examples; you can also type
    "cubical_complexes."  and hit TAB to get a list::

        Sphere
        Torus
        RealProjectivePlane
        KleinBottle
        SurfaceOfGenus
        Cube

    EXAMPLES::

        sage: cubical_complexes.Torus()  # indirect doctest
        Cubical complex with 16 vertices and 64 cubes
        sage: cubical_complexes.Cube(7)
        Cubical complex with 128 vertices and 2187 cubes
        sage: cubical_complexes.Sphere(7)
        Cubical complex with 256 vertices and 6560 cubes
    """

    def Sphere(self,n):
        r"""
        A cubical complex representation of the `n`-dimensional sphere,
        formed by taking the boundary of an `(n+1)`-dimensional cube.

        :param n: the dimension of the sphere
        :type n: non-negative integer

        EXAMPLES::

            sage: cubical_complexes.Sphere(7)
            Cubical complex with 256 vertices and 6560 cubes
        """
        return CubicalComplex(Cube([[0,1]]*(n+1)).faces())

    def Torus(self):
        r"""
        A cubical complex representation of the torus, obtained by
        taking the product of the circle with itself.

        EXAMPLES::

            sage: cubical_complexes.Torus()
            Cubical complex with 16 vertices and 64 cubes
        """
        S1 = cubical_complexes.Sphere(1)
        return S1.product(S1)

    def RealProjectivePlane(self):
        r"""
        A cubical complex representation of the real projective plane.
        This is taken from the examples from CHomP, the Computational
        Homology Project: http://chomp.rutgers.edu/.

        EXAMPLES::

            sage: cubical_complexes.RealProjectivePlane()
            Cubical complex with 21 vertices and 81 cubes
        """
        return CubicalComplex([
            ([0, 1], [0], [0], [0, 1], [0]),
            ([0, 1], [0], [0], [0], [0, 1]),
            ([0], [0, 1], [0, 1], [0], [0]),
            ([0], [0, 1], [0], [0, 1], [0]),
            ([0], [0], [0, 1], [0], [0, 1]),
            ([0, 1], [0, 1], [1], [0], [0]),
            ([0, 1], [1], [0, 1], [0], [0]),
            ([1], [0, 1], [0, 1], [0], [0]),
            ([0, 1], [0, 1], [0], [0], [1]),
            ([0, 1], [1], [0], [0], [0, 1]),
            ([1], [0, 1], [0], [0], [0, 1]),
            ([0, 1], [0], [0, 1], [1], [0]),
            ([0, 1], [0], [1], [0, 1], [0]),
            ([1], [0], [0, 1], [0, 1], [0]),
            ([0], [0, 1], [0], [0, 1], [1]),
            ([0], [0, 1], [0], [1], [0, 1]),
            ([0], [1], [0], [0, 1], [0, 1]),
            ([0], [0], [0, 1], [0, 1], [1]),
            ([0], [0], [0, 1], [1], [0, 1]),
            ([0], [0], [1], [0, 1], [0, 1])])

    def KleinBottle(self):
        r"""
        A cubical complex representation of the Klein bottle, formed
        by taking the connected sum of the real projective plane with
        itself.

        EXAMPLES::

            sage: cubical_complexes.KleinBottle()
            Cubical complex with 42 vertices and 168 cubes
        """
        RP2 = cubical_complexes.RealProjectivePlane()
        return RP2.connected_sum(RP2)

    def SurfaceOfGenus(self, g, orientable=True):
        """
        A surface of genus g as a cubical complex.

        :param g: the genus
        :type g: non-negative integer
        :param orientable: whether the surface should be orientable
        :type orientable: bool, optional, default True

        In the orientable case, return a sphere if `g` is zero, and
        otherwise return a `g`-fold connected sum of a torus with
        itself.

        In the non-orientable case, raise an error if `g` is zero.  If
        `g` is positive, return a `g`-fold connected sum of a
        real projective plane with itself.

        EXAMPLES::

            sage: cubical_complexes.SurfaceOfGenus(2)
            Cubical complex with 32 vertices and 134 cubes
            sage: cubical_complexes.SurfaceOfGenus(1, orientable=False)
            Cubical complex with 21 vertices and 81 cubes
        """
        try:
            g = Integer(g)
        except TypeError:
            raise ValueError("genus must be a non-negative integer")
        if g < 0:
            raise ValueError("genus must be a non-negative integer")
        if g == 0:
            if not orientable:
                raise ValueError("no non-orientable surface of genus zero")
            else:
                return cubical_complexes.Sphere(2)
        if orientable:
            T = cubical_complexes.Torus()
        else:
            T = cubical_complexes.RealProjectivePlane()
        S = T
        for i in range(g-1):
            S = S.connected_sum(T)
        return S

    def Cube(self, n):
        r"""
        A cubical complex representation of an `n`-dimensional cube.

        :param n: the dimension
        :type n: non-negative integer

        EXAMPLES::

            sage: cubical_complexes.Cube(0)
            Cubical complex with 1 vertex and 1 cube
            sage: cubical_complexes.Cube(3)
            Cubical complex with 8 vertices and 27 cubes
        """
        if n == 0:
            return CubicalComplex([Cube([[0]])])
        else:
            return CubicalComplex([Cube([[0,1]]*n)])

cubical_complexes = CubicalComplexExamples()
