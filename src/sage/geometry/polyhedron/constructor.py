r"""
Polyhedra

In this module, a polyhedron is a convex (possibly unbounded) set in
Euclidean space cut out by a finite set of linear inequalities and
linear equations. Note that the dimension of the polyhedron can be
less than the dimension of the ambient space. There are two
complementary representations of the same data:

**H(alf-space/Hyperplane)-representation**
    This describes a polyhedron as the common solution set of a
    finite number of

    * linear **inequalities** `A \vec{x} + b \geq 0`, and

    * linear **equations**  `C \vec{x} + d = 0`.


**V(ertex)-representation**
    The other representation is as the convex hull of vertices (and
    rays and lines to all for unbounded polyhedra) as generators. The
    polyhedron is then the Minkowski sum

    .. MATH::

        P = \text{conv}\{v_1,\dots,v_k\} +
        \sum_{i=1}^m \RR_+ r_i +
        \sum_{j=1}^n \RR \ell_j

    where

    * **vertices** `v_1`, `\dots`, `v_k` are a finite number of
      points. Each vertex is specified by an arbitrary vector, and two
      points are equal if and only if the vector is the same.

    * **rays** `r_1`, `\dots`, `r_m` are a finite number of directions
      (directions of infinity). Each ray is specified by a non-zero
      vector, and two rays are equal if and only if the vectors are
      the same up to rescaling with a positive constant.

    * **lines** `\ell_1`, `\dots`, `\ell_n` are a finite number of
      unoriented directions. In other words, a line is equivalent to
      the set `\{r, -r\}` for a ray `r`. Each line is specified by a
      non-zero vector, and two lines are equivalent if and only if the
      vectors are the same up to rescaling with a non-zero (possibly
      negative) constant.

When specifying a polyhedron, you can input a non-minimal set of
inequalities/equations or generating vertices/rays/lines. The
non-minimal generators are usually called points, non-extremal rays,
and non-extremal lines, but for our purposes it is more convenient to
always talk about vertices/rays/lines. Sage will remove any
superfluous representation objects and always return a minimal
representation. For example, `(0,0)` is a superfluous vertex here::

    sage: triangle = Polyhedron(vertices=[(0,2), (-1,0), (1,0), (0,0)])
    sage: triangle.vertices()
    (A vertex at (-1, 0), A vertex at (1, 0), A vertex at (0, 2))


Unbounded Polyhedra
-------------------

A polytope is defined as a bounded polyhedron. In this case, the
minimal representation is unique and a vertex of the minimal
representation is equivalent to a 0-dimensional face of the
polytope. This is why one generally does not distinguish vertices and
0-dimensional faces. But for non-bounded polyhedra we have to allow
for a more general notion of "vertex" in order to make sense of the
Minkowsi sum presentation::

    sage: half_plane = Polyhedron(ieqs=[(0,1,0)])
    sage: half_plane.Hrepresentation()
    (An inequality (1, 0) x + 0 >= 0,)
    sage: half_plane.Vrepresentation()
    (A line in the direction (0, 1), A ray in the direction (1, 0), A vertex at (0, 0))

Note how we need a point in the above example to anchor the ray and
line. But any point on the boundary of the half-plane would serve the
purpose just as well. Sage picked the origin here, but this choice is
not unique. Similarly, the choice of ray is arbitrary but necessary to
generate the half-plane.

Finally, note that while rays and lines generate unbounded edges of
the polyhedron they are not in a one-to-one correspondence with
them. For example, the infinite strip has two infinite edges (1-faces)
but only one generating line::

    sage: strip = Polyhedron(vertices=[(1,0),(-1,0)], lines=[(0,1)])
    sage: strip.Hrepresentation()
    (An inequality (1, 0) x + 1 >= 0, An inequality (-1, 0) x + 1 >= 0)
    sage: strip.lines()
    (A line in the direction (0, 1),)
    sage: strip.faces(1)
    (<0,1>, <0,2>)
    sage: for face in strip.faces(1):
    ....:      print("{} = {}".format(face, face.as_polyhedron().Vrepresentation()))
    <0,1> = (A line in the direction (0, 1), A vertex at (-1, 0))
    <0,2> = (A line in the direction (0, 1), A vertex at (1, 0))

EXAMPLES::

    sage: trunc_quadr = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
    sage: trunc_quadr
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 2 rays
    sage: v = next(trunc_quadr.vertex_generator())  # the first vertex in the internal enumeration
    sage: v
    A vertex at (0, 1)
    sage: v.vector()
    (0, 1)
    sage: list(v)
    [0, 1]
    sage: len(v)
    2
    sage: v[0] + v[1]
    1
    sage: v.is_vertex()
    True
    sage: type(v)
    <class 'sage.geometry.polyhedron.representation.Vertex'>
    sage: type( v() )
    <type 'sage.modules.vector_integer_dense.Vector_integer_dense'>
    sage: v.polyhedron()
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 2 rays
    sage: r = next(trunc_quadr.ray_generator())
    sage: r
    A ray in the direction (0, 1)
    sage: r.vector()
    (0, 1)
    sage: list( v.neighbors() )
    [A ray in the direction (0, 1), A vertex at (1, 0)]

Inequalities `A \vec{x} + b \geq 0` (and, similarly, equations) are
specified by a list ``[b, A]``::

    sage: Polyhedron(ieqs=[(0,1,0),(0,0,1),(1,-1,-1)]).Hrepresentation()
    (An inequality (-1, -1) x + 1 >= 0,
     An inequality (1, 0) x + 0 >= 0,
     An inequality (0, 1) x + 0 >= 0)

See :func:`Polyhedron` for a detailed description of all possible ways
to construct a polyhedron.


Base Rings
----------

The base ring of the polyhedron can be specified by the ``base_ring``
optional keyword argument. If not specified, a suitable common base
ring for all coordinates/coefficients will be chosen
automatically. Important cases are:

* ``base_ring=QQ`` uses a fast implementation for exact rational
  numbers.

* ``base_ring=ZZ`` is similar to ``QQ``, but the resulting polyhedron
  object will have extra methods for lattice polyhedra.

* ``base_ring=RDF`` uses floating point numbers, this is fast but
  susceptible to numerical errors.

Polyhedra with symmetries often are defined over some algebraic field
extension of the rationals. As a simple example, consider the
equilateral triangle whose vertex coordinates involve `\sqrt{3}`. An
exact way to work with roots in Sage is the :mod:`Algebraic Real Field
<sage.rings.qqbar>` ::

    sage: triangle = Polyhedron([(0,0), (1,0), (1/2, sqrt(3)/2)], base_ring=AA)
    sage: triangle.Hrepresentation()
    (An inequality (-1, -0.5773502691896258?) x + 1 >= 0, 
     An inequality (1, -0.5773502691896258?) x + 0 >= 0, 
     An inequality (0, 1.154700538379252?) x + 0 >= 0)

Without specifying the ``base_ring``, the ``sqrt(3)`` would be a
symbolic ring element and, therefore, the polyhedron defined over the
symbolic ring. This is possible as well, but rather slow::

    sage: Polyhedron([(0,0), (1,0), (1/2, sqrt(3)/2)])
    A 2-dimensional polyhedron in (Symbolic Ring)^2 defined as the convex 
    hull of 3 vertices

Even faster than all algebraic real numbers (the field ``AA``) is
to take the smallest extension field. For the equilateral
triangle, that would be::

    sage: K.<sqrt3> = NumberField(x^2-3)
    sage: Polyhedron([(0,0), (1,0), (1/2, sqrt3/2)])
    A 2-dimensional polyhedron in (Number Field in sqrt3 with defining 
    polynomial x^2 - 3)^2 defined as the convex hull of 3 vertices


Appendix
--------

REFERENCES:

    Komei Fukuda's `FAQ in Polyhedral Computation
    <http://www.ifor.math.ethz.ch/~fukuda/polyfaq/polyfaq.html>`_

AUTHORS:

    - Marshall Hampton: first version, bug fixes, and various improvements, 2008 and 2009
    - Arnaud Bergeron: improvements to triangulation and rendering, 2008
    - Sebastien Barthelemy: documentation improvements, 2008
    - Volker Braun: refactoring, handle non-compact case, 2009 and 2010
    - Andrey Novoseltsev: added Hasse_diagram_from_incidences, 2010
    - Volker Braun: rewrite to use PPL instead of cddlib, 2011
    - Volker Braun: Add support for arbitrary subfields of the reals
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################
from __future__ import print_function
from __future__ import absolute_import

from sage.rings.all import QQ, ZZ, RDF, RR
from sage.misc.decorators import rename_keyword

from .misc import _make_listlist, _common_length_of


#########################################################################
@rename_keyword(deprecation=11634, field='base_ring')
def Polyhedron(vertices=None, rays=None, lines=None,
               ieqs=None, eqns=None,
               ambient_dim=None, base_ring=None, minimize=True, verbose=False,
               backend=None):
    """
    Construct a polyhedron object.

    You may either define it with vertex/ray/line or
    inequalities/equations data, but not both. Redundant data will
    automatically be removed (unless ``minimize=False``), and the
    complementary representation will be computed.

    INPUT:

    - ``vertices`` -- list of point. Each point can be specified as
      any iterable container of ``base_ring`` elements. If ``rays`` or
      ``lines`` are specified but no ``vertices``, the origin is
      taken to be the single vertex.

    - ``rays`` -- list of rays. Each ray can be specified as any
      iterable container of ``base_ring`` elements.

    - ``lines`` -- list of lines. Each line can be specified as any
      iterable container of ``base_ring`` elements.

    - ``ieqs`` -- list of inequalities. Each line can be specified as any
      iterable container of ``base_ring`` elements. An entry equal to
      ``[-1,7,3,4]`` represents the inequality `7x_1+3x_2+4x_3\geq 1`.

    - ``eqns`` -- list of equalities. Each line can be specified as
      any iterable container of ``base_ring`` elements. An entry equal to
      ``[-1,7,3,4]`` represents the equality `7x_1+3x_2+4x_3= 1`.

    - ``base_ring`` -- a sub-field of the reals implemented in
      Sage. The field over which the polyhedron will be defined. For
      ``QQ`` and algebraic extensions, exact arithmetic will be
      used. For ``RDF``, floating point numbers will be used. Floating
      point arithmetic is faster but might give the wrong result for
      degenerate input.

    - ``ambient_dim`` -- integer. The ambient space dimension. Usually
      can be figured out automatically from the H/Vrepresentation
      dimensions.

    - ``backend`` -- string or ``None`` (default). The backend to use. Valid choices are

      * ``'cdd'``: use cdd
        (:mod:`~sage.geometry.polyhedron.backend_cdd`) with `\QQ` or
        `\RDF` coefficients depending on ``base_ring``.

      * ``'ppl'``: use ppl
        (:mod:`~sage.geometry.polyhedron.backend_ppl`) with `\ZZ` or
        `\QQ` coefficients depending on ``base_ring``.

      * ``'field'``: use python implementation
        (:mod:`~sage.geometry.polyhedron.backend_field`) for any field

    Some backends support further optional arguments:

    - ``minimize`` -- boolean (default: ``True``). Whether to
      immediately remove redundant H/V-representation data. Currently
      not used.

    - ``verbose`` -- boolean (default: ``False``). Whether to print
      verbose output for debugging purposes. Only supported by the cdd
      backends.

    OUTPUT:

    The polyhedron defined by the input data.

    EXAMPLES:

    Construct some polyhedra::

        sage: square_from_vertices = Polyhedron(vertices = [[1, 1], [1, -1], [-1, 1], [-1, -1]])
        sage: square_from_ieqs = Polyhedron(ieqs = [[1, 0, 1], [1, 1, 0], [1, 0, -1], [1, -1, 0]])
        sage: list(square_from_ieqs.vertex_generator())
        [A vertex at (1, -1),
         A vertex at (1, 1),
         A vertex at (-1, 1),
         A vertex at (-1, -1)]
        sage: list(square_from_vertices.inequality_generator())
        [An inequality (1, 0) x + 1 >= 0,
         An inequality (0, 1) x + 1 >= 0,
         An inequality (-1, 0) x + 1 >= 0,
         An inequality (0, -1) x + 1 >= 0]
        sage: p = Polyhedron(vertices = [[1.1, 2.2], [3.3, 4.4]], base_ring=RDF)
        sage: p.n_inequalities()
        2

    The same polyhedron given in two ways::

        sage: p = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0]])
        sage: p.Vrepresentation()
        (A line in the direction (0, 0, 1),
         A ray in the direction (1, 0, 0),
         A ray in the direction (0, 1, 0),
         A vertex at (0, 0, 0))
        sage: q = Polyhedron(vertices=[[0,0,0]], rays=[[1,0,0],[0,1,0]], lines=[[0,0,1]])
        sage: q.Hrepresentation()
        (An inequality (1, 0, 0) x + 0 >= 0,
         An inequality (0, 1, 0) x + 0 >= 0)

    Finally, a more complicated example. Take `\mathbb{R}_{\geq 0}^6` with
    coordinates `a, b, \dots, f` and

      * The inequality `e+b \geq c+d`
      * The inequality `e+c \geq b+d`
      * The equation `a+b+c+d+e+f = 31`

    ::

        sage: positive_coords = Polyhedron(ieqs=[
        ...       [0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0],
        ...       [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])
        sage: P = Polyhedron(ieqs=positive_coords.inequalities() + (
        ...       [0,0,1,-1,-1,1,0], [0,0,-1,1,-1,1,0]), eqns=[[-31,1,1,1,1,1,1]])
        sage: P
        A 5-dimensional polyhedron in QQ^6 defined as the convex hull of 7 vertices
        sage: P.dim()
        5
        sage: P.Vrepresentation()
        (A vertex at (31, 0, 0, 0, 0, 0), A vertex at (0, 0, 0, 0, 0, 31),
         A vertex at (0, 0, 0, 0, 31, 0), A vertex at (0, 0, 31/2, 0, 31/2, 0),
         A vertex at (0, 31/2, 31/2, 0, 0, 0), A vertex at (0, 31/2, 0, 0, 31/2, 0),
         A vertex at (0, 0, 0, 31/2, 31/2, 0))

    .. NOTE::

      * Once constructed, a ``Polyhedron`` object is immutable.

      * Although the option ``field=RDF`` allows numerical data to
        be used, it might not give the right answer for degenerate
        input data - the results can depend upon the tolerance
        setting of cdd.
    """
    # Clean up the arguments
    vertices = _make_listlist(vertices)
    rays     = _make_listlist(rays)
    lines    = _make_listlist(lines)
    ieqs     = _make_listlist(ieqs)
    eqns     = _make_listlist(eqns)

    got_Vrep = (len(vertices+rays+lines) > 0)
    got_Hrep = (len(ieqs+eqns) > 0)

    if got_Vrep and got_Hrep:
        raise ValueError('You cannot specify both H- and V-representation.')
    elif got_Vrep:
        deduced_ambient_dim = _common_length_of(vertices, rays, lines)[1]
    elif got_Hrep:
        deduced_ambient_dim = _common_length_of(ieqs, eqns)[1] - 1
    else:
        if ambient_dim is None:
            deduced_ambient_dim = 0
        else:
            deduced_ambient_dim = ambient_dim
        if base_ring is None:
            base_ring = ZZ

    # set ambient_dim
    if ambient_dim is not None and deduced_ambient_dim != ambient_dim:
        raise ValueError('Ambient space dimension mismatch. Try removing the "ambient_dim" parameter.')
    ambient_dim = deduced_ambient_dim

    # figure out base_ring
    from sage.misc.flatten import flatten
    values = flatten(vertices+rays+lines+ieqs+eqns)
    if base_ring is not None:
        try:
            convert = not all(x.parent() is base_ring for x in values)
        except AttributeError:   # No x.parent() method?
            convert = True
    else:
        from sage.rings.integer import is_Integer
        from sage.rings.rational import is_Rational
        from sage.rings.real_double import is_RealDoubleElement
        if all(is_Integer(x) for x in values):
            if got_Vrep:
                base_ring = ZZ
            else:   # integral inequalities usually do not determine a lattice polytope!
                base_ring = QQ
            convert = False
        elif all(is_Rational(x) for x in values):
            base_ring = QQ
            convert = False
        elif all(is_RealDoubleElement(x) for x in values):
            base_ring = RDF
            convert = False
        else:
            try:
                for v in values:
                    ZZ(v)
                if got_Vrep:
                    base_ring = ZZ
                else:
                    base_ring = QQ
                convert = True
            except (TypeError, ValueError):
                from sage.structure.sequence import Sequence
                values = Sequence(values)
                common_ring = values.universe()
                if QQ.has_coerce_map_from(common_ring):
                    base_ring = QQ
                    convert = True
                elif common_ring is RR:   # DWIM: replace with RDF
                    base_ring = RDF
                    convert = True
                else:
                    base_ring = common_ring
                    convert = True

    # Add the origin if necesarry
    if got_Vrep and len(vertices)==0:
        vertices = [ [0]*ambient_dim ]

    # Specific backends can override the base_ring
    from sage.geometry.polyhedron.parent import Polyhedra
    parent = Polyhedra(base_ring, ambient_dim, backend=backend)
    base_ring = parent.base_ring()


    # finally, construct the Polyhedron
    Hrep = Vrep = None
    if got_Hrep:
        Hrep = [ieqs, eqns]
    if got_Vrep:
        Vrep = [vertices, rays, lines]
    return parent(Vrep, Hrep, convert=convert, verbose=verbose)
