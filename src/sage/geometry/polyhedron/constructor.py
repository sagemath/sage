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

.. SEEALSO::

    If one only needs to keep track of a system of linear system of
    inequalities, one should also consider the class for mixed integer linear
    programming.

    - :mod:`Mixed Integer Linear Programming <sage.numerical.mip>`


Unbounded Polyhedra
-------------------

A polytope is defined as a bounded polyhedron. In this case, the
minimal representation is unique and a vertex of the minimal
representation is equivalent to a 0-dimensional face of the
polytope. This is why one generally does not distinguish vertices and
0-dimensional faces. But for non-bounded polyhedra we have to allow
for a more general notion of "vertex" in order to make sense of the
Minkowski sum presentation::

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
    sage: [f.ambient_V_indices() for f in strip.faces(1)]
    [(0, 2), (0, 1)]
    sage: for face in strip.faces(1):
    ....:      print(face.ambient_V_indices())
    (0, 2)
    (0, 1)
    sage: for face in strip.faces(1):
    ....:      print("{} = {}".format(face.ambient_V_indices(), face.as_polyhedron().Vrepresentation()))
    (0, 2) = (A line in the direction (0, 1), A vertex at (1, 0))
    (0, 1) = (A line in the direction (0, 1), A vertex at (-1, 0))

EXAMPLES::

    sage: trunc_quadr = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
    sage: trunc_quadr
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 2 rays
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
    <class 'sage.modules.vector_rational_dense.Vector_rational_dense'>
    sage: v.polyhedron()
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 2 rays
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

    sage: triangle = Polyhedron([(0,0), (1,0), (1/2, sqrt(3)/2)], base_ring=AA) # optional - sage.rings.number_field  # optional - sage.symbolic
    sage: triangle.Hrepresentation()                                            # optional - sage.rings.number_field  # optional - sage.symbolic
    (An inequality (-1, -0.5773502691896258?) x + 1 >= 0,
     An inequality (1, -0.5773502691896258?) x + 0 >= 0,
     An inequality (0, 1.154700538379252?) x + 0 >= 0)

Without specifying the ``base_ring``, the ``sqrt(3)`` would be a
symbolic ring element and, therefore, the polyhedron defined over the
symbolic ring. This is currently not supported as SR is not exact::

    sage: Polyhedron([(0,0), (1,0), (1/2, sqrt(3)/2)])                          # optional - sage.symbolic
    Traceback (most recent call last):
    ...
    ValueError: no default backend for computations with Symbolic Ring

    sage: SR.is_exact()                                                         # optional - sage.symbolic
    False

Even faster than all algebraic real numbers (the field ``AA``) is
to take the smallest extension field. For the equilateral
triangle, that would be::

    sage: K.<sqrt3> = NumberField(x^2 - 3, embedding=AA(3)**(1/2))              # optional - sage.rings.number_field
    sage: Polyhedron([(0,0), (1,0), (1/2, sqrt3/2)])                            # optional - sage.rings.number_field
    A 2-dimensional polyhedron in (Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?)^2 defined as the convex hull of 3 vertices

.. WARNING::

    Be careful when you construct polyhedra with floating point numbers. The only
    available backend for such computation is ``cdd`` which uses machine floating
    point numbers which have have limited precision. If the input consists of
    floating point numbers and the ``base_ring`` is not specified, the base ring is
    set to be the ``RealField`` with the precision given by the minimal bit precision
    of the input. Then, if the obtained minimum is 53 bits of precision, the
    constructor converts automatically the base ring to ``RDF``. Otherwise,
    it returns an error::

        sage: Polyhedron(vertices = [[1.12345678901234, 2.12345678901234]])
        A 0-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex
        sage: Polyhedron(vertices = [[1.12345678901234, 2.123456789012345]])
        A 0-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex
        sage: Polyhedron(vertices = [[1.123456789012345, 2.123456789012345]])
        Traceback (most recent call last):
        ...
        ValueError: the only allowed inexact ring is 'RDF' with backend 'cdd'

    The strongly suggested method to input floating point numbers is to specify the
    ``base_ring`` to be ``RDF``::

        sage: Polyhedron(vertices = [[1.123456789012345, 2.123456789012345]], base_ring=RDF)
        A 0-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex

.. SEEALSO::

    :mod:`Parents for polyhedra <sage.geometry.polyhedron.parent.Polyhedra>`

Base classes
------------

Depending on the chosen base ring, a specific class is used to represent the polyhedron object.

.. SEEALSO::

    - :mod:`Base class for polyhedra <sage.geometry.polyhedron.base.Polyhedron_base>`
    - :mod:`Base class for polyhedra over integers <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`
    - :mod:`Base class for polyhedra over rationals <sage.geometry.polyhedron.base_QQ.Polyhedron_QQ>`
    - :mod:`Base class for polyhedra over RDF <sage.geometry.polyhedron.base_RDF.Polyhedron_RDF>`

The most important base class is **Base class for polyhedra** from which other base classes and backends inherit.

Backends
--------

There are different backends available to deal with polyhedron objects.

.. SEEALSO::

    - :mod:`cdd backend for polyhedra <sage.geometry.polyhedron.backend_cdd.Polyhedron_cdd>`
    - :mod:`field backend for polyhedra <sage.geometry.polyhedron.backend_field.Polyhedron_field>`
    - :mod:`normaliz backend for polyhedra <sage.geometry.polyhedron.backend_normaliz.Polyhedron_normaliz>`
    - :mod:`ppl backend for polyhedra <sage.geometry.polyhedron.backend_ppl.Polyhedron_ppl>`

.. NOTE::

    Depending on the backend used, it may occur that different methods are
    available or not.

Appendix
--------

REFERENCES:

    Komei Fukuda's `FAQ in Polyhedral Computation
    <https://www.inf.ethz.ch/personal/fukudak/polyfaq/polyfaq.html>`_

AUTHORS:

    - Marshall Hampton: first version, bug fixes, and various improvements, 2008 and 2009
    - Arnaud Bergeron: improvements to triangulation and rendering, 2008
    - Sebastien Barthelemy: documentation improvements, 2008
    - Volker Braun: refactoring, handle non-compact case, 2009 and 2010
    - Andrey Novoseltsev: added lattice_from_incidences, 2010
    - Volker Braun: rewrite to use PPL instead of cddlib, 2011
    - Volker Braun: Add support for arbitrary subfields of the reals
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RR

from .misc import _make_listlist, _common_length_of


#########################################################################
def Polyhedron(vertices=None, rays=None, lines=None,
               ieqs=None, eqns=None,
               ambient_dim=None, base_ring=None, minimize=True, verbose=False,
               backend=None, mutable=False):
    r"""
    Construct a polyhedron object.

    You may either define it with vertex/ray/line or
    inequalities/equations data, but not both. Redundant data will
    automatically be removed (unless ``minimize=False``), and the
    complementary representation will be computed.

    INPUT:

    - ``vertices`` -- list of points. Each point can be specified as
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
        `\RDF` coefficients depending on ``base_ring``

      * ``'normaliz'``: use normaliz
        (:mod:`~sage.geometry.polyhedron.backend_normaliz`) with `\ZZ` or
        `\QQ` coefficients depending on ``base_ring``

      * ``'polymake'``: use polymake
        (:mod:`~sage.geometry.polyhedron.backend_polymake`) with `\QQ`, `\RDF` or
        ``QuadraticField`` coefficients depending on ``base_ring``

      * ``'ppl'``: use ppl
        (:mod:`~sage.geometry.polyhedron.backend_ppl`) with `\ZZ` or
        `\QQ` coefficients depending on ``base_ring``

      * ``'field'``: use python implementation
        (:mod:`~sage.geometry.polyhedron.backend_field`) for any field

    Some backends support further optional arguments:

    - ``minimize`` -- boolean (default: ``True``); whether to
      immediately remove redundant H/V-representation data; currently
      not used.

    - ``verbose`` -- boolean (default: ``False``); whether to print
      verbose output for debugging purposes; only supported by the cdd and
      normaliz backends

    - ``mutable`` -- boolean (default: ``False``); whether the polyhedron
      is mutable

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
        ....:     [0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0],
        ....:     [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])
        sage: P = Polyhedron(ieqs=positive_coords.inequalities() + (
        ....:     [0,0,1,-1,-1,1,0], [0,0,-1,1,-1,1,0]), eqns=[[-31,1,1,1,1,1,1]])
        sage: P
        A 5-dimensional polyhedron in QQ^6 defined as the convex hull of 7 vertices
        sage: P.dim()
        5
        sage: P.Vrepresentation()
        (A vertex at (31, 0, 0, 0, 0, 0), A vertex at (0, 0, 0, 0, 0, 31),
         A vertex at (0, 0, 0, 0, 31, 0), A vertex at (0, 0, 31/2, 0, 31/2, 0),
         A vertex at (0, 31/2, 31/2, 0, 0, 0), A vertex at (0, 31/2, 0, 0, 31/2, 0),
         A vertex at (0, 0, 0, 31/2, 31/2, 0))

    Regular icosahedron, centered at `0` with edge length `2`, with vertices given
    by the cyclic shifts of `(0, \pm 1, \pm (1+\sqrt(5))/2)`, cf.
    :wikipedia:`Regular_icosahedron`. It needs a number field::

        sage: R0.<r0> = QQ[]                                                    # optional - sage.rings.number_field
        sage: R1.<r1> = NumberField(r0^2-5, embedding=AA(5)**(1/2))             # optional - sage.rings.number_field
        sage: gold = (1+r1)/2                                                   # optional - sage.rings.number_field
        sage: v = [[0, 1, gold], [0, 1, -gold], [0, -1, gold], [0, -1, -gold]]  # optional - sage.rings.number_field
        sage: pp = Permutation((1, 2, 3))          # optional - sage.combinat   # optional - sage.rings.number_field
        sage: icosah = Polyhedron(                 # optional - sage.combinat   # optional - sage.rings.number_field
        ....:    [(pp^2).action(w) for w in v] + [pp.action(w) for w in v] + v,
        ....:    base_ring=R1)
        sage: len(icosah.faces(2))                 # optional - sage.combinat   # optional - sage.rings.number_field
        20

    When the input contains elements of a Number Field, they require an
    embedding::

        sage: K = NumberField(x^2-2,'s')                                        # optional - sage.rings.number_field
        sage: s = K.0                                                           # optional - sage.rings.number_field
        sage: L = NumberField(x^3-2,'t')                                        # optional - sage.rings.number_field
        sage: t = L.0                                                           # optional - sage.rings.number_field
        sage: P = Polyhedron(vertices = [[0,s],[t,0]])                          # optional - sage.rings.number_field
        Traceback (most recent call last):
        ...
        ValueError: invalid base ring

    Create a mutable polyhedron::

        sage: P = Polyhedron(vertices=[[0, 1], [1, 0]], mutable=True)
        sage: P.is_mutable()
        True
        sage: hasattr(P, "_Vrepresentation")
        False
        sage: P.Vrepresentation()
        (A vertex at (0, 1), A vertex at (1, 0))
        sage: hasattr(P, "_Vrepresentation")
        True

    .. NOTE::

      * Once constructed, a ``Polyhedron`` object is immutable.

      * Although the option ``base_ring=RDF`` allows numerical data to
        be used, it might not give the right answer for degenerate
        input data - the results can depend upon the tolerance
        setting of cdd.


    TESTS:

    Check that giving ``float`` input gets converted to ``RDF`` (see :trac:`22605`)::

        sage: f = float(1.1)
        sage: Polyhedron(vertices=[[f]])
        A 0-dimensional polyhedron in RDF^1 defined as the convex hull of 1 vertex

    Check that giving ``int`` input gets converted to ``ZZ`` (see :trac:`22605`)::

        sage: Polyhedron(vertices=[[int(42)]])
        A 0-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex

    Check that giving ``Fraction`` input gets converted to ``QQ`` (see :trac:`22605`)::

        sage: from fractions import Fraction
        sage: f = Fraction(int(6), int(8))
        sage: Polyhedron(vertices=[[f]])
        A 0-dimensional polyhedron in QQ^1 defined as the convex hull of 1 vertex

    Check that non-compact polyhedra given by V-representation have base ring ``QQ``,
    not ``ZZ`` (see :trac:`27840`)::

        sage: Q = Polyhedron(vertices=[(1, 2, 3), (1, 3, 2), (2, 1, 3),
        ....:                          (2, 3, 1), (3, 1, 2), (3, 2, 1)],
        ....:                rays=[[1, 1, 1]], lines=[[1, 2, 3]], backend='ppl')
        sage: Q.base_ring()
        Rational Field

    Check that enforcing base ring `ZZ` for this example gives an error::

        sage: Q = Polyhedron(vertices=[(1, 2, 3), (1, 3, 2), (2, 1, 3),
        ....:                          (2, 3, 1), (3, 1, 2), (3, 2, 1)],
        ....:                rays=[[1, 1, 1]], lines=[[1, 2, 3]], backend='ppl',
        ....:                base_ring=ZZ)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    Check that input with too many bits of precision returns an error (see
    :trac:`22552`)::

        sage: Polyhedron(vertices=[(8.3319544851638732, 7.0567045956967727), (6.4876921900819049, 4.8435898415984129)])
        Traceback (most recent call last):
        ...
        ValueError: the only allowed inexact ring is 'RDF' with backend 'cdd'

    Check that setting ``base_ring`` to a ``RealField`` returns an error (see :trac:`22552`)::

        sage: Polyhedron(vertices =[(8.3, 7.0), (6.4, 4.8)], base_ring=RealField(40))
        Traceback (most recent call last):
        ...
        ValueError: no default backend for computations with Real Field with 40 bits of precision
        sage: Polyhedron(vertices =[(8.3, 7.0), (6.4, 4.8)], base_ring=RealField(53))
        Traceback (most recent call last):
        ...
        ValueError: no default backend for computations with Real Field with 53 bits of precision

    Check that :trac:`17339` is fixed::

        sage: Polyhedron(ambient_dim=0, ieqs=[], eqns=[[1]], base_ring=QQ)
        The empty polyhedron in QQ^0
        sage: P = Polyhedron(ambient_dim=0, ieqs=[], eqns=[], base_ring=QQ); P
        A 0-dimensional polyhedron in QQ^0 defined as the convex hull of 1 vertex
        sage: P.Vrepresentation()
        (A vertex at (),)
        sage: Polyhedron(ambient_dim=2, ieqs=[], eqns=[], base_ring=QQ)
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull
         of 1 vertex and 2 lines
        sage: Polyhedron(ambient_dim=2, ieqs=[], eqns=[], base_ring=QQ, backend='field')
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull
         of 1 vertex and 2 lines
        sage: Polyhedron(ambient_dim=0, ieqs=[], eqns=[[1]], base_ring=QQ, backend="cdd")
        The empty polyhedron in QQ^0
        sage: Polyhedron(ambient_dim=0, ieqs=[], eqns=[[1]], base_ring=QQ, backend="ppl")
        The empty polyhedron in QQ^0
        sage: Polyhedron(ambient_dim=0, ieqs=[], eqns=[[1]], base_ring=QQ, backend="field")
        The empty polyhedron in QQ^0

        sage: Polyhedron(ambient_dim=2, vertices=[], rays=[], lines=[], base_ring=QQ)
        The empty polyhedron in QQ^2

    .. SEEALSO::

        :mod:`Library of polytopes <sage.geometry.polyhedron.library>`
    """
    got_Vrep = not ((vertices is None) and (rays is None) and (lines is None))
    got_Hrep = not ((ieqs is None) and (eqns is None))

    # Clean up the arguments
    vertices = _make_listlist(vertices)
    rays = _make_listlist(rays)
    lines = _make_listlist(lines)
    ieqs = _make_listlist(ieqs)
    eqns = _make_listlist(eqns)

    if got_Vrep and got_Hrep:
        raise ValueError('cannot specify both H- and V-representation.')
    elif got_Vrep:
        deduced_ambient_dim = _common_length_of(vertices, rays, lines)[1]
        if deduced_ambient_dim is None:
            if ambient_dim is not None:
                deduced_ambient_dim = ambient_dim
            else:
                deduced_ambient_dim = 0
    elif got_Hrep:
        deduced_ambient_dim = _common_length_of(ieqs, eqns)[1]
        if deduced_ambient_dim is None:
            if ambient_dim is not None:
                deduced_ambient_dim = ambient_dim
            else:
                deduced_ambient_dim = 0
        else:
            deduced_ambient_dim -= 1
    else:
        if ambient_dim is None:
            deduced_ambient_dim = 0
        else:
            deduced_ambient_dim = ambient_dim
        if base_ring is None:
            base_ring = ZZ

    # set ambient_dim
    if ambient_dim is not None and deduced_ambient_dim != ambient_dim:
        raise ValueError('ambient space dimension mismatch. Try removing the "ambient_dim" parameter.')
    ambient_dim = deduced_ambient_dim

    # figure out base_ring
    from sage.misc.flatten import flatten
    from sage.structure.element import parent
    from sage.categories.fields import Fields
    from sage.categories.rings import Rings

    values = flatten(vertices + rays + lines + ieqs + eqns)
    if base_ring is not None:
        convert = any(parent(x) is not base_ring for x in values)
    elif not values:
        base_ring = ZZ
        convert = False
    else:
        P = parent(values[0])
        if any(parent(x) is not P for x in values):
            from sage.structure.sequence import Sequence
            P = Sequence(values).universe()
            convert = True
        else:
            convert = False

        from sage.structure.coerce import py_scalar_parent
        if isinstance(P, type):
            base_ring = py_scalar_parent(P)
            convert = convert or P is not base_ring
        else:
            base_ring = P

        if base_ring not in Fields():
            got_compact_Vrep = got_Vrep and not rays and not lines
            got_cone_Vrep = got_Vrep and all(all(x == 0 for x in v) for v in vertices)
            if not got_compact_Vrep and not got_cone_Vrep:
                base_ring = base_ring.fraction_field()
                convert = True

        if base_ring not in Rings():
            raise ValueError('invalid base ring')

        try:
            from sage.symbolic.ring import SR
        except ImportError:
            SR = None
        if base_ring is not SR and not base_ring.is_exact():
            # TODO: remove this hack?
            if base_ring is RR:
                base_ring = RDF
                convert = True
            elif base_ring is not RDF:
                raise ValueError("the only allowed inexact ring is 'RDF' with backend 'cdd'")

    # Add the origin if necessary
    if got_Vrep and len(vertices) == 0 and len(rays + lines) > 0:
        vertices = [[0] * ambient_dim]

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
    return parent(Vrep, Hrep, convert=convert, verbose=verbose, mutable=mutable)
