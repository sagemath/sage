r"""
Convex rational polyhedral cones

This module was designed as a part of framework for toric varieties
(:mod:`~sage.schemes.toric.variety`,
:mod:`~sage.schemes.toric.fano_variety`). While the emphasis is on
strictly convex cones, non-strictly convex cones are supported as well. Work
with distinct lattices (in the sense of discrete subgroups spanning vector
spaces) is supported. The default lattice is :class:`ToricLattice
<sage.geometry.toric_lattice.ToricLatticeFactory>` `N` of the appropriate
dimension. The only case when you must specify lattice explicitly is creation
of a 0-dimensional cone, where dimension of the ambient space cannot be
guessed.

AUTHORS:

- Andrey Novoseltsev (2010-05-13): initial version.

- Andrey Novoseltsev (2010-06-17): substantial improvement during review by
  Volker Braun.

- Volker Braun (2010-06-21): various spanned/quotient/dual lattice
  computations added.

- Volker Braun (2010-12-28): Hilbert basis for cones.

- Andrey Novoseltsev (2012-02-23): switch to PointCollection container.

EXAMPLES:

Use :func:`Cone` to construct cones::

    sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
    sage: halfspace = Cone([(1,0,0), (0,1,0), (-1,-1,0), (0,0,1)])
    sage: positive_xy = Cone([(1,0,0), (0,1,0)])
    sage: four_rays = Cone([(1,1,1), (1,-1,1), (-1,-1,1), (-1,1,1)])

For all of the cones above we have provided primitive generating rays, but in
fact this is not necessary - a cone can be constructed from any collection of
rays (from the same space, of course). If there are non-primitive (or even
non-integral) rays, they will be replaced with primitive ones. If there are
extra rays, they will be discarded. Of course, this means that :func:`Cone`
has to do some work before actually constructing the cone and sometimes it is
not desirable, if you know for sure that your input is already "good". In this
case you can use options ``check=False`` to force :func:`Cone` to use
exactly the directions that you have specified and ``normalize=False`` to
force it to use exactly the rays that you have specified. However, it is
better not to use these possibilities without necessity, since cones are
assumed to be represented by a minimal set of primitive generating rays.
See :func:`Cone` for further documentation on construction.

Once you have a cone, you can perform numerous operations on it. The most
important ones are, probably, ray accessing methods::

    sage: rays = halfspace.rays()
    sage: rays
    N( 0,  0, 1),
    N( 0,  1, 0),
    N( 0, -1, 0),
    N( 1,  0, 0),
    N(-1,  0, 0)
    in 3-d lattice N
    sage: rays.set()
    frozenset({N(-1, 0, 0), N(0, -1, 0), N(0, 0, 1), N(0, 1, 0), N(1, 0, 0)})
    sage: rays.matrix()
    [ 0  0  1]
    [ 0  1  0]
    [ 0 -1  0]
    [ 1  0  0]
    [-1  0  0]
    sage: rays.column_matrix()
    [ 0  0  0  1 -1]
    [ 0  1 -1  0  0]
    [ 1  0  0  0  0]
    sage: rays(3)
    N(1, 0, 0)
    in 3-d lattice N
    sage: rays[3]
    N(1, 0, 0)
    sage: halfspace.ray(3)
    N(1, 0, 0)

The method :meth:`~IntegralRayCollection.rays` returns a
:class:`~sage.geometry.point_collection.PointCollection` with the
`i`-th element being the primitive integral generator of the `i`-th
ray. It is possible to convert this collection to a matrix with either
rows or columns corresponding to these generators. You may also change
the default
:meth:`~sage.geometry.point_collection.PointCollection.output_format`
of all point collections to be such a matrix.

If you want to do something with each ray of a cone, you can write ::

    sage: for ray in positive_xy: print(ray)
    N(1, 0, 0)
    N(0, 1, 0)

There are two dimensions associated to each cone - the dimension of the
subspace spanned by the cone and the dimension of the space where it lives::

    sage: positive_xy.dim()
    2
    sage: positive_xy.lattice_dim()
    3

You also may be interested in this dimension::

    sage: dim(positive_xy.linear_subspace())
    0
    sage: dim(halfspace.linear_subspace())
    2

Or, perhaps, all you care about is whether it is zero or not::

    sage: positive_xy.is_strictly_convex()
    True
    sage: halfspace.is_strictly_convex()
    False

You can also perform these checks::

    sage: positive_xy.is_simplicial()
    True
    sage: four_rays.is_simplicial()
    False
    sage: positive_xy.is_smooth()
    True

You can work with subcones that form faces of other cones::

    sage: face = four_rays.faces(dim=2)[0]
    sage: face
    2-d face of 3-d cone in 3-d lattice N
    sage: face.rays()
    N(-1, -1, 1),
    N(-1,  1, 1)
    in 3-d lattice N
    sage: face.ambient_ray_indices()
    (2, 3)
    sage: four_rays.rays(face.ambient_ray_indices())
    N(-1, -1, 1),
    N(-1,  1, 1)
    in 3-d lattice N

If you need to know inclusion relations between faces, you can use ::

    sage: L = four_rays.face_lattice()
    sage: [len(s) for s in L.level_sets()]
    [1, 4, 4, 1]
    sage: face = L.level_sets()[2][0]
    sage: face.rays()
    N(1,  1, 1),
    N(1, -1, 1)
    in 3-d lattice N
    sage: L.hasse_diagram().neighbors_in(face)
    [1-d face of 3-d cone in 3-d lattice N,
     1-d face of 3-d cone in 3-d lattice N]

.. WARNING::

    The order of faces in level sets of
    the face lattice may differ from the order of faces returned by
    :meth:`~ConvexRationalPolyhedralCone.faces`. While the first order is
    random, the latter one ensures that one-dimensional faces are listed in
    the same order as generating rays.

When all the functionality provided by cones is not enough, you may want to
check if you can do necessary things using polyhedra corresponding to cones::

    sage: four_rays.polyhedron()
    A 3-dimensional polyhedron in ZZ^3 defined as
    the convex hull of 1 vertex and 4 rays

And of course you are always welcome to suggest new features that should be
added to cones!

REFERENCES:

- [Ful1993]_
"""

# ****************************************************************************
#       Copyright (C) 2010-2014 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010-2018 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#       Copyright (C) 2012      Christian Stump
#       Copyright (C) 2014-2018 Frédéric Chapoton
#       Copyright (C) 2014      Peter Bruin
#       Copyright (C) 2015-2017 Jori Mäntysalo
#       Copyright (C) 2015-2020 Michael Orlitzky
#       Copyright (C) 2016-2020 John H. Palmieri
#       Copyright (C) 2018      David Coudert
#       Copyright (C) 2019-2020 Jonathan Kliem
#       Copyright (C) 2020-2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Hashable, Iterable, Container
from copy import copy
from warnings import warn

from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.posets.posets', 'FinitePoset')
from sage.arith.misc import GCD as gcd
from sage.arith.functions import lcm
from sage.geometry.point_collection import PointCollection
from sage.geometry.polyhedron.constructor import Polyhedron
lazy_import('sage.geometry.hasse_diagram', 'lattice_from_incidences')
from sage.geometry.toric_lattice import (ToricLattice, ToricLattice_generic,
                                         ToricLattice_quotient)
lazy_import('sage.geometry.toric_plotter', ['ToricPlotter', 'label_list'])
from sage.geometry.relative_interior import RelativeInterior
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.special import column_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.latex import latex
from sage.modules.free_module import span, VectorSpace
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.structure.element import parent
from sage.structure.richcmp import richcmp_method, richcmp
lazy_import('sage.geometry.integral_points', 'parallelotope_points')
from sage.geometry.convex_set import ConvexSet_closed
import sage.geometry.abc

from sage.features import PythonModule
lazy_import('ppl', ['C_Polyhedron', 'Generator_System', 'Constraint_System',
                    'Linear_Expression', 'Poly_Con_Relation'],
                    feature=PythonModule("ppl", spkg='pplpy', type='standard'))
lazy_import('ppl', ['ray', 'point'], as_=['PPL_ray', 'PPL_point'],
                    feature=PythonModule("ppl", spkg='pplpy', type='standard'))


def Cone(rays, lattice=None, check=True, normalize=True):
    r"""
    Construct a (not necessarily strictly) convex rational polyhedral cone.

    INPUT:

    - ``rays`` -- list of rays; each ray should be given as a list
      or a vector convertible to the rational extension of the given
      ``lattice``. May also be specified by a
      :class:`~sage.geometry.polyhedron.base.Polyhedron_base` object;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically;

    - ``check`` -- by default the input data will be checked for
      correctness (e.g. that all rays have the same number of
      components) and generating rays will be constructed from
      ``rays``. If you know that the input is a minimal set of
      generators of a valid cone, you may significantly decrease
      construction time using ``check=False`` option;

    - ``normalize`` -- you can further speed up construction using
      ``normalize=False`` option. In this case ``rays`` must be a list of
      immutable primitive rays in ``lattice``. In general, you should not use
      this option, it is designed for code optimization and does not give as
      drastic improvement in speed as the previous one.

    OUTPUT: convex rational polyhedral cone determined by ``rays``

    EXAMPLES:

    Let's define a cone corresponding to the first quadrant of the plane
    (note, you can even mix objects of different types to represent rays, as
    long as you let this function to perform all the checks and necessary
    conversions!)::

        sage: quadrant = Cone([(1,0), [0,1]])
        sage: quadrant
        2-d cone in 2-d lattice N
        sage: quadrant.rays()
        N(1, 0),
        N(0, 1)
        in 2-d lattice N

    If you give more rays than necessary, the extra ones will be discarded::

        sage: Cone([(1,0), (0,1), (1,1), (0,1)]).rays()
        N(0, 1),
        N(1, 0)
        in 2-d lattice N

    However, this work is not done with ``check=False`` option, so use it
    carefully! ::

        sage: Cone([(1,0), (0,1), (1,1), (0,1)], check=False).rays()
        N(1, 0),
        N(0, 1),
        N(1, 1),
        N(0, 1)
        in 2-d lattice N

    Even worse things can happen with ``normalize=False`` option::

        sage: Cone([(1,0), (0,1)], check=False, normalize=False)
        Traceback (most recent call last):
        ...
        AttributeError: 'tuple' object has no attribute 'parent'...

    You can construct different "not" cones: not full-dimensional, not
    strictly convex, not containing any rays::

        sage: one_dimensional_cone = Cone([(1,0)])
        sage: one_dimensional_cone.dim()
        1
        sage: half_plane = Cone([(1,0), (0,1), (-1,0)])
        sage: half_plane.rays()
        N( 0, 1),
        N( 1, 0),
        N(-1, 0)
        in 2-d lattice N
        sage: half_plane.is_strictly_convex()
        False
        sage: origin = Cone([(0,0)])
        sage: origin.rays()
        Empty collection
        in 2-d lattice N
        sage: origin.dim()
        0
        sage: origin.lattice_dim()
        2

    You may construct the cone above without giving any rays, but in this case
    you must provide ``lattice`` explicitly::

        sage: origin = Cone([])
        Traceback (most recent call last):
        ...
        ValueError: lattice must be given explicitly if there are no rays!
        sage: origin = Cone([], lattice=ToricLattice(2))
        sage: origin.dim()
        0
        sage: origin.lattice_dim()
        2
        sage: origin.lattice()
        2-d lattice N

    However, the trivial cone in ``n`` dimensions has a predefined
    constructor for you to use::

        sage: origin = cones.trivial(2)
        sage: origin.rays()
        Empty collection
        in 2-d lattice N

    Of course, you can also provide ``lattice`` in other cases::

        sage: L = ToricLattice(3, "L")
        sage: c1 = Cone([(1,0,0),(1,1,1)], lattice=L)
        sage: c1.rays()
        L(1, 0, 0),
        L(1, 1, 1)
        in 3-d lattice L

    Or you can construct cones from rays of a particular lattice::

        sage: ray1 = L(1,0,0)
        sage: ray2 = L(1,1,1)
        sage: c2 = Cone([ray1, ray2])
        sage: c2.rays()
        L(1, 0, 0),
        L(1, 1, 1)
        in 3-d lattice L
        sage: c1 == c2
        True

    When the cone in question is not strictly convex, the standard form for
    the "generating rays" of the linear subspace is "basis vectors and their
    negatives", as in the following example::

        sage: plane = Cone([(1,0), (0,1), (-1,-1)])
        sage: plane.rays()
        N( 0,  1),
        N( 0, -1),
        N( 1,  0),
        N(-1,  0)
        in 2-d lattice N

    The cone can also be specified by a
    :class:`~sage.geometry.polyhedron.base.Polyhedron_base`::

        sage: p = plane.polyhedron()
        sage: Cone(p)
        2-d cone in 2-d lattice N
        sage: Cone(p) == plane
        True

    TESTS::

        sage: N = ToricLattice(2)
        sage: Nsub = N.span([ N(1,2) ])
        sage: Cone(Nsub.basis())
        1-d cone in Sublattice <N(1, 2)>
        sage: Cone([N(0)])
        0-d cone in 2-d lattice N
    """
    # Cone from Polyhedron
    if isinstance(rays, sage.geometry.abc.Polyhedron):
        polyhedron = rays
        if lattice is None:
            lattice = ToricLattice(polyhedron.ambient_dim())
        if polyhedron.n_vertices() > 1:
            raise ValueError("%s is not a cone!" % polyhedron)
        apex = polyhedron.vertices()[0]
        if apex.count(0) != len(apex):
            raise ValueError("the apex of %s is not at the origin!"
                             % polyhedron)
        rays = normalize_rays(polyhedron.rays(), lattice)
        for line in normalize_rays(polyhedron.lines(), lattice):
            rays.append(line)
            rays.append(-line)
            rays[-1].set_immutable()
        return ConvexRationalPolyhedralCone(rays, lattice)
    # Cone from rays
    if check or normalize:
        rays = normalize_rays(rays, lattice)
    if lattice is None:
        if rays:
            lattice = rays[0].parent()
        else:
            raise ValueError(
                "lattice must be given explicitly if there are no rays!")
    if not check or not rays:
        return ConvexRationalPolyhedralCone(rays, lattice)
    # Any set of rays forms a cone, but we want to keep only generators
    if isinstance(lattice, ToricLattice_quotient):
        gs = Generator_System(
                        PPL_point(Linear_Expression(lattice(0).vector(), 0)))
        for r in rays:
            if not r.is_zero():
                gs.insert(PPL_ray(Linear_Expression(r.vector(), 0)))
    else:
        gs = Generator_System( PPL_point(Linear_Expression(lattice(0),0)) )
        for r in rays:
            if not r.is_zero():
                gs.insert( PPL_ray(Linear_Expression(r,0)) )
    cone = C_Polyhedron(gs)
    return _Cone_from_PPL(cone, lattice, rays)


def _Cone_from_PPL(cone, lattice, original_rays=None):
    r"""
    Construct a cone from a :class:`~ppl.polyhedron.Polyhedron`.

    This is a private function and not intended to be exposed to the
    end user. It is used internally by :func:`Cone` and in
    :meth:`ConvexRationalPolyhedralCone.intersection`.

    INPUT:

    - ``cone`` -- a :class:`~ppl.polyhedron.Polyhedron` having the
      origin as its single point

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these.

    - ``original_rays`` -- (default: ``None``) if given, must be a minimal list
      of normalized generating rays of ``cone``. If ``cone`` is strictly convex
      and ``original_rays`` were given, they will be used as internal rays of
      the constructed cone, in the given order.

    OUTPUT: a :class:`ConvexRationalPolyhedralCone`

    TESTS::

        sage: Cone([(1,0), (0,1), (1,1), (0,1)]).rays()   # indirect doctest
        N(0, 1),
        N(1, 0)
        in 2-d lattice N
    """
    rays = []
    lines = []
    for g in cone.minimized_generators():
        if g.is_ray():
            rays.append(g)
        if g.is_line():
            lines.append(g)
    if (original_rays is not None and not lines and
        len(rays) == len(original_rays)):
        return ConvexRationalPolyhedralCone(original_rays, lattice, PPL=cone)
    rays = [ray.coefficients() for ray in rays]
    for line in lines:
        rays.append(line.coefficients())
        rays.append(-vector(ZZ, rays[-1]))
    try:
        for i, ray in enumerate(rays):
            rays[i] = lattice(ray)
            rays[i].set_immutable()
    except TypeError:
        rays = normalize_rays(rays, lattice)
    return ConvexRationalPolyhedralCone(rays, lattice, PPL=cone)


def _ambient_space_point(body, data):
    r"""
    Try to convert ``data`` to a point of the ambient space of ``body``.

    INPUT:

    - ``body`` -- a cone, fan, or lattice polytope with ``lattice()`` method

    - ``data`` -- anything

    OUTPUT:

    An integral, rational, real algebraic, or numeric point of the
    ambient space of ``body`` is returned if ``data`` were
    successfully interpreted in such a way. A :exc:`TypeError` is raised
    otherwise.

    TESTS::

        sage: from sage.geometry.cone import _ambient_space_point
        sage: c = Cone([(1,0), (0,1)])
        sage: _ambient_space_point(c, [1,1])
        N(1, 1)
        sage: _ambient_space_point(c, vector(ZZ,[1,1]))
        N(1, 1)
        sage: _ambient_space_point(c, c.dual_lattice()([1,1]))
        Traceback (most recent call last):
        ...
        TypeError: the point M(1, 1) and
         2-d cone in 2-d lattice N have incompatible lattices
        sage: _ambient_space_point(c, [1,1/3])
        (1, 1/3)
        sage: _ambient_space_point(c, vector(QQ,[1,1/3]))
        (1, 1/3)
        sage: _ambient_space_point(c, [1/2, 1/sqrt(3)])                                 # needs sage.rings.number_field sage.symbolic
        (1/2, 0.5773502691896258?)
        sage: _ambient_space_point(c, vector(AA, [1/2, 1/sqrt(3)]))                     # needs sage.rings.number_field sage.symbolic
        (1/2, 0.5773502691896258?)
        sage: _ambient_space_point(c, [1,1,3])
        Traceback (most recent call last):
        ...
        TypeError: [1, 1, 3] does not represent a valid point
        in the ambient space of 2-d cone in 2-d lattice N
        sage: _ambient_space_point(c, vector(ZZ,[1,1,3]))
        Traceback (most recent call last):
        ...
        TypeError: (1, 1, 3) does not represent a valid point
        in the ambient space of 2-d cone in 2-d lattice N

    Ensure that transcendental elements can, at the very least, be
    represented numerically::

        sage: from sage.geometry.cone import _ambient_space_point
        sage: c = Cone([(1,0), (0,1)])
        sage: _ambient_space_point(c, [1, pi])                                          # needs sage.rings.number_field sage.symbolic
        (1.00000000000000, 3.14159265358979)
        sage: _ambient_space_point(c, vector(SR,[1, pi]))                               # needs sage.rings.number_field sage.symbolic
        (1.00000000000000, 3.14159265358979)
    """
    L = body.lattice()

    def try_base_extend(ring):
        # Factor out the "try this ring..." code that's repeated four
        # times.
        try:
            return L.base_extend(ring)(data)
        except TypeError:
            pass
        except ValueError as ex:
            if str(ex).startswith("Cannot coerce"):
                pass

    # Special treatment for toric lattice elements
    p = try_base_extend(ZZ)
    if p is not None:
        return p
    if isinstance(parent(data), ToricLattice_generic):
        raise TypeError("the point %s and %s have incompatible "
                        "lattices" % (data, body))

    # If we don't have a lattice element, try successively
    # less-desirable ambient spaces until (as a last resort) we
    # attempt a numerical representation.
    rings = [QQ]
    try:
        from sage.rings.qqbar import AA
    except ImportError:
        pass
    else:
        rings.append(AA)
    try:
        from sage.rings.real_mpfr import RR
    except ImportError:
        pass
    else:
        rings.append(RR)

    for ring in rings:
        p = try_base_extend(ring)
        if p is not None:
            return p

    # Raise TypeError with our own message
    raise TypeError("%s does not represent a valid point in the ambient "
                    "space of %s" % (data, body))


def integral_length(v):
    """
    Compute the integral length of a given rational vector.

    INPUT:

    - ``v`` -- any object which can be converted to a list of rationals

    OUTPUT:

    Rational number ``r`` such that ``v = r * u``, where ``u`` is the
    primitive integral vector in the direction of ``v``.

    EXAMPLES::

        sage: from sage.geometry.cone import integral_length
        sage: integral_length([1, 2, 4])
        1
        sage: integral_length([2, 2, 4])
        2
        sage: integral_length([2/3, 2, 4])
        2/3
    """
    data = [QQ(e) for e in list(v)]
    ns = [e.numerator() for e in data]
    ds = [e.denominator() for e in data]
    return gcd(ns) / lcm(ds)


def normalize_rays(rays, lattice):
    r"""
    Normalize a list of rational rays: make them primitive and immutable.

    INPUT:

    - ``rays`` -- list of rays which can be converted to the rational
      extension of ``lattice``;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If ``None``, an attempt will
      be made to determine an appropriate toric lattice automatically.

    OUTPUT:

    - list of immutable primitive vectors of the ``lattice`` in the same
      directions as original ``rays``.

    EXAMPLES::

        sage: from sage.geometry.cone import normalize_rays
        sage: normalize_rays([(0, 1), (0, 2), (3, 2), (5/7, 10/3)], None)
        [N(0, 1), N(0, 1), N(3, 2), N(3, 14)]
        sage: L = ToricLattice(2, "L")
        sage: normalize_rays([(0, 1), (0, 2), (3, 2), (5/7, 10/3)], L.dual())
        [L*(0, 1), L*(0, 1), L*(3, 2), L*(3, 14)]
        sage: ray_in_L = L(0,1)
        sage: normalize_rays([ray_in_L, (0, 2), (3, 2), (5/7, 10/3)], None)
        [L(0, 1), L(0, 1), L(3, 2), L(3, 14)]
        sage: normalize_rays([(0, 1), (0, 2), (3, 2), (5/7, 10/3)], ZZ^2)
        [(0, 1), (0, 1), (3, 2), (3, 14)]
        sage: normalize_rays([(0, 1), (0, 2), (3, 2), (5/7, 10/3)], ZZ^3)
        Traceback (most recent call last):
        ...
        TypeError: cannot convert (0, 1) to
        Vector space of dimension 3 over Rational Field!
        sage: normalize_rays([], ZZ^3)
        []
    """
    if rays is None:
        rays = []
    try:
        rays = list(rays)
    except TypeError:
        raise TypeError(
                    "rays must be given as a list or a compatible structure!"
                    "\nGot: %s" % rays)
    if rays:
        if lattice is None:
            ray_parent = parent(rays[0])
            lattice = (ray_parent if isinstance(ray_parent, ToricLattice_generic)
                                  else ToricLattice(len(rays[0])))
        if lattice.base_ring() is not ZZ:
            raise TypeError("lattice must be a free module over ZZ")
        # Are we dealing with a quotient lattice?
        try:
            if not lattice.is_torsion_free():
                raise ValueError("cannot normalize rays of torsion quotients!")
        except AttributeError:
            pass
        V = None
        try:
            if lattice.is_ambient():
                # Handle the most common case efficiently.
                V = lattice.base_extend(QQ)
                length = integral_length
        except AttributeError:
            pass
        if V is None:
            # Use a more general, but slower way.
            V = lattice.vector_space_span_of_basis(lattice.basis())
            length = lambda ray: integral_length(V.coordinate_vector(ray))
        for n, ray in enumerate(rays):
            try:
                if isinstance(ray, (list, tuple, V.element_class)):
                    ray = V(ray)
                else:
                    ray = V(list(ray))
            except TypeError:
                raise TypeError("cannot convert %s to %s!" % (ray, V))
            if ray.is_zero():
                ray = lattice(0)
            else:
                ray = lattice(ray / length(ray))
            ray.set_immutable()
            rays[n] = ray
    return rays


@richcmp_method
class IntegralRayCollection(SageObject, Hashable, Iterable):
    r"""
    Create a collection of integral rays.

    .. WARNING::

        No correctness check or normalization is performed on the input data.
        This class is designed for internal operations and you probably should
        not use it directly.

    This is a base class for :class:`convex rational polyhedral cones
    <ConvexRationalPolyhedralCone>` and :class:`fans
    <sage.geometry.fan.RationalPolyhedralFan>`.

    Ray collections are immutable, but they cache most of the returned values.

    INPUT:

    - ``rays`` -- list of immutable vectors in ``lattice``;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If ``None``, it will be determined
      as :func:`parent` of the first ray. Of course, this cannot be done if
      there are no rays, so in this case you must give an appropriate
      ``lattice`` directly. Note that ``None`` is *not* the default value -
      you always *must* give this argument explicitly, even if it is ``None``.

    OUTPUT:

    - collection of given integral rays.
    """

    def __init__(self, rays, lattice):
        r"""
        See :class:`IntegralRayCollection` for documentation.

        TESTS::

            sage: from sage.geometry.cone import (
            ....:           IntegralRayCollection)
            sage: v = vector([1,0])
            sage: v.set_immutable()
            sage: c = IntegralRayCollection([v], ZZ^2)
            sage: c = IntegralRayCollection([v], None)
            sage: c.lattice()  # Determined automatically
            Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: c.rays()
            (1, 0)
            in Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: TestSuite(c).run()
        """
        if lattice is None:
            lattice = rays[0].parent()
        self._rays = PointCollection(rays, lattice)
        self._lattice = lattice

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- anything

        OUTPUT: boolean

        There is equality if ``other`` is of the same type as
        ``self``, they have the same ambient lattices, and their
        rays are the same and listed in the same order.

        TESTS::

            sage: c1 = Cone([(1,0), (0,1)])
            sage: c2 = Cone([(0,1), (1,0)])
            sage: c3 = Cone([(0,1), (1,0)])
            sage: c1 > c2
            True
            sage: c2 < c1
            True
            sage: c2 == c3
            True
            sage: c2 is c3
            False
        """
        if type(self) is not type(other):
            return NotImplemented

        # We probably do need to have explicit comparison of lattices here
        # since if one of the collections does not live in a toric lattice,
        # comparison of rays may miss the difference.
        return richcmp((self.lattice(), self.rays()),
                       (other.lattice(), other.rays()), op)

    def __hash__(self):
        r"""
        Return the hash of ``self`` computed from rays.

        OUTPUT: integer

        TESTS::

            sage: c = Cone([(1,0), (0,1)])
            sage: hash(c) == hash(c)
            True
        """
        if "_hash" not in self.__dict__:
            self._hash = hash(self._rays)
        return self._hash

    def __iter__(self):
        r"""
        Return an iterator over rays of ``self``.

        OUTPUT: iterator

        TESTS::

            sage: c = Cone([(1,0), (0,1)])
            sage: for ray in c: print(ray)
            N(1, 0)
            N(0, 1)
        """
        return iter(self._rays)

    def cartesian_product(self, other, lattice=None):
        r"""
        Return the Cartesian product of ``self`` with ``other``.

        INPUT:

        - ``other`` -- an :class:`IntegralRayCollection`;

        - ``lattice`` -- (optional) the ambient lattice for the result. By
          default, the direct sum of the ambient lattices of ``self`` and
          ``other`` is constructed.

        OUTPUT: an :class:`IntegralRayCollection`

        By the Cartesian product of ray collections `(r_0, \dots, r_{n-1})` and
        `(s_0, \dots, s_{m-1})` we understand the ray collection of the form
        `((r_0, 0), \dots, (r_{n-1}, 0), (0, s_0), \dots, (0, s_{m-1}))`, which
        is suitable for Cartesian products of cones and fans. The ray order is
        guaranteed to be as described.

        EXAMPLES::

            sage: c = Cone([(1,)])
            sage: c.cartesian_product(c)    # indirect doctest
            2-d cone in 2-d lattice N+N
            sage: _.rays()
            N+N(1, 0),
            N+N(0, 1)
            in 2-d lattice N+N
        """
        assert isinstance(other, IntegralRayCollection)
        if lattice is None:
            lattice = self.lattice().direct_sum(other.lattice())
        suffix = [0] * other.lattice_dim()
        rays = [lattice(list(r1) + suffix) for r1 in self.rays()]
        prefix = [0] * self.lattice_dim()
        rays.extend(lattice(prefix + list(r2)) for r2 in other.rays())
        for r in rays:
            r.set_immutable()
        return IntegralRayCollection(rays, lattice)

    def __neg__(self):
        """
        Return the collection with opposite rays.

        EXAMPLES::

            sage: c = Cone([(1,1),(0,1)]); c
            2-d cone in 2-d lattice N
            sage: d = -c  # indirect doctest
            sage: d.rays()
            N(-1, -1),
            N( 0, -1)
            in 2-d lattice N
        """
        lattice = self.lattice()
        rays = [-r1 for r1 in self.rays()]
        for r in rays:
            r.set_immutable()
        return IntegralRayCollection(rays, lattice)

    def dim(self):
        r"""
        Return the dimension of the subspace spanned by rays of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: c = Cone([(1,0)])
            sage: c.lattice_dim()
            2
            sage: c.dim()
            1
        """
        if "_dim" not in self.__dict__:
            self._dim = self.rays().matrix().rank()
        return self._dim

    def lattice(self):
        r"""
        Return the ambient lattice of ``self``.

        OUTPUT: lattice

        EXAMPLES::

            sage: c = Cone([(1,0)])
            sage: c.lattice()
            2-d lattice N
            sage: Cone([], ZZ^3).lattice()
            Ambient free module of rank 3
            over the principal ideal domain Integer Ring
        """
        return self._lattice

    def ambient_vector_space(self, base_field=None):
        r"""
        Return the ambient vector space.

        It is the ambient lattice (:meth:`lattice`) tensored with a field.

        INPUT:

        - ``base_field`` -- (default: the rationals) a field

        EXAMPLES::

            sage: c = Cone([(1,0)])
            sage: c.ambient_vector_space()
            Vector space of dimension 2 over Rational Field
            sage: c.ambient_vector_space(AA)                                            # needs sage.rings.number_field
            Vector space of dimension 2 over Algebraic Real Field
        """
        return self.lattice().vector_space(base_field=base_field)

    @cached_method
    def dual_lattice(self):
        r"""
        Return the dual of the ambient lattice of ``self``.

        OUTPUT:

        - lattice. If possible (that is, if :meth:`lattice` has a
          ``dual()`` method), the dual lattice is returned. Otherwise,
          `\ZZ^n` is returned, where `n` is the dimension of :meth:`lattice`.

        EXAMPLES::

            sage: c = Cone([(1,0)])
            sage: c.dual_lattice()
            2-d lattice M
            sage: Cone([], ZZ^3).dual_lattice()
            Ambient free module of rank 3
            over the principal ideal domain Integer Ring

        """
        try:
            return self.lattice().dual()
        except AttributeError:
            return ZZ**self.lattice_dim()

    def lattice_dim(self):
        r"""
        Return the dimension of the ambient lattice of ``self``.

        An alias is :meth:`ambient_dim`.

        OUTPUT: integer

        EXAMPLES::

            sage: c = Cone([(1,0)])
            sage: c.lattice_dim()
            2
            sage: c.dim()
            1
        """
        return self.lattice().dimension()

    ambient_dim = lattice_dim

    def n_rays(self) -> int:
        r"""
        Return the number of rays of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: c = Cone([(1,0), (0,1)])
            sage: c.n_rays()
            2

        TESTS:

        The old method name is kept as an alias::

            sage: c = Cone([(1,0), (0,1)])
            sage: c.nrays()
            2
        """
        return len(self._rays)

    nrays = n_rays

    def plot(self, **options):
        r"""
        Plot ``self``.

        INPUT:

        - any options for toric plots (see :func:`toric_plotter.options
          <sage.geometry.toric_plotter.options>`), none are mandatory.

        OUTPUT: a plot

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant.plot()                                                       # needs sage.plot sage.symbolic
            Graphics object consisting of 9 graphics primitives
        """
        tp = ToricPlotter(options, self.lattice().degree(), self.rays())
        return tp.plot_lattice() + tp.plot_rays() + tp.plot_generators()

    def ray(self, n):
        r"""
        Return the ``n``-th ray of ``self``.

        INPUT:

        - ``n`` -- integer; an index of a ray of ``self``. Enumeration of rays
          starts with zero

        OUTPUT: ray; an element of the lattice of ``self``

        EXAMPLES::

            sage: c = Cone([(1,0), (0,1)])
            sage: c.ray(0)
            N(1, 0)
        """
        return self._rays[n]

    def rays(self, *args):
        r"""
        Return (some of the) rays of ``self``.

        INPUT:

        - ``ray_list`` -- list of integers, the indices of the requested
          rays. If not specified, all rays of ``self`` will be returned

        OUTPUT: a :class:`~sage.geometry.point_collection.PointCollection`
        of primitive integral ray generators

        EXAMPLES::

            sage: c = Cone([(1,0), (0,1), (-1, 0)])
            sage: c.rays()
            N( 0, 1),
            N( 1, 0),
            N(-1, 0)
            in 2-d lattice N
            sage: c.rays([0, 2])
            N( 0, 1),
            N(-1, 0)
            in 2-d lattice N

        You can also give ray indices directly, without packing them into a
        list::

            sage: c.rays(0, 2)
            N( 0, 1),
            N(-1, 0)
            in 2-d lattice N
        """
        return self._rays if not args else self._rays(*args)

    def codim(self):
        r"""
        Return the codimension of ``self``.

        The codimension of a collection of rays (of a cone/fan) is the
        difference between the dimension of the ambient space and the
        dimension of the subspace spanned by those rays (of the cone/fan).

        OUTPUT: nonnegative integer representing the codimension of ``self``

        .. SEEALSO::

            :meth:`dim`, :meth:`lattice_dim`

        EXAMPLES:

        The codimension of the nonnegative orthant is zero, since the
        span of its generators equals the entire ambient space::

            sage: K = cones.nonnegative_orthant(3)
            sage: K.codim()
            0

        However, if we remove a ray so that the entire cone is contained
        within the `x`-`y` plane, then the resulting cone will have
        codimension one, because the `z`-axis is perpendicular to every
        element of the cone::

            sage: K = Cone([(1,0,0), (0,1,0)])
            sage: K.codim()
            1

        If our cone is all of `\mathbb{R}^{2}`, then its codimension is
        zero::

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: K.codim()
            0

        And if the cone is trivial in any space, then its codimension is
        equal to the dimension of the ambient space::

            sage: K = cones.trivial(0)
            sage: K.lattice_dim()
            0
            sage: K.codim()
            0

            sage: K = cones.trivial(1)
            sage: K.lattice_dim()
            1
            sage: K.codim()
            1

            sage: K = cones.trivial(2)
            sage: K.lattice_dim()
            2
            sage: K.codim()
            2

        """
        # same as ConvexSet_base.codim; the main point is the much more detailed
        # docstring.
        return (self.lattice_dim() - self.dim())

    codimension = codim

    def span(self, base_ring=None):
        r"""
        Return the span of ``self``.

        INPUT:

        - ``base_ring`` -- (default: from lattice) the base ring to use
          for the generated module

        OUTPUT: a module spanned by the generators of ``self``

        EXAMPLES:

        The span of a single ray is a one-dimensional sublattice::

            sage: K1 = Cone([(1,)])
            sage: K1.span()
            Sublattice <N(1)>
            sage: K2 = Cone([(1,0)])
            sage: K2.span()
            Sublattice <N(1, 0)>

        The span of the nonnegative orthant is the entire ambient lattice::

            sage: K = cones.nonnegative_orthant(3)
            sage: K.span() == K.lattice()
            True

        By specifying a ``base_ring``, we can obtain a vector space::

            sage: K = Cone([(1,0,0),(0,1,0),(0,0,1)])
            sage: K.span(base_ring=QQ)
            Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]

        TESTS:

        We can take the span of the trivial cone::

            sage: cones.trivial(0).span()
            Sublattice <>

        """
        L = self.lattice()

        if base_ring is None:
            base_ring = L.base_ring()

        return L.span(self, base_ring)

    def _macaulay2_init_(self, macaulay2=None):
        """
        Conversion to Macaulay2.

        EXAMPLES::

            sage: # optional - macaulay2
            sage: m2 = macaulay2
            sage: c = Cone([(3,4), (0,1)])
            sage: m2(c).rays()  # indirect doctest
            | 0 3 |
            | 1 4 |
            sage: m2(c) == c._macaulay2_init_()
            True
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default

        return macaulay2(matrix([tuple(v) for v in self._rays]).transpose())


def classify_cone_2d(ray0, ray1, check=True):
    r"""
    Return `(d,k)` classifying the lattice cone spanned by the two rays.

    INPUT:

    - ``ray0``, ``ray1`` -- two primitive integer vectors; the
      generators of the two rays generating the two-dimensional cone

    - ``check`` -- boolean (default: ``True``); whether to check the
      input rays for consistency

    OUTPUT:

    A pair `(d,k)` of integers classifying the cone up to `GL(2, \ZZ)`
    equivalence. See Proposition 10.1.1 of [CLS2011]_ for the
    definition. We return the unique `(d,k)` with minimal `k`, see
    Proposition 10.1.3 of [CLS2011]_.

    EXAMPLES::

        sage: ray0 = vector([1,0])
        sage: ray1 = vector([2,3])
        sage: from sage.geometry.cone import classify_cone_2d
        sage: classify_cone_2d(ray0, ray1)
        (3, 2)

        sage: ray0 = vector([2,4,5])
        sage: ray1 = vector([5,19,11])
        sage: classify_cone_2d(ray0, ray1)
        (3, 1)

        sage: m = matrix(ZZ, [(19, -14, -115), (-2, 5, 25), (43, -42, -298)])
        sage: m.det()   # check that it is in GL(3,ZZ)
        -1
        sage: classify_cone_2d(m*ray0, m*ray1)
        (3, 1)

    TESTS:

    Check using the connection between the Hilbert basis of the cone
    spanned by the two rays (in arbitrary dimension) and the
    Hirzebruch-Jung continued fraction expansion, see Chapter 10 of
    [CLS2011]_ ::

        sage: from sage.geometry.cone import normalize_rays
        sage: for i in range(10):
        ....:     ray0 = random_vector(ZZ, 3)
        ....:     ray1 = random_vector(ZZ, 3)
        ....:     if ray0.is_zero() or ray1.is_zero(): continue
        ....:     ray0, ray1 = normalize_rays([ray0, ray1], ZZ^3)
        ....:     d, k = classify_cone_2d(ray0, ray1, check=True)
        ....:     assert (d,k) == classify_cone_2d(ray1, ray0)
        ....:     if d == 0: continue
        ....:     frac = (k/d).continued_fraction_list("hj")
        ....:     if len(frac)>100: continue   # avoid expensive computation
        ....:     hilb = Cone([ray0, ray1]).Hilbert_basis()
        ....:     assert len(hilb) == len(frac) + 1
    """
    if check:
        assert ray0.parent() is ray1.parent()
        assert ray0.base_ring() is ZZ
        assert gcd(ray0) == 1
        assert gcd(ray1) == 1
        assert not ray0.is_zero() and not ray1.is_zero()

    m = matrix([ray0, ray1])              # dim(ray) x 2 matrix
    basis = m.saturation().solve_left(m)  # 2-d basis for the span of the cone
    basis = basis.change_ring(ZZ).transpose()
    if basis.nrows() < 2:
        d = 0
        k = basis[0,1]
    else:
        basis.echelonize()                    # columns are the "cone normal form"
        d = basis[1,1]
        k = basis[0,1]

    if check:
        if d == 0:  # degenerate cone
            assert basis[0,0] == 1
            assert k == -1 or k == +1
        else:       # non-degenerate cone
            assert basis[0,0] == 1 and basis[1,0] == 0
            assert d > 0
            assert 0 <= k < d
            assert gcd(d,k) == 1

    # compute unique k, see Proposition 10.1.3 of [CLS2011]
    if d > 0:
        for ktilde in range(k):
            if (k*ktilde) % d == 1:
                k = ktilde
                break
    return (d,k)


# Derived classes MUST allow construction of their objects using ``ambient``
# and ``ambient_ray_indices`` keyword parameters. See ``intersection`` method
# for an example why this is needed.
@richcmp_method
class ConvexRationalPolyhedralCone(IntegralRayCollection, Container, ConvexSet_closed, sage.geometry.abc.ConvexRationalPolyhedralCone):
    r"""
    Create a convex rational polyhedral cone.

    .. WARNING::

        This class does not perform any checks of correctness of input nor
        does it convert input into the standard representation. Use
        :func:`Cone` to construct cones.

    Cones are immutable, but they cache most of the returned values.

    INPUT:

    The input can be either:

    - ``rays`` -- list of immutable primitive vectors in ``lattice``;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If ``None``, it will be determined
      as :func:`parent` of the first ray. Of course, this cannot be done if
      there are no rays, so in this case you must give an appropriate
      ``lattice`` directly.

    or (these parameters must be given as keywords):

    - ``ambient`` -- ambient structure of this cone, a bigger :class:`cone
      <ConvexRationalPolyhedralCone>` or a :class:`fan
      <sage.geometry.fan.RationalPolyhedralFan>`, this cone *must be a face
      of* ``ambient``;

    - ``ambient_ray_indices`` -- increasing list or tuple of integers, indices
      of rays of ``ambient`` generating this cone

    In both cases, the following keyword parameter may be specified in addition:

    - ``PPL`` -- either ``None`` (default) or a
      :class:`~ppl.polyhedron.C_Polyhedron` representing the cone. This
      serves only to cache the polyhedral data if you know it
      already. The constructor does not make a copy so the ``PPL`` object
      should not be modified afterwards.

    OUTPUT: convex rational polyhedral cone

    .. NOTE::

        Every cone has its ambient structure. If it was not specified, it is
        this cone itself.
    """

    def __init__(self, rays=None, lattice=None,
                 ambient=None, ambient_ray_indices=None, PPL=None):
        r"""
        See :class:`ConvexRationalPolyhedralCone` for documentation.

        TESTS::

            sage: from sage.geometry.cone import (
            ....:       ConvexRationalPolyhedralCone)
            sage: v1 = vector([1,0])
            sage: v2 = vector([0,1])
            sage: v1.set_immutable()
            sage: v2.set_immutable()
            sage: ac = ConvexRationalPolyhedralCone([v1, v2], ZZ^2)
            sage: ac = ConvexRationalPolyhedralCone([v1, v2], None)
            sage: ac.lattice()  # Determined automatically
            Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: ac.rays()
            (1, 0),
            (0, 1)
            in Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: ac.ambient() is ac
            True
            sage: TestSuite(ac).run()
            sage: sc = ConvexRationalPolyhedralCone(ambient=ac,
            ....:                       ambient_ray_indices=[1])
            sage: sc.rays()
            (0, 1)
            in Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: sc.ambient() is ac
            True
            sage: TestSuite(sc).run()
        """
        superinit = super().__init__
        if ambient is None:
            superinit(rays, lattice)
            self._ambient = self
            self._ambient_ray_indices = tuple(range(self.n_rays()))
        else:
            self._ambient = ambient
            self._ambient_ray_indices = tuple(ambient_ray_indices)
            superinit(ambient.rays(self._ambient_ray_indices),
                      ambient.lattice())
        if PPL is not None:
            self._PPL_C_Polyhedron = PPL

    def _sage_input_(self, sib, coerced):
        """
        Return Sage command to reconstruct ``self``.

        See :mod:`sage.misc.sage_input` for details.

        EXAMPLES::

            sage: cone = Cone([(1,0), (1,1)])
            sage: sage_input(cone)
            Cone([(1, 0), (1, 1)])
        """
        return sib.name('Cone')([sib(tuple(r)) for r in self.rays()])

    def _PPL_cone(self):
        r"""
        Return the Parma Polyhedra Library (PPL) representation of the cone.

        OUTPUT: a :class:`~ppl.polyhedron.C_Polyhedron` representing the cone

        EXAMPLES::

            sage: c = Cone([(1,0), (1,1), (0,1)])
            sage: c._PPL_cone()
            A 2-dimensional polyhedron in QQ^2 defined as
            the convex hull of 1 point, 2 rays
            sage: c._PPL_cone().minimized_generators()
            Generator_System {point(0/1, 0/1), ray(0, 1), ray(1, 0)}
            sage: c._PPL_cone().minimized_constraints()
            Constraint_System {x1>=0, x0>=0}

        TESTS:

        There are no empty cones, the origin always belongs to them::

            sage: Cone([(0,0)])._PPL_cone()
            A 0-dimensional polyhedron in QQ^2
            defined as the convex hull of 1 point
            sage: cones.trivial(2)._PPL_cone()
            A 0-dimensional polyhedron in QQ^2
            defined as the convex hull of 1 point
        """
        if "_PPL_C_Polyhedron" not in self.__dict__:
            gs = Generator_System(
                            PPL_point(Linear_Expression(self._lattice(0), 0)))
            for r in self.rays():
                gs.insert( PPL_ray(Linear_Expression(r,0)) )
            self._PPL_C_Polyhedron = C_Polyhedron(gs)
        return self._PPL_C_Polyhedron

    def _macaulay2_init_(self, macaulay2=None):
        """
        Conversion to Macaulay2.

        EXAMPLES::

            sage: # optional - macaulay2
            sage: C = Cone([[1,1],[1,2]])
            sage: m2 = macaulay2
            sage: c = m2(C); c.rays()  # indirect doctest
            | 1 1 |
            | 1 2 |
            sage: c.dualCone().rays()  # random
            | 2  -1 |
            | -1 1  |
            sage: c == C._macaulay2_init_()
            True
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default

        return super()._macaulay2_init_(macaulay2).coneFromVData()

    def __contains__(self, point) -> bool:
        r"""
        Check if ``point`` is contained in ``self``.

        See :meth:`_contains` (which is called by this function) for
        documentation.

        TESTS::

            sage: c = Cone([(1,0), (0,1)])
            sage: (1,1) in c
            True
            sage: [1,1] in c
            True
            sage: (-1,0) in c
            False
        """
        return self._contains(point)

    def __getstate__(self):
        r"""
        Return the dictionary that should be pickled.

        OUTPUT: :class:`dict`

        TESTS::

            sage: C = Cone([(1,0)])
            sage: C.face_lattice()                                                      # needs sage.combinat sage.graphs
            Finite lattice containing 2 elements with distinguished linear extension
            sage: C._test_pickling()
            sage: C2 = loads(dumps(C)); C2
            1-d cone in 2-d lattice N
            sage: C2 == C
            True
            sage: C2 is C      # Is this desirable?
            False
        """
        state = copy(self.__dict__)
        state.pop("_PPL_C_Polyhedron", None)  # PPL is not picklable.

        # TODO: do we want to keep the face lattice in the pickle?
        # Currently there is an unpickling loop if do:
        # Unpickling a cone C requires first to unpickle its face lattice.
        # The latter is a Poset which takes C among its arguments. Due
        # to UniqueRepresentation, this triggers a call to hash(C) which
        # itself depends on the attribute C._rays which have not yet
        # been unpickled.  See ``explain_pickle(dumps(C))``.
        state.pop("_face_lattice", None)
        return state

    def _contains(self, point, region='whole cone') -> bool:
        r"""
        Check if ``point`` is contained in ``self``.

        This function is called by :meth:`__contains__` and :meth:`contains`
        to ensure the same call depth for warning messages.

        By default, a point on the boundary of the cone is considered
        part of the cone. If you want to test whether the
        **interior** of the cone contains the point, you need to pass
        the optional argument ``'interior'``.  If you want to test
        whether the **relative interior** of the cone contains the
        point, you need to pass the optional argument
        ``'relative_interior'``.

        .. WARNING::

            The boundary of a closed convex cone is determined by a
            set of inequalities. If your ``point`` has entries in an
            inexact ring, it will sometimes be impossible to say (with
            confidence) if that point lies on the boundary of the cone
            or slightly inside it.

        INPUT:

        - ``point`` -- anything; an attempt will be made to convert it
          into an element compatible with the ambient space of ``self``

        - ``region`` -- string (default: 'whole cone'); can be
          either 'whole cone', 'interior', or 'relative interior'

        OUTPUT:

        ``True`` is returned if ``point`` is contained in the
        specified ``region`` of ``self``. ``False`` is returned
        otherwise, in particular when ``point`` is incompatible with
        the ambient space.

        A :exc:`ValueError` is raised if ``region`` is not one of the
        three allowed values.

        TESTS::

            sage: c = Cone([(1,0), (0,1)])
            sage: c._contains((1,1))
            True

        We can test vectors with irrational components::

            sage: c = Cone([(1,0), (0,1)])
            sage: c._contains((1, sqrt(2)))                                             # needs sage.symbolic
            True
            sage: c._contains(vector(SR, [1, pi]))                                      # needs sage.symbolic
            True

        Ensure that complex vectors are not contained in a real cone::

            sage: c = Cone([(1,0), (0,1)])
            sage: c._contains((1,I))                                                    # needs sage.symbolic
            False
            sage: c._contains(vector(QQbar, [1,I]))                                     # needs sage.rings.number_field sage.symbolic
            False

        And we refuse to coerce elements of another lattice into ours::

            sage: c = Cone([(1,0), (0,1)])
            sage: c._contains(c.dual().ray(0))
            False
        """
        try:
            point = _ambient_space_point(self, point)
        except TypeError as ex:
            if str(ex).endswith("have incompatible lattices!"):
                warn("you have checked if a cone contains a point "
                     "from an incompatible lattice, this is False!",
                     stacklevel=3)
            return False

        if region not in ("whole cone", "relative interior", "interior"):
            raise ValueError("%s is an unknown region of the cone!" % region)
        if region == "interior" and self.dim() < self.lattice_dim():
            return False
        if self.is_trivial():
            # We always have whole cone = relative interior = {0} for
            # the trivial cone, and in addition we have interior = {0}
            # when the ambient space is trivial (the preceding "if"
            # ensures this).
            return point.is_zero()

        need_strict = region.endswith("interior")
        M = self.dual_lattice()
        for c in self._PPL_cone().minimized_constraints():
            pr = M(c.coefficients()) * point
            if c.is_equality():
                if pr != 0:
                    return False
            elif pr < 0 or need_strict and pr == 0:
                return False
        return True

    def interior_contains(self, *args):
        r"""
        Check if a given point is contained in the interior of ``self``.

        For a cone of strictly lower-dimension than the ambient space,
        the interior is always empty. You probably want to use
        :meth:`relative_interior_contains` in this case.

        INPUT:

        - anything. An attempt will be made to convert all arguments into a
          single element of the ambient space of ``self``. If it fails,
          ``False`` will be returned.

        OUTPUT:

        - ``True`` if the given point is contained in the interior of
          ``self``, ``False`` otherwise.

        EXAMPLES::

            sage: c = Cone([(1,0), (0,1)])
            sage: c.contains((1,1))
            True
            sage: c.interior_contains((1,1))
            True
            sage: c.contains((1,0))
            True
            sage: c.interior_contains((1,0))
            False
        """
        point = flatten(args)
        if len(point) == 1:
            point = point[0]
        return self._contains(point, 'interior')

    @cached_method
    def interior(self):
        r"""
        Return the interior of ``self``.

        OUTPUT:

        - either ``self``, an empty polyhedron, or an instance of
          :class:`~sage.geometry.relative_interior.RelativeInterior`.

        EXAMPLES::

            sage: c = Cone([(1,0,0), (0,1,0)]); c
            2-d cone in 3-d lattice N
            sage: c.interior()
            The empty polyhedron in ZZ^3

            sage: origin = cones.trivial(2); origin
            0-d cone in 2-d lattice N
            sage: origin.interior()
            The empty polyhedron in ZZ^2

            sage: K = cones.nonnegative_orthant(2); K
            2-d cone in 2-d lattice N
            sage: K.interior()
            Relative interior of 2-d cone in 2-d lattice N

            sage: K2 = Cone([(1,0),(-1,0),(0,1),(0,-1)]); K2
            2-d cone in 2-d lattice N
            sage: K2.interior() is K2
            True
        """
        if self.is_solid():
            return self.relative_interior()
        return Polyhedron(ambient_dim=self.lattice_dim())

    def relative_interior_contains(self, *args):
        r"""
        Check if a given point is contained in the relative interior of ``self``.

        For a full-dimensional cone the relative interior is simply
        the interior, so this method will do the same check as
        :meth:`interior_contains`. For a strictly lower-dimensional cone, the
        relative interior is the cone without its facets.

        INPUT:

        - anything. An attempt will be made to convert all arguments into a
          single element of the ambient space of ``self``. If it fails,
          ``False`` will be returned.

        OUTPUT:

        - ``True`` if the given point is contained in the relative
          interior of ``self``, ``False`` otherwise.

        EXAMPLES::

            sage: c = Cone([(1,0,0), (0,1,0)])
            sage: c.contains((1,1,0))
            True
            sage: c.relative_interior_contains((1,1,0))
            True
            sage: c.interior_contains((1,1,0))
            False
            sage: c.contains((1,0,0))
            True
            sage: c.relative_interior_contains((1,0,0))
            False
            sage: c.interior_contains((1,0,0))
            False
        """
        point = flatten(args)
        if len(point) == 1:
            point = point[0]
        return self._contains(point, 'relative interior')

    @cached_method
    def relative_interior(self):
        r"""
        Return the relative interior of ``self``.

        OUTPUT:

        - either ``self`` or an instance of
          :class:`~sage.geometry.relative_interior.RelativeInterior`.

        EXAMPLES::

            sage: c = Cone([(1,0,0), (0,1,0)]); c
            2-d cone in 3-d lattice N
            sage: c.relative_interior()
            Relative interior of 2-d cone in 3-d lattice N

            sage: origin = cones.trivial(2); origin
            0-d cone in 2-d lattice N
            sage: origin.relative_interior() is origin
            True

            sage: K1 = Cone([(1,0), (-1,0)]); K1
            1-d cone in 2-d lattice N
            sage: K1.relative_interior() is K1
            True

            sage: K2 = Cone([(1,0),(-1,0),(0,1),(0,-1)]); K2
            2-d cone in 2-d lattice N
            sage: K2.relative_interior() is K2
            True
        """
        if self.is_relatively_open():
            return self
        return RelativeInterior(self)

    def cartesian_product(self, other, lattice=None):
        r"""
        Return the Cartesian product of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a :class:`cone <ConvexRationalPolyhedralCone>`

        - ``lattice`` -- (optional) the ambient lattice for the
          Cartesian product cone; by default, the direct sum of the
          ambient lattices of ``self`` and ``other`` is constructed

        OUTPUT: a :class:`cone <ConvexRationalPolyhedralCone>`

        EXAMPLES::

            sage: c = Cone([(1,)])
            sage: c.cartesian_product(c)
            2-d cone in 2-d lattice N+N
            sage: _.rays()
            N+N(1, 0),
            N+N(0, 1)
            in 2-d lattice N+N
        """
        assert isinstance(other, sage.geometry.abc.ConvexRationalPolyhedralCone)
        rc = super().cartesian_product(other, lattice)
        return ConvexRationalPolyhedralCone(rc.rays(), rc.lattice())

    def __neg__(self):
        """
        Return the cone with opposite rays.

        OUTPUT: a :class:`cone <ConvexRationalPolyhedralCone>`

        EXAMPLES::

            sage: c = Cone([(1,1),(0,1)]); c
            2-d cone in 2-d lattice N
            sage: d = -c; d  # indirect doctest
            2-d cone in 2-d lattice N
            sage: -d == c
            True
            sage: d.rays()
            N(-1, -1),
            N( 0, -1)
            in 2-d lattice N
        """
        rc = super().__neg__()
        return ConvexRationalPolyhedralCone(rc.rays(), rc.lattice())

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- anything

        OUTPUT: boolean

        There is equality if ``self`` and ``other`` are cones of any
        kind in the same lattice with the same rays listed in the
        same order.

        TESTS::

            sage: c1 = Cone([(1,0), (0,1)])
            sage: c2 = Cone([(0,1), (1,0)])
            sage: c3 = Cone([(0,1), (1,0)])
            sage: c1 > c2
            True
            sage: c2 < c1
            True
            sage: c2 == c3
            True
            sage: c2 is c3
            False
        """
        if isinstance(other, sage.geometry.abc.ConvexRationalPolyhedralCone):
            # We don't care about particular type of other in this case
            return richcmp((self.lattice(), self.rays()),
                           (other.lattice(), other.rays()), op)
        return NotImplemented

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant._latex_()
            '\\sigma^{2}'
            sage: quadrant.facets()[0]._latex_()                                        # needs sage.graphs
            '\\sigma^{1} \\subset \\sigma^{2}'
        """
        if self.ambient() is self:
            return r"\sigma^{%d}" % self.dim()
        return r"\sigma^{%d} \subset %s" % (self.dim(),
                                            latex(self.ambient()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant._repr_()
            '2-d cone in 2-d lattice N'
            sage: quadrant
            2-d cone in 2-d lattice N
            sage: quadrant.facets()[0]                                                  # needs sage.graphs
            1-d face of 2-d cone in 2-d lattice N
        """
        result = "%d-d" % self.dim()
        if self.ambient() is self:
            result += " cone in"
            if isinstance(self.lattice(), ToricLattice_generic):
                result += " %s" % self.lattice()
            else:
                result += " %d-d lattice" % self.lattice_dim()
        else:
            result += " face of %s" % self.ambient()
        return result

    def _some_elements_(self):
        r"""
        Generate some points of ``self``.

        EXAMPLES::

            sage: K = cones.nonnegative_orthant(3)
            sage: K.some_elements()  # indirect doctest
            [(0, 0, 0), (1/2, 0, 0), (1/4, 1/2, 0), (1/8, 1/4, 1/2)]
        """
        V = self.ambient_vector_space()
        r_iter = iter(self._rays)
        p = V.zero()
        yield p
        for i in range(5):
            try:
                p = (p + next(r_iter)) / 2
            except StopIteration:
                return
            yield p

    def _sort_faces(self, faces):
        r"""
        Return sorted (if necessary) ``faces`` as a tuple.

        This function ensures that one-dimensional faces are listed in
        agreement with the order of corresponding rays and facets with
        facet normals.

        INPUT:

        - ``faces`` -- iterable of :class:`cones
          <ConvexRationalPolyhedralCone>`

        OUTPUT: :class:`tuple` of :class:`cones <ConvexRationalPolyhedralCone>`

        TESTS::

            sage: octant = Cone(identity_matrix(3).columns())
            sage: # indirect doctest
            sage: for i, face in enumerate(octant.faces(1)):                            # needs sage.graphs
            ....:     if face.ray(0) != octant.ray(i):
            ....:         print("Wrong order!")
        """
        faces = tuple(faces)
        if len(faces) > 1: # Otherwise there is nothing to sort
            if faces[0].n_rays() == 1:
                faces = tuple(sorted(faces,
                                     key=lambda f: f._ambient_ray_indices))
            elif faces[0].dim() == self.dim() - 1 and \
                    self.facet_normals.is_in_cache():
                # If we already have facet normals, sort according to them
                faces = set(faces)
                sorted_faces = [None] * len(faces)
                for i, n in enumerate(self.facet_normals()):
                    for f in faces:
                        if n*f.rays() == 0:
                            sorted_faces[i] = f
                            faces.remove(f)
                            break
                faces = tuple(sorted_faces)
        return faces

    @cached_method
    def adjacent(self):
        r"""
        Return faces adjacent to ``self`` in the ambient face lattice.

        Two *distinct* faces `F_1` and `F_2` of the same face lattice are
        **adjacent** if all of the following conditions hold:

        * `F_1` and `F_2` have the same dimension `d`;

        * `F_1` and `F_2` share a facet of dimension `d-1`;

        * `F_1` and `F_2` are facets of some face of dimension `d+1`, unless
          `d` is the dimension of the ambient structure.

        OUTPUT: :class:`tuple` of :class:`cones <ConvexRationalPolyhedralCone>`

        EXAMPLES::

            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: octant.adjacent()
            ()
            sage: one_face = octant.faces(1)[0]
            sage: len(one_face.adjacent())
            2
            sage: one_face.adjacent()[1]
            1-d face of 3-d cone in 3-d lattice N

        Things are a little bit subtle with fans, as we illustrate below.

        First, we create a fan from two cones in the plane::

            sage: fan = Fan(cones=[(0,1), (1,2)],
            ....:           rays=[(1,0), (0,1), (-1,0)])
            sage: cone = fan.generating_cone(0)
            sage: len(cone.adjacent())                                                  # needs sage.graphs
            1

        The second generating cone is adjacent to this one. Now we create the
        same fan, but embedded into the 3-dimensional space::

            sage: fan = Fan(cones=[(0,1), (1,2)],
            ....:           rays=[(1,0,0), (0,1,0), (-1,0,0)])
            sage: cone = fan.generating_cone(0)
            sage: len(cone.adjacent())                                                  # needs sage.graphs
            1

        The result is as before, since we still have::

            sage: fan.dim()
            2

        Now we add another cone to make the fan 3-dimensional::

            sage: fan = Fan(cones=[(0,1), (1,2), (3,)],
            ....:           rays=[(1,0,0), (0,1,0), (-1,0,0), (0,0,1)])
            sage: cone = fan.generating_cone(0)
            sage: len(cone.adjacent())                                                  # needs sage.graphs
            0

        Since now ``cone`` has smaller dimension than ``fan``, it and its
        adjacent cones must be facets of a bigger one, but since ``cone``
        in this example is generating, it is not contained in any other.
        """
        L = self._ambient._face_lattice_function()
        adjacent = set()
        facets = self.facets()
        superfaces = self.facet_of()
        if superfaces:
            for superface in superfaces:
                for facet in facets:
                    adjacent.update(L.open_interval(facet, superface))
            if adjacent:
                adjacent.remove(L(self))
            return self._sort_faces(adjacent)
        if self.dim() == self._ambient.dim():
            # Special treatment relevant for fans
            for facet in facets:
                adjacent.update(facet.facet_of())
            if adjacent:
                adjacent.remove(self)
            return self._sort_faces(adjacent)
        return ()

    def ambient(self):
        r"""
        Return the ambient structure of ``self``.

        OUTPUT: cone or fan containing ``self`` as a face

        EXAMPLES::

            sage: cone = Cone([(1,2,3), (4,6,5), (9,8,7)])
            sage: cone.ambient()
            3-d cone in 3-d lattice N
            sage: cone.ambient() is cone
            True

            sage: face = cone.faces(1)[0]
            sage: face
            1-d face of 3-d cone in 3-d lattice N
            sage: face.ambient()
            3-d cone in 3-d lattice N
            sage: face.ambient() is cone
            True
        """
        return self._ambient

    def ambient_ray_indices(self):
        r"""
        Return indices of rays of the ambient structure generating ``self``.

        OUTPUT: increasing :class:`tuple` of integers

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant.ambient_ray_indices()
            (0, 1)
            sage: quadrant.facets()[1].ambient_ray_indices()                            # needs sage.graphs
            (1,)
        """
        return self._ambient_ray_indices

    def contains(self, *args):
        r"""
        Check if a given point is contained in ``self``.

        INPUT:

        - anything. An attempt will be made to convert all arguments into a
          single element of the ambient space of ``self``. If it fails,
          ``False`` will be returned.

        OUTPUT:

        - ``True`` if the given point is contained in ``self``, ``False``
          otherwise.

        EXAMPLES::

            sage: c = Cone([(1,0), (0,1)])
            sage: c.contains(c.lattice()(1,0))
            True
            sage: c.contains((1,0))
            True
            sage: c.contains((1,1))
            True
            sage: c.contains(1,1)
            True
            sage: c.contains((-1,0))
            False
            sage: c.contains(c.dual_lattice()(1,0))  # random output (warning)
            False
            sage: c.contains(c.dual_lattice()(1,0))
            False
            sage: c.contains(1)
            False
            sage: c.contains(1/2, sqrt(3))                                              # needs sage.symbolic
            True
            sage: c.contains(-1/2, sqrt(3))                                             # needs sage.symbolic
            False
        """
        point = flatten(args)
        if len(point) == 1:
            point = point[0]
        return self._contains(point)

    def dual(self):
        r"""
        Return the dual cone of ``self``.

        OUTPUT: :class:`cone <ConvexRationalPolyhedralCone>`

        EXAMPLES::

            sage: cone = Cone([(1,0), (-1,3)])
            sage: cone.dual().rays()
            M(0, 1),
            M(3, 1)
            in 2-d lattice M

        Now let's look at a more complicated case::

            sage: cone = Cone([(-2,-1,2), (4,1,0), (-4,-1,-5), (4,1,5)])
            sage: cone.is_strictly_convex()
            False
            sage: cone.dim()
            3
            sage: cone.dual().rays()
            M(7, -18, -2),
            M(1,  -4,  0)
            in 3-d lattice M
            sage: cone.dual().dual() is cone
            True

        We correctly handle the degenerate cases::

            sage: N = ToricLattice(2)
            sage: Cone([], lattice=N).dual().rays()  # empty cone
            M( 1,  0),
            M(-1,  0),
            M( 0,  1),
            M( 0, -1)
            in 2-d lattice M
            sage: Cone([(1,0)], lattice=N).dual().rays()  # ray in 2d
            M(1,  0),
            M(0,  1),
            M(0, -1)
            in 2-d lattice M
            sage: Cone([(1,0),(-1,0)], lattice=N).dual().rays()  # line in 2d
            M(0,  1),
            M(0, -1)
            in 2-d lattice M
            sage: Cone([(1,0),(0,1)], lattice=N).dual().rays()  # strictly convex cone
            M(0, 1),
            M(1, 0)
            in 2-d lattice M
            sage: Cone([(1,0),(-1,0),(0,1)], lattice=N).dual().rays()  # half space
            M(0, 1)
            in 2-d lattice M
            sage: Cone([(1,0),(0,1),(-1,-1)], lattice=N).dual().rays()  # whole space
            Empty collection
            in 2-d lattice M

        """
        if "_dual" not in self.__dict__:
            rays = list(self.facet_normals())
            for ray in self.orthogonal_sublattice().gens():
                rays.append(ray)
                rays.append(-ray)
            self._dual = Cone(rays, lattice=self.dual_lattice(), check=False)
            self._dual._dual = self
        return self._dual

    def embed(self, cone):
        r"""
        Return the cone equivalent to the given one, but sitting in ``self`` as
        a face.

        You may need to use this method before calling methods of ``cone`` that
        depend on the ambient structure, such as
        :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.ambient_ray_indices`
        or
        :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.facet_of`. The
        cone returned by this method will have ``self`` as ambient. If ``cone``
        does not represent a valid cone of ``self``, :exc:`ValueError`
        exception is raised.

        .. NOTE::

            This method is very quick if ``self`` is already the ambient
            structure of ``cone``, so you can use without extra checks and
            performance hit even if ``cone`` is likely to sit in ``self`` but
            in principle may not.

        INPUT:

        - ``cone`` -- a :class:`cone
          <sage.geometry.cone.ConvexRationalPolyhedralCone>`

        OUTPUT:

        - a :class:`cone <sage.geometry.cone.ConvexRationalPolyhedralCone>`,
          equivalent to ``cone`` but sitting inside ``self``.

        EXAMPLES:

        Let's take a 3-d cone on 4 rays::

            sage: c = Cone([(1,0,1), (0,1,1), (-1,0,1), (0,-1,1)])

        Then any ray generates a 1-d face of this cone, but if you construct
        such a face directly, it will not "sit" inside the cone::

            sage: ray = Cone([(0,-1,1)])
            sage: ray
            1-d cone in 3-d lattice N
            sage: ray.ambient_ray_indices()
            (0,)
            sage: ray.adjacent()                                                        # needs sage.graphs
            ()
            sage: ray.ambient()
            1-d cone in 3-d lattice N

        If we want to operate with this ray as a face of the cone, we need to
        embed it first::

            sage: e_ray = c.embed(ray)
            sage: e_ray
            1-d face of 3-d cone in 3-d lattice N
            sage: e_ray.rays()
            N(0, -1, 1)
            in 3-d lattice N
            sage: e_ray is ray
            False
            sage: e_ray.is_equivalent(ray)
            True
            sage: e_ray.ambient_ray_indices()
            (3,)
            sage: e_ray.adjacent()
            (1-d face of 3-d cone in 3-d lattice N,
             1-d face of 3-d cone in 3-d lattice N)
            sage: e_ray.ambient()
            3-d cone in 3-d lattice N

        Not every cone can be embedded into a fixed ambient cone::

            sage: c.embed(Cone([(0,0,1)]))
            Traceback (most recent call last):
            ...
            ValueError: 1-d cone in 3-d lattice N is not a face
            of 3-d cone in 3-d lattice N!
            sage: c.embed(Cone([(1,0,1), (-1,0,1)]))                                    # needs sage.graphs
            Traceback (most recent call last):
            ...
            ValueError: 2-d cone in 3-d lattice N is not a face
            of 3-d cone in 3-d lattice N!
        """
        assert isinstance(cone, sage.geometry.abc.ConvexRationalPolyhedralCone)
        if cone.ambient() is self:
            return cone
        if self.is_strictly_convex():
            rays = self.rays()
            try:
                ray_indices = tuple(sorted(rays.index(ray)
                                           for ray in cone.rays()))
                for face in self.faces(cone.dim()):
                    if face.ambient_ray_indices() == ray_indices:
                        return face
            except ValueError:
                pass
        else:
            # We cannot use the trick with indices since rays are not unique.
            for face in self.faces(cone.dim()):
                if cone.is_equivalent(face):
                    return face
        # If we are here, then either ValueError was raised or we went through
        # all faces and didn't find the matching one.
        raise ValueError("%s is not a face of %s!" % (cone, self))

    def face_lattice(self):
        r"""
        Return the face lattice of ``self``.

        This lattice will have the origin as the bottom (we do not include the
        empty set as a face) and this cone itself as the top.

        OUTPUT:

        - :class:`finite poset <sage.combinat.posets.posets.FinitePoset>` of
          :class:`cones <ConvexRationalPolyhedralCone>`.

        EXAMPLES:

        Let's take a look at the face lattice of the first quadrant::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: L = quadrant.face_lattice()                                           # needs sage.combinat sage.graphs
            sage: L                                                                     # needs sage.combinat sage.graphs
            Finite lattice containing 4 elements with distinguished linear extension

        To see all faces arranged by dimension, you can do this::

            sage: for level in L.level_sets(): print(level)                             # needs sage.combinat sage.graphs
            [0-d face of 2-d cone in 2-d lattice N]
            [1-d face of 2-d cone in 2-d lattice N,
             1-d face of 2-d cone in 2-d lattice N]
            [2-d cone in 2-d lattice N]

        For a particular face you can look at its actual rays... ::

            sage: face = L.level_sets()[1][0]                                           # needs sage.combinat sage.graphs
            sage: face.rays()                                                           # needs sage.combinat sage.graphs
            N(1, 0)
            in 2-d lattice N

        ... or you can see the index of the ray of the original cone that
        corresponds to the above one::

            sage: face.ambient_ray_indices()                                            # needs sage.combinat sage.graphs
            (0,)
            sage: quadrant.ray(0)
            N(1, 0)

        An alternative to extracting faces from the face lattice is to use
        :meth:`faces` method::

            sage: face is quadrant.faces(dim=1)[0]                                      # needs sage.combinat sage.graphs
            True

        The advantage of working with the face lattice directly is that you
        can (relatively easily) get faces that are related to the given one::

            sage: face = L.level_sets()[1][0]                                           # needs sage.combinat sage.graphs
            sage: D = L.hasse_diagram()                                                 # needs sage.combinat sage.graphs
            sage: sorted(D.neighbors(face))                                             # needs sage.combinat sage.graphs
            [0-d face of 2-d cone in 2-d lattice N,
             2-d cone in 2-d lattice N]

        However, you can achieve some of this functionality using
        :meth:`facets`, :meth:`facet_of`, and :meth:`adjacent` methods::

            sage: face = quadrant.faces(1)[0]
            sage: face
            1-d face of 2-d cone in 2-d lattice N
            sage: face.rays()
            N(1, 0)
            in 2-d lattice N
            sage: face.facets()
            (0-d face of 2-d cone in 2-d lattice N,)
            sage: face.facet_of()
            (2-d cone in 2-d lattice N,)
            sage: face.adjacent()
            (1-d face of 2-d cone in 2-d lattice N,)
            sage: face.adjacent()[0].rays()
            N(0, 1)
            in 2-d lattice N

        Note that if ``cone`` is a face of ``supercone``, then the face
        lattice of ``cone`` consists of (appropriate) faces of ``supercone``::

            sage: supercone = Cone([(1,2,3,4), (5,6,7,8),
            ....:                   (1,2,4,8), (1,3,9,7)])
            sage: supercone.face_lattice()
            Finite lattice containing 16 elements with distinguished linear extension
            sage: supercone.face_lattice().top()
            4-d cone in 4-d lattice N
            sage: cone = supercone.facets()[0]
            sage: cone
            3-d face of 4-d cone in 4-d lattice N
            sage: cone.face_lattice()
            Finite poset containing 8 elements with distinguished linear extension
            sage: cone.face_lattice().bottom()
            0-d face of 4-d cone in 4-d lattice N
            sage: cone.face_lattice().top()
            3-d face of 4-d cone in 4-d lattice N
            sage: cone.face_lattice().top() == cone
            True

        TESTS::

            sage: C1 = Cone([(0,1)])
            sage: C2 = Cone([(0,1)])
            sage: C1 == C2
            True
            sage: C1 is C2
            False

        C1 and C2 are equal, but not identical. We currently want them
        to have non identical face lattices, even if the faces
        themselves are equal (see :issue:`10998`)::

            sage: C1.face_lattice() is C2.face_lattice()                                # needs sage.combinat sage.graphs
            False

            sage: C1.facets()[0]                                                        # needs sage.graphs
            0-d face of 1-d cone in 2-d lattice N
            sage: C2.facets()[0]                                                        # needs sage.graphs
            0-d face of 1-d cone in 2-d lattice N

            sage: C1.facets()[0].ambient() is C1                                        # needs sage.graphs
            True
            sage: C2.facets()[0].ambient() is C1                                        # needs sage.graphs
            False
            sage: C2.facets()[0].ambient() is C2                                        # needs sage.graphs
            True
        """
        if "_face_lattice" not in self.__dict__:
            if self._ambient is self:
                # We need to compute face lattice on our own. To accommodate
                # non-strictly convex cones we split rays (or rather their
                # indices) into those in the linear subspace and others, which
                # we refer to as atoms.
                S = self.linear_subspace()
                subspace_rays = []
                atom_to_ray = []
                for i, ray in enumerate(self):
                    # This try...except tests whether ray lies in S;
                    # "ray in S" does not work because ray lies in a
                    # toric lattice and S is a "plain" vector space,
                    # and there is only a conversion (no coercion)
                    # between them as of Github issue #10513.
                    try:
                        S(ray)
                        subspace_rays.append(i)
                    except (TypeError, ValueError):
                        atom_to_ray.append(i)

                def ConeFace(atoms, facets):
                    if facets:
                        rays = sorted([atom_to_ray[a] for a in atoms]
                                      + subspace_rays)
                        face = ConvexRationalPolyhedralCone(
                                    ambient=self, ambient_ray_indices=rays)
                        # It may be nice if this functionality is exposed,
                        # however it makes sense only for cones which are
                        # thought of as faces of a single cone, not of a fan.
                        face._containing_cone_facets = facets
                        return face
                    return self

                # Obtain a modified version of the incidence matrix,
                # with rows corresponding to rays in subspace removed.
                mod_incidence_matrix = self.incidence_matrix()[atom_to_ray]

                atom_to_facets = [row.nonzero_positions()
                                  for row in mod_incidence_matrix.rows()]
                facet_to_atoms = [column.nonzero_positions()
                                  for column in mod_incidence_matrix.columns()]

                self._face_lattice = lattice_from_incidences(
                                    atom_to_facets, facet_to_atoms, ConeFace,
                                    key=id(self))
            else:
                # Get face lattice as a sublattice of the ambient one
                allowed_indices = frozenset(self._ambient_ray_indices)
                from sage.graphs.digraph import DiGraph
                L = DiGraph()
                origin = \
                    self._ambient._face_lattice_function().bottom()
                L.add_vertex(0) # In case it is the only one
                dfaces = [origin]
                faces = [origin]
                face_to_index = {origin:0}
                next_index = 1
                next_d = 1 # Dimension of faces to be considered next.
                while next_d < self.dim():
                    ndfaces = []
                    for face in dfaces:
                        face_index = face_to_index[face]
                        for new_face in face.facet_of():
                            if not allowed_indices.issuperset(
                                            new_face._ambient_ray_indices):
                                continue
                            if new_face in ndfaces:
                                new_face_index = face_to_index[new_face]
                            else:
                                ndfaces.append(new_face)
                                face_to_index[new_face] = next_index
                                new_face_index = next_index
                                next_index += 1
                            L.add_edge(face_index, new_face_index)
                    faces.extend(ndfaces)
                    dfaces = ndfaces
                    next_d += 1
                if self.dim() > 0:
                    # Last level is very easy to build, so we do it separately
                    # even though the above cycle could do it too.
                    faces.append(self)
                    for face in dfaces:
                        L.add_edge(face_to_index[face], next_index)
                D = dict(enumerate(faces))
                L.relabel(D)
                self._face_lattice = FinitePoset(L, faces, key=id(self))
        return self._face_lattice

    # Internally we use this name for a uniform behaviour of cones and fans.
    _face_lattice_function = face_lattice

    def faces(self, dim=None, codim=None):
        r"""
        Return faces of ``self`` of specified (co)dimension.

        INPUT:

        - ``dim`` -- integer; dimension of the requested faces

        - ``codim`` -- integer; codimension of the requested faces

        .. NOTE::

            You can specify at most one parameter. If you don't give any, then
            all faces will be returned.

        OUTPUT:

        - if either ``dim`` or ``codim`` is given, the output will be a
          :class:`tuple` of :class:`cones <ConvexRationalPolyhedralCone>`;

        - if neither ``dim`` nor ``codim`` is given, the output will be the
          :class:`tuple` of tuples as above, giving faces of all existing
          dimensions. If you care about inclusion relations between faces,
          consider using :meth:`face_lattice` or :meth:`adjacent`,
          :meth:`facet_of`, and :meth:`facets`.

        EXAMPLES:

        Let's take a look at the faces of the first quadrant::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant.faces()                                                      # needs sage.graphs
            ((0-d face of 2-d cone in 2-d lattice N,),
             (1-d face of 2-d cone in 2-d lattice N,
              1-d face of 2-d cone in 2-d lattice N),
             (2-d cone in 2-d lattice N,))
            sage: quadrant.faces(dim=1)                                                 # needs sage.graphs
            (1-d face of 2-d cone in 2-d lattice N,
             1-d face of 2-d cone in 2-d lattice N)
            sage: face = quadrant.faces(dim=1)[0]                                       # needs sage.graphs

        Now you can look at the actual rays of this face... ::

            sage: face.rays()                                                           # needs sage.graphs
            N(1, 0)
            in 2-d lattice N

        ... or you can see indices of the rays of the original cone that
        correspond to the above ray::

            sage: face.ambient_ray_indices()                                            # needs sage.graphs
            (0,)
            sage: quadrant.ray(0)
            N(1, 0)

        Note that it is OK to ask for faces of too small or high dimension::

            sage: quadrant.faces(-1)                                                    # needs sage.graphs
            ()
            sage: quadrant.faces(3)                                                     # needs sage.graphs
            ()

        In the case of non-strictly convex cones even faces of small
        nonnegative dimension may be missing::

            sage: halfplane = Cone([(1,0), (0,1), (-1,0)])
            sage: halfplane.faces(0)
            ()
            sage: halfplane.faces()
            ((1-d face of 2-d cone in 2-d lattice N,),
             (2-d cone in 2-d lattice N,))
            sage: plane = Cone([(1,0), (0,1), (-1,-1)])
            sage: plane.faces(1)
            ()
            sage: plane.faces()
            ((2-d cone in 2-d lattice N,),)

        TESTS:

        Now we check that "general" cones whose dimension is smaller than the
        dimension of the ambient space work as expected (see :issue:`9188`)::

            sage: c = Cone([(1,1,1,3),(1,-1,1,3),(-1,-1,1,3)])
            sage: c.faces()                                                             # needs sage.graphs
            ((0-d face of 3-d cone in 4-d lattice N,),
             (1-d face of 3-d cone in 4-d lattice N,
              1-d face of 3-d cone in 4-d lattice N,
              1-d face of 3-d cone in 4-d lattice N),
             (2-d face of 3-d cone in 4-d lattice N,
              2-d face of 3-d cone in 4-d lattice N,
              2-d face of 3-d cone in 4-d lattice N),
             (3-d cone in 4-d lattice N,))

        We also ensure that a call to this function does not break
        :meth:`facets` method (see :issue:`9780`)::

            sage: # needs palp sage.graphs
            sage: cone = toric_varieties.dP8().fan().generating_cone(0); cone
            2-d cone of Rational polyhedral fan in 2-d lattice N
            sage: for f in cone.facets(): print(f.rays())
            N(1, 1)
            in 2-d lattice N
            N(0, 1)
            in 2-d lattice N
            sage: len(cone.faces())
            3
            sage: for f in cone.facets(): print(f.rays())
            N(1, 1)
            in 2-d lattice N
            N(0, 1)
            in 2-d lattice N
        """
        if dim is not None and codim is not None:
            raise ValueError(
                    "dimension and codimension cannot be specified together!")
        dim = self.dim() - codim if codim is not None else dim
        if "_faces" not in self.__dict__:
            self._faces = tuple(map(self._sort_faces,
                                    self.face_lattice().level_sets()))
        if dim is None:
            return self._faces
        lsd = self.linear_subspace().dimension()
        return self._faces[dim - lsd] if lsd <= dim <= self.dim() else ()

    @cached_method
    def facet_normals(self):
        r"""
        Return inward normals to facets of ``self``.

        .. NOTE::

            #. For a not full-dimensional cone facet normals will specify
               hyperplanes whose intersections with the space spanned by
               ``self`` give facets of ``self``.

            #. For a not strictly convex cone facet normals will be orthogonal
               to the linear subspace of ``self``, i.e. they always will be
               elements of the dual cone of ``self``.

            #. The order of normals is random, but consistent with
               :meth:`facets`.

        OUTPUT: a :class:`~sage.geometry.point_collection.PointCollection`

        If the ambient :meth:`~IntegralRayCollection.lattice` of ``self`` is a
        :class:`toric lattice
        <sage.geometry.toric_lattice.ToricLatticeFactory>`, the facet normals
        will be elements of the dual lattice. If it is a general lattice (like
        ``ZZ^n``) that does not have a ``dual()`` method, the facet normals
        will be returned as integral vectors.

        EXAMPLES::

            sage: cone = Cone([(1,0), (-1,3)])
            sage: cone.facet_normals()
            M(0, 1),
            M(3, 1)
            in 2-d lattice M

        Now let's look at a more complicated case::

            sage: cone = Cone([(-2,-1,2), (4,1,0), (-4,-1,-5), (4,1,5)])
            sage: cone.is_strictly_convex()
            False
            sage: cone.dim()
            3
            sage: cone.linear_subspace().dimension()
            1
            sage: lsg = (QQ^3)(cone.linear_subspace().gen(0)); lsg
            (1, 1/4, 5/4)
            sage: cone.facet_normals()
            M(7, -18, -2),
            M(1,  -4,  0)
            in 3-d lattice M
            sage: [lsg*normal for normal in cone.facet_normals()]
            [0, 0]

        A lattice that does not have a ``dual()`` method::

            sage: Cone([(1,1),(0,1)], lattice=ZZ^2).facet_normals()
            (-1, 1),
            ( 1, 0)
            in Ambient free module of rank 2
            over the principal ideal domain Integer Ring

        We correctly handle the degenerate cases::

            sage: N = ToricLattice(2)
            sage: Cone([], lattice=N).facet_normals()  # empty cone
            Empty collection
            in 2-d lattice M
            sage: Cone([(1,0)], lattice=N).facet_normals()  # ray in 2d
            M(1, 0)
            in 2-d lattice M
            sage: Cone([(1,0),(-1,0)], lattice=N).facet_normals()  # line in 2d
            Empty collection
            in 2-d lattice M
            sage: Cone([(1,0),(0,1)], lattice=N).facet_normals()  # strictly convex cone
            M(0, 1),
            M(1, 0)
            in 2-d lattice M
            sage: Cone([(1,0),(-1,0),(0,1)], lattice=N).facet_normals()  # half space
            M(0, 1)
            in 2-d lattice M
            sage: Cone([(1,0),(0,1),(-1,-1)], lattice=N).facet_normals()  # whole space
            Empty collection
            in 2-d lattice M
        """
        cone = self._PPL_cone()
        normals = []
        for c in cone.minimized_constraints():
            assert c.inhomogeneous_term() == 0
            if c.is_inequality():
                normals.append(c.coefficients())
        M = self.dual_lattice()
        normals = tuple(map(M, normals))
        for n in normals:
            n.set_immutable()
        if len(normals) > 1:
            # Sort normals if they are rays
            if self.dim() == 2 and normals[0]*self.ray(0) != 0:
                normals = (normals[1], normals[0])
            else:
                try:    # or if we have combinatorial faces already
                    facets = self._faces[-2]
                    normals = set(normals)
                    sorted_normals = [None] * len(normals)
                    for i, f in enumerate(facets):
                        for n in normals:
                            if n*f.rays() == 0:
                                sorted_normals[i] = n
                                normals.remove(n)
                                break
                    normals = tuple(sorted_normals)
                except AttributeError:
                    pass
        return PointCollection(normals, M)

    @cached_method
    def facet_of(self):
        r"""
        Return *cones* of the ambient face lattice having ``self`` as a facet.

        OUTPUT: :class:`tuple` of :class:`cones <ConvexRationalPolyhedralCone>`

        EXAMPLES::

            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: octant.facet_of()
            ()
            sage: one_face = octant.faces(1)[0]
            sage: len(one_face.facet_of())
            2
            sage: one_face.facet_of()[1]
            2-d face of 3-d cone in 3-d lattice N

        While fan is the top element of its own cone lattice, which is a
        variant of a face lattice, we do not refer to cones as its facets::

            sage: fan = Fan([octant])                                                   # needs sage.graphs
            sage: fan.generating_cone(0).facet_of()                                     # needs sage.graphs
            ()

        Subcones of generating cones work as before::

            sage: one_cone = fan(1)[0]                                                  # needs sage.graphs
            sage: len(one_cone.facet_of())                                              # needs sage.graphs
            2
        """
        L = self._ambient._face_lattice_function()
        H = L.hasse_diagram()
        return self._sort_faces(
            f for f in H.neighbors_out(L(self)) if isinstance(f, sage.geometry.abc.ConvexRationalPolyhedralCone))

    def facets(self):
        r"""
        Return facets (faces of codimension 1) of ``self``.

        OUTPUT: :class:`tuple` of :class:`cones <ConvexRationalPolyhedralCone>`

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant.facets()                                                     # needs sage.graphs
            (1-d face of 2-d cone in 2-d lattice N,
             1-d face of 2-d cone in 2-d lattice N)
        """
        return self.faces(codim=1)

    @cached_method
    def incidence_matrix(self):
        r"""
        Return the incidence matrix.

        .. NOTE::

           The columns correspond to facets/facet normals
           in the order of :meth:`facet_normals`, the rows
           correspond to the rays in the order of
           :meth:`~sage.geometry.cone.IntegralRayCollection.rays`.

        EXAMPLES::

            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: octant.incidence_matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]

            sage: halfspace = Cone([(1,0,0), (0,1,0), (-1,-1,0), (0,0,1)])
            sage: halfspace.incidence_matrix()
            [0]
            [1]
            [1]
            [1]
            [1]

        TESTS::

            sage: halfspace.incidence_matrix().is_immutable()
            True

        Check that the base ring is ``ZZ``, see :issue:`29840`::

            sage: halfspace.incidence_matrix().base_ring()
            Integer Ring
        """
        normals = self.facet_normals()
        incidence_matrix = matrix(ZZ, self.n_rays(),
                                  len(normals), 0)

        for Hindex, normal in enumerate(self.facet_normals()):
            for Vindex, ray in enumerate(self.rays()):
                if normal*ray == 0:
                    incidence_matrix[Vindex, Hindex] = 1

        incidence_matrix.set_immutable()
        return incidence_matrix

    def intersection(self, other):
        r"""
        Compute the intersection of two cones.

        INPUT:

        - ``other`` -- :class:`cone <ConvexRationalPolyhedralCone>`

        OUTPUT: :class:`cone <ConvexRationalPolyhedralCone>`

        This raises :exc:`ValueError` if the ambient space dimensions
        are not compatible.

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (-1, 3)])
            sage: cone2 = Cone([(-1,0), (2, 5)])
            sage: cone1.intersection(cone2).rays()
            N(-1, 3),
            N( 2, 5)
            in 2-d lattice N

        The intersection can also be expressed using the operator ``&``::

            sage: (cone1 & cone2).rays()
            N(-1, 3),
            N( 2, 5)
            in 2-d lattice N

        It is OK to intersect cones living in sublattices of the same ambient
        lattice::

            sage: N = cone1.lattice()
            sage: Ns = N.submodule([(1,1)])
            sage: cone3 = Cone([(1,1)], lattice=Ns)
            sage: I = cone1.intersection(cone3)
            sage: I.rays()
            N(1, 1)
            in Sublattice <N(1, 1)>
            sage: I.lattice()
            Sublattice <N(1, 1)>

        But you cannot intersect cones from incompatible lattices without
        explicit conversion::

            sage: cone1.intersection(cone1.dual())
            Traceback (most recent call last):
            ...
            ValueError: 2-d lattice N and 2-d lattice M
            have different ambient lattices!
            sage: cone1.intersection(Cone(cone1.dual().rays(), N)).rays()
            N(3, 1),
            N(0, 1)
            in 2-d lattice N

        An intersection with a polyhedron returns a polyhedron::

            sage: cone = Cone([(1,0), (-1,0), (0,1)])
            sage: p = polytopes.hypercube(2)
            sage: cone & p
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: sorted(_.vertices_list())
            [[-1, 0], [-1, 1], [1, 0], [1, 1]]
        """
        if not isinstance(other, ConvexRationalPolyhedralCone):
            return self.polyhedron().intersection(other)
        if self._ambient is other._ambient:
            # Cones of the same ambient cone or fan intersect nicely/quickly.
            # Can we maybe even return an element of the cone lattice?..
            # But currently it can be done only for strictly convex cones.
            ambient_ray_indices = tuple(r for r in self._ambient_ray_indices
                                          if r in other._ambient_ray_indices)
            # type(self) allows this code to work nicely for derived classes,
            # although it forces all of them to accept such input
            return type(self)(ambient=self._ambient,
                              ambient_ray_indices=ambient_ray_indices)
        # Generic (slow) intersection, returning a generic cone.
        p = C_Polyhedron(self._PPL_cone())
        p.add_constraints(other._PPL_cone().constraints())
        return _Cone_from_PPL(p, self.lattice().intersection(other.lattice()))

    __and__ = intersection

    def is_equivalent(self, other):
        r"""
        Check if ``self`` is "mathematically" the same as ``other``.

        INPUT:

        - ``other`` -- cone

        OUTPUT:

        - ``True`` if ``self`` and ``other`` define the same cones as sets of
          points in the same lattice, ``False`` otherwise.

        There are three different equivalences between cones `C_1` and `C_2`
        in the same lattice:

        #. They have the same generating rays in the same order.
           This is tested by ``C1 == C2``.
        #. They describe the same sets of points.
           This is tested by ``C1.is_equivalent(C2)``.
        #. They are in the same orbit of `GL(n,\ZZ)` (and, therefore,
           correspond to isomorphic affine toric varieties).
           This is tested by ``C1.is_isomorphic(C2)``.

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (-1, 3)])
            sage: cone2 = Cone([(-1,3), (1, 0)])
            sage: cone1.rays()
            N( 1, 0),
            N(-1, 3)
            in 2-d lattice N
            sage: cone2.rays()
            N(-1, 3),
            N( 1, 0)
            in 2-d lattice N
            sage: cone1 == cone2
            False
            sage: cone1.is_equivalent(cone2)
            True

        """
        if self is other:
            return True
        # TODO: Next check is pointless if cones and fans are made to be unique
        if self.ambient() is other.ambient() and self.is_strictly_convex():
            return self.ambient_ray_indices() == other.ambient_ray_indices()
        if self.lattice() != other.lattice():
            return False
        return self._PPL_cone() == other._PPL_cone()

    def is_face_of(self, cone):
        r"""
        Check if ``self`` forms a face of another ``cone``.

        INPUT:

        - ``cone`` -- cone

        OUTPUT: ``True`` if ``self`` is a face of ``cone``, ``False`` otherwise

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: cone1 = Cone([(1,0)])
            sage: cone2 = Cone([(1,2)])
            sage: quadrant.is_face_of(quadrant)
            True
            sage: cone1.is_face_of(quadrant)
            True
            sage: cone2.is_face_of(quadrant)
            False

        Being a face means more than just saturating a facet
        inequality::

            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: cone = Cone([(2,1,0),(1,2,0)])
            sage: cone.is_face_of(octant)
            False

        """
        if self.lattice() != cone.lattice():
            return False
        if self._ambient is cone._ambient:
            # Cones are always faces of their ambient structure, so
            return self.rays().set().issubset(cone.rays().set())
        if self.is_equivalent(cone):
            return True
        # Obviously False case
        if self.dim() >= cone.dim(): # if == and face, we return True above
            return False

        # It remains to test whether self is a proper face of cone:
        # 1) self must saturate at least one facet inequality
        saturates = Poly_Con_Relation.saturates()
        supporting_hyperplanes = Constraint_System()
        for c in cone._PPL_cone().minimized_constraints():
            rel = self._PPL_cone().relation_with(c)
            if c.is_equality() and not rel.implies(saturates):
                return False
            if c.is_inequality() and rel.implies(saturates):
                c_eq = (Linear_Expression(c.coefficients(),
                                          c.inhomogeneous_term()) == 0)
                supporting_hyperplanes.insert(c_eq)
        if supporting_hyperplanes.empty():
            return False
        # 2) self must be a whole face, and not just a part of one
        cone_face = C_Polyhedron(cone._PPL_cone())
        cone_face.add_constraints(supporting_hyperplanes)
        return cone_face == self._PPL_cone()

    def is_isomorphic(self, other):
        r"""
        Check if ``self`` is in the same `GL(n, \ZZ)`-orbit as ``other``.

        INPUT:

        - ``other`` -- cone

        OUTPUT: ``True`` if ``self`` and ``other`` are in the same
        `GL(n, \ZZ)`-orbit, ``False`` otherwise

        There are three different equivalences between cones `C_1` and `C_2`
        in the same lattice:

        #. They have the same generating rays in the same order.
           This is tested by ``C1 == C2``.
        #. They describe the same sets of points.
           This is tested by ``C1.is_equivalent(C2)``.
        #. They are in the same orbit of `GL(n,\ZZ)` (and, therefore,
           correspond to isomorphic affine toric varieties).
           This is tested by ``C1.is_isomorphic(C2)``.

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (0, 3)])
            sage: m = matrix(ZZ, [(1, -5), (-1, 4)]) # a GL(2,ZZ)-matrix
            sage: cone2 = Cone( m*r for r in cone1.rays() )
            sage: cone1.is_isomorphic(cone2)
            True

            sage: cone1 = Cone([(1,0), (0, 3)])
            sage: cone2 = Cone([(-1,3), (1, 0)])
            sage: cone1.is_isomorphic(cone2)
            False

        TESTS::

            sage: from sage.geometry.cone import classify_cone_2d
            sage: classify_cone_2d(*cone1.rays())
            (1, 0)
            sage: classify_cone_2d(*cone2.rays())
            (3, 2)

        We check that :issue:`18613` is fixed::

            sage: K = cones.trivial(0)
            sage: K.is_isomorphic(K)
            True
            sage: K = cones.trivial(1)
            sage: K.is_isomorphic(K)
            True
            sage: K = cones.trivial(2)
            sage: K.is_isomorphic(K)
            True

        """
        if self.is_strictly_convex() and other.is_strictly_convex():
            from sage.geometry.fan import Fan
            return Fan([self]).is_isomorphic(Fan([other]))
        if self.is_strictly_convex() ^ other.is_strictly_convex():
            return False
        raise NotImplementedError("isomorphism check for not strictly convex "
                                  "cones is not implemented")

    def is_simplicial(self) -> bool:
        r"""
        Check if ``self`` is simplicial.

        A cone is called **simplicial** if primitive vectors along its
        generating rays form a part of a *rational* basis of the ambient
        space.

        OUTPUT: ``True`` if ``self`` is simplicial, ``False`` otherwise

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (0, 3)])
            sage: cone2 = Cone([(1,0), (0, 3), (-1,-1)])
            sage: cone1.is_simplicial()
            True
            sage: cone2.is_simplicial()
            False
        """
        return self.n_rays() == self.dim()

    @cached_method
    def is_smooth(self) -> bool:
        r"""
        Check if ``self`` is smooth.

        A cone is called **smooth** if primitive vectors along its
        generating rays form a part of an *integral* basis of the
        ambient space. Equivalently, they generate the whole lattice
        on the linear subspace spanned by the rays.

        OUTPUT: ``True`` if ``self`` is smooth, ``False`` otherwise

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (0, 1)])
            sage: cone2 = Cone([(1,0), (-1, 3)])
            sage: cone1.is_smooth()
            True
            sage: cone2.is_smooth()
            False

        The following cones are the same up to a `SL(2,\ZZ)`
        coordinate transformation::

            sage: Cone([(1,0,0), (2,1,-1)]).is_smooth()
            True
            sage: Cone([(1,0,0), (2,1,1)]).is_smooth()
            True
            sage: Cone([(1,0,0), (2,1,2)]).is_smooth()
            True
        """
        if not self.is_simplicial():
            return False
        return self.rays().matrix().elementary_divisors() == [1] * self.n_rays()

    def is_empty(self):
        """
        Return whether ``self`` is the empty set.

        Because a cone always contains the origin, this method returns ``False``.

        EXAMPLES::

            sage: trivial_cone = cones.trivial(3)
            sage: trivial_cone.is_empty()
            False
        """
        return False

    def is_trivial(self) -> bool:
        """
        Check if the cone has no rays.

        OUTPUT: ``True`` if the cone has no rays, ``False`` otherwise

        EXAMPLES::

            sage: c0 = cones.trivial(3)
            sage: c0.is_trivial()
            True
            sage: c0.n_rays()
            0
        """
        return self.n_rays() == 0

    is_compact = is_trivial

    @cached_method
    def is_strictly_convex(self) -> bool:
        r"""
        Check if ``self`` is strictly convex.

        A cone is called **strictly convex** if it does not contain any lines.

        OUTPUT: ``True`` if ``self`` is strictly convex, ``False`` otherwise

        .. SEEALSO::

            :meth:`is_pointed`

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (0, 1)])
            sage: cone2 = Cone([(1,0), (-1, 0)])
            sage: cone1.is_strictly_convex()
            True
            sage: cone2.is_strictly_convex()
            False
        """
        return all(not gs.is_line()
                   for gs in self._PPL_cone().minimized_generators())

    def is_pointed(self) -> bool:
        r"""
        Check if ``self`` is pointed.

        A convex cone is said to be *pointed* if (and only if) it
        contains no lines. For a closed convex cone, this is
        equivalent to saying that the origin is the only point
        contained in both the cone and its negation.

        This is an alias for :meth:`is_strictly_convex`. The "pointed"
        terminology is however standard in many settings, and for
        example is consistent with the ``isPointed()`` function in
        Macaulay2.

        OUTPUT:

        ``True`` if :meth:`lineality` is zero, and ``False`` otherwise.

        .. SEEALSO::

            :meth:`is_strictly_convex`,
            :meth:`is_solid`,
            :meth:`is_proper`

        EXAMPLES:

        The Barker-Foran cone must be pointed, because it is self-dual
        under the canonical dual space identification::

            sage: cones.barker_foran().is_pointed()
            True

        The nonnegative orthant and trivial cone are always pointed::

            sage: n = ZZ.random_element(10)
            sage: L = ToricLattice(n)
            sage: cones.nonnegative_orthant(n, lattice=L).is_pointed()
            True
            sage: cones.trivial(lattice=L).is_pointed()
            True

        """
        return self.is_strictly_convex()

    @cached_method
    def linear_subspace(self):
        r"""
        Return the largest linear subspace contained inside of ``self``.

        OUTPUT: subspace of the ambient space of ``self``

        EXAMPLES::

            sage: halfplane = Cone([(1,0), (0,1), (-1,0)])
            sage: halfplane.linear_subspace()
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]

        """
        if self.is_strictly_convex():
            return span([vector(QQ, self.lattice_dim())], QQ)
        return span(self.lines(), QQ)

    @cached_method
    def lines(self):
        r"""
        Return lines generating the linear subspace of ``self``.

        OUTPUT:

        - :class:`tuple` of primitive vectors in the lattice of ``self``
          giving directions of lines that span the linear subspace of
          ``self``.

        EXAMPLES::

            sage: halfplane = Cone([(1,0), (0,1), (-1,0)])
            sage: halfplane.lines()
            N(1, 0)
            in 2-d lattice N
            sage: fullplane = Cone([(1,0), (0,1), (-1,-1)])
            sage: fullplane.lines()
            N(0, 1),
            N(1, 0)
            in 2-d lattice N

        TESTS:

        The line generators are not necessary orthogonal to one
        another::

            sage: K = Cone([(-1, 0, 0),
            ....:           (1, 0, 1),
            ....:           (-1, 0, -1),
            ....:           (1, 1, 0),
            ....:           (-1, -1, 0)])
            sage: K.lines()
            N(1, 0, 1),
            N(1, 1, 0)
            in 3-d lattice N

        """
        lines = []
        for g in self._PPL_cone().minimized_generators():
            if g.is_line():
                lines.append(g.coefficients())
        N = self.lattice()
        lines = tuple(map(N, lines))
        for l in lines:
            l.set_immutable()
        return PointCollection(lines, N)

    def plot(self, **options):
        r"""
        Plot ``self``.

        INPUT:

        - any options for toric plots (see :func:`toric_plotter.options
          <sage.geometry.toric_plotter.options>`), none are mandatory

        OUTPUT: a plot

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant.plot()                                                       # needs sage.plot sage.symbolic
            Graphics object consisting of 9 graphics primitives
        """
        # What to do with 3-d cones in 5-d? Use some projection method?
        deg = self.lattice().degree()
        tp = ToricPlotter(options, deg, self.rays())
        # Modify ray labels to match the ambient cone or fan.
        tp.ray_label = label_list(tp.ray_label, self.n_rays(), deg <= 2,
                                   self.ambient_ray_indices())
        result = tp.plot_lattice() + tp.plot_generators()
        # To deal with non-strictly convex cones we separate rays and labels.
        result += tp.plot_ray_labels()
        tp.ray_label = None
        lsd = self.linear_subspace().dimension()
        if lsd == 1:
            # Plot only rays of the line
            v = self.lines()[0]
            tp.set_rays([v, -v])
        if lsd <= 1:
            result += tp.plot_rays()
        # Modify wall labels to match the ambient cone or fan too.
        walls = self.faces(2)
        try:
            ambient_walls = self.ambient().faces(2)
        except AttributeError:
            ambient_walls = self.ambient().cones(2)
        tp.wall_label = label_list(tp.wall_label, len(walls), deg <= 2,
                            [ambient_walls.index(wall) for wall in walls])
        tp.set_rays(self.ambient().rays())
        result += tp.plot_walls(walls)
        return result

    def polyhedron(self, **kwds):
        r"""
        Return the polyhedron associated to ``self``.

        Mathematically this polyhedron is the same as ``self``.

        OUTPUT: :class:`~sage.geometry.polyhedron.base.Polyhedron_base`

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant.polyhedron()
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull
            of 1 vertex and 2 rays
            sage: line = Cone([(1,0), (-1,0)])
            sage: line.polyhedron()
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull
            of 1 vertex and 1 line

        Here is an example of a trivial cone (see :issue:`10237`)::

            sage: origin = Cone([], lattice=ZZ^2)
            sage: origin.polyhedron()
            A 0-dimensional polyhedron in ZZ^2 defined as the convex hull
            of 1 vertex
        """
        return Polyhedron(rays=self.rays(), vertices=[self.lattice()(0)], **kwds)

    def an_affine_basis(self):
        r"""
        Return points in ``self`` that form a basis for the affine hull.

        EXAMPLES::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: quadrant.an_affine_basis()
            [(0, 0), (1, 0), (0, 1)]
            sage: ray = Cone([(1, 1)])
            sage: ray.an_affine_basis()
            [(0, 0), (1, 1)]
            sage: line = Cone([(1,0), (-1,0)])
            sage: line.an_affine_basis()
            [(1, 0), (0, 0)]
        """
        return [vector(v) for v in self.polyhedron().an_affine_basis()]

    @cached_method
    def strict_quotient(self):
        r"""
        Return the quotient of ``self`` by the linear subspace.

        We define the **strict quotient** of a cone to be the image of this
        cone in the quotient of the ambient space by the linear subspace of
        the cone, i.e. it is the "complementary part" to the linear subspace.

        OUTPUT: cone

        EXAMPLES::

            sage: halfplane = Cone([(1,0), (0,1), (-1,0)])
            sage: ssc = halfplane.strict_quotient()
            sage: ssc
            1-d cone in 1-d lattice N
            sage: ssc.rays()
            N(1)
            in 1-d lattice N
            sage: line = Cone([(1,0), (-1,0)])
            sage: ssc = line.strict_quotient()
            sage: ssc
            0-d cone in 1-d lattice N
            sage: ssc.rays()
            Empty collection
            in 1-d lattice N

        The quotient of the trivial cone is trivial::

            sage: K = cones.trivial(0)
            sage: K.strict_quotient()
            0-d cone in 0-d lattice N
            sage: K = Cone([(0,0,0,0)])
            sage: K.strict_quotient()
            0-d cone in 4-d lattice N

        """
        if self.is_strictly_convex():
            return self
        L = self.lattice()
        Q = L.base_extend(QQ) / self.linear_subspace()
        # Maybe we can improve this one if we create something special
        # for sublattices. But it seems to be the most natural choice
        # for names. If many subcones land in the same lattice -
        # that's just how it goes.
        if isinstance(L, ToricLattice_generic):
            S = ToricLattice(Q.dimension(), L._name, L._dual_name,
                             L._latex_name, L._latex_dual_name)
        else:
            S = ZZ**Q.dimension()
        rays = ( Q(ray) for ray in self if not Q(ray).is_zero() )
        quotient = Cone(rays, S, check=False)
        quotient._is_strictly_convex = True
        return quotient

    @cached_method
    def solid_restriction(self):
        r"""
        Return a solid representation of this cone in terms of a basis
        of its :meth:`sublattice`.

        We define the **solid restriction** of a cone to be a
        representation of that cone in a basis of its own
        sublattice. Since a cone's sublattice is just large enough to
        hold the cone (by definition), the resulting solid restriction
        :meth:`is_solid`. For convenience, the solid restriction lives
        in a new lattice (of the appropriate dimension) and not actually
        in the sublattice object returned by :meth:`sublattice`.

        OUTPUT:

        A solid cone in a new lattice having the same dimension as this
        cone's :meth:`sublattice`.

        EXAMPLES:

        The nonnegative quadrant in the plane is left after we take its
        solid restriction in space::

            sage: K = Cone([(1,0,0), (0,1,0)])
            sage: K.solid_restriction().rays()
            N(0, 1),
            N(1, 0)
            in 2-d lattice N

        The solid restriction of a single ray has the same
        representation regardless of the ambient space::

            sage: K = Cone([(1,0)])
            sage: K.solid_restriction().rays()
            N(1)
            in 1-d lattice N
            sage: K = Cone([(1,1,1)])
            sage: K.solid_restriction().rays()
            N(1)
            in 1-d lattice N

        The solid restriction of the trivial cone lives in a trivial space::

            sage: K = cones.trivial(0)
            sage: K.solid_restriction()
            0-d cone in 0-d lattice N
            sage: K = cones.trivial(4)
            sage: K.solid_restriction()
            0-d cone in 0-d lattice N

        The solid restriction of a solid cone is itself::

            sage: K = Cone([(1,1),(1,2)])
            sage: K.solid_restriction() is K
            True

        """
        if self.is_solid():
            return self
        # Construct a NEW lattice ``S`` (of the appropriate dimension)
        # to use. This works around the fact that it's difficult to
        # work with sublattice objects. There are naming issues here
        # similar to those in the strict_quotient() method.
        L = self.lattice()
        subL = self.sublattice()
        S = ToricLattice(subL.dimension(), L._name,
                         L._dual_name, L._latex_name, L._latex_dual_name)

        # We don't need to check if these rays are zero: they will all
        # have at least one nonzero coordinate; otherwise they would
        # lie outside of the span of our cone. And they don't, because
        # they generate the cone.
        rays = ( S(subL.coordinates(ray)) for ray in self )
        return Cone(rays, lattice=S, check=False)

    def _split_ambient_lattice(self):
        r"""
        Compute a decomposition of the ``N``-lattice into `N_\sigma`
        and its complement isomorphic to `N(\sigma)`.

        You should not call this function directly, but call
        :meth:`sublattice` and :meth:`sublattice_complement` instead.

        EXAMPLES::

            sage: c = Cone([ (1,2) ])
            sage: c._split_ambient_lattice()
            sage: c._sublattice
            Sublattice <N(1, 2)>
            sage: c._sublattice_complement
            Sublattice <N(0, 1)>

        Degenerate cases::

            sage: C2_Z2 = Cone([(1,0),(1,2)])
            sage: C2_Z2._split_ambient_lattice()
            sage: C2_Z2._sublattice
            Sublattice <N(1, 0), N(0, 1)>

        Trivial cone::

            sage: trivial_cone = cones.trivial(3)
            sage: trivial_cone._split_ambient_lattice()
            sage: trivial_cone._sublattice
            Sublattice <>
            sage: trivial_cone._sublattice_complement
            Sublattice <N(1, 0, 0), N(0, 1, 0), N(0, 0, 1)>
        """
        N = self.lattice()
        n = N.dimension()
        basis = self.rays().basis()
        r = len(basis)
        Nsigma = matrix(ZZ, r, n, ( N.coordinates(v) for v in basis ))
        D, U, V = Nsigma.smith_form()  # D = U*N*V <=> N = Uinv*D*Vinv
        basis = (V.inverse() * N.basis_matrix()).rows()
        # spanned lattice N_sigma
        self._sublattice = N.submodule_with_basis(basis[:r])
        # complement to the spanned lattice, isomorphic to N(sigma)
        self._sublattice_complement = N.submodule_with_basis(basis[r:])

    def sublattice(self, *args, **kwds):
        r"""
        The sublattice spanned by the cone.

        Let `\sigma` be the given cone and `N=` ``self.lattice()`` the
        ambient lattice. Then, in the notation of [Ful1993]_, this
        method returns the sublattice

        .. MATH::

            N_\sigma \stackrel{\text{def}}{=} \mathop{span}( N\cap \sigma )

        INPUT:

        - either nothing or something that can be turned into an element of
          this lattice.

        OUTPUT:

        - if no arguments were given, a :class:`toric sublattice
          <sage.geometry.toric_lattice.ToricLattice_sublattice_with_basis>`,
          otherwise the corresponding element of it.

        .. NOTE::

            * The sublattice spanned by the cone is the saturation of
              the sublattice generated by the rays of the cone.

            * If you only need a `\QQ`-basis, you may want to try the
              :meth:`~sage.geometry.point_collection.PointCollection.basis`
              method on the result of :meth:`~IntegralRayCollection.rays`.

            * The returned lattice points are usually not rays of the
              cone. In fact, for a non-smooth cone the rays do not
              generate the sublattice `N_\sigma`, but only a finite
              index sublattice.

        EXAMPLES::

            sage: cone = Cone([(1, 1, 1), (1, -1, 1), (-1, -1, 1), (-1, 1, 1)])
            sage: cone.rays().basis()
            N( 1,  1, 1),
            N( 1, -1, 1),
            N(-1, -1, 1)
            in 3-d lattice N
            sage: cone.rays().basis().matrix().det()
            -4
            sage: cone.sublattice()
            Sublattice <N(1, 1, 1), N(0, -1, 0), N(-1, -1, 0)>
            sage: matrix( cone.sublattice().gens() ).det()
            -1

        Another example::

            sage: c = Cone([(1,2,3), (4,-5,1)])
            sage: c
            2-d cone in 3-d lattice N
            sage: c.rays()
            N(1,  2, 3),
            N(4, -5, 1)
            in 3-d lattice N
            sage: c.sublattice()
            Sublattice <N(4, -5, 1), N(1, 2, 3)>
            sage: c.sublattice(5, -3, 4)
            N(5, -3, 4)
            sage: c.sublattice(1, 0, 0)
            Traceback (most recent call last):
            ...
            TypeError: element [1, 0, 0] is not in free module
        """
        if "_sublattice" not in self.__dict__:
            self._split_ambient_lattice()
        if args or kwds:
            return self._sublattice(*args, **kwds)
        return self._sublattice

    def sublattice_quotient(self, *args, **kwds):
        r"""
        The quotient of the ambient lattice by the sublattice spanned
        by the cone.

        INPUT:

        - either nothing or something that can be turned into an element of
          this lattice.

        OUTPUT:

        - if no arguments were given, a :class:`quotient of a toric lattice
          <sage.geometry.toric_lattice.ToricLattice_quotient>`,
          otherwise the corresponding element of it.

        EXAMPLES::

            sage: C2_Z2 = Cone([(1,0), (1,2)])     # C^2/Z_2
            sage: c1, c2 = C2_Z2.facets()
            sage: c2.sublattice_quotient()
            1-d lattice, quotient of 2-d lattice N by Sublattice <N(1, 2)>
            sage: N = C2_Z2.lattice()
            sage: n = N(1,1)
            sage: n_bar = c2.sublattice_quotient(n); n_bar
            N[1, 1]
            sage: n_bar.lift()
            N(1, 1)
            sage: vector(n_bar)
            (-1)
        """
        if "_sublattice_quotient" not in self.__dict__:
            self._sublattice_quotient = self.lattice() / self.sublattice()
        if args or kwds:
            return self._sublattice_quotient(*args, **kwds)
        return self._sublattice_quotient

    def sublattice_complement(self, *args, **kwds):
        r"""
        A complement of the sublattice spanned by the cone.

        In other words, :meth:`sublattice` and
        :meth:`sublattice_complement` together form a
        `\ZZ`-basis for the ambient :meth:`lattice()
        <sage.geometry.cone.IntegralRayCollection.lattice>`.

        In the notation of [Ful1993]_, let `\sigma` be the given cone
        and `N=` ``self.lattice()`` the ambient lattice. Then this
        method returns

        .. MATH::

            N(\sigma) \stackrel{\text{def}}{=} N / N_\sigma

        lifted (non-canonically) to a sublattice of `N`.

        INPUT:

        - either nothing or something that can be turned into an element of
          this lattice.

        OUTPUT:

        - if no arguments were given, a :class:`toric sublattice
          <sage.geometry.toric_lattice.ToricLattice_sublattice_with_basis>`,
          otherwise the corresponding element of it.

        EXAMPLES::

            sage: C2_Z2 = Cone([(1,0), (1,2)])     # C^2/Z_2
            sage: c1, c2 = C2_Z2.facets()                                               # needs sage.graphs
            sage: c2.sublattice()                                                       # needs sage.graphs
            Sublattice <N(1, 2)>
            sage: c2.sublattice_complement()                                            # needs sage.graphs
            Sublattice <N(0, 1)>

        A more complicated example::

            sage: c = Cone([(1,2,3), (4,-5,1)])
            sage: c.sublattice()
            Sublattice <N(4, -5, 1), N(1, 2, 3)>
            sage: c.sublattice_complement()
            Sublattice <N(2, -3, 0)>
            sage: m = matrix( c.sublattice().gens() + c.sublattice_complement().gens() )
            sage: m
            [ 4 -5  1]
            [ 1  2  3]
            [ 2 -3  0]
            sage: m.det()
            -1
        """
        if "_sublattice_complement" not in self.__dict__:
            self._split_ambient_lattice()
        if args or kwds:
            return self._sublattice_complement(*args, **kwds)
        return self._sublattice_complement

    def orthogonal_sublattice(self, *args, **kwds):
        r"""
        The sublattice (in the dual lattice) orthogonal to the
        sublattice spanned by the cone.

        Let `M=` ``self.dual_lattice()`` be the lattice dual to the
        ambient lattice of the given cone `\sigma`. Then, in the
        notation of [Ful1993]_, this method returns the sublattice

        .. MATH::

            M(\sigma) \stackrel{\text{def}}{=}
            \sigma^\perp \cap M
            \subset M

        INPUT:

        - either nothing or something that can be turned into an element of
          this lattice.

        OUTPUT:

        - if no arguments were given, a :class:`toric sublattice
          <sage.geometry.toric_lattice.ToricLattice_sublattice_with_basis>`,
          otherwise the corresponding element of it.

        EXAMPLES::

            sage: c = Cone([(1,1,1), (1,-1,1), (-1,-1,1), (-1,1,1)])
            sage: c.orthogonal_sublattice()
            Sublattice <>
            sage: c12 = Cone([(1,1,1), (1,-1,1)])
            sage: c12.sublattice()
            Sublattice <N(1, 1, 1), N(0, -1, 0)>
            sage: c12.orthogonal_sublattice()
            Sublattice <M(1, 0, -1)>

        TESTS:

        We check that :issue:`24541` is fixed::

            sage: c = Cone([(1,0)], lattice=ZZ^2)
            sage: c.orthogonal_sublattice()
            Free module of degree 2 and rank 1 over Integer Ring
            User basis matrix:
            [0 1]
            sage: c.dual()
            2-d cone in 2-d lattice
        """
        if "_orthogonal_sublattice" not in self.__dict__:
            try:
                self._orthogonal_sublattice = self.sublattice_quotient().dual()
            except AttributeError:
                N = self.lattice()
                basis = self.rays().basis()
                Nsigma = column_matrix(ZZ, (N.coordinates(v) for v in basis))
                D, U, V = Nsigma.smith_form()  # D = U * Nsigma * V
                M = self.dual_lattice()
                self._orthogonal_sublattice = M.submodule_with_basis(
                    U.rows()[len(basis):])
        if args or kwds:
            return self._orthogonal_sublattice(*args, **kwds)
        return self._orthogonal_sublattice

    def relative_quotient(self, subcone):
        r"""
        The quotient of the spanned lattice by the lattice spanned by
        a subcone.

        In the notation of [Ful1993]_, let `N` be the ambient lattice
        and `N_\sigma` the sublattice spanned by the given cone
        `\sigma`. If `\rho < \sigma` is a subcone, then `N_\rho` =
        ``rho.sublattice()`` is a saturated sublattice of `N_\sigma` =
        ``self.sublattice()``. This method returns the quotient
        lattice. The lifts of the quotient generators are
        `\dim(\sigma)-\dim(\rho)` linearly independent primitive
        lattice points that, together with `N_\rho`, generate
        `N_\sigma`.

        OUTPUT:

        - :class:`toric lattice quotient
          <sage.geometry.toric_lattice.ToricLattice_quotient>`.

        .. NOTE::

            * The quotient `N_\sigma / N_\rho` of spanned sublattices
              has no torsion since the sublattice `N_\rho` is saturated.

            * In the codimension one case, the generator of
              `N_\sigma / N_\rho` is chosen to be in the same direction as the
              image `\sigma / N_\rho`

        EXAMPLES::

            sage: sigma = Cone([(1,1,1,3),(1,-1,1,3),(-1,-1,1,3),(-1,1,1,3)])
            sage: rho   = Cone([(-1, -1, 1, 3), (-1, 1, 1, 3)])
            sage: sigma.sublattice()
            Sublattice <N(1, 1, 1, 3), N(0, -1, 0, 0), N(-1, -1, 0, 0)>
            sage: rho.sublattice()
            Sublattice <N(-1, -1, 1, 3), N(0, 1, 0, 0)>
            sage: sigma.relative_quotient(rho)
            1-d lattice, quotient
            of Sublattice <N(1, 1, 1, 3), N(0, -1, 0, 0), N(-1, -1, 0, 0)>
            by Sublattice <N(1, 0, -1, -3), N(0, 1, 0, 0)>
            sage: sigma.relative_quotient(rho).gens()
            (N[1, 0, 0, 0],)

        More complicated example::

            sage: rho = Cone([(1, 2, 3), (1, -1, 1)])
            sage: sigma = Cone([(1, 2, 3), (1, -1, 1), (-1, 1, 1), (-1, -1, 1)])
            sage: N_sigma = sigma.sublattice()
            sage: N_sigma
            Sublattice <N(1, 2, 3), N(1, -1, 1), N(-1, -1, -2)>
            sage: N_rho = rho.sublattice()
            sage: N_rho
            Sublattice <N(1, -1, 1), N(1, 2, 3)>
            sage: sigma.relative_quotient(rho).gens()
            (N[-1, -1, -2],)
            sage: N = rho.lattice()
            sage: N_sigma == N.span(N_rho.gens() + tuple(q.lift()
            ....:            for q in sigma.relative_quotient(rho).gens()))
            True

        Sign choice in the codimension one case::

            sage: sigma1 = Cone([(1, 2, 3), (1, -1, 1), (-1, 1, 1), (-1, -1, 1)])  # 3d
            sage: sigma2 = Cone([(1, 1, -1), (1, 2, 3), (1, -1, 1), (1, -1, -1)])  # 3d
            sage: rho = sigma1.intersection(sigma2)
            sage: rho.sublattice()
            Sublattice <N(1, -1, 1), N(1, 2, 3)>
            sage: sigma1.relative_quotient(rho)
            1-d lattice, quotient
            of Sublattice <N(1, 2, 3), N(1, -1, 1), N(-1, -1, -2)>
            by Sublattice <N(1, 2, 3), N(0, 3, 2)>
            sage: sigma1.relative_quotient(rho).gens()
            (N[-1, -1, -2],)
            sage: sigma2.relative_quotient(rho).gens()
            (N[0, 2, 1],)
        """
        try:
            cached_values = self._relative_quotient
        except AttributeError:
            self._relative_quotient = {}
            cached_values = self._relative_quotient

        try:
            return cached_values[subcone]
        except KeyError:
            pass

        Ncone = self.sublattice()
        Nsubcone = subcone.sublattice()

        extra_ray = None
        if Ncone.dimension()-Nsubcone.dimension() == 1:
            extra_ray = set(self.rays().set() - subcone.rays().set()).pop()

        Q = Ncone.quotient(Nsubcone, positive_point=extra_ray)
        assert Q.is_torsion_free()
        cached_values[subcone] = Q
        return Q

    def relative_orthogonal_quotient(self, supercone):
        r"""
        The quotient of the dual spanned lattice by the dual of the
        supercone's spanned lattice.

        In the notation of [Ful1993]_, if ``supercone`` = `\rho >
        \sigma` = ``self`` is a cone that contains `\sigma` as a face,
        then `M(\rho)` = ``supercone.orthogonal_sublattice()`` is a
        saturated sublattice of `M(\sigma)` =
        ``self.orthogonal_sublattice()``. This method returns the
        quotient lattice. The lifts of the quotient generators are
        `\dim(\rho)-\dim(\sigma)` linearly independent M-lattice
        lattice points that, together with `M(\rho)`, generate
        `M(\sigma)`.

        OUTPUT:

        - :class:`toric lattice quotient
          <sage.geometry.toric_lattice.ToricLattice_quotient>`.

        If we call the output ``Mrho``, then

            - ``Mrho.cover() == self.orthogonal_sublattice()``, and

            - ``Mrho.relations() == supercone.orthogonal_sublattice()``.

        .. NOTE::

            * `M(\sigma) / M(\rho)` has no torsion since the sublattice
              `M(\rho)` is saturated.

            * In the codimension one case, (a lift of) the generator of
              `M(\sigma) / M(\rho)` is chosen to be positive on `\sigma`.

        EXAMPLES::

            sage: rho = Cone([(1,1,1,3), (1,-1,1,3), (-1,-1,1,3), (-1,1,1,3)])
            sage: rho.orthogonal_sublattice()
            Sublattice <M(0, 0, 3, -1)>
            sage: sigma = rho.facets()[1]
            sage: sigma.orthogonal_sublattice()
            Sublattice <M(0, 1, 1, 0), M(0, 0, 3, -1)>
            sage: sigma.is_face_of(rho)
            True
            sage: Q = sigma.relative_orthogonal_quotient(rho); Q
            1-d lattice, quotient
            of Sublattice <M(0, 1, 1, 0), M(0, 0, 3, -1)>
            by Sublattice <M(0, 0, 3, -1)>
            sage: Q.gens()
            (M[0, 1, 1, 0],)

        Different codimension::

            sage: rho = Cone([[1,-1,1,3],[-1,-1,1,3]])
            sage: sigma = rho.facets()[0]
            sage: sigma.orthogonal_sublattice()
            Sublattice <M(1, 0, 2, -1), M(0, 1, 1, 0), M(0, 0, 3, -1)>
            sage: rho.orthogonal_sublattice()
            Sublattice <M(0, 1, 1, 0), M(0, 0, 3, -1)>
            sage: sigma.relative_orthogonal_quotient(rho).gens()
            (M[-1, 0, -2, 1],)

        Sign choice in the codimension one case::

            sage: sigma1 = Cone([(1, 2, 3), (1, -1, 1), (-1, 1, 1), (-1, -1, 1)])  # 3d
            sage: sigma2 = Cone([(1, 1, -1), (1, 2, 3), (1, -1, 1), (1, -1, -1)])  # 3d
            sage: rho = sigma1.intersection(sigma2)
            sage: rho.relative_orthogonal_quotient(sigma1).gens()
            (M[-5, -2, 3],)
            sage: rho.relative_orthogonal_quotient(sigma2).gens()
            (M[5, 2, -3],)
        """
        try:
            cached_values = self._relative_orthogonal_quotient
        except AttributeError:
            self._relative_orthogonal_quotient = {}
            cached_values = self._relative_orthogonal_quotient

        try:
            return cached_values[supercone]
        except KeyError:
            pass

        Mcone = self.orthogonal_sublattice()
        Msupercone = supercone.orthogonal_sublattice()

        extra_ray = None
        if Mcone.dimension()-Msupercone.dimension() == 1:
            extra_ray = set(supercone.rays().set() - self.rays().set()).pop()

        Q = Mcone.quotient(Msupercone, positive_dual_point=extra_ray)
        assert Q.is_torsion_free()
        cached_values[supercone] = Q
        return Q

    def semigroup_generators(self):
        r"""
        Return generators for the semigroup of lattice points of ``self``.

        OUTPUT:

        - a :class:`~sage.geometry.point_collection.PointCollection`
          of lattice points generating the semigroup of lattice points
          contained in ``self``.

        .. NOTE::

            No attempt is made to return a minimal set of generators, see
            :meth:`Hilbert_basis` for that.

        EXAMPLES:

        The following command ensures that the output ordering in the examples
        below is independent of TOPCOM, you don't have to use it::

            sage: PointConfiguration.set_engine('internal')

        We start with a simple case of a non-smooth 2-dimensional cone::

            sage: Cone([(1,0), (1,2)]).semigroup_generators()
            N(1, 1),
            N(1, 0),
            N(1, 2)
            in 2-d lattice N

        A non-simplicial cone works, too::

            sage: cone = Cone([(3,0,-1), (1,-1,0), (0,1,0), (0,0,1)])
            sage: sorted(cone.semigroup_generators())
            [N(0, 0, 1), N(0, 1, 0), N(1, -1, 0), N(1, 0, 0), N(3, 0, -1)]

        GAP's toric package thinks this is challenging::

            sage: cone = Cone([[1,2,3,4], [0,1,0,7], [3,1,0,2], [0,0,1,0]]).dual()
            sage: len(cone.semigroup_generators())
            2806

        The cone need not be strictly convex::

            sage: halfplane = Cone([(1,0), (2,1), (-1,0)])
            sage: sorted(halfplane.semigroup_generators())
            [N(-1, 0), N(0, 1), N(1, 0)]
            sage: line = Cone([(1,1,1), (-1,-1,-1)])
            sage: sorted(line.semigroup_generators())
            [N(-1, -1, -1), N(1, 1, 1)]
            sage: wedge = Cone([(1,0,0), (1,2,0), (0,0,1), (0,0,-1)])
            sage: sorted(wedge.semigroup_generators())
            [N(0, 0, -1), N(0, 0, 1), N(1, 0, 0), N(1, 1, 0), N(1, 2, 0)]

        Nor does it have to be full-dimensional (see :issue:`11312`)::

            sage: Cone([(1,1,0), (-1,1,0)]).semigroup_generators()
            N( 0, 1, 0),
            N( 1, 1, 0),
            N(-1, 1, 0)
            in 3-d lattice N

        Neither full-dimensional nor simplicial::

            sage: A = matrix([(1, 3, 0), (-1, 0, 1), (1, 1, -2), (15, -2, 0)])
            sage: A.elementary_divisors()
            [1, 1, 1, 0]
            sage: cone3d = Cone([(3,0,-1), (1,-1,0), (0,1,0), (0,0,1)])
            sage: rays = (A*vector(v) for v in cone3d.rays())
            sage: gens = Cone(rays).semigroup_generators(); sorted(gens)
            [N(-2, -1, 0, 17),
             N(0, 1, -2, 0),
             N(1, -1, 1, 15),
             N(3, -4, 5, 45),
             N(3, 0, 1, -2)]
            sage: set(map(tuple,gens)) == set(tuple(A*r) for r in cone3d.semigroup_generators())
            True

        TESTS::

            sage: len(Cone(identity_matrix(10).rows()).semigroup_generators())
            10

            sage: trivial_cone = cones.trivial(3)
            sage: trivial_cone.semigroup_generators()
            Empty collection
            in 3-d lattice N

        ALGORITHM:

        If the cone is not simplicial, it is first triangulated. Each
        simplicial subcone has the integral points of the spaned
        parallelotope as generators. This is the first step of the
        primal Normaliz algorithm, see [Normaliz]_. For each
        simplicial cone (of dimension `d`), the integral points of the
        open parallelotope

        .. MATH::

            par \langle x_1, \dots, x_d \rangle =
            \ZZ^n \cap
            \left\{
            q_1 x_1 + \cdots +q_d x_d
            :~
            0 \leq q_i < 1
            \right\}

        are then computed [BK2001]_.

        Finally, the union of the generators of all simplicial
        subcones is returned.
        """
        # if the cone is not simplicial, triangulate and run
        # recursively
        N = self.lattice()
        if not self.is_simplicial():
            from sage.geometry.triangulation.point_configuration \
                    import PointConfiguration
            origin = self.n_rays() # last one in pc
            pc = PointConfiguration(tuple(self.rays()) + (N(0),), star=origin)
            triangulation = pc.triangulate()
            subcones = ( Cone(( self.ray(i) for i in simplex if i != origin ),
                              lattice=N, check=False)
                         for simplex in triangulation )
            gens = set()
            for cone in subcones:
                gens.update(cone.semigroup_generators())
            return tuple(gens)

        gens = list(parallelotope_points(self.rays(), N)) + list(self.rays())
        gens = ( v for v in gens if gcd(v) == 1 )
        return PointCollection(gens, N)

    @cached_method
    def Hilbert_basis(self):
        r"""
        Return the Hilbert basis of the cone.

        Given a strictly convex cone `C\subset \RR^d`, the Hilbert
        basis of `C` is the set of all irreducible elements in the
        semigroup `C\cap \ZZ^d`. It is the unique minimal generating
        set over `\ZZ` for the integral points `C\cap \ZZ^d`.

        If the cone `C` is not strictly convex, this method finds the
        (unique) minimal set of lattice points that need to be added
        to the defining rays of the cone to generate the whole
        semigroup `C\cap \ZZ^d`. But because the rays of the cone are
        not unique nor necessarily minimal in this case, neither is
        the returned generating set (consisting of the rays plus
        additional generators).

        See also :meth:`semigroup_generators` if you are not
        interested in a minimal set of generators.

        OUTPUT:

        - a
          :class:`~sage.geometry.point_collection.PointCollection`. The
          rays of ``self`` are the first ``self.n_rays()`` entries.

        EXAMPLES:

        The following command ensures that the output ordering in the examples
        below is independent of TOPCOM, you don't have to use it::

            sage: PointConfiguration.set_engine('internal')

        We start with a simple case of a non-smooth 2-dimensional cone::

            sage: Cone([(1,0), (1,2)]).Hilbert_basis()
            N(1, 0),
            N(1, 2),
            N(1, 1)
            in 2-d lattice N

        A more complicated example from GAP/toric::

            sage: Cone([[1,0], [3,4]]).dual().Hilbert_basis()
            M(0,  1),
            M(4, -3),
            M(1,  0),
            M(2, -1),
            M(3, -2)
            in 2-d lattice M

        Examples verified independently using [Normaliz]_::

            sage: cones.rearrangement(3,4).Hilbert_basis()
            N(-2,  1,  1,  1),
            N( 1,  1,  1, -2),
            N( 1,  1, -2,  1),
            N( 1, -2,  1,  1),
            N( 1,  0,  0,  0),
            N(-1,  1,  1,  0),
            N(-1,  1,  0,  1),
            N( 0, -1,  1,  1),
            N( 1,  0, -1,  1),
            N( 1,  0,  1, -1),
            N(-1,  0,  1,  1),
            N( 0,  1, -1,  1),
            N( 1, -1,  0,  1),
            N( 0,  1,  1, -1),
            N( 1, -1,  1,  0),
            N( 0,  1,  0,  0),
            N( 1,  1,  0, -1),
            N( 1,  1, -1,  0),
            N( 0,  0,  1,  0),
            N( 0,  0,  0,  1)
            in 4-d lattice N

            sage: # long time
            sage: K = Cone([(1,0,1), (-1,0,1), (0,1,1), (0,-1,1)])
            sage: K.Hilbert_basis()
            N( 1,  0, 1),
            N(-1,  0, 1),
            N( 0,  1, 1),
            N( 0, -1, 1),
            N( 0,  0, 1)
            in 3-d lattice N

            sage: # long time
            sage: K = Cone([(1,0,1,0), (-1,0,1,0), (0,1,1,0), (0,-1,1,0)])
            sage: K.Hilbert_basis()
            N( 1,  0, 1, 0),
            N(-1,  0, 1, 0),
            N( 0,  1, 1, 0),
            N( 0, -1, 1, 0),
            N( 0,  0, 1, 0)
            in 4-d lattice N

        Not a strictly convex cone::

            sage: wedge = Cone([(1,0,0), (1,2,0), (0,0,1), (0,0,-1)])
            sage: sorted(wedge.semigroup_generators())
            [N(0, 0, -1), N(0, 0, 1), N(1, 0, 0), N(1, 1, 0), N(1, 2, 0)]
            sage: wedge.Hilbert_basis()
            N(1, 2,  0),
            N(1, 0,  0),
            N(0, 0,  1),
            N(0, 0, -1),
            N(1, 1,  0)
            in 3-d lattice N

        Not full-dimensional cones are ok, too (see :issue:`11312`)::

            sage: Cone([(1,1,0), (-1,1,0)]).Hilbert_basis()
            N( 1, 1, 0),
            N(-1, 1, 0),
            N( 0, 1, 0)
            in 3-d lattice N

        ALGORITHM:

        The primal Normaliz algorithm, see [Normaliz]_.
        """

        # Used in lieu of our linear_subspace() for containment
        # testing when this cone is strictly convex. None of the
        # points we test can be zero, so checking the empty tuple is
        # equivalent to checking a trivial subspace, and the empty
        # tuple is faster.
        L = ()

        if not self.is_strictly_convex():
            # Our linear_subspace(), but as a cone, so that
            # containment testing using "in" works properly.
            L = Cone((c*r for c in (1, -1) for r in self.lines()),
                     self.lattice(),
                     check=False)

        irreducible = list(self.rays())  # these are irreducible for sure
        irr_modified = False  # have we appended to "irreducible"?
        gens = [x for x in self.semigroup_generators()
                if x not in irreducible]

        from itertools import chain
        while gens:
            x = gens.pop()
            if all((y in L) or (x-y not in self)
                   for y in chain(irreducible, gens)):
                irreducible.append(x)
                irr_modified = True

        if irr_modified:
            return PointCollection(irreducible, self.lattice())

        # Avoid the PointCollection overhead if nothing was
        # added to the irreducible list beyond self.rays().
        return self.rays()

    def Hilbert_coefficients(self, point, solver=None, verbose=0,
                             *, integrality_tolerance=1e-3):
        r"""
        Return the expansion coefficients of ``point`` with respect to
        :meth:`Hilbert_basis`.

        INPUT:

        - ``point`` -- a :meth:`~IntegralRayCollection.lattice` point
          in the cone, or something that can be converted to a
          point. For example, a list or tuple of integers.

        - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear Programming
          (MILP) solver to be used. If set to ``None``, the default one is used. For
          more information on MILP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: 0); sets the level of verbosity
          of the LP solver. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- parameter for use with MILP solvers over an
          inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`

        OUTPUT:

        A `\ZZ`-vector of length ``len(self.Hilbert_basis())`` with nonnegative
        components.

        .. NOTE::

            Since the Hilbert basis elements are not necessarily linearly
            independent, the expansion coefficients are not unique. However,
            this method will always return the same expansion coefficients when
            invoked with the same argument.

        EXAMPLES::

            sage: cone = Cone([(1,0), (0,1)])
            sage: cone.rays()
            N(1, 0),
            N(0, 1)
            in 2-d lattice N
            sage: cone.Hilbert_coefficients([3,2])
            (3, 2)

        A more complicated example::

            sage: N = ToricLattice(2)
            sage: cone = Cone([N(1,0), N(1,2)])
            sage: cone.Hilbert_basis()
            N(1, 0),
            N(1, 2),
            N(1, 1)
            in 2-d lattice N
            sage: cone.Hilbert_coefficients(N(1,1))
            (0, 0, 1)

        The cone need not be strictly convex::

            sage: N = ToricLattice(3)
            sage: cone = Cone([N(1,0,0), N(1,2,0), N(0,0,1), N(0,0,-1)])
            sage: cone.Hilbert_basis()
            N(1, 2,  0),
            N(1, 0,  0),
            N(0, 0,  1),
            N(0, 0, -1),
            N(1, 1,  0)
            in 3-d lattice N
            sage: cone.Hilbert_coefficients(N(1,1,3))
            (0, 0, 3, 0, 1)
        """
        point = self.lattice()(point)
        if point not in self:
            raise ValueError('The given point is not in the cone!')
        basis = self.Hilbert_basis()

        from sage.numerical.mip import MixedIntegerLinearProgram
        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        p.set_objective(None)
        x = p.new_variable(integer=True, nonnegative=True)
        for i in range(self.lattice_dim()):
            p.add_constraint(p.sum(b[i]*x[j] for j,b in enumerate(basis)) == point[i])
        p.solve(log=verbose)

        return vector(ZZ, p.get_values(x, convert=ZZ, tolerance=integrality_tolerance))

    def is_solid(self):
        r"""
        Check if this cone is solid.

        A cone is said to be solid if it has nonempty interior. That
        is, if its extreme rays span the entire ambient space.

        An alias is :meth:`is_full_dimensional`.

        OUTPUT: ``True`` if this cone is solid, and ``False`` otherwise

        .. SEEALSO::

            :meth:`is_proper`

        EXAMPLES:

        The nonnegative orthant is always solid::

            sage: quadrant = cones.nonnegative_orthant(2)
            sage: quadrant.is_solid()
            True
            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: octant.is_solid()
            True

        However, if we embed the two-dimensional nonnegative quadrant
        into three-dimensional space, then the resulting cone no longer
        has interior, so it is not solid::

            sage: quadrant = Cone([(1,0,0), (0,1,0)])
            sage: quadrant.is_solid()
            False

        """
        return (self.dim() == self.lattice_dim())

    is_full_dimensional = is_solid

    def is_proper(self):
        r"""
        Check if this cone is proper.

        A cone is said to be proper if it is closed, convex, solid,
        and contains no lines. This cone is assumed to be closed and
        convex; therefore it is proper if it is solid and contains no
        lines.

        OUTPUT: ``True`` if this cone is proper, and ``False`` otherwise

        .. SEEALSO::

            :meth:`is_strictly_convex`,
            :meth:`is_pointed`,
            :meth:`is_solid`

        EXAMPLES:

        The nonnegative orthant is always proper::

            sage: quadrant = cones.nonnegative_orthant(2)
            sage: quadrant.is_proper()
            True
            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: octant.is_proper()
            True

        However, if we embed the two-dimensional nonnegative quadrant
        into three-dimensional space, then the resulting cone no longer
        has interior, so it is not solid, and thus not proper::

            sage: quadrant = Cone([(1,0,0), (0,1,0)])
            sage: quadrant.is_proper()
            False

        Likewise, a half-space contains at least one line, so it is not
        proper::

            sage: halfspace = Cone([(1,0), (0,1), (-1,0)])
            sage: halfspace.is_proper()
            False
        """
        return (self.is_strictly_convex() and self.is_solid())

    def is_full_space(self):
        r"""
        Check if this cone is equal to its ambient vector space.

        An alias is :meth:`is_universe`.

        OUTPUT:

        ``True`` if this cone equals its entire ambient vector
        space and ``False`` otherwise.

        EXAMPLES:

        A single ray in two dimensions is not equal to the entire
        space::

            sage: K = Cone([(1,0)])
            sage: K.is_full_space()
            False

        Neither is the nonnegative orthant::

            sage: K = cones.nonnegative_orthant(2)
            sage: K.is_full_space()
            False

        The right half-space contains a vector subspace, but it is
        still not equal to the entire space::

            sage: K = Cone([(1,0), (-1,0), (0,1)])
            sage: K.is_full_space()
            False

        However, if we allow conic combinations of both axes, then
        the resulting cone is the entire two-dimensional space::

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
        """
        return self.linear_subspace() == self.lattice().vector_space()

    is_universe = is_full_space

    def lineality(self):
        r"""
        Return the lineality of this cone.

        The lineality of a cone is the dimension of the largest linear
        subspace contained in that cone.

        OUTPUT:

        A nonnegative integer; the dimension of the largest subspace
        contained within this cone.

        REFERENCES:

        - [Roc1970]_

        EXAMPLES:

        The lineality of the nonnegative orthant is zero, since it clearly
        contains no lines::

            sage: K = cones.nonnegative_orthant(3)
            sage: K.lineality()
            0

        However, if we add another ray so that the entire `x`-axis belongs
        to the cone, then the resulting cone will have lineality one::

            sage: K = Cone([(1,0,0), (-1,0,0), (0,1,0), (0,0,1)])
            sage: K.lineality()
            1

        If our cone is all of `\mathbb{R}^{2}`, then its lineality is equal
        to the dimension of the ambient space (i.e. two)::

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: K.lineality()
            2
            sage: K.lattice_dim()
            2

        Per the definition, the lineality of the trivial cone in a trivial
        space is zero::

            sage: K = cones.trivial(0)
            sage: K.lineality()
            0

        """
        return self.linear_subspace().dimension()

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is relatively open.

        OUTPUT: boolean

        EXAMPLES::

            sage: K = cones.nonnegative_orthant(3)
            sage: K.is_relatively_open()
            False

            sage: K1 = Cone([(1,0), (-1,0)]); K1
            1-d cone in 2-d lattice N
            sage: K1.is_relatively_open()
            True
        """
        return self.lineality() == self.dim()

    @cached_method
    def discrete_complementarity_set(self):
        r"""
        Compute a discrete complementarity set of this cone.

        A discrete complementarity set of a cone is the set of all
        orthogonal pairs `(x,s)` where `x` is in some fixed generating
        set of the cone, and `s` is in some fixed generating set of its
        dual. The generators chosen for this cone and its dual are
        simply their :meth:`~IntegralRayCollection.rays`.

        OUTPUT:

        A tuple of pairs `(x,s)` such that,

          * `x` and `s` are nonzero.
          * `s(x)` is zero.
          * `x` is one of this cone's :meth:`~IntegralRayCollection.rays`.
          * `s` is one of the :meth:`~IntegralRayCollection.rays` of this
            cone's :meth:`dual`.

        REFERENCES:

        - [Or2017]_

        EXAMPLES:

        Pairs of standard basis elements form a discrete complementarity
        set for the nonnegative orthant::

            sage: K = cones.nonnegative_orthant(2)
            sage: K.discrete_complementarity_set()
            ((N(1, 0), M(0, 1)), (N(0, 1), M(1, 0)))

        If a cone consists of a single ray, then the second components
        of a discrete complementarity set for that cone should generate
        the orthogonal complement of the ray::

            sage: K = Cone([(1,0)])
            sage: K.discrete_complementarity_set()
            ((N(1, 0), M(0, 1)), (N(1, 0), M(0, -1)))
            sage: K = Cone([(1,0,0)])
            sage: K.discrete_complementarity_set()
            ((N(1, 0, 0), M(0, 1, 0)),
             (N(1, 0, 0), M(0, -1, 0)),
             (N(1, 0, 0), M(0, 0, 1)),
             (N(1, 0, 0), M(0, 0, -1)))

        When a cone is the entire space, its dual is the trivial cone,
        so the only discrete complementarity set for it is empty::

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: K.discrete_complementarity_set()
            ()

        Likewise for trivial cones, whose duals are the entire space::

            sage: cones.trivial(0).discrete_complementarity_set()
            ()

        """
        # Return an immutable tuple instead of a mutable list because
        # the result will be cached.
        return tuple( (x,s) for x in self
                            for s in self.dual()
                            if s*x == 0 )

    def lyapunov_like_basis(self):
        r"""
        Compute a basis of Lyapunov-like transformations on this cone.

        A linear transformation `L` is said to be Lyapunov-like on this
        cone if `L(x)` and `s` are orthogonal for every pair `(x,s)` in
        its :meth:`discrete_complementarity_set`. The set of all such
        transformations forms a vector space, namely the Lie algebra of
        the automorphism group of this cone.

        OUTPUT:

        A list of matrices forming a basis for the space of all
        Lyapunov-like transformations on this cone.

        .. SEEALSO::

            :meth:`cross_positive_operators_gens`,
            :meth:`positive_operators_gens`,
            :meth:`Z_operators_gens`

        REFERENCES:

        - [Or2017]_
        - [RNPA2011]_

        EXAMPLES:

        Every transformation is Lyapunov-like on the trivial cone::

            sage: K = cones.trivial(2)
            sage: M = MatrixSpace(K.lattice().base_field(), K.lattice_dim())
            sage: list(M.basis()) == K.lyapunov_like_basis()
            True

        And by duality, every transformation is Lyapunov-like on the
        ambient space::

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: M = MatrixSpace(K.lattice().base_field(), K.lattice_dim())
            sage: list(M.basis()) == K.lyapunov_like_basis()
            True

        However, in a trivial space, there are no non-trivial linear maps,
        so there can be no Lyapunov-like basis::

            sage: K = cones.trivial(0)
            sage: K.lyapunov_like_basis()
            []

        The Lyapunov-like transformations on the nonnegative orthant are
        diagonal matrices::

            sage: K = cones.nonnegative_orthant(1)
            sage: K.lyapunov_like_basis()
            [[1]]

            sage: K = cones.nonnegative_orthant(2)
            sage: K.lyapunov_like_basis()
            [
            [1 0]  [0 0]
            [0 0], [0 1]
            ]

            sage: K = cones.nonnegative_orthant(3)
            sage: K.lyapunov_like_basis()
            [
            [1 0 0]  [0 0 0]  [0 0 0]
            [0 0 0]  [0 1 0]  [0 0 0]
            [0 0 0], [0 0 0], [0 0 1]
            ]

        Only the identity matrix is Lyapunov-like on the pyramids
        defined by the one- and infinity-norms [RNPA2011]_::

            sage: l31 = Cone([(1,0,1), (0,-1,1), (-1,0,1), (0,1,1)])
            sage: l31.lyapunov_like_basis()
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]

            sage: l3infty = Cone([(0,1,1), (1,0,1), (0,-1,1), (-1,0,1)])
            sage: l3infty.lyapunov_like_basis()
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]

        """
        # Matrices are not vectors in Sage, so we have to convert them
        # to vectors explicitly before we can find a basis. We need these
        # two values to construct the appropriate "long vector" space.
        F = self.lattice().base_field()
        n = self.lattice_dim()

        # These tensor products contain a basis for the orthogonal
        # complement of the Lyapunov-like transformations on this cone.
        tensor_products = ( s.tensor_product(x)
                            for (x,s) in self.discrete_complementarity_set() )

        # Convert those tensor products to long vectors.
        W = VectorSpace(F, n**2)
        perp_vectors = ( W(tp.list()) for tp in tensor_products )

        # Now find the Lyapunov-like transformations (as long vectors).
        LL_vectors = W.span(perp_vectors).complement()

        # And finally convert the long vectors back to matrices.
        M = MatrixSpace(F, n, n)
        return [ M(v.list()) for v in LL_vectors.basis() ]

    def lyapunov_rank(self):
        r"""
        Compute the Lyapunov rank of this cone.

        The Lyapunov rank of a cone is the dimension of the space of its
        Lyapunov-like transformations --- that is, the length of a
        :meth:`lyapunov_like_basis`. Equivalently, the Lyapunov rank is
        the dimension of the Lie algebra of the automorphism group of
        the cone.

        OUTPUT: nonnegative integer representing the Lyapunov rank of this cone

        If the ambient space is trivial, then the Lyapunov rank will be
        zero. On the other hand, if the dimension of the ambient vector
        space is `n > 0`, then the resulting Lyapunov rank will be
        between `1` and `n^2` inclusive. If this cone :meth:`is_proper`,
        then that upper bound reduces from `n^2` to `n`. A Lyapunov rank
        of `n-1` is not possible (by Lemma 6 [Or2017]_) in either case.

        ALGORITHM:

        Algorithm 3 [Or2017]_ is used. Every closed convex cone is
        isomorphic to a Cartesian product of a proper cone, a subspace,
        and a trivial cone. The Lyapunov ranks of the subspace and
        trivial cone are easy to compute. Essentially, we "peel off"
        those easy parts of the cone and compute their Lyapunov ranks
        separately. We then compute the rank of the proper cone by
        counting a :meth:`lyapunov_like_basis` for it. Summing the
        individual ranks gives the Lyapunov rank of the original cone.

        REFERENCES:

        - [GT2014]_
        - [Or2017]_
        - [RNPA2011]_

        EXAMPLES:

        The Lyapunov rank of the nonnegative orthant is the same as the
        dimension of the ambient space [RNPA2011]_::

            sage: positives = cones.nonnegative_orthant(1)
            sage: positives.lyapunov_rank()
            1
            sage: quadrant = cones.nonnegative_orthant(2)
            sage: quadrant.lyapunov_rank()
            2
            sage: octant = cones.nonnegative_orthant(3)
            sage: octant.lyapunov_rank()
            3

        A vector space of dimension `n` has Lyapunov rank `n^{2}`
        [Or2017]_::

            sage: Q5 = VectorSpace(QQ, 5)
            sage: gs = Q5.basis() + [-r for r in Q5.basis()]
            sage: K = Cone(gs)
            sage: K.lyapunov_rank()
            25

        A pyramid in three dimensions has Lyapunov rank one [RNPA2011]_::

            sage: l31 = Cone([(1,0,1), (0,-1,1), (-1,0,1), (0,1,1)])
            sage: l31.lyapunov_rank()
            1
            sage: l3infty = Cone([(0,1,1), (1,0,1), (0,-1,1), (-1,0,1)])
            sage: l3infty.lyapunov_rank()
            1

        A ray in `n` dimensions has Lyapunov rank `n^{2} - n + 1`
        [Or2017]_::

            sage: K = Cone([(1,0,0,0,0)])
            sage: K.lyapunov_rank()
            21
            sage: K.lattice_dim()**2 - K.lattice_dim() + 1
            21

        A subspace of dimension `m` in an `n`-dimensional ambient space
        has Lyapunov rank `n^{2} - m(n - m)` [Or2017]_::

            sage: e1 = vector(QQ, [1,0,0,0,0])
            sage: e2 = vector(QQ, [0,1,0,0,0])
            sage: z = (0,0,0,0,0)
            sage: K = Cone([e1, -e1, e2, -e2, z, z, z])
            sage: K.lyapunov_rank()
            19
            sage: K.lattice_dim()**2 - K.dim()*K.codim()
            19

        Lyapunov rank is additive on a product of proper cones [RNPA2011]_::

            sage: l31 = Cone([(1,0,1), (0,-1,1), (-1,0,1), (0,1,1)])
            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: K = l31.cartesian_product(octant)
            sage: K.lyapunov_rank()
            4
            sage: l31.lyapunov_rank() + octant.lyapunov_rank()
            4

        Two linearly-isomorphic cones have the same Lyapunov rank
        [RNPA2011]_. A cone linearly-isomorphic to the nonnegative octant
        will have Lyapunov rank ``3``::

            sage: K = Cone([(1,2,3), (-1,1,0), (1,0,6)])
            sage: K.lyapunov_rank()
            3

        """
        # The solid_restriction() and strict_quotient() methods
        # already check if the cone is solid or strictly convex, so we
        # can't save any additional time here by seeing if those
        # methods would be no-ops.
        #
        # The call to solid_restriction() restricts K to its own span,
        # resulting in the cone K_S from the paper. The call to
        # strict_quotient() then restricts K_S to the span of its dual.
        K_SP = self.solid_restriction().strict_quotient()

        # K_SP is proper, so we have to compute its Lyapunov rank the
        # hard way -- by counting a Lyapunov-like basis for it.
        m = self.dim()
        n = self.lattice_dim()
        l = self.lineality()

        # cf. Theorem 2
        return len(K_SP.lyapunov_like_basis()) + l*m + (n - m)*n

    def random_element(self, ring=ZZ):
        r"""
        Return a random element of this cone.

        All elements of a convex cone can be represented as a
        nonnegative linear combination of its generators. A random
        element is thus constructed by assigning random nonnegative
        weights to the generators of this cone.  By default, these
        weights are integral and the resulting random element will live
        in the same lattice as the cone.

        The random nonnegative weights are chosen from ``ring`` which
        defaults to ``ZZ``. When ``ring`` is not ``ZZ``, the random
        element returned will be a vector. Only the rings ``ZZ`` and
        ``QQ`` are currently supported.

        INPUT:

          - ``ring`` -- (default: ``ZZ``) the ring from which the random
            generator weights are chosen; either ``ZZ`` or ``QQ``

        OUTPUT:

        Either a lattice element or vector contained in both this cone
        and its ambient vector space. If ``ring`` is ``ZZ``, a lattice
        element is returned; otherwise a vector is returned. If ``ring``
        is neither ``ZZ`` nor ``QQ``, then a :exc:`NotImplementedError` is
        raised.

        EXAMPLES:

        The trivial element ``()`` is always returned in a trivial space::

            sage: K = cones.trivial(0)
            sage: K.random_element()
            N()
            sage: K.random_element(ring=QQ)
            ()

        A random element of the trivial cone in a nontrivial space is zero::

            sage: K = cones.trivial(3)
            sage: K.random_element()
            N(0, 0, 0)
            sage: K.random_element(ring=QQ)
            (0, 0, 0)

        A random element of the nonnegative orthant should have all
        components nonnegative::

            sage: K = cones.nonnegative_orthant(3)
            sage: all(x >= 0 for x in K.random_element())
            True
            sage: all(x >= 0 for x in K.random_element(ring=QQ))
            True

        If ``ring`` is not ``ZZ`` or ``QQ``, an error is raised::

            sage: K = Cone([(1,0), (0,1)])
            sage: K.random_element(ring=RR)
            Traceback (most recent call last):
            ...
            NotImplementedError: ring must be either ZZ or QQ.

        """
        if ring not in [ZZ, QQ]:
            # This cone theoretically lives in a real vector space,
            # but in Sage, we work over the rationals to avoid
            # numerical issues. Thus ``ring`` must consist of
            # rationals so that the ambient vector space will contain
            # the resulting random element.
            raise NotImplementedError('ring must be either ZZ or QQ.')

        # The lattice or vector space in which the return value will live.
        L = self.lattice()
        if ring is not ZZ:
            L = L.vector_space()

        # Scale each generator by a random nonnegative factor.
        terms = ( ring.random_element().abs()*L(g) for g in self )

        # Make sure we return a lattice element or vector. Without the
        # explicit conversion, we return ``0`` when we have no rays.
        return L(sum(terms))

    def _positive_operators_dual(self, K2):
        r"""
        Return the dual cone of the positive operators from
        ``self`` to ``K2`` under the matrix <-> long-vector
        isometry.

        In :meth:`positive_operators_gens` we compute the cone of
        positive operators from ``self`` to ``K2`` by taking the
        dual of a dual. This method computes the latter dual. Matrices
        cannot generate cones in Sage, so long vectors are used as
        proxies in this step. The matrix <-> long-vector map is an
        isometry, so everything works out the same in the end.

        The intermediate dual cone of long vectors turns out to be
        useful in at least one doctest, so it has been factored out
        for efficiency.

        REFERENCES:

        - [Or2018b]_

        EXAMPLES:

        On the nonnegative orthant, the positive operators are the
        (self-dual) cone of nonnegative matrices, which in long-vector
        form is just a (bigger) nonnegative orthant::

            sage: K = cones.nonnegative_orthant(3)
            sage: J = K._positive_operators_dual(K)
            sage: J.is_equivalent(cones.nonnegative_orthant(9))
            True

        """
        # Matrices are not vectors in Sage, so we have to convert them
        # to vectors explicitly before we can find a basis. We need these
        # two values to construct the appropriate "long vector" space.
        F = self.lattice().base_field()
        n = self.lattice_dim()
        m = K2.lattice_dim()

        tensor_products = ( s.tensor_product(x) for x in self
                                                for s in K2.dual() )

        # Convert those tensor products to long vectors.
        W = VectorSpace(F, n*m)
        vectors = ( W(tp.list()) for tp in tensor_products )

        check = True
        if self.is_proper() and K2.is_proper():
            # All of the generators involved are extreme vectors and
            # therefore minimal. If this cone is neither solid nor
            # strictly convex, then the tensor product of ``s`` and ``x``
            # is the same as that of ``-s`` and ``-x``. However, as a
            # /set/, ``tensor_products`` may still be minimal.
            check = False

        # Create the dual cone of the positive operators, expressed as
        # long vectors.
        return Cone(vectors, ToricLattice(W.dimension()), check=check)

    def positive_operators_gens(self, K2=None):
        r"""
        Compute minimal generators of the positive operators on this cone.

        A linear operator on a cone is positive if the image of
        the cone under the operator is a subset of the cone. This
        concept can be extended to two cones: the image of the
        first cone under a positive operator is a subset of the
        second cone, which may live in a different space.

        The positive operators (on one or two fixed cones) themselves
        form a closed convex cone. This method computes and returns
        the generators of that cone as a list of matrices.

        INPUT:

        - ``K2`` -- (default: ``self``) the codomain cone; the image of
          this cone under the returned generators is a subset of ``K2``

        OUTPUT:

        A list of `m`-by-`n` matrices where `m` is the ambient dimension
        of ``K2`` and `n` is the ambient dimension of this cone. Each
        matrix `P` in the list has the property that `P(x)` is an
        element of ``K2`` whenever `x` is an element of this cone.

        The returned matrices generate the cone of positive operators
        from this cone to ``K2``; that is,

        - Any nonnegative linear combination of the returned matrices
          sends elements of this cone to ``K2``.

        - Every positive operator on this cone (with respect to ``K2``)
          is some nonnegative linear combination of the returned matrices.

        ALGORITHM:

        Computing positive operators directly is difficult, but
        computing their dual is straightforward using the generators of
        Berman and Gaiha. We construct the dual of the positive
        operators, and then return the dual of that, which is guaranteed
        to be the desired positive operators because everything is
        closed, convex, and polyhedral.

        .. SEEALSO::

            :meth:`cross_positive_operators_gens`,
            :meth:`lyapunov_like_basis`,
            :meth:`Z_operators_gens`

        REFERENCES:

        - [BG1972]_
        - [BP1994]_
        - [Or2018b]_

        EXAMPLES:

        Positive operators on the nonnegative orthant are nonnegative
        matrices::

            sage: K = Cone([(1,)])
            sage: K.positive_operators_gens()
            [[1]]

            sage: K = Cone([(1,0), (0,1)])
            sage: K.positive_operators_gens()
            [
            [1 0]  [0 1]  [0 0]  [0 0]
            [0 0], [0 0], [1 0], [0 1]
            ]

        The trivial cone in a trivial space has no positive operators::

            sage: K = cones.trivial(0)
            sage: K.positive_operators_gens()
            []

        Every operator is positive on the trivial cone::

            sage: K = cones.trivial(1)
            sage: K.positive_operators_gens()
            [[1], [-1]]

            sage: K = cones.trivial(2)
            sage: K.is_trivial()
            True
            sage: K.positive_operators_gens()
            [
            [1 0]  [-1  0]  [0 1]  [ 0 -1]  [0 0]  [ 0  0]  [0 0]  [ 0  0]
            [0 0], [ 0  0], [0 0], [ 0  0], [1 0], [-1  0], [0 1], [ 0 -1]
            ]

        Every operator is positive on the ambient vector space::

            sage: K = Cone([(1,), (-1,)])
            sage: K.is_full_space()
            True
            sage: K.positive_operators_gens()
            [[1], [-1]]

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: K.positive_operators_gens()
            [
            [1 0]  [-1  0]  [0 1]  [ 0 -1]  [0 0]  [ 0  0]  [0 0]  [ 0  0]
            [0 0], [ 0  0], [0 0], [ 0  0], [1 0], [-1  0], [0 1], [ 0 -1]
            ]

        A non-obvious application is to find the positive operators on the
        right half-plane [Or2018b]_::

            sage: K = Cone([(1,0), (0,1), (0,-1)])
            sage: K.positive_operators_gens()
            [
            [1 0]  [0 0]  [ 0  0]  [0 0]  [ 0  0]
            [0 0], [1 0], [-1  0], [0 1], [ 0 -1]
            ]

        TESTS:

        The trivial cone, full space, and half-plane all give rise to the
        expected dimensions [Or2018b]_::

            sage: n = ZZ.random_element(5)
            sage: K = cones.trivial(n)
            sage: L = ToricLattice(n^2)
            sage: pi_gens = K.positive_operators_gens()
            sage: pi_cone = Cone((g.list() for g in pi_gens),
            ....:                lattice=L,
            ....:                check=False)
            sage: pi_cone.dim() == n^2
            True

            sage: K = K.dual()
            sage: K.is_full_space()
            True
            sage: pi_gens = K.positive_operators_gens()
            sage: pi_cone = Cone((g.list() for g in pi_gens),
            ....:                lattice=L,
            ....:                check=False)
            sage: pi_cone.dim() == n^2
            True

            sage: K = Cone([(1,0), (0,1), (0,-1)])
            sage: pi_gens = K.positive_operators_gens()
            sage: pi_cone = Cone((g.list() for g in pi_gens),
            ....:                check=False)
            sage: pi_cone.dim() == 3
            True

        The trivial cone, full space, and half-plane all give rise to
        the expected linealities [Or2018b]_::

            sage: n = ZZ.random_element(5)
            sage: K = cones.trivial(n)
            sage: L = ToricLattice(n^2)
            sage: pi_gens = K.positive_operators_gens()
            sage: pi_cone = Cone((g.list() for g in pi_gens),
            ....:                lattice=L,
            ....:                check=False)
            sage: pi_cone.lineality() == n^2
            True

            sage: K = K.dual()
            sage: K.is_full_space()
            True
            sage: pi_gens = K.positive_operators_gens()
            sage: pi_cone = Cone((g.list() for g in pi_gens),
            ....:                lattice=L,
            ....:                check=False)
            sage: pi_cone.lineality() == n^2
            True

            sage: K = Cone([(1,0), (0,1), (0,-1)])
            sage: pi_gens = K.positive_operators_gens()
            sage: pi_cone = Cone((g.list() for g in pi_gens), check=False)
            sage: pi_cone.lineality() == 2
            True

        """
        if K2 is None:
            K2 = self

        # Create cone of positive operators, expressed as long
        # vectors.
        pi_cone = self._positive_operators_dual(K2).dual()

        # And finally convert its rays back to matrix representations.
        F = self.lattice().base_field()
        n = self.lattice_dim()
        m = K2.lattice_dim()
        M = MatrixSpace(F, m, n)
        return [ M(v.list()) for v in pi_cone ]

    def _cross_positive_operators_dual(self):
        r"""
        Return the dual cone of the cross-positive operators on
        ``self`` under the matrix <-> long-vector isometry.

        In :meth:`cross_positive_operators_gens` we compute the cone
        of cross positive operators on ``self`` by taking the dual of
        its dual. This method computes that dual. Matrices cannot
        generate cones in Sage, so long vectors are used as proxies in
        this step. The matrix <-> long-vector map is an isometry, so
        everything works out the same in the end.

        The intermediate dual cone of long vectors turns out to be
        useful in one of the tests for :meth:`is_reducible`, so it has
        been factored out for convenience.

        REFERENCES:

        - [Or2018b]_

        EXAMPLES:

        For the nonnegative orthant generated by the standard basis
        vectors `e_{i}`, we should return the cone generated by the
        long-vector counterparts of all `e_{i}e_{j}^{T}` such that
        `e_{i}` and `e_{j}` are orthogonal. There are `n^{2} - n` such
        pairs::

            sage: K = cones.nonnegative_orthant(3)
            sage: K._cross_positive_operators_dual().rays()
            N(0, 0, 0, 1, 0, 0, 0, 0, 0),
            N(0, 0, 0, 0, 0, 0, 1, 0, 0),
            N(0, 1, 0, 0, 0, 0, 0, 0, 0),
            N(0, 0, 0, 0, 0, 0, 0, 1, 0),
            N(0, 0, 1, 0, 0, 0, 0, 0, 0),
            N(0, 0, 0, 0, 0, 1, 0, 0, 0)
            in 9-d lattice N
        """
        # These tensor products contain generators for the dual cone of
        # the cross-positive operators.
        tensor_products = ( s.tensor_product(x)
                            for (x,s) in self.discrete_complementarity_set() )

        # Turn our matrices into long vectors (lists of coordinates).
        vectors = ( m.list() for m in tensor_products )

        check = True
        if self.is_proper():
            # All of the generators involved are extreme vectors and
            # therefore minimal. If this cone is neither solid nor
            # strictly convex, then the tensor product of ``s`` and
            # ``x`` is the same as that of ``-s`` and ``-x``. However,
            # as a /set/, ``tensor_products`` may still be minimal.
            check = False

        # Create the dual cone of the cross-positive operators,
        # expressed as long vectors. Specify a lattice in case
        # I'm the trivial cone.
        return Cone(vectors,
                    lattice=ToricLattice(self.lattice_dim()**2),
                    check=check)

    def cross_positive_operators_gens(self):
        r"""
        Compute minimal generators of the cross-positive operators on this
        cone.

        Any positive operator `P` on this cone will have `s(P(x)) \ge 0`
        whenever `x` is an element of this cone and `s` is an element of
        its dual. By contrast, the cross-positive operators need only
        satisfy that property on the :meth:`discrete_complementarity_set`;
        that is, when `x` and `s` are "cross" (orthogonal).

        The cross-positive operators (on some fixed cone) themselves
        form a closed convex cone. This method computes and returns
        the generators of that cone as a list of matrices.

        Cross-positive operators are also called exponentially-positive,
        since they become positive operators when exponentiated. Other
        equivalent names are resolvent-positive, essentially-positive,
        and quasimonotone.

        OUTPUT:

        A list of `n`-by-`n` matrices where `n` is the ambient dimension
        of this cone. Each matrix `L` in the list has the property that
        `s(L(x)) \ge 0` whenever `(x,s)` is an element of this cone's
        :meth:`discrete_complementarity_set`.

        The returned matrices generate the cone of cross-positive operators
        on this cone; that is,

        - Any nonnegative linear combination of the returned matrices
          is cross-positive on this cone.

        - Every cross-positive operator on this cone is some nonnegative
          linear combination of the returned matrices.

        .. SEEALSO::

           :meth:`lyapunov_like_basis`,
           :meth:`positive_operators_gens`,
           :meth:`Z_operators_gens`

        REFERENCES:

        - [SV1970]_
        - [Or2018b]_

        EXAMPLES:

        Cross-positive operators on the nonnegative orthant are
        negations of Z-matrices; that is, matrices whose off-diagonal
        elements are nonnegative::

            sage: K = cones.nonnegative_orthant(2)
            sage: K.cross_positive_operators_gens()
            [
            [0 1]  [0 0]  [1 0]  [-1  0]  [0 0]  [ 0  0]
            [0 0], [1 0], [0 0], [ 0  0], [0 1], [ 0 -1]
            ]
            sage: K = Cone([(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)])
            sage: all(c[i][j] >= 0 for c in K.cross_positive_operators_gens()
            ....:                  for i in range(c.nrows())
            ....:                  for j in range(c.ncols())
            ....:                  if i != j)
            True

        The trivial cone in a trivial space has no cross-positive
        operators::

            sage: K = cones.trivial(0)
            sage: K.cross_positive_operators_gens()
            []

        Every operator is a cross-positive operator on the ambient
        vector space::

            sage: K = Cone([(1,), (-1,)])
            sage: K.is_full_space()
            True
            sage: K.cross_positive_operators_gens()
            [[1], [-1]]

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: K.cross_positive_operators_gens()
            [
            [1 0]  [-1  0]  [0 1]  [ 0 -1]  [0 0]  [ 0  0]  [0 0]  [ 0  0]
            [0 0], [ 0  0], [0 0], [ 0  0], [1 0], [-1  0], [0 1], [ 0 -1]
            ]

        A non-obvious application is to find the cross-positive
        operators on the right half-plane [Or2018b]_::

            sage: K = Cone([(1,0), (0,1), (0,-1)])
            sage: K.cross_positive_operators_gens()
            [
            [1 0]  [-1  0]  [0 0]  [ 0  0]  [0 0]  [ 0  0]
            [0 0], [ 0  0], [1 0], [-1  0], [0 1], [ 0 -1]
            ]

        Cross-positive operators on a subspace are Lyapunov-like and
        vice-versa::

            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: lls = span(vector(l.list())
            ....:            for l in K.lyapunov_like_basis())
            sage: cs  = span(vector(c.list())
            ....:            for c in K.cross_positive_operators_gens())
            sage: cs == lls
            True

        TESTS:

        The trivial cone, full space, and half-plane all give rise to
        the expected dimensions [Or2018b]_::

            sage: n = ZZ.random_element(5)
            sage: K = cones.trivial(n)
            sage: L = ToricLattice(n^2)
            sage: cp_gens = K.cross_positive_operators_gens()
            sage: cp_cone = Cone((g.list() for g in cp_gens),
            ....:                lattice=L,
            ....:                check=False)
            sage: cp_cone.dim() == n^2
            True

            sage: K = K.dual()
            sage: K.is_full_space()
            True
            sage: cp_gens = K.cross_positive_operators_gens()
            sage: cp_cone = Cone((g.list() for g in cp_gens),
            ....:                lattice=L,
            ....:                check=False)
            sage: cp_cone.dim() == n^2
            True

            sage: K = Cone([(1,0),(0,1),(0,-1)])
            sage: cp_gens = K.cross_positive_operators_gens()
            sage: cp_cone = Cone(( g.list() for g in cp_gens ), check=False)
            sage: cp_cone.dim() == 3
            True

        """
        # Compute the desired cone from its dual...
        cp_cone = self._cross_positive_operators_dual().dual()

        # Convert its rays from long vectors back to matrices.
        MS = MatrixSpace(self.lattice().base_field(), self.lattice_dim())
        return [MS(v.list()) for v in cp_cone]

    def Z_operators_gens(self):
        r"""
        Compute minimal generators of the Z-operators on this cone.

        The Z-operators on a cone generalize the Z-matrices over the
        nonnegative orthant. They are simply negations of the
        :meth:`cross_positive_operators_gens`.

        OUTPUT:

        A list of `n`-by-`n` matrices where `n` is the ambient dimension
        of this cone. Each matrix `L` in the list has the property that
        `s(L(x)) \le 0` whenever `(x,s)` is an element of this cone's
        :meth:`discrete_complementarity_set`.

        The returned matrices generate the cone of Z-operators on this
        cone; that is,

        - Any nonnegative linear combination of the returned matrices
          is a Z-operator on this cone.

        - Every Z-operator on this cone is some nonnegative linear
          combination of the returned matrices.

        .. SEEALSO::

           :meth:`cross_positive_operators_gens`,
           :meth:`lyapunov_like_basis`,
           :meth:`positive_operators_gens`

        REFERENCES:

        - [BP1994]_
        - [Or2018b]_

        """
        return [ -cp for cp in self.cross_positive_operators_gens() ]

    def max_angle(self, other=None, exact=True, epsilon=0):
        r"""
        Return the maximal angle between ``self`` and ``other``.

        The maximal angle between two closed convex cones is the
        unique largest angle formed by any two unit-norm vectors in
        those cones. In pathological cases, this computation can fail.

        If it fails when ``exact`` is ``True`` and if each of the
        cones :meth:`is_strictly_convex`, then a second attempt will
        be made using inexact arithmetic. (This sometimes avoids the
        problem noted in [Or2024]_). If the computation fails when the
        cones are not strictly convex or when ``exact`` is ``False``,
        a :exc:`ValueError` is raised.

        INPUT:

        - ``other`` -- (default: ``None``) a rational, polyhedral
          convex cone

        - ``exact`` -- boolean (default: ``True``); whether or not to use exact
          rational arithmetic instead of floating point computations;
          beware that ``True`` is not guaranteed to avoid floating
          point computations if the algorithm runs into trouble in
          rational arithmetic

        - ``epsilon`` -- (default: ``0``) the tolerance to use when
          making comparisons

        .. WARNING::

            Using inexact arithmetic (``exact=False``) is faster, but
            this computation is only known to be stable when both of
            the cones are strictly convex (or when one of them is the
            entire space, but the maximal angle is obviously `\pi` in
            that case).

        OUTPUT:

        A triple `\left( \theta_{\max}, u, v \right)`
        containing:

        - the maximal angle `\theta_{\max}` between ``self`` and
          ``other``

        - a vector `u` in ``self`` that achieves the maximal angle

        - a vector `v` in ``other`` that achieves the maximal angle

        If ``other`` is ``None`` (the default), then the maximal angle
        within this cone (between this cone and itself) is returned.

        If an eigenspace of dimension greater than one is encountered
        and if the corresponding angle cannot be ruled out as a
        maximum, the behavior of this function depends on
        ``exact``:

        - If ``exact`` is ``True`` and if both ``self`` and ``other``
          are strictly convex, then the algorithm will fall back to
          inexact arithmetic. In that case, the returned angle and
          vectors will be over ``RDF``.

        - If ``exact`` is ``False`` or if either cone is not strictly
          convex, then a :exc:`ValueError` is raised to indicate
          that we have failed; i.e. we cannot say with certainty what
          the maximal angle is.

        REFERENCES:

        - [IS2005]_
        - [SS2016]_
        - [Or2020]_
        - [Or2024]_

        ALGORITHM:

        Algorithm 3 in [Or2020]_ is used. If a potentially-maximal
        angle corresponds to an eigenspace of dimension two or more,
        we sometimes fall back to inexact arithmetic which has the
        effect of perturbing the cones. That this will not affect the
        answer too much is one conclusion of [Or2024]_.

        EXAMPLES:

        The maximal angle in a single ray is zero::

            sage: K = random_cone(min_rays=1, max_rays=1, max_ambient_dim=5)
            sage: K.max_angle()[0]
            0

        The maximal angle in the nonnegative octant is `\pi/2`::

            sage: K = cones.nonnegative_orthant(3)
            sage: K.max_angle()[0]
            1/2*pi

        The maximal angle between the nonnegative quintant and the
        Schur cone of dimension 5 is about `0.8524 \pi`. The same
        result can be obtained faster using inexact arithmetic, but
        only confidently so because we already know the answer::

            sage: # long time
            sage: P = cones.nonnegative_orthant(5)
            sage: Q = cones.schur(5)
            sage: actual = P.max_angle(Q)[0]
            sage: expected = 0.8524*pi
            sage: bool( (actual - expected).abs() < 0.0001 )
            True
            sage: actual = P.max_angle(Q,exact=False)[0]
            sage: bool( (actual - expected).abs() < 0.0001 )
            True

        The maximal angle within the Schur cone is known explicitly
        via Gourion and Seeger's Proposition 2 [GS2010]_::

            sage: n = 3
            sage: K = cones.schur(n)
            sage: bool(K.max_angle()[0] == ((n-1)/n)*pi)
            True

        Sage can't prove that the actual and expected results are
        equal in the next two cases without a little nudge in the
        right direction, and, moreover, it's crashy about it::

            sage: n = 4
            sage: K = cones.schur(n)
            sage: actual = K.max_angle()[0].simplify()._sympy_()._sage_()
            sage: expected = ((n-1)/n)*pi
            sage: bool( actual == expected )
            True
            sage: n = 5
            sage: K = cones.schur(n)
            sage: actual = K.max_angle()[0].simplify()._sympy_()._sage_()
            sage: expected = ((n-1)/n)*pi
            sage: bool( actual == expected )
            True

        When there's a unit norm vector in ``self`` whose negation is
        in ``other``, they form a maximal angle of `\pi`::

            sage: P = Cone([(5,1), (1,-1)])
            sage: Q = Cone([(-1,0), (-1,0)])
            sage: P.max_angle(Q)[0]
            pi

        TESTS:

        When ``self`` and ``-other`` have nontrivial intersection, we
        expect the maximal angle to be `\pi`::

            sage: # long time
            sage: from sage.geometry.cone_critical_angles import (
            ....:   _random_admissible_cone )
            sage: n = ZZ.random_element(1,3)
            sage: P = _random_admissible_cone(ambient_dim=n)
            sage: Q = _random_admissible_cone(ambient_dim=n)
            sage: expected = P.intersection(-Q).is_trivial()
            sage: actual = bool(P.max_angle(Q)[0] != pi)
            sage: actual == expected
            True

        In particular, this should happen when either cone is the full
        space::

            sage: from sage.geometry.cone_critical_angles import (
            ....:   _random_admissible_cone )
            sage: n = ZZ.random_element(1,5)
            sage: P = _random_admissible_cone(ambient_dim=n)
            sage: Q = cones.trivial(n, P.dual().lattice()).dual()
            sage: Q.is_full_space()
            True
            sage: P.max_angle(Q)[0]
            pi
            sage: Q.max_angle(P)[0]
            pi
        """
        # We do the argument checking here, in the public cone method,
        # because the error message should say something like "this
        # cone" and not reference the argument of a function the user
        # has never heard of in some internal module.
        if self.is_trivial():
            raise ValueError("cone should not be trivial")

        if other is None:
            other = self
        else:
            if (other.lattice_dim() != self.lattice_dim()):
                raise ValueError("lattice dimensions of self and other "
                                 "must agree")
            if other.is_trivial():
                raise ValueError("other cone cannot be trivial")

        from sage.geometry.cone_critical_angles import max_angle
        return max_angle(self, other, exact, epsilon)

    def irreducible_factors(self):
        r"""
        Decompose a strictly convex (AKA pointed) convex cone
        into nontrivial irreducible factors.

        A convex cone is said to be *reducible* if it can be expressed as
        the direct sum of two subcones. An *irreducible* cone is a cone
        that is not reducible. Every cone that :meth:`is_proper` can be
        uniquely decomposed -- up to isomorphism and the order of the
        factors -- as a direct sum of irreducible, pointed, nontrivial
        factors (subcones).

        In the literature, "reducible" and "decomposable" are
        interchangeable.

        The cone being strictly convex / pointed is essential to the
        uniqueness of the decomposition: the plane is a cone, and it can
        be decomposed into the direct sum of the x-axis and y-axis, but
        any other pair of (unequal) lines through the origin would work
        just as well.

        Being solid is less essential. The definition of "direct sum"
        requires that the ambient space be fully decomposed into a set of
        subspaces, and if the cone is not solid, we must (to satisfy the
        definition) manufacture trivial cones to use as the factors in the
        subspaces orthogonal to the given cone. This destroys the
        uniqueness of the direct sum decomposition in the sense that the
        vector subspaces corresponding to the trivial factors are
        arbitrary. If the given cone is one-dimensional and the ambient
        space three-dimensional, then a full decomposition would include
        two trivial cones in any two one-dimensional subspaces of the
        plane. As in the preceding paragraph, those one-dimensional
        subspaces can be chosen arbitrarily.

        Considering that this method does not return the ambient vector
        space decomposition, adding trivial cones to the list of
        irreducible factors is not helpful: any cone is equal to the sum
        of itself with a trivial cone, or two trivial cones, etc. This is
        all to say: this method will accept non-solid cones, but it will
        not return any trivial factors. The result does not technically
        correspond to a direct sum, but by omitting the trivial factors,
        we recover uniqueness of the factors in the non-solid case.

        OUTPUT:

        A set of ``sage.geometry.cone.ConvexRationalPolyhedralCone``,
        each of which is nontrivial, strictly convex (pointed), and
        irreducible. Each factor lives in the same ambient space as
        ``K`` to avoid confusing isomorphic factors that live in
        different subspaces. As a consequence, factors of solid cones
        will not in general be solid.

        ALGORITHM:

        If a strictly convex / pointed cone ``K`` is reducible to ``K1 +
        K2`` in the vector space ``V = V1 + V2``, then the generators
        (extreme rays) of ``K`` can be split into subsets ``G1`` and
        ``G2`` of ``V1`` and ``V2`` that generate ``K1`` and ``K2``,
        respectively. Following [BSPRS2014]_, we find a basis for the
        span of ``K`` consisting of extreme rays of ``K``, and then
        express each extreme ray in terms of that basis. By looking for
        groups of rays that require a common subset of basis elements, we
        determine which, if any, generators can be split accordingly.

        REFERENCES:

        - [HausGul2002]_
        - [BSPRS2014]_

        .. SEEALSO::

            :meth:`is_reducible`

        EXAMPLES:

        The nonnegative orthant is a direct sum of its extreme rays::

            sage: K = cones.nonnegative_orthant(3)
            sage: expected = {Cone([r]) for r in K.rays()}
            sage: K.irreducible_factors() == expected
            True

        An irreducible example (the l1-norm cone)::

            sage: K = Cone([(1,0,1), (0,1,1), (-1,0,1), (0,-1,1)])
            sage: K.irreducible_factors() == {K}
            True

        We can decompose reducible cones that are not solid::

            sage: K = Cone([(1,0,0), (0,1,0)])
            sage: expected = {Cone([r]) for r in K.rays()}
            sage: K.irreducible_factors() == expected
            True

        And nothing bad happens with irreducible ones::

            sage: K = Cone([(1,0,1,0), (0,1,1,0), (-1,0,1,0), (0,-1,1,0)])
            sage: K.irreducible_factors() == {K}
            True

        The rays of the Schur cone are linearly-independent::

            sage: K = cones.schur(5)
            sage: expected = {Cone([r]) for r in K.rays()}
            sage: K.irreducible_factors() == expected
            True

        The Cartesian product can be used to construct larger
        examples. Here, two of the factors are irreducible, but the
        nonnegative orthant should reduce to two rays::

            sage: J1 = Cone([(1,0,1,0), (0,1,1,0), (-1,0,1,0), (0,-1,1,0)])
            sage: J2 = cones.nonnegative_orthant(2)
            sage: J3 = cones.rearrangement(2, 5)
            sage: J1
            3-d cone in 4-d lattice N
            sage: J3
            5-d cone in 5-d lattice N
            sage: K = J1.cartesian_product(J2).cartesian_product(J3)
            sage: sorted(K.irreducible_factors())
            [5-d cone in 11-d lattice N+N+N,
             1-d cone in 11-d lattice N+N+N,
             1-d cone in 11-d lattice N+N+N,
             3-d cone in 11-d lattice N+N+N]

        All proper cones in two dimensions are isomorphic to the
        nonnegative orthant and should decompose::

            sage: K = Cone([(1,2), (-3,4)])
            sage: sorted(K.irreducible_factors())
            [1-d cone in 2-d lattice N, 1-d cone in 2-d lattice N]

        In three dimensions, we can hit the nonnegative orthant with an
        invertible map, and it will still decompose::

            sage: A = matrix(QQ, [[1, -1/2, 0], [1, 1, -2], [-2, 0, -1]])
            sage: K = cones.nonnegative_orthant(3)
            sage: AK = Cone([ r*A for r in K.rays() ])
            sage: expected = {Cone([r]) for r in AK.rays()}
            sage: AK.irreducible_factors() == expected
            True

        This example decomposes into a ray and a cone with four
        generators in a three-dimensional subspace that trivially
        intersects the span of the ray::

            sage: K = Cone([ (-1, -1, 0, -1),
            ....:            (-1,  0, 1,  0),
            ....:            ( 0,  1, 4, -2),
            ....:            ( 1,  0, 0, -1),
            ....:            ( 1,  1, 0, -1) ])
            sage: K1, K2 = sorted(K.irreducible_factors())
            sage: K1
            3-d cone in 4-d lattice N
            sage: K1.rays()
            N(-1, -1, 0, -1),
            N( 1,  0, 0, -1),
            N( 0,  1, 4, -2),
            N(-1,  0, 1,  0)
            in 4-d lattice N
            sage: K2
            1-d cone in 4-d lattice N
            sage: K2.rays()
            N(1, 1, 0, -1)
            in 4-d lattice N
            sage: K1.span().intersection(K2.span())
            Sublattice <>

        TESTS:

        The error message looks right::

            sage: K = Cone([(1,0),(-1,0)])
            sage: K.irreducible_factors()
            Traceback (most recent call last):
            ...
            ValueError: cone must be strictly convex (AKA pointed) for its
            irreducible factors to be well-defined

        Trivial cones are handled correctly::

            sage: K = cones.trivial(0)
            sage: K.irreducible_factors() == {K}
            True
            sage: K = cones.trivial(4)
            sage: K.irreducible_factors() == {K}
            True
        """
        if not self.is_strictly_convex():
            raise ValueError("cone must be strictly convex (AKA pointed) for"
                             " its irreducible factors to be well-defined")

        if self.is_trivial():
            # Trivial cones are valid inputs, but they have no generators
            # so a special case is required.
            return {self}

        V = self.ambient_vector_space()

        # Create a column matrix from the generators of K, and then
        # perform Gaussian elimination so that the resulting pivots
        # specify a linearly-independent subset of the original columns.
        M = self.rays().column_matrix().change_ring(V.base_ring())
        pivots = M.echelon_form(algorithm="classical").pivots()
        M = M.matrix_from_columns(pivots)

        # A basis for span(self)
        B = M.columns(copy=False)

        # The span of self, but now with a basis of extreme rays. If
        # self is not solid, B will not be a basis of the entire space
        # V, only of a subspace.
        W = V.subspace_with_basis(B)

        # Make a graph with ray -> ray edges where the coefficients in
        # the W-representation nonzero. Beware: while they ultimately
        # represent rays of self, the W-coordinates have been
        # renumbered by pivots().
        vertices = list(range(self.nrays()))
        edges = [(i, pivots[j])
                 for i in vertices
                 for j in W.coordinate_vector(self.ray(i)).nonzero_positions()
                 if pivots[j] != i]
        from sage.graphs.graph import Graph
        G = Graph([vertices, edges], format='vertices_and_edges')

        from sage.geometry.cone import Cone
        if G.connected_components_number() == 1:
            # Special case where we don't want to pointlessly return an
            # equivalent but unequal copy of the input cone.
            return {self}
        return {Cone(self.rays(c)) for c in G.connected_components()}

    def is_reducible(self) -> bool:
        r"""
        Return whether or not this cone is reducible.

        A pointed convex cone is reducible if some other cone appears
        in its :meth:`irreducible_factors`.

        .. SEEALSO::

            :meth:`irreducible_factors`

        EXAMPLES:

        The nonnegative orthant is always reducible in dimension two or
        more::

            sage: cones.nonnegative_orthant(1).is_reducible()
            False
            sage: cones.nonnegative_orthant(2).is_reducible()
            True
            sage: cones.nonnegative_orthant(3).is_reducible()
            True

        Theorem 4.1 in [JG2018]_ says that the Lyapunov rank of a
        permutation-invariant cone is either ``n`` or ``1``, depending
        on whether or not it is reducible. Corollary 5.2.4 of
        [Jeong2017]_ then implies that the ``(p,n)`` rearrangement
        cone is reducible if and only if ``p`` is either ``1`` or
        ``n-1``. We exclude the possibility of ``p == n`` since that
        returns a (not pointed) half-space::

            sage: n = ZZ.random_element(10) + 2
            sage: p = ZZ.random_element(1, n)
            sage: K = cones.rearrangement(p, n)
            sage: K.is_reducible() == (p in [1, n-1])
            True

        """
        return len(self.irreducible_factors()) > 1


def random_cone(lattice=None, min_ambient_dim=0, max_ambient_dim=8,
                min_rays=0, max_rays=16, strictly_convex=None, solid=None):
    r"""
    Generate a random convex rational polyhedral cone.

    Lower and upper bounds may be provided for both the dimension of
    the ambient space and the number of generating rays that the cone
    will possess. If lower bounds are not provided, they default to
    zero. Unspecified upper bounds default to reasonable values. You
    can additionally request that the cone be strictly convex (or
    not) or solid (or not).

    As an alternative to the dimension bounds, you can provide a
    ``lattice`` for the returned cone to live in. In that case, the
    ``min_ambient_dim`` and ``max_ambient_dim`` parameters are
    ignored.

    .. WARNING::

        If you request a large number of rays in a low-dimensional
        space, you may be waiting for a while.

        In three dimensions, and ignoring for the moment the need for
        rational coordinates, it is possible to obtain an "octagon"
        raised up to height one (all z-coordinates equal to one). But
        in practice we usually generate the entire three-dimensional
        ambient vector space (having six rays) before we obtain the
        eight that we want. At that point it is futile to add more
        rays; we have to throw the cone out and start from
        scratch. This process repeats until we get lucky.

    INPUT:

    - ``lattice`` -- (default: random) a ``ToricLattice`` object in
      which the returned cone will live. By default a new lattice will
      be constructed with a randomly-chosen rank (subject to
      ``min_ambient_dim`` and ``max_ambient_dim``).

    - ``min_ambient_dim`` -- (default: zero) a nonnegative integer
      representing the minimum dimension of the ambient lattice;
      ignored if ``lattice`` is given explicitly

    - ``max_ambient_dim`` -- (default: 8) a nonnegative integer
      representing the maximum dimension of the ambient lattice;
      ignored if ``lattice`` is given explicitly

    - ``min_rays`` -- (default: zero) a nonnegative integer representing
      the minimum number of generating rays of the cone

    - ``max_rays`` -- (default: 16) a nonnegative integer representing
      the maximum number of generating rays of the cone.

    - ``strictly_convex`` -- (default: random) whether or not to make the
      returned cone strictly convex. Specify ``True`` for a strictly convex
      cone, ``False`` for a non-strictly-convex cone, or ``None`` if you
      don't care.

    - ``solid`` -- (default: random) whether or not to make the returned
      cone solid. Specify ``True`` for a solid cone, ``False`` for a
      non-solid cone, or ``None`` if you don't care.

    The ``max_ambient_dim`` and ``max_rays`` should be chosen
    carefully. While in theory it is possible for `n+1` rays to
    generate an `n`-dimensional vector space, if you actually try it,
    they will normalized to a set containing `2n` rays (plus/minus a
    basis). If the ambient space is to be generated randomly,
    ``max_rays`` must therefore be at least twice ``max_ambient_dim``.

    On the other hand, the higher ``max_rays`` is, the less likely we
    are in a given iteration (see the ALGORITHM) to find a satisfatory
    cone. As we attempt to find more than `2n` independent rays in an
    `n`-dimensional space, it becomes increasingly likely that we will
    generate the entire vector space, fail, and be forced start
    over. While there do exist pointed cones (solid or not) in `n \ge
    3` dimensions with more than `2n` rays, they are unusually hard to
    find, and not comparably interesting.

    OUTPUT:

    In most cases, a new, randomly generated
    :class:`ConvexRationalPolyhedralCone` will be returned. A
    :exc:`ValueError` will however be raised under the following
    conditions:

    * Any of ``min_ambient_dim``, ``max_ambient_dim``, ``min_rays``,
      or ``max_rays`` are negative.

    * ``max_ambient_dim`` is less than ``min_ambient_dim``.

    * ``max_rays`` is less than ``min_rays``.

    * A non-strictly-convex cone is requested but ``max_rays`` is less
      than two.

    * ``max_ambient_dim`` is two or less, and ``min_rays`` is
      greater than ``2*max_ambient_dim``.

    * ``lattice`` has dimension two or less, and ``min_rays`` is
      larger than twice that dimension.

    * A strictly-convex cone is requested with ``max_ambient_dim <=
      2``, and ``min_rays`` is greater than ``max_ambient_dim``.

    * A strictly-convex cone is requested in a ``lattice`` that has
      dimension two or less, and ``min_rays`` is greater than that
      dimension.

    * A non-strictly-convex cone is requested and ``max_ambient_dim``
      is zero.

    * A non-strictly-convex cone is requested in a trivial lattice.

    * A non-solid cone is requested but ``max_ambient_dim`` is zero.

    * A non-solid cone is requested but ``lattice`` has dimension
      zero.

    * A solid cone is requested but ``max_rays`` is less than
      ``min_ambient_dim``.

    * A solid cone is requested but ``max_rays`` is less than the
      dimension of ``lattice``.

    ALGORITHM:

    Loop until a valid cone is found.

    The first step inside the loop is to decide on a lattice. If a
    ``lattice`` was supplied, it is used. Otherwise, we construct a new
    lattice whose dimension is sampled randomly from the interval
    ``[min_ambient_dim, max_ambient_dim]``.

    Next we decide the number of rays to seek in this iteration. This
    number is chosen randomly from ``[min_rays, max_rays]``. In both
    cases (dimension and number of rays), minor adjustments may be
    made to the interval to avoid impossible combinations.

    With the lattice dimension and number of rays fixed, random
    lattice elements (rays) are generated in batches and added to a
    cone. This continues until the cone either has enough rays or
    becomes the entire ambient space, at which point it is futile to
    append more.

    If the cone is equal to the ambient space, we check to see if it
    meets the user's requirements. If it does, it is returned; if not,
    we throw it away and start over with a new (random, admissible)
    dimension and number of rays in mind.

    If the cone is *not* equal to the ambient space, and if ``solid``
    or ``strictly_convex`` is not ``None``, we attempt some
    adjustments:

    * make a non-solid cone solid by adding more rays,

    * make a solid cone non-solid by embedding it in a space
      of higher dimension,

    * make a pointed cone non-pointed by appending negative
      multiples of existing rays (if there's room for more rays),

    * make a pointed cone non-pointed by replacing existing rays
      with negative multiples of others (if we are at max_rays).

    (We do *not* use
    :meth:`ConvexRationalPolyhedralCone.solid_restriction` to make a
    non-solid cone solid. The solid restriction does not affect the
    ray count, but it does involve a change of basis that turns your
    existing rays into boring standard basis vectors.)

    Afterwards the cone is checked against the requirements. If it is
    valid, we return it; otherwise we throw it out and begin anew.

    Building a cone that meets the user's critereia is inherently a
    bottom-up process. We cannot know how many rays in a given set are
    extreme until we reduce them with :func:`Cone`, but in the
    current implementation, reducing the rays also happens to
    normalize them. If we desire ``r`` rays and generate ``r`` at
    random, we are likely to discover that some are redundant and that
    we need to append more. On the other hand, if generate more than
    ``r`` with the expectation that some will be eliminated, we
    increase the likelihood that :func:`Cone` normalizes them to
    boring standard basis vectors.

    This bottom-up approach motivates fixing the number of rays ``r``
    attempted in each iteration. If ``r`` were not fixed -- if we were
    willing to accept any number of rays between the min/max in any
    iteration -- then we would almost always stop as soon as we hit
    ``min_rays``. That would be legal, but perhaps unexpected: it is
    reasonable to assume that the result will have a number of rays
    distributed throughout ``[min_rays, max_rays]``. That said, this
    approach greatly increases the failure rate of the loop.

    The default upper bounds on the ambient dimension and number of
    rays ensure that failed iterations remain fast. Adding rays to the
    cone consists of appending vectors to the ray set, and calling
    :func:`Cone` on the result. This can be quite slow, and we may
    have to do it many times. As the speed of :func:`Cone` depends
    mainly on the degree and number of rays given to it, those are the
    critical factors in the worst-case performance.

    For the time being, we forego the following adjustments because
    they remove the cone's
    :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.lines`,
    thereby skewing the number of rays towards the lower end of
    ``[min_rays, max_rays]``:

    * Using :meth:`ConvexRationalPolyhedralCone.strict_quotient` to
      make a non-pointed cone pointed.

    * Projecting a non-pointed cone onto the orthogonal complement
      of its :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.lines`
      to obtain a
      pointed cone.

    EXAMPLES:

    Generate a trivial cone in a trivial space::

        sage: random_cone(max_ambient_dim=0, max_rays=0)
        0-d cone in 0-d lattice N

    We can predict the ambient dimension when
    ``min_ambient_dim == max_ambient_dim``::

        sage: K = random_cone(min_ambient_dim=4, max_ambient_dim=4)
        sage: K.lattice_dim()
        4

    Likewise for the number of rays when ``min_rays == max_rays``::

        sage: K = random_cone(min_rays=3, max_rays=3)
        sage: K.n_rays()
        3

    If we specify a lattice, then the returned cone will live in it::

        sage: L = ToricLattice(5, "L")
        sage: K = random_cone(lattice=L)
        sage: K.lattice() is L
        True

    We can also request a strictly convex cone::

        sage: K = random_cone(max_ambient_dim=8, max_rays=10,
        ....:                 strictly_convex=True)
        sage: K.is_strictly_convex()
        True

    Or one that isn't strictly convex::

        sage: K = random_cone(min_ambient_dim=5, min_rays=2,
        ....:                 max_ambient_dim=8, max_rays=10,
        ....:                 strictly_convex=False)
        sage: K.is_strictly_convex()
        False

    An example with all parameters set::

        sage: K = random_cone(min_ambient_dim=4, max_ambient_dim=7,
        ....:                 min_rays=2, max_rays=10,
        ....:                 strictly_convex=False, solid=True)
        sage: 4 <= K.lattice_dim() and K.lattice_dim() <= 7
        True
        sage: 2 <= K.n_rays() and K.n_rays() <= 10
        True
        sage: K.is_strictly_convex()
        False
        sage: K.is_solid()
        True

    TESTS:

    It's hard to test the output of a random process, but we can at
    least make sure that we get a cone back::

        sage: K = random_cone()
        sage: isinstance(K, sage.geometry.abc.ConvexRationalPolyhedralCone)
        True

    The upper/lower bounds are respected::

        sage: K = random_cone(min_ambient_dim=5, max_ambient_dim=8,
        ....:                 min_rays=3, max_rays=4)
        sage: 5 <= K.lattice_dim() and K.lattice_dim() <= 8
        True
        sage: 3 <= K.n_rays() and K.n_rays() <= 4
        True

    Ensure that an exception is raised when either lower bound is greater
    than its respective upper bound::

        sage: random_cone(min_ambient_dim=5, max_ambient_dim=2)
        Traceback (most recent call last):
        ...
        ValueError: max_ambient_dim cannot be less than min_ambient_dim.

        sage: random_cone(min_rays=5, max_rays=2)
        Traceback (most recent call last):
        ...
        ValueError: max_rays cannot be less than min_rays.

    If the user requests too many rays in zero, one, or two dimensions,
    a :exc:`ValueError` is thrown::

        sage: random_cone(max_ambient_dim=0, min_rays=1)
        Traceback (most recent call last):
        ...
        ValueError: all cones in dimension d <= 2 have 2d or fewer rays.
        Please increase max_ambient_dim or decrease min_rays.

        sage: random_cone(max_ambient_dim=1, min_rays=3)
        Traceback (most recent call last):
        ...
        ValueError: all cones in dimension d <= 2 have 2d or fewer rays.
        Please increase max_ambient_dim or decrease min_rays.

        sage: random_cone(max_ambient_dim=2, min_rays=5)
        Traceback (most recent call last):
        ...
        ValueError: all cones in dimension d <= 2 have 2d or fewer rays.
        Please increase max_ambient_dim or decrease min_rays.

    This happens with a lattice, too::

        sage: L = ToricLattice(0)
        sage: random_cone(lattice=L, min_rays=1)
        Traceback (most recent call last):
        ...
        ValueError: all cones in this lattice have 0 or fewer rays.
        Please decrease min_rays.

        sage: L = ToricLattice(1)
        sage: random_cone(lattice=L, min_rays=3)
        Traceback (most recent call last):
        ...
        ValueError: all cones in this lattice have 2 or fewer rays.
        Please decrease min_rays.

        sage: L = ToricLattice(2)
        sage: random_cone(lattice=L, min_rays=5)
        Traceback (most recent call last):
        ...
        ValueError: all cones in this lattice have 4 or fewer
        rays. Please decrease min_rays.

    Ensure that we can obtain a cone in three dimensions with a large
    number (in particular, more than 2*dim) rays. Note that at least
    dimension three is required for seven rays; we should throw out
    cones in zero/one/two dimensions until we get lucky. We set the
    random seed manually to ensure that the test completes in a
    reasonable amount of time::

        sage: # long time
        sage: _initial_seed = initial_seed()
        sage: set_random_seed(12)
        sage: K = random_cone(max_ambient_dim=3, min_rays=7)
        sage: K.n_rays() >= 7
        True
        sage: K.lattice_dim()
        3
        sage: set_random_seed(_initial_seed)

    It is an error to request a non-strictly-convex cone in a trivial
    lattice::

        sage: random_cone(max_ambient_dim=0, strictly_convex=False)
        Traceback (most recent call last):
        ...
        ValueError: all cones are strictly convex when max_ambient_dim
        is zero.

        sage: L = ToricLattice(0, "L")
        sage: random_cone(lattice=L, strictly_convex=False)
        Traceback (most recent call last):
        ...
        ValueError: all cones in the trivial lattice are strictly
        convex (trivial).

    Or a non-strictly-convex cone with fewer than two rays::

        sage: random_cone(max_rays=1, strictly_convex=False)
        Traceback (most recent call last):
        ...
        ValueError: all cones are strictly convex when max_rays is
        less than two.

        sage: L = ToricLattice(5)
        sage: random_cone(lattice=L, max_rays=1, strictly_convex=False)
        Traceback (most recent call last):
        ...
        ValueError: all cones are strictly convex when max_rays is
        less than two.

    But fine to ask for a strictly convex trivial cone::

        sage: L = ToricLattice(0, "L")
        sage: random_cone(lattice=L, strictly_convex=True)
        0-d cone in 0-d lattice L

    A :exc:`ValueError` is thrown if a non-solid cone is requested in a
    zero-dimensional lattice::

        sage: L = ToricLattice(0)
        sage: random_cone(lattice=L, solid=False)
        Traceback (most recent call last):
        ...
        ValueError: all cones in the trivial lattice are solid.

        sage: random_cone(max_ambient_dim=0, solid=False)
        Traceback (most recent call last):
        ...
        ValueError: all cones are solid when max_ambient_dim is zero.

    A :exc:`ValueError` is thrown if a solid cone is requested but the
    maximum number of rays is too few::

        sage: random_cone(min_ambient_dim=4, max_rays=3, solid=True)
        Traceback (most recent call last):
        ...
        ValueError: max_rays must be at least min_ambient_dim for a
        solid cone.

        sage: L = ToricLattice(5)
        sage: random_cone(lattice=L, max_rays=3, solid=True)
        Traceback (most recent call last):
        ...
        ValueError: max_rays must be at least 5 for a solid cone in this
        lattice.

    A :exc:`ValueError` is thrown if a cone with too many rays is
    requested in dimension less than three::

        sage: K = random_cone(max_ambient_dim=1,
        ....:                 min_rays=3,
        ....:                 max_rays=3)
        Traceback (most recent call last):
        ...
        ValueError: all cones in dimension d <= 2 have 2d or fewer
        rays. Please increase max_ambient_dim or decrease min_rays.

    A :exc:`ValueError` is thrown if a strictly convex cone is
    requested with too many rays in dimension less than three (the
    bound is lower than if we did not specify
    ``strictly_convex=True``)::

        sage: K = random_cone(max_ambient_dim=2,
        ....:                 min_rays=3,
        ....:                 max_rays=3,
        ....:                 strictly_convex=True)
        Traceback (most recent call last):
        ...
        ValueError: in dimension d <= 2, all strictly convex cones
        have d rays. Please increase max_ambient_dim, or decrease
        min_rays.

    Ensure that we can produce a cone whose
    :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.lines` are not
    orthogonal::

        sage: set_random_seed(1)
        sage: K = random_cone(min_ambient_dim=2, strictly_convex=False)
        sage: V = K.ambient_vector_space()
        sage: all( V(x).inner_product(V(y)).is_zero()
        ....:      for x in K.lines()
        ....:      for y in K.lines()
        ....:      if not x == y )
        False

    """
    from sage.misc.prandom import randint, sample

    # Catch obvious mistakes so that we can generate clear error
    # messages. We sanity check the min_rays and max_rays parameters
    # first, since they do not depend on whether or not a lattice
    # was specified.

    if min_rays < 0:
        raise ValueError("min_rays must be nonnegative.")

    if max_rays < 0:
        raise ValueError("max_rays must be nonnegative.")

    if max_rays < min_rays:
        raise ValueError("max_rays cannot be less than min_rays.")

    if strictly_convex is False and max_rays < 2:
        raise ValueError("all cones are strictly convex when "
                         "max_rays is less than two.")

    if lattice is None:
        if min_ambient_dim < 0:
            raise ValueError("min_ambient_dim must be nonnegative.")

        if max_ambient_dim < 0:
            raise ValueError("max_ambient_dim must be nonnegative.")

        if max_ambient_dim < min_ambient_dim:
            raise ValueError("max_ambient_dim cannot be less than "
                             "min_ambient_dim.")

        # The next check prevents an infinite loop (a futile search
        # for more rays) in zero, one, or two dimensions.
        if max_ambient_dim <= 2:
            if min_rays > 2*max_ambient_dim:
                raise ValueError("all cones in dimension d <= 2 have "
                                 "2d or fewer rays. Please increase "
                                 "max_ambient_dim or decrease "
                                 "min_rays.")

            elif strictly_convex and min_rays > max_ambient_dim:
                raise ValueError("in dimension d <= 2, all strictly "
                                 "convex cones have d rays. Please "
                                 "increase max_ambient_dim, or "
                                 "decrease min_rays.")

            if max_ambient_dim == 0:
                if strictly_convex is False:
                    raise ValueError("all cones are strictly convex "
                                     "when max_ambient_dim is zero.")
                if solid is False:
                    raise ValueError("all cones are solid when "
                                     "max_ambient_dim is zero.")

        if solid and max_rays < min_ambient_dim:
            raise ValueError("max_rays must be at least "
                             "min_ambient_dim for a solid cone.")
    else:
        # Also perform the "futile search" checks when a lattice is
        # given, using its dimension rather than max_ambient_dim.
        d = lattice.dimension()

        if d <= 2:
            if min_rays > 2*d:
                raise ValueError("all cones in this lattice have "
                                 f"{2*d} or fewer rays. Please "
                                 "decrease min_rays.")
            elif strictly_convex and min_rays > d:
                raise ValueError("all strictly convex cones in "
                                 f"this lattice have {d} or fewer "
                                 "rays. Please decrease min_rays.")

            if d.is_zero():
                if strictly_convex is False:
                    raise ValueError("all cones in the trivial lattice "
                                     "are strictly convex (trivial).")
                if solid is False:
                    raise ValueError("all cones in the trivial lattice "
                                     "are solid.")

        if solid and max_rays < d:
            raise ValueError(f"max_rays must be at least {d} for a "
                             "solid cone in this lattice.")


    # Parameter adjustment. In some cases are are able to adjust the
    # dim/ray bounds to eliminate impossible combinations. We are
    # going to attempt those combinations at random, so avoiding the
    # ones that are doomed to fail should improve the average runtime.

    if solid and max_rays < max_ambient_dim:
        # If the user wants a solid cone, we need at least as many
        # rays as the dimension. If max_rays was specified, then
        # above we have already verified that
        #
        # (a) max_rays >= min_ambient_dim
        # (b) max_ambient_dim >= min_ambient_dim, if max_ambient_dim
        #     is set
        #
        # As a result we can fudge the max_ambient_dim down to
        # max_rays, ensuring that we don't perform any pointless
        # iterations in a space that's too big for max_rays to
        # generate a solid cone.
        max_ambient_dim = max_rays

    if strictly_convex is False and min_rays < 2:
        # We've already verified that max_rays >= 2. Increase min_rays
        # to avoid pointless loops with r=0 or r=1
        min_rays = 2

    def is_valid(K):
        r"""
        Check if the given cone is valid; that is, if its ambient
        dimension and number of rays meet the upper and lower bounds
        provided by the user.
        """
        return all((
          lattice is None or K.lattice() is lattice,
          lattice is not None
            or min_ambient_dim <= K.lattice_dim() <= max_ambient_dim,
          min_rays <= K.n_rays() <= max_rays,
          solid is None or K.is_solid() == solid,
          strictly_convex is None
            or K.is_strictly_convex() == strictly_convex
        ))

    # Now we actually compute the thing. To avoid recursion (and the
    # associated "maximum recursion depth exceeded" error), we loop
    # until we have a valid cone and occasionally throw everything out
    # and start over from scratch.
    while True:
        # d = lattice.dimension() is already set above if
        # the lattice is not None
        L = lattice

        if lattice is None:
            # No lattice given, make our own.
            d = randint(min_ambient_dim, max_ambient_dim)
            L = ToricLattice(d)

        _min_rays_this_iter = min_rays
        if solid:
            # If a solid cone was requested, then before the loop, we
            # set max_ambient_dim to max_rays as a performance hint.
            # Without more information, that was the best we could do,
            # but now we know the dimension that we are seeking for
            # this iteration. If min_rays it too low to achieve it,
            # we can use a higher lower bound (but still within
            # max_rays). Note: the aforementioned hint guarantees
            # that d <= max_rays.
            _min_rays_this_iter = d

        # The number of rays that we will try to attain in this iteration.
        r = randint(_min_rays_this_iter, max_rays)

        # The rays are trickier to generate, since we could generate v and
        # 2*v as our "two rays." In that case, the resulting cone would
        # have only one generating ray -- not what we want.
        #
        # Let's begin with an easier question: how many rays should we
        # start with? If we want to attain r rays in this iteration,
        # then surely r is a good number to start with?
        rays = [L.random_element() for _ in range(r)]

        # The lattice parameter is required when no rays are given, so
        # we pass it in case r == 0 or d == 0 (or d == 1 but we're
        # making a strictly convex cone).
        K = Cone(rays, lattice=L)

        # Now, some of the rays that we generated were probably
        # redundant, so we need to come up with more. We can obviously
        # stop if K becomes the entire ambient vector space.
        #
        # We're still not guaranteed to have the correct number of
        # rays after this! Since we normalize the generators in the
        # loop above, we can jump from two to four generators by
        # adding e.g. (1,1) to [(0,1), (0,-1)]. Rather than trying to
        # mangle what we have, we just start over if we get a cone
        # that won't work.
        while True:
            # How many more rays do we need? *At least* this many.
            ray_gap = r - K.n_rays()

            if solid:
                # But if we need a solid cone, and if the
                # corresponding "dimension gap" is larger than the
                # ray_gap, we should use the dimension gap instead
                # because each additional ray can occupy at most one
                # vacant dimension.
                ray_gap = max(ray_gap, d - K.dim())

            # ray_gap can be negative if e.g. r=3 rays that generate
            # R^2 get normalized to four plus/minus basis vectors.
            if ray_gap <= 0 or K.is_full_space():
                break

            rays = list(K.rays())
            rays.extend(L.random_element() for _ in range(ray_gap))
            K = Cone(rays, lattice=L)

        if K.is_full_space():
            # When K is the full space, Cone() normalizes its rays to
            # plus/minus the standard basis. If the full space is an
            # acceptable result, OK, whatever, return it. But
            # otherwise, do not try to tweak the world's most
            # uninteresting set of generators.
            if is_valid(K):
                return K
            continue

        if solid and not K.is_solid():
            # Just add rays until we get a solid cone? Adding rays is
            # allowed. The solid_restriction() would also give us a
            # solid cone with the same number of rays, but it does a
            # change of basis that results in extremely boring
            # generators.
            while K.n_rays() <= max_rays and not K.is_solid():
                dim_gap = d - K.dim()
                rays = list(K.rays())
                rays.extend(L.random_element() for _ in range(dim_gap))
                K = Cone(rays, lattice=L)
            if not K.is_solid():
                # We had to stop because we hit max_rays
                continue

        if solid is False and K.is_solid():
            # To go from not-solid to solid, we add rays. But to go
            # from solid to not-solid? Removing rays sounds like a
            # good idea at first, but thanks to the normalization that
            # happens in Cone(), you are likely to wind up with more
            # boring rays than you would have otherwise. (The ambient
            # space is a solid cone, but every cone obtained by
            # deleting generators from it has only orthogonal standard
            # basis vectors as generators.) In any case, that's why
            # there is no fallback here for when K is already of
            # maximum dimension or a lattice was specified.
            if lattice is None and K.dim() < max_ambient_dim:
                # Embed what we've got into a space of one greater
                # dimension. Afterward we hit the cone with a random
                # unitary matrix, so that the form of the result is
                # not so predictable.
                K = Cone([r.list() + [0] for r in K.rays()],
                         check=False)
                A = matrix.random(K.lattice().base_field(),
                                  K.lattice_dim(),
                                  algorithm="unitary")
                K = Cone([A*r.dense_vector() for r in K.rays()],
                         check=False)
                L = K.lattice()

        if strictly_convex is False and K.is_strictly_convex():
            # To make a pointed cone not-pointed is easy: add a
            # negative multiple of an existing ray. Or if you're at
            # max_rays already, replace an existing ray. Below we do
            # this with a random number of new rays or replacements.
            ray_slots = max_rays - K.n_rays()
            if ray_slots == 0:
                # Replace some rays with flipped copies of other rays.
                if K.n_rays() >= 2:
                    # Divide by four because highly un-pointed cones
                    # are also un-interesting. (We have to divide by
                    # at least two: if we flip half, we need the other
                    # half for storage.)
                    num_to_flip = randint(1, max(1, K.n_rays()//4))
                    rays = [
                        -K.ray(K.n_rays() - i - 1) if i <= num_to_flip
                        else K.ray(i)
                        for i in range(K.n_rays())
                    ]
                    K = Cone(rays, lattice=L)
            elif K.n_rays() >= 1:
                # Append negative multiples of some existing rays to
                # guarantee that the cone contains lines.
                num_to_flip = randint(1, min(K.n_rays(), ray_slots))
                rays = list(K.rays())
                rays.extend( -r for r in sample(rays, num_to_flip) )
                K = Cone(rays, lattice=L)

        if is_valid(K):
            # Loop if we don't have a valid cone.
            return K
