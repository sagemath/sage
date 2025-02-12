# sage.doctest: needs sage.graphs sage.combinat
r"""
Rational polyhedral fans

This module was designed as a part of the framework for toric varieties
(:mod:`~sage.schemes.toric.variety`,
:mod:`~sage.schemes.toric.fano_variety`). While the emphasis is on
complete full-dimensional fans, arbitrary fans are supported. Work
with distinct lattices. The default lattice is :class:`ToricLattice
<sage.geometry.toric_lattice.ToricLatticeFactory>` `N` of the appropriate
dimension. The only case when you must specify lattice explicitly is creation
of a 0-dimensional fan, where dimension of the ambient space cannot be
guessed.

A **rational polyhedral fan** is a *finite* collection of *strictly* convex
rational polyhedral cones, such that the intersection of any two cones of the
fan is a face of each of them and each face of each cone is also a cone of the
fan.

AUTHORS:

- Andrey Novoseltsev (2010-05-15): initial version.

- Andrey Novoseltsev (2010-06-17): substantial improvement during review by
  Volker Braun.

EXAMPLES:

Use :func:`Fan` to construct fans "explicitly"::

    sage: fan = Fan(cones=[(0,1), (1,2)],
    ....:           rays=[(1,0), (0,1), (-1,0)])
    sage: fan
    Rational polyhedral fan in 2-d lattice N

In addition to giving such lists of cones and rays you can also create cones
first using :func:`~sage.geometry.cone.Cone` and then combine them into a fan.
See the documentation of :func:`Fan` for details.

In 2 dimensions there is a unique maximal fan determined by rays, and
you can use :func:`Fan2d` to construct it::

    sage: fan2d = Fan2d(rays=[(1,0), (0,1), (-1,0)])
    sage: fan2d.is_equivalent(fan)
    True

But keep in mind that in higher dimensions the cone data is essential
and cannot be omitted. Instead of building a fan from scratch, for
this tutorial we will use an easy way to get two fans associated to
:class:`lattice polytopes
<sage.geometry.lattice_polytope.LatticePolytopeClass>`:
:func:`FaceFan` and :func:`NormalFan`::

    sage: fan1 = FaceFan(lattice_polytope.cross_polytope(3))
    sage: fan2 = NormalFan(lattice_polytope.cross_polytope(3))

Given such "automatic" fans, you may wonder what are their rays and cones::

    sage: fan1.rays()
    M( 1,  0,  0),
    M( 0,  1,  0),
    M( 0,  0,  1),
    M(-1,  0,  0),
    M( 0, -1,  0),
    M( 0,  0, -1)
    in 3-d lattice M
    sage: fan1.generating_cones()
    (3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M)

The last output is not very illuminating. Let's try to improve it::

    sage: for cone in fan1: print(cone.rays())
    M( 0, 1, 0),
    M( 0, 0, 1),
    M(-1, 0, 0)
    in 3-d lattice M
    M( 0,  0, 1),
    M(-1,  0, 0),
    M( 0, -1, 0)
    in 3-d lattice M
    M(-1,  0,  0),
    M( 0, -1,  0),
    M( 0,  0, -1)
    in 3-d lattice M
    M( 0, 1,  0),
    M(-1, 0,  0),
    M( 0, 0, -1)
    in 3-d lattice M
    M(1, 0,  0),
    M(0, 1,  0),
    M(0, 0, -1)
    in 3-d lattice M
    M(1, 0, 0),
    M(0, 1, 0),
    M(0, 0, 1)
    in 3-d lattice M
    M(1,  0, 0),
    M(0,  0, 1),
    M(0, -1, 0)
    in 3-d lattice M
    M(1,  0,  0),
    M(0, -1,  0),
    M(0,  0, -1)
    in 3-d lattice M

You can also do ::

    sage: for cone in fan1: print(cone.ambient_ray_indices())
    (1, 2, 3)
    (2, 3, 4)
    (3, 4, 5)
    (1, 3, 5)
    (0, 1, 5)
    (0, 1, 2)
    (0, 2, 4)
    (0, 4, 5)

to see indices of rays of the fan corresponding to each cone.

While the above cycles were over "cones in fan", it is obvious that we did not
get ALL the cones: every face of every cone in a fan must also be in the fan,
but all of the above cones were of dimension three. The reason for this
behaviour is that in many cases it is enough to work with generating cones of
the fan, i.e. cones which are not faces of bigger cones. When you do need to
work with lower dimensional cones, you can easily get access to them using
:meth:`~sage.geometry.fan.RationalPolyhedralFan.cones`::

    sage: [cone.ambient_ray_indices() for cone in fan1.cones(2)]
    [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4),
     (2, 4), (3, 4), (0, 5), (1, 5), (3, 5), (4, 5)]

In fact, you do not have to type ``.cones``::

    sage: [cone.ambient_ray_indices() for cone in fan1(2)]
    [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4),
     (2, 4), (3, 4), (0, 5), (1, 5), (3, 5), (4, 5)]

You may also need to know the inclusion relations between all of the cones of
the fan. In this case check out
:meth:`~sage.geometry.fan.RationalPolyhedralFan.cone_lattice`::

    sage: L = fan1.cone_lattice()
    sage: L
    Finite lattice containing 28 elements with distinguished linear extension
    sage: L.bottom()
    0-d cone of Rational polyhedral fan in 3-d lattice M
    sage: L.top()
    Rational polyhedral fan in 3-d lattice M
    sage: cone = L.level_sets()[2][0]
    sage: cone
    2-d cone of Rational polyhedral fan in 3-d lattice M
    sage: sorted(L.hasse_diagram().neighbors(cone))
    [1-d cone of Rational polyhedral fan in 3-d lattice M,
     1-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M,
     3-d cone of Rational polyhedral fan in 3-d lattice M]

You can check how "good" a fan is::

    sage: fan1.is_complete()
    True
    sage: fan1.is_simplicial()
    True
    sage: fan1.is_smooth()
    True

The face fan of the octahedron is really good! Time to remember that we have
also constructed its normal fan::

    sage: fan2.is_complete()
    True
    sage: fan2.is_simplicial()
    False
    sage: fan2.is_smooth()
    False

This one does have some "problems," but we can fix them::

    sage: fan3 = fan2.make_simplicial()
    sage: fan3.is_simplicial()
    True
    sage: fan3.is_smooth()
    False

Note that we had to save the result of
:meth:`~sage.geometry.fan.RationalPolyhedralFan.make_simplicial` in a new fan.
Fans in Sage are immutable, so any operation that does change them constructs
a new fan.

We can also make ``fan3`` smooth, but it will take a bit more work::

    sage: # needs palp
    sage: cube = lattice_polytope.cross_polytope(3).polar()
    sage: sk = cube.skeleton_points(2)
    sage: rays = [cube.point(p) for p in sk]
    sage: fan4 = fan3.subdivide(new_rays=rays)
    sage: fan4.is_smooth()
    True

Let's see how "different" are ``fan2`` and ``fan4``::

    sage: fan2.ngenerating_cones()
    6
    sage: fan2.nrays()
    8
    sage: fan4.ngenerating_cones()                                                      # needs palp
    48
    sage: fan4.nrays()                                                                  # needs palp
    26

Smoothness does not come for free!

Please take a look at the rest of the available functions below and their
complete descriptions. If you need any features that are missing, feel free to
suggest them. (Or implement them on your own and submit a patch to Sage for
inclusion!)
"""

# ****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Callable, Container
from copy import copy
from warnings import warn

import sage.geometry.abc

from sage.structure.richcmp import richcmp_method, richcmp
from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.combination', 'Combinations')
lazy_import('sage.combinat.posets.posets', 'FinitePoset')
from sage.geometry.cone import (_ambient_space_point,
                                Cone,
                                ConvexRationalPolyhedralCone,
                                IntegralRayCollection,
                                normalize_rays)
lazy_import('sage.geometry.hasse_diagram', 'lattice_from_incidences')
from sage.geometry.point_collection import PointCollection
from sage.geometry.toric_lattice import ToricLattice, ToricLattice_generic
lazy_import('sage.geometry.toric_plotter', 'ToricPlotter')
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.timing import walltime
from sage.misc.misc_c import prod
from sage.modules.free_module import span
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


def is_Fan(x) -> bool:
    r"""
    Check if ``x`` is a Fan.

    INPUT:

    - ``x`` -- anything

    OUTPUT: ``True`` if ``x`` is a fan and ``False`` otherwise

    EXAMPLES::

        sage: from sage.geometry.fan import is_Fan
        sage: is_Fan(1)
        doctest:warning...
        DeprecationWarning: The function is_Fan is deprecated; use 'isinstance(..., RationalPolyhedralFan)' instead.
        See https://github.com/sagemath/sage/issues/38126 for details.
        False
        sage: fan = toric_varieties.P2().fan(); fan                                     # needs palp
        Rational polyhedral fan in 2-d lattice N
        sage: is_Fan(fan)                                                               # needs palp
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38126,
                "The function is_Fan is deprecated; "
                "use 'isinstance(..., RationalPolyhedralFan)' instead.")
    return isinstance(x, RationalPolyhedralFan)


def Fan(cones, rays=None, lattice=None, check=True, normalize=True,
        is_complete=None, virtual_rays=None, discard_faces=False,
        allow_arrangement=False):
    r"""
    Construct a rational polyhedral fan.

    .. NOTE::

        Approximate time to construct a fan consisting of `n` cones is `n^2/5`
        seconds. That is half an hour for 100 cones. This time can be
        significantly reduced in the future, but it is still likely to be
        `\sim n^2` (with, say, `/500` instead of `/5`). If you know that your
        input does form a valid fan, use ``check=False`` option to skip
        consistency checks.

    INPUT:

    - ``cones`` -- list of either
      :class:`Cone<sage.geometry.cone.ConvexRationalPolyhedralCone>` objects
      or lists of integers interpreted as indices of generating rays in
      ``rays``. These must be only **maximal** cones of the fan, unless
      ``discard_faces=True`` or ``allow_arrangement=True`` option is specified;

    - ``rays`` -- list of rays given as list or vectors convertible to the
      rational extension of ``lattice``. If ``cones`` are given by
      :class:`Cone<sage.geometry.cone.ConvexRationalPolyhedralCone>` objects
      ``rays`` may be determined automatically. You still may give them
      explicitly to ensure a particular order of rays in the fan. In this case
      you must list all rays that appear in ``cones``. You can give "extra"
      ones if it is convenient (e.g. if you have a big list of rays for
      several fans), but all "extra" rays will be discarded;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically;

    - ``check`` -- by default the input data will be checked for correctness
      (e.g. that intersection of any two given cones is a face of each),
      unless ``allow_arrangement=True`` option is specified. If you
      know for sure that the input is correct, you may significantly decrease
      construction time using ``check=False`` option;

    - ``normalize`` -- you can further speed up construction using
      ``normalize=False`` option. In this case ``cones`` must be a list of
      **sorted** :class:`tuples` and ``rays`` must be immutable primitive
      vectors in ``lattice``. In general, you should not use this option, it
      is designed for code optimization and does not give as drastic
      improvement in speed as the previous one;

    - ``is_complete`` -- every fan can determine on its own if it is complete
      or not, however it can take quite a bit of time for "big" fans with many
      generating cones. On the other hand, in some situations it is known in
      advance that a certain fan is complete. In this case you can pass
      ``is_complete=True`` option to speed up some computations. You may also
      pass ``is_complete=False`` option, although it is less likely to be
      beneficial. Of course, passing a wrong value can compromise the
      integrity of data structures of the fan and lead to wrong results, so
      you should be very careful if you decide to use this option;

    - ``virtual_rays`` -- (optional, computed automatically if needed) a list of
      ray generators to be used for :meth:`virtual_rays`;

    - ``discard_faces`` -- by default, the fan constructor expects the list of
      **maximal** cones, unless ``allow_arrangement=True`` option is specified.
      If you provide "extra" ones and leave ``allow_arrangement=False`` (default)
      and ``check=True`` (default), an exception will be raised.
      If you provide "extra" cones and set ``allow_arrangement=False`` (default)
      and ``check=False``, you may get wrong results as assumptions on internal
      data structures will be invalid. If you want the fan constructor to
      select the maximal cones from the given input, you may provide
      ``discard_faces=True`` option (it works both for ``check=True`` and
      ``check=False``).

    - ``allow_arrangement`` -- by default (``allow_arrangement=False``),
      the fan constructor expects that the intersection of any two given cones is
      a face of each. If ``allow_arrangement=True`` option is specified, then
      construct a rational polyhedralfan from the cone arrangement, so that the
      union of the cones in the polyhedral fan equals to the union of the given
      cones, and each given cone is the union of some cones in the polyhedral fan.

    OUTPUT: a :class:`fan <RationalPolyhedralFan>`

    .. SEEALSO::

        In 2 dimensions you can cyclically order the rays. Hence the
        rays determine a unique maximal fan without having to specify
        the cones, and you can use :func:`Fan2d` to construct this
        fan from just the rays.

    EXAMPLES:

    Let's construct a fan corresponding to the projective plane in several
    ways::

        sage: cone1 = Cone([(1,0), (0,1)])
        sage: cone2 = Cone([(0,1), (-1,-1)])
        sage: cone3 = Cone([(-1,-1), (1,0)])
        sage: P2 = Fan([cone1, cone2, cone2])
        Traceback (most recent call last):
        ...
        ValueError: you have provided 3 cones, but only 2 of them are maximal!
        Use discard_faces=True if you indeed need to construct a fan from
        these cones.

    Oops! There was a typo and ``cone2`` was listed twice as a generating cone
    of the fan. If it was intentional (e.g. the list of cones was generated
    automatically and it is possible that it contains repetitions or faces of
    other cones), use ``discard_faces=True`` option::

        sage: P2 = Fan([cone1, cone2, cone2], discard_faces=True)
        sage: P2.ngenerating_cones()
        2

    However, in this case it was definitely a typo, since the fan of
    `\mathbb{P}^2` has 3 maximal cones::

        sage: P2 = Fan([cone1, cone2, cone3])
        sage: P2.ngenerating_cones()
        3

    Looks better. An alternative way is ::

        sage: rays = [(1,0), (0,1), (-1,-1)]
        sage: cones = [(0,1), (1,2), (2,0)]
        sage: P2a = Fan(cones, rays)
        sage: P2a.ngenerating_cones()
        3
        sage: P2 == P2a
        False

    That may seem wrong, but it is not::

        sage: P2.is_equivalent(P2a)
        True

    See :meth:`~RationalPolyhedralFan.is_equivalent` for details.

    Yet another way to construct this fan is ::

        sage: P2b = Fan(cones, rays, check=False)
        sage: P2b.ngenerating_cones()
        3
        sage: P2a == P2b
        True

    If you try the above examples, you are likely to notice the difference in
    speed, so when you are sure that everything is correct, it is a good idea
    to use ``check=False`` option. On the other hand, it is usually **NOT** a
    good idea to use ``normalize=False`` option::

        sage: P2c = Fan(cones, rays, check=False, normalize=False)
        Traceback (most recent call last):
        ...
        AttributeError: 'tuple' object has no attribute 'parent'...

    Yet another way is to use functions :func:`FaceFan` and :func:`NormalFan`
    to construct fans from :class:`lattice polytopes
    <sage.geometry.lattice_polytope.LatticePolytopeClass>`.

    We have not yet used ``lattice`` argument, since if was determined
    automatically::

        sage: P2.lattice()
        2-d lattice N
        sage: P2b.lattice()
        2-d lattice N

    However, it is necessary to specify it explicitly if you want to construct
    a fan without rays or cones::

        sage: Fan([], [])
        Traceback (most recent call last):
        ...
        ValueError: you must specify the lattice
        when you construct a fan without rays and cones!
        sage: F = Fan([], [], lattice=ToricLattice(2, "L"))
        sage: F
        Rational polyhedral fan in 2-d lattice L
        sage: F.lattice_dim()
        2
        sage: F.dim()
        0

    In the following examples, we test the ``allow_arrangement=True`` option.
    See :issue:`25122`.

    The intersection of the two cones is not a face of each. Therefore,
    they do not belong to the same rational polyhedral fan::

        sage: c1 = Cone([(-2,-1,1), (-2,1,1), (2,1,1), (2,-1,1)])
        sage: c2 = Cone([(-1,-2,1), (-1,2,1), (1,2,1), (1,-2,1)])
        sage: c1.intersection(c2).is_face_of(c1)
        False
        sage: c1.intersection(c2).is_face_of(c2)
        False
        sage: Fan([c1, c2])
        Traceback (most recent call last):
        ...
        ValueError: these cones cannot belong to the same fan!
        ...

    Let's construct the fan using ``allow_arrangement=True`` option::

        sage: fan = Fan([c1, c2], allow_arrangement=True)
        sage: fan.ngenerating_cones()
        5

    Another example where cone c2 is inside cone c1::

        sage: c1 = Cone([(4, 0, 0), (0, 4, 0), (0, 0, 4)])
        sage: c2 = Cone([(2, 1, 1), (1, 2, 1), (1, 1, 2)])
        sage: fan = Fan([c1, c2], allow_arrangement=True)
        sage: fan.ngenerating_cones()
        7
        sage: fan.plot()                                                                # needs sage.plot
        Graphics3d Object

    Cones of different dimension::

        sage: c1 = Cone([(1,0), (0,1)])
        sage: c2 = Cone([(2,1)])
        sage: c3 = Cone([(-1,-2)])
        sage: fan = Fan([c1, c2, c3], allow_arrangement=True)
        sage: for cone in sorted(fan.generating_cones()): print(sorted(cone.rays()))
        [N(-1, -2)]
        [N(0, 1), N(1, 2)]
        [N(1, 0), N(2, 1)]
        [N(1, 2), N(2, 1)]

    A 3-d cone and a 1-d cone::

        sage: c3 = Cone([[0, 1, 1], [1, 0, 1], [0, -1, 1], [-1, 0, 1]])
        sage: c1 = Cone([[0, 0, 1]])
        sage: fan1 = Fan([c1, c3], allow_arrangement=True)
        sage: fan1.plot()                                                               # needs sage.plot
        Graphics3d Object

    A 3-d cone and two 2-d cones::

        sage: c2v = Cone([[0, 1, 1], [0, -1, 1]])
        sage: c2h = Cone([[1, 0, 1], [-1, 0, 1]])
        sage: fan2 = Fan([c2v, c2h, c3], allow_arrangement=True)
        sage: fan2.is_simplicial()
        True
        sage: fan2.is_equivalent(fan1)
        True
    """
    def result():
        # "global" does not work here...
        R, V = rays, virtual_rays
        if V is not None:
            if normalize:
                V = normalize_rays(V, lattice)
            if check:
                R = PointCollection(V, lattice)
                V = PointCollection(V, lattice)
                d = lattice.dimension()
                if len(V) != d - R.dim() or (R + V).dim() != d:
                    raise ValueError("virtual rays must be linearly "
                    "independent and with other rays span the ambient space.")
        return RationalPolyhedralFan(cones, R, lattice, is_complete, V)

    if not check and not normalize and not discard_faces and not allow_arrangement:
        return result()
    if not isinstance(cones, list):
        try:
            cones = list(cones)
        except TypeError:
            raise TypeError(
                "cones must be given as an iterable!"
                "\nGot: %s" % cones)
    if not cones:
        if lattice is None:
            if rays is not None and rays:
                lattice = normalize_rays(rays, lattice)[0].parent()
            else:
                raise ValueError("you must specify the lattice when you "
                                 "construct a fan without rays and cones!")
        cones = ((), )
        rays = ()
        return result()
    if isinstance(cones[0], sage.geometry.abc.ConvexRationalPolyhedralCone):
        # Construct the fan from Cone objects
        if lattice is None:
            lattice = cones[0].lattice()
            # If we determine the lattice automatically, we do not
            # want to force any conversion. TODO: take into account
            # coercions?
            if check:
                for cone in cones:
                    if cone.lattice() != lattice:
                        raise ValueError("cones belong to different lattices "
                            "(%s and %s), cannot determine the lattice of the "
                            "fan!" % (lattice, cone.lattice()))
        for i, cone in enumerate(cones):
            if cone.lattice() != lattice:
                cones[i] = Cone(cone.rays(), lattice, check=False)
        if check:
            for cone in cones:
                if not cone.is_strictly_convex():
                    raise ValueError(
                                    "cones of a fan must be strictly convex!")
        # Optimization for fans generated by a single cone
        if len(cones) == 1 and rays is None:
            cone = cones[0]
            cones = (tuple(range(cone.nrays())), )
            rays = cone.rays()
            is_complete = lattice.dimension() == 0
            return result()
        if allow_arrangement:
            cones = _refine_arrangement_to_fan(cones)
            cones = _discard_faces(cones)
        elif check:
            # Maybe we should compute all faces of all cones and save them for
            # later if we are doing this check?
            generating_cones = []
            for cone in sorted(cones, key=lambda cone: cone.dim(),
                               reverse=True):
                is_generating = True
                for g_cone in generating_cones:
                    i_cone = cone.intersection(g_cone)
                    if i_cone.is_face_of(cone) and i_cone.is_face_of(g_cone):
                        if i_cone.dim() == cone.dim():
                            is_generating = False  # cone is a face of g_cone
                            break
                    else:
                        raise ValueError(
                                "these cones cannot belong to the same fan!"
                                "\nCone 1 rays: %s\nCone 2 rays: %s"
                                % (g_cone.rays(), cone.rays()))
                if is_generating:
                    generating_cones.append(cone)
            if len(cones) > len(generating_cones):
                if discard_faces:
                    cones = generating_cones
                else:
                    raise ValueError("you have provided %d cones, but only %d "
                        "of them are maximal! Use discard_faces=True if you "
                        "indeed need to construct a fan from these cones." %
                        (len(cones), len(generating_cones)))
        elif discard_faces:
            cones = _discard_faces(cones)
        ray_set = set([])
        for cone in cones:
            ray_set.update(cone.rays())
        if rays:    # Preserve the initial order of rays, if they were given
            rays = normalize_rays(rays, lattice)
            new_rays = []
            for ray in rays:
                if ray in ray_set and ray not in new_rays:
                    new_rays.append(ray)
            if len(new_rays) != len(ray_set):
                raise ValueError(
                  "if rays are given, they must include all rays of the fan!")
            rays = new_rays
        else:
            rays = tuple(sorted(ray_set))
        ray_to_index = {ray: i for i, ray in enumerate(rays)}
        cones = (tuple(sorted(ray_to_index[ray] for ray in cone.rays()))
                 for cone in cones)
        return result()
    # Construct the fan from rays and "tuple cones"
    rays = normalize_rays(rays, lattice)
    for n, cone in enumerate(cones):
        try:
            cones[n] = sorted(cone)
        except TypeError:
            raise TypeError("cannot interpret %s as a cone!" % cone)
    if not check and not discard_faces and not allow_arrangement:
        return result()
    # If we do need to make all the check, build explicit cone objects first
    return Fan((Cone([rays[n] for n in cone], lattice) for cone in cones),
               rays, lattice, is_complete=is_complete,
               virtual_rays=virtual_rays, discard_faces=discard_faces,
               allow_arrangement=allow_arrangement)


def FaceFan(polytope, lattice=None):
    r"""
    Construct the face fan of the given rational ``polytope``.

    INPUT:

    - ``polytope`` -- a :func:`polytope
      <sage.geometry.polyhedron.constructor.Polyhedron>` over `\QQ` or
      a :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`. A (not
      necessarily full-dimensional) polytope containing the origin in
      its :meth:`relative interior
      <sage.geometry.polyhedron.base.Polyhedron_base.relative_interior_contains>`.

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically.

    OUTPUT: :class:`rational polyhedral fan <RationalPolyhedralFan>`

    See also :func:`NormalFan`.

    EXAMPLES:

    Let's construct the fan corresponding to the product of two projective
    lines::

        sage: diamond = lattice_polytope.cross_polytope(2)
        sage: P1xP1 = FaceFan(diamond)
        sage: P1xP1.rays()
        M( 1,  0),
        M( 0,  1),
        M(-1,  0),
        M( 0, -1)
        in 2-d lattice M
        sage: for cone in P1xP1: print(cone.rays())
        M(-1,  0),
        M( 0, -1)
        in 2-d lattice M
        M( 0, 1),
        M(-1, 0)
        in 2-d lattice M
        M(1, 0),
        M(0, 1)
        in 2-d lattice M
        M(1,  0),
        M(0, -1)
        in 2-d lattice M

    TESTS::

        sage: cuboctahed = polytopes.cuboctahedron()
        sage: FaceFan(cuboctahed)
        Rational polyhedral fan in 3-d lattice M
        sage: cuboctahed.is_lattice_polytope(), cuboctahed.dilation(1/2).is_lattice_polytope()
        (True, False)
        sage: fan1 = FaceFan(cuboctahed)
        sage: fan2 = FaceFan(cuboctahed.dilation(2).lattice_polytope())
        sage: fan1.is_equivalent(fan2)
        True

        sage: ray = Polyhedron(vertices=[(-1,-1)], rays=[(1,1)])
        sage: FaceFan(ray)
        Traceback (most recent call last):
        ...
        ValueError: face fans are defined only for
        polytopes containing the origin as an interior point!

        sage: interval_in_QQ2 = Polyhedron([(0,-1), (0,+1)])
        sage: FaceFan(interval_in_QQ2).generating_cones()
        (1-d cone of Rational polyhedral fan in 2-d lattice M,
         1-d cone of Rational polyhedral fan in 2-d lattice M)

        sage: FaceFan(Polyhedron([(-1,0), (1,0), (0,1)])) # origin on facet
        Traceback (most recent call last):
        ...
        ValueError: face fans are defined only for
        polytopes containing the origin as an interior point!
    """
    interior_point_error = ValueError(
        "face fans are defined only for polytopes containing "
        "the origin as an interior point!")
    if isinstance(polytope, sage.geometry.abc.LatticePolytope):
        if any(d <= 0 for d in polytope.distances([0] * polytope.dim())):
            raise interior_point_error
        cones = (f.ambient_vertex_indices() for f in polytope.facets())
        rays = polytope.vertices()
        is_complete = polytope.dim() == polytope.lattice_dim()
    else:
        origin = polytope.ambient_space().zero()
        if not (polytope.is_compact() and
                polytope.relative_interior_contains(origin)):
            raise interior_point_error
        cones = [[v.index() for v in facet.incident()]
                  for facet in polytope.inequalities()]
        rays = [vector(_) for _ in polytope.vertices()]
        is_complete = polytope.dim() == polytope.ambient_dim()
        if lattice is None:
            # Since default lattice polytopes are in the M lattice,
            # treat polyhedra as being there as well.
            lattice = ToricLattice(len(origin)).dual()
    return Fan(cones, rays, lattice=lattice, check=False,
               is_complete=is_complete)


def NormalFan(polytope, lattice=None):
    r"""
    Construct the normal fan of the given rational ``polytope``.

    This returns the inner normal fan. For the outer normal fan, use
    ``NormalFan(-P)``.

    INPUT:

    - ``polytope`` -- a full-dimensional :func:`polytope
      <sage.geometry.polyhedron.constructor.Polyhedron>` over `\QQ`
      or:class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`.

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically.

    OUTPUT: :class:`rational polyhedral fan <RationalPolyhedralFan>`

    See also :func:`FaceFan`.

    EXAMPLES:

    Let's construct the fan corresponding to the product of two projective
    lines::

        sage: square = LatticePolytope([(1,1), (-1,1), (-1,-1), (1,-1)])
        sage: P1xP1 = NormalFan(square)
        sage: P1xP1.rays()
        N( 1,  0),
        N( 0,  1),
        N(-1,  0),
        N( 0, -1)
        in 2-d lattice N
        sage: for cone in P1xP1: print(cone.rays())
        N(-1,  0),
        N( 0, -1)
        in 2-d lattice N
        N(1,  0),
        N(0, -1)
        in 2-d lattice N
        N(1, 0),
        N(0, 1)
        in 2-d lattice N
        N( 0, 1),
        N(-1, 0)
        in 2-d lattice N

        sage: cuboctahed = polytopes.cuboctahedron()
        sage: NormalFan(cuboctahed)
        Rational polyhedral fan in 3-d lattice N

    TESTS::

        sage: cuboctahed.is_lattice_polytope(), cuboctahed.dilation(1/2).is_lattice_polytope()
        (True, False)
        sage: fan1 = NormalFan(cuboctahed)
        sage: fan2 = NormalFan(cuboctahed.dilation(2).lattice_polytope())
        sage: fan1.is_equivalent(fan2)
        True
    """
    dimension_error = ValueError(
        'the normal fan is only defined for full-dimensional polytopes')
    if isinstance(polytope, sage.geometry.abc.LatticePolytope):
        if polytope.dim() != polytope.lattice_dim():
            raise dimension_error
        rays = polytope.facet_normals()
        cones = (v.ambient_facet_indices() for v in polytope.faces(dim=0))
    else:
        if polytope.dim() != polytope.ambient_dim():
            raise dimension_error
        if not polytope.is_compact():
            raise NotImplementedError('the normal fan is only supported for polytopes (compact polyhedra).')
        cones = [[ieq.index() for ieq in vertex.incident()]
                 for vertex in polytope.vertices()]
        rays = [ieq.A() for ieq in polytope.inequalities()]
    return Fan(cones, rays, lattice=lattice, check=False, is_complete=True)


def Fan2d(rays, lattice=None):
    r"""
    Construct the maximal 2-d fan with given ``rays``.

    In two dimensions we can uniquely construct a fan from just rays,
    just by cyclically ordering the rays and constructing as many
    cones as possible. This is why we implement a special constructor
    for this case.

    INPUT:

    - ``rays`` -- list of rays given as list or vectors convertible to
      the rational extension of ``lattice``. Duplicate rays are
      removed without changing the ordering of the remaining rays.

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically.

    EXAMPLES::

        sage: Fan2d([(0,1), (1,0)])
        Rational polyhedral fan in 2-d lattice N
        sage: Fan2d([], lattice=ToricLattice(2, 'myN'))
        Rational polyhedral fan in 2-d lattice myN

    The ray order is as specified, even if it is not the cyclic order::

        sage: fan1 = Fan2d([(0,1), (1,0)])
        sage: fan1.rays()
        N(0, 1),
        N(1, 0)
        in 2-d lattice N

        sage: fan2 = Fan2d([(1,0), (0,1)])
        sage: fan2.rays()
        N(1, 0),
        N(0, 1)
        in 2-d lattice N

        sage: fan1 == fan2, fan1.is_equivalent(fan2)
        (False, True)

        sage: fan = Fan2d([(1,1), (-1,-1), (1,-1), (-1,1)])
        sage: [cone.ambient_ray_indices() for cone in fan]
        [(2, 1), (1, 3), (3, 0), (0, 2)]
        sage: fan.is_complete()
        True

    TESTS::

        sage: Fan2d([(0,1), (0,1)]).generating_cones()
        (1-d cone of Rational polyhedral fan in 2-d lattice N,)

        sage: Fan2d([(1,1), (-1,-1)]).generating_cones()
        (1-d cone of Rational polyhedral fan in 2-d lattice N,
         1-d cone of Rational polyhedral fan in 2-d lattice N)

        sage: Fan2d([])
        Traceback (most recent call last):
        ...
        ValueError: you must specify a 2-dimensional lattice
        when you construct a fan without rays.

        sage: Fan2d([(3,4)]).rays()
        N(3, 4)
        in 2-d lattice N

        sage: Fan2d([(0,1,0)])
        Traceback (most recent call last):
        ...
        ValueError: the lattice must be 2-dimensional.

        sage: Fan2d([(0,1), (1,0), (0,0)])
        Traceback (most recent call last):
        ...
        ValueError: only nonzero vectors define rays

        sage: Fan2d([(0, -2), (2, -10), (1, -3), (2, -9), (2, -12), (1, 1),
        ....:        (2, 1), (1, -5), (0, -6), (1, -7), (0, 1), (2, -4),
        ....:        (2, -2), (1, -9), (1, -8), (2, -6), (0, -1), (0, -3),
        ....:        (2, -11), (2, -8), (1, 0), (0, -5), (1, -4), (2, 0),
        ....:        (1, -6), (2, -7), (2, -5), (-1, -3), (1, -1), (1, -2),
        ....:        (0, -4), (2, -3), (2, -1)]).cone_lattice()
        Finite lattice containing 44 elements with distinguished linear extension

        sage: Fan2d([(1,1)]).is_complete()
        False
        sage: Fan2d([(1,1), (-1,-1)]).is_complete()
        False
        sage: Fan2d([(1,0), (0,1)]).is_complete()
        False
    """
    if not rays:
        if lattice is None or lattice.dimension() != 2:
            raise ValueError('you must specify a 2-dimensional lattice when '
                             'you construct a fan without rays.')
        return RationalPolyhedralFan(cones=((), ), rays=(), lattice=lattice)

    # remove multiple rays without changing order
    rays = normalize_rays(rays, lattice)
    rays = sorted((r, i) for i, r in enumerate(rays))
    distinct_rays = [rays[i] for i in range(len(rays))
                     if rays[i][0] != rays[i-1][0]]
    if distinct_rays:
        rays = sorted((i, r) for r, i in distinct_rays)
        rays = [r[1] for r in rays]
    else:  # all given rays were the same
        rays = [rays[0][0]]
    lattice = rays[0].parent()
    if lattice.dimension() != 2:
        raise ValueError('the lattice must be 2-dimensional.')
    n = len(rays)
    if n == 1 or n == 2 and rays[0] == -rays[1]:
        cones = [(i, ) for i in range(n)]
        return RationalPolyhedralFan(cones, rays, lattice, False)

    import math
    # each sorted_rays entry = (angle, ray, original_ray_index)
    sorted_rays = sorted((math.atan2(r[0], r[1]), r, i)
                         for i, r in enumerate(rays))
    cones = []
    is_complete = True
    for i in range(n):
        r0 = sorted_rays[i-1][1]
        r1 = sorted_rays[i][1]
        if r1.is_zero():
            raise ValueError('only nonzero vectors define rays')
        assert r0 != r1
        cross_prod = r0[0] * r1[1] - r0[1] * r1[0]
        if cross_prod < 0:
            r0_index = (i-1) % len(sorted_rays)
            r1_index = i
            cones.append((sorted_rays[r0_index][2], sorted_rays[r1_index][2]))
        else:
            is_complete = False
    return RationalPolyhedralFan(cones, rays, lattice, is_complete)


class Cone_of_fan(ConvexRationalPolyhedralCone):
    r"""
    Construct a cone belonging to a fan.

    .. WARNING::

        This class does not check that the input defines a valid cone of a
        fan. You must not construct objects of this class directly.

    In addition to all of the properties of "regular" :class:`cones
    <sage.geometry.cone.ConvexRationalPolyhedralCone>`, such cones know their
    relation to the fan.

    INPUT:

    - ``ambient`` -- fan whose cone is constructed

    - ``ambient_ray_indices`` -- increasing list or tuple of integers, indices
      of rays of ``ambient`` generating this cone

    OUTPUT: cone of ``ambient``

    EXAMPLES:

    The intended way to get objects of this class is the following::

        sage: # needs palp
        sage: fan = toric_varieties.P1xP1().fan()
        sage: cone = fan.generating_cone(0); cone
        2-d cone of Rational polyhedral fan in 2-d lattice N
        sage: cone.ambient_ray_indices()
        (0, 2)
        sage: cone.star_generator_indices()
        (0,)
    """

    def __init__(self, ambient, ambient_ray_indices):
        r"""
        See :class:`Cone_of_Fan` for documentation.

        TESTS:

        The following code is likely to construct an invalid object, we just
        test that creation of cones of fans is working::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: cone = sage.geometry.fan.Cone_of_fan(fan, (0,)); cone                 # needs palp
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: TestSuite(cone).run()                                                 # needs palp
        """
        super().__init__(ambient=ambient,
                         ambient_ray_indices=ambient_ray_indices)
        self._is_strictly_convex = True
        # Because if not, this cone should not have been constructed

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: # needs palp
            sage: P1xP1 = toric_varieties.P1xP1()
            sage: cone = P1xP1.fan().generating_cone(0)
            sage: cone._repr_()
            '2-d cone of Rational polyhedral fan in 2-d lattice N'
            sage: cone.facets()[0]._repr_()
            '1-d cone of Rational polyhedral fan in 2-d lattice N'
        """
        # The base class would print "face of" instead of  "cone of"
        return "%d-d cone of %s" % (self.dim(), self.ambient())

    def star_generator_indices(self):
        r"""
        Return indices of generating cones of the "ambient fan" containing
        ``self``.

        OUTPUT: increasing :class:`tuple` of integers

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()                                       # needs palp
            sage: cone = P1xP1.fan().generating_cone(0)                                 # needs palp
            sage: cone.star_generator_indices()                                         # needs palp
            (0,)

        TESTS:

        A mistake in this function used to cause the problem reported in
        :issue:`9782`. We check that now everything is working smoothly::

            sage: f = Fan([(0, 2, 4),
            ....:          (0, 4, 5),
            ....:          (0, 3, 5),
            ....:          (0, 1, 3),
            ....:          (0, 1, 2),
            ....:          (2, 4, 6),
            ....:          (4, 5, 6),
            ....:          (3, 5, 6),
            ....:          (1, 3, 6),
            ....:          (1, 2, 6)],
            ....:         [(-1, 0, 0),
            ....:          (0, -1, 0),
            ....:          (0, 0, -1),
            ....:          (0, 0, 1),
            ....:          (0, 1, 2),
            ....:          (0, 1, 3),
            ....:          (1, 0, 4)])
            sage: f.is_complete()
            True
            sage: X = ToricVariety(f)
            sage: X.fan().is_complete()
            True
        """
        if "_star_generator_indices" not in self.__dict__:
            fan = self.ambient()
            sgi = set(range(fan.ngenerating_cones()))
            for ray in self.ambient_ray_indices():
                sgi.intersection_update(fan._ray_to_cones(ray))
            self._star_generator_indices = tuple(sorted(sgi))
        return self._star_generator_indices

    def star_generators(self):
        r"""
        Return indices of generating cones of the "ambient fan" containing
        ``self``.

        OUTPUT: increasing :class:`tuple` of integers

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()                                       # needs palp
            sage: cone = P1xP1.fan().generating_cone(0)                                 # needs palp
            sage: cone.star_generators()                                                # needs palp
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)
        """
        if "_star_generators" not in self.__dict__:
            self._star_generators = tuple(self.ambient().generating_cone(i)
                                    for i in self.star_generator_indices())
        return self._star_generators


@richcmp_method
class RationalPolyhedralFan(IntegralRayCollection, Callable, Container):
    r"""
    Create a rational polyhedral fan.

    .. WARNING::

        This class does not perform any checks of correctness of input nor
        does it convert input into the standard representation. Use
        :func:`Fan` to construct fans from "raw data" or :func:`FaceFan` and
        :func:`NormalFan` to get fans associated to polytopes.

    Fans are immutable, but they cache most of the returned values.

    INPUT:

    - ``cones`` -- list of generating cones of the fan, each cone given as a
      list of indices of its generating rays in ``rays``;

    - ``rays`` -- list of immutable primitive vectors in ``lattice``
      consisting of exactly the rays of the fan (i.e. no "extra" ones);

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If ``None``, it will be determined
      as :func:`parent` of the first ray. Of course, this cannot be done if
      there are no rays, so in this case you must give an appropriate
      ``lattice`` directly;

    - ``is_complete`` -- if given, must be ``True`` or ``False`` depending on
      whether this fan is complete or not. By default, it will be determined
      automatically if necessary;

    - ``virtual_rays`` -- if given, must be a list of immutable primitive
      vectors in ``lattice``, see :meth:`virtual_rays` for details. By default,
      it will be determined automatically if necessary.

    OUTPUT:

    rational polyhedral fan
    """

    def __init__(self, cones, rays, lattice,
                 is_complete=None, virtual_rays=None):
        r"""
        See :class:`RationalPolyhedralFan` for documentation.

        TESTS::

            sage: v = vector([0,1])
            sage: v.set_immutable()
            sage: f = sage.geometry.fan.RationalPolyhedralFan(
            ....:                       [(0,)], [v], None)
            sage: f.rays()
            (0, 1)
            in Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: TestSuite(f).run()
            sage: f = Fan([(0,)], [(0,1)])
            sage: TestSuite(f).run()
        """
        super().__init__(rays, lattice)
        self._generating_cones = tuple(Cone_of_fan(self, c) for c in cones)
        for i, cone in enumerate(self._generating_cones):
            cone._star_generator_indices = (i,)
        # Knowing completeness drastically affects the speed of cone lattice
        # computation and containment check, so we have a special way to
        # optimize it.
        if is_complete is not None:
            self._is_complete = is_complete
        # Computing virtual rays is fast, but it may be convenient to choose
        # them based on relation to other cones and fans.
        if virtual_rays is not None:
            self._virtual_rays = PointCollection(virtual_rays, self.lattice())

    def _sage_input_(self, sib, coerced):
        """
        Return Sage command to reconstruct ``self``.

        See :mod:`sage.misc.sage_input` for details.

        EXAMPLES::

            sage: fan = Fan([Cone([(1,0), (1,1)]), Cone([(-1,-1)])])
            sage: sage_input(fan)
            Fan(cones=[[1, 2], [0]], rays=[(-1, -1), (1, 0), (1, 1)])
        """
        cones = [[ZZ(_) for _ in c.ambient_ray_indices()] for c in self.generating_cones()]
        rays = [sib(tuple(r)) for r in self.rays()]
        return sib.name('Fan')(cones=cones, rays=rays)

    def __call__(self, dim=None, codim=None):
        r"""
        Return the specified cones of ``self``.

        .. NOTE::

            "Direct call" syntax is a synonym for :meth:`cones` method except
            that in the case of no input parameters this function returns
            just ``self``.

        INPUT:

        - ``dim`` -- dimension of the requested cones

        - ``codim`` -- codimension of the requested cones

        OUTPUT:

        cones of ``self`` of the specified (co)dimension if it was given,
        otherwise ``self``

        TESTS::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan(1)
            (1-d cone of Rational polyhedral fan in 2-d lattice N,
             1-d cone of Rational polyhedral fan in 2-d lattice N,
             1-d cone of Rational polyhedral fan in 2-d lattice N)
            sage: fan(2)
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(dim=2)
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(codim=2)
            (0-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(dim=1, codim=1)
            Traceback (most recent call last):
            ...
            ValueError: dimension and codimension
            cannot be specified together!
            sage: fan() is fan
            True
        """
        if dim is None and codim is None:
            # "self.cones()" returns all cones, but for the call syntax
            # "self()" we return just "self", which seems to be more natural
            # and convenient for ToricVariety.fan() method.
            return self
        else:
            return self.cones(dim, codim)

    def __richcmp__(self, right, op):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything

        OUTPUT: boolean

        There is equality if ``right`` is also a fan, their rays are
        the same and stored in the same order, and their generating
        cones are the same and stored in the same order.

        TESTS::

            sage: f1 = Fan(cones=[(0,1), (1,2)],
            ....:          rays=[(1,0), (0,1), (-1, 0)],
            ....:          check=False)
            sage: f2 = Fan(cones=[(1,2), (0,1)],
            ....:          rays=[(1,0), (0,1), (-1, 0)],
            ....:          check=False)
            sage: f3 = Fan(cones=[(1,2), (0,1)],
            ....:          rays=[(1,0), (0,1), (-1, 0)],
            ....:          check=False)
            sage: f1 > f2
            True
            sage: f2 < f1
            True
            sage: f2 == f3
            True
            sage: f2 is f3
            False
        """
        if isinstance(right, RationalPolyhedralFan):
            return richcmp([self.rays(), self.virtual_rays(),
                            self.generating_cones()],
                           [right.rays(), right.virtual_rays(),
                            right.generating_cones()], op)
        else:
            return NotImplemented

    def __contains__(self, cone):
        r"""
        Check if ``cone`` is equivalent to a cone of the fan.

        See :meth:`_contains` (which is called by this function) for
        documentation.

        TESTS::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])
            sage: f.generating_cone(0) in f
            True
            sage: cone1 in f
            True
            sage: (1,1) in f    # not a cone
            False
            sage: "Ceci n'est pas un cone" in f
            False
        """
        return self._contains(cone)

    def __iter__(self):
        r"""
        Return an iterator over generating cones of ``self``.

        OUTPUT: iterator

        TESTS::

            sage: f = Fan(cones=[(0,1), (1,2)],
            ....:         rays=[(1,0), (0,1), (-1, 0)],
            ....:         check=False)
            sage: for cone in f: print(cone.rays())
            N(1, 0),
            N(0, 1)
            in 2-d lattice N
            N( 0, 1),
            N(-1, 0)
            in 2-d lattice N
        """
        return iter(self.generating_cones())

    def _compute_cone_lattice(self):
        r"""
        Compute the cone lattice of ``self``.

        See :meth:`cone_lattice` for documentation.

        TESTS:

        We use different algorithms depending on available information. One of
        the common cases is a fan which is KNOWN to be complete, i.e. we do
        not even need to check if it is complete::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.cone_lattice()  # indirect doctest                                # needs palp
            Finite lattice containing 10 elements with distinguished linear extension

        These 10 elements are: 1 origin, 4 rays, 4 generating cones, 1 fan.

        Another common case is the fan of faces of a single cone::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: fan = Fan([quadrant])
            sage: fan.cone_lattice() # indirect doctest
            Finite poset containing 5 elements with distinguished linear extension

        These 5 elements are: 1 origin, 2 rays, 1 generating cone, 1 fan.

        A subcase of this common case is treatment of fans consisting of the
        origin only, which used to be handled incorrectly :issue:`18613`::

            sage: fan = Fan([Cone([], ToricLattice(0))])
            sage: list(fan.cone_lattice())
            [0-d cone of Rational polyhedral fan in 0-d lattice N,
             Rational polyhedral fan in 0-d lattice N]
            sage: fan = Fan([Cone([], ToricLattice(1))])
            sage: list(fan.cone_lattice())
            [0-d cone of Rational polyhedral fan in 1-d lattice N,
             Rational polyhedral fan in 1-d lattice N]

        Finally, we have "intermediate" fans which are incomplete but are
        generated by more than one cone::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.rays()
            N(-1, 0),
            N( 0, 1),
            N( 1, 0)
            in 2-d lattice N
            sage: for cone in fan: print(cone.ambient_ray_indices())
            (1, 2)
            (0,)
            sage: L = fan.cone_lattice() # indirect doctest
            sage: L
            Finite poset containing 6 elements with distinguished linear extension

        Here we got 1 origin, 3 rays (one is a generating cone),
        1 2-dimensional cone (a generating one), and 1 fan.
        """
        # Define a face constructor
        def FanFace(rays, cones):
            if not cones:       # The top face, fan itself
                return self
            if len(cones) == 1:  # MAY be a generating cone or NOT!!!
                g_cone = self.generating_cone(cones[0])
                if g_cone.ambient_ray_indices() == rays:
                    return g_cone
            face = Cone_of_fan(ambient=self, ambient_ray_indices=rays)
            face._star_generator_indices = cones
            return face
        # Check directly if we know completeness already, since *determining*
        # completeness relies on this function
        if "_is_complete" in self.__dict__ and self._is_complete:
            # We can use a fast way for complete fans
            self._cone_lattice = lattice_from_incidences(
                                # When there are no rays, fan is the only atom
                                self._ray_to_cones() if self.rays() else [()],
                                (cone.ambient_ray_indices() for cone in self),
                                FanFace, key=id(self))
        else:
            # For general fans we will "merge" face lattices of generating
            # cones.
            from sage.graphs.digraph import DiGraph
            L = DiGraph()
            face_to_rays = {}  # face |---> (indices of fan rays)
            rays_to_index = {}  # (indices of fan rays) |---> face index
            # face index |---> (indices of containing generating cones)
            index_to_cones = []
            # During construction index 0 will correspond to the fan
            # We think of the fan not being in the cone even when there is
            # only one cone
            index_to_cones.append(())
            next_index = 1
            for i, cone in enumerate(self):
                # Set up translation of faces of cone to rays and indices
                # We make a standalone cone to compute its standalone face
                # lattice, since cones of fans get their lattices from fans
                L_cone = Cone(cone.rays(), lattice=self.lattice(),
                              check=False, normalize=False).face_lattice()
                for f in L_cone:
                    f_rays = tuple(cone.ambient_ray_indices()[ray]
                                   for ray in f.ambient_ray_indices())
                    face_to_rays[f] = f_rays
                    try:
                        f_index = rays_to_index[f_rays]
                        index_to_cones[f_index].append(i)
                    except KeyError:        # Did not see f before
                        f_index = next_index
                        next_index += 1
                        rays_to_index[f_rays] = f_index
                        index_to_cones.append([i])
                # Add all relations between faces of cone to L
                for f, g in L_cone.cover_relations_iterator():
                    L.add_edge(rays_to_index[face_to_rays[f]],
                               rays_to_index[face_to_rays[g]])
                # Add the inclusion of cone into the fan itself
                L.add_edge(
                        rays_to_index[face_to_rays[L_cone.top()]], 0)

            # Enumeration of graph vertices must be a linear extension of the
            # poset
            new_order = L.topological_sort()
            # Make sure that generating cones are in the end in proper order
            tail = [rays_to_index[gc.ambient_ray_indices()] for gc in self]
            tail.append(0)  # We know that the fan itself has index 0
            new_order = [n for n in new_order if n not in tail] + tail
            # Make sure that rays are in the beginning in proper order
            head = [rays_to_index[()]]  # Empty face
            head.extend(rays_to_index[(n,)] for n in range(self.nrays()))
            new_order = head + [n for n in new_order if n not in head]
            # "Invert" this list to a dictionary
            labels = {old: new for new, old in enumerate(new_order)}
            L.relabel(labels)

            elements = [None] * next_index
            for rays, index in rays_to_index.items():
                elements[labels[index]] = FanFace(
                                           rays, tuple(index_to_cones[index]))
            # We need "special treatment" for the whole fan. If we added its
            # ray incidence information to the total list, it would be
            # confused with the generating cone in the case of a single cone.
            elements[labels[0]] = FanFace(tuple(range(self.nrays())), ())
            D = dict(enumerate(elements))
            L.relabel(D)
            self._cone_lattice = FinitePoset(L, elements, key=id(self))

    def _contains(self, cone) -> bool:
        r"""
        Check if ``cone`` is equivalent to a cone of the fan.

        This function is called by :meth:`__contains__` and :meth:`contains`
        to ensure the same call depth for warning messages.

        INPUT:

        - ``cone`` -- anything

        OUTPUT:

        ``False`` if ``cone`` is not a cone or if ``cone`` is not
        equivalent to a cone of the fan, ``True`` otherwise

        TESTS::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])
            sage: f._contains(cone1)
            True
            sage: f._contains((1,1))  # this is not a cone!
            False

        Note that the ambient fan of the cone does not matter::

            sage: cone1_f = f.generating_cone(0)
            sage: cone1_f is cone1
            False
            sage: cone1_f.is_equivalent(cone1)
            True
            sage: cone1 in Fan([cone1, cone2])  # not a cone of any particular fan
            True
            sage: cone1_f in Fan([cone1, cone2])  # belongs to different fan, but equivalent cone
            True
        """
        try:
            self.embed(cone)    # Fails if cone is not in self.
            return True
        except TypeError:   # cone is not a cone
            return False
        except ValueError:  # cone is a cone, but wrong
            if not cone.lattice().is_submodule(self.lattice()):
                warn("you have checked if a fan contains a cone "
                     "from another lattice, this is always False!",
                     stacklevel=3)
            return False

    def support_contains(self, *args):
        r"""
        Check if a point is contained in the support of the fan.

        The support of a fan is the union of all cones of the fan. If
        you want to know whether the fan contains a given cone, you
        should use :meth:`contains` instead.

        INPUT:

        - ``*args`` -- an element of ``self.lattice()`` or something
          that can be converted to it (for example, a list of
          coordinates).

        OUTPUT:

        ``True`` if ``point`` is contained in the support of the
        fan, ``False`` otherwise

        TESTS::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])

        We check if some points are in this fan::

            sage: f.support_contains(f.lattice()(1,0))
            True
            sage: f.support_contains(cone1)    # a cone is not a point of the lattice
            False
            sage: f.support_contains((1,0))
            True
            sage: f.support_contains(1,1)
            True
            sage: f.support_contains((-1,0))
            False
            sage: f.support_contains(f.lattice().dual()(1,0))  # random output (warning)
            False
            sage: f.support_contains(f.lattice().dual()(1,0))
            False
            sage: f.support_contains(1)
            False
            sage: f.support_contains(0)   # 0 converts to the origin in the lattice
            True
            sage: f.support_contains(1/2, sqrt(3))                                      # needs sage.symbolic
            True
            sage: f.support_contains(-1/2, sqrt(3))                                     # needs sage.symbolic
            False
        """
        if len(args) == 1:
            point = args[0]
        else:
            point = args

        try:
            point = _ambient_space_point(self, point)
        except TypeError as ex:
            if str(ex).endswith("have incompatible lattices!"):
                warn("you have checked if a fan contains a point "
                     "from an incompatible lattice, this is False!",
                     stacklevel=3)
            return False
        if self.is_complete():
            return True
        return any(point in cone for cone in self)

    def cartesian_product(self, other, lattice=None):
        r"""
        Return the Cartesian product of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a :class:`rational polyhedral fan
          <sage.geometry.fan.RationalPolyhedralFan>`;

        - ``lattice`` -- (optional) the ambient lattice for the
          Cartesian product fan. By default, the direct sum of the
          ambient lattices of ``self`` and ``other`` is constructed.

        OUTPUT:

        a :class:`fan <RationalPolyhedralFan>` whose cones are all pairwise
        Cartesian products of the cones of ``self`` and ``other``

        EXAMPLES::

            sage: K = ToricLattice(1, 'K')
            sage: fan1 = Fan([[0],[1]], [(1,),(-1,)], lattice=K)
            sage: L = ToricLattice(2, 'L')
            sage: fan2 = Fan(rays=[(1,0), (0,1), (-1,-1)],
            ....:            cones=[[0,1], [1,2], [2,0]], lattice=L)
            sage: fan1.cartesian_product(fan2)
            Rational polyhedral fan in 3-d lattice K+L
            sage: _.ngenerating_cones()
            6
        """
        assert isinstance(other, RationalPolyhedralFan)
        rc = super().cartesian_product(other, lattice)
        self_cones = [cone.ambient_ray_indices() for cone in self]
        n = self.nrays()
        other_cones = [tuple(n + i for i in cone.ambient_ray_indices())
                       for cone in other]
        new_cones = [c1 + c2 for c1 in self_cones for c2 in other_cones]
        try:    # Is completeness of the result obvious?
            return RationalPolyhedralFan(new_cones, rc.rays(), rc.lattice(),
                                    self._is_complete and other._is_complete)
        except AttributeError:  # The result is either incomplete or unknown.
            return RationalPolyhedralFan(new_cones, rc.rays(), rc.lattice())

    def __neg__(self):
        """
        Return the fan where each cone is replaced by the opposite cone.

        EXAMPLES::

            sage: c0 = Cone([(1,1),(0,1)])
            sage: c1 = Cone([(1,1),(1,0)])
            sage: F = Fan([c0, c1]); F
            Rational polyhedral fan in 2-d lattice N
            sage: G = -F; G  # indirect doctest
            Rational polyhedral fan in 2-d lattice N
            sage: -G==F
            True
            sage: G.rays()
            N( 0, -1),
            N(-1,  0),
            N(-1, -1)
            in 2-d lattice N
        """
        new_rays = [-r1 for r1 in self.rays()]
        for r in new_rays:
            r.set_immutable()
        self_cones = [cone.ambient_ray_indices() for cone in self]
        return RationalPolyhedralFan(self_cones, new_rays, self.lattice())

    def common_refinement(self, other):
        """
        Return the common refinement of this fan and ``other``.

        INPUT:

        - ``other`` -- a :class:`fan <RationalPolyhedralFan>` in the same
          :meth:`lattice` and with the same support as this fan

        OUTPUT: a :class:`fan <RationalPolyhedralFan>`

        EXAMPLES:

        Refining a fan with itself gives itself::

            sage: F0 = Fan2d([(1,0), (0,1), (-1,0), (0,-1)])
            sage: F0.common_refinement(F0) == F0
            True

        A more complex example with complete fans::

            sage: F1 = Fan([[0],[1]], [(1,),(-1,)])
            sage: F2 = Fan2d([(1,0), (1,1), (0,1), (-1,0), (0,-1)])
            sage: F3 = F2.cartesian_product(F1)
            sage: F4 = F1.cartesian_product(F2)
            sage: FF = F3.common_refinement(F4)
            sage: F3.ngenerating_cones()
            10
            sage: F4.ngenerating_cones()
            10
            sage: FF.ngenerating_cones()
            13

        An example with two non-complete fans with the same support::

            sage: F5 = Fan2d([(1,0), (1,2), (0,1)])
            sage: F6 = Fan2d([(1,0), (2,1), (0,1)])
            sage: F5.common_refinement(F6).ngenerating_cones()
            3

        Both fans must live in the same lattice::

            sage: F0.common_refinement(F1)
            Traceback (most recent call last):
            ...
            ValueError: the fans are not in the same lattice
        """
        from sage.categories.homset import End
        from sage.geometry.fan_morphism import FanMorphism
        N = self.lattice()
        if other.lattice() is not N:
            raise ValueError('the fans are not in the same lattice')
        id = End(N).identity()
        subdivision = FanMorphism(id, self, other, subdivide=True).domain_fan()
        if not self.is_complete():
            # Construct the opposite morphism to ensure support equality
            FanMorphism(id, other, self, subdivide=True)
        return subdivision

    def _latex_(self) -> str:
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: f = Fan(cones=[(0,1), (1,2)],
            ....:         rays=[(1,0), (0,1), (-1, 0)],
            ....:         check=False)
            sage: f._latex_()
            '\\Sigma^{2}'
        """
        return r"\Sigma^{%s}" % self.lattice_dim()

    def _ray_to_cones(self, i=None):
        r"""
        Return the set of generating cones containing the ``i``-th ray.

        INPUT:

        - ``i`` -- integer; index of a ray of ``self``

        OUTPUT:

        :class:`frozenset` of indices of generating cones of ``self``
        containing the ``i``-th ray if ``i`` was given, :class:`tuple` of
        these sets for all rays otherwise.

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan._ray_to_cones(0)                                                  # needs palp
            frozenset({0, 3})
            sage: fan._ray_to_cones()                                                   # needs palp
            (frozenset({0, 3}), frozenset({1, 2}), frozenset({0, 1}), frozenset({2, 3}))
        """
        # This function is close to self(1)[i].star_generator_indices(), but
        # it does not require computation of the cone lattice and is
        # convenient for internal purposes.
        if "_ray_to_cones_tuple" not in self.__dict__:
            ray_to_cones = []
            for _ in self.rays():
                ray_to_cones.append([])
            for k, cone in enumerate(self):
                for j in cone.ambient_ray_indices():
                    ray_to_cones[j].append(k)
            self._ray_to_cones_tuple = tuple(frozenset(rtc)
                                             for rtc in ray_to_cones)
        if i is None:
            return self._ray_to_cones_tuple
        else:
            return self._ray_to_cones_tuple[i]

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: f = Fan(cones=[(0,1), (1,2)],
            ....:         rays=[(1,0), (0,1), (-1, 0)],
            ....:         check=False)
            sage: f._repr_()
            'Rational polyhedral fan in 2-d lattice N'
            sage: f = Fan(cones=[(0,1), (1,2)],
            ....:         rays=[(1,0), (0,1), (-1, 0)],
            ....:         lattice=ZZ^2,
            ....:         check=False)
            sage: f._repr_()
            'Rational polyhedral fan in 2-d lattice'
        """
        result = "Rational polyhedral fan in"
        if isinstance(self.lattice(), ToricLattice_generic):
            result += " %s" % self.lattice()
        else:
            result += " %d-d lattice" % self.lattice_dim()
        return result

    def _subdivide_stellar(self, new_rays, verbose):
        r"""
        Return iterative stellar subdivision of ``self`` via ``new_rays``.

        INPUT:

        - ``new_rays`` -- immutable primitive vectors in the lattice of
          ``self``

        - ``verbose`` -- if ``True``, some timing information will be printed

        OUTPUT: rational polyhedral fan

        TESTS::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: new_rays = sage.geometry.cone.normalize_rays([(1,1)], None)
            sage: fan = Fan([cone1, cone2])
            sage: fan._subdivide_stellar(new_rays, False)
            Rational polyhedral fan in 2-d lattice N
            sage: fan = Fan([cone1])
            sage: new_fan = fan._subdivide_stellar(new_rays, True)
            R:1/1  C:2  T:...(ms)  T/new:...(ms)  T/all:...(ms)
            sage: new_fan.rays()
            N(1, 0),
            N(0, 1),
            N(1, 1)
            in 2-d lattice N
            sage: for cone in new_fan: print(cone.ambient_ray_indices())
            (0, 2)
            (1, 2)

        We make sure that this function constructs cones with ordered ambient
        ray indices (see :issue:`9812`)::

            sage: C = Cone([(1,0,0), (0,1,0), (1,0,1), (0,1,1)])
            sage: F = Fan([C]).make_simplicial()
            sage: [cone.ambient_ray_indices() for cone in F]
            [(0, 1, 3), (0, 2, 3)]
        """
        cones = self.generating_cones()
        for n, ray in enumerate(new_rays):
            if verbose:
                start = walltime()
            new = []
            for cone in cones:
                if ray in cone:
                    new.extend(Cone(tuple(facet.rays())+(ray,), check=False)
                               for facet in cone.facets() if ray not in facet)
                else:
                    new.append(cone)
            if verbose:
                t = walltime(start)
                added = len(new) - len(cones)
                T_new = "%d" % (t / added * 1000) if added else "-"
                print("R:%d/%d  C:%d  T:%d(ms)  T/new:%s(ms)  T/all:%d(ms)"
                      % (n + 1, len(new_rays), len(new), t * 1000,
                         T_new, t / len(new) * 1000))
            cones = new
        new_fan_rays = list(self.rays())
        new_fan_rays.extend(ray for ray in new_rays
                                if ray not in self.rays().set())
        ray_to_index = {ray: i for i, ray in enumerate(new_fan_rays)}
        cones = tuple(tuple(sorted(ray_to_index[ray] for ray in cone))
                      for cone in cones)
        fan = Fan(cones, new_fan_rays, check=False, normalize=False)
        return fan

    def cone_containing(self, *points):
        r"""
        Return the smallest cone of ``self`` containing all given points.

        INPUT:

        - either one or more indices of rays of ``self``, or one or more
          objects representing points of the ambient space of ``self``, or a
          list of such objects (you CANNOT give a list of indices).

        OUTPUT:

        A :class:`cone of fan <Cone_of_fan>` whose ambient fan is
        ``self``

        .. NOTE::

            We think of the origin as of the smallest cone containing no rays
            at all. If there is no ray in ``self`` that contains all ``rays``,
            a :exc:`ValueError` exception will be raised.

        EXAMPLES::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])
            sage: f.rays()
            N(0, -1),
            N(0,  1),
            N(1,  0)
            in 2-d lattice N
            sage: f.cone_containing(0)  # ray index
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing(0, 1) # ray indices
            Traceback (most recent call last):
            ...
            ValueError: there is no cone in
            Rational polyhedral fan in 2-d lattice N
            containing all of the given rays! Ray indices: [0, 1]
            sage: f.cone_containing(0, 2) # ray indices
            2-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((0,1))  # point
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing([(0,1)]) # point
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((1,1))
            2-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((1,1), (1,0))
            2-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing()
            0-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((0,0))
            0-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((-1,1))
            Traceback (most recent call last):
            ...
            ValueError: there is no cone in
            Rational polyhedral fan in 2-d lattice N
            containing all of the given points! Points: [N(-1, 1)]

        TESTS::

            sage: fan = Fan(cones=[(0,1,2,3), (0,1,4)],
            ....:           rays=[(1,1,1), (1,-1,1), (1,-1,-1), (1,1,-1), (0,0,1)])
            sage: fan.cone_containing(0).rays()
            N(1, 1, 1)
            in 3-d lattice N
        """
        if not points:
            return self.cones(dim=0)[0]
        try:
            rays = [int(_) for _ in points]
            # Got ray indices
            generating_cones = set(range(self.ngenerating_cones()))
            for ray in rays:
                generating_cones.intersection_update(self._ray_to_cones(ray))
            if not generating_cones:
                raise ValueError("there is no cone in %s containing all of "
                        "the given rays! Ray indices: %s" % (self, rays))
            containing_cone = self.generating_cone(generating_cones.pop())
            for cone in generating_cones:
                containing_cone = containing_cone.intersection(
                                                self.generating_cone(cone))
            if not self.is_complete():
                # This cone may be too big in the case of incomplete fans
                rays = frozenset(rays)
                facets = containing_cone.facets()
                for facet in facets:
                    if rays.issubset(facet._ambient_ray_indices):
                        containing_cone = containing_cone.intersection(facet)
            return containing_cone
        except TypeError:
            # Got points (hopefully)
            try:
                points = [_ambient_space_point(self, p) for p in points]
            except TypeError:
                if len(points) == 1:
                    points = [_ambient_space_point(self, p) for p in points[0]]
                else:
                    raise
            # If we are still here, points are good
            # First we try to find a generating cone containing all points
            containing_cone = None
            for cone in self:
                contains_all = True
                for point in points:
                    if point not in cone:
                        contains_all = False
                        break
                if contains_all:
                    containing_cone = cone
                    break
            if containing_cone is None:
                raise ValueError("there is no cone in %s containing all of "
                            "the given points! Points: %s" % (self, points))
            # Now we take the intersection of facets that contain all points
            facets = containing_cone.facets()
            for facet in facets:
                contains_all = True
                for point in points:
                    if point not in facet:
                        contains_all = False
                        break
                if contains_all:
                    containing_cone = containing_cone.intersection(facet)
            return containing_cone

    def cone_lattice(self):
        r"""
        Return the cone lattice of ``self``.

        This lattice will have the origin as the bottom (we do not include the
        empty set as a cone) and the fan itself as the top.

        OUTPUT:

        :class:`finite poset <sage.combinat.posets.posets.FinitePoset` of
        :class:`cones of fan<Cone_of_fan>`, behaving like "regular" cones,
        but also containing the information about their relation to this
        fan, namely, the contained rays and containing generating cones. The
        top of the lattice will be this fan itself (*which is not a*
        :class:`cone of fan<Cone_of_fan>`).

        See also :meth:`cones`.

        EXAMPLES:

        Cone lattices can be computed for arbitrary fans::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.rays()
            N(-1, 0),
            N( 0, 1),
            N( 1, 0)
            in 2-d lattice N
            sage: for cone in fan: print(cone.ambient_ray_indices())
            (1, 2)
            (0,)
            sage: L = fan.cone_lattice()
            sage: L
            Finite poset containing 6 elements with distinguished linear extension

        These 6 elements are the origin, three rays, one two-dimensional
        cone, and the fan itself\ . Since we do add the fan itself as the
        largest face, you should be a little bit careful with this last
        element::

            sage: for face in L: print(face.ambient_ray_indices())
            Traceback (most recent call last):
            ...
            AttributeError: 'RationalPolyhedralFan'
            object has no attribute 'ambient_ray_indices'
            sage: L.top()
            Rational polyhedral fan in 2-d lattice N

        For example, you can do ::

            sage: for l in L.level_sets()[:-1]:
            ....:     print([f.ambient_ray_indices() for f in l])
            [()]
            [(0,), (1,), (2,)]
            [(1, 2)]

        If the fan is complete, its cone lattice is atomic and coatomic and
        can (and will!) be computed in a much more efficient way, but the
        interface is exactly the same::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: L = fan.cone_lattice()                                                # needs palp
            sage: for l in L.level_sets()[:-1]:                                         # needs palp
            ....:     print([f.ambient_ray_indices() for f in l])
            [()]
            [(0,), (1,), (2,), (3,)]
            [(0, 2), (1, 2), (1, 3), (0, 3)]

        Let's also consider the cone lattice of a fan generated by a single
        cone::

            sage: fan = Fan([cone1])
            sage: L = fan.cone_lattice()
            sage: L
            Finite poset containing 5 elements with distinguished linear extension

        Here these 5 elements correspond to the origin, two rays, one
        generating cone of dimension two, and the whole fan. While this single
        cone "is" the whole fan, it is consistent and convenient to
        distinguish them in the cone lattice.
        """
        if "_cone_lattice" not in self.__dict__:
            self._compute_cone_lattice()
        return self._cone_lattice

    def f_vector(self) -> tuple:
        r"""
        Return the f-vector of the fan.

        This is the tuple `(f_0, f_1, \ldots, f_d)`
        where `f_i` is the number of cones of dimension `i`.

        EXAMPLES::

            sage: F = ClusterAlgebra(['A',2]).cluster_fan()
            sage: F.f_vector()
            (1, 5, 5)
        """
        return tuple(len(d) for d in self.cones())

    # Internally we use this name for a uniform behaviour of cones and fans.
    _face_lattice_function = cone_lattice

    def __getstate__(self):
        r"""
        Return the dictionary that should be pickled.

        OUTPUT: :class:`dict`

        TESTS::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.cone_lattice()
            Finite poset containing 6 elements with distinguished linear extension
            sage: fan._test_pickling()
        """
        state = copy(self.__dict__)
        # TODO: do we want to keep the cone lattice in the pickle?
        # Currently there is an unpickling loop if do.
        # See Cone.__getstate__ for a similar problem and discussion.
        state.pop("_cone_lattice", None)
        return state

    def cones(self, dim=None, codim=None):
        r"""
        Return the specified cones of ``self``.

        INPUT:

        - ``dim`` -- dimension of the requested cones

        - ``codim`` -- codimension of the requested cones

        .. NOTE::

            You can specify at most one input parameter.

        OUTPUT:

        :class:`tuple` of cones of ``self`` of the specified (co)dimension,
        if either ``dim`` or ``codim`` is given. Otherwise :class:`tuple` of
        such tuples for all existing dimensions.

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan(dim=0)
            (0-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(codim=2)
            (0-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: for cone in fan.cones(1): cone.ray(0)
            N(-1, 0)
            N(0, 1)
            N(1, 0)
            sage: fan.cones(2)
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)

        You cannot specify both dimension and codimension, even if they
        "agree"::

            sage: fan(dim=1, codim=1)
            Traceback (most recent call last):
            ...
            ValueError: dimension and codimension
            cannot be specified together!

        But it is OK to ask for cones of too high or low (co)dimension::

            sage: fan(-1)
            ()
            sage: fan(3)
            ()
            sage: fan(codim=4)
            ()
        """
        if "_cones" not in self.__dict__:
            levels = self.cone_lattice().level_sets()
            levels.pop()  # The very last level is this FAN, not cone.
            # It seems that there is no reason to believe that the order of
            # faces in level sets has anything to do with the order of
            # vertices in the Hasse diagram of FinitePoset. So, while
            # lattice_from_incidences tried to ensure a "good order,"
            # we will sort faces corresponding to rays, as well as faces
            # corresponding to generating cones, if they are all of the same
            # dimension (otherwise it is not very useful).
            if len(levels) >= 3:  # There are cones of dimension higher than 1
                top_cones = list(levels[-1])
                if len(top_cones) == self.ngenerating_cones():
                    top_cones.sort(key=lambda cone:
                                            cone.star_generator_indices()[0])
                levels[-1] = top_cones
            if len(levels) >= 2:  # We have rays
                rays = list(levels[1])
                rays.sort(key=lambda cone: cone.ambient_ray_indices()[0])
                levels[1] = rays
            self._cones = tuple(tuple(level) for level in levels)
        if dim is None:
            if codim is None:
                return self._cones
            dim = self.dim() - codim
        elif codim is not None:
            raise ValueError(
                    "dimension and codimension cannot be specified together!")
        return self._cones[dim] if 0 <= dim < len(self._cones) else ()

    def contains(self, cone) -> bool:
        r"""
        Check if a given ``cone`` is equivalent to a cone of the fan.

        INPUT:

        - ``cone`` -- anything

        OUTPUT:

        ``False`` if ``cone`` is not a cone or if ``cone`` is not
        equivalent to a cone of the fan, ``True`` otherwise

        .. NOTE::

            Recall that a fan is a (finite) collection of cones. A
            cone is contained in a fan if it is equivalent to one of
            the cones of the fan. In particular, it is possible that
            all rays of the cone are in the fan, but the cone itself
            is not.

            If you want to know whether a point is in the support of
            the fan, you should use :meth:`support_contains`.

        EXAMPLES:

        We first construct a simple fan::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])

        Now we check if some cones are in this fan. First, we make sure that
        the order of rays of the input cone does not matter (``check=False``
        option ensures that rays of these cones will be listed exactly as they
        are given)::

            sage: f.contains(Cone([(1,0), (0,1)], check=False))
            True
            sage: f.contains(Cone([(0,1), (1,0)], check=False))
            True

        Now we check that a non-generating cone is in our fan::

            sage: f.contains(Cone([(1,0)]))
            True
            sage: Cone([(1,0)]) in f   # equivalent to the previous command
            True

        Finally, we test some cones which are not in this fan::

            sage: f.contains(Cone([(1,1)]))
            False
            sage: f.contains(Cone([(1,0), (-0,1)]))
            True

        A point is not a cone::

            sage: n = f.lattice()(1,1); n
            N(1, 1)
            sage: f.contains(n)
            False
        """
        return self._contains(cone)

    def embed(self, cone):
        r"""
        Return the cone equivalent to the given one, but sitting in ``self``.

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

        a :class:`cone of fan <Cone_of_fan>`, equivalent to ``cone`` but
        sitting inside ``self``

        EXAMPLES:

        Let's take a 3-d fan generated by a cone on 4 rays::

            sage: f = Fan([Cone([(1,0,1), (0,1,1), (-1,0,1), (0,-1,1)])])

        Then any ray generates a 1-d cone of this fan, but if you construct
        such a cone directly, it will not "sit" inside the fan::

            sage: ray = Cone([(0,-1,1)])
            sage: ray
            1-d cone in 3-d lattice N
            sage: ray.ambient_ray_indices()
            (0,)
            sage: ray.adjacent()
            ()
            sage: ray.ambient()
            1-d cone in 3-d lattice N

        If we want to operate with this ray as a part of the fan, we need to
        embed it first::

            sage: e_ray = f.embed(ray)
            sage: e_ray
            1-d cone of Rational polyhedral fan in 3-d lattice N
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
            (1-d cone of Rational polyhedral fan in 3-d lattice N,
             1-d cone of Rational polyhedral fan in 3-d lattice N)
            sage: e_ray.ambient()
            Rational polyhedral fan in 3-d lattice N

        Not every cone can be embedded into a fixed fan::

            sage: f.embed(Cone([(0,0,1)]))
            Traceback (most recent call last):
            ...
            ValueError: 1-d cone in 3-d lattice N does not belong
            to Rational polyhedral fan in 3-d lattice N!
            sage: f.embed(Cone([(1,0,1), (-1,0,1)]))
            Traceback (most recent call last):
            ...
            ValueError: 2-d cone in 3-d lattice N does not belong
            to Rational polyhedral fan in 3-d lattice N!
        """
        if not isinstance(cone, sage.geometry.abc.ConvexRationalPolyhedralCone):
            raise TypeError("%s is not a cone!" % cone)
        if cone.ambient() is self:
            return cone
        rays = self.rays()
        try:
            # Compute ray indices.
            ray_indices = [rays.index(ray) for ray in cone.rays()]
            # Get the smallest cone containing them
            result = self.cone_containing(*ray_indices)
            # If there is a cone containing all of the rays of the given cone,
            # they must be among its generating rays and we only need to worry
            # if there are any extra ones.
            if cone.nrays() != result.nrays():
                raise ValueError
        except ValueError:
            raise ValueError("%s does not belong to %s!" % (cone, self))
        return result

    @cached_method
    def Gale_transform(self):
        r"""
        Return the Gale transform of ``self``.

        OUTPUT: a matrix over `ZZ`

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.Gale_transform()                                                  # needs palp
            [ 1  1  0  0 -2]
            [ 0  0  1  1 -2]
            sage: _.base_ring()                                                         # needs palp
            Integer Ring
        """
        m = self.rays().matrix().stack(matrix(ZZ, 1, self.lattice_dim()))
        m = m.augment(matrix(ZZ, m.nrows(), 1, [1] * m.nrows()))
        return matrix(ZZ, m.integer_kernel().matrix())

    def is_polytopal(self):
        r"""
        Check if ``self`` is the normal fan of a polytope.

        A rational polyhedral fan is *polytopal* if it is the normal fan of a
        polytope. This is also called *regular*, or provide a *coherent*
        subdivision or leads to a *projective* toric variety.

        OUTPUT: ``True`` if ``self`` is polytopal and ``False`` otherwise

        EXAMPLES:

        This is the mother of all examples (see Section 7.1.1 in
        [DLRS2010]_)::

            sage: def mother(epsilon=0):
            ....:     rays = [(4-epsilon,epsilon,0),(0,4-epsilon,epsilon),(epsilon,0,4-epsilon),(2,1,1),(1,2,1),(1,1,2),(-1,-1,-1)]
            ....:     L = [(0,1,4),(0,3,4),(1,2,5),(1,4,5),(0,2,3),(2,3,5),(3,4,5),(6,0,1),(6,1,2),(6,2,0)]
            ....:     S1 = [Cone([rays[i] for i in indices]) for indices in L]
            ....:     return Fan(S1)

        When epsilon=0, it is not polytopal::

            sage: epsilon = 0
            sage: mother(epsilon).is_polytopal()
            False

        Doing a slight perturbation makes the same subdivision polytopal::

            sage: epsilon = 1/2
            sage: mother(epsilon).is_polytopal()
            True

        .. SEEALSO::

            :meth:`is_projective`.
        """
        if not self.is_complete():
            raise ValueError('To be polytopal, the fan should be complete.')
        from sage.geometry.triangulation.point_configuration import PointConfiguration
        from sage.geometry.polyhedron.constructor import Polyhedron
        pc = PointConfiguration(self.rays())
        v_pc = [tuple(p) for p in pc]
        pc_to_indices = {tuple(p):i for (i,p) in enumerate(pc)}
        indices_to_vr = [tuple(r) for r in self.rays()]
        cone_indices = [cone.ambient_ray_indices() for cone in self.generating_cones()]
        translator = [pc_to_indices[t] for t in indices_to_vr]
        translated_cone_indices = [[translator[i] for i in ci] for ci in cone_indices]
        dc_pc = pc.deformation_cone(translated_cone_indices)
        lift = dc_pc.an_element()
        ieqs = [(lift_i,) + v for (lift_i, v) in zip(lift, v_pc)]
        poly = Polyhedron(ieqs=ieqs)
        return self.is_equivalent(poly.normal_fan())

    def generating_cone(self, n):
        r"""
        Return the ``n``-th generating cone of ``self``.

        INPUT:

        - ``n`` -- integer; the index of a generating cone

        OUTPUT: :class:`cone of fan<Cone_of_fan>`

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.generating_cone(0)                                                # needs palp
            2-d cone of Rational polyhedral fan in 2-d lattice N
        """
        return self._generating_cones[n]

    def generating_cones(self):
        r"""
        Return generating cones of ``self``.

        OUTPUT: :class:`tuple` of :class:`cones of fan<Cone_of_fan>`

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.generating_cones()                                                # needs palp
            (2-d cone of Rational polyhedral fan in 2-d lattice N,
             2-d cone of Rational polyhedral fan in 2-d lattice N,
             2-d cone of Rational polyhedral fan in 2-d lattice N,
             2-d cone of Rational polyhedral fan in 2-d lattice N)
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.generating_cones()
            (2-d cone of Rational polyhedral fan in 2-d lattice N,
             1-d cone of Rational polyhedral fan in 2-d lattice N)
        """
        return self._generating_cones

    @cached_method
    def vertex_graph(self):
        r"""
        Return the graph of 1- and 2-cones.

        OUTPUT:

        An edge-colored graph. The vertices correspond to the 1-cones
        (i.e. rays) of
        the fan. Two vertices are joined by an edge iff the rays span
        a 2-cone of the fan. The edges are colored by pairs of
        integers that classify the 2-cones up to `GL(2,\ZZ)`
        transformation, see
        :func:`~sage.geometry.cone.classify_cone_2d`.

        EXAMPLES::

            sage: # needs palp
            sage: dP8 = toric_varieties.dP8()
            sage: g = dP8.fan().vertex_graph(); g
            Graph on 4 vertices
            sage: set(dP8.fan(1)) == set(g.vertices(sort=False))
            True
            sage: g.edge_labels()  # all edge labels the same since every cone is smooth
            [(1, 0), (1, 0), (1, 0), (1, 0)]

            sage: g = toric_varieties.Cube_deformation(10).fan().vertex_graph()
            sage: g.automorphism_group().order()                                        # needs sage.groups
            48
            sage: g.automorphism_group(edge_labels=True).order()                        # needs sage.groups
            4
        """
        from sage.geometry.cone import classify_cone_2d
        graph = {}
        cones_1d = list(self(1))
        while cones_1d:
            c0 = cones_1d.pop()
            c0_edges = {}
            for c1 in c0.adjacent():
                if c1 not in cones_1d:
                    continue
                label = classify_cone_2d(c0.ray(0), c1.ray(0), check=False)
                c0_edges[c1] = label
            graph[c0] = c0_edges
        from sage.graphs.graph import Graph
        return Graph(graph)

    def is_complete(self) -> bool:
        r"""
        Check if ``self`` is complete.

        A rational polyhedral fan is *complete* if its cones fill the whole
        space.

        OUTPUT: ``True`` if ``self`` is complete and ``False`` otherwise

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.is_complete()                                                     # needs palp
            True
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.is_complete()
            False
        """
        if "_is_complete" in self.__dict__:
            return self._is_complete
        d = self.lattice_dim()
        if self.dim() != d:
            self._is_complete = False
            return False
        for cone in self:
            if cone.dim() != d:
                self._is_complete = False
                return False
        # Now we know that all generating cones are full-dimensional.
        # Then boundary cones are (d-1)-dimensional.
        for cone in self(codim=1):
            if len(cone.star_generator_indices()) != 2:
                self._is_complete = False
                return False
        self._is_complete = True
        return True

    def is_equivalent(self, other) -> bool:
        r"""
        Check if ``self`` is "mathematically" the same as ``other``.

        INPUT:

        - ``other`` -- fan

        OUTPUT:

        ``True`` if ``self`` and ``other`` define the same fans as
        collections of equivalent cones in the same lattice, ``False``
        otherwise.

        There are three different equivalences between fans `F_1` and `F_2`
        in the same lattice:

        #. They have the same rays in the same order and the same generating
           cones in the same order.
           This is tested by ``F1 == F2``.
        #. They have the same rays and the same generating cones without
           taking into account any order.
           This is tested by ``F1.is_equivalent(F2)``.
        #. They are in the same orbit of `GL(n,\ZZ)` (and, therefore,
           correspond to isomorphic toric varieties).
           This is tested by ``F1.is_isomorphic(F2)``.

        Note that :meth:`virtual_rays` are included into consideration for all
        of the above equivalences.

        EXAMPLES::

            sage: fan1 = Fan(cones=[(0,1), (1,2)],
            ....:            rays=[(1,0), (0,1), (-1,-1)],
            ....:            check=False)
            sage: fan2 = Fan(cones=[(2,1), (0,2)],
            ....:            rays=[(1,0), (-1,-1), (0,1)],
            ....:            check=False)
            sage: fan3 = Fan(cones=[(0,1), (1,2)],
            ....:            rays=[(1,0), (0,1), (-1,1)],
            ....:            check=False)
            sage: fan1 == fan2
            False
            sage: fan1.is_equivalent(fan2)
            True
            sage: fan1 == fan3
            False
            sage: fan1.is_equivalent(fan3)
            False
        """
        if (self.lattice() != other.lattice()
              or self.dim() != other.dim()
              or self.ngenerating_cones() != other.ngenerating_cones()
              or self.rays().set() != other.rays().set()
              or self.virtual_rays().set() != other.virtual_rays().set()):
            return False
        # Now we need to really compare cones, which can take a while
        return sorted(sorted(cone.rays()) for cone in self) \
               == sorted(sorted(cone.rays()) for cone in other)

    def is_isomorphic(self, other) -> bool:
        r"""
        Check if ``self`` is in the same `GL(n, \ZZ)`-orbit as ``other``.

        There are three different equivalences between fans `F_1` and `F_2`
        in the same lattice:

        #. They have the same rays in the same order and the same generating
           cones in the same order.
           This is tested by ``F1 == F2``.
        #. They have the same rays and the same generating cones without
           taking into account any order.
           This is tested by ``F1.is_equivalent(F2)``.
        #. They are in the same orbit of `GL(n,\ZZ)` (and, therefore,
           correspond to isomorphic toric varieties).
           This is tested by ``F1.is_isomorphic(F2)``.

        Note that :meth:`virtual_rays` are included into consideration for all
        of the above equivalences.

        INPUT:

        - ``other`` -- a :class:`fan <RationalPolyhedralFan>`

        OUTPUT:

        ``True`` if ``self`` and ``other`` are in the same
        `GL(n, \ZZ)`-orbit, ``False`` otherwise

        .. SEEALSO::

            If you want to obtain the actual fan isomorphism, use
            :meth:`isomorphism`.

        EXAMPLES:

        Here we pick an `SL(2,\ZZ)` matrix ``m`` and then verify that
        the image fan is isomorphic::

            sage: rays = ((1, 1), (0, 1), (-1, -1), (1, 0))
            sage: cones = [(0,1), (1,2), (2,3), (3,0)]
            sage: fan1 = Fan(cones, rays)
            sage: m = matrix([[-2,3], [1,-1]])
            sage: fan2 = Fan(cones, [vector(r)*m for r in rays])
            sage: fan1.is_isomorphic(fan2)
            True
            sage: fan1.is_equivalent(fan2)
            False
            sage: fan1 == fan2
            False

        These fans are "mirrors" of each other::

            sage: fan1 = Fan(cones=[(0,1), (1,2)],
            ....:            rays=[(1,0), (0,1), (-1,-1)],
            ....:            check=False)
            sage: fan2 = Fan(cones=[(0,1), (1,2)],
            ....:            rays=[(1,0), (0,-1), (-1,1)],
            ....:            check=False)
            sage: fan1 == fan2
            False
            sage: fan1.is_equivalent(fan2)
            False
            sage: fan1.is_isomorphic(fan2)
            True
            sage: fan1.is_isomorphic(fan1)
            True
        """
        from sage.geometry.fan_isomorphism import \
            fan_isomorphic_necessary_conditions, fan_isomorphism_generator
        if not fan_isomorphic_necessary_conditions(self, other):
            return False
        if self.lattice_dim() == 2:
            if self._2d_echelon_forms.cache is None:
                return self._2d_echelon_form() in other._2d_echelon_forms()
            else:
                return other._2d_echelon_form() in self._2d_echelon_forms()
        generator = fan_isomorphism_generator(self, other)
        try:
            next(generator)
            return True
        except StopIteration:
            return False

    @cached_method
    def _2d_echelon_forms(self):
        """
        Return all echelon forms of the cyclically ordered rays of a 2-d fan.

        OUTPUT: a set of integer matrices

        EXAMPLES::

            sage: fan = toric_varieties.dP8().fan()                                     # needs palp
            sage: fan._2d_echelon_forms()                                               # needs palp
            frozenset({[ 1  0 -1 -1]
                       [ 0  1  0 -1], [ 1  0 -1  0]
                       [ 0  1 -1 -1], [ 1  0 -1  0]
                       [ 0  1  1 -1], [ 1  0 -1  1]
                       [ 0  1  0 -1]})
        """
        from sage.geometry.fan_isomorphism import fan_2d_echelon_forms
        return fan_2d_echelon_forms(self)

    @cached_method
    def _2d_echelon_form(self):
        """
        Return the echelon form of one particular cyclic order of rays of a 2-d fan.

        OUTPUT: integer matrix whose columns are the rays in the echelon form

        EXAMPLES::

            sage: fan = toric_varieties.dP8().fan()                                     # needs palp
            sage: fan._2d_echelon_form()                                                # needs palp
            [ 1  0 -1 -1]
            [ 0  1  0 -1]
        """
        from sage.geometry.fan_isomorphism import fan_2d_echelon_form
        return fan_2d_echelon_form(self)

    def isomorphism(self, other):
        r"""
        Return a fan isomorphism from ``self`` to ``other``.

        INPUT:

        - ``other`` -- fan

        OUTPUT:

        A fan isomorphism. If no such isomorphism exists, a
        :class:`~sage.geometry.fan_isomorphism.FanNotIsomorphicError`
        is raised.

        EXAMPLES::

            sage: rays = ((1, 1), (0, 1), (-1, -1), (3, 1))
            sage: cones = [(0,1), (1,2), (2,3), (3,0)]
            sage: fan1 = Fan(cones, rays)
            sage: m = matrix([[-2,3], [1,-1]])
            sage: fan2 = Fan(cones, [vector(r)*m for r in rays])

            sage: fan1.isomorphism(fan2)
            Fan morphism defined by the matrix
            [-2  3]
            [ 1 -1]
            Domain fan: Rational polyhedral fan in 2-d lattice N
            Codomain fan: Rational polyhedral fan in 2-d lattice N

            sage: fan2.isomorphism(fan1)
            Fan morphism defined by the matrix
            [1 3]
            [1 2]
            Domain fan: Rational polyhedral fan in 2-d lattice N
            Codomain fan: Rational polyhedral fan in 2-d lattice N

            sage: fan1.isomorphism(toric_varieties.P2().fan())                          # needs palp
            Traceback (most recent call last):
            ...
            FanNotIsomorphicError
        """
        from sage.geometry.fan_isomorphism import find_isomorphism
        return find_isomorphism(self, other, check=False)

    def is_simplicial(self) -> bool:
        r"""
        Check if ``self`` is simplicial.

        A rational polyhedral fan is **simplicial** if all of its cones are,
        i.e. primitive vectors along generating rays of every cone form a part
        of a *rational* basis of the ambient space.

        OUTPUT: ``True`` if ``self`` is simplicial and ``False`` otherwise

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.is_simplicial()                                                   # needs palp
            True
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.is_simplicial()
            True

        In fact, any fan in a two-dimensional ambient space is simplicial.
        This is no longer the case in dimension three::

            sage: fan = NormalFan(lattice_polytope.cross_polytope(3))
            sage: fan.is_simplicial()
            False
            sage: fan.generating_cone(0).nrays()
            4
        """
        if "is_simplicial" not in self.__dict__:
            self._is_simplicial = all(cone.is_simplicial() for cone in self)
        return self._is_simplicial

    @cached_method
    def is_smooth(self, codim=None) -> bool:
        r"""
        Check if ``self`` is smooth.

        A rational polyhedral fan is **smooth** if all of its cones
        are, i.e. primitive vectors along generating rays of every
        cone form a part of an *integral* basis of the ambient
        space. In this case the corresponding toric variety is smooth.

        A fan in an `n`-dimensional lattice is smooth up to codimension `c`
        if all cones of codimension greater than or equal to `c` are smooth,
        i.e. if all cones of dimension less than or equal to `n-c` are smooth.
        In this case the singular set of the corresponding toric variety is of
        dimension less than `c`.

        INPUT:

        - ``codim`` -- codimension in which smoothness has to be checked, by
          default complete smoothness will be checked

        OUTPUT:

        ``True`` if ``self`` is smooth (in codimension ``codim``, if it was
        given) and ``False`` otherwise.

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.is_smooth()                                                       # needs palp
            True
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.is_smooth()
            True
            sage: fan = NormalFan(lattice_polytope.cross_polytope(2))
            sage: fan.is_smooth()
            False
            sage: fan.is_smooth(codim=1)
            True
            sage: fan.generating_cone(0).rays()
            N(-1, -1),
            N(-1,  1)
            in 2-d lattice N
            sage: fan.generating_cone(0).rays().matrix().det()
            -2
        """
        if codim is None or codim < 0:
            codim = 0
        if codim > self.lattice_dim() - 2:
            return True
        return all(cone.is_smooth() for cone in self(codim=codim)) and \
               self.is_smooth(codim + 1)

    def make_simplicial(self, **kwds):
        r"""
        Construct a simplicial fan subdividing ``self``.

        It is a synonym for :meth:`subdivide` with ``make_simplicial=True``
        option.

        INPUT:

        - this functions accepts only keyword arguments. See :meth:`subdivide`
          for documentation.

        OUTPUT:

        :class:`rational polyhedral fan
        <sage.geometry.fan.RationalPolyhedralFan>`

        EXAMPLES::

            sage: fan = NormalFan(lattice_polytope.cross_polytope(3))
            sage: fan.is_simplicial()
            False
            sage: fan.ngenerating_cones()
            6
            sage: new_fan = fan.make_simplicial()
            sage: new_fan.is_simplicial()
            True
            sage: new_fan.ngenerating_cones()
            12
        """
        return self.subdivide(make_simplicial=True, **kwds)

    def ngenerating_cones(self):
        r"""
        Return the number of generating cones of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: fan = toric_varieties.P1xP1().fan()                                   # needs palp
            sage: fan.ngenerating_cones()                                               # needs palp
            4
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.ngenerating_cones()
            2
        """
        return len(self.generating_cones())

    def plot(self, **options):
        r"""
        Plot ``self``.

        INPUT:

        - any options for toric plots (see :func:`toric_plotter.options
          <sage.geometry.toric_plotter.options>`), none are mandatory.

        OUTPUT: a plot

        EXAMPLES::

            sage: fan = toric_varieties.dP6().fan()                                     # needs palp
            sage: fan.plot()                                                            # needs palp sage.plot
            Graphics object consisting of 31 graphics primitives
        """
        tp = ToricPlotter(options, self.lattice().degree(), self.rays())
        result = tp.plot_lattice() + tp.plot_rays() + tp.plot_generators()
        if self.dim() >= 2:
            result += tp.plot_walls(self(2))
        return result

    def subdivide(self, new_rays=None, make_simplicial=False,
                  algorithm='default', verbose=False):
        r"""
        Construct a new fan subdividing ``self``.

        INPUT:

        - ``new_rays`` -- list of new rays to be added during subdivision, each
          ray must be a list or a vector. May be empty or ``None`` (default);

        - ``make_simplicial`` -- if ``True``, the returned fan is guaranteed to
          be simplicial, default is ``False``;

        - ``algorithm`` -- string with the name of the algorithm used for
          subdivision. Currently there is only one available algorithm called
          "default";

        - ``verbose`` -- if ``True``, some timing information may be printed
          during the process of subdivision

        OUTPUT:

        :class:`rational polyhedral fan
        <sage.geometry.fan.RationalPolyhedralFan>`

        Currently the "default" algorithm corresponds to iterative stellar
        subdivision for each ray in ``new_rays``.

        EXAMPLES::

            sage: fan = NormalFan(lattice_polytope.cross_polytope(3))
            sage: fan.is_simplicial()
            False
            sage: fan.ngenerating_cones()
            6
            sage: fan.nrays()
            8
            sage: new_fan = fan.subdivide(new_rays=[(1,0,0)])
            sage: new_fan.is_simplicial()
            False
            sage: new_fan.ngenerating_cones()
            9
            sage: new_fan.nrays()
            9

        TESTS:

        We check that :issue:`11902` is fixed::

            sage: fan = toric_varieties.P2().fan()                                      # needs palp
            sage: fan.subdivide(new_rays=[(0,0)])                                       # needs palp
            Traceback (most recent call last):
            ...
            ValueError: the origin cannot be used for fan subdivision!
        """
        # Maybe these decisions should be done inside the algorithms
        # We can figure it out once we have at least two of them.
        if make_simplicial and not self.is_simplicial():
            rays = list(self.rays())
        else:
            rays = []
        rays.extend(ray for ray in normalize_rays(new_rays, self.lattice())
                        if ray not in self.rays().set())
        if not rays:
            return self  # Nothing has to be done
        if self.lattice().zero() in rays:
            raise ValueError("the origin cannot be used for fan subdivision!")
        if algorithm == "default":
            algorithm = "stellar"
        method_name = "_subdivide_" + algorithm
        if not hasattr(self, method_name):
            raise ValueError('"%s" is an unknown subdivision algorithm!'
                             % algorithm)
        return getattr(self, method_name)(rays, verbose)

    def virtual_rays(self, *args):
        r"""
        Return (some of the) virtual rays of ``self``.

        Let `N` be the `D`-dimensional
        :meth:`~sage.geometry.cone.IntegralRayCollection.lattice`
        of a `d`-dimensional fan `\Sigma` in `N_\RR`. Then the corresponding
        toric variety is of the form `X \times (\CC^*)^{D-d}`. The actual
        :meth:`~sage.geometry.cone.IntegralRayCollection.rays` of `\Sigma`
        give a canonical choice of homogeneous coordinates on `X`. This function
        returns an arbitrary but fixed choice of virtual rays corresponding to a
        (non-canonical) choice of homogeneous coordinates on the torus factor.
        Combinatorially primitive integral generators of virtual rays span the
        `D-d` dimensions of `N_\QQ` "missed" by the actual rays. (In general
        addition of virtual rays is not sufficient to span `N` over `\ZZ`.)

        .. NOTE::

            You may use a particular choice of virtual rays by passing optional
            argument ``virtual_rays`` to the :func:`Fan` constructor.

        INPUT:

        - ``ray_list`` -- list of integers; the indices of the
          requested virtual rays. If not specified, all virtual rays of ``self``
          will be returned.

        OUTPUT:

        a :class:`~sage.geometry.point_collection.PointCollection` of
        primitive integral ray generators. Usually (if the fan is
        full-dimensional) this will be empty.

        EXAMPLES::

            sage: f = Fan([Cone([(1,0,1,0), (0,1,1,0)])])
            sage: f.virtual_rays()
            N(1, 0, 0, 0),
            N(0, 0, 0, 1)
            in 4-d lattice N

            sage: f.rays()
            N(1, 0, 1, 0),
            N(0, 1, 1, 0)
            in 4-d lattice N

            sage: f.virtual_rays([0])
            N(1, 0, 0, 0)
            in 4-d lattice N

        You can also give virtual ray indices directly, without
        packing them into a list::

            sage: f.virtual_rays(0)
            N(1, 0, 0, 0)
            in 4-d lattice N

        Make sure that :issue:`16344` is fixed and one can compute
        the virtual rays of fans in non-saturated lattices::

            sage: N = ToricLattice(1)
            sage: B = N.submodule([(2,)]).basis()
            sage: f = Fan([Cone([B[0]])])
            sage: len(f.virtual_rays())
            0

        TESTS::

            sage: N = ToricLattice(4)
            sage: for i in range(10):
            ....:      c = Cone([N.random_element() for j in range(i//2)], lattice=N)
            ....:      if not c.is_strictly_convex():
            ....:          continue
            ....:      f = Fan([c])
            ....:      assert matrix(f.rays() + f.virtual_rays()).rank() == 4
            ....:      assert f.dim() + len(f.virtual_rays()) == 4
        """
        try:
            virtual = self._virtual_rays
        except AttributeError:
            N = self.lattice()
            Np = N.ambient_module()
            qp = Np.quotient(self.rays().matrix().saturation().rows())
            quotient = qp.submodule(N.gens())
            virtual = [gen.lift() for gen in quotient.gens()]
            for v in virtual:
                v.set_immutable()
            virtual = PointCollection(virtual, N)
            self._virtual_rays = virtual
        if args:
            return virtual(*args)
        else:
            return virtual

    def primitive_collections(self):
        r"""
        Return the primitive collections.

        OUTPUT:

        Return the subsets `\{i_1,\dots,i_k\} \subset \{ 1,\dots,n\}`
        such that

        * The points `\{p_{i_1},\dots,p_{i_k}\}` do not span a cone of
          the fan.

        * If you remove any one `p_{i_j}` from the set, then they do
          span a cone of the fan.

        .. NOTE::

            By replacing the multiindices `\{i_1,\dots,i_k\}` of each
            primitive collection with the monomials `x_{i_1}\cdots
            x_{i_k}` one generates the Stanley-Reisner ideal in
            `\ZZ[x_1,\dots]`.

        REFERENCES:

        - [Bat1991]_

        EXAMPLES::

            sage: fan = Fan([[0,1,3], [3,4], [2,0], [1,2,4]],
            ....:           [(-3, -2, 1), (0, 0, 1), (3, -2, 1), (-1, -1, 1), (1, -1, 1)])
            sage: fan.primitive_collections()
            [frozenset({0, 4}),
             frozenset({2, 3}),
             frozenset({0, 1, 2}),
             frozenset({1, 3, 4})]
        """
        try:
            return self._primitive_collections
        except AttributeError:
            pass

        def is_not_facet(I):
            return all(not (I <= f) for f in facets)

        def is_in_SR(I):
            return all(not (I >= sr) for sr in SR)

        # Generators of SR are index sets I = {i1, ..., ik}
        # called "primitive collections" such that
        # 1) I is not contained in a face
        # 2) if you remove any one entry j, then I-{j} is contained in a facet
        facets = [frozenset(c.ambient_ray_indices())
                  for c in self.generating_cones()]
        all_points = frozenset(range(self.nrays()))
        d_max = max(map(len, facets)) + 1
        SR = []
        for d in range(1, d_max):
            checked = set()
            for facet in facets:
                for I_minus_j_list in Combinations(facet, d):
                    I_minus_j = frozenset(I_minus_j_list)
                    for j in all_points - I_minus_j:
                        I = I_minus_j.union(frozenset([j]))

                        if I in checked:
                            continue
                        else:
                            checked.add(I)

                        if is_not_facet(I) and is_in_SR(I):
                            SR.append(I)

        self._primitive_collections = SR
        return self._primitive_collections

    def Stanley_Reisner_ideal(self, ring):
        """
        Return the Stanley-Reisner ideal.

        INPUT:

        - A polynomial ring in ``self.nrays()`` variables.

        OUTPUT: the Stanley-Reisner ideal in the given polynomial ring

        EXAMPLES::

            sage: fan = Fan([[0,1,3], [3,4], [2,0], [1,2,4]],
            ....:           [(-3, -2, 1), (0, 0, 1), (3, -2, 1), (-1, -1, 1), (1, -1, 1)])
            sage: fan.Stanley_Reisner_ideal(PolynomialRing(QQ, 5, 'A, B, C, D, E'))
            Ideal (A*E, C*D, A*B*C, B*D*E) of
             Multivariate Polynomial Ring in A, B, C, D, E over Rational Field
        """
        generators_indices = self.primitive_collections()
        return ring.ideal([prod([ring.gen(i) for i in sr])
                           for sr in generators_indices])

    def linear_equivalence_ideal(self, ring):
        """
        Return the ideal generated by linear relations.

        INPUT:

        - A polynomial ring in ``self.nrays()`` variables.

        OUTPUT:

        Return the ideal, in the given ``ring``, generated by the
        linear relations of the rays. In toric geometry, this
        corresponds to rational equivalence of divisors.

        EXAMPLES::

            sage: fan = Fan([[0,1,3],[3,4],[2,0],[1,2,4]],
            ....:           [(-3, -2, 1), (0, 0, 1), (3, -2, 1), (-1, -1, 1), (1, -1, 1)])
            sage: fan.linear_equivalence_ideal(PolynomialRing(QQ, 5, 'A, B, C, D, E'))
            Ideal (-3*A + 3*C - D + E, -2*A - 2*C - D - E, A + B + C + D + E) of
             Multivariate Polynomial Ring in A, B, C, D, E over Rational Field
        """
        gens = []
        for d in range(self.dim()):
            gens.append(sum([self.ray(i)[d] * ring.gen(i)
                             for i in range(self.nrays())]))
        return ring.ideal(gens)

    def oriented_boundary(self, cone):
        r"""
        Return the facets bounding ``cone`` with their induced
        orientation.

        INPUT:

        - ``cone`` -- a cone of the fan or the whole fan

        OUTPUT:

        The boundary cones of ``cone`` as a formal linear combination
        of cones with coefficients `\pm 1`. Each summand is a facet of
        ``cone`` and the coefficient indicates whether their (chosen)
        orientation agrees or disagrees with the "outward normal
        first" boundary orientation. Note that the orientation of any
        individual cone is arbitrary. This method once and for all
        picks orientations for all cones and then computes the
        boundaries relative to that chosen orientation.

        If ``cone`` is the fan itself, the generating cones with their
        orientation relative to the ambient space are returned.

        See :meth:`complex` for the associated chain complex. If you
        do not require the orientation, use :meth:`cone.facets()
        <sage.geometry.cone.ConvexRationalPolyhedralCone.facets>`
        instead.

        EXAMPLES::

            sage: # needs palp
            sage: fan = toric_varieties.P(3).fan()
            sage: cone = fan(2)[0]
            sage: bdry = fan.oriented_boundary(cone);  bdry
            -1-d cone of Rational polyhedral fan in 3-d lattice N
            + 1-d cone of Rational polyhedral fan in 3-d lattice N
            sage: bdry[0]
            (-1, 1-d cone of Rational polyhedral fan in 3-d lattice N)
            sage: bdry[1]
            (1, 1-d cone of Rational polyhedral fan in 3-d lattice N)
            sage: fan.oriented_boundary(bdry[0][1])
            -0-d cone of Rational polyhedral fan in 3-d lattice N
            sage: fan.oriented_boundary(bdry[1][1])
            -0-d cone of Rational polyhedral fan in 3-d lattice N

        If you pass the fan itself, this method returns the
        orientation of the generating cones which is determined by the
        order of the rays in :meth:`cone.ray_basis()
        <sage.geometry.cone.IntegralRayCollection.ray_basis>` ::

            sage: fan.oriented_boundary(fan)                                            # needs palp
            -3-d cone of Rational polyhedral fan in 3-d lattice N
            + 3-d cone of Rational polyhedral fan in 3-d lattice N
            - 3-d cone of Rational polyhedral fan in 3-d lattice N
            + 3-d cone of Rational polyhedral fan in 3-d lattice N
            sage: [cone.rays().basis().matrix().det()                                   # needs palp
            ....:  for cone in fan.generating_cones()]
            [-1, 1, -1, 1]

        A non-full dimensional fan::

            sage: cone = Cone([(4,5)])
            sage: fan = Fan([cone])
            sage: fan.oriented_boundary(cone)
            0-d cone of Rational polyhedral fan in 2-d lattice N
            sage: fan.oriented_boundary(fan)
            1-d cone of Rational polyhedral fan in 2-d lattice N

        TESTS::

            sage: fan = toric_varieties.P2().fan()                                      # needs palp
            sage: trivial_cone = fan(0)[0]                                              # needs palp
            sage: fan.oriented_boundary(trivial_cone)                                   # needs palp
            0
        """
        if cone is not self:
            cone = self.embed(cone)
        if '_oriented_boundary' in self.__dict__:
            return self._oriented_boundary[cone]

        # Fix (arbitrary) orientations of the generating cones. Induced
        # by ambient space orientation for full-dimensional cones
        from sage.structure.formal_sum import FormalSum

        def sign(x):
            assert x != 0
            if x > 0:
                return 1
            else:
                return -1
        N_QQ = self.lattice().base_extend(QQ)
        dim = self.lattice_dim()
        outward_vectors = {}
        generating_cones = []
        for c in self.generating_cones():
            if c.dim() == dim:
                outward_v = []
            else:
                Q = N_QQ.quotient(c.rays())
                outward_v = [Q.lift(q) for q in Q.gens()]

            outward_vectors[c] = outward_v
            orientation = sign(matrix(outward_v + list(c.rays().basis())).det())
            generating_cones.append(tuple([orientation, c]))
        boundaries = {self: FormalSum(generating_cones)}

        # The orientation of each facet is arbitrary, but the
        # partition of the boundary in positively and negatively
        # oriented facets is not.
        for d in range(dim, -1, -1):
            for c in self(d):
                c_boundary = []
                c_matrix = matrix(outward_vectors[c] + list(c.rays().basis()))
                c_matrix_inv = c_matrix.inverse()
                for facet in c.facets():
                    outward_ray_indices = set(c.ambient_ray_indices()) \
                              .difference(set(facet.ambient_ray_indices()))
                    outward_vector = - sum(self.ray(i) for i in outward_ray_indices)
                    outward_vectors[facet] = [outward_vector] + outward_vectors[c]
                    facet_matrix = matrix(outward_vectors[facet] + list(facet.rays().basis()))
                    orientation = sign((c_matrix_inv * facet_matrix).det())
                    c_boundary.append(tuple([orientation, facet]))
                boundaries[c] = FormalSum(c_boundary)

        self._oriented_boundary = boundaries
        return boundaries[cone]

    def toric_variety(self, *args, **kwds):
        """
        Return the associated toric variety.

        INPUT:

        Same arguments as :func:`~sage.schemes.toric.variety.ToricVariety`.

        OUTPUT: a toric variety

        This is equivalent to the command ``ToricVariety(self)`` and
        is provided only as a convenient alternative method to go from the
        fan to the associated toric variety.

        EXAMPLES::

            sage: Fan([Cone([(1,0)]), Cone([(0,1)])]).toric_variety()
            2-d toric variety covered by 2 affine patches
        """
        from sage.schemes.toric.variety import ToricVariety
        return ToricVariety(self, *args, **kwds)

    def complex(self, base_ring=ZZ, extended=False):
        r"""
        Return the chain complex of the fan.

        To a `d`-dimensional fan `\Sigma`, one can canonically
        associate a chain complex `K^\bullet`

        .. MATH::

            0 \longrightarrow
            \ZZ^{\Sigma(d)} \longrightarrow
            \ZZ^{\Sigma(d-1)} \longrightarrow
            \cdots \longrightarrow
            \ZZ^{\Sigma(0)} \longrightarrow
            0

        where the leftmost nonzero entry is in degree `0` and the
        rightmost entry in degree `d`. See [Kly1990]_, eq. (3.2). This
        complex computes the homology of `|\Sigma|\subset N_\RR` with
        arbitrary support,

        .. MATH::

            H_i(K) = H_{d-i}(|\Sigma|, \ZZ)_{\text{non-cpct}}

        For a complete fan, this is just the non-compactly supported
        homology of `\RR^d`. In this case, `H_0(K)=\ZZ` and `0` in all
        nonzero degrees.

        For a complete fan, there is an extended chain complex

        .. MATH::

            0 \longrightarrow
            \ZZ \longrightarrow
            \ZZ^{\Sigma(d)} \longrightarrow
            \ZZ^{\Sigma(d-1)} \longrightarrow
            \cdots \longrightarrow
            \ZZ^{\Sigma(0)} \longrightarrow
            0

        where we take the first `\ZZ` term to be in degree -1. This
        complex is an exact sequence, that is, all homology groups
        vanish.

        The orientation of each cone is chosen as in
        :meth:`oriented_boundary`.

        INPUT:

        - ``extended`` -- boolean (default: ``False``); whether to
          construct the extended complex, that is, including the
          `\ZZ`-term at degree -1 or not

        - ``base_ring`` -- a ring (default: ``ZZ``); the ring to use
          instead of `\ZZ`

        OUTPUT:

        The complex associated to the fan as a :class:`ChainComplex
        <sage.homology.chain_complex.ChainComplex>`. This raises a
        :exc:`ValueError` if the extended complex is requested for a
        non-complete fan.

        EXAMPLES::

            sage: # needs palp
            sage: fan = toric_varieties.P(3).fan()
            sage: K_normal = fan.complex(); K_normal
            Chain complex with at most 4 nonzero terms over Integer Ring
            sage: K_normal.homology()
            {0: Z, 1: 0, 2: 0, 3: 0}
            sage: K_extended = fan.complex(extended=True); K_extended
            Chain complex with at most 5 nonzero terms over Integer Ring
            sage: K_extended.homology()
            {-1: 0, 0: 0, 1: 0, 2: 0, 3: 0}

        Homology computations are much faster over `\QQ` if you do not
        care about the torsion coefficients::

            sage: toric_varieties.P2_123().fan().complex(extended=True,                 # needs palp
            ....:                                        base_ring=QQ)
            Chain complex with at most 4 nonzero terms over Rational Field
            sage: _.homology()                                                          # needs palp
            {-1: Vector space of dimension 0 over Rational Field,
             0: Vector space of dimension 0 over Rational Field,
             1: Vector space of dimension 0 over Rational Field,
             2: Vector space of dimension 0 over Rational Field}

        The extended complex is only defined for complete fans::

            sage: fan = Fan([Cone([(1,0)])])
            sage: fan.is_complete()
            False
            sage: fan.complex(extended=True)
            Traceback (most recent call last):
            ...
            ValueError: The extended complex is only defined for complete fans!

        The definition of the complex does not refer to the ambient
        space of the fan, so it does not distinguish a fan from the
        same fan embedded in a subspace::

            sage: K1 = Fan([Cone([(-1,)]), Cone([(1,)])]).complex()
            sage: K2 = Fan([Cone([(-1,0,0)]), Cone([(1,0,0)])]).complex()
            sage: K1 == K2
            True

        Things get more complicated for non-complete fans::

            sage: fan = Fan([Cone([(1,1,1)]),
            ....:            Cone([(1,0,0), (0,1,0)]),
            ....:            Cone([(-1,0,0), (0,-1,0), (0,0,-1)])])
            sage: fan.complex().homology()
            {0: 0, 1: 0, 2: Z x Z, 3: 0}
            sage: fan = Fan([Cone([(1,0,0), (0,1,0)]),
            ....:            Cone([(-1,0,0), (0,-1,0), (0,0,-1)])])
            sage: fan.complex().homology()
            {0: 0, 1: 0, 2: Z, 3: 0}
            sage: fan = Fan([Cone([(-1,0,0), (0,-1,0), (0,0,-1)])])
            sage: fan.complex().homology()
            {0: 0, 1: 0, 2: 0, 3: 0}
        """
        dim = self.dim()
        delta = {}
        for degree in range(1, dim + 1):
            m = matrix(base_ring, len(self(degree - 1)), len(self(degree)),
                       base_ring.zero())
            for i, cone in enumerate(self(degree)):
                boundary = self.oriented_boundary(cone)
                for orientation, d_cone in boundary:
                    m[self(degree - 1).index(d_cone), i] = orientation
            delta[dim - degree] = m

        from sage.homology.chain_complex import ChainComplex
        if not extended:
            return ChainComplex(delta, base_ring=base_ring)

        # add the extra entry for the extended complex
        if not self.is_complete():
            raise ValueError('The extended complex is only defined for complete fans!')
        extension = matrix(base_ring, len(self(dim)), 1, base_ring.zero())
        generating_cones = self.oriented_boundary(self)
        for orientation, d_cone in generating_cones:
            extension[self(dim).index(d_cone), 0] = orientation
        delta[-1] = extension
        return ChainComplex(delta, base_ring=base_ring)


def discard_faces(cones):
    r"""
    Return the cones of the given list which are not faces of each other.

    INPUT:

    - ``cones`` -- list of
      :class:`cones <sage.geometry.cone.ConvexRationalPolyhedralCone>`

    OUTPUT:

    a list of
    :class:`cones <sage.geometry.cone.ConvexRationalPolyhedralCone>`,
    sorted by dimension in decreasing order

    EXAMPLES:

    Consider all cones of a fan::

        sage: Sigma = toric_varieties.P2().fan()                                        # needs palp
        sage: cones = flatten(Sigma.cones())                                            # needs palp
        sage: len(cones)                                                                # needs palp
        7

    Most of them are not necessary to generate this fan::

        sage: from sage.geometry.fan import discard_faces
        sage: len(discard_faces(cones))                                                 # needs palp
        3
        sage: Sigma.ngenerating_cones()                                                 # needs palp
        3
    """
    # Convert to a list or make a copy, so that the input is unchanged.
    cones = list(cones)
    cones.sort(key=lambda cone: cone.dim(), reverse=True)
    generators = []
    for cone in cones:
        if not any(cone.is_face_of(other) for other in generators):
            generators.append(cone)
    return generators


_discard_faces = discard_faces  # Due to a name conflict in Fan constructor


def _refine_arrangement_to_fan(cones):
    """
    Refine the cones of the given list so that they can belong to the same fan.

    INPUT:

    - ``cones`` -- list of rational cones that are possibly overlapping

    OUTPUT: list of refined cones

    EXAMPLES::

        sage: from sage.geometry.fan import _refine_arrangement_to_fan
        sage: c1 = Cone([(-2,-1,1), (-2,1,1), (2,1,1), (2,-1,1)])
        sage: c2 = Cone([(-1,-2,1), (-1,2,1), (1,2,1), (1,-2,1)])
        sage: refined_cones = _refine_arrangement_to_fan([c1, c2])
        sage: for cone in refined_cones: print(cone.rays())
        N(-1,  1, 1),
        N(-1, -1, 1),
        N( 1, -1, 1),
        N( 1,  1, 1)
        in 3-d lattice N
        N(1, -1, 1),
        N(1,  1, 1),
        N(2, -1, 1),
        N(2,  1, 1)
        in 3-d lattice N
        N(-2,  1, 1),
        N(-1, -1, 1),
        N(-1,  1, 1),
        N(-2, -1, 1)
        in 3-d lattice N
        N(-1, 1, 1),
        N(-1, 2, 1),
        N( 1, 1, 1),
        N( 1, 2, 1)
        in 3-d lattice N
        N(-1, -1, 1),
        N(-1, -2, 1),
        N( 1, -2, 1),
        N( 1, -1, 1)
        in 3-d lattice N
    """
    dual_lattice = cones[0].dual_lattice()
    is_face_to_face = True
    for i in range(len(cones)):
        ci = cones[i]
        for j in range(i):
            cj = cones[j]
            c = ci.intersection(cj)
            if not (c.is_face_of(ci)) or not (c.is_face_of(cj)):
                is_face_to_face = False
                break
        if not is_face_to_face:
            break
    if is_face_to_face:
        return cones
    facet_normal_vectors = []
    for c in cones:
        for l in c.polyhedron().Hrepresentation():
            v = l[1::]
            is_new = True
            for fnv in facet_normal_vectors:
                if span([v, fnv]).dimension() < 2:
                    is_new = False
                    break
            if is_new:
                facet_normal_vectors.append(v)
    for v in facet_normal_vectors:
        halfspace1 = Cone([v], lattice=dual_lattice).dual()
        halfspace2 = Cone([-v], lattice=dual_lattice).dual()
        subcones = []
        for c in cones:
            subc1 = c.intersection(halfspace1)
            subc2 = c.intersection(halfspace2)
            for subc in [subc1, subc2]:
                if subc.dim() == c.dim():
                    is_new = True
                    for subcone in subcones:
                        if subc.dim() == subcone.dim() and subc.is_equivalent(subcone):
                            is_new = False
                            break
                    if is_new:
                        subcones.append(subc)
        cones = subcones
    return cones
