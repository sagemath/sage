# sage.doctest: needs sage.geometry.polyhedron sage.graphs
r"""
Morphisms of toric varieties

There are three "obvious" ways to map toric varieties to toric
varieties:

1. Polynomial maps in local coordinates, the usual morphisms in
   algebraic geometry.

2. Polynomial maps in the (global) homogeneous coordinates.

3. Toric morphisms, that is, algebraic morphisms equivariant with
   respect to the torus action on the toric variety.

Both 2 and 3 are special cases of 1, which is just to say that we
always remain within the realm of algebraic geometry. But apart from
that, none is included in one of the other cases. In the examples
below, we will explore some algebraic maps that can or can not be
written as a toric morphism. Often a toric morphism can be written
with polynomial maps in homogeneous coordinates, but sometimes it
cannot.

The toric morphisms are perhaps the most mysterious at the
beginning. Let us quickly review their definition (See Definition
3.3.3 of [CLS2011]_). Let `\Sigma_1` be a fan in `N_{1,\RR}` and `\Sigma_2` be a
fan in `N_{2,\RR}`. A morphism `\phi: X_{\Sigma_1} \to X_{\Sigma_2}`
of the associated toric varieties is toric if `\phi` maps the maximal
torus `T_{N_1} \subseteq X_{\Sigma_1}` into `T_{N_2} \subseteq
X_{\Sigma_2}` and `\phi|_{T_N}` is a group homomorphism.

The data defining a toric morphism is precisely what defines a fan
morphism (see :mod:`~sage.geometry.fan_morphism`), extending the more
familiar dictionary between toric varieties and fans. Toric geometry
is a functor from the category of fans and fan morphisms to the
category of toric varieties and toric morphisms.

.. NOTE::

    Do not create the toric morphisms (or any morphism of schemes)
    directly from the ``SchemeMorphism...`` classes. Instead, use the
    :meth:`~sage.schemes.generic.scheme.hom` method common to all
    algebraic schemes to create new homomorphisms.

EXAMPLES:

First, consider the following embedding of `\mathbb{P}^1` into
`\mathbb{P}^2` ::

    sage: P2.<x,y,z> = toric_varieties.P2()
    sage: P1.<u,v> = toric_varieties.P1()
    sage: P1.hom([0, u^2 + v^2, u*v], P2)
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending [u : v] to [0 : u^2 + v^2 : u*v]

This is a well-defined morphism of algebraic varieties because
homogeneously rescaled coordinates of a point of `\mathbb{P}^1` map to the same
point in `\mathbb{P}^2` up to its homogeneous rescalings. It is not
equivariant with respect to the torus actions

.. MATH::

    \CC^\times \times \mathbb{P}^1,
    (\mu,[u:v]) \mapsto [u:\mu v]
    \quad\text{and}\quad
    \left(\CC^\times\right)^2 \times \mathbb{P}^2,
    ((\alpha,\beta),[x:y:z]) \mapsto [x:\alpha y:\beta z]
    ,

hence it is not a toric morphism. Clearly, the problem is that
the map in homogeneous coordinates contains summands that transform
differently under the torus action. However, this is not the only
difficulty. For example, consider ::

    sage: phi = P1.hom([0,u,v], P2);  phi
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending [u : v] to [0 : u : v]

This map is actually the embedding of the
:meth:`~sage.schemes.toric.variety.ToricVariety_field.orbit_closure`
associated to one of the rays of the fan of `\mathbb{P}^2`. Now the
morphism is equivariant with respect to **some** map `\CC^\times \to
(\CC^\times)^2` of the maximal tori of `\mathbb{P}^1` and
`\mathbb{P}^2`. But this map of the maximal tori cannot be the same as
``phi`` defined above. Indeed, the image of ``phi`` completely misses
the maximal torus `T_{\mathbb{P}^2} = \{ [x:y:z] | x\not=0, y\not=0,
z\not=0 \}` of `\mathbb{P}^2`.

Consider instead the following morphism of fans::

    sage: fm = FanMorphism(matrix(ZZ, [[1,0]]), P1.fan(), P2.fan());  fm
    Fan morphism defined by the matrix
    [1 0]
    Domain fan: Rational polyhedral fan in 1-d lattice N
    Codomain fan: Rational polyhedral fan in 2-d lattice N

which also defines a morphism of toric varieties::

    sage: P1.hom(fm, P2)
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined by sending Rational polyhedral fan in 1-d lattice N
            to Rational polyhedral fan in 2-d lattice N.

The fan morphism map is equivalent to the following polynomial map::

    sage: _.as_polynomial_map()
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending [u : v] to [u : v : v]

Finally, here is an example of a fan morphism that cannot be written
using homogeneous polynomials. Consider the blowup `O_{\mathbb{P}^1}(2)
\to \CC^2/\ZZ_2`. In terms of toric data, this blowup is::

    sage: A2_Z2 = toric_varieties.A2_Z2()
    sage: A2_Z2.fan().rays()
    N(1, 0),
    N(1, 2)
    in 2-d lattice N
    sage: O2_P1 = A2_Z2.resolve(new_rays=[(1,1)])
    sage: blowup = O2_P1.hom(identity_matrix(2), A2_Z2)
    sage: blowup.as_polynomial_map()
    Traceback (most recent call last):
    ...
    TypeError: the fan morphism cannot be written in homogeneous polynomials

If we denote the homogeneous coordinates of `O_{\mathbb{P}^1}(2)` by
`x`, `t`, `y` corresponding to the rays `(1,2)`, `(1,1)`, and `(1,0)`
then the blow-up map is [BB2013]_:

.. MATH::

    f: O_{\mathbb{P}^1}(2) \to \CC^2/\ZZ_2, \quad
    (x,t,y) \mapsto \left( x\sqrt{t}, y\sqrt{t} \right)

which requires square roots.


Fibrations
----------

If a toric morphism is :meth:`dominant
<SchemeMorphism_fan_toric_variety.is_dominant>`, then all fibers over
a fixed torus orbit in the base are isomorphic. Hence, studying the
fibers is again a combinatorial question and Sage implements
additional methods to study such fibrations that are not available
otherwise (however, note that you can always
:meth:`~SchemeMorphism_fan_toric_variety.factor` to pick out the part
that is dominant over the image or its closure).

For example, consider the blow-up restricted to one of the two
coordinate charts of `O_{\mathbb{P}^1}(2)` ::


    sage: O2_P1_chart = ToricVariety(Fan([O2_P1.fan().generating_cones()[0]]))
    sage: single_chart = O2_P1_chart.hom(identity_matrix(2), A2_Z2)
    sage: single_chart.is_dominant()
    True
    sage: single_chart.is_surjective()
    False

    sage: fiber = single_chart.fiber_generic();  fiber
    (0-d affine toric variety, 1)
    sage: fiber[0].embedding_morphism().as_polynomial_map()
    Scheme morphism:
      From: 0-d affine toric variety
      To:   2-d affine toric variety
      Defn: Defined on coordinates by sending [] to [1 : 1]

The fibers are labeled by torus orbits in the base, that is, cones of
the codomain fan. In this case, the fibers over lower-dimensional
torus orbits are::

    sage: A2_Z2_cones = flatten(A2_Z2.fan().cones())
    sage: table([('cone', 'dim')] +
    ....:       [(cone.ambient_ray_indices(), single_chart.fiber_dimension(cone))
    ....:        for cone in A2_Z2_cones], header_row=True)
        cone     dim
      ├────────┼─────┤
        ()       0
        (0,)     0
        (1,)     -1
        (0, 1)   1

Lets look closer at the one-dimensional fiber. Although not the case
in this example, connected components of fibers over higher-dimensional cones
(corresponding
to lower-dimensional torus orbits) of the base are often not
irreducible. The irreducible components are labeled by the
:meth:`~sage.geometry.fan_morphism.FanMorphism.primitive_preimage_cones`,
which are certain cones of the domain fan that map to the cone in the
base that defines the torus orbit::

    sage: table([('base cone', 'primitive preimage cones')] +
    ....:       [(cone.ambient_ray_indices(),
    ....:         single_chart.fan_morphism().primitive_preimage_cones(cone))
    ....:        for cone in A2_Z2_cones], header_row=True)
        base cone   primitive preimage cones
      ├───────────┼─────────────────────────────────────────────────────────┤
        ()          (0-d cone of Rational polyhedral fan in 2-d lattice N,)
        (0,)        (1-d cone of Rational polyhedral fan in 2-d lattice N,)
        (1,)        ()
        (0, 1)      (1-d cone of Rational polyhedral fan in 2-d lattice N,)

The fiber over the trivial cone is the generic fiber that we have
already encountered. The interesting fiber is the one over the
2-dimensional cone, which represents the exceptional set of the
blow-up in this single coordinate chart. Lets investigate further::

    sage: fm = single_chart.fan_morphism()
    sage: exceptional_cones = fm.primitive_preimage_cones(A2_Z2.fan(2)[0])
    sage: exceptional_set = single_chart.fiber_component(exceptional_cones[0])
    sage: exceptional_set
    1-d affine toric variety
    sage: exceptional_set.embedding_morphism().as_polynomial_map()
    Scheme morphism:
      From: 1-d affine toric variety
      To:   2-d affine toric variety
      Defn: Defined on coordinates by sending [z0] to [z0 : 0]

So we see that the fiber over this point is an affine line. Together
with another affine line in the other coordinate patch, this covers
the exceptional `\mathbb{P}^1` of the blowup `O_{\mathbb{P}^1}(2) \to
\CC^2/\ZZ_2`.

Here is an example with higher dimensional varieties involved::

    sage: A3 = toric_varieties.A(3)
    sage: P3 = toric_varieties.P(3)
    sage: m = matrix([(2,0,0), (1,1,0), (3,1,0)])
    sage: phi = A3.hom(m, P3)
    sage: phi.as_polynomial_map()
    Scheme morphism:
      From: 3-d affine toric variety
      To:   3-d CPR-Fano toric variety covered by 4 affine patches
      Defn: Defined on coordinates by sending [z0 : z1 : z2] to
            [z0^2*z1*z2^3 : z1*z2 : 1 : 1]
    sage: phi.fiber_generic()
    Traceback (most recent call last):
    ...
    AttributeError: 'SchemeMorphism_fan_toric_variety' object
    has no attribute 'fiber_generic'...

Let's use factorization mentioned above::

    sage: phi_i, phi_b, phi_s = phi.factor()

It is possible to study fibers of the last two morphisms or their composition::

    sage: phi_d = phi_b * phi_s
    sage: phi_d
    Scheme morphism:
      From: 3-d affine toric variety
      To:   2-d toric variety covered by 3 affine patches
      Defn: Defined by sending Rational polyhedral fan in 3-d lattice N to
            Rational polyhedral fan in Sublattice <N(1, 0, 0), N(0, 1, 0)>.
    sage: phi_d.as_polynomial_map()
    Scheme morphism:
      From: 3-d affine toric variety
      To:   2-d toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending [z0 : z1 : z2] to
            [1 : z1*z2 : z0^2*z1*z2^3]
    sage: phi_d.codomain().fan().rays()
    N(-1, -1, 0),
    N( 0,  1, 0),
    N( 1,  0, 0)
    in Sublattice <N(1, 0, 0), N(0, 1, 0)>
    sage: for c in phi_d.codomain().fan():
    ....:     c.ambient_ray_indices()
    (0, 1)
    (0, 2)
    (1, 2)

We see that codomain fan of this morphism is a projective plane, which can be
verified by ::

    sage: phi_d.codomain().fan().is_isomorphic(toric_varieties.P2().fan()) # known bug
    True

(Unfortunately it cannot be verified correctly until :issue:`16012` is fixed.)

We now have access to fiber methods::

    sage: fiber = phi_d.fiber_generic()
    sage: fiber
    (1-d affine toric variety, 2)
    sage: fiber[0].embedding_morphism()
    Scheme morphism:
      From: 1-d affine toric variety
      To:   3-d affine toric variety
      Defn: Defined by sending
            Rational polyhedral fan in Sublattice <N(1, 1, -1)> to
            Rational polyhedral fan in 3-d lattice N.
    sage: fiber[0].embedding_morphism().as_polynomial_map()
    Traceback (most recent call last):
    ...
    NotImplementedError: polynomial representations for
    fans with virtual rays are not implemented yet
    sage: fiber[0].fan().rays()
    Empty collection
    in Sublattice <N(1, 1, -1)>

We see that generic fibers of this morphism consist of 2 one-dimensional tori
each. To see what happens over boundary points we can look at fiber components
corresponding to the cones of the domain fan::

    sage: fm = phi_d.fan_morphism()
    sage: for c in flatten(phi_d.domain().fan().cones()):
    ....:     fc, m = phi_d.fiber_component(c, multiplicity=True)
    ....:     print("{} |-> {} ({} rays, multiplicity {}) over {}".format(
    ....:       c.ambient_ray_indices(), fc, fc.fan().nrays(),
    ....:       m, fm.image_cone(c).ambient_ray_indices()))
    () |-> 1-d affine toric variety (0 rays, multiplicity 2) over ()
    (0,) |-> 1-d affine toric variety (0 rays, multiplicity 1) over (2,)
    (1,) |-> 2-d affine toric variety (2 rays, multiplicity 1) over (1, 2)
    (2,) |-> 2-d affine toric variety (2 rays, multiplicity 1) over (1, 2)
    (0, 1) |-> 1-d affine toric variety (1 rays, multiplicity 1) over (1, 2)
    (1, 2) |-> 1-d affine toric variety (1 rays, multiplicity 1) over (1, 2)
    (0, 2) |-> 1-d affine toric variety (1 rays, multiplicity 1) over (1, 2)
    (0, 1, 2) |-> 0-d affine toric variety (0 rays, multiplicity 1) over (1, 2)

Now we see that over one of the coordinate lines of the projective plane we also
have one-dimensional tori (but only one in each fiber), while over one of the
points fixed by torus action we have two affine planes intersecting along an
affine line. An alternative perspective is provided by cones of the codomain
fan::

    sage: for c in flatten(phi_d.codomain().fan().cones()):
    ....:     print("{} connected components over {}, each with {} irreducible components.".format(
    ....:       fm.index(c), c.ambient_ray_indices(),
    ....:       len(fm.primitive_preimage_cones(c))))
    2 connected components over (), each with 1 irreducible components.
    None connected components over (0,), each with 0 irreducible components.
    None connected components over (1,), each with 0 irreducible components.
    1 connected components over (2,), each with 1 irreducible components.
    None connected components over (0, 1), each with 0 irreducible components.
    None connected components over (0, 2), each with 0 irreducible components.
    1 connected components over (1, 2), each with 2 irreducible components.
"""

# ****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# For now, the scheme morphism base class cannot derive from Morphism
# since this would clash with elliptic curves. So we derive only on
# the toric varieties level from Morphism. See
# https://groups.google.com/d/msg/sage-devel/qF4yU6Vdmao/wQlNrneSmWAJ
from sage.categories.morphism import Morphism
from sage.structure.richcmp import richcmp_not_equal, richcmp
from sage.structure.sequence import Sequence
from sage.rings.integer_ring import ZZ
from sage.arith.misc import GCD as gcd
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix, identity_matrix
from sage.modules.free_module_element import vector
from sage.geometry.cone import Cone
from sage.geometry.fan import Fan

from sage.schemes.generic.scheme import Scheme
from sage.schemes.generic.morphism import (
    SchemeMorphism, SchemeMorphism_point, SchemeMorphism_polynomial
)


############################################################################
# A points on a toric variety determined by homogeneous coordinates.
class SchemeMorphism_point_toric_field(SchemeMorphism_point, Morphism):
    """
    A point of a toric variety determined by homogeneous coordinates
    in a field.

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.generic.scheme.hom` method of
        :class:`toric varieties
        <sage.schemes.toric.variety.ToricVariety_field>`
        instead.

    INPUT:

    - ``X`` -- toric variety or subscheme of a toric variety

    - ``coordinates`` -- list of coordinates in the base field of ``X``

    - ``check`` -- if ``True`` (default), the input will be checked for
      correctness

    OUTPUT: a :class:`SchemeMorphism_point_toric_field`

    TESTS::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1(1,2,3,4)
        [1 : 2 : 3 : 4]
    """
    # Mimicking affine/projective classes
    def __init__(self, X, coordinates, check=True):
        r"""
        See :class:`SchemeMorphism_point_toric_field` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1(1,2,3,4)
            [1 : 2 : 3 : 4]
        """
        # Convert scheme to its set of points over the base ring
        if isinstance(X, Scheme):
            X = X(X.base_ring())
        super().__init__(X)
        if check:
            # Verify that there are the right number of coords
            # Why is it not done in the parent?
            if isinstance(coordinates, SchemeMorphism):
                coordinates = list(coordinates)
            if not isinstance(coordinates, (list, tuple)):
                raise TypeError("coordinates must be a scheme point, list, "
                                "or tuple; got %s" % coordinates)
            d = X.codomain().ambient_space().ngens()
            if len(coordinates) != d:
                raise ValueError("there must be %d coordinates; got only %d: "
                                 "%s" % (d, len(coordinates), coordinates))
            # Make sure the coordinates all lie in the appropriate ring
            coordinates = Sequence(coordinates, X.value_ring())
            # Verify that the point satisfies the equations of X.
            X.codomain()._check_satisfies_equations(coordinates)
        self._coords = coordinates


############################################################################
# A morphism of toric varieties determined by homogeneous polynomials.
class SchemeMorphism_polynomial_toric_variety(SchemeMorphism_polynomial, Morphism):
    """
    A morphism determined by homogeneous polynomials.

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.generic.scheme.hom` method of
        :class:`toric varieties
        <sage.schemes.toric.variety.ToricVariety_field>`
        instead.

    INPUT:

    Same as for
    :class:`~sage.schemes.toric.morphism.SchemeMorphism_polynomial`.

    OUTPUT: a :class:`~sage.schemes.toric.morphism.SchemeMorphism_polynomial_toric_variety`

    TESTS::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1.inject_variables()
        Defining s, t, x, y
        sage: P1 = P1xP1.subscheme(s - t)
        sage: H = P1xP1.Hom(P1)
        sage: import sage.schemes.toric.morphism as MOR
        sage: MOR.SchemeMorphism_polynomial_toric_variety(H, [s, s, x, y])
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   Closed subscheme of 2-d CPR-Fano toric variety
                covered by 4 affine patches defined by:
          s - t
          Defn: Defined on coordinates by sending [s : t : x : y] to
                [s : s : x : y]
    """

    def __init__(self, parent, polynomials, check=True):
        r"""
        See :class:`SchemeMorphism_polynomial_toric_variety` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: P1 = P1xP1.subscheme(s - t)
            sage: H = P1xP1.Hom(P1)
            sage: import sage.schemes.toric.morphism as MOR
            sage: MOR.SchemeMorphism_polynomial_toric_variety(H, [s, s, x, y])
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   Closed subscheme of 2-d CPR-Fano toric variety
                    covered by 4 affine patches defined by:
              s - t
              Defn: Defined on coordinates by sending [s : t : x : y] to
                    [s : s : x : y]
        """
        SchemeMorphism_polynomial.__init__(self, parent, polynomials, check)
        if check:
            # Check that defining polynomials are homogeneous (degrees can be
            # different if the target uses weighted coordinates)
            for p in self.defining_polynomials():
                if not self.domain().ambient_space().is_homogeneous(p):
                    raise ValueError("%s is not homogeneous" % p)

    def as_fan_morphism(self):
        """
        Express the morphism as a map defined by a fan morphism.

        OUTPUT: a :class:`SchemeMorphism_polynomial_toric_variety`

        This raises a :exc:`TypeError` if the morphism cannot be written
        in such a way.

        EXAMPLES::

            sage: A1.<z> = toric_varieties.A1()
            sage: P1 = toric_varieties.P1()
            sage: patch = A1.hom([1,z], P1)
            sage: patch.as_fan_morphism()
            Traceback (most recent call last):
            ...
            NotImplementedError: expressing toric morphisms as fan morphisms is
            not implemented yet
        """
        raise NotImplementedError("expressing toric morphisms as fan "
                                  "morphisms is not implemented yet")


############################################################################
# The embedding morphism of an orbit closure
class SchemeMorphism_orbit_closure_toric_variety(SchemeMorphism, Morphism):
    """
    The embedding of an orbit closure.

    INPUT:

    - ``parent`` -- the parent homset

    - ``defining_cone`` -- the defining cone

    - ``ray_map`` -- dictionary ``{ambient ray generator: orbit ray
      generator}``. Note that the image of the ambient ray generator
      is not necessarily primitive.

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.toric.variety.ToricVariety_field.orbit_closure`
        method of :class:`toric varieties
        <sage.schemes.toric.variety.ToricVariety_field>`
        instead.

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: H = P1xP1.fan(1)[0]
        sage: V = P1xP1.orbit_closure(H)
        sage: V.embedding_morphism()
        Scheme morphism:
          From: 1-d toric variety covered by 2 affine patches
          To:   2-d CPR-Fano toric variety covered by 4 affine patches
          Defn: Defined by embedding the torus closure associated to the 1-d
                cone of Rational polyhedral fan in 2-d lattice N.

    TESTS::

        sage: V.embedding_morphism()._reverse_ray_map()
        {N(-1): 3, N(1): 2}
        sage: V.embedding_morphism()._defining_cone
        1-d cone of Rational polyhedral fan in 2-d lattice N
    """
    def __init__(self, parent, defining_cone, ray_map):
        """
        The Python constructor.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: P1 = P2.orbit_closure(P2.fan(1)[0])
            sage: P1.embedding_morphism()
            Scheme morphism:
              From: 1-d toric variety covered by 2 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by embedding the torus closure associated to the 1-d cone
                    of Rational polyhedral fan in 2-d lattice N.
        """
        SchemeMorphism.__init__(self, parent)
        self._defining_cone = defining_cone
        self._ray_map = ray_map

    def defining_cone(self):
        r"""
        Return the cone corresponding to the torus orbit.

        OUTPUT: a cone of the fan of the ambient toric variety

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: cone = P2.fan(1)[0]
            sage: P1 = P2.orbit_closure(cone)
            sage: P1.embedding_morphism().defining_cone()
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: _ is cone
            True
        """
        return self._defining_cone

    @cached_method
    def _reverse_ray_map(self):
        """
        Reverse ``self._ray_map``.

        OUTPUT:

        A dictionary ``{orbit ray generator: preimage ray
        index}``. Note that the orbit ray generator need not be
        primitive. Also, the preimage ray is not necessarily unique.

        EXAMPLES::

            sage: P2_112 = toric_varieties.P2_112()
            sage: P1 = P2_112.orbit_closure(Cone([(1,0)]))
            sage: f = P1.embedding_morphism()
            sage: f._ray_map
            {N(-1, -2): (-2), N(0, 1): (1), N(1, 0): (0)}
            sage: f._reverse_ray_map()
            {N(-2): 2, N(1): 1}
        """
        orbit = self.parent().domain()
        codomain_fan = self.parent().codomain().fan()
        reverse_ray_dict = {}
        for n1, n2 in self._ray_map.items():
            ray_index = codomain_fan.rays().index(n1)
            if n2.is_zero():
                assert ray_index in self._defining_cone.ambient_ray_indices()
                continue
            n2 = orbit.fan().lattice()(n2)
            n2.set_immutable()
            reverse_ray_dict[n2] = ray_index
        return reverse_ray_dict

    def _repr_defn(self):
        """
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: V = P2.orbit_closure(P2.fan(1)[0]);  V
            1-d toric variety covered by 2 affine patches
            sage: V.embedding_morphism()._repr_defn()
            'Defined by embedding the torus closure associated to the 1-d cone of
             Rational polyhedral fan in 2-d lattice N.'
        """
        s = 'Defined by embedding the torus closure associated to the '
        s += str(self._defining_cone)
        s += '.'
        return s

    def as_polynomial_map(self):
        """
        Express the morphism via homogeneous polynomials.

        OUTPUT: a :class:`SchemeMorphism_polynomial_toric_variety`

        This raises a :exc:`TypeError` if the morphism cannot be
        written in terms of homogeneous polynomials.

        The defining polynomials are not necessarily unique. There are
        choices if multiple ambient space ray generators project to
        the same orbit ray generator, and one such choice is made
        implicitly. The orbit embedding can be written as a polynomial
        map if and only if each primitive orbit ray generator is the
        image of at least one primitive ray generator of the ambient
        toric variety.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: V = P2.orbit_closure(P2.fan(1)[0]);  V
            1-d toric variety covered by 2 affine patches
            sage: V.embedding_morphism().as_polynomial_map()
            Scheme morphism:
              From: 1-d toric variety covered by 2 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined on coordinates by sending [z0 : z1] to [0 : z1 : z0]

        If the toric variety is singular, then some orbit closure
        embeddings cannot be written with homogeneous polynomials::

            sage: P2_112 = toric_varieties.P2_112()
            sage: P1 = P2_112.orbit_closure(Cone([(1,0)]))
            sage: P1.embedding_morphism().as_polynomial_map()
            Traceback (most recent call last):
            ...
            TypeError: the embedding cannot be written with homogeneous polynomials
        """
        orbit = self.domain()
        codomain_fan = self.codomain().fan()
        R = orbit.coordinate_ring()
        polys = [R.one()] * codomain_fan.nrays()
        for i in self._defining_cone.ambient_ray_indices():
            polys[i] = R.zero()
        ray_index_map = self._reverse_ray_map()
        for i, ray in enumerate(orbit.fan().rays()):
            try:
                ray_index = ray_index_map[ray]
            except KeyError:
                raise TypeError('the embedding cannot be written with homogeneous polynomials')
            polys[ray_index] = R.gen(i)
        return SchemeMorphism_polynomial_toric_variety(self.parent(), polys)

    def pullback_divisor(self, divisor):
        r"""
        Pull back a toric divisor.

        INPUT:

        - ``divisor`` -- a torus-invariant `\QQ`-Cartier divisor on the
          codomain of the embedding map

        OUTPUT:

        A divisor on the domain of the embedding map (the orbit
        closure) that is isomorphic to the pull-back divisor `f^*(D)`
        but with possibly different linearization.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: P1 = P2.orbit_closure(P2.fan(1)[0])
            sage: f = P1.embedding_morphism()
            sage: D = P2.divisor([1,2,3]); D
            V(x) + 2*V(y) + 3*V(z)
            sage: f.pullback_divisor(D)
            4*V(z0) + 2*V(z1)
        """
        from sage.schemes.toric.divisor import ToricDivisor_generic
        if not (isinstance(divisor, ToricDivisor_generic) and divisor.is_QQ_Cartier()):
            raise ValueError('the divisor must be torus-invariant and QQ-Cartier')
        m = divisor.m(self._defining_cone)
        values = []
        codomain_rays = self.codomain().fan().rays()
        for ray in self.domain().fan().rays():
            ray = codomain_rays[self._reverse_ray_map()[ray]]
            value = divisor.function_value(ray) - m * ray
            values.append(value)
        return self.domain().divisor(values)


############################################################################
# A morphism of toric varieties determined by a fan morphism
class SchemeMorphism_fan_toric_variety(SchemeMorphism, Morphism):
    """
    Construct a morphism determined by a fan morphism.

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.generic.scheme.hom` method of
        :class:`toric varieties
        <sage.schemes.toric.variety.ToricVariety_field>`
        instead.

    INPUT:

    - ``parent`` -- Hom-set whose domain and codomain are toric varieties

    - ``fan_morphism`` -- a morphism of fans whose domain and codomain
      fans equal the fans of the domain and codomain in the ``parent``
      Hom-set.

    - ``check`` -- boolean (default: ``True``); whether to
      check the input for consistency

    .. WARNING::

        A fibration is a dominant morphism; if you are interested in
        these then you have to make sure that your fan morphism is
        dominant. For example, this can be achieved by
        :meth:`factoring the morphism
        <sage.schemes.toric.morphism.SchemeMorphism_fan_toric_variety.factor>`. See
        :class:`SchemeMorphism_fan_toric_variety_dominant` for
        additional functionality for fibrations.

    OUTPUT: a :class:`~sage.schemes.toric.morphism.SchemeMorphism_fan_toric_variety`

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1 = toric_varieties.P1()
        sage: f = P1.hom(matrix([[1,0]]), P1xP1);  f
        Scheme morphism:
          From: 1-d CPR-Fano toric variety covered by 2 affine patches
          To:   2-d CPR-Fano toric variety covered by 4 affine patches
          Defn: Defined by sending Rational polyhedral fan in 1-d lattice N
                to Rational polyhedral fan in 2-d lattice N.
        sage: type(f)
        <class 'sage.schemes.toric.morphism.SchemeMorphism_fan_toric_variety'>

    Slightly more explicit construction::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1 = toric_varieties.P1()
        sage: hom_set = P1xP1.Hom(P1)
        sage: fm = FanMorphism(matrix(ZZ, [[1],[0]]), P1xP1.fan(), P1.fan())
        sage: hom_set(fm)
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
          Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                to Rational polyhedral fan in 1-d lattice N.

        sage: P1xP1.hom(fm, P1)
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
          Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                to Rational polyhedral fan in 1-d lattice N.
    """

    def __init__(self, parent, fan_morphism, check=True):
        r"""
        See :class:`SchemeMorphism_polynomial_toric_variety` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: hom_set = P1xP1.Hom(P1)
            sage: fan_morphism = FanMorphism(matrix(ZZ, [[1],[0]]), P1xP1.fan(), P1.fan())
            sage: from sage.schemes.toric.morphism import SchemeMorphism_fan_toric_variety
            sage: SchemeMorphism_fan_toric_variety(hom_set, fan_morphism)
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 1-d lattice N.
        """
        SchemeMorphism.__init__(self, parent)
        if check and self.domain().fan() != fan_morphism.domain_fan():
            raise ValueError('the fan morphism domain must be the fan of the domain')
        if check and self.codomain().fan() != fan_morphism.codomain_fan():
            raise ValueError('the fan morphism codomain must be the fan of the codomain')
        self._fan_morphism = fan_morphism

    def _richcmp_(self, right, op):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- another toric morphism

        OUTPUT: boolean

        Comparison is done first by domain, then by codomain, then by
        fan morphism.

        TESTS::

            sage: A2 = toric_varieties.A2()
            sage: P3 = toric_varieties.P(3)
            sage: m = matrix([(2,0,0), (1,1,0)])
            sage: phi = A2.hom(m, P3)
            sage: phi == phi
            True
            sage: phi == prod(phi.factor())
            True
            sage: phi == phi.factor()[0]
            False
        """
        if not isinstance(right, SchemeMorphism_fan_toric_variety):
            return NotImplemented

        lx = self.domain()
        rx = right.domain()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.codomain()
        rx = right.codomain()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        return richcmp(self.fan_morphism(), right.fan_morphism(), op)

    def _composition_(self, right, homset):
        """
        Return the composition of ``self`` and ``right``.

        INPUT:

        - ``right`` -- a toric morphism defined by a fan morphism

        OUTPUT: a toric morphism

        EXAMPLES::

            sage: A2 = toric_varieties.A2()
            sage: P3 = toric_varieties.P(3)
            sage: m = matrix([(2,0,0), (1,1,0)])
            sage: phi = A2.hom(m, P3)
            sage: phi1, phi2, phi3 = phi.factor()
            sage: phi1 * phi2
            Scheme morphism:
              From: 2-d affine toric variety
              To:   3-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined by sending Rational polyhedral fan in Sublattice
                    <N(1, 0, 0), N(0, 1, 0)> to Rational polyhedral fan in 3-d lattice N.
            sage: phi1 * phi2 * phi3 == phi
            True
        """
        f = self.fan_morphism() * right.fan_morphism()
        return homset(f, self.codomain())

    def _repr_defn(self):
        """
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: f = P1xP1.hom(matrix([[1],[0]]), P1)
            sage: f._repr_defn()
            'Defined by sending Rational polyhedral fan in 2-d lattice N to Rational polyhedral fan in 1-d lattice N.'
        """
        s = 'Defined by sending '
        s += str(self.domain().fan())
        s += ' to '
        s += str(self.codomain().fan())
        s += '.'
        return s

    def factor(self):
        r"""
        Factor ``self`` into injective * birational * surjective morphisms.

        OUTPUT:

        - a triple of toric morphisms `(\phi_i, \phi_b, \phi_s)`, such that
          `\phi_s` is surjective, `\phi_b` is birational, `\phi_i` is injective,
          and ``self`` is equal to `\phi_i \circ \phi_b \circ \phi_s`.

        The intermediate varieties are universal in the following sense. Let
        ``self`` map `X` to `X'` and let `X_s`, `X_i` sit in between, that is,

        .. MATH::

            X
            \twoheadrightarrow
            X_s
            \to
            X_i
            \hookrightarrow
            X'.

        Then any toric morphism from `X` coinciding with ``self`` on the maximal
        torus factors through `X_s` and any toric morphism into `X'` coinciding
        with ``self`` on the maximal torus factors through `X_i`. In particular,
        `X_i` is the closure of the image of ``self`` in `X'`.

        See
        :meth:`~sage.geometry.fan_morphism.FanMorphism.factor`
        for a description of the toric algorithm.

        EXAMPLES:

        We map an affine plane into a projective 3-space in such a way, that it
        becomes "a double cover of a chart of the blow up of one of the
        coordinate planes"::

            sage: A2 = toric_varieties.A2()
            sage: P3 = toric_varieties.P(3)
            sage: m = matrix([(2,0,0), (1,1,0)])
            sage: phi = A2.hom(m, P3)
            sage: phi.as_polynomial_map()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   3-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [x : y] to
                    [x^2*y : y : 1 : 1]

            sage: phi.is_surjective(), phi.is_birational(), phi.is_injective()
            (False, False, False)
            sage: phi_i, phi_b, phi_s = phi.factor()
            sage: phi_s.is_surjective(), phi_b.is_birational(), phi_i.is_injective()
            (True, True, True)
            sage: prod(phi.factor()) == phi
            True

        Double cover (surjective)::

            sage: phi_s.as_polynomial_map()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d affine toric variety
              Defn: Defined on coordinates by sending [x : y] to [x^2 : y]

        Blowup chart (birational)::

            sage: phi_b.as_polynomial_map()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d toric variety covered by 3 affine patches
              Defn: Defined on coordinates by sending [z0 : z1] to [1 : z1 : z0*z1]

        Coordinate plane inclusion (injective)::

            sage: phi_i.as_polynomial_map()
            Scheme morphism:
              From: 2-d toric variety covered by 3 affine patches
              To:   3-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [z0 : z1 : z2] to [z2 : z1 : z0 : z0]
        """
        phi_i, phi_b, phi_s = self.fan_morphism().factor()
        from sage.schemes.toric.variety import ToricVariety
        X = self.domain()
        X_s = ToricVariety(phi_s.codomain_fan())
        X_i = ToricVariety(phi_i.domain_fan())
        X_prime = self.codomain()
        return X_i.hom(phi_i, X_prime), X_s.hom(phi_b, X_i), X.hom(phi_s, X_s)

    def fan_morphism(self):
        """
        Return the defining fan morphism.

        OUTPUT: a :class:`~sage.geometry.fan_morphism.FanMorphism`

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: f = P1xP1.hom(matrix([[1],[0]]), P1)
            sage: f.fan_morphism()
            Fan morphism defined by the matrix
            [1]
            [0]
            Domain fan: Rational polyhedral fan in 2-d lattice N
            Codomain fan: Rational polyhedral fan in 1-d lattice N
        """
        return self._fan_morphism

    def as_polynomial_map(self):
        """
        Express the morphism via homogeneous polynomials.

        OUTPUT: a :class:`SchemeMorphism_polynomial_toric_variety`

        This raises a :exc:`TypeError` if the morphism cannot be written
        in terms of homogeneous polynomials.

        EXAMPLES::

            sage: A1 = toric_varieties.A1()
            sage: square = A1.hom(matrix([[2]]), A1)
            sage: square.as_polynomial_map()
            Scheme endomorphism of 1-d affine toric variety
              Defn: Defined on coordinates by sending [z] to [z^2]

            sage: P1 = toric_varieties.P1()
            sage: patch = A1.hom(matrix([[1]]), P1)
            sage: patch.as_polynomial_map()
            Scheme morphism:
              From: 1-d affine toric variety
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
              Defn: Defined on coordinates by sending [z] to [z : 1]
        """
        R = self.domain().coordinate_ring()
        phi = self.fan_morphism()
        polys = [R.one()] * self.codomain().ngens()
        for rho, x in zip(phi.domain_fan(1), R.gens()):
            ray = rho.ray(0)
            sigma = phi.image_cone(rho)
            degrees = sigma.rays().matrix().solve_left(phi(ray))
            for i, d in zip(sigma.ambient_ray_indices(), degrees):
                try:
                    d = ZZ(d)
                except TypeError:
                    raise TypeError('the fan morphism cannot be written in '
                                    'homogeneous polynomials')
                polys[i] *= x**d
        if phi.domain_fan().virtual_rays():
            raise NotImplementedError("polynomial representations for fans with"
                                      " virtual rays are not implemented yet")
        return SchemeMorphism_polynomial_toric_variety(self.parent(), polys)

    def is_bundle(self):
        r"""
        Check if ``self`` is a bundle.

        See :meth:`~sage.geometry.fan_morphism.FanMorphism.is_bundle`
        for fan morphisms for details.

        OUTPUT: ``True`` if ``self`` is a bundle, ``False`` otherwise

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: P1xP1.hom(matrix([[1],[0]]), P1).is_bundle()
            True
        """
        return self.fan_morphism().is_bundle()

    def is_fibration(self):
        r"""
        Check if ``self`` is a fibration.

        See
        :meth:`~sage.geometry.fan_morphism.FanMorphism.is_fibration`
        for fan morphisms for details.

        OUTPUT: ``True`` if ``self`` is a fibration, ``False`` otherwise

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: P1xP1.hom(matrix([[1],[0]]), P1).is_fibration()
            True
        """
        return self.fan_morphism().is_fibration()

    def is_injective(self):
        r"""
        Check if ``self`` is injective.

        See
        :meth:`~sage.geometry.fan_morphism.FanMorphism.is_injective`
        for fan morphisms for a description of the toric algorithm.

        OUTPUT: boolean; whether ``self`` is injective

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: P1xP1.hom(matrix([[1],[0]]), P1).is_injective()
            False

            sage: X = toric_varieties.A(2)
            sage: m = identity_matrix(2)
            sage: f = X.hom(m, X)
            sage: f.is_injective()
            True

            sage: Y = ToricVariety(Fan([Cone([(1,0), (1,1)])]))
            sage: f = Y.hom(m, X)
            sage: f.is_injective()
            False
        """
        return self.fan_morphism().is_injective()

    def is_surjective(self):
        r"""
        Check if ``self`` is surjective.

        See
        :meth:`~sage.geometry.fan_morphism.FanMorphism.is_surjective`
        for fan morphisms for a description of the toric algorithm.

        OUTPUT: boolean; whether ``self`` is surjective

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: P1xP1.hom(matrix([[1],[0]]), P1).is_surjective()
            True

            sage: X = toric_varieties.A(2)
            sage: m = identity_matrix(2)
            sage: f = X.hom(m, X)
            sage: f.is_surjective()
            True

            sage: Y = ToricVariety(Fan([Cone([(1,0), (1,1)])]))
            sage: f = Y.hom(m, X)
            sage: f.is_surjective()
            False
        """
        return self.fan_morphism().is_surjective()

    def is_birational(self):
        r"""
        Check if ``self`` is birational.

        See
        :meth:`~sage.geometry.fan_morphism.FanMorphism.is_birational`
        for fan morphisms for a description of the toric algorithm.

        OUTPUT: boolean; whether ``self`` is birational

        EXAMPLES::

            sage: dP8 = toric_varieties.dP8()
            sage: P2 = toric_varieties.P2()
            sage: dP8.hom(identity_matrix(2), P2).is_birational()
            True

            sage: X = toric_varieties.A(2)
            sage: Y = ToricVariety(Fan([Cone([(1,0), (1,1)])]))
            sage: m = identity_matrix(2)
            sage: f = Y.hom(m, X)
            sage: f.is_birational()
            True
        """
        return self.fan_morphism().is_birational()

    def is_dominant(self):
        r"""
        Return whether ``self`` is dominant.

        See
        :meth:`~sage.geometry.fan_morphism.FanMorphism.is_dominant`
        for fan morphisms for a description of the toric algorithm.

        OUTPUT: boolean; whether ``self`` is a dominant scheme morphism

        EXAMPLES::

            sage: P1 = toric_varieties.P1()
            sage: A1 = toric_varieties.A1()
            sage: phi = A1.hom(identity_matrix(1), P1);  phi
            Scheme morphism:
              From: 1-d affine toric variety
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
              Defn: Defined by sending Rational polyhedral fan in 1-d lattice N
                    to Rational polyhedral fan in 1-d lattice N.
            sage: phi.is_dominant()
            True
            sage: phi.is_surjective()
            False
        """
        return self.fan_morphism().is_dominant()

    def pullback_divisor(self, divisor):
        r"""
        Pull back a toric divisor.

        INPUT:

        - ``divisor`` -- a torus-invariant `\QQ`-Cartier divisor on the
          codomain of ``self``

        OUTPUT: the pull-back divisor `f^*(D)`

        EXAMPLES::

            sage: A2_Z2 = toric_varieties.A2_Z2()
            sage: A2 = toric_varieties.A2()
            sage: f = A2.hom(matrix([[1,0], [1,2]]), A2_Z2)
            sage: f.pullback_divisor(A2_Z2.divisor(0))
            V(x)

            sage: A1 = toric_varieties.A1()
            sage: square = A1.hom(matrix([[2]]), A1)
            sage: D = A1.divisor(0);  D
            V(z)
            sage: square.pullback_divisor(D)
            2*V(z)
        """
        from sage.schemes.toric.divisor import ToricDivisor_generic
        if not (isinstance(divisor, ToricDivisor_generic) and divisor.is_QQ_Cartier()):
            raise ValueError('the divisor must be torus-invariant and QQ-Cartier')
        fm = self.fan_morphism()
        values = []
        for ray in self.domain().fan().rays():
            value = divisor.function_value(fm(ray))
            values.append(value)
        return self.domain().divisor(values)


############################################################################
# A morphism of toric varieties determined by a dominant fan morphism
class SchemeMorphism_fan_toric_variety_dominant(SchemeMorphism_fan_toric_variety):
    """
    Construct a morphism determined by a dominant fan morphism.

    A dominant morphism is one that is surjective onto a dense
    subset. In the context of toric morphisms, this means that it is
    onto the big torus orbit.

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.generic.scheme.hom` method of
        :class:`toric varieties
        <sage.schemes.toric.variety.ToricVariety_field>`
        instead.

    INPUT:

    See :class:`SchemeMorphism_fan_toric_variety`. The given fan
    morphism :meth:`must be dominant
    <sage.geometry.fan_morphism.FanMorphism.is_dominant>`.

    OUTPUT: a :class:`~sage.schemes.toric.morphism.SchemeMorphism_fan_toric_variety_dominant`

    EXAMPLES::

        sage: P2 = toric_varieties.P2()
        sage: dP8 = toric_varieties.dP8()
        sage: f = dP8.hom(identity_matrix(2), P2);  f
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   2-d CPR-Fano toric variety covered by 3 affine patches
          Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                to Rational polyhedral fan in 2-d lattice N.
        sage: type(f)
        <class 'sage.schemes.toric.morphism.SchemeMorphism_fan_toric_variety_dominant'>
    """

    @cached_method
    def fiber_generic(self):
        r"""
        Return the generic fiber.

        OUTPUT:

        - a tuple `(X, n)`, where `X` is a :class:`toric variety
          <sage.schemes.toric.variety.ToricVariety_field>` with the
          embedding morphism into domain of ``self`` and `n` is an integer.

        The fiber over the base point with homogeneous coordinates
        `[1:1:\cdots:1]` consists of `n` disjoint toric varieties isomorphic to
        `X`. Note that fibers of a dominant toric morphism are isomorphic over
        all points of a fixed torus orbit of its codomain, in particular over
        all points of the maximal torus, so it makes sense to talk about "the
        generic" fiber.

        The embedding of `X` is a toric morphism with
        the :meth:`~sage.geometry.fan_morphism.FanMorphism.domain_fan`
        being the
        :meth:`~sage.geometry.fan_morphism.FanMorphism.kernel_fan` of
        the defining fan morphism. By contrast, embeddings of fiber components
        over lower-dimensional torus orbits of the image are not toric
        morphisms. Use :meth:`fiber_component` for the latter
        (non-generic) fibers.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: fiber = P1xP1.hom(matrix([[1],[0]]), P1).fiber_generic()
            sage: fiber
            (1-d toric variety covered by 2 affine patches, 1)
            sage: f = fiber[0].embedding_morphism();  f
            Scheme morphism:
              From: 1-d toric variety covered by 2 affine patches
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined by sending Rational polyhedral fan in Sublattice <N(0, 1)> to
                    Rational polyhedral fan in 2-d lattice N.
            sage: f.as_polynomial_map()
            Scheme morphism:
              From: 1-d toric variety covered by 2 affine patches
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [z0 : z1] to [1 : 1 : z0 : z1]

            sage: A1 = toric_varieties.A1()
            sage: fan = Fan([(0,1,2)], [(1,1,0), (1,0,1), (1,-1,-1)])
            sage: fan = fan.subdivide(new_rays=[(1,0,0)])
            sage: f = ToricVariety(fan).hom(matrix([[1],[0],[0]]), A1)
            sage: f.fiber_generic()
            (2-d affine toric variety, 1)
            sage: _[0].fan().generating_cones()
            (0-d cone of Rational polyhedral fan in Sublattice <N(0, 1, 0), N(0, 0, 1)>,)
        """
        from sage.schemes.toric.variety import ToricVariety
        fm = self.fan_morphism()
        X = ToricVariety(fm.kernel_fan())
        m = X.fan().lattice().echelonized_basis_matrix()
        N = fm.domain()     # May be a sublattice as well
        m *= N.basis_matrix().solve_right(identity_matrix(N.dimension()))
        X._embedding_morphism = X.hom(m, self.domain())
        return X, fm.index()

    def fiber_component(self, domain_cone, multiplicity=False):
        r"""
        Return a fiber component corresponding to ``domain_cone``.

        INPUT:

        - ``domain_cone`` -- a cone of the domain fan of ``self``

        - ``multiplicity`` -- boolean (default: ``False``); whether to return
          the number of fiber components corresponding to ``domain_cone`` as well

        OUTPUT:

        - either `X` or a tuple `(X, n)`, where `X` is a :class:`toric variety
          <sage.schemes.toric.variety.ToricVariety_field>` with the
          embedding morphism into domain of ``self`` and `n` is an integer.

        Let `\phi: \Sigma \to \Sigma'` be the :class:`fan morphism
        <sage.geometry.fan_morphism.FanMorphism>` corresponding to
        ``self``. Let `\sigma \in \Sigma` and `\sigma' \in \Sigma'` be
        the :meth:`~sage.geometry.fan_morphism.FanMorphism.image_cone`
        of `\sigma`.  The fiber over any point of the torus orbit corresponding
        to `\sigma'` consists of `n` isomorphic connected components with each
        component being a union of toric varieties intersecting along
        their torus invariant subvarieties. The latter correspond to
        :meth:`~sage.geometry.fan_morphism.FanMorphism.preimage_cones` of
        `\sigma'` and `X` is one of the `n` components corresponding to
        `\sigma`. The irreducible components correspond to
        :meth:`~sage.geometry.fan_morphism.FanMorphism.primitive_preimage_cones`.

        EXAMPLES::

            sage: polytope = LatticePolytope(
            ....:     [(-3,0,-1,-1),(-1,2,-1,-1),(0,-1,0,0),(0,0,0,1),(0,0,1,0),
            ....:      (0,1,0,0),(0,2,-1,-1),(1,0,0,0),(2,0,-1,-1)])
            sage: coarse_fan = FaceFan(polytope)
            sage: P2 = toric_varieties.P2()
            sage: proj24 = matrix([[0,0], [1,0], [0,0], [0,1]])
            sage: fm = FanMorphism(proj24, coarse_fan, P2.fan(), subdivide=True)
            sage: fibration = ToricVariety(fm.domain_fan()).hom(fm, P2)
            sage: ffm = fibration.fan_morphism()
            sage: primitive_cones = ffm.primitive_preimage_cones(P2.fan(1)[0])
            sage: primitive_cone = primitive_cones[0]
            sage: fibration.fiber_component(primitive_cone)
            2-d toric variety covered by 4 affine patches
            sage: fibration.fiber_component(primitive_cone, True)
            (2-d toric variety covered by 4 affine patches, 1)

            sage: for primitive_cone in primitive_cones:
            ....:     print(fibration.fiber_component(primitive_cone))
            2-d toric variety covered by 4 affine patches
            2-d toric variety covered by 3 affine patches
            2-d toric variety covered by 3 affine patches
        """
        domain_cone = self.domain().fan().embed(domain_cone)
        if domain_cone.is_trivial():
            if multiplicity:
                return self.fiber_generic()
            else:
                return self.fiber_generic()[0]
        embedding = SchemeMorphism_fan_fiber_component_toric_variety(self, domain_cone)
        if multiplicity:
            return embedding.domain(), \
                self.fan_morphism().index(embedding.base_cone())
        else:
            return embedding.domain()

    @cached_method
    def fiber_dimension(self, codomain_cone):
        r"""
        Return the dimension of the fiber over a particular torus
        orbit in the base.

        INPUT:

        - ``codomain_cone`` -- a cone `\sigma` of the codomain,
          specifying a torus orbit `O(\sigma)`

        OUTPUT:

        An integer. The dimension of the fiber over the torus orbit
        corresponding to ``codomain_cone``. If the fiber is the empty
        set, ``-1`` is returned. Note that all fibers over this torus
        orbit are isomorphic, and therefore have the same dimension.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: f = P1xP1.hom(matrix([[1],[0]]), P1)
            sage: f.fiber_dimension(P1.fan(0)[0])
            1
            sage: f.fiber_dimension(P1.fan(1)[0])
            1
            sage: f.fiber_dimension(P1.fan(1)[1])
            1

        Here is a more complicated example that is not a flat fibration::

            sage: A2_Z2 = toric_varieties.A2_Z2()
            sage: O2_P1 = A2_Z2.resolve(new_rays=[(1,1)])
            sage: blowup = O2_P1.hom(identity_matrix(2), A2_Z2)
            sage: blowup.fiber_dimension(A2_Z2.fan(0)[0])
            0
            sage: blowup.fiber_dimension(A2_Z2.fan(1)[0])
            0
            sage: blowup.fiber_dimension(A2_Z2.fan(2)[0])
            1

        This corresponds to the three different fibers::

            sage: blowup.fiber_generic()
            (0-d affine toric variety, 1)
            sage: blowup.fiber_component(Cone([(1,0)]))
            0-d affine toric variety
            sage: blowup.fiber_component(Cone([(1,1)]))
            1-d toric variety covered by 2 affine patches
        """
        dim = []
        fm = self.fan_morphism()
        base_dim = codomain_cone.dim()
        dim.extend(base_dim - c.dim()
                   for c in fm.primitive_preimage_cones(codomain_cone))
        if dim:
            return max(dim) + self.domain().dimension() - self.codomain().dimension()
        return ZZ(-1)

    def fiber_graph(self, codomain_cone):
        r"""
        Return the fiber over a given torus orbit in the codomain.

        INPUT:

        - ``codomain_cone`` -- a cone `\sigma` of the codomain,
          specifying a torus orbit `O(\sigma)`

        OUTPUT:

        A graph whose nodes are the irreducible components of a connected
        component of the fiber over a point of `O(\sigma)`. If two irreducible
        components intersect, the
        corresponding nodes of the graph are joined by an edge. Note that
        irreducible components do not have to be of the same dimension.

        .. SEEALSO::

            :meth:`~SchemeMorphism_fan_toric_variety_dominant.fiber_component`.

        EXAMPLES::

            sage: polytope = Polyhedron(
            ....:     [(-3,0,-1,-1),(-1,2,-1,-1),(0,-1,0,0),(0,0,0,1),(0,0,1,0),
            ....:      (0,1,0,0),(0,2,-1,-1),(1,0,0,0),(2,0,-1,-1)])
            sage: coarse_fan = FaceFan(polytope, lattice=ToricLattice(4))

            sage: P2 = toric_varieties.P2()
            sage: proj34 = block_matrix(2, 1, [zero_matrix(2,2),
            ....:                              identity_matrix(2)])
            sage: fm = FanMorphism(proj34, coarse_fan, P2.fan(), subdivide=True)
            sage: fibration = ToricVariety(fm.domain_fan()).hom(fm, P2)

            sage: fibration.fiber_graph(P2.fan(0)[0])
            Graph on 1 vertex
            sage: for c1 in P2.fan(1):
            ....:     fibration.fiber_graph(c1)
            Graph on 1 vertex
            Graph on 1 vertex
            Graph on 4 vertices

            sage: fibration.fiber_graph(P2.fan(1)[2]).get_vertices()
            {0: 2-d toric variety covered by 4 affine patches,
             1: 2-d toric variety covered by 3 affine patches,
             2: 2-d toric variety covered by 3 affine patches,
             3: 2-d toric variety covered by 4 affine patches}

            sage: fibration
            Scheme morphism:
              From: 4-d toric variety covered by 18 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by sending Rational polyhedral fan in 4-d lattice N
                    to Rational polyhedral fan in 2-d lattice N.
        """
        fm = self.fan_morphism()
        prim = fm.primitive_preimage_cones(codomain_cone)
        n = len(prim)

        def is_union_in_fan(self, c0, c1):
            indices = c0.ambient_ray_indices() + c1.ambient_ray_indices()
            try:
                fm.domain_fan().cone_containing(*indices)
                return True
            except ValueError:
                return False

        m = matrix(ZZ, n, n, lambda i, j: is_union_in_fan(self, prim[i], prim[j]))

        for i in range(n):
            m[i, i] = 0
        from sage.graphs.graph import Graph
        graph = Graph(m, loops=False, multiedges=False)
        for i in range(n):
            graph.set_vertex(i, self.fiber_component(prim[i]))
        return graph


############################################################################
# The embedding morphism of a fiber component
class SchemeMorphism_fan_fiber_component_toric_variety(SchemeMorphism):
    """
    The embedding of a fiber component of a toric morphism.

    Note that the embedding map of a fiber component of a toric morphism is
    itself not a toric morphism!

    INPUT:

    - ``toric_morphism`` -- a toric morphism; the toric morphism whose
      fiber component we are describing

    - ``defining_cone`` -- a cone of the fan of the domain of
      ``toric_morphism``; see
      :meth:`~SchemeMorphism_fan_toric_variety_dominant.fiber_component` for
      details

    EXAMPLES::

        sage: polytope = Polyhedron(
        ....:     [(-3,0,-1,-1),(-1,2,-1,-1),(0,-1,0,0),(0,0,0,1),(0,0,1,0),
        ....:      (0,1,0,0),(0,2,-1,-1),(1,0,0,0),(2,0,-1,-1)])
        sage: coarse_fan = FaceFan(polytope, lattice=ToricLattice(4))
        sage: P2 = toric_varieties.P2()
        sage: proj24 = matrix([[0,0],[1,0],[0,0],[0,1]])
        sage: fm = FanMorphism(proj24, coarse_fan, P2.fan(), subdivide=True)
        sage: fibration = ToricVariety(fm.domain_fan()).hom(fm, P2)
        sage: ffm = fibration.fan_morphism()
        sage: primitive_cones = ffm.primitive_preimage_cones(P2.fan(1)[0])
        sage: primitive_cone = primitive_cones[0]
        sage: fiber_component = fibration.fiber_component(primitive_cone)
        sage: fiber_component
        2-d toric variety covered by 4 affine patches
        sage: fiber_component.embedding_morphism()
        Scheme morphism:
          From: 2-d toric variety covered by 4 affine patches
          To:   4-d toric variety covered by 23 affine patches
          Defn: Defined by embedding a fiber component corresponding to
                1-d cone of Rational polyhedral fan in 4-d lattice N.
        sage: fiber_component.embedding_morphism().as_polynomial_map()
        Scheme morphism:
          From: 2-d toric variety covered by 4 affine patches
          To:   4-d toric variety covered by 23 affine patches
          Defn: Defined on coordinates by sending [z0 : z1 : z2 : z3] to
                [1 : 1 : 1 : 1 : z2 : 0 : 1 : z3 : 1 : 1 : 1 : z1 : z0 : 1 : 1]
        sage: type(fiber_component.embedding_morphism())
        <class 'sage.schemes.toric.morphism.SchemeMorphism_fan_fiber_component_toric_variety'>
    """

    def __init__(self, toric_morphism, defining_cone):
        """
        The Python constructor.

        TESTS::

            sage: polytope = Polyhedron(
            ....:     [(-3,0,-1,-1),(-1,2,-1,-1),(0,-1,0,0),(0,0,0,1),(0,0,1,0),
            ....:      (0,1,0,0),(0,2,-1,-1),(1,0,0,0),(2,0,-1,-1)])
            sage: coarse_fan = FaceFan(polytope, lattice=ToricLattice(4))
            sage: P2 = toric_varieties.P2()
            sage: proj24 = matrix([[0,0], [1,0], [0,0], [0,1]])
            sage: fm = FanMorphism(proj24, coarse_fan, P2.fan(), subdivide=True)
            sage: fibration = ToricVariety(fm.domain_fan()).hom(fm, P2)
            sage: primitive_cone = Cone([(-1, 2, -1, 0)])
            sage: fibration.fiber_component(primitive_cone).embedding_morphism()
            Scheme morphism:
              From: 2-d toric variety covered by 3 affine patches
              To:   4-d toric variety covered by 23 affine patches
              Defn: Defined by embedding a fiber component corresponding to
                    1-d cone of Rational polyhedral fan in 4-d lattice N.
        """
        fm = toric_morphism.fan_morphism()
        self._fan_morphism = fm
        defining_cone = fm.domain_fan().embed(defining_cone)
        self._defining_cone = defining_cone
        self._base_cone = fm.image_cone(defining_cone)
        fc = self._make_fiber_component()
        fc._embedding_morphism = self
        parent = fc.Hom(toric_morphism.domain())
        SchemeMorphism.__init__(self, parent)

    def _repr_defn(self):
        """
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: fc = P1xP1.hom(matrix([[1],[0]]), P1).fiber_component(Cone([(1,0)]))
            sage: fc.embedding_morphism()._repr_defn()
            'Defined by embedding a fiber component corresponding to 1-d cone of Rational polyhedral fan in 2-d lattice N.'
        """
        return 'Defined by embedding a fiber component corresponding to {}.'.format(self.defining_cone())

    def as_polynomial_map(self):
        """
        Express the embedding morphism via homogeneous polynomials.

        OUTPUT: a :class:`SchemeMorphism_polynomial_toric_variety`

        This raises a :exc:`ValueError` if the morphism cannot be
        written in terms of homogeneous polynomials.

        EXAMPLES::

            sage: polytope = Polyhedron(
            ....:     [(-3,0,-1,-1),(-1,2,-1,-1),(0,-1,0,0),(0,0,0,1),(0,0,1,0),
            ....:      (0,1,0,0),(0,2,-1,-1),(1,0,0,0),(2,0,-1,-1)])
            sage: coarse_fan = FaceFan(polytope, lattice=ToricLattice(4))
            sage: P2 = toric_varieties.P2()
            sage: proj24 = matrix([[0,0], [1,0], [0,0], [0,1]])
            sage: fm = FanMorphism(proj24, coarse_fan, P2.fan(), subdivide=True)
            sage: fibration = ToricVariety(fm.domain_fan()).hom(fm, P2)

            sage: primitive_cone = Cone([(0, 1, 0, 0)])
            sage: f = fibration.fiber_component(primitive_cone).embedding_morphism()
            sage: f.as_polynomial_map()
            Scheme morphism:
              From: 2-d toric variety covered by 4 affine patches
              To:   4-d toric variety covered by 23 affine patches
              Defn: Defined on coordinates by sending [z0 : z1 : z2 : z3] to
                    [1 : 1 : 1 : 1 : z2 : 0 : 1 : z3 : 1 : 1 : 1 : z1 : z0 : 1 : 1]

            sage: primitive_cone = Cone([(-1, 2, -1, 0)])
            sage: f = fibration.fiber_component(primitive_cone).embedding_morphism()
            sage: f.as_polynomial_map()
            Traceback (most recent call last):
            ...
            ValueError: the morphism cannot be written using homogeneous polynomials
        """
        fc = self.domain()
        toric_variety = self.codomain()
        R = fc.coordinate_ring()
        polys = [R.one()] * toric_variety.fan().nrays()
        for i in self.defining_cone().ambient_ray_indices():
            polys[i] = R.zero()
        for ray, x in zip(fc.fan().rays(), R.gens()):
            try:
                ray_index = self._ray_index_map[ray]
            except KeyError:
                raise ValueError('the morphism cannot be written using homogeneous polynomials')
            polys[ray_index] = x
        return SchemeMorphism_polynomial_toric_variety(self.parent(), polys)

    def _make_fiber_component(self):
        """
        Construct the fiber component as a toric variety.

        OUTPUT: the fiber component as a toric variety

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: fc = P1xP1.hom(matrix([[1],[0]]), P1).fiber_component(Cone([(1,0)]))
            sage: f = fc.embedding_morphism()
            sage: f._ray_index_map  # indirect doctest
            {N(-1): 3, N(1): 2}

        TESTS::

            sage: A2_Z2 = toric_varieties.A2_Z2()
            sage: O2_P1 = A2_Z2.resolve(new_rays=[(1,1)])
            sage: blowup = O2_P1.hom(identity_matrix(2), A2_Z2)
            sage: blowup.fiber_generic()
            (0-d affine toric variety, 1)
            sage: blowup.fiber_component(Cone([(1,0)]))
            0-d affine toric variety
            sage: blowup.fiber_component(Cone([(1,1)]))
            1-d toric variety covered by 2 affine patches

            sage: P1 = toric_varieties.P1()
            sage: f = P1.hom(matrix([2]), P1)
            sage: f.fiber_component(P1.fan(1)[0])
            0-d affine toric variety
            sage: f.fan_morphism().index(P1.fan(1)[0])
            1
            sage: f.fiber_generic()
            (0-d affine toric variety, 2)
        """
        fm = self._fan_morphism
        defining_cone = self._defining_cone
        base_cone = self._base_cone

        ker = fm.kernel().basis()
        m = fm.matrix() * base_cone.lattice().basis_matrix()
        base_cone_preimg = [m.solve_left(r) for r in base_cone.rays()]
        L = fm.domain_fan().lattice().span(ker + base_cone_preimg).saturation()

        cone_L = Cone([L.coordinates(r) for r in defining_cone.rays()])
        L_quotient = cone_L.sublattice_quotient()

        def projection(ray):
            ray_L = L.coordinates(ray)
            return vector(ZZ, L_quotient(ray_L))

        cones = []
        star_rays = set()
        for cone in fm.relative_star_generators(defining_cone):
            star_rays.update(cone.rays())
            projected_rays = [projection(r) for r in cone.rays()]
            cones.append(Cone(projected_rays))
        fiber_fan = Fan(cones)

        ray_index_map = {}
        for ray in star_rays:
            ray_index = fm.domain_fan().rays().index(ray)
            projected_ray = fiber_fan.lattice()(projection(ray))
            if projected_ray.is_zero():
                assert ray in defining_cone.rays()
                continue
            projected_ray.set_immutable()
            ray_index_map[projected_ray] = ray_index
        self._ray_index_map = ray_index_map

        from sage.schemes.toric.variety import ToricVariety
        return ToricVariety(fiber_fan)

    def defining_cone(self):
        r"""
        Return the cone corresponding to the fiber torus orbit.

        OUTPUT: a cone of the fan of the total space of the toric fibration

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: fc = P1xP1.hom(matrix([[1],[0]]), P1).fiber_component(Cone([(1,0)]))
            sage: f = fc.embedding_morphism()
            sage: f.defining_cone().rays()
            N(1, 0)
            in 2-d lattice N
            sage: f.base_cone().rays()
            N(1)
            in 1-d lattice N
        """
        return self._defining_cone

    def base_cone(self):
        r"""
        Return the base cone `\sigma`.

        The fiber is constant over the base orbit closure `V(\sigma)`.

        OUTPUT: a cone of the base of the toric fibration

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: fc = P1xP1.hom(matrix([[1],[0]]), P1).fiber_component(Cone([(1,0)]))
            sage: f = fc.embedding_morphism()
            sage: f.defining_cone().rays()
            N(1, 0)
            in 2-d lattice N
            sage: f.base_cone().rays()
            N(1)
            in 1-d lattice N
        """
        return self._base_cone

    def _image_ray_multiplicity(self, fiber_ray):
        """
        Find the image ray of ``fiber_ray`` with multiplicity in the relative star.

        INPUT:

        - ``fiber_ray`` -- a ray of the domain fan (the fiber component)

        OUTPUT: a pair ``(codomain ray index, multiplicity)``

        EXAMPLES::

            sage: polytope = Polyhedron(
            ....:     [(-3,0,-1,-1),(-1,2,-1,-1),(0,-1,0,0),(0,0,0,1),(0,0,1,0),
            ....:      (0,1,0,0),(0,2,-1,-1),(1,0,0,0),(2,0,-1,-1)])
            sage: coarse_fan = FaceFan(polytope, lattice=ToricLattice(4))
            sage: P2 = toric_varieties.P2()
            sage: proj24 = matrix([[0,0], [1,0], [0,0], [0,1]])
            sage: fm = FanMorphism(proj24, coarse_fan, P2.fan(), subdivide=True)
            sage: fibration = ToricVariety(fm.domain_fan()).hom(fm, P2)
            sage: primitive_cone = Cone([(-1, 2, -1, 0)])
            sage: fc = fibration.fiber_component(primitive_cone)
            sage: f = fc.embedding_morphism()
            sage: for r in fc.fan().rays():
            ....:     print("{} {}".format(r, f._image_ray_multiplicity(r)))
            N(-1, -1) (9, 2)
            N(0, 1) (5, 1)
            N(1, 0) (11, 1)
            sage: f._ray_index_map
            {N(-2, -2): 9, N(-1, 2): 4, N(0, 1): 5, N(1, 0): 11, N(3, -2): 10}
        """
        try:
            image_ray_index = self._ray_index_map[fiber_ray]
            return (image_ray_index, 1)
        except KeyError:
            pass
        multiplicity = None
        image_ray_index = None
        for ray, index in self._ray_index_map.items():
            d = gcd(ray)
            if d * fiber_ray != ray:
                continue
            if multiplicity is not None and d > multiplicity:
                continue
            multiplicity = d
            image_ray_index = index
        return (image_ray_index, multiplicity)

    def pullback_divisor(self, divisor):
        r"""
        Pull back a toric divisor.

        INPUT:

        - ``divisor`` -- a torus-invariant `\QQ`-Cartier divisor on the
          codomain of the embedding map

        OUTPUT:

        A divisor on the domain of the embedding map (irreducible
        component of a fiber of a toric morphism) that is isomorphic
        to the pull-back divisor `f^*(D)` but with possibly different
        linearization.

        EXAMPLES::

            sage: A1 = toric_varieties.A1()
            sage: fan = Fan([(0,1,2)], [(1,1,0),(1,0,1),(1,-1,-1)]).subdivide(new_rays=[(1,0,0)])
            sage: f = ToricVariety(fan).hom(matrix([[1],[0],[0]]), A1)
            sage: D = f.domain().divisor([1,1,3,4]); D
            V(z0) + V(z1) + 3*V(z2) + 4*V(z3)
            sage: fc = f.fiber_component(Cone([(1,1,0)]))
            sage: fc.embedding_morphism().pullback_divisor(D)
            4*V(z0) + V(z1) + 4*V(z2)
            sage: fc = f.fiber_component(Cone([(1,0,0)]))
            sage: fc.embedding_morphism().pullback_divisor(D)
            -V(z0) - 3*V(z1) - 3*V(z2)
        """
        from sage.schemes.toric.divisor import ToricDivisor_generic
        if not (isinstance(divisor, ToricDivisor_generic) and divisor.is_QQ_Cartier()):
            raise ValueError('the divisor must be torus-invariant and QQ-Cartier')
        m = divisor.m(self.defining_cone())
        values = []
        codomain_rays = self.codomain().fan().rays()
        for ray in self.domain().fan().rays():
            image_ray_index, multiplicity = self._image_ray_multiplicity(ray)
            image_ray = codomain_rays[image_ray_index]
            value = divisor.function_value(image_ray) - m * image_ray
            value /= multiplicity
            values.append(value)
        return self.domain().divisor(values)
