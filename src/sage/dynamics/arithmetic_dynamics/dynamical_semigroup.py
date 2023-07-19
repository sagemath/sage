r"""
Dynamical semigroups

A dynamical semigroup is a finitely generated subsemigroup of
the endomorphism ring of a subscheme of projective or affine space.

AUTHORS:

 - Dang Phan (July 12th, 2023): initial implementation
"""

#*****************************************************************************
# Dang Phan <dang8phan@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from collections.abc import Collection

from sage.structure.parent import Parent

from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sage.categories.number_fields import NumberFields
from sage.rings.rational_field import QQ
from sage.categories.semigroups import Semigroups

class DynamicalSemigroup(Parent):
    r"""
    A dynamical semigroup defined by a multiple dynamical systems on projective or affine space.

    INPUT:

    - ``ds_data`` -- list or tuple of dynamical systems or objects that define dynamical systems

    OUTPUT:

    :class:`DynamicalSemigroup_affine` if ``ds_data`` only contains or defines dynamical systems
    over affine space; and :class:`DynamicalSemigroup_projective` otherwise.

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: DynamicalSemigroup(([x, y], [x^2, y^2]))
        Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([x^2, y^2], P)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    ::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(x^2, A)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 2 dynamical systems:
        Dynamical System of Affine Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x) to
                (x)
        Dynamical System of Affine Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x) to
                (x^2)

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: X = P.subscheme(x - y)
        sage: f = DynamicalSystem_projective([x, y], X)
        sage: g = DynamicalSystem_projective([x^2, y^2], X)
        sage: DynamicalSemigroup_projective([f, g])
        Dynamical semigroup over Closed subscheme of Projective Space of dimension 1 over Rational Field defined by:
          x - y defined by 2 dynamical systems:
        Dynamical System of Closed subscheme of Projective Space of dimension 1 over Rational Field defined by:
          x - y
        Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Closed subscheme of Projective Space of dimension 1 over Rational Field defined by:
          x - y
        Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    If a dynamical semigroup is built from dynamical systems with different base rings, all systems will be coerced
    to the largest base ring::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: Q.<z,w> = ProjectiveSpace(RR, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Real Field with 53 bits of precision defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Real Field with 53 bits of precision
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Real Field with 53 bits of precision
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    ::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: B.<y> = AffineSpace(RR, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(y^2, B)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over Real Field with 53 bits of precision defined by 2 dynamical systems:
        Dynamical System of Affine Space of dimension 1 over Real Field with 53 bits of precision
          Defn: Defined on coordinates by sending (x) to
                (x)
        Dynamical System of Affine Space of dimension 1 over Real Field with 53 bits of precision
          Defn: Defined on coordinates by sending (x) to
                (x^2)

    If a dynamical semigroup is built from dynamical systems over number fields, a composite number field is created
    and all systems will be coerced to it. This composite number field contains all of the initial number fields::

        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: Q.<x,y> = ProjectiveSpace(K, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Number Field in k with defining polynomial r^2 - 2 defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Number Field in k with defining polynomial r^2 - 2
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Number Field in k with defining polynomial r^2 - 2
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    ::

        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: L.<l> = NumberField(r^2 - 3)
        sage: P.<x,y> = ProjectiveSpace(K, 1)
        sage: Q.<z,w> = ProjectiveSpace(L, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Number Field in kl with defining polynomial r^4 - 10*r^2 + 1 defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    ::

        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: L.<l> = NumberField(r^2 - 3)
        sage: P.<x> = AffineSpace(K, 1)
        sage: Q.<y> = AffineSpace(L, 1)
        sage: f = DynamicalSystem(x, P)
        sage: g = DynamicalSystem(y^2, Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over Number Field in kl with defining polynomial r^4 - 10*r^2 + 1 defined by 2 dynamical systems:
        Dynamical System of Affine Space of dimension 1 over Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
          Defn: Defined on coordinates by sending (x) to
                (x)
        Dynamical System of Affine Space of dimension 1 over Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
          Defn: Defined on coordinates by sending (x) to
                (x^2)

    If a dynamical semigroup is built from dynamical systems over both projective and affine spaces, all systems
    will be homogenized to dynamical systems over projective space::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: A.<z> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem(z^2, A)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

    TESTS::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: DynamicalSemigroup(1)
        ValueError: 1 does not define a 'DynamicalSystem' object

    ::

        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: P.<x,y> = ProjectiveSpace(RR, 1)
        sage: Q.<z,w> = ProjectiveSpace(K, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        ValueError: given dynamical systems are not automorphic under global composition

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: Q.<u,v,w> = ProjectiveSpace(QQ, 2)
        sage: f = DynamicalSystem([x, y])
        sage: g = DynamicalSystem([u^2, v^2, w^2])
        sage: DynamicalSemigroup((f, g))
        ValueError: domains of 'DynamicalSystem' objects must be of the same dimension
    """

    def __init__(self, ds_data):
        r"""
        The Python constructor.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: DynamicalSemigroup(([x, y], [x^2, y^2]))
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x : y)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)
        """
        systems = []

        if isinstance(ds_data, Collection):
            for ds_datum in ds_data:
                if isinstance(ds_datum, DynamicalSystem):
                    systems.append(ds_datum)
                else:
                    try:
                        systems.append(DynamicalSystem(ds_datum))
                    except ValueError:
                        raise ValueError(str(ds_datum) + " does not define a 'DynamicalSystem' object")
        else:
            if isinstance(ds_data, DynamicalSystem):
                    systems.append(ds_data)
            else:
                try:
                    systems.append(DynamicalSystem(ds_data))
                except ValueError:
                    raise ValueError(str(ds_data) + " does not define a 'DynamicalSystem' object")

        exists_projective_ds = any(isinstance(element, DynamicalSystem_projective) for element in systems)

        if exists_projective_ds:
            dimension = None
            for ds in systems:
                if isinstance(ds, DynamicalSystem_projective):
                    dimension = ds.domain().dimension()
                    break
            for i in range(len(systems)):
                if isinstance(systems[i], DynamicalSystem_affine):
                    systems[i] = systems[i].homogenize(dimension)

        identical_domains = True
        for ds in systems:
            if ds.domain() != systems[0].domain():
                identical_domains = False
                break

        over_number_fields = True
        all_over_QQ = True
        for ds in systems:
            if ds.base_ring() not in NumberFields():
                over_number_fields = False
            if ds.base_ring() is not QQ:
                all_over_QQ = False

        biggest_ring = None

        if over_number_fields and not all_over_QQ:
            number_fields = []
            for ds in systems:
                number_fields.append(ds.base_ring())

            minimal_composite_field = None
            for field in number_fields:
                if field is not QQ:
                    if minimal_composite_field is None:
                        minimal_composite_field = field
                    else:
                        minimal_composite_field = minimal_composite_field.composite_fields(field)[0]

            biggest_ring = minimal_composite_field
        else:
            for ds in systems:
                if biggest_ring is None:
                    biggest_ring = ds.base_ring()
                elif ds.base_ring().has_coerce_map_from(biggest_ring):
                    biggest_ring = ds.base_ring()
                elif biggest_ring.has_coerce_map_from(ds.base_ring()):
                    pass
                else:
                    raise ValueError("given dynamical systems are not automorphic"
                                    " under global composition")

        for i in range(len(systems)):
            if systems[i].base_ring() != biggest_ring:
                    systems[i] = systems[i].change_ring(biggest_ring)

        identical_domains = all(ds.domain() == systems[0].domain() for ds in systems)
        if not identical_domains:
            for ds in systems:
                if ds.domain().dimension() != systems[0].domain().dimension():
                    raise ValueError("domains of 'DynamicalSystem' objects must be of the same dimension")

        gens = systems[0].domain().ambient_space().gens()
        for i in range(len(systems)):
            if systems[i].domain().coordinate_ring() != systems[0].domain().coordinate_ring():
                sub_dict = {}
                old_gens = systems[i].domain().ambient_space().gens()
                for j in range(len(old_gens)):
                    sub_dict[old_gens[j]] = gens[j]
                new_polys = []
                for poly in systems[i].defining_polynomials():
                    new_polys.append(poly.subs(sub_dict))
                if exists_projective_ds:
                    systems[i] = DynamicalSystem_projective(new_polys)
                else:
                    systems[i] = DynamicalSystem_affine(new_polys)

        self._dynamical_systems = systems
        self._domain = systems[0].domain()
        self._codomain = systems[0].codomain()
        self._dimension = self._domain.dimension()

        Parent.__init__(self, category=Semigroups().Finite().FinitelyGenerated())

    def __call__(self, input):
        r"""
        The result after evaluating this dynamical semigroup on a list or tuple of values.

        INPUT:

        - ``input`` -- one element, a list, or a tuple of values that can be evaluated
          with the generators of this dynamical semigroup.

        OUTPUT: A tuple of the resulting values after applying all of this dynamical semigroup's generators to ``input``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f(2)
            ((2 : 1), (4 : 1))

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f([2, 1])
            ((2 : 1), (4 : 1))

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f(f(2))
            TypeError: unable to convert (2 : 1) to an element of Rational Field
        """
        result = []
        for ds in self._dynamical_systems:
            result.append(ds(self.domain()(input)))
        return tuple(result)

    def base_ring(self):
        r"""
        The base ring of this dynamical semigroup. This is identical
        to the base ring of all of its defining dynamical system.

        OUTPUT: A ring.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.base_ring()
            Rational Field
        """
        return self._dynamical_systems[0].base_ring()

    def change_ring(self, new_ring):
        r"""
        Return a new :class:`DynamicalSemigroup` whose generators
        are the initial dynamical systems coerced to ``new_ring``.

        INPUT:

        - ``new_ring`` -- a ring.

        OUTPUT:

        A :class:`DynamicalSemigroup` defined by this dynamical
        semigroup's generators, but coerced to ``new_ring``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.change_ring(RR)
            Dynamical semigroup over Projective Space of dimension 1 over Real Field with 53 bits of precision defined by 2 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Real Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y) to
                    (x : y)
            Dynamical System of Projective Space of dimension 1 over Real Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)
        """
        new_systems = []
        for ds in self._dynamical_systems:
            new_systems.append(ds.change_ring(new_ring))
        return DynamicalSemigroup_projective(new_systems)

    def domain(self):
        r"""
        Return the domain of the generators of this dynamical semigroup.

        OUTPUT: A subscheme of a projective space or affine space.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.domain()
            Projective Space of dimension 1 over Rational Field
        """
        return self._domain

    def codomain(self):
        r"""
        Return the codomain of the generators of this dynamical semigroup.

        OUTPUT: A subscheme of a projective space or affine space.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.codomain()
            Projective Space of dimension 1 over Rational Field
        """
        return self._codomain

    def defining_polynomials(self):
        r"""
        Return the tuple of polynomials that define the generators of this dynamical semigroup.

        OUTPUT: A tuple of polynomials.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.defining_polynomials()
            ((x, y), (x^2, y^2))
        """
        result = []
        for ds in self._dynamical_systems:
            result.append(ds.defining_polynomials())
        return tuple(result)

    def defining_systems(self):
        r"""
        Return the generators of this dynamical semigroup.

        OUTPUT: A tuple of dynamical systems.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.defining_systems()
            (Dynamical System of Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y) to
                     (x : y),
             Dynamical System of Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y) to
                     (x^2 : y^2))
        """
        return tuple(self._dynamical_systems)

    def _repr_(self):
        r"""
        Return the :class:`String` representation of this dynamical semigroup.

        OUTPUT: A :class:`String` displaying information about this dynamical semigroup.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: DynamicalSemigroup(([x, y], [x^2, y^2]))
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x : y)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)
        """
        header = "Dynamical semigroup over %s defined by %d dynamical system"
        if (len(self._dynamical_systems) > 1):
            header += "s"
        header += ":"
        header = header % (str(self.domain()), len(self._dynamical_systems))
        systems = []
        for ds in self._dynamical_systems:
            systems.append(str(ds))
        systems = '\n'.join(systems)
        return header + "\n" + systems

    def __eq__(self, other):
        r"""
        Return whether two dynamical semigroups are equal.

        OUTPUT:

        A boolean that is True if and only if the two dynamical semigroups contain
        the same dynamical systems in the same order.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: g = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f == g
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: g = DynamicalSemigroup(([x^2, y^2], [x, y]))
            sage: f == g
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: g = DynamicalSemigroup(([x^2, y^2], [y, x]))
            sage: f == g
            False
        """
        if isinstance(other, DynamicalSemigroup):
            return len(self._dynamical_systems) == len(other._dynamical_systems) and \
                all(systems in other._dynamical_systems for systems in self._dynamical_systems)
        return False

class DynamicalSemigroup_projective(DynamicalSemigroup):
    r"""
    A dynamical semigroup defined by a multiple dynamical systems on projective space.

    INPUT:

    - ``ds_data`` -- list or tuple of dynamical systems or objects that define dynamical systems
      over projective space.

    OUTPUT: :class:`DynamicalSemigroup_projective`

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: DynamicalSemigroup_projective(([x, y], [x^2, y^2]))
        Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)
    """
    pass

class DynamicalSemigroup_affine(DynamicalSemigroup):
    r"""
    A dynamical semigroup defined by a multiple dynamical systems on affine space.

    INPUT:

    - ``ds_data`` -- list or tuple of dynamical systems or objects that define dynamical systems
      over affine space.

    OUTPUT: :class:`DynamicalSemigroup_affine`

    EXAMPLES::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(x^2, A)
        sage: DynamicalSemigroup_affine((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 2 dynamical systems:
        Dynamical System of Affine Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x) to
                (x)
        Dynamical System of Affine Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x) to
                (x^2)
    """

    def __init__(self, ds_data):
        r"""
        The Python constructor.

        EXAMPLES::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(x^2, A)
        sage: DynamicalSemigroup_affine((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 2 dynamical systems:
        Dynamical System of Affine Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x) to
                (x)
        Dynamical System of Affine Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x) to
                (x^2)
        """
        super().__init__(ds_data)
        for i in range(len(self._dynamical_systems)):
            if isinstance(self._dynamical_systems[i], DynamicalSystem_projective):
                self._dynamical_systems[i] = self._dynamical_systems[i].dehomogenize(self._dimension)