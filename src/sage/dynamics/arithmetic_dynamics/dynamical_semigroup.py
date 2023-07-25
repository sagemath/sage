r"""
Dynamical semigroups

A dynamical semigroup is a finitely generated subsemigroup of
the endomorphism ring of a subscheme of projective or affine space.

AUTHORS:

 - Dang Phan (July 21st, 2023): initial implementation
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
from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.categories.semigroups import Semigroups
from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
from sage.misc.classcall_metaclass import typecall
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent

class DynamicalSemigroup(Parent, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A dynamical semigroup defined by a multiple dynamical systems on projective or affine space.

    INPUT:

    - ``ds_data`` -- list or tuple of dynamical systems or objects that define dynamical systems

    OUTPUT:

    :class:`DynamicalSemigroup_affine` if ``ds_data`` only contains dynamical systems
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

    A dynamical semigroup may contain dynamical systems over finite fields::

        sage: F = FiniteField(5)
        sage: P.<x,y> = ProjectiveSpace(F, 1)
        sage: DynamicalSemigroup(([x, y], [x^2, y^2]))
        Dynamical semigroup over Projective Space of dimension 1 over Finite Field of size 5 defined by 2 dynamical systems:
        Dynamical System of Projective Space of dimension 1 over Finite Field of size 5
          Defn: Defined on coordinates by sending (x : y) to
                (x : y)
        Dynamical System of Projective Space of dimension 1 over Finite Field of size 5
          Defn: Defined on coordinates by sending (x : y) to
                (x^2 : y^2)

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
        Traceback (most recent call last):
        ...
        TypeError: 1 does not define a 'DynamicalSemigroup' object

    ::

        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: P.<x,y> = ProjectiveSpace(RR, 1)
        sage: Q.<z,w> = ProjectiveSpace(K, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Traceback (most recent call last):
        ...
        ValueError: given dynamical systems are not automorphic under global composition

    ::

        sage: F = FiniteField(5)
        sage: P.<x,y> = ProjectiveSpace(F, 1)
        sage: Q.<z,w> = ProjectiveSpace(QQ, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Traceback (most recent call last):
        ...
        ValueError: given dynamical systems are not automorphic under global composition

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: Q.<u,v,w> = ProjectiveSpace(QQ, 2)
        sage: f = DynamicalSystem([x, y])
        sage: g = DynamicalSystem([u^2, v^2, w^2])
        sage: DynamicalSemigroup((f, g))
        Traceback (most recent call last):
        ...
        ValueError: domains of 'DynamicalSystem' objects must be of the same dimension
    """

    @staticmethod
    def __classcall_private__(cls, ds_data):
        if isinstance(ds_data, Collection):
            all_affine_systems = all(isinstance(ds_datum, DynamicalSystem_affine) for ds_datum in ds_data)
            if all_affine_systems:
                return DynamicalSemigroup_affine(ds_data)
        elif not isinstance(ds_data, DynamicalSystem):
            raise TypeError(str(ds_data) + " does not define a 'DynamicalSemigroup' object")
        return DynamicalSemigroup_projective(ds_data)

    def __init__(self, systems):
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

        self._domain = systems[0].domain()
        self._codomain = systems[0].codomain()
        self._dynamical_systems = systems
        Parent.__init__(self, category=Semigroups().Finite().FinitelyGeneratedAsMagma())

    def __call__(self, input):
        r"""
        The result after evaluating this dynamical semigroup on a value.

        INPUT:

        - ``input`` -- one value that can be evaluated
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
            Traceback (most recent call last):
            ...
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

    def nth_iterate(self, p, n):
        r"""
        Return a tuple of values that results from evaluating this dynamical semigroup
        on the value ``p`` a total of ``n`` times.

        INPUT:

        - ``p`` -- a value on which dynamical systems can evaluate
        - ``n`` -- a positive integer

        OUTPUT: a tuple of values

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: f.nth_iterate(2, 3)
            ((2 : 1), (9 : 1), (5/3 : 1), (16 : 1))
        """
        if not isinstance(n, Integer):
            raise TypeError(str(n) + " must be an integer")
        if not n > 0:
            raise ValueError(str(n) + " must be a positive integer")
        result = self(p)
        if n == 1:
            return result
        for i in range(n - 2):
            next_iteration = []
            for point in result:
                next_iteration.extend(self(point))
            result = next_iteration
        return tuple(result)

    def __pow__(self, n):
        r"""
        Return a new dynamical semigroup that is the product of this dynamical semigroup and itself ``n`` times.

        INPUT:

        - ``n`` -- a positive integer

        OUTPUT: :class:`DynamicalSemigroup`

        EXAMPLES:

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f^2
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 8 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x : y)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 : y^4)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x : y)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 : y^4)
        """
        if not isinstance(n, Integer):
            raise TypeError(str(n) + " must be an integer")
        if not n > 0:
            raise ValueError(str(n) + " must be a positive integer")
        result = self
        for i in range(n - 1):
            result = result * self
        return result

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

    @staticmethod
    def __classcall_private__(cls, ds_data):
        systems = []

        if isinstance(ds_data, Collection):
            for ds_datum in ds_data:
                if isinstance(ds_datum, DynamicalSystem_projective):
                    systems.append(ds_datum)
                elif isinstance(ds_datum, DynamicalSystem_affine):
                    dimension = ds_datum.domain().dimension()
                    systems.append(ds_datum.homogenize(dimension))
                else:
                    try:
                        systems.append(DynamicalSystem_projective(ds_datum))
                    except ValueError:
                        raise ValueError(str(ds_datum) + " does not define a 'DynamicalSystem_projective' object")
        else:
            if isinstance(ds_data, DynamicalSystem_projective):
                    systems.append(ds_data)
            else:
                try:
                    systems.append(DynamicalSystem_projective(ds_data))
                except ValueError:
                    raise ValueError(str(ds_data) + " does not define a 'DynamicalSystem_projective' object")

        systems = _standardize_domains_of_(systems)
        if systems[0].base_ring() not in Fields():
            return typecall(cls, systems)
        if isinstance(systems[0].base_ring(), FiniteField):
            return DynamicalSemigroup_projective_finite_field(systems)
        return DynamicalSemigroup_projective_field(systems)

    def dehomogenize(self, n):
        r"""
        Return a new :class:`DynamicalSemigroup_projective` with the dehomogenization at ``n`` of
        the generators of this dynamical semigroup.

        INPUT:

        - ``n`` -- a tuple of nonnegative integers. If ``n`` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: :class:`DynamicalSemigroup_affine`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x, y], P)
            sage: g = DynamicalSystem([x^2, y^2], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.dehomogenize(0)
            Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 2 dynamical systems:
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (y) to
                    (y)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (y) to
                    (y^2)

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x, y], P)
            sage: g = DynamicalSystem([x^2, y^2], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.dehomogenize((1, 0))
            Traceback (most recent call last):
            ...
            ValueError: Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x : y) dehomogenized at (1, 0) is not a `DynamicalSystem_affine` object
        """
        new_systems = []
        for ds in self._dynamical_systems:
            new_system = ds.dehomogenize(n)
            if not isinstance(new_system, DynamicalSystem_affine):
                raise ValueError(str(ds) + " dehomogenized at " + str(n) + " is not a `DynamicalSystem_affine` object")
            new_systems.append(new_system)
        return DynamicalSemigroup_affine(new_systems)

    def _mul_(self, other_dynamical_semigroup):
        r"""
        Return a new :class:`DynamicalSystem_projective` with???

        INPUT:

        - ``other_dynamical_semigroup`` -- a dynamical semigroup over projective space

        OUTPUT: :class:`DynamicalSemigroup_projective`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: f2 = DynamicalSystem([x^2, y^2], P)
            sage: g1 = DynamicalSystem([x^3, y^3], P)
            sage: g2 = DynamicalSystem([x^4, y^4], P)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 8 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 : y^3)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 : y^4)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^6 : y^6)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^8 : y^8)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 : y^3)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^6 : y^6)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 : y^4)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^8 : y^8)

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: f2 = DynamicalSystem([x^2, y^2], P)
            sage: g1 = DynamicalSystem([x^3, y^3], P)
            sage: g2 = DynamicalSystem([x^4, y^4], P)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g == g*f
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: f2 = DynamicalSystem([x^2, y^2], P)
            sage: g1 = DynamicalSystem([x^3, y^3], P)
            sage: g2 = DynamicalSystem([x^4, y^4], P)
            sage: h1 = DynamicalSystem([x^5, y^5], P)
            sage: h2 = DynamicalSystem([x^6, y^6], P)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: h = DynamicalSemigroup((h1, h2))
            sage: f*(g*h) == (f*g)*h
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: f2 = DynamicalSystem([x^2, y^2], P)
            sage: g1 = DynamicalSystem(z^3, A)
            sage: g2 = DynamicalSystem(z^4, A)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g
            Traceback (most recent call last):
            ...
            TypeError: can only multiply `DynamicalSemigroup_projective` objects

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: Q.<u, v, w> = ProjectiveSpace(QQ, 2)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: f2 = DynamicalSystem([x^2, y^2], P)
            sage: g1 = DynamicalSystem([u^3, v^3, w^3], Q)
            sage: g2 = DynamicalSystem([u^4, v^4, w^4], Q)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g
            Traceback (most recent call last):
            ...
            ValueError: cannot multiply dynamical semigroups of different dimensions
        """
        if not isinstance(other_dynamical_semigroup, DynamicalSemigroup_projective):
            raise TypeError("can only multiply `DynamicalSemigroup_projective` objects")
        if self._dynamical_systems[0].domain().dimension() != other_dynamical_semigroup._dynamical_systems[0].domain().dimension():
            raise ValueError("cannot multiply dynamical semigroups of different dimensions")
        composite_systems = []
        my_polys = self.defining_polynomials()
        other_polys = other_dynamical_semigroup.defining_polynomials()
        for my_poly in my_polys:
            for other_poly in other_polys:
                composite_poly = []
                for coordinate_poly in my_poly:
                    composite_poly.append(coordinate_poly(other_poly))
                composite_system = DynamicalSystem_projective(composite_poly)
                composite_systems.append(composite_system)
        for other_poly in other_polys:
            for my_poly in my_polys:
                composite_poly = []
                for coordinate_poly in other_poly:
                    composite_poly.append(coordinate_poly(my_poly))
                composite_system = DynamicalSystem_projective(composite_poly)
                composite_systems.append(composite_system)
        return DynamicalSemigroup_projective(composite_systems)

class DynamicalSemigroup_projective_field(DynamicalSemigroup_projective):
    pass

class DynamicalSemigroup_projective_finite_field(DynamicalSemigroup_projective_field):
    pass

class DynamicalSemigroup_affine(DynamicalSemigroup):
    r"""
    A dynamical semigroup defined by multiple dynamical systems on affine space.

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

    @staticmethod
    def __classcall_private__(cls, ds_data):
        systems = []

        if isinstance(ds_data, Collection):
            for ds_datum in ds_data:
                if isinstance(ds_datum, DynamicalSystem_affine):
                    systems.append(ds_datum)
                else:
                    try:
                        systems.append(DynamicalSystem_affine(ds_datum))
                    except ValueError:
                        raise ValueError(str(ds_datum) + " does not define a 'DynamicalSystem_affine' object")
        else:
            if isinstance(ds_data, DynamicalSystem_affine):
                    systems.append(ds_data)
            else:
                try:
                    systems.append(DynamicalSystem_affine(ds_data))
                except ValueError:
                    raise ValueError(str(ds_data) + " does not define a 'DynamicalSystem_affine' object")

        systems = _standardize_domains_of_(systems)
        if systems[0].base_ring() not in Fields():
            return typecall(cls, systems)
        if isinstance(systems[0].base_ring(), FiniteField):
            return DynamicalSemigroup_affine_finite_field(systems)
        return DynamicalSemigroup_affine_field(systems)

    def homogenize(self, n):
        r"""
        Return a new :class:`DynamicalSemigroup_projective` with the homogenization at ``n`` of
        the generators of this dynamical semigroup.

        INPUT:

        - ``n`` -- a tuple of nonnegative integers. If ``n`` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: :class:`DynamicalSemigroup_projective`

        EXAMPLES::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem(x + 1, A)
            sage: g = DynamicalSystem(x^2, A)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.homogenize(1)
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x0 + x1 : x1)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x0^2 : x1^2)

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem(x + 1, A)
            sage: g = DynamicalSystem(x^2, A)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.homogenize((1, 0))
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x1 : x0 + x1)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x1^2 : x0^2)
        """
        new_systems = []
        for ds in self._dynamical_systems:
            new_system = ds.homogenize(n)
            if not isinstance(new_system, DynamicalSystem_projective):
                raise ValueError(str(ds) + " homogenized at " + str(n) + " is not a `DynamicalSystem_projective` object")
            new_systems.append(new_system)
        return DynamicalSemigroup_projective(new_systems)

    def _mul_(self, other_dynamical_semigroup):
        r"""
        Return a new :class:`DynamicalSystem_affine` with???

        INPUT:

        - ``other_dynamical_semigroup`` -- a dynamical semigroup over affine space

        OUTPUT: :class:`DynamicalSemigroup_affine`

        EXAMPLES::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x, A)
            sage: f2 = DynamicalSystem(x^2, A)
            sage: g1 = DynamicalSystem(x^3, A)
            sage: g2 = DynamicalSystem(x^4, A)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g
            Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 8 dynamical systems:
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^3)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^4)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^6)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^8)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^3)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^6)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^4)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^8)

        TESTS::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x, A)
            sage: f2 = DynamicalSystem(x^2, A)
            sage: g1 = DynamicalSystem(x^3, A)
            sage: g2 = DynamicalSystem(x^4, A)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g == g*f
            True

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x, A)
            sage: f2 = DynamicalSystem(x^2, A)
            sage: g1 = DynamicalSystem(x^3, A)
            sage: g2 = DynamicalSystem(x^4, A)
            sage: h1 = DynamicalSystem(x^5, A)
            sage: h2 = DynamicalSystem(x^6, A)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: h = DynamicalSemigroup((h1, h2))
            sage: f*(g*h) == (f*g)*h
            True

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: P.<y,z> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x, A)
            sage: f2 = DynamicalSystem(x^2, A)
            sage: g1 = DynamicalSystem([y^3, z^3], P)
            sage: g2 = DynamicalSystem([y^4, z^4], P)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g
            Traceback (most recent call last):
            ...
            TypeError: can only multiply `DynamicalSemigroup_affine` objects

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: B.<y,z> = AffineSpace(QQ, 2)
            sage: f1 = DynamicalSystem(x, A)
            sage: f2 = DynamicalSystem(x^2, A)
            sage: g1 = DynamicalSystem([y^3, z^3], B)
            sage: g2 = DynamicalSystem([y^4, z^4], B)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g
            Traceback (most recent call last):
            ...
            ValueError: cannot multiply dynamical semigroups of different dimensions
        """
        if not isinstance(other_dynamical_semigroup, DynamicalSemigroup_affine):
            raise TypeError("can only multiply `DynamicalSemigroup_affine` objects")
        if self._dynamical_systems[0].domain().dimension() != other_dynamical_semigroup._dynamical_systems[0].domain().dimension():
            raise ValueError("cannot multiply dynamical semigroups of different dimensions")
        composite_systems = []
        my_polys = self.defining_polynomials()
        other_polys = other_dynamical_semigroup.defining_polynomials()
        for my_poly in my_polys:
            for other_poly in other_polys:
                composite_poly = []
                for coordinate_poly in my_poly:
                    composite_poly.append(coordinate_poly(other_poly))
                composite_system = DynamicalSystem_affine(composite_poly)
                composite_systems.append(composite_system)
        for other_poly in other_polys:
            for my_poly in my_polys:
                composite_poly = []
                for coordinate_poly in other_poly:
                    composite_poly.append(coordinate_poly(my_poly))
                composite_system = DynamicalSystem_affine(composite_poly)
                composite_systems.append(composite_system)
        return DynamicalSemigroup_affine(composite_systems)

class DynamicalSemigroup_affine_field(DynamicalSemigroup_affine):
    pass

class DynamicalSemigroup_affine_finite_field(DynamicalSemigroup_affine_field):
    pass

def _standardize_domains_of_(systems):
    r"""
    Coerces dynamical systems to the same domain and have the same generators.

    INPUT:

    - ``systems`` -- list of dynamical systems

    OUTPUT: list of dynamical systems from ``systems`` coerced to the same domain with the same generators

    EXAMPLES::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: B.<y> = AffineSpace(RR, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(x^2, B)
        sage: sage.dynamics.arithmetic_dynamics.dynamical_semigroup._standardize_domains_of_([f, g])
        [Dynamical System of Affine Space of dimension 1 over Real Field with 53 bits of precision
           Defn: Defined on coordinates by sending (x) to
             (x),
         Dynamical System of Affine Space of dimension 1 over Real Field with 53 bits of precision
           Defn: Defined on coordinates by sending (x) to
             (x^2)]
    """
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

    domain = systems[0].domain()

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
            systems[i] = DynamicalSystem(new_polys, domain)

    return systems