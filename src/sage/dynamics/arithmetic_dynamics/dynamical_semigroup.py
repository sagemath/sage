# sage_setup: distribution = sagemath-schemes
r"""
Dynamical semigroups

A dynamical semigroup is a finitely generated subsemigroup of
the endomorphism ring of a subscheme of projective or affine space.

AUTHORS:

 - Dang Phan (August 6th, 2023): initial implementation
"""

# ****************************************************************************
# Dang Phan <dang8phan@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

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
from sage.rings.integer_ring import ZZ
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
        Dynamical semigroup over Projective Space of dimension 1 over Rational Field
         defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([x^2, y^2], P)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Rational Field
         defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    ::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem_affine(x, A)
        sage: DynamicalSemigroup(f)
        Dynamical semigroup over Affine Space of dimension 1 over Rational Field
         defined by 1 dynamical system:
          Dynamical System of Affine Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x) to (x)

    ::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(x^2, A)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over Rational Field
         defined by 2 dynamical systems:
          Dynamical System of Affine Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x) to (x)
          Dynamical System of Affine Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x) to (x^2)

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: X = P.subscheme(x - y)
        sage: f = DynamicalSystem_projective([x, y], X)
        sage: g = DynamicalSystem_projective([x^2, y^2], X)
        sage: DynamicalSemigroup_projective([f, g])
        Dynamical semigroup over Closed subscheme of Projective Space of dimension 1
         over Rational Field defined by: x - y
         defined by 2 dynamical systems:
          Dynamical System of Closed subscheme of Projective Space of dimension 1
           over Rational Field defined by: x - y
             Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Closed subscheme of Projective Space of dimension 1
           over Rational Field defined by: x - y
             Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    If a dynamical semigroup is built from dynamical systems with different base rings,
    all systems will be coerced to the largest base ring::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: Q.<z,w> = ProjectiveSpace(RR, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over
         Real Field with 53 bits of precision defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1
           over Real Field with 53 bits of precision
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1
           over Real Field with 53 bits of precision
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    ::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: B.<y> = AffineSpace(RR, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(y^2, B)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over
         Real Field with 53 bits of precision defined by 2 dynamical systems:
          Dynamical System of Affine Space of dimension 1 over
           Real Field with 53 bits of precision
            Defn: Defined on coordinates by sending (x) to (x)
          Dynamical System of Affine Space of dimension 1 over
           Real Field with 53 bits of precision
            Defn: Defined on coordinates by sending (x) to (x^2)

    If a dynamical semigroup is built from dynamical systems over number fields, a composite number field is created
    and all systems will be coerced to it. This composite number field contains all of the initial number fields::

        sage: # needs sage.rings.number_field
        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: Q.<x,y> = ProjectiveSpace(K, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over
         Number Field in k with defining polynomial r^2 - 2 defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over
           Number Field in k with defining polynomial r^2 - 2
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1 over
           Number Field in k with defining polynomial r^2 - 2
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    ::

        sage: # needs sage.rings.number_field
        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: L.<l> = NumberField(r^2 - 3)
        sage: P.<x,y> = ProjectiveSpace(K, 1)
        sage: Q.<z,w> = ProjectiveSpace(L, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem([z^2, w^2], Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over
         Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
         defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over
           Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1 over
           Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    ::

        sage: # needs sage.rings.number_field
        sage: R.<r> = QQ[]
        sage: K.<k> = NumberField(r^2 - 2)
        sage: L.<l> = NumberField(r^2 - 3)
        sage: P.<x> = AffineSpace(K, 1)
        sage: Q.<y> = AffineSpace(L, 1)
        sage: f = DynamicalSystem(x, P)
        sage: g = DynamicalSystem(y^2, Q)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over
         Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
         defined by 2 dynamical systems:
          Dynamical System of Affine Space of dimension 1 over
           Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
            Defn: Defined on coordinates by sending (x) to (x)
          Dynamical System of Affine Space of dimension 1 over
           Number Field in kl with defining polynomial r^4 - 10*r^2 + 1
            Defn: Defined on coordinates by sending (x) to (x^2)

    A dynamical semigroup may contain dynamical systems over function fields::

        sage: R.<r> = QQ[]
        sage: P.<x,y> = ProjectiveSpace(R, 1)
        sage: f = DynamicalSystem([r * x, y], P)
        sage: g = DynamicalSystem([x, r * y], P)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Univariate
         Polynomial Ring in r over Rational Field defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over
           Univariate Polynomial Ring in r over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (r*x : y)
          Dynamical System of Projective Space of dimension 1 over
           Univariate Polynomial Ring in r over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x : r*y)

    ::

        sage: R.<r> = QQ[]
        sage: P.<x,y> = ProjectiveSpace(R, 1)
        sage: f = DynamicalSystem([r * x, y], P)
        sage: g = DynamicalSystem([x, y], P)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Univariate
         Polynomial Ring in r over Rational Field defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over
           Univariate Polynomial Ring in r over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (r*x : y)
          Dynamical System of Projective Space of dimension 1 over
           Univariate Polynomial Ring in r over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x : y)

    ::

        sage: R.<r,s> = QQ[]
        sage: P.<x,y> = ProjectiveSpace(R, 1)
        sage: f = DynamicalSystem([r * x, y], P)
        sage: g = DynamicalSystem([s * x, y], P)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Multivariate
         Polynomial Ring in r, s over Rational Field defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over
           Multivariate Polynomial Ring in r, s over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (r*x : y)
          Dynamical System of Projective Space of dimension 1 over
           Multivariate Polynomial Ring in r, s over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (s*x : y)

    ::

        sage: R.<r,s> = QQ[]
        sage: P.<x,y> = ProjectiveSpace(R, 1)
        sage: f = DynamicalSystem([r * x, s * y], P)
        sage: g = DynamicalSystem([s * x, r * y], P)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over
         Multivariate Polynomial Ring in r, s over Rational Field
         defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over
           Multivariate Polynomial Ring in r, s over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (r*x : s*y)
          Dynamical System of Projective Space of dimension 1 over
           Multivariate Polynomial Ring in r, s over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (s*x : r*y)

    A dynamical semigroup may contain dynamical systems over finite fields::

        sage: F = FiniteField(5)
        sage: P.<x,y> = ProjectiveSpace(F, 1)
        sage: DynamicalSemigroup(([x, y], [x^2, y^2]))
        Dynamical semigroup over Projective Space of dimension 1 over
         Finite Field of size 5 defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over Finite Field of size 5
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1 over Finite Field of size 5
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    If a dynamical semigroup is built from dynamical systems over both projective and
    affine spaces, all systems will be homogenized to dynamical systems over projective space::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: A.<z> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem([x, y], P)
        sage: g = DynamicalSystem(z^2, A)
        sage: DynamicalSemigroup((f, g))
        Dynamical semigroup over Projective Space of dimension 1 over Rational Field
         defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)

    TESTS::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: DynamicalSemigroup(1)
        Traceback (most recent call last):
        ...
        TypeError: 1 does not define a 'DynamicalSemigroup' object

    ::

        sage: # needs sage.rings.number_field
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
        elif isinstance(ds_data, DynamicalSystem_affine):
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

        self._dynamical_systems = []
        for ds in systems:
            if ds not in self._dynamical_systems:
                self._dynamical_systems.append(ds)
        Parent.__init__(self, category=Semigroups().FinitelyGeneratedAsMagma())

    def __call__(self, input):
        r"""
        The result after evaluating this dynamical semigroup on a value.

        INPUT:

        - ``input`` -- one value that can be evaluated
          with the generators of this dynamical semigroup

        OUTPUT: a set of the resulting values after applying all of this
        dynamical semigroup's generators to ``input``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f(2)
            {(2 : 1), (4 : 1)}

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f([2, 1])
            {(2 : 1), (4 : 1)}

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f(f(2))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert {(4 : 1), (2 : 1)} to an element of Rational Field
        """
        result = set()
        for ds in self.defining_systems():
            result.add(ds(self.domain()(input)))
        return result

    def base_ring(self):
        r"""
        The base ring of this dynamical semigroup. This is identical
        to the base ring of all of its defining dynamical system.

        OUTPUT: a ring

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.base_ring()
            Rational Field
        """
        return self.defining_systems()[0].base_ring()

    def change_ring(self, new_ring):
        r"""
        Return a new :class:`DynamicalSemigroup` whose generators
        are the initial dynamical systems coerced to ``new_ring``.

        INPUT:

        - ``new_ring`` -- a ring

        OUTPUT:

        A :class:`DynamicalSemigroup` defined by this dynamical
        semigroup's generators, but coerced to ``new_ring``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.change_ring(RR)
            Dynamical semigroup over Projective Space of dimension 1 over
             Real Field with 53 bits of precision defined by 2 dynamical systems:
              Dynamical System of Projective Space of dimension 1 over
               Real Field with 53 bits of precision
                Defn: Defined on coordinates by sending (x : y) to (x : y)
              Dynamical System of Projective Space of dimension 1 over
               Real Field with 53 bits of precision
                Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)
        """
        new_systems = []
        for ds in self.defining_systems():
            new_systems.append(ds.change_ring(new_ring))
        return DynamicalSemigroup_projective(new_systems)

    def domain(self):
        r"""
        Return the domain of the generators of this dynamical semigroup.

        OUTPUT: a subscheme of a projective space or affine space

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.domain()
            Projective Space of dimension 1 over Rational Field
        """
        return self.defining_systems()[0].domain()

    def codomain(self):
        r"""
        Return the codomain of the generators of this dynamical semigroup.

        OUTPUT: a subscheme of a projective space or affine space

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.codomain()
            Projective Space of dimension 1 over Rational Field
        """
        return self.defining_systems()[0].codomain()

    def defining_polynomials(self):
        r"""
        Return the set of polynomials that define the generators of this dynamical semigroup.

        OUTPUT: a set of polynomials

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.defining_polynomials()
            {(x, y), (x^2, y^2)}
        """
        result = set()
        for ds in self.defining_systems():
            result.add(ds.defining_polynomials())
        return result

    def defining_systems(self):
        r"""
        Return the generators of this dynamical semigroup.

        OUTPUT: a tuple of dynamical systems

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f.defining_systems()
            (Dynamical System of Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y) to (x : y),
             Dynamical System of Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2))
        """
        return tuple(self._dynamical_systems)

    def nth_iterate(self, p, n):
        r"""
        Return a set of values that results from evaluating this dynamical semigroup
        on the value ``p`` a total of ``n`` times.

        INPUT:

        - ``p`` -- a value on which dynamical systems can evaluate
        - ``n`` -- nonnegative integer

        OUTPUT: a set of values

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x^2, y^2],))
            sage: f.nth_iterate(2, 0)
            {(2 : 1)}

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x^2, y^2],))
            sage: f.nth_iterate(2, 1)
            {(4 : 1)}

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x^2, y^2],))
            sage: f.nth_iterate(2, 2)
            {(16 : 1)}

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: f.nth_iterate(2, 0)
            {(2 : 1)}

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: f.nth_iterate(2, 1)
            {(3 : 1), (4 : 1)}

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: f.nth_iterate(2, 2)
            {(5/3 : 1), (2 : 1), (9 : 1), (16 : 1)}

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: f.nth_iterate(2, 3.5)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: f.nth_iterate(2, -3)
            Traceback (most recent call last):
            ...
            ValueError: -3 must be a nonnegative integer

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: f.nth_iterate(3, 2) == (f * f)(3)
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x + y, x - y], [x^2, y^2]))
            sage: one = QQ(1)
            sage: f.nth_iterate(2, one)
            {(3 : 1), (4 : 1)}
        """
        n = ZZ(n)
        if n < 0:
            raise ValueError(str(n) + " must be a nonnegative integer")
        result = {self.domain()(p)}
        for i in range(1, n + 1):
            next_iteration = set()
            for point in result:
                next_iteration.update(self(point))
            result = next_iteration
        return result

    def orbit(self, p, n):
        r"""
        If ``n`` is an integer, return `(p, f(p), f^2(p), \dots, f^n(p))`. If ``n`` is a list or tuple in interval
        notation `[a, b]`, return `(f^a(p), \dots, f^b(p))`.

        INPUT:

        - ``p`` -- value on which this dynamical semigroup can be evaluated
        - ``n`` -- nonnegative integer or a list or tuple of length 2 describing an
          interval of the number line containing entirely nonnegative integers

        OUTPUT: a tuple of sets of values on the domain of this dynamical semigroup

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, 0)
            ({(2 : 1)},)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, 1)
            ({(2 : 1)}, {(2 : 1), (4 : 1)})

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, 2)
            ({(2 : 1)}, {(2 : 1), (4 : 1)}, {(2 : 1), (4 : 1), (16 : 1)})

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, [1, 2])
            ({(2 : 1), (4 : 1)}, {(2 : 1), (4 : 1), (16 : 1)})

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: one = QQ(1)
            sage: d.orbit(2, one)
            ({(2 : 1)}, {(2 : 1), (4 : 1)})

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, -2)
            Traceback (most recent call last):
            ...
            ValueError: -2 must be a nonnegative integer

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, x)
            Traceback (most recent call last):
            ...
            TypeError: x is not a constant polynomial

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, [1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: [1, 2, 3] must be an integer or list or tuple of two integers

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, [-2, 1])
            Traceback (most recent call last):
            ...
            ValueError: [-2, 1] must contain exactly two nonnegative integers

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: d = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: d.orbit(2, [2, 1])
            Traceback (most recent call last):
            ...
            ValueError: [2, 1] cannot be in descending order
        """

        if not isinstance(n, Collection):
            n = ZZ(n)
            if n < 0:
                raise ValueError(str(n) + " must be a nonnegative integer")
            return self.orbit(p, [0, n])

        if not len(n) == 2:
            raise ValueError(str(n) + " must be an integer or list or tuple of two integers")
        if ZZ(n[0]) < 0 or ZZ(n[1]) < 0:
            raise ValueError(str(n) + " must contain exactly two nonnegative integers")
        if ZZ(n[0]) > ZZ(n[1]):
            raise ValueError(str(n) + " cannot be in descending order")

        result = []
        current_iterate = self.nth_iterate(p, n[0])
        result.append(current_iterate)
        for i in range(n[0] + 1, n[1] + 1):
            next_iterate = set()
            for value in current_iterate:
                next_iterate.update(self(value))
            result.append(next_iterate)
            current_iterate = next_iterate
        return tuple(result)

    def specialization(self, assignments):
        r"""
        Return the specialization of the generators of this dynamical semigroup.

        INPUT:

        - ``assignments`` -- argument for specialization of the generators of
          this dynamical semigroup

        OUTPUT: a dynamical semigroup with the specialization of the generators
        of this dynamical semigroup

        EXAMPLES::

            sage: R.<r> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem([r * x, y], P)
            sage: g = DynamicalSystem([x, r * y], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.specialization({r:2})
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field
             defined by 2 dynamical systems:
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (2*x : y)
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (x : 2*y)

        ::

            sage: R.<r> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem([r * x, y], P)
            sage: g = DynamicalSystem([x, y], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.specialization({r:2})
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field
             defined by 2 dynamical systems:
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (2*x : y)
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (x : y)

        ::

            sage: R.<r,s> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem([r * x, y], P)
            sage: g = DynamicalSystem([s * x, y], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.specialization({r:2, s:3})
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field
             defined by 2 dynamical systems:
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (2*x : y)
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (3*x : y)

        ::

            sage: R.<r,s> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem([r * x, s * y], P)
            sage: g = DynamicalSystem([s * x, r * y], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.specialization({s:3})
            Dynamical semigroup over Projective Space of dimension 1 over
             Univariate Polynomial Ring in r over Rational Field
             defined by 2 dynamical systems:
              Dynamical System of Projective Space of dimension 1 over
               Univariate Polynomial Ring in r over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (r*x : 3*y)
              Dynamical System of Projective Space of dimension 1 over
               Univariate Polynomial Ring in r over Rational Field
                Defn: Defined on coordinates by sending (x : y) to (3*x : r*y)
        """
        specialized_systems = []
        for ds in self.defining_systems():
            specialized_systems.append(ds.specialization(assignments))
        return DynamicalSemigroup(specialized_systems)

    def __mul__(self, other_dynamical_semigroup):
        r"""
        Return a new :class:`DynamicalSemigroup` that is the result of multiplying
        this dynamical semigroup with another dynamical semigroup of the same type
        using the * operator.

        Let `f` be a dynamical semigroup with generators `\{ f_1, f_2, \dots, f_m \}`
        and `g` be a dynamical semigroup with generators `\{ g_1, g_2, \dots, g_n \}`.
        The product `f * g` has generators
        `\{ f_i(g_j) : 1 \leq i \leq m, 1 \leq j \leq n \}`.

        INPUT:

        - ``other_dynamical_semigroup`` -- a dynamical semigroup

        OUTPUT: :class:`DynamicalSemigroup`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x^10, y^10], P)
            sage: f = DynamicalSemigroup(f1)
            sage: f*f
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 1 dynamical system:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^100 : y^100)

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem_affine(x^10, A)
            sage: f = DynamicalSemigroup(f1)
            sage: f*f
            Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 1 dynamical system:
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^100)

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x^2 + x * y, y^2 + x * y], P)
            sage: g1 = DynamicalSystem([x^3 + x^2 * y, y^3 + x * y^2], P)
            sage: f = DynamicalSemigroup(f1)
            sage: g = DynamicalSemigroup(g1)
            sage: f*g
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 1 dynamical system:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^6 + 2*x^5*y + 2*x^4*y^2 + 2*x^3*y^3 + x^2*y^4 : x^4*y^2 + 2*x^3*y^3 + 2*x^2*y^4 + 2*x*y^5 + y^6)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x^2 + x * y, y^2 + x * y], P)
            sage: f2 = DynamicalSystem([x^2 - x * y, y^2 - x * y], P)
            sage: g1 = DynamicalSystem([x^3 + x^2 * y, y^3 + x * y^2], P)
            sage: g2 = DynamicalSystem([x^3 - x^2 * y, y^3 - x * y^2], P)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f*g
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field defined by 2 dynamical systems:
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^6 + 2*x^5*y + 2*x^4*y^2 + 2*x^3*y^3 + x^2*y^4 : x^4*y^2 + 2*x^3*y^3 + 2*x^2*y^4 + 2*x*y^5 + y^6)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^6 - 2*x^5*y + 2*x^3*y^3 - x^2*y^4 : -x^4*y^2 + 2*x^3*y^3 - 2*x*y^5 + y^6)

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x^2 + x * y, y^2 - x * y], P)
            sage: g1 = DynamicalSystem([x^3 - x^2 * y, y^3 + x * y^2], P)
            sage: f = DynamicalSemigroup(f1)
            sage: g = DynamicalSemigroup(g1)
            sage: f*g == g*f
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x^2 + x * y, y^2 + x * y], P)
            sage: g1 = DynamicalSystem([x^3 + x^2 * y, y^3 + x * y^2], P)
            sage: h1 = DynamicalSystem([x^4 + x^3 * y, y^4 + x * y^3], P)
            sage: f = DynamicalSemigroup(f1)
            sage: g = DynamicalSemigroup(g1)
            sage: h = DynamicalSemigroup(h1)
            sage: (f*g)*h == f*(g*h)
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: f = DynamicalSemigroup(f1)
            sage: f*1
            Traceback (most recent call last):
            ...
            TypeError: can only multiply dynamical semigroups with other dynamical semigroups of the same type

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem_projective([x, y], P)
            sage: g1 = DynamicalSystem_affine(z, A)
            sage: f = DynamicalSemigroup_projective(f1)
            sage: g = DynamicalSemigroup_affine(g1)
            sage: f*g
            Traceback (most recent call last):
            ...
            TypeError: can only multiply dynamical semigroups with other dynamical semigroups of the same type

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: Q.<z,w> = ProjectiveSpace(RR, 1)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: g1 = DynamicalSystem([z, w], Q)
            sage: f = DynamicalSemigroup(f1)
            sage: g = DynamicalSemigroup(g1)
            sage: f*g
            Traceback (most recent call last):
            ...
            ValueError: left dynamical semigroup's domain must equal right dynamical semigroup's codomain

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: Q.<z,w> = ProjectiveSpace(QQ, 1)
            sage: f1 = DynamicalSystem([x, y], P)
            sage: g1 = DynamicalSystem([z, w], Q)
            sage: f = DynamicalSemigroup(f1)
            sage: g = DynamicalSemigroup(g1)
            sage: f*g
            Traceback (most recent call last):
            ...
            ValueError: left dynamical semigroup's domain must equal right dynamical semigroup's codomain
        """
        if type(self) is not type(other_dynamical_semigroup):
            raise TypeError("can only multiply dynamical semigroups with other dynamical semigroups of the same type")
        if self.domain() != other_dynamical_semigroup.codomain():
            raise ValueError("left dynamical semigroup's domain must equal right dynamical semigroup's codomain")
        composite_systems = []
        for f in self.defining_systems():
            for g in other_dynamical_semigroup.defining_systems():
                composite_systems.append(DynamicalSystem(f * g))
        return DynamicalSemigroup(composite_systems)

    def __pow__(self, n):
        r"""
        Return a new dynamical semigroup that is this dynamical semigroup with itself ``n`` times.
        If ``n`` is zero, return a dynamical semigroup with the identity map.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: :class:`DynamicalSemigroup`

        EXAMPLES::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem(x^2, A)
            sage: d = DynamicalSemigroup(f)
            sage: d^2
            Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 1 dynamical system:
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^4)

        ::

            sage: A.<x, y, z, w, t, u ,v> = AffineSpace(QQ, 7)
            sage: f = DynamicalSystem([t + u, v - w, x + y, z^2, u * t, v^2 - w^2, x * y * z], A)
            sage: d = DynamicalSemigroup(f)
            sage: d^0
            Dynamical semigroup over Affine Space of dimension 7 over Rational Field defined by 1 dynamical system:
            Dynamical System of Affine Space of dimension 7 over Rational Field
              Defn: Defined on coordinates by sending (x, y, z, w, t, u, v) to
                    (x, y, z, w, t, u, v)

        TESTS::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: d = DynamicalSemigroup(f)
            sage: d*d == d^2
            True

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: d = DynamicalSemigroup(f)
            sage: d^3 * d^2 == d^(3 + 2)
            True

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: d = DynamicalSemigroup(f)
            sage: (d^3)^2 == d^(3 * 2)
            True

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: g1 = DynamicalSystem(x^3 + x - 1, A)
            sage: g2 = DynamicalSystem(x^2 - x + 1, A)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: (f * g) ^ 2 == f^2 * g^2
            False

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: g1 = DynamicalSystem(x^3 + x - 1, A)
            sage: g2 = DynamicalSystem(x^2 - x + 1, A)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: f * g^0 == f
            True

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: g1 = DynamicalSystem(x^3 + x - 1, A)
            sage: g2 = DynamicalSystem(x^2 - x + 1, A)
            sage: f = DynamicalSemigroup((f1, f2))
            sage: g = DynamicalSemigroup((g1, g2))
            sage: g * f^0 == g
            True

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: d = DynamicalSemigroup((f1, f2))
            sage: d^1.5
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f1 = DynamicalSystem(x^2 + 1, A)
            sage: f2 = DynamicalSystem(x^3 - 1, A)
            sage: d = DynamicalSemigroup((f1, f2))
            sage: d^(-1)
            Traceback (most recent call last):
            ...
            ValueError: -1 must be a nonnegative integer

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem(x^2, A)
            sage: d = DynamicalSemigroup(f)
            sage: two = RR(2)
            sage: d^two
            Dynamical semigroup over Affine Space of dimension 1 over Rational Field defined by 1 dynamical system:
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    (x^4)
        """
        n = ZZ(n)
        if n < 0:
            raise ValueError(str(n) + " must be a nonnegative integer")
        if n == 0:
            return DynamicalSemigroup(DynamicalSystem(self.defining_systems()[0] ** 0))
        result = self
        for i in range(n - 1):
            result = result * self
        return result

    def _repr_(self):
        r"""
        Return the :class:`String` representation of this dynamical semigroup.

        OUTPUT: a :class:`String` displaying information about this dynamical semigroup

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
        if (len(self.defining_systems()) > 1):
            header += "s"
        header += ":"
        header = header % (str(self.domain()), len(self.defining_systems()))
        systems = []
        for ds in self.defining_systems():
            systems.append(str(ds))
        systems = '\n'.join(systems)
        return header + "\n" + systems

    def __eq__(self, other):
        r"""
        Return whether two dynamical semigroups are equal.

        OUTPUT:

        A boolean that is ``True`` if and only if the generators of the two
        dynamical semigroups are equal as sets and no generator is of degree 1.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x^2, y^2], [x^3, y^3]))
            sage: g = DynamicalSemigroup(([x^2, y^2], [x^3, y^3]))
            sage: f == g
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x^2, y^2], [x^2, y^2]))
            sage: g = DynamicalSemigroup(([x^2, y^2], [x^2, y^2], [x^2, y^2]))
            sage: f == g
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x^3, y^3], [x^2, y^2]))
            sage: g = DynamicalSemigroup(([x^2, y^2], [x^3, y^3]))
            sage: f == g
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x^3, y^3], [x^2, y^2]))
            sage: g = DynamicalSemigroup(([x^2, y^2], [y^3, x^3]))
            sage: f == g
            False

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: f == 1
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSemigroup(([x, y], [x^2, y^2]))
            sage: g = DynamicalSemigroup(([x^2, y^2], [x^3, y^3]))
            sage: f == g
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot compare dynamical semigroups with at least one generator of degree 1
        """
        if isinstance(other, DynamicalSemigroup):
            if any(ds.degree() == 1 for ds in self.defining_systems()) or \
                    any(ds.degree() == 1 for ds in other.defining_systems()):
                raise NotImplementedError("cannot compare dynamical semigroups with at least one generator of degree 1")
            return all(ds in other.defining_systems() for ds in self.defining_systems()) and \
                all(ds in self.defining_systems() for ds in other.defining_systems())
        return False


class DynamicalSemigroup_projective(DynamicalSemigroup):
    r"""
    A dynamical semigroup defined by a multiple dynamical systems on projective space.

    INPUT:

    - ``ds_data`` -- list or tuple of dynamical systems or objects that define
      dynamical systems over projective space

    OUTPUT: :class:`DynamicalSemigroup_projective`

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: DynamicalSemigroup_projective(([x, y], [x^2, y^2]))
        Dynamical semigroup over Projective Space of dimension 1 over
         Rational Field defined by 2 dynamical systems:
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x : y)
          Dynamical System of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to (x^2 : y^2)
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

        - ``n`` -- tuple of nonnegative integers; if `n` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: :class:`DynamicalSemigroup_affine`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x, y], P)
            sage: g = DynamicalSystem([x^2, y^2], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.dehomogenize(0)
            Dynamical semigroup over Affine Space of dimension 1 over
             Rational Field defined by 2 dynamical systems:
              Dynamical System of Affine Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (y) to (y)
              Dynamical System of Affine Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (y) to (y^2)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x, y], P)
            sage: g = DynamicalSystem([x^2, y^2], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.dehomogenize(1)
            Dynamical semigroup over Affine Space of dimension 1 over
             Rational Field defined by 2 dynamical systems:
              Dynamical System of Affine Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x) to (x)
              Dynamical System of Affine Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x) to (x^2)

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x, y], P)
            sage: g = DynamicalSystem([x^2, y^2], P)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.dehomogenize((1, 0))
            Traceback (most recent call last):
            ...
            ValueError: Scheme morphism:
              From: Affine Space of dimension 1 over Rational Field
              To:   Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to (1/x)
            is not a `DynamicalSystem_affine` object
        """
        new_systems = []
        for ds in self.defining_systems():
            new_system = ds.dehomogenize(n)
            if not isinstance(new_system, DynamicalSystem_affine):
                raise ValueError(str(new_system) + " is not a `DynamicalSystem_affine` object")
            new_systems.append(new_system)
        return DynamicalSemigroup_affine(new_systems)


class DynamicalSemigroup_projective_field(DynamicalSemigroup_projective):
    pass


class DynamicalSemigroup_projective_finite_field(DynamicalSemigroup_projective_field):
    pass


class DynamicalSemigroup_affine(DynamicalSemigroup):
    r"""
    A dynamical semigroup defined by multiple dynamical systems on affine space.

    INPUT:

    - ``ds_data`` -- list or tuple of dynamical systems or objects that define dynamical systems
      over affine space

    OUTPUT: :class:`DynamicalSemigroup_affine`

    EXAMPLES::

        sage: A.<x> = AffineSpace(QQ, 1)
        sage: f = DynamicalSystem(x, A)
        sage: g = DynamicalSystem(x^2, A)
        sage: DynamicalSemigroup_affine((f, g))
        Dynamical semigroup over Affine Space of dimension 1 over
         Rational Field defined by 2 dynamical systems:
          Dynamical System of Affine Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x) to (x)
          Dynamical System of Affine Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x) to (x^2)
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

        - ``n`` -- tuple of nonnegative integers; if `n` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: :class:`DynamicalSemigroup_projective`

        EXAMPLES::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem(x + 1, A)
            sage: g = DynamicalSystem(x^2, A)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.homogenize(1)
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field
             defined by 2 dynamical systems:
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x0 : x1) to (x0 + x1 : x1)
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x0 : x1) to (x0^2 : x1^2)

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem(x + 1, A)
            sage: g = DynamicalSystem(x^2, A)
            sage: d = DynamicalSemigroup((f, g))
            sage: d.homogenize((1, 0))
            Dynamical semigroup over Projective Space of dimension 1 over Rational Field
             defined by 2 dynamical systems:
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x0 : x1) to (x1 : x0 + x1)
              Dynamical System of Projective Space of dimension 1 over Rational Field
                Defn: Defined on coordinates by sending (x0 : x1) to (x1^2 : x0^2)
        """
        new_systems = []
        for ds in self.defining_systems():
            new_systems.append(ds.homogenize(n))
        return DynamicalSemigroup_projective(new_systems)


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
           Defn: Defined on coordinates by sending (x) to (x),
         Dynamical System of Affine Space of dimension 1 over Real Field with 53 bits of precision
           Defn: Defined on coordinates by sending (x) to (x^2)]
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
                raise ValueError("given dynamical systems are not automorphic \
                                under global composition")

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
