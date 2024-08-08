r"""
Derivations

Let `A` be a ring and `B` be a bimodule over `A`.
A derivation `d : A \to B` is an additive map that satisfies
the Leibniz rule

.. MATH::

    d(xy) = x d(y) + d(x) y.

If `B` is an algebra over `A` and if we are given in addition a
ring homomorphism `\theta : A \to B`, a twisted derivation with respect
to `\theta` (or a `\theta`-derivation) is an additive map `d : A \to B`
such that

.. MATH::

    d(xy) = \theta(x) d(y) + d(x) y.

When `\theta` is the morphism defining the structure of `A`-algebra
on `B`, a `\theta`-derivation is nothing but a derivation.
In general, if `\iota : A \to B` denotes the defining morphism above,
one easily checks that `\theta - \iota` is a `\theta`-derivation.

This file provides support for derivations and twisted derivations
over commutative rings with values in algebras (i.e. we require
that `B` is a commutative `A`-algebra).
In this case, the set of derivations (resp. `\theta`-derivations)
is a module over `B`.

Given a ring `A`, the module of derivations over `A` can be created
as follows::

    sage: A.<x,y,z> = QQ[]
    sage: M = A.derivation_module()
    sage: M
    Module of derivations over
     Multivariate Polynomial Ring in x, y, z over Rational Field

The method :meth:`~sage.rings.derivation.RingDerivationModule.gens`
returns the generators of this module::

    sage: A.<x,y,z> = QQ[]
    sage: M = A.derivation_module()
    sage: M.gens()
    (d/dx, d/dy, d/dz)

We can combine them in order to create all derivations::

    sage: d = 2*M.gen(0) + z*M.gen(1) + (x^2 + y^2)*M.gen(2)
    sage: d
    2*d/dx + z*d/dy + (x^2 + y^2)*d/dz

and now play with them::

    sage: d(x + y + z)
    x^2 + y^2 + z + 2
    sage: P = A.random_element()
    sage: Q = A.random_element()
    sage: d(P*Q) == P*d(Q) + d(P)*Q
    True

Alternatively we can use the method
:meth:`~sage.rings.ring.CommutativeRing.derivation`
of the ring `A` to create derivations::

    sage: Dx = A.derivation(x); Dx
    d/dx
    sage: Dy = A.derivation(y); Dy
    d/dy
    sage: Dz = A.derivation(z); Dz
    d/dz
    sage: A.derivation([2, z, x^2+y^2])
    2*d/dx + z*d/dy + (x^2 + y^2)*d/dz

Sage knows moreover that `M` is a Lie algebra::

    sage: M.category()
    Join of
     Category of Lie algebras with basis over Rational Field and
     Category of modules with basis over
      Multivariate Polynomial Ring in x, y, z over Rational Field

Computations of Lie brackets are implemented as well::

    sage: Dx.bracket(Dy)
    0
    sage: d.bracket(Dx)
    -2*x*d/dz

At the creation of a module of derivations, a codomain can be specified::

    sage: B = A.fraction_field()
    sage: A.derivation_module(B)
    Module of derivations from Multivariate Polynomial Ring in x, y, z over Rational Field
     to Fraction Field of Multivariate Polynomial Ring in x, y, z over Rational Field

Alternatively, one can specify a morphism `f` with domain `A`.
In this case, the codomain of the derivations is the codomain of
`f` but the latter is viewed as an algebra over `A` through the
homomorphism `f`.
This construction is useful, for example, if we want to work with
derivations on `A` at a certain point, e.g. `(0,1,2)`. Indeed,
in order to achieve this, we first define the evaluation map at
this point::

    sage: ev = A.hom([QQ(0), QQ(1), QQ(2)])
    sage: ev
    Ring morphism:
      From: Multivariate Polynomial Ring in x, y, z over Rational Field
      To:   Rational Field
      Defn: x |--> 0
            y |--> 1
            z |--> 2

Now we use this ring homomorphism to define a structure of `A`-algebra
on `\QQ` and then build the following module of derivations::

    sage: M = A.derivation_module(ev)
    sage: M
    Module of derivations
     from Multivariate Polynomial Ring in x, y, z over Rational Field
       to Rational Field
    sage: M.gens()
    (d/dx, d/dy, d/dz)

Elements in `M` then acts as derivations at `(0,1,2)`::

    sage: Dx = M.gen(0)
    sage: Dy = M.gen(1)
    sage: Dz = M.gen(2)
    sage: f = x^2 + y^2 + z^2
    sage: Dx(f)  # = 2*x evaluated at (0,1,2)
    0
    sage: Dy(f)  # = 2*y evaluated at (0,1,2)
    2
    sage: Dz(f)  # = 2*z evaluated at (0,1,2)
    4


Twisted derivations are handled similarly::

    sage: theta = B.hom([B(y),B(z),B(x)])
    sage: theta
    Ring endomorphism of Fraction Field of
     Multivariate Polynomial Ring in x, y, z over Rational Field
      Defn: x |--> y
            y |--> z
            z |--> x

    sage: M = B.derivation_module(twist=theta)
    sage: M
    Module of twisted derivations over Fraction Field of Multivariate Polynomial Ring
     in x, y, z over Rational Field (twisting morphism: x |--> y, y |--> z, z |--> x)

Over a field, one proves that every `\theta`-derivation is a multiple
of `\theta - id`, so that::

    sage: d = M.gen(); d
    [x |--> y, y |--> z, z |--> x] - id

and then::

    sage: d(x)
    -x + y
    sage: d(y)
    -y + z
    sage: d(z)
    x - z
    sage: d(x + y + z)
    0

AUTHOR:

- Xavier Caruso (2018-09)
"""

# ***************************************************************************
#    Copyright (C) 2018 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family
from sage.modules.module import Module
from sage.structure.element import ModuleElement
from sage.rings.integer_ring import ZZ

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.power_series_ring import PowerSeriesRing_generic
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.fraction_field import FractionField_generic
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.function_field.function_field import FunctionField
from sage.rings.function_field.function_field_rational import RationalFunctionField
from sage.categories.number_fields import NumberFields
from sage.categories.finite_fields import FiniteFields
from sage.categories.modules import Modules
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.lie_algebras import LieAlgebras

from sage.categories.map import Map
from sage.categories.rings import Rings

from sage.misc.latex import latex


class RingDerivationModule(Module, UniqueRepresentation):
    """
    A class for modules of derivations over a commutative ring.
    """
    def __init__(self, domain, codomain, twist=None):
        """
        Initialize this module of derivation.

        TESTS::

            sage: A.<x,y> = QQ[]
            sage: M = A.derivation_module()
            sage: TestSuite(M).run()

            sage: from sage.rings.derivation import RingDerivationModule
            sage: R5.<x> = GF(5)[]
            sage: R25.<x> = GF(25)[]                                                    # needs sage.rings.finite_rings
            sage: R7.<x> = GF(7)[]

            sage: RingDerivationModule(R5, R25)                                         # needs sage.rings.finite_rings
            Module of derivations
             from Univariate Polynomial Ring in x over Finite Field of size 5
               to Univariate Polynomial Ring in x over Finite Field in z2 of size 5^2
            sage: RingDerivationModule(R5, R5^2)
            Traceback (most recent call last):
            ...
            TypeError: the codomain must be an algebra over the domain
             or a morphism with the correct domain
            sage: RingDerivationModule(R5, R7)                                          # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: the codomain must be an algebra over the domain
             or a morphism with the correct domain

            sage: theta = R5.hom([R5.gen()^2])
            sage: RingDerivationModule(R5, R25, twist=theta)                            # needs sage.rings.finite_rings
            Module of twisted derivations
             from Univariate Polynomial Ring in x over Finite Field of size 5
               to Univariate Polynomial Ring in x over Finite Field in z2 of size 5^2
             (twisting morphism: x |--> x^2)
            sage: RingDerivationModule(R7, R7, twist=theta)                             # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: the domain of the derivation must coerce to the domain of the twisting homomorphism
        """
        if domain not in Rings().Commutative():
            raise TypeError("the domain must be a commutative ring")

        if codomain in Rings().Commutative() and codomain.has_coerce_map_from(domain):
            defining_morphism = codomain.coerce_map_from(domain)
        elif (isinstance(codomain,Map)
              and codomain.category_for().is_subcategory(Rings())
              and codomain.domain().has_coerce_map_from(domain)):
            if codomain.domain() is domain:
                defining_morphism = codomain
            else:
                defining_morphism = codomain * codomain.domain().coerce_map_from(domain)
            codomain = defining_morphism.codomain()
        else:
            raise TypeError("the codomain must be an algebra over the domain"
                            " or a morphism with the correct domain")

        if twist is not None:
            if not (isinstance(twist, Map) and twist.category_for().is_subcategory(Rings())):
                raise TypeError("the twisting homomorphism must be a homomorphism of rings")
            if twist.domain() is not domain:
                map = twist.domain().coerce_map_from(domain)
                if map is None:
                    raise TypeError("the domain of the derivation must coerce"
                                    " to the domain of the twisting homomorphism")
                twist = twist * map
            if twist.codomain() is not codomain:
                map = codomain.coerce_map_from(twist.codomain())
                if map is None:
                    raise TypeError("the codomain of the twisting homomorphism"
                                    " must coerce to the codomain of the derivation")
                twist = map * twist
            # We check if the twisting morphism is the defining morphism
            try:
                if twist == defining_morphism:
                    twist = None
                else:
                    for g in domain.gens():
                        if twist(g) != defining_morphism(g):
                            break
                    else:
                        twist = None
            except (AttributeError, NotImplementedError):
                pass

        self._domain = domain
        self._codomain = codomain
        self._defining_morphism = defining_morphism
        self._twist = twist
        self._base_derivation = None
        self._gens = None
        self._basis = self._dual_basis = None
        # Currently basis and gens play exactly the same role because
        # the only rings that are supported lead to free modules of derivations
        # So the code is a bit redundant but we except to be able to cover more
        # rings (with non free modules of derivations) in a near future
        self._constants = (ZZ, False)
        if twist is not None:
            self.Element = RingDerivationWithTwist_generic
            if domain.is_field():
                self._gens = [ 1 ]
                self._basis = [ 1 ]
        elif (domain is ZZ or domain in NumberFields() or domain in FiniteFields()
              or isinstance(domain, IntegerModRing_generic)
              or (isinstance(domain, pAdicGeneric) and (domain.is_field() or domain.absolute_e() == 1))):
            self.Element = RingDerivationWithoutTwist_zero
            self._gens = [ ]
            self._basis = [ ]
            self._dual_basis = [ ]
            self._constants = (domain, True)
        elif (isinstance(domain, (PolynomialRing_general, MPolynomialRing_base, PowerSeriesRing_generic, LaurentSeriesRing))
              or (isinstance(domain, FractionField_generic)
                  and isinstance(domain.ring(), (PolynomialRing_general, MPolynomialRing_base)))):
            self._base_derivation = RingDerivationModule(domain.base_ring(), defining_morphism)
            self.Element = RingDerivationWithoutTwist_function
            try:
                self._gens = self._base_derivation.gens() + domain.gens()
            except NotImplementedError:
                pass
            try:
                self._basis = tuple(self._base_derivation.basis()) + domain.gens()
                self._dual_basis = tuple(self._base_derivation.dual_basis()) + domain.gens()
            except NotImplementedError:
                pass
            constants, sharp = self._base_derivation._constants
            if domain.characteristic() == 0:
                self._constants = (constants, sharp)
            else:
                # in this case, the constants are polynomials in x^p
                # TODO: implement this
                self._constants = (constants, False)
        elif isinstance(domain, FractionField_generic):
            self._base_derivation = RingDerivationModule(domain.ring(), defining_morphism)
            self.Element = RingDerivationWithoutTwist_fraction_field
            try:
                self._gens = self._base_derivation.gens()
            except NotImplementedError:
                pass
            try:
                self._basis = self._base_derivation.basis()
                self._dual_basis = self._base_derivation.dual_basis()
            except NotImplementedError:
                pass
            constants, sharp = self._base_derivation._constants
            self._constants = (constants.fraction_field(), False)
        elif isinstance(domain, PolynomialQuotientRing_generic):
            self._base_derivation = RingDerivationModule(domain.base(), defining_morphism)
            modulus = domain.modulus()
            for der in self._base_derivation.gens():
                if der(modulus) != 0:
                    raise NotImplementedError("derivations over quotient rings"
                                              " are not fully supported")
            self.Element = RingDerivationWithoutTwist_quotient
            try:
                self._gens = self._base_derivation.gens()
            except NotImplementedError:
                pass
            try:
                self._basis = self._base_derivation.basis()
                self._dual_basis = self._base_derivation.dual_basis()
            except NotImplementedError:
                pass
            constants, sharp = self._base_derivation._constants
            self._constants = (constants, False)  # can we do better?
        elif isinstance(domain, QuotientRing_generic):
            self._base_derivation = RingDerivationModule(domain.cover_ring(), defining_morphism)
            if any(der(modulus) != 0 for modulus in domain.defining_ideal().gens()
                   for der in self._base_derivation.gens()):
                raise NotImplementedError("derivations over quotient rings"
                                          " are not fully supported")
            self.Element = RingDerivationWithoutTwist_quotient
            try:
                self._gens = self._base_derivation.gens()
            except NotImplementedError:
                pass
            try:
                self._basis = self._base_derivation.basis()
                self._dual_basis = self._base_derivation.dual_basis()
            except NotImplementedError:
                pass
            constants, sharp = self._base_derivation._constants
            self._constants = (constants, False)  # can we do better?
        elif isinstance(domain, RationalFunctionField):
            from sage.rings.function_field.derivations_rational import FunctionFieldDerivation_rational
            self.Element = FunctionFieldDerivation_rational
            self._gens = self._basis = [ None ]
            self._dual_basis = [ domain.gen() ]
        elif isinstance(domain, FunctionField):
            if domain.is_separable():
                from sage.rings.function_field.derivations_polymod import FunctionFieldDerivation_separable
                self._base_derivation = RingDerivationModule(domain.base_ring(), defining_morphism)
                self.Element = FunctionFieldDerivation_separable
                try:
                    self._gens = self._base_derivation.gens()
                except NotImplementedError:
                    pass
                try:
                    self._basis = self._base_derivation.basis()
                    self._dual_basis = self._base_derivation.dual_basis()
                except NotImplementedError:
                    pass
            else:
                from sage.rings.function_field.derivations_polymod import FunctionFieldDerivation_inseparable
                M, f, self._t = domain.separable_model()
                self._base_derivation = RingDerivationModule(M, defining_morphism * f)
                self._d = self._base_derivation(None)
                self.Element = FunctionFieldDerivation_inseparable
                self._gens = self._basis = [ None ]
                self._dual_basis = [ f(M.base_ring().gen()) ]
        else:
            raise NotImplementedError("derivations over this ring is not implemented")
        if self._basis is None:
            category = Modules(codomain)
        else:
            category = ModulesWithBasis(codomain)
        if self._twist is None and domain is codomain:
            category &= LieAlgebras(self._constants[0])
        Module.__init__(self, codomain, category=category)
        if self._gens is not None:
            self._gens = [self.element_class(self, x) for x in self._gens]
        if self._basis is not None:
            self._basis = [self.element_class(self, x) for x in self._basis]
        if self._dual_basis is not None:
            self._dual_basis = [domain(x) for x in self._dual_basis]

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module()
            sage: hash(M) == hash((M.domain(), M.codomain(), M.twisting_morphism()))
            True
        """
        return hash((self._domain, self._codomain, self._twist))

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion map from ``R``
        to this module.

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: B.<y> = A[]
            sage: M1 = A.derivation_module(); M1
            Module of derivations over Univariate Polynomial Ring in x over Rational Field
            sage: M2 = A.derivation_module(B); M2
            Module of derivations
             from Univariate Polynomial Ring in x over Rational Field
               to Univariate Polynomial Ring in y over
                   Univariate Polynomial Ring in x over Rational Field
            sage: M1._coerce_map_from_(M2) is None
            True
            sage: M1.has_coerce_map_from(M2)
            False
            sage: M2.has_coerce_map_from(M1)
            True

            sage: M1.has_coerce_map_from(ZZ)
            False
            sage: M1.has_coerce_map_from(QQ)
            False
            sage: M1.has_coerce_map_from(A)
            False
        """
        if isinstance(R, RingDerivationModule):
            if R.domain().has_coerce_map_from(self._domain) and self._codomain.has_coerce_map_from(R.codomain()):
                morR = R.defining_morphism()
                morS = self._defining_morphism
                try:
                    # this test is not perfect
                    for g in self._domain.gens():
                        if morR(g) != morS(g):
                            return False
                    return True
                except (AttributeError, NotImplementedError):
                    pass
        return super()._coerce_map_from_(R)

    def _repr_(self):
        """
        Return a string representation of this module of derivations.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: R.derivation_module()
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring

            sage: theta = R.hom([y,x])
            sage: R.derivation_module(twist=theta)
            Module of twisted derivations over Multivariate Polynomial Ring in x, y
             over Integer Ring (twisting morphism: x |--> y, y |--> x)
        """
        t = ""
        if self._twist is None:
            s = "Module of derivations"
        else:
            s = "Module of twisted derivations"
            try:
                t = " (twisting morphism: %s)" % self._twist._repr_short()
            except AttributeError:
                pass
        if self._domain is self._codomain:
            s += " over %s" % self._domain
        else:
            s += " from %s to %s" % (self._domain, self._codomain)
        return s + t

    def domain(self):
        """
        Return the domain of the derivations in this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.domain()
            Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of the derivations in this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.codomain()
            Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self._codomain

    def defining_morphism(self):
        """
        Return the morphism defining the structure of algebra
        of the codomain over the domain.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: M = R.derivation_module()
            sage: M.defining_morphism()
            Identity endomorphism of Univariate Polynomial Ring in x over Rational Field

            sage: S.<y> = R[]
            sage: M = R.derivation_module(S)
            sage: M.defining_morphism()
            Polynomial base injection morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Univariate Polynomial Ring in y over
                     Univariate Polynomial Ring in x over Rational Field

            sage: ev = R.hom([QQ(0)])
            sage: M = R.derivation_module(ev)
            sage: M.defining_morphism()
            Ring morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Rational Field
              Defn: x |--> 0
        """
        return self._defining_morphism

    def twisting_morphism(self):
        r"""
        Return the twisting homomorphism of the derivations in this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: M = R.derivation_module(twist=theta); M
            Module of twisted derivations over Multivariate Polynomial Ring in x, y
             over Integer Ring (twisting morphism: x |--> y, y |--> x)
            sage: M.twisting_morphism()
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Integer Ring
              Defn: x |--> y
                    y |--> x

        When the derivations are untwisted, this method returns nothing::

            sage: M = R.derivation_module()
            sage: M.twisting_morphism()
        """
        return self._twist

    def ngens(self):
        r"""
        Return the number of generators of this module of derivations.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.ngens()
            2

        Indeed, generators are::

            sage: M.gens()
            (d/dx, d/dy)

        We check that, for a nontrivial twist over a field, the module of
        twisted derivation is a vector space of dimension 1 generated by
        ``twist - id``::

            sage: K = R.fraction_field()
            sage: theta = K.hom([K(y),K(x)])
            sage: M = K.derivation_module(twist=theta); M
            Module of twisted derivations over Fraction Field of Multivariate Polynomial
             Ring in x, y over Integer Ring (twisting morphism: x |--> y, y |--> x)
            sage: M.ngens()
            1
            sage: M.gen()
            [x |--> y, y |--> x] - id
        """
        if self._gens is None:
            raise NotImplementedError("generators are not implemented for this derivation module")
        return len(self._gens)

    def gens(self):
        r"""
        Return the generators of this module of derivations.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.gens()
            (d/dx, d/dy)

        We check that, for a nontrivial twist over a field, the module of
        twisted derivation is a vector space of dimension 1 generated by
        ``twist - id``::

            sage: K = R.fraction_field()
            sage: theta = K.hom([K(y),K(x)])
            sage: M = K.derivation_module(twist=theta); M
            Module of twisted derivations over Fraction Field of Multivariate Polynomial
             Ring in x, y over Integer Ring (twisting morphism: x |--> y, y |--> x)
            sage: M.gens()
            ([x |--> y, y |--> x] - id,)
        """
        if self._gens is None:
            raise NotImplementedError("generators are not implemented for this derivation module")
        return tuple(self._gens)

    def gen(self, n=0):
        r"""
        Return the `n`-th generator of this module of derivations.

        INPUT:

        - ``n`` -- integer (default: `0`)

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.gen()
            d/dx
            sage: M.gen(1)
            d/dy
        """
        if self._gens is None:
            raise NotImplementedError("generators are not implemented for this derivation module")
        try:
            return self._gens[n]
        except IndexError:
            raise ValueError("generator not defined")

    def basis(self):
        r"""
        Return a basis of this module of derivations.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module()
            sage: M.basis()
            Family (d/dx, d/dy)
        """
        if self._basis is None:
            raise NotImplementedError("basis is not implemented for this derivation module")
        return Family(self._basis)

    def dual_basis(self):
        r"""
        Return the dual basis of the canonical basis of this module of
        derivations (which is that returned by the method :meth:`basis`).

        .. NOTE::

            The dual basis of `(d_1, \dots, d_n)` is a family
            `(x_1, \ldots, x_n)` of elements in the domain such
            that `d_i(x_i) = 1` and `d_i(x_j) = 0` if `i \neq j`.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module()
            sage: M.basis()
            Family (d/dx, d/dy)
            sage: M.dual_basis()
            Family (x, y)
        """
        if self._dual_basis is None:
            raise NotImplementedError("basis is not implemented for this derivation module")
        return Family(self._dual_basis)

    def ring_of_constants(self):
        r"""
        Return the subring of the domain consisting of elements
        `x` such that `d(x) = 0` for all derivation `d` in this module.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: M = R.derivation_module()
            sage: M.basis()
            Family (d/dx, d/dy)
            sage: M.ring_of_constants()
            Rational Field
        """
        if not self._constants[1]:
            raise NotImplementedError("the computation of the ring of constants"
                                      " is not implemented for this derivation module")
        return self._constants[0]

    def random_element(self, *args, **kwds):
        r"""
        Return a random derivation in this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module()
            sage: M.random_element()  # random
            (x^2 + x*y - 3*y^2 + x + 1)*d/dx + (-2*x^2 + 3*x*y + 10*y^2 + 2*x + 8)*d/dy
        """
        if self._gens is None:
            raise NotImplementedError("generators are not implemented for this derivation module")
        return self([ self._codomain.random_element(*args, **kwds) for _ in range(len(self._gens)) ])

    def some_elements(self):
        r"""
        Return a list of elements of this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module()
            sage: M.some_elements()
            [d/dx, d/dy, x*d/dx, x*d/dy, y*d/dx, y*d/dy]
        """
        if self._gens is None:
            return self.an_element()
        if self._dual_basis is None:
            return self._gens
        return self._gens + [f * D for f in self._dual_basis for D in self._gens]


# The class RingDerivation does not derive from Map (or RingMap)
# because we don't want to see derivations as morphisms in some
# category since they are not stable by composition.
class RingDerivation(ModuleElement):
    r"""
    An abstract class for twisted and untwisted derivations over
    commutative rings.

    TESTS::

        sage: R.<x,y> = ZZ[]
        sage: f = R.derivation(x) + 2*R.derivation(y); f
        d/dx + 2*d/dy
        sage: f(x*y)
        2*x + y
    """
    def __call__(self, x):
        """
        Return the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x*R.derivation(x) + y*R.derivation(y)
            sage: f(x^2 + 3*x*y - y^2)
            2*x^2 + 6*x*y - 2*y^2
        """
        arg = self.parent().domain()(x)
        return self._call_(arg)

    def domain(self):
        """
        Return the domain of this derivation.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.derivation(y); f
            d/dy
            sage: f.domain()
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: f.domain() is R
            True
        """
        return self.parent().domain()

    def codomain(self):
        """
        Return the codomain of this derivation.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = R.derivation(); f
            d/dx
            sage: f.codomain()
            Univariate Polynomial Ring in x over Rational Field
            sage: f.codomain() is R
            True

        ::

            sage: S.<y> = R[]
            sage: M = R.derivation_module(S)
            sage: M.random_element().codomain()
            Univariate Polynomial Ring in y over
             Univariate Polynomial Ring in x over Rational Field
            sage: M.random_element().codomain() is S
            True
        """
        return self.parent().codomain()


class RingDerivationWithoutTwist(RingDerivation):
    """
    An abstract class for untwisted derivations.
    """
    def _repr_(self):
        r"""
        Return a string representation of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: R.derivation(x)
            d/dx
            sage: R.derivation(y)
            d/dy
        """
        parent = self.parent()
        try:
            dual_basis = parent.dual_basis()
        except NotImplementedError:
            return "A derivation on %s" % parent.domain()
        coeffs = self.list()
        s = ""
        for i in range(len(dual_basis)):
            c = coeffs[i]
            sc = str(c)
            if sc == "0":
                continue
            name = str(dual_basis[i])
            if all(name.find(c) == -1 for c in "+-*/ ()"):
                ddx = "d/d%s" % name
            else:
                ddx = "d_{%s}" % name
            if sc == "1":
                s += " + " + ddx
            elif sc == "-1":
                s += " - " + ddx
            elif c._is_atomic() and sc[0] != "-":
                s += " + %s*%s" % (sc, ddx)
            elif (-c)._is_atomic():
                s += " - %s*%s" % (-c, ddx)
            else:
                s += " + (%s)*%s" % (sc, ddx)
        if s[:3] == " + ":
            return s[3:]
        elif s[:3] == " - ":
            return "-" + s[3:]
        elif s == "":
            return "0"
        else:
            return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: ddx = R.derivation(x)
            sage: ddy = R.derivation(y)

            sage: latex(ddx)
            \frac{d}{dx}
            sage: latex(ddy)
            \frac{d}{dy}
            sage: latex(ddx + ddy)
            \frac{d}{dx} + \frac{d}{dy}
        """
        parent = self.parent()
        try:
            dual_basis = parent.dual_basis()
        except NotImplementedError:
            return "\\text{A derivation on } %s" % latex(parent.domain())
        coeffs = self.list()
        s = ""
        for i in range(len(dual_basis)):
            c = coeffs[i]
            sc = str(c)
            if sc == "0":
                continue
            name = str(dual_basis[i])
            if all(name.find(c) == -1 for c in "+-*/ ()"):
                ddx = "\\frac{d}{d%s}" % latex(dual_basis[i])
            else:
                ddx = "d_{%s}" % latex(dual_basis[i])
            if sc == "1":
                s += " + " + ddx
            elif sc == "-1":
                s += " - " + ddx
            elif c._is_atomic() and sc[0] != "-":
                s += " + %s %s" % (sc, ddx)
            elif (-c)._is_atomic():
                s += " - %s %s" % (-c, ddx)
            else:
                s += " + \\left(%s\\right) %s" % (sc, ddx)
        if s[:3] == " + ":
            return s[3:]
        elif s[:3] == " - ":
            return "-" + s[3:]
        elif s == "":
            return "0"
        else:
            return s

    def list(self):
        """
        Return the list of coefficient of this derivation
        on the canonical basis.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: M = R.derivation_module()
            sage: M.basis()
            Family (d/dx, d/dy)

            sage: R.derivation(x).list()
            [1, 0]
            sage: R.derivation(y).list()
            [0, 1]

            sage: f = x*R.derivation(x) + y*R.derivation(y); f
            x*d/dx + y*d/dy
            sage: f.list()
            [x, y]
        """
        parent = self.parent()
        return [self(x) for x in parent.dual_basis()]

    def monomial_coefficients(self):
        r"""
        Return dictionary of nonzero coordinates (on the canonical
        basis) of this derivation.

        More precisely, this returns a dictionary whose keys are indices
        of basis elements and whose values are the corresponding coefficients.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: M = R.derivation_module()
            sage: M.basis()
            Family (d/dx, d/dy)

            sage: R.derivation(x).monomial_coefficients()
            {0: 1}
            sage: R.derivation(y).monomial_coefficients()
            {1: 1}

            sage: f = x*R.derivation(x) + y*R.derivation(y); f
            x*d/dx + y*d/dy
            sage: f.monomial_coefficients()
            {0: x, 1: y}
        """
        dual_basis = self.parent().dual_basis()
        dict = { }
        for i in range(len(dual_basis)):
            c = self(dual_basis[i])
            if c != 0:
                dict[i] = c
        return dict

    def is_zero(self):
        """
        Return ``True`` if this derivation is zero.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(); f
            d/dx
            sage: f.is_zero()
            False

            sage: (f-f).is_zero()
            True
        """
        for c in self.list():
            if not c.is_zero():
                return False
        return True

    def _richcmp_(self, other, op):
        """
        Compare this derivation with ``other`` according
        to the comparison operator ``op``.

        EXAMPLES::

            sage: R.<x,y,z> = GF(5)[]
            sage: D = sum(v*R.derivation(v) for v in R.gens()); D
            x*d/dx + y*d/dy + z*d/dz
            sage: D.pth_power() == D
            True
        """
        if op == op_EQ:
            if isinstance(other, RingDerivationWithoutTwist):
                return self.list() == other.list()
            else:
                return False
        if op == op_NE:
            if isinstance(other, RingDerivationWithoutTwist):
                return self.list() != other.list()
            else:
                return True
        return NotImplemented

    def _bracket_(self, other):
        """
        Return the Lie bracket (that is the commutator) of
        this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: Dx = R.derivation(x)
            sage: Dy = R.derivation(y)
            sage: Dx._bracket_(Dy)
            0

            sage: Dx.bracket(x*Dy)
            d/dy

        TESTS::

            sage: M = R.derivation_module()
            sage: X = M.random_element()
            sage: X.bracket(X)
            0

        We check the Jacobi identity::

            sage: Y = M.random_element()
            sage: Z = M.random_element()
            sage: X.bracket(Y.bracket(Z)) + Y.bracket(Z.bracket(X)) + Z.bracket(X.bracket(Y))
            0

        and the product rule::

            sage: f = R.random_element()
            sage: X.bracket(f*Y) == X(f)*Y + f*X.bracket(Y)
            True
        """
        parent = self.parent()
        if parent.domain() is not parent.codomain():
            raise TypeError("the bracket is only defined for derivations with same domain and codomain")
        arg = [self(other(x)) - other(self(x)) for x in parent.dual_basis()]
        return parent(arg)

    def pth_power(self):
        r"""
        Return the `p`-th power of this derivation where `p`
        is the characteristic of the domain.

        .. NOTE::

            Leibniz rule implies that this is again a derivation.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = GF(5)[]
            sage: Dx = R.derivation(x)
            sage: Dx.pth_power()
            0
            sage: (x*Dx).pth_power()
            x*d/dx
            sage: (x^6*Dx).pth_power()
            x^26*d/dx

            sage: Dy = R.derivation(y)                                                  # needs sage.rings.finite_rings
            sage: (x*Dx + y*Dy).pth_power()                                             # needs sage.rings.finite_rings
            x*d/dx + y*d/dy

        An error is raised if the domain has characteristic zero::

            sage: R.<x,y> = QQ[]
            sage: Dx = R.derivation(x)
            sage: Dx.pth_power()
            Traceback (most recent call last):
            ...
            TypeError: the domain of the derivation must have positive and prime characteristic

        or if the characteristic is not a prime number::

            sage: R.<x,y> = Integers(10)[]
            sage: Dx = R.derivation(x)
            sage: Dx.pth_power()
            Traceback (most recent call last):
            ...
            TypeError: the domain of the derivation must have positive and prime characteristic

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = GF(3)[]
            sage: D = R.derivation_module().random_element()
            sage: Dp = D.pth_power()
            sage: f = R.random_element()
            sage: Dp(f) == D(D(D(f)))
            True

            sage: D.bracket(Dp)                                                         # needs sage.rings.finite_rings
            0
        """
        parent = self.parent()
        if parent.domain() is not parent.codomain():
            raise TypeError("the derivation must have the same domain and codomain")
        p = parent.domain().characteristic()
        if not p.is_prime():
            raise TypeError("the domain of the derivation must have positive and prime characteristic")
        arg = [ ]
        for x in parent.dual_basis():
            res = x
            for _ in range(p):
                res = self(res)
            arg.append(res)
        return parent(arg)

    def precompose(self, morphism):
        r"""
        Return the derivation obtained by applying first
        ``morphism`` and then this derivation.

        INPUT:

        - ``morphism`` -- a homomorphism of rings whose codomain is
          the domain of this derivation or a ring that coerces to
          the domain of this derivation

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: B.<x,y> = QQ[]
            sage: D = B.derivation(x) - 2*x*B.derivation(y); D
            d/dx - 2*x*d/dy

        When restricting to ``A``, the term ``d/dy`` disappears
        (since it vanishes on ``A``)::

            sage: D.precompose(A)
            d/dx

        If we restrict to another well chosen subring, the derivation vanishes::

            sage: C.<t> = QQ[]
            sage: f = C.hom([x^2 + y]); f
            Ring morphism:
              From: Univariate Polynomial Ring in t over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn: t |--> x^2 + y
            sage: D.precompose(f)
            0

        Note that this method cannot be used to compose derivations::

            sage: D.precompose(D)
            Traceback (most recent call last):
            ...
            TypeError: you must give a homomorphism of rings

        TESTS::

            sage: D.precompose(C)
            Traceback (most recent call last):
            ...
            TypeError: the given ring does not coerce to the domain of the derivation
        """
        parent = self.parent()
        if morphism in Rings().Commutative():
            if parent.domain().has_coerce_map_from(morphism):
                morphism = parent.domain().coerce_map_from(morphism)
            else:
                raise TypeError("the given ring does not coerce to the domain of the derivation")
        elif not (isinstance(morphism, Map) and morphism.category_for().is_subcategory(Rings())):
            raise TypeError("you must give a homomorphism of rings")
        M = RingDerivationModule(morphism.domain(), parent.defining_morphism() * morphism)
        arg = [self(morphism(x)) for x in M.dual_basis()]
        return M(arg)

    def postcompose(self, morphism):
        """
        Return the derivation obtained by applying first
        this derivation and then ``morphism``.

        INPUT:

        - ``morphism`` -- a homomorphism of rings whose domain is
          the codomain of this derivation or a ring into which the
          codomain of this derivation coerces

        EXAMPLES::

            sage: A.<x,y>= QQ[]
            sage: ev = A.hom([QQ(0), QQ(1)])
            sage: Dx = A.derivation(x)
            sage: Dy = A.derivation(y)

        We can define the derivation at `(0,1)` just by postcomposing
        with ``ev``::

            sage: dx = Dx.postcompose(ev)
            sage: dy = Dy.postcompose(ev)
            sage: f = x^2 + y^2
            sage: dx(f)
            0
            sage: dy(f)
            2

        Note that we cannot avoid the creation of the evaluation morphism:
        if we pass in ``QQ`` instead, an error is raised since there is
        no coercion morphism from ``A`` to ``QQ``::

            sage: Dx.postcompose(QQ)
            Traceback (most recent call last):
            ...
            TypeError: the codomain of the derivation does not coerce to the given ring

        Note that this method cannot be used to compose derivations::

            sage: Dx.precompose(Dy)
            Traceback (most recent call last):
            ...
            TypeError: you must give a homomorphism of rings
        """
        parent = self.parent()
        if morphism in Rings().Commutative():
            if morphism.has_coerce_map_from(parent.codomain()):
                morphism = morphism.coerce_map_from(parent.codomain())
            else:
                raise TypeError("the codomain of the derivation does not coerce to the given ring")
        elif not (isinstance(morphism, Map) and morphism.category_for().is_subcategory(Rings())):
            raise TypeError("you must give a homomorphism of rings")
        M = RingDerivationModule(parent.domain(), morphism * parent.defining_morphism())
        arg = [morphism(self(x)) for x in M.dual_basis()]
        return M(arg)

    def extend_to_fraction_field(self):
        r"""
        Return the extension of this derivation to fraction fields of
        the domain and the codomain.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: d = S.derivation()
            sage: d
            d/dx

            sage: D = d.extend_to_fraction_field()
            sage: D
            d/dx
            sage: D.domain()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field

            sage: D(1/x)
            -1/x^2
        """
        parent = self.parent()
        domain = parent.domain().fraction_field()
        codomain = parent.codomain().fraction_field()
        M = RingDerivationModule(domain, codomain)
        try:
            return M(self)
        except (ValueError, NotImplementedError):
            return M(self.list())


class RingDerivationWithoutTwist_zero(RingDerivationWithoutTwist):
    """
    This class can only represent the zero derivation.

    It is used when the parent is the zero derivation module
    (e.g., when its domain is ``ZZ``, ``QQ``, a finite field, etc.)
    """
    def __init__(self, parent, arg=None):
        """
        Initialize this derivation.

        TESTS::

            sage: M = ZZ.derivation_module()
            sage: der = M(); der
            0

            sage: from sage.rings.derivation import RingDerivationWithoutTwist_zero
            sage: isinstance(der, RingDerivationWithoutTwist_zero)
            True

            sage: TestSuite(der).run()
        """
        if isinstance(arg, list) and len(arg) == 1 and isinstance(arg[0], RingDerivation):
            arg = arg[0]
        if arg and not (isinstance(arg, RingDerivation) and arg.is_zero()):
            raise ValueError("unable to create the derivation")
        RingDerivation.__init__(self, parent)

    def _repr_(self):
        """
        Return a string representation of this derivation.

        EXAMPLES::

            sage: M = ZZ.derivation_module()
            sage: M()
            0
        """
        return "0"

    def _latex_(self):
        """
        Return a string representation of this derivation.

        EXAMPLES::

            sage: M = ZZ.derivation_module()
            sage: latex(M())
            0
        """
        return "0"

    def __hash__(self):
        """
        Return a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(x)
            sage: hash(f)  # random
            3713081631936575706
        """
        return hash(tuple(self.list()))

    def _add_(self, other):
        """
        Return the sum of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: Dy = R.derivation(y)
            sage: Dx + Dy
            d/dx + d/dy
        """
        return other

    def _sub_(self, other):
        """
        Return the difference of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: Dy = R.derivation(y)
            sage: Dx - Dy
            d/dx - d/dy
        """
        return -other

    def _neg_(self):
        """
        Return the opposite of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: -Dx
            -d/dx
        """
        return self

    def _lmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: Dx * 2
            2*d/dx
            sage: Dx * x^2
            x^2*d/dx
        """
        return self

    def _rmul_(self, left):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: 2 * Dx
            2*d/dx
            sage: x^2 * Dx
            x^2*d/dx
        """
        return self

    def _call_(self, x):
        """
        Return the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x*R.derivation(x) + y*R.derivation(y)
            sage: f(x^2 + 3*x*y - y^2)
            2*x^2 + 6*x*y - 2*y^2
        """
        return self.parent().codomain().zero()

    def _bracket_(self, other):
        """
        Return the Lie bracket (that is the commutator) of
        this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: Dx = R.derivation(x)
            sage: Dy = R.derivation(y)
            sage: Dx._bracket_(Dy)
            0
        """
        return self

    def is_zero(self):
        """
        Return ``True`` if this derivation vanishes.

        EXAMPLES::

            sage: M = QQ.derivation_module()
            sage: M().is_zero()
            True
        """
        return True

    def list(self):
        """
        Return the list of coefficient of this derivation
        on the canonical basis.

        EXAMPLES::

            sage: M = QQ.derivation_module()
            sage: M().list()
            []
        """
        return []


class RingDerivationWithoutTwist_wrapper(RingDerivationWithoutTwist):
    """
    This class is a wrapper for derivation.

    It is useful for changing the parent without changing the
    computation rules for derivations. It is used for derivations
    over fraction fields and quotient rings.
    """
    def __init__(self, parent, arg=None):
        """
        Initialize this derivation.

        TESTS::

            sage: # needs sage.libs.singular
            sage: from sage.rings.derivation import RingDerivationWithoutTwist_wrapper
            sage: R.<x,y> = GF(5)[]
            sage: S = R.quo([x^5, y^5])
            sage: M = S.derivation_module()
            sage: der = M.random_element()
            sage: isinstance(der, RingDerivationWithoutTwist_wrapper)
            True
            sage: TestSuite(der).run()
        """
        if isinstance(arg, list) and len(arg) == 1 and isinstance(arg[0], RingDerivation):
            arg = arg[0]
        if isinstance(arg, RingDerivationWithoutTwist_wrapper):
            self._base_derivation = arg._base_derivation
        else:
            self._base_derivation = parent._base_derivation(arg)
        RingDerivation.__init__(self, parent)

    def __hash__(self):
        """
        Return a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(x)
            sage: hash(f)  # random
            3713081631936575706
        """
        return hash(tuple(self.list()))

    def _add_(self, other):
        """
        Return the sum of this derivation and ``other``.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<X,Y> = GF(5)[]
            sage: S.<x,y> = R.quo([X^5, Y^5])
            sage: Dx = S.derivation(x)
            sage: Dy = S.derivation(y)
            sage: Dx + Dy
            d/dx + d/dy
        """
        return type(self)(self.parent(), self._base_derivation + other._base_derivation)

    def _sub_(self, other):
        """
        Return the difference of this derivation and ``other``.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<X,Y> = GF(5)[]
            sage: S.<x,y> = R.quo([X^5, Y^5])
            sage: Dx = S.derivation(x)
            sage: Dy = S.derivation(y)
            sage: Dx - Dy
            d/dx - d/dy
        """
        return type(self)(self.parent(), self._base_derivation - other._base_derivation)

    def _neg_(self):
        """
        Return the opposite of this derivation.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<X,Y> = GF(5)[]
            sage: S.<x,y> = R.quo([X^5, Y^5])
            sage: Dx = S.derivation(x)
            sage: -Dx
            -d/dx
        """
        return type(self)(self.parent(), -self._base_derivation)

    def _lmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<X,Y> = GF(5)[]
            sage: S.<x,y> = R.quo([X^5, Y^5])
            sage: Dx = S.derivation(x)
            sage: Dx * 2
            2*d/dx
            sage: Dx * x^2
            x^2*d/dx
        """
        return type(self)(self.parent(), self._base_derivation * factor)

    def _rmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<X,Y> = GF(5)[]
            sage: S.<x,y> = R.quo([X^5, Y^5])
            sage: Dx = S.derivation(x)
            sage: 2 * Dx
            2*d/dx
            sage: x^2 * Dx
            x^2*d/dx
        """
        return type(self)(self.parent(), factor * self._base_derivation)

    def list(self):
        """
        Return the list of coefficient of this derivation
        on the canonical basis.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<X,Y> = GF(5)[]
            sage: S.<x,y> = R.quo([X^5, Y^5])
            sage: M = S.derivation_module()
            sage: M.basis()
            Family (d/dx, d/dy)
            sage: S.derivation(x).list()
            [1, 0]
            sage: S.derivation(y).list()
            [0, 1]
            sage: f = x*S.derivation(x) + y*S.derivation(y); f
            x*d/dx + y*d/dy
            sage: f.list()
            [x, y]
        """
        return self._base_derivation.list()


class RingDerivationWithoutTwist_function(RingDerivationWithoutTwist):
    """
    A class for untwisted derivations over rings whose elements
    are either polynomials, rational fractions, power series or
    Laurent series.
    """
    def __init__(self, parent, arg=None):
        """
        Initialize this derivation.

        TESTS::

            sage: R.<x,y> = ZZ[]
            sage: R.derivation(x)
            d/dx
            sage: der = R.derivation([1,2])
            sage: der
            d/dx + 2*d/dy

            sage: TestSuite(der).run()
        """
        domain = parent.domain()
        codomain = parent.codomain()
        ngens = domain.ngens()
        self._base_derivation = parent._base_derivation()
        self._images = [codomain.zero() for _ in range(ngens)]
        if arg is None:
            arg = domain.gen()
        if isinstance(arg, list) and len(arg) == 1 and isinstance(arg[0], RingDerivation):
            arg = arg[0]
        if not arg:
            pass
        elif (isinstance(arg, RingDerivationWithoutTwist_function)
              and parent.has_coerce_map_from(arg.parent())):
            self._base_derivation = parent._base_derivation(arg._base_derivation)
            self._images = [codomain(x) for x in arg._images]
        elif isinstance(arg, (tuple, list)):
            if len(arg) < ngens:
                raise ValueError("the number of images is incorrect")
            self._base_derivation = parent._base_derivation(arg[:-ngens])
            self._images = [codomain(x) for x in arg[-ngens:]]
        else:
            for i in range(ngens):
                if arg == domain.gen(i):
                    self._base_derivation = parent._base_derivation()
                    self._images[i] = codomain.one()
                    break
            else:
                self._base_derivation = parent._base_derivation(arg)
        RingDerivation.__init__(self, parent)

    def __hash__(self):
        """
        Return a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(x)
            sage: hash(f)  # random
            3713081631936575706
        """
        return hash(tuple(self.list()))

    def _add_(self, other):
        """
        Return the sum of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: Dy = R.derivation(y)
            sage: Dx + Dy
            d/dx + d/dy
        """
        base_derivation = self._base_derivation + other._base_derivation
        im = [ self._images[i] + other._images[i] for i in range(self.parent().domain().ngens()) ]
        return type(self)(self.parent(), [base_derivation] + im)

    def _sub_(self, other):
        """
        Return the subtraction of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: Dy = R.derivation(y)
            sage: Dx - Dy
            d/dx - d/dy
        """
        base_derivation = self._base_derivation - other._base_derivation
        im = [ self._images[i] - other._images[i] for i in range(self.parent().domain().ngens()) ]
        return type(self)(self.parent(), [base_derivation] + im)

    def _rmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: 2 * Dx
            2*d/dx
            sage: x^2 * Dx
            x^2*d/dx
        """
        factor = self.parent().codomain()(factor)
        base_derivation = factor * self._base_derivation
        im = [ factor*x for x in self._images ]
        return type(self)(self.parent(), [base_derivation] + im)

    def _lmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Dx = R.derivation(x)
            sage: Dx * 2
            2*d/dx
            sage: Dx * x^2
            x^2*d/dx
        """
        return self._rmul_(factor)

    def _call_(self, x):
        """
        Return the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: D = x*R.derivation(x) + y*R.derivation(y)
            sage: D(x^2 + 3*x*y - y^2)
            2*x^2 + 6*x*y - 2*y^2
        """
        parent = self.parent()
        domain = parent.domain()
        codomain = parent.codomain()
        defining_morphism = parent.defining_morphism()
        if isinstance(domain, FractionField_generic):
            num = x.numerator()
            den = x.denominator()
            u = defining_morphism(num)
            v = defining_morphism(den)
            up = num.map_coefficients(self._base_derivation, codomain)(*domain.gens())
            vp = den.map_coefficients(self._base_derivation, codomain)(*domain.gens())
            res = (up*v - u*vp) / (v*v)
        else:
            res = x.map_coefficients(self._base_derivation, codomain)(*domain.gens())
        for i in range(len(self._images)):
            res += defining_morphism(x.derivative(domain.gen(i))) * self._images[i]
        return res

    def is_zero(self):
        """
        Return ``True`` if this derivation is zero.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(); f
            d/dx
            sage: f.is_zero()
            False

            sage: (f-f).is_zero()
            True
        """
        if not self._base_derivation.is_zero():
            return False
        return all(im == 0 for im in self._images)

    def list(self):
        """
        Return the list of coefficient of this derivation
        on the canonical basis.

        EXAMPLES::

            sage: R.<x,y> = GF(5)[[]]
            sage: M = R.derivation_module()
            sage: M.basis()
            Family (d/dx, d/dy)

            sage: R.derivation(x).list()
            [1, 0]
            sage: R.derivation(y).list()
            [0, 1]

            sage: f = x*R.derivation(x) + y*R.derivation(y); f
            x*d/dx + y*d/dy
            sage: f.list()
            [x, y]
        """
        return self._base_derivation.list() + self._images


class RingDerivationWithoutTwist_fraction_field(RingDerivationWithoutTwist_wrapper):
    """
    This class handles derivations over fraction fields.
    """
    def __hash__(self):
        """
        Return a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(x)
            sage: hash(f)  # random
            3713081631936575706
        """
        return hash(tuple(self.list()))

    def _call_(self, x):
        """
        Return the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: K = R.fraction_field()
            sage: f = K.derivation(); f
            d/dx
            sage: f(1/x)
            (-1)/x^2
        """
        defining_morphism = self.parent().defining_morphism()
        num = x.numerator()
        den = x.denominator()
        u = defining_morphism(num)
        v = defining_morphism(den)
        up = self._base_derivation(u)
        vp = self._base_derivation(v)
        return (up*v - u*vp) / (v*v)


class RingDerivationWithoutTwist_quotient(RingDerivationWithoutTwist_wrapper):
    """
    This class handles derivations over quotient rings.
    """
    def __hash__(self):
        """
        Return a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(x)
            sage: hash(f)  # random
            3713081631936575706
        """
        return hash(tuple(self.list()))

    def _call_(self, x):
        """
        Return the image of ``x`` under this derivation.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<X,Y> = GF(5)[]
            sage: S.<x,y> = R.quo([X^5, Y^5])
            sage: f = x^3*S.derivation(); f
            x^3*d/dx
            sage: f(x^3)
            0
        """
        return self._base_derivation(x.lift())


class RingDerivationWithTwist_generic(RingDerivation):
    r"""
    The class handles `\theta`-derivations of the form
    `\lambda (\theta - \iota)` (where `\iota` is the defining
    morphism of the codomain over the domain) for a scalar
    `\lambda` varying in the codomain.
    """
    def __init__(self, parent, scalar=0):
        """
        Initialize this derivation.

        TESTS::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: R.derivation(twist=theta)
            0
            sage: R.derivation(1, twist=theta)
            [x |--> y, y |--> x] - id

            sage: der = R.derivation(x, twist=theta)
            sage: TestSuite(der).run()
        """
        codomain = parent.codomain()
        self._scalar = codomain(scalar)
        RingDerivation.__init__(self, parent)

    def __hash__(self):
        """
        Return a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: f = R.derivation(1, twist=theta)
            sage: hash(f)  # random
            -6511057926760520014
        """
        return hash(self._scalar)

    def _repr_(self):
        r"""
        Return a string representation of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: R.derivation(1, twist=theta)
            [x |--> y, y |--> x] - id
        """
        scalar = self._scalar
        sc = str(scalar)
        if sc == "0":
            return "0"
        defining_morphism = self.parent().defining_morphism()
        twisting_morphism = self.parent().twisting_morphism()
        try:
            if defining_morphism.is_identity():
                sdef = "id"
            else:
                sdef = "[%s]" % defining_morphism._repr_short()
        except (AttributeError, NotImplementedError):
            sdef = "defining_morphism"
        try:
            stwi = "[%s]" % twisting_morphism._repr_short()
        except AttributeError:
            stwi = "twisting_morphism"
        if sc == "1":
            return "%s - %s" % (stwi, sdef)
        elif sc == "-1":
            s = "-"
        elif scalar._is_atomic():
            s = "%s*" % sc
        elif (-scalar)._is_atomic():
            s = "-%s*" % (-scalar)
        else:
            s = "(%s)*" % sc
        return "%s(%s - %s)" % (s, stwi, sdef)

    def _latex_(self):
        r"""
        Return a LaTeX representation of this derivation.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a + 1, twist=Frob)
            sage: latex(der)
            \left(a + 1\right)  \left(\left[a \mapsto a^{5}\right] - \text{id}\right)
        """
        scalar = self._scalar
        sc = str(scalar)
        if sc == "0":
            return "0"
        defining_morphism = self.parent().defining_morphism()
        twisting_morphism = self.parent().twisting_morphism()
        try:
            if defining_morphism.is_identity():
                sdef = "\\text{id}"
            else:
                sdef = "\\left[%s\\right]" % latex(defining_morphism)
        except (AttributeError, NotImplementedError):
            sdef = "\\text{defining morphism}"
        try:
            stwi = "\\left[%s\\right]" % latex(twisting_morphism)
        except AttributeError:
            stwi = "\\text{twisting morphism}"
        if sc == "1":
            return "%s - %s" % (stwi, sdef)
        elif sc == "-1":
            s = "-"
        elif scalar._is_atomic():
            s = "%s " % sc
        elif (-scalar)._is_atomic():
            s = "-%s " % (-scalar)
        else:
            s = "\\left(%s\\right) " % sc
        return "%s \\left(%s - %s\\right)" % (s, stwi, sdef)

    def _add_(self, other):
        """
        Return the sum of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: der1 = R.derivation(x, twist=theta); der1
            x*([x |--> y, y |--> x] - id)
            sage: der2 = R.derivation(y, twist=theta); der2
            y*([x |--> y, y |--> x] - id)
            sage: der1 + der2
            (x + y)*([x |--> y, y |--> x] - id)
        """
        return type(self)(self.parent(), self._scalar + other._scalar)

    def _sub_(self, other):
        """
        Return the subtraction of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: der1 = R.derivation(x, twist=theta); der1
            x*([x |--> y, y |--> x] - id)
            sage: der2 = R.derivation(y, twist=theta); der2
            y*([x |--> y, y |--> x] - id)
            sage: der1 - der2
            (x - y)*([x |--> y, y |--> x] - id)

        TESTS::

            sage: der1 - der1
            0
            sage: der2 - der2
            0
        """
        return type(self)(self.parent(), self._scalar - other._scalar)

    def _rmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: D = R.derivation(x, twist=theta); D
            x*([x |--> y, y |--> x] - id)
            sage: y * D
            x*y*([x |--> y, y |--> x] - id)
        """
        return type(self)(self.parent(), factor * self._scalar)

    def _lmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: D = R.derivation(x, twist=theta); D
            x*([x |--> y, y |--> x] - id)
            sage: D * y
            x*y*([x |--> y, y |--> x] - id)
        """
        return self._rmul_(factor)

    def _call_(self, x):
        """
        Return the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: f = R.derivation(1, twist=theta); f
            [x |--> y, y |--> x] - id
            sage: f(x)
            -x + y
        """
        parent = self.parent()
        return self._scalar * (parent.twisting_morphism()(x) - parent.defining_morphism()(x))

    def list(self):
        """
        Return the list of coefficient of this twisted derivation
        on the canonical basis.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: K = R.fraction_field()
            sage: theta = K.hom([y,x])
            sage: M = K.derivation_module(twist=theta)
            sage: M.basis()
            Family (twisting_morphism - id,)
            sage: f = (x+y) * M.gen()
            sage: f
            (x + y)*(twisting_morphism - id)
            sage: f.list()
            [x + y]
        """
        return [ self._scalar ]

    def precompose(self, morphism):
        r"""
        Return the twisted derivation obtained by applying first
        ``morphism`` and then this twisted derivation.

        INPUT:

        - ``morphism`` -- a homomorphism of rings whose codomain is
          the domain of this derivation or a ring that coerces to
          the domain of this derivation

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: D = R.derivation(x, twist=theta); D
            x*([x |--> y, y |--> x] - id)

            sage: f = R.hom([x^2, y^3])
            sage: g = D.postcompose(f); g
            x^2*([x |--> y^3, y |--> x^2] - [x |--> x^2, y |--> y^3])

        Observe that the `g` is no longer a `\theta`-derivation but
        a `(f \circ \theta)`-derivation::

            sage: g.parent().twisting_morphism()
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Integer Ring
              Defn: x |--> y^3
                    y |--> x^2
        """
        parent = self.parent()
        if morphism in Rings().Commutative():
            if parent.domain().has_coerce_map_from(morphism):
                morphism = parent.domain().coerce_map_from(morphism)
            else:
                raise TypeError("the given ring does not coerce to the domain of the derivation")
        elif not (isinstance(morphism, Map) and morphism.category_for().is_subcategory(Rings())):
            raise TypeError("you must give a homomorphism of rings")
        M = RingDerivationModule(morphism.domain(), parent.defining_morphism() * morphism,
                                 parent.twisting_morphism() * morphism)
        return M(self._scalar)

    def postcompose(self, morphism):
        r"""
        Return the twisted derivation obtained by applying first
        this twisted derivation and then ``morphism``.

        INPUT:

        - ``morphism`` -- a homomorphism of rings whose domain is
          the codomain of this derivation or a ring into which the
          codomain of this derivation

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: D = R.derivation(x, twist=theta); D
            x*([x |--> y, y |--> x] - id)

            sage: f = R.hom([x^2, y^3])
            sage: g = D.precompose(f); g
            x*([x |--> y^2, y |--> x^3] - [x |--> x^2, y |--> y^3])

        Observe that the `g` is no longer a `\theta`-derivation but
        a `(\theta \circ f)`-derivation::

            sage: g.parent().twisting_morphism()
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Integer Ring
              Defn: x |--> y^2
                    y |--> x^3
        """
        parent = self.parent()
        if morphism in Rings().Commutative():
            if morphism.has_coerce_map_from(parent.codomain()):
                morphism = morphism.coerce_map_from(parent.codomain())
            else:
                raise TypeError("the codomain of the derivation does not coerce to the given ring")
        elif not (isinstance(morphism, Map) and morphism.category_for().is_subcategory(Rings())):
            raise TypeError("you must give a homomorphism of rings")
        M = RingDerivationModule(parent.domain(), morphism * parent.defining_morphism(),
                                 morphism * parent.twisting_morphism())
        return M(morphism(self._scalar))

    def _richcmp_(self, other, op):
        """
        Compare this derivation with ``other`` according
        to the comparison operator ``op``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: Dx = R.derivation(x, twist=theta); Dx
            x*([x |--> y, y |--> x] - id)
            sage: Dy = R.derivation(y, twist=theta); Dy
            y*([x |--> y, y |--> x] - id)
            sage: D = R.derivation(x+y, twist=theta); D
            (x + y)*([x |--> y, y |--> x] - id)

            sage: Dx == Dy
            False
            sage: D == Dx + Dy
            True

            sage: D != Dy
            True
        """
        if op == op_EQ:
            if isinstance(other, RingDerivationWithTwist_generic):
                return self._scalar == other._scalar
            else:
                return False
        if op == op_NE:
            if isinstance(other, RingDerivationWithTwist_generic):
                return self._scalar != other._scalar
            else:
                return True
        return NotImplemented

    def extend_to_fraction_field(self):
        r"""
        Return the extension of this derivation to fraction fields of
        the domain and the codomain.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: d = R.derivation(x, twist=theta)
            sage: d
            x*([x |--> y, y |--> x] - id)

            sage: D = d.extend_to_fraction_field(); D                                   # needs sage.libs.singular
            x*([x |--> y, y |--> x] - id)
            sage: D.domain()                                                            # needs sage.libs.singular
            Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring

            sage: D(1/x)                                                                # needs sage.libs.singular
            (x - y)/y
        """
        parent = self.parent()
        domain = parent.domain().fraction_field()
        codomain = parent.codomain().fraction_field()
        twist = parent.twisting_morphism().extend_to_fraction_field()
        M = RingDerivationModule(domain, codomain, twist)
        return M(codomain(self._scalar))
