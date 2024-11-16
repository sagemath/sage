# sage.doctest: needs sage.rings.finite_rings
r"""
Drinfeld module morphisms

This module provides the class
:class:`sage.rings.function_fields.drinfeld_module.morphism.DrinfeldModuleMorphism`.

AUTHORS:
- Antoine Leudière (2022-04)
"""

# *****************************************************************************
#        Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.latex import latex
from sage.categories.morphism import Morphism
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import coerce_binop

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix


class DrinfeldModuleMorphism(Morphism, UniqueRepresentation,
                             metaclass=InheritComparisonClasscallMetaclass):
    r"""
    This class represents Drinfeld `\mathbb{F}_q[T]`-module morphisms.

    Let `\phi` and `\psi` be two Drinfeld `\mathbb{F}_q[T]`-modules over
    a field `K`. A *morphism of Drinfeld modules* `\phi \to \psi` is an
    Ore polynomial `f \in K\{\tau\}` such that `f \phi_a = \psi_a f` for
    every `a \in \mathbb{F}_q[T]`. In our case, this is equivalent to `f
    \phi_T = \psi_T f`. An *isogeny* is a nonzero morphism.

    To create a morphism object, the user should never explicitly
    instantiate :class:`DrinfeldModuleMorphism`, but rather call the
    parent homset with the defining Ore polynomial::

        sage: Fq = GF(4)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(3)
        sage: phi = DrinfeldModule(A, [z, z^2 + z, z^2 + z])
        sage: t = phi.ore_polring().gen()
        sage: ore_pol = t + z^5 + z^3 + z + 1
        sage: psi = phi.velu(ore_pol)
        sage: morphism = Hom(phi, psi)(ore_pol)
        sage: morphism
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> (z^2 + z)*t^2 + (z^2 + z)*t + z
          To:   Drinfeld module defined by T |--> (z^5 + z^2 + z + 1)*t^2 + (z^4 + z + 1)*t + z
          Defn: t + z^5 + z^3 + z + 1


    The given Ore polynomial must indeed define a morphism::

        sage: morphism = Hom(phi, psi)(1)
        Traceback (most recent call last):
        ...
        ValueError: Ore polynomial does not define a morphism

    One can get basic data on the morphism::

        sage: morphism.domain()
        Drinfeld module defined by T |--> (z^2 + z)*t^2 + (z^2 + z)*t + z
        sage: morphism.domain() is phi
        True

        sage: morphism.codomain()
        Drinfeld module defined by T |--> (z^5 + z^2 + z + 1)*t^2 + (z^4 + z + 1)*t + z
        sage: morphism.codomain() is psi
        True

    ::

        sage: morphism.ore_polynomial()
        t + z^5 + z^3 + z + 1
        sage: morphism.ore_polynomial() is ore_pol
        True

    One can check various properties::

        sage: morphism.is_zero()
        False
        sage: morphism.is_isogeny()
        True
        sage: morphism.is_endomorphism()
        False
        sage: morphism.is_isomorphism()
        False

    TESTS::

        sage: morphism.parent() is Hom(phi, psi)
        True
        sage: Hom(phi, psi)(morphism) == morphism
        True

    ::

        sage: phi.frobenius_endomorphism().parent() == End(phi)
        True
        sage: End(phi)(0).parent() == End(phi)
        True

    For the sake of completeness, we explain how the user can directly
    instantiate the class, even though this should never be explicitly
    done::

        sage: from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
        sage: DrinfeldModuleMorphism(Hom(phi, psi), ore_pol)
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> (z^2 + z)*t^2 + (z^2 + z)*t + z
          To:   Drinfeld module defined by T |--> (z^5 + z^2 + z + 1)*t^2 + (z^4 + z + 1)*t + z
          Defn: t + z^5 + z^3 + z + 1
        sage: DrinfeldModuleMorphism(Hom(phi, psi), ore_pol) is morphism
        True
    """
    @staticmethod
    def __classcall_private__(cls, parent, x):
        """
        Create the morphism.

        INPUT:

        - ``cls`` -- the class ``DrinfeldModuleMorphism``

        - ``parent`` -- the Drinfeld module homset

        - ``x`` -- a morphism of Drinfeld modules or an element
          (either an Ore polynomial or an element in the base
          ring) defining it

        TESTS::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: End(phi)(T + 1)
            Endomorphism of Drinfeld module defined by T |--> t^2 + t + z6
              Defn: t^2 + t + z6 + 1

        ::

            sage: t = phi.ore_polring().gen()
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism is Hom(phi, psi)(morphism)
            True

        ::

            sage: from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
            sage: morphism = DrinfeldModuleMorphism(Sets(), t + 1)
            Traceback (most recent call last):
            ...
            TypeError: parent should be a DrinfeldModuleHomset
        """
        from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        if not isinstance(parent, DrinfeldModuleHomset):
            raise TypeError('parent should be a DrinfeldModuleHomset')
        domain = parent.domain()
        codomain = parent.codomain()
        A = domain.category().function_ring()
        # NOTE: it used to be x.parent() is parent, but this was False.
        # DrinfeldModuleHomset inherits Homset, which does NOT inherit
        # UniqueRepresentation
        if x.parent() == parent:
            # x is a DrinfeldModuleMorphism
            ore_pol = x.ore_polynomial()
        elif domain is codomain and A.has_coerce_map_from(x.parent()):
            # x is in the function field; we return the endomorphism phi_x
            ore_pol = domain(A(x))
        else:
            # x is an Ore polynomial
            ore_pol = domain.ore_polring()(x)
        if ore_pol * domain.gen() != codomain.gen() * ore_pol:
            raise ValueError('Ore polynomial does not define a morphism')
        return cls.__classcall__(cls, parent, ore_pol)

    def __init__(self, parent, ore_pol):
        r"""
        Initialize ``self``.

        INPUT:

        - ``parent`` -- the Drinfeld module homset

        - ``ore_pol`` -- the Ore polynomial that defines the morphism

        TESTS::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism._domain is phi
            True
            sage: morphism._codomain is psi
            True
            sage: morphism._ore_polynomial == t + z6^5 + z6^2 + 1
            True
        """
        super().__init__(parent)
        self._domain = parent.domain()
        self._codomain = parent.codomain()
        self._ore_polynomial = ore_pol

    def _latex_(self):
        r"""
        Return a LaTeX representation of the morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: latex(morphism)
            t + z_{6}^{5} + z_{6}^{2} + 1
        """
        return f'{latex(self._ore_polynomial)}'

    def _repr_(self):
        r"""
        Return a string representation of the morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> t^2 + t + z6
              To:   Drinfeld module defined by T |--> t^2 + (z6^4 + z6^2 + 1)*t + z6
              Defn: t + z6^5 + z6^2 + 1
        """
        if self.is_identity():
            return f'Identity morphism of {self._domain}'
        elif self.is_endomorphism():
            return f'Endomorphism of {self._domain}\n' \
                   f'  Defn: {self._ore_polynomial}'
        else:
            return f'Drinfeld Module morphism:\n' \
                   f'  From: {self._domain}\n'  \
                   f'  To:   {self._codomain}\n' \
                   f'  Defn: {self._ore_polynomial}'

    def __hash__(self):
        r"""
        Return a hash of ``self``.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: hash(morphism)  # random
            -4214883752078138009
        """
        return hash((self._domain, self._codomain, self._ore_polynomial))

    def is_zero(self):
        r"""
        Return ``True`` whether the morphism is the zero morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_zero()
            False

        ::

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_zero()
            True
        """
        return self._ore_polynomial.is_zero()

    def is_identity(self):
        r"""
        Return ``True`` whether the morphism is the identity morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: morphism = End(phi)(1)
            sage: morphism.is_identity()
            True

        ::

            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_identity()
            False
        """
        return self._ore_polynomial == 1

    def is_isogeny(self):
        r"""
        Return ``True`` whether the morphism is an isogeny.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_isogeny()
            True

        ::

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_isogeny()
            False

        ::

            sage: identity_morphism = End(phi)(1)
            sage: identity_morphism.is_isogeny()
            True

        ::

            sage: frobenius_endomorphism = phi.frobenius_endomorphism()
            sage: frobenius_endomorphism.is_isogeny()
            True
        """
        return not self.is_zero()

    def is_isomorphism(self):
        r"""
        Return ``True`` whether the morphism is an isomorphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: morphism.is_isomorphism()
            False

        ::

            sage: zero_morphism = End(phi)(0)
            sage: zero_morphism.is_isomorphism()
            False

        ::

            sage: identity_morphism = End(phi)(1)
            sage: identity_morphism.is_isomorphism()
            True

        ::

            sage: frobenius_endomorphism = phi.frobenius_endomorphism()
            sage: frobenius_endomorphism.is_isomorphism()
            False
        """
        return self.is_isogeny() and self._ore_polynomial.degree() == 0

    def ore_polynomial(self):
        r"""
        Return the Ore polynomial that defines the morphism.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z6, 1, 1])
            sage: psi = DrinfeldModule(A, [z6, z6^4 + z6^2 + 1, 1])
            sage: t = phi.ore_polring().gen()
            sage: morphism = Hom(phi, psi)(t + z6^5 + z6^2 + 1)
            sage: ore_pol = morphism.ore_polynomial()
            sage: ore_pol
            t + z6^5 + z6^2 + 1

        ::

            sage: ore_pol * phi(T) == psi(T) * ore_pol
            True
        """
        return self._ore_polynomial

    @coerce_binop
    def __add__(self, other):
        r"""
        Return the sum of this morphism and ``other``.

        INPUT:

        - ``other`` -- a morphism of Drinfeld modules in the
          same category

        TESTS::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: t = phi.ore_variable()
            sage: f = phi.hom(t + 1)
            sage: g = T * f
            sage: f + g  # indirect doctest
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> z*t^3 + t^2 + z
              To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^3 + (3*z^2 + 2*z + 2)*t^2 + (2*z^2 + 3*z + 4)*t + z
              Defn: (2*z^2 + 4*z + 4)*t^4 + (z + 1)*t^3 + t^2 + (2*z^2 + 4*z)*t + z + 1

        We check that coercion from the function ring works::

            sage: F = phi.frobenius_endomorphism()
            sage: F + T
            Endomorphism of Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: (z + 1)*t^3 + t^2 + z
        """
        return self.parent()(self.ore_polynomial() + other.ore_polynomial())

    def _composition_(self, other, H):
        r"""
        Return the composite of this morphism and ``other``.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 1, z, z^2])
            sage: f = phi.frobenius_endomorphism()
            sage: f
            Endomorphism of Drinfeld module defined by T |--> z^2*t^3 + z*t^2 + t + z
              Defn: t^3
            sage: f * f  # indirect doctest
            Endomorphism of Drinfeld module defined by T |--> z^2*t^3 + z*t^2 + t + z
              Defn: t^6
        """
        return H(self.ore_polynomial() * other.ore_polynomial())

    def inverse(self):
        r"""
        Return the inverse of this morphism.

        Only morphisms defined by constant nonzero Ore
        polynomials are invertible.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 1, z, z^2])
            sage: f = phi.hom(2); f
            Endomorphism of Drinfeld module defined by T |--> z^2*t^3 + z*t^2 + t + z
              Defn: 2
            sage: f.inverse()
            Endomorphism of Drinfeld module defined by T |--> z^2*t^3 + z*t^2 + t + z
              Defn: 3

        Inversion of general isomorphisms between different Drinfeld modules
        also works::

            sage: g = phi.hom(z); g
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> z^2*t^3 + z*t^2 + t + z
              To:   Drinfeld module defined by T |--> z^2*t^3 + (z^2 + 2*z + 3)*t^2 + (z^2 + 3*z)*t + z
              Defn: z
            sage: g.inverse()
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> z^2*t^3 + (z^2 + 2*z + 3)*t^2 + (z^2 + 3*z)*t + z
              To:   Drinfeld module defined by T |--> z^2*t^3 + z*t^2 + t + z
              Defn: 3*z^2 + 4

        When the morphism is not invertible, an error is raised::

            sage: F = phi.frobenius_endomorphism()
            sage: F.inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: this morphism is not invertible
        """
        return self.__invert__()

    def __invert__(self):
        r"""
        Return the inverse of this morphism.

        TESTS::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: coeffs = [z] + [K.random_element() for _ in range(10)]
            sage: phi = DrinfeldModule(A, coeffs)
            sage: a = K.random_element()
            sage: while a.is_zero():
            ....:     a = K.random_element()
            sage: f = phi.hom(a)
            sage: g = ~f
            sage: (f*g).is_identity()
            True
            sage: (g*f).is_identity()
            True
        """
        if not self.is_isomorphism():
            raise ZeroDivisionError("this morphism is not invertible")
        H = self.codomain().Hom(self.domain())
        return H(~(self.ore_polynomial()[0]))

    def _motive_matrix(self):
        r"""
        Return the matrix giving the action of this morphism
        on the motives of the underlying Drinfeld modules.

        For internal use. Do not call this method directly.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1])
            sage: t = phi.ore_variable()
            sage: u = t^2 + (2*z^2 + 3*z + 3)*t + (2*z + 3)
            sage: f = phi.hom(u)
            sage: f._motive_matrix()
            [                      T + 3 + z                 3 + 3*z + 2*z^2]
            [(1 + z + z^2)*T + 3 + 2*z - z^2               T + 2 - z + 2*z^2]
        """
        phi = self.domain()
        phiT = phi.gen()
        r = phiT.degree()
        K = phi.base_over_constants_field()
        S = phi.ore_polring()
        Frob = S.twisting_morphism()
        KT = PolynomialRing(K, name='T')

        # The first row:
        # we write u = u0 + u1*phiT + u2*phiT^2 + ...
        u = self.ore_polynomial()
        us = [ ]
        while not u.is_zero():
            u, ui = u.right_quo_rem(phiT)
            us.append(ui)
        l = len(us)
        row = [KT([us[i][j] for i in range(l)]) for j in range(r)]
        rows = [row]

        # The next rows:
        # each row is obtained from the previous one by
        # applying the semi-linear transformation f |-> t*f
        inv = K(phiT[r]).inverse()
        B = inv * phiT
        T = KT.gen()
        for i in range(1, r):
            twist = [c.map_coefficients(Frob) for c in row]
            row = [(inv*T - B[0]) * twist[-1]]
            row += [twist[j-1] - B[j]*twist[-1] for j in range(1, r)]
            rows.append(row)

        return matrix(KT, rows)

    def norm(self, ideal=True):
        r"""
        Return the norm of this isogeny.

        INPUT:

        - ``ideal`` -- boolean (default: ``True``); if ``True``, return the
          norm as an ideal in the function ring of the Drinfeld modules; if
          ``False``, return the norm as an element in this function ring (only
          relevant for endomorphisms)

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: t = phi.ore_variable()
            sage: f = phi.hom(t + 1)
            sage: f.norm()
            Principal ideal (T + 4) of Univariate Polynomial Ring in T over Finite Field of size 5

        The norm of the Frobenius endomorphism is equal to the characteristic::

            sage: F = phi.frobenius_endomorphism()
            sage: F.norm()
            Principal ideal (T^3 + 3*T + 3) of Univariate Polynomial Ring in T over Finite Field of size 5
            sage: phi.characteristic()
            T^3 + 3*T + 3

        For `a` in the underlying function ring, the norm of the
        endomorphism given by `\phi_a` is `a^r` where `r` is the rank::

            sage: g = phi.hom(T)
            sage: g.norm()
            Principal ideal (T^3) of Univariate Polynomial Ring in T over Finite Field of size 5

            sage: h = phi.hom(T+1)
            sage: h.norm()
            Principal ideal (T^3 + 3*T^2 + 3*T + 1) of Univariate Polynomial Ring in T over Finite Field of size 5

        For endomorphisms, the norm is not an ideal of `A` but it makes
        sense as an actual element of `A`. We can get this element by passing
        in the argument ``ideal=False``::

            sage: phi.hom(2*T).norm(ideal=False)
            3*T^3

            sage: f.norm(ideal=False)
            Traceback (most recent call last):
            ...
            ValueError: norm is defined as an actual element only for endomorphisms
        """
        nu = self._motive_matrix().det()
        # We cast to A
        A = self.domain().function_ring()
        if ideal:
            nu = A([c.in_base() for c in nu.monic().list()])
            return A.ideal(nu)
        elif self.domain() is self.codomain():
            return A([c.in_base() for c in nu.list()])
        else:
            raise ValueError("norm is defined as an actual element only for endomorphisms")

    def dual_isogeny(self):
        r"""
        Return a dual isogeny to this morphism.

        By definition, a dual isogeny of `f : \phi \to \psi` is an
        isogeny `g : \psi \to \phi` such that the composite `g \circ f`
        is the multiplication by a generator of the norm of `f`.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: t = phi.ore_variable()
            sage: f = phi.hom(t + 1)
            sage: f
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> z*t^3 + t^2 + z
              To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^3 + (3*z^2 + 2*z + 2)*t^2 + (2*z^2 + 3*z + 4)*t + z
              Defn: t + 1
            sage: g = f.dual_isogeny()
            sage: g
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^3 + (3*z^2 + 2*z + 2)*t^2 + (2*z^2 + 3*z + 4)*t + z
              To:   Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: z*t^2 + (4*z + 1)*t + z + 4

        We check that `f \circ g` (resp. `g \circ f`) is the multiplication
        by the norm of `f`::

            sage: a = f.norm().gen(); a
            T + 4
            sage: g * f == phi.hom(a)
            True

            sage: psi = f.codomain()
            sage: f * g == psi.hom(a)
            True

        TESTS::

            sage: zero = phi.hom(0)
            sage: zero.dual_isogeny()
            Traceback (most recent call last):
            ...
            ValueError: the dual isogeny of the zero morphism is not defined
        """
        if not self.is_isogeny():
            raise ValueError("the dual isogeny of the zero morphism is not defined")
        nu = self._motive_matrix().det().monic()
        A = self.domain().function_ring()
        nu = A([c.in_base() for c in nu.list()])
        dual = self.domain()(nu) // self.ore_polynomial()
        return self.codomain().hom(dual, codomain=self.domain())

    def characteristic_polynomial(self, var='X'):
        r"""
        Return the characteristic polynomial of this endomorphism.

        INPUT:

        - ``var`` -- string (default: ``X``), the name of the
          variable of the characteristic polynomial

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])

            sage: f = phi.frobenius_endomorphism()
            sage: f.characteristic_polynomial()
            X^3 + (T + 1)*X^2 + (2*T + 3)*X + 2*T^3 + T + 1

        We verify, on an example, that the caracteristic polynomial
        of a morphism corresponding to `\phi_a` is `(X-a)^r` where `r`
        is the rank::

            sage: g = phi.hom(T^2 + 1)
            sage: chi = g.characteristic_polynomial()
            sage: chi.factor()
            (X + 4*T^2 + 4)^3

        An example with another variable name::

            sage: f.characteristic_polynomial(var='Y')
            Y^3 + (T + 1)*Y^2 + (2*T + 3)*Y + 2*T^3 + T + 1

        TESTS::

            sage: t = phi.ore_variable()
            sage: isog = phi.hom(t + 1)
            sage: isog.characteristic_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: characteristic polynomial is only defined for endomorphisms
        """
        if self.domain() is not self.codomain():
            raise ValueError("characteristic polynomial is only defined for endomorphisms")
        P = self._motive_matrix().charpoly()
        # We cast to the correct parent
        A = self.domain().function_ring()
        parent = PolynomialRing(A, name=var)
        return parent([A([c.in_base() for c in co.list()]) for co in P.list()])

    def charpoly(self, var='X'):
        r"""
        Return the characteristic polynomial of this endomorphism.

        INPUT:

        - ``var`` -- string (default: ``'X'``); the name of the
          variable of the characteristic polynomial

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])

            sage: f = phi.frobenius_endomorphism()
            sage: chi = f.charpoly()
            sage: chi
            X^3 + (T + 1)*X^2 + (2*T + 3)*X + 2*T^3 + T + 1

        We check that the characteristic polynomial annihilates the
        morphism (Cayley-Hamilton's theorem)::

            sage: chi(f)
            Endomorphism of Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: 0

        We verify, on an example, that the caracteristic polynomial
        of the morphism corresponding to `\phi_a` is `(X-a)^r` where `r`
        is the rank::

            sage: g = phi.hom(T^2 + 1)
            sage: g.charpoly().factor()
            (X + 4*T^2 + 4)^3

        An example with another variable name::

            sage: f.charpoly(var='Y')
            Y^3 + (T + 1)*Y^2 + (2*T + 3)*Y + 2*T^3 + T + 1
        """
        return self.characteristic_polynomial(var)
