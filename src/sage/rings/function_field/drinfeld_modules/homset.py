# sage.doctest: needs sage.rings.finite_rings
r"""
Set of morphisms between two Drinfeld modules

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.homset.DrinfeldModuleHomset`.

AUTHORS:

- Antoine Leudière (2022-04)
- Xavier Caruso, Yossef Musleh (2025-08): added computation of bases of homsets
"""

# *****************************************************************************
#        Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   https://www.gnu.org/licenses/
# *****************************************************************************

import operator

from sage.categories.finite_fields import FiniteFields
from sage.categories.drinfeld_modules import DrinfeldModules
from sage.categories.homset import Homset
from sage.categories.action import Action
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import Matrix
from sage.matrix.special import identity_matrix
from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
from sage.structure.parent import Parent
from sage.functions.log import logb


class DrinfeldModuleMorphismAction(Action):
    r"""
    Action of the function ring on the homset of a Drinfeld module.

    EXAMPLES::

        sage: Fq = GF(5)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(3)
        sage: phi = DrinfeldModule(A, [z, 1, z])
        sage: psi = DrinfeldModule(A, [z, z^2 + 4*z + 3, 2*z^2 + 4*z + 4])
        sage: H = Hom(phi, psi)
        sage: tau = phi.ore_variable()
        sage: f = H(tau + 2)

    Left action::

        sage: (T + 1) * f
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> z*τ^2 + τ + z
          To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*τ^2 + (z^2 + 4*z + 3)*τ + z
          Defn: (2*z^2 + 4*z + 4)*τ^3 + (2*z + 1)*τ^2 + (2*z^2 + 4*z + 2)*τ + 2*z + 2

    Right action currently does not work (it is a known bug, due to an
    incompatibility between multiplication of morphisms and the coercion
    system)::

        sage: f * (T + 1)
        Traceback (most recent call last):
        ...
        TypeError: right (=T + 1) must be a map to multiply it by Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> z*τ^2 + τ + z
          To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*τ^2 + (z^2 + 4*z + 3)*τ + z
          Defn: τ + 2
    """
    def __init__(self, A, H, is_left, op) -> None:
        r"""
        Initialize this action.

        INPUT:

        - ``A`` -- the function ring of the underlying Drinfeld module

        - ``H`` -- a homset between Drinfeld modules

        - ``is_left`` -- boolean

        - ``op`` -- an operator

        TESTS::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)

            sage: from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleMorphismAction
            sage: left_action = DrinfeldModuleMorphismAction(A, H, True, operator.mul)
            sage: TestSuite(left_action).run(skip='_test_pickling')

            sage: right_action = DrinfeldModuleMorphismAction(A, H, False, operator.mul)
            sage: TestSuite(right_action).run(skip='_test_pickling')
        """
        Action.__init__(self, A, H, is_left, op)
        if is_left:
            self._phi = H.codomain()
        else:
            self._phi = H.domain()

    def _act_(self, a, f):
        r"""
        Return the action of `a` on `f`.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: t = phi.ore_variable()
            sage: f = phi.hom(t + 1)
            sage: T*f  # indirect doctest
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
              To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*τ^3 + (3*z^2 + 2*z + 2)*τ^2 + (2*z^2 + 3*z + 4)*τ + z
              Defn: (2*z^2 + 4*z + 4)*τ^4 + (z + 1)*τ^3 + τ^2 + (2*z^2 + 4*z + 4)*τ + z
        """
        u = f.ore_polynomial()
        if self._is_left:
            u = self._phi(a) * u
        else:
            u = u * self._phi(a)
        return f.parent()(u)


class DrinfeldModuleHomset(Homset):
    r"""
    This class implements the set of morphisms between two Drinfeld
    `\GF{q}[T]`-modules.

    INPUT:

    - ``X`` -- the domain

    - ``Y`` -- the codomain

    EXAMPLES::

        sage: Fq = GF(27)
        sage: A.<T> = Fq[]
        sage: K.<z6> = Fq.extension(2)
        sage: phi = DrinfeldModule(A, [z6, z6, 2])
        sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
        sage: H = Hom(phi, psi)
        sage: H
        Set of Drinfeld module morphisms
         from (gen) 2*τ^2 + z6*τ + z6
           to (gen) 2*τ^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*τ + z6

    ::

        sage: from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        sage: isinstance(H, DrinfeldModuleHomset)
        True

    There is a simpler syntax for endomorphisms sets::

        sage: E = End(phi)
        sage: E
        Set of Drinfeld module morphisms from (gen) 2*τ^2 + z6*τ + z6 to (gen) 2*τ^2 + z6*τ + z6
        sage: E is Hom(phi, phi)
        True

    The domain and codomain must have the same Drinfeld modules
    category::

        sage: rho = DrinfeldModule(A, [T, 1])
        sage: Hom(phi, rho)
        Traceback (most recent call last):
        ...
        ValueError: Drinfeld modules must be in the same category

    ::

        sage: sigma = DrinfeldModule(A, [1, z6, 2])
        sage: Hom(phi, sigma)
        Traceback (most recent call last):
        ...
        ValueError: Drinfeld modules must be in the same category

    One can create morphism objects by calling the homset::

        sage: identity_morphism = E(1)
        sage: identity_morphism
        Identity morphism of Drinfeld module defined by T |--> 2*τ^2 + z6*τ + z6

    ::

        sage: tau = phi.ore_variable()
        sage: frobenius_endomorphism = E(tau^6)
        sage: frobenius_endomorphism
        Endomorphism of Drinfeld module defined by T |--> 2*τ^2 + z6*τ + z6
          Defn: τ^6

    ::

        sage: isogeny = H(tau + 1)
        sage: isogeny
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> 2*τ^2 + z6*τ + z6
          To:   Drinfeld module defined by T |--> 2*τ^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*τ + z6
          Defn: τ + 1

    And one can test if an Ore polynomial defines a morphism using the
    ``in`` syntax::

        sage: 1 in H
        False
        sage: tau^6 in H
        False
        sage: tau + 1 in H
        True
        sage: 1 in E
        True
        sage: tau^6 in E
        True
        sage: tau + 1 in E
        False

    This also works if the candidate is a morphism object::

        sage: isogeny in H
        True
        sage: E(0) in E
        True
        sage: identity_morphism in H
        False
        sage: frobenius_endomorphism in H
        False
    """
    Element = DrinfeldModuleMorphism

    def __init__(self, X, Y, category=None, check=True) -> None:
        """
        Initialize ``self``.

        INPUT:

        - ``X`` -- the domain of the homset

        - ``Y`` -- the codomain of the homset

        - ``category`` -- (default: ``None``) the Drinfeld modules category of
          the domain and codomain

        - ``check`` -- boolean (default: ``True``); check the validity of the
          category

        TESTS::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: H.domain() is phi
            True
            sage: H.codomain() is psi
            True
        """
        if category is None:
            category = X.category()
        if check:
            if X.category() != Y.category() \
                    or not isinstance(X.category(), DrinfeldModules):
                raise ValueError('Drinfeld modules must be in the same category')
            if category != X.category():
                raise ValueError('category should be DrinfeldModules')
        base = category.base()
        super().__init__(X, Y, category=category, base=base, check=check)
        A = X.function_ring()
        self.register_action(DrinfeldModuleMorphismAction(A, self, True, operator.mul))
        # ARGH: the next line does not work
        # because Map.__mul__ does not call the coercion system
        self.register_action(DrinfeldModuleMorphismAction(A, self, False, operator.mul))
        if X is Y:
            self.register_coercion(A)
        self._basis = None

    def _latex_(self) -> str:
        r"""
        Return a LaTeX representation of the homset.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: latex(H)
            \text{Set{ }of{ }Drinfeld{ }module{ }morphisms{ }from{ }(gen){ }}2 τ^{2} + z_{6} τ + z_{6}\text{{ }to{ }(gen){ }}2 τ^{2} + \left(2 z_{6}^{5} + 2 z_{6}^{4} + 2 z_{6} + 1\right) τ + z_{6}
        """
        return f'\\text{{Set{{ }}of{{ }}Drinfeld{{ }}module{{ }}morphisms' \
               f'{{ }}from{{ }}(gen){{ }}}}{latex(self.domain().gen())}' \
               f'\\text{{{{ }}to{{ }}(gen){{ }}}}'\
               f'{latex(self.codomain().gen())}'

    def _repr_(self) -> str:
        r"""
        Return a string representation of the homset.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: H
            Set of Drinfeld module morphisms from (gen) 2*τ^2 + z6*τ + z6 to (gen) 2*τ^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*τ + z6
        """
        return f'Set of Drinfeld module morphisms from (gen) '\
               f'{self.domain().gen()} to (gen) {self.codomain().gen()}'

    def __contains__(self, x) -> bool:
        r"""
        Return ``True`` if the input defines a morphism in the homset.

        INPUT:

        - ``x`` -- an Ore polynomial or a Drinfeld module morphism

        EXAMPLES:

        In the next examples, the input is an Ore polynomial::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: E = End(phi)
            sage: tau = phi.ore_variable()
            sage: 1 in H
            False
            sage: tau^6 in H
            False
            sage: tau + 1 in H
            True
            sage: 1 in E
            True
            sage: tau^6 in E
            True
            sage: tau + 1 in E
            False

        Whereas the input is now a Drinfeld module morphism::

            sage: isogeny = H(tau + 1)
            sage: isogeny in H
            True
            sage: E(0) in E
            True
            sage: identity_morphism = E(1)
            sage: identity_morphism in H
            False
            sage: frobenius_endomorphism = phi.frobenius_endomorphism()
            sage: frobenius_endomorphism in H
            False
        """
        try:
            x = self(x)
            return True
        except (AttributeError, ValueError, TypeError):
            return False

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return the Drinfeld module morphism defined by the given Ore
        polynomial.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: E = End(phi)
            sage: tau = phi.ore_variable()
            sage: identity_morphism = E(1)
            sage: identity_morphism
            Identity morphism of Drinfeld module defined by T |--> 2*τ^2 + z6*τ + z6

        ::

            sage: scalar_multiplication = E(T)
            sage: scalar_multiplication
            Endomorphism of Drinfeld module defined by T |--> 2*τ^2 + z6*τ + z6
              Defn: 2*τ^2 + z6*τ + z6

        ::

            sage: frobenius_endomorphism = E(tau^6)
            sage: frobenius_endomorphism
            Endomorphism of Drinfeld module defined by T |--> 2*τ^2 + z6*τ + z6
              Defn: τ^6

        ::

            sage: isogeny = H(tau + 1)
            sage: isogeny
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> 2*τ^2 + z6*τ + z6
              To:   Drinfeld module defined by T |--> 2*τ^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*τ + z6
              Defn: τ + 1
        """
        # NOTE: This used to be self.element_class(self, ...), but this
        # would call __init__ instead of __classcall_private__. This
        # seems to work, but I don't know what I'm doing.
        return DrinfeldModuleMorphism(self, *args, **kwds)

    def an_element(self, degree=None):
        r"""
        Return an element in this homset.

        INPUT:

        - ``degree`` (default: ``None``) -- an integer; if given,
          return an isogeny of this degree

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: H = Hom(phi, psi)
            sage: H.an_element()
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
              To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
              Defn: z^2*τ^3

        We can also ask for an isogeny with a required degree::

            sage: H.an_element(degree=2)
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
              To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
              Defn: (z^2 + 1)*τ^2 + τ + z + 1

        If there is no isogeny with the required degree, an error is raised::

            sage: H.an_element(degree=1)
            Traceback (most recent call last):
            ...
            ValueError: no isogeny of given degree

        Below, `\phi` and `\psi` are not isogenous, so :meth:`an_element`
        returns the zero morphism (which is the unique element in the
        homset)::

            sage: psi = DrinfeldModule(A, [z, z^2, z^3])
            sage: H = Hom(phi, psi)
            sage: H.an_element()
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
              To:   Drinfeld module defined by T |--> (z + 1)*τ^2 + z^2*τ + z
              Defn: 0

        When the base is not a finite field, computing isogenies is not
        implemented so far, so this method always returns the zero morphism::

            sage: phi = DrinfeldModule(A, [T, 1])
            sage: End(phi).an_element()
            Endomorphism of Drinfeld module defined by T |--> τ + T
              Defn: 0
            sage: End(phi).an_element(degree=1)
            Traceback (most recent call last):
            ...
            NotImplementedError: computing isogenies are currently only implemented over finite fields
        """
        if self.base() not in FiniteFields():
            if degree is None:
                return self.zero()
            else:
                raise NotImplementedError("computing isogenies are currently only implemented over finite fields")
        if degree is None:
            basis = self._A_basis()
            if len(basis) == 0:
                return self.zero()
            return basis[0]
        else:
            basis = self._Fq_basis(degree=degree)
            for isogeny in basis:
                if isogeny.degree() == degree:
                    return isogeny
            raise ValueError("no isogeny of given degree")

    def zero(self):
        r"""
        Return the zero morphism is this homset.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: H = Hom(phi, psi)
            sage: H.zero()
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
              To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
              Defn: 0
        """
        return self(self.domain().ore_polring().zero())

    @cached_method
    def _A_basis(self):
        r"""
        Return a basis of this homset over the underlying
        function ring.

        Currently, it only works over finite fields.

        This method should not be called directly; call
        :meth:`basis` instead.

        ALGORITHM:

        The isogenies `u : \phi \to \psi' correspond to
        the vectors `u \in M(\phi)` such that

        .. MATH::

            \sum_{k=0}^r g_k \tau^k(u) = T u

        where the `g_k` are the coefficients of `\psi_T`

        The algorithm consists in solving the above system
        viewed as a linear system over `\GF{q}[T]`.

        We refer to [Mus2023]_, Section 7.3 for more details.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: End(phi).basis()
            Traceback (most recent call last):
            ...
            NotImplementedError: computing basis of homsets are currently only implemented over finite fields
        """
        phi = self.domain()
        psi = self.codomain()
        r = phi.rank()
        if psi.rank() != r:
            return []
        Fo = phi.base_over_constants_field()
        Fq = Fo.base()
        F = Fo.backend()
        d = Fo.degree(Fq)
        q = Fq.cardinality()
        A = phi.function_ring()
        T = A.gen()

        AF = PolynomialRing(F, name='T')
        TF = AF.gen()

        Frob = lambda x: x**q
        FrobT = lambda P: P.map_coefficients(Frob)

        phiT = phi.gen()
        psiT = psi.gen()

        # We compute the tau^i in M(phi)
        lc = ~(phiT[r])
        xT = [-AF(lc * phiT[i]) for i in range(r)]
        xT[0] += lc * TF
        taus = []
        for i in range(r):
            taui = r * [AF.zero()]
            taui[i] = AF.one()
            taus.append(taui)
        for i in range(r):
            s = FrobT(taui[-1])
            taui = [s * xT[0]] + [FrobT(taui[j - 1]) + s * xT[j] for j in range(1, r)]
            taus.append(taui)

        # We precompute the Frob^k(z^i)
        z = F(Fo.gen())
        zs = []
        zq = z
        for k in range(r + 1):
            x = F.one()
            for i in range(d):
                zs.append(x)
                x *= zq
            zq = zq ** q

        # We compute the linear system to solve
        rows = []
        for i in range(d):
            for j in range(r):
                # For x = z^i * tau^j, we compute
                #    sum(g_k*tau^k(x), k=0..r) - T*x
                #  = sum(g_k*Frob^k(z^i)*tau^(k+j), k=0..r) - T*x
                row = r * [AF.zero()]
                for k in range(r + 1):
                    s = psiT[k] * zs[k * d + i]
                    for l in range(r):
                        row[l] += s * taus[k + j][l]
                row[j] -= zs[i] * TF
                # We write it in the A-basis
                rowA = []
                for c in row:
                    c0 = Fo(c[0]).vector()
                    c1 = Fo(c[1]).vector()
                    rowA += [c0[k] + T * c1[k] for k in range(d)]
                rows.append(rowA)
        M = Matrix(rows)

        # We solve the linear system
        ker = M.minimal_kernel_basis()

        # We reconstruct the isogenies
        isogenies = []
        S = phi.ore_polring()
        t = S.gen()
        for row in range(ker.nrows()):
            u = S.zero()
            for i in range(d):
                for j in range(r):
                    a = ker[row, i * r + j]
                    u += zs[i] * t**j * sum(a[k] * phiT**k for k in range(a.degree() + 1))
            isogenies.append(self(u))

        return isogenies

    def _Fq_basis(self, degree):
        r"""
        Return a `\GF{q}`-basis of the space of morphisms
        in this homset of degree at most ``degree``.

        Currently, it only works over finite fields.

        This method should not be called directly; call
        :meth:`basis` instead.

        INPUT:

        - ``degree`` -- an integer

        ALGORITHM:

        We return the basis of the kernel of a matrix derived from the
        constraint that `f \phi_T = \psi_T f`. See [Wes2022]_ for
        details on this algorithm.

        TESTS::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: H = Hom(phi, psi)
            sage: B = H._Fq_basis(3)
            sage: M = B[0]
            sage: M_poly = M.ore_polynomial()
            sage: M_poly*phi.gen() - psi.gen()*M_poly
            0

        ::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z, 1])
            sage: psi = DrinfeldModule(A, [z, 1, 1])
            sage: hom = Hom(phi, psi)
            sage: hom._Fq_basis(4)
            []
        """
        domain, codomain = self.domain(), self.codomain()
        Fq = domain._Fq
        K = domain.base_over_constants_field()
        q = Fq.cardinality()
        r = domain.rank()
        if codomain.rank() != r:
            return []
        n = K.degree(Fq)
        # shorten name for readability
        d = degree
        K_basis = K.basis_over(Fq)
        phiT = domain.coefficients(sparse=False)
        psiT = codomain.coefficients(sparse=False)

        # We precompute the matrices of the iterates of
        # the Frobenius of K/Fq
        frob_matrices = [identity_matrix(Fq, n)] + [Matrix(Fq, n) for _ in range(d + r)]
        for i, elem in enumerate(K_basis):
            for k in range(1, d + r + 1):
                elem = elem ** q
                v = elem.vector()
                for j in range(n):
                    frob_matrices[k][i, j] = v[j]

        # We write the linear system and solve it
        sys = Matrix(Fq, (d + r + 1) * n, (d + 1) * n)
        for k in range(0, d + r + 1):
            for i in range(max(0, k - r), min(k, d) + 1):
                # We represent multiplication and Frobenius
                # as operators acting on K as a vector space
                # over Fq
                oper = K(phiT[k - i] ** (q**i)).matrix() \
                     - frob_matrices[k - i] * K(psiT[k - i]).matrix()
                for j in range(n):
                    for l in range(n):
                        sys[k * n + j, i * n + l] = oper[l, j]
        sol = sys.right_kernel().basis()

        # Reconstruct the Ore polynomial from the coefficients
        basis = []
        tau = domain.ore_polring().gen()
        for basis_elem in sol:
            ore_poly = sum([sum([K_basis[j].backend() * basis_elem[i * n + j]
                               for j in range(n)]) * (tau**i)
                               for i in range(d + 1)])
            basis.append(self(ore_poly))

        return basis

    def basis(self, degree=None):
        r"""
        Return a basis of this homset.

        INPUT:

        - ``degree`` -- an integer or ``None`` (default: ``None``)

        If ``degree`` is ``None``, a basis over the underlying
        function ring is returned.
        Otherwise, a `\GF{q}`-basis of the set of morphisms of
        degree at most ``degree`` is returned.

        ALGORITHM:

        We reformulate the problem as a linear system over `\GF{q}`
        or `A = \GF{q}[T]` and solve it.
        We refer to [Wes2022]_ and [Mus2023]_, Section 7.3 for more details.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: End(phi).basis()
            [Identity morphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z,
             Endomorphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               Defn: τ^2,
             Endomorphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               Defn: 2*τ^4 + z*τ^3 + z]

        If we specify a degree, a basis over `\GF{q}` is computed::

            sage: End(phi).basis(degree=5)
            [Identity morphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z,
             Endomorphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               Defn: z*τ^3 + z,
             Endomorphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               Defn: τ^2,
             Endomorphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               Defn: z*τ^5 + z*τ^2,
             Endomorphism of Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               Defn: τ^4]

        Here is another example where the domain and the codomain differ::

            sage: psi = DrinfeldModule(A, [z, 3*z + 1, 2*z, 4*z + 1])
            sage: H = Hom(phi, psi)
            sage: H.basis()
            [Drinfeld Module morphism:
               From: Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               To:   Drinfeld module defined by T |--> (4*z + 1)*τ^3 + 2*z*τ^2 + (3*z + 1)*τ + z
               Defn: τ + 1,
             Drinfeld Module morphism:
               From: Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               To:   Drinfeld module defined by T |--> (4*z + 1)*τ^3 + 2*z*τ^2 + (3*z + 1)*τ + z
               Defn: 3*τ^3 + (z + 2)*τ^2 + 4*z*τ + z + 4,
             Drinfeld Module morphism:
               From: Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               To:   Drinfeld module defined by T |--> (4*z + 1)*τ^3 + 2*z*τ^2 + (3*z + 1)*τ + z
               Defn: (z + 4)*τ^2 + 4*z*τ + z + 4]

            sage: H.basis(degree=2)
            [Drinfeld Module morphism:
             From: Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               To:   Drinfeld module defined by T |--> (4*z + 1)*τ^3 + 2*z*τ^2 + (3*z + 1)*τ + z
             Defn: τ + 1,
             Drinfeld Module morphism:
             From: Drinfeld module defined by T |--> z*τ^3 + τ^2 + z
               To:   Drinfeld module defined by T |--> (4*z + 1)*τ^3 + 2*z*τ^2 + (3*z + 1)*τ + z
             Defn: (z + 4)*τ^2 + (4*z + 1)*τ + z]

        We can check that `\phi` and `\psi` are not isomorphic by checking
        that there is no isogeny of degree `0` between them::

            sage: H.basis(degree=0)
            []

        When `\phi` and `\psi` are not isogenous, an empty list is returned::

            sage: psi = DrinfeldModule(A, [z, 3*z, 2*z, 4*z])
            sage: Hom(phi, psi).basis()
            []

        Currently, this method only works over finite fields::

            sage: A.<T> = GF(5)[]
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: End(phi).basis()
            Traceback (most recent call last):
            ...
            NotImplementedError: computing basis of homsets are currently only implemented over finite fields
        """
        if self.base() not in FiniteFields():
            raise NotImplementedError("computing basis of homsets are currently only implemented over finite fields")
        if degree is None:
            return self._A_basis()
        else:
            return self._Fq_basis(degree)

    def basis_over_frobenius(self):
        r"""
        Return a basis of this homser over `\GF{q}[\tau^n]` where
        `n = [K:\GF{q}]` (and thus `\tau^n` is to the Frobenius endomorphism).

        ALGORITHM:

        We return the basis of the kernel of a matrix derived from the
        constraint that `\iota \phi_T = \psi_T \iota` for any morphism
        `iota: \phi \to \psi`.
        We refer to [Mus2023]_, Section 7.3 for more details.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: H = Hom(phi, psi)
            sage: H.basis_over_frobenius()
            [Drinfeld Module morphism:
               From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
               To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
               Defn: (z^2 + 1)*τ^2 + τ + z + 1,
             Drinfeld Module morphism:
               From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
               To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
               Defn: (z^2 + 1)*τ^5 + (z + 1)*τ^4 + z*τ^3 + τ^2 + (z^2 + z)*τ + z,
             Drinfeld Module morphism:
               From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
               To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
               Defn: z^2]

        ::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 3*z, 4*z])
            sage: chi = DrinfeldModule(A, [z, 2*z^2 + 3, 4*z^2 + 4*z])
            sage: H = Hom(phi, chi)
            sage: H.basis_over_frobenius()
            []

        TESTS::

            sage: Fq = GF(4)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z^5 + z^3 + 1, 1])
            sage: phi = DrinfeldModule(A, [z, z^4 + z^3 + 1, 1])
            sage: H = Hom(phi, psi)
            sage: basis = H.basis_over_frobenius()
            sage: basis
            [... Defn: τ^2 + (z^5 + z^4 + z^3 + z^2)*τ + z^5 + z^4 + z^3 + 1,
             ... Defn: (z^3 + z^2 + z + 1)*τ^2 + (z^5 + z^4 + z^3 + 1)*τ + z^4 + z^3 + z,
             ... Defn: (z^3 + z^2 + z + 1)*τ^5 + (z^5 + z^4 + z^3)*τ^4 + z^2*τ^3 + (z^3 + 1)*τ^2 + (z^5 + z^4 + z^2 + 1)*τ + z^2]
            sage: basis[2].ore_polynomial()*phi.gen() - psi.gen()*basis[2].ore_polynomial()
            0

        ::

            sage: A.<T> = GF(5)[]
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: End(phi).basis_over_frobenius()
            Traceback (most recent call last):
            ...
            ValueError: basis over Frobenius only makes sense for Drinfeld module defined over finite fields
        """
        if self.base() not in FiniteFields():
            raise ValueError("basis over Frobenius only makes sense for Drinfeld module defined over finite fields")
        Fq = self.domain()._Fq
        K = self.domain().base_over_constants_field()
        r = self.domain().rank()
        n = K.degree(Fq)
        q = Fq.cardinality()
        K_basis = K.basis_over(Fq)
        phiT = self.domain().coefficients(sparse=False)
        psiT = self.codomain().coefficients(sparse=False)

        # The commutative polynomial ring in tau^n.
        poly_taun = PolynomialRing(Fq, 'taun')
        taun = poly_taun.gen()

        sys = Matrix(poly_taun, n**2, n**2)

        # Build a linear system over the commutative polynomial ring
        # Fq[tau^n]. The kernel of this system consists of all
        # morphisms from domain to codomain.
        for j in range(n):
            for k in range(n):
                for i in range(r + 1):
                    # Coefficients of tau^{i + k} coming from the
                    # relation defining morphisms of Drinfeld modules
                    # These are elements of K, expanded in terms of
                    # K_basis.
                    poly = K(phiT[i]**(q**k) * K_basis[j]
                           - psiT[i] * K_basis[j]**(q**i)).polynomial()
                    deg = (i + k) // n
                    row = n * (i + k - n * deg)
                    col = k * n + j
                    for b in range(poly.degree() + 1):
                        sys[row + b, col] += poly[b] * taun**deg

        sol = sys.right_kernel().basis()

        # Reconstruct basis as skew polynomials.
        basis = []
        tau = self.domain().ore_variable()
        for basis_vector in sol:
            basis_poly = 0
            for i in range(n):
                for j in range(n):
                    basis_poly += basis_vector[n * i + j].subs(tau**n) * K_basis[j].backend() * tau**i
            basis.append(self(basis_poly))

        return basis

    def random_element(self, degree=None):
        r"""
        Return a random morphism in this homset.

        INPUT:

        - ``degree`` (default: ``None``) -- the maximum degree of
          the morphism

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: H = Hom(phi, psi)
            sage: H.random_element()  # random
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
              To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
              Defn: z*τ^7 + (z^2 + 1)*τ^6 + τ^5 + z^2*τ^4 + (z^2 + z + 1)*τ^3 + τ^2 + (z^2 + z)*τ + z

        When ``degree`` is given, a uniformly distributed random isogeny
        of degree *at most* the given value is outputted::

            sage: H.random_element(3)  # random
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
              To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
              Defn: (z^2 + 1)*τ^2 + τ + z + 1

        For producing a random isogeny with accurate degree, we can
        proceed as follows::

            sage: H.an_element(3) + H.random_element(2)  # random
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> (z^2 + z)*τ^2 + (z^2 + z + 1)*τ + z
              To:   Drinfeld module defined by T |--> (z^2 + z + 1)*τ^2 + (z + 1)*τ + z
              Defn: z^2*τ^3 + (z^2 + 1)*τ^2 + τ + z^2 + z + 1

        Currently, this method only works over finite fields::

            sage: phi = DrinfeldModule(A, [T, 1])
            sage: End(phi).random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError: computing isogenies are currently only implemented over finite fields
        """
        if self.base() not in FiniteFields():
            raise NotImplementedError("computing isogenies are currently only implemented over finite fields")
        domain = self.domain()
        if degree is None:
            scalars = domain._function_ring
            basis = self._A_basis()
        else:
            scalars = domain._Fq
            basis = self._Fq_basis(degree)
        isog = self.zero()
        for u in basis:
            isog += scalars.random_element() * u
        return isog
