# sage.doctest: optional - sage.rings.finite_rings
r"""
Set of morphisms between two Drinfeld modules

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.homset.DrinfeldModuleHomset`.

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

import operator

from sage.categories.drinfeld_modules import DrinfeldModules
from sage.categories.homset import Homset
from sage.categories.action import Action
from sage.misc.latex import latex
from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
from sage.structure.parent import Parent
from sage.functions.log import logb
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.randstate import set_random_seed


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
        sage: t = phi.ore_variable()
        sage: f = H(t + 2)

    Left action::

        sage: (T + 1) * f
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> z*t^2 + t + z
          To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^2 + (z^2 + 4*z + 3)*t + z
          Defn: (2*z^2 + 4*z + 4)*t^3 + (2*z + 1)*t^2 + (2*z^2 + 4*z + 2)*t + 2*z + 2

    Right action currently does not work (it is a known bug, due to an
    incompatibility between multiplication of morphisms and the coercion
    system)::

        sage: f * (T + 1)
        Traceback (most recent call last):
        ...
        TypeError: right (=T + 1) must be a map to multiply it by Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> z*t^2 + t + z
          To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^2 + (z^2 + 4*z + 3)*t + z
          Defn: t + 2

    """
    def __init__(self, A, H, is_left, op):
        r"""
        Initialize this action.

        INPUT:

        - ``A`` -- the function ring of the underlying Drinfeld module

        - ``H`` -- a homset between Drinfeld modules

        - ``is_left`` -- a boolean

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
              From: Drinfeld module defined by T |--> z*t^3 + t^2 + z
              To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^3 + (3*z^2 + 2*z + 2)*t^2 + (2*z^2 + 3*z + 4)*t + z
              Defn: (2*z^2 + 4*z + 4)*t^4 + (z + 1)*t^3 + t^2 + (2*z^2 + 4*z + 4)*t + z

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
    `\mathbb{F}_q[T]`-modules.

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
         from (gen) 2*t^2 + z6*t + z6
           to (gen) 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6

    ::

        sage: from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        sage: isinstance(H, DrinfeldModuleHomset)
        True

    There is a simpler syntax for endomorphisms sets::

        sage: E = End(phi)
        sage: E
        Set of Drinfeld module morphisms from (gen) 2*t^2 + z6*t + z6 to (gen) 2*t^2 + z6*t + z6
        sage: E is Hom(phi, phi)
        True

    The domain and codomain must have the same Drinfeld modules
    category::

        sage: rho = DrinfeldModule(A, [Frac(A)(T), 1])
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
        Identity morphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6

    ::

        sage: t = phi.ore_polring().gen()
        sage: frobenius_endomorphism = E(t^6)
        sage: frobenius_endomorphism
        Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
          Defn: t^6

    ::

        sage: isogeny = H(t + 1)
        sage: isogeny
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
          To:   Drinfeld module defined by T |--> 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
          Defn: t + 1

    And one can test if an Ore polynomial defines a morphism using the
    ``in`` syntax::

        sage: 1 in H
        False
        sage: t^6 in H
        False
        sage: t + 1 in H
        True
        sage: 1 in E
        True
        sage: t^6 in E
        True
        sage: t + 1 in E
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

    def __init__(self, X, Y, category=None, check=True):
        """
        Initialize ``self``.

        INPUT:

        - ``X`` -- the domain of the homset

        - ``Y`` -- the codomain of the homset

        - ``category`` (default: ``None``) -- the Drinfeld modules category of
          the domain and codomain

        - ``check`` (default: ``True``) -- check the validity of the category

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

    def _latex_(self):
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
            \text{Set{ }of{ }Drinfeld{ }module{ }morphisms{ }from{ }(gen){ }}2 t^{2} + z_{6} t + z_{6}\text{{ }to{ }(gen){ }}2 t^{2} + \left(2 z_{6}^{5} + 2 z_{6}^{4} + 2 z_{6} + 1\right) t + z_{6}
        """
        return f'\\text{{Set{{ }}of{{ }}Drinfeld{{ }}module{{ }}morphisms' \
               f'{{ }}from{{ }}(gen){{ }}}}{latex(self.domain().gen())}' \
               f'\\text{{{{ }}to{{ }}(gen){{ }}}}'\
               f'{latex(self.codomain().gen())}'

    def _repr_(self):
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
            Set of Drinfeld module morphisms from (gen) 2*t^2 + z6*t + z6 to (gen) 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
        """
        return f'Set of Drinfeld module morphisms from (gen) '\
               f'{self.domain().gen()} to (gen) {self.codomain().gen()}'

    def __contains__(self, x):
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
            sage: t = phi.ore_polring().gen()
            sage: 1 in H
            False
            sage: t^6 in H
            False
            sage: t + 1 in H
            True
            sage: 1 in E
            True
            sage: t^6 in E
            True
            sage: t + 1 in E
            False

        Whereas the input is now a Drinfeld module morphism::

            sage: isogeny = H(t + 1)
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
            sage: t = phi.ore_polring().gen()
            sage: identity_morphism = E(1)
            sage: identity_morphism
            Identity morphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6

        ::

            sage: scalar_multiplication = E(T)
            sage: scalar_multiplication
            Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              Defn: 2*t^2 + z6*t + z6

        ::

            sage: frobenius_endomorphism = E(t^6)
            sage: frobenius_endomorphism
            Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              Defn: t^6

        ::

            sage: isogeny = H(t + 1)
            sage: isogeny
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              To:   Drinfeld module defined by T |--> 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
              Defn: t + 1
        """
        # NOTE: This used to be self.element_class(self, ...), but this
        # would call __init__ instead of __classcall_private__. This
        # seems to work, but I don't know what I'm doing.
        return DrinfeldModuleMorphism(self, *args, **kwds)

    def element(self, degree):
        r"""
        Return an element of the space of morphisms between the domain and
        codomain. By default, chooses an element of largest degree less than
        or equal to the parameter `degree`.

        INPUT:

        - ``degree`` -- the maximum degree of the morphism

        OUTPUT: a univariate ore polynomials with coefficients in `K`

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: H = Hom(phi, psi)
            sage: M = H.element(3)
            sage: M_poly = M.ore_polynomial()
            sage: M_poly*phi.gen() - psi.gen()*M_poly
            0

        ALGORITHM:

            We scan the basis for the first element of maximal degree
            and return it.
        """
        basis = self.Fq_basis(degree)
        elem = basis[0]
        max_deg = elem.ore_polynomial().degree()
        for basis_elem in basis:
            if basis_elem.ore_polynomial().degree() > max_deg:
                elem = basis_elem
                max_deg = elem.ore_polynomial().degree()
        return elem

    def Fq_basis(self, degree):
        r"""
        Return a basis for the `\mathbb{F}_q`-space of morphisms from `phi` to
        a Drinfeld module `\psi` of degree at most `degree`. A morphism
        `\iota: \phi \to psi` is an element `\iota \in K\{\tau\}` such that
        `iota \phi_T = \psi_T \iota`. The degree of a morphism is the
        `\tau`-degree of `\iota`.

        INPUT:

        - ``degree`` -- the maximum degree of the morphisms in the span.

        OUTPUT: a list of Drinfeld module morphisms.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: H = Hom(phi, psi)
            sage: B = H.Fq_basis(3)
            sage: M = B[0]
            sage: M_poly = M.ore_polynomial()
            sage: M_poly*phi.gen() - psi.gen()*M_poly
            0

        ALGORITHM:

            We return the basis of the kernel of a matrix derived from the
            constraint that `\iota \phi_T = \psi_T \iota`. See [Wes2022]_ for
            details on this algorithm.
        """
        domain, codomain = self.domain(), self.codomain()
        Fq = domain._Fq
        K = domain.base_over_constants_field()
        q = Fq.cardinality()
        char = Fq.characteristic()
        r = domain.rank()
        n = K.degree(Fq)
        # shorten name for readability
        d = degree
        qorder = logb(q, char)
        K_basis = K.basis_over(Fq)
        dom_coeffs = domain.coefficients(sparse=False)
        cod_coeffs = codomain.coefficients(sparse=False)

        sys = Matrix(Fq, (d + r + 1)*n, (d + 1)*n)
        for k in range(0, d + r + 1):
            for i in range(max(0, k - r), min(k, d) + 1):
                # We represent multiplication and Frobenius
                # as operators acting on K as a vector space
                # over Fq
                # Require matrices act on the right, so we
                # take a transpose of operators here
                oper = K(dom_coeffs[k-i]
                       .frobenius(qorder*i)).matrix().transpose() \
                       - K(cod_coeffs[k-i]).matrix().transpose() \
                       * self._frobenius_matrix(k - i)
                for j in range(n):
                    for l in range(n):
                        sys[k*n + j, i*n + l] = oper[j, l]
        sol = sys.right_kernel().basis()
        # Reconstruct the Ore polynomial from the coefficients
        basis = []
        tau = domain.ore_polring().gen()
        for basis_elem in sol:
            basis.append(self(sum([sum([K_basis[j]*basis_elem[i*n + j]
                               for j in range(n)])*(tau**i)
                               for i in range(d + 1)])))
        return basis

    def basis(self):
        r"""
        Return a basis for the `\mathbb{F}_q[\tau^n]`-module of morphisms from
        the domain to the codomain.

        OUTPUT: a list of Drinfeld module morphisms.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: H = Hom(phi, psi)
            sage: basis = H.basis()
            sage: [b.ore_polynomial() for b in basis]
            [(z^2 + 1)*t^2 + t + z + 1, (z^2 + 1)*t^5 + (z + 1)*t^4 + z*t^3 + t^2 + (z^2 + z)*t + z, z^2]

            ::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 3*z, 4*z])
            sage: chi = DrinfeldModule(A, [z, 2*z^2 + 3, 4*z^2 + 4*z])
            sage: H = Hom(phi, chi)
            sage: H.basis()
            []

        ALGORITHM:

            We return the basis of the kernel of a matrix derived from the
            constraint that `\iota \phi_T = \psi_T \iota` for any morphism
            `iota`. 
        """
        domain, codomain = self.domain(), self.codomain()
        Fq = domain._Fq
        K = domain.base_over_constants_field()
        q = Fq.cardinality()
        char = Fq.characteristic()
        r = domain.rank()
        n = K.degree(Fq)
        qorder = logb(q, char)
        K_basis = [K.gen()**i for i in range(n)]
        dom_coeffs = domain.coefficients(sparse=False)
        cod_coeffs = codomain.coefficients(sparse=False)
        # The commutative polynomial ring in tau^n.
        poly_taun = PolynomialRing(Fq, 'taun')

        sys = Matrix(poly_taun, n**2, n**2)

        # Build a linear system over the commutative polynomial ring
        # Fq[tau^n]. The kernel of this system consists of all
        # morphisms from domain -> codomain.
        for j in range(n):
            for k in range(n):
                for i in range(r+1):
                    # Coefficients of tau^{i + k} coming from the
                    # relation defining morphisms of Drinfeld modules
                    # These are elements of K, expanded in terms of
                    # K_basis.
                    c_tik = (dom_coeffs[i].frobenius(qorder*k)*K_basis[j] \
                            - cod_coeffs[i]*K_basis[j].frobenius(qorder*i)) \
                            .polynomial().coefficients(sparse=False)
                    c_tik += [0 for _ in range(n - len(c_tik))]
                    colpos = k*n + j
                    taudeg = i + k
                    for b in range(n):
                        sys[(taudeg % n)*n + b, colpos] += c_tik[b] * \
                                poly_taun.gen()**(taudeg // n)

        sol = sys.right_kernel().basis()
        # Reconstruct basis as skew polynomials.
        basis = []
        tau = domain.ore_polring().gen()
        for basis_vector in sol:
            basis_poly = 0
            for i in range(n):
                for j in range(n):
                    basis_poly += basis_vector[n*i + j].subs(tau**n)*K_basis[j]*tau**i
            basis.append(self(basis_poly))
        return basis
    

    def motive_power_decomposition(self, upper):
        r"""
	Computes coefficients a_0, .., a_{r-1} such that
	\tau^t = a_{r-1}\tau^{r-1} + .... + a_{0}
	with each a_i \in Fq[x] for each t < upper..
	"""
        domain = self.domain()
        r = domain.rank()
        A = domain.function_ring()
        x = A.gen()
        Fq = domain._Fq
        char, q = Fq.characteristic(), Fq.cardinality()
        qord = char.log(q)
        K = domain.base_over_constants_field()
        coeff = domain.coefficients(sparse=False)
        recurrence_relation = [K(coeff[i]/coeff[r]) for i in range(r)]
#	base_expansions = matrix.identity(r)
        expansion_list = [[K(1) if i == j else K(0) for j in range(r)] for i in range(r)]
        for i in range(0, upper - r):
            #print(f'x: {x}')
            z = K.gen()
            #print(f'genr: {z}')
            #print(f'multi: {x*z}')
            #print(f'{[expansion_list[i][k] for k in range(r)]}')
            #print(f'expansions: {[K(expansion_list[i][k])*x for k in range(r)]}')
            next_expansion = [-1*sum([recurrence_relation[j].frobenius(i*qord)*expansion_list[i+j][k] \
                                for j in range(r)]) + expansion_list[i][k]*x for k in range(r)]
            expansion_list.append(next_expansion)
        return expansion_list
				

    def motive_basis(self):
        r"""
	
        """
        domain, codomain = self.domain(), self.codomain()
        Fq = domain._Fq
        K = domain.base_over_constants_field()
        q = Fq.cardinality()
        char = Fq.characteristic()
        r = domain.rank()
        n = K.degree(Fq)
        qorder = logb(q, char)


    def _frobenius_matrix(self, order=1, K_basis=None):
        r"""
        Internal method for computing the matrix of the Frobenius endomorphism
        for K/Fq. This is a useful method for computing morphism ring bases so
        we provide a helper method here. This should probably be part of the
        Finite field implementation.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: H = Hom(phi, psi)
            sage: m = H._frobenius_matrix()
            sage: e = K(z^2 + z + 1)
            sage: frob = K.frobenius_endomorphism()
            sage: frob(e) == K(m*vector(e.polynomial().coefficients(sparse=False)))
            True
        """
        Fq = self.domain()._Fq
        K = self.domain().base_over_constants_field()
        n = K.degree(Fq)
        frob = K.frobenius_endomorphism(order)
        if K_basis is None:
            K_basis = [K.gen()**i for i in range(n)]
        pol_var = K_basis[0].polynomial().parent().gen()
        pol_ring = PolynomialRing(Fq, str(pol_var))
        frob_matrix = Matrix(Fq, n, n)
        for i, elem in enumerate(K_basis):
            col = pol_ring(frob(elem).polynomial()).coefficients(sparse=False)
            col += [0 for _ in range(n - len(col))]
            for j in range(n):
                frob_matrix[j, i] = col[j]
        return frob_matrix

    def random_element(self, degree, seed=None):
        r"""
        Return a random morphism chosen uniformly from the space of morphisms
        of degree at most `degree`.

        INPUT:

        - ``degree`` -- the maximum degree of the morphism

        OUTPUT: a univariate ore polynomials with coefficients in `K`

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: H = Hom(phi, psi)
            sage: M = H.random_element(3, seed=25)
            sage: M_poly = M.ore_polynomial()
            sage: M_poly*phi.gen() - psi.gen()*M_poly
            0
        """
        set_random_seed(seed)
        return self(sum([self.domain()._Fq.random_element()
                * elem.ore_polynomial() for elem in self.Fq_basis(degree)]))
