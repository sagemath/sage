# sage.doctest: needs sage.rings.finite_rings
r"""
Finite Drinfeld modules

This module provides the class
:class:`sage.rings.function_fields.drinfeld_module.finite_drinfeld_module.DrinfeldModule_finite`,
which inherits
:class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

AUTHORS:

- Antoine Leudière (2022-04)
- Yossef Musleh (2023-02): added characteristic polynomial methods
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

from sage.functions.log import logb
from sage.functions.other import ceil, sqrt
from sage.matrix.constructor import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.special import companion_matrix
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class DrinfeldModule_finite(DrinfeldModule):
    r"""
    This class implements finite Drinfeld `\mathbb{F}_q[T]`-modules.

    A *finite Drinfeld module* is a Drinfeld module whose base field is
    finite. In this case, the function field characteristic is a prime
    ideal.

    For general definitions and help on Drinfeld modules, see class
    :class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

    .. RUBRIC:: Construction:

    The user does not ever need to directly call
    ``DrinfeldModule_finite`` --- the metaclass ``DrinfeldModule`` is
    responsible for instantiating ``DrinfeldModule`` or
    ``DrinfeldModule_finite`` depending on the input::

        sage: Fq = GF(343)
        sage: A.<T> = Fq[]
        sage: K.<z6> = Fq.extension(2)
        sage: phi = DrinfeldModule(A, [z6, 0, 5])
        sage: phi
        Drinfeld module defined by T |--> 5*t^2 + z6

    ::

        sage: isinstance(phi, DrinfeldModule)
        True
        sage: from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import DrinfeldModule_finite
        sage: isinstance(phi, DrinfeldModule_finite)
        True

    The user should never use ``DrinfeldModule_finite`` to test if a
    Drinfeld module is finite, but rather the ``is_finite`` method::

        sage: phi.is_finite()
        True

    .. RUBRIC:: Complex multiplication of rank two finite Drinfeld modules

    We can handle some aspects of the theory of complex multiplication
    of finite Drinfeld modules. Apart from the method
    ``frobenius_endomorphism``, we only handle rank two Drinfeld
    modules.

    First of all, it is easy to create the Frobenius endomorphism::

        sage: frobenius_endomorphism = phi.frobenius_endomorphism()
        sage: frobenius_endomorphism
        Endomorphism of Drinfeld module defined by T |--> 5*t^2 + z6
          Defn: t^2

    Its characteristic polynomial can be computed::

        sage: chi = phi.frobenius_charpoly()
        sage: chi
        X^2 + (T + 2*z3^2 + 2*z3 + 1)*X + 2*T^2 + (z3^2 + z3 + 4)*T + 2*z3
        sage: frob_pol = frobenius_endomorphism.ore_polynomial()
        sage: chi(frob_pol, phi(T))
        0

    as well as its trace and norm::

        sage: phi.frobenius_trace()
        6*T + 5*z3^2 + 5*z3 + 6
        sage: phi.frobenius_trace() == -chi[1]
        True
        sage: phi.frobenius_norm()
        2*T^2 + (z3^2 + z3 + 4)*T + 2*z3

    We can decide if a Drinfeld module is ordinary or supersingular::

        sage: phi.is_ordinary()
        True
        sage: phi.is_supersingular()
        False

    .. RUBRIC:: Inverting the Drinfeld module

    The morphism that defines a Drinfeld module is injective (see
    [Gos1998]_, cor. 4.5.2). If the Drinfeld module is finite, one can
    retrieve preimages::

        sage: a = A.random_element()
        sage: phi.invert(phi(a)) == a
        True
    """

    def __init__(self, gen, category):
        """
        Initialize ``self``.

        Validity of the input is checked in `__classcall_private__`. The
        `__init__` just saves attributes.

        INPUT:

        - ``function_ring`` -- a univariate polynomial ring whose base
          is a finite field

        - ``gen`` -- the generator of the Drinfeld module as a list of
          coefficients or an Ore polynomial

        - ``name`` -- (default: ``'t'``) the name of the Ore polynomial
          ring gen

        TESTS::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: gen = [p_root, z12^3, z12^5]
            sage: phi = DrinfeldModule(A, gen)
            sage: ore_polring = phi.ore_polring()
            sage: phi._gen == ore_polring(gen)
            True
        """
        # NOTE: There used to be no __init__ here (which was fine). I
        # added one to ensure that DrinfeldModule_finite would always
        # have _frobenius_norm and _frobenius_trace attributes.
        super().__init__(gen, category)
        self._base_degree_over_constants = self.base_over_constants_field().degree(self._Fq)
        self._frobenius_norm = None
        self._frobenius_trace = None
        self._frobenius_charpoly = None

    def _frobenius_crystalline_matrix(self):
        r"""
        Return the matrix representing the Frobenius endomorphism on the
        crystalline cohomology of the Drinfeld module. This is done up to
        precision [K:Fq] + 1, which is enough to ensure the characteristic
        polynomial can be recovered.

        For the formal construction of the crystalline cohomology, see
        [Ang1997]_. It is a free module of rank r over the ring of Witt
        vectors in `T - \mathfrak{p}`.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi._frobenius_crystalline_matrix()
            [...((z2 + 3) + z12 + (4*z2 + 1)*z12^2 + (z2 + 4)*z12^3 + (z2 + 4)*z12^4 + (2*z2 + 3)*z12^5)*T^3 + (2*z2 + 2*z2*z12 + (2*z2 + 3)*z12^2 + ... + (3*z2 + 4)*z12^3 + (2*z2 + 3)*z12^4 + 3*z2*z12^5]

        ALGORITHM:

        Compute the action of the Frobenius on an explicit basis for the
        cohomology using a recurrence relation, and return the resulting
        matrix.
        """
        Fq = self._Fq
        K = self.base_over_constants_field()
        A = self.function_ring()
        char, q = Fq.characteristic(), Fq.cardinality()
        qdeg = logb(q, char)
        r, n = self.rank(), self._base_degree_over_constants
        nstar = ceil(sqrt(n))
        nquo, nrem = divmod(n, nstar)
        drin_coeffs = self.coefficients(sparse=False)
        poly_K = PolynomialRing(K, name=str(A.gen()))
        matrix_poly_K = MatrixSpace(poly_K, r, r)
        mu_coeffs = ((poly_K.gen() - drin_coeffs[0])**(n+1)) \
                    .coefficients(sparse=False)

        def companion(order):
            # + [1] is required to satisfy formatting for companion matrix
            M = matrix_poly_K(companion_matrix([(drin_coeffs[i]/drin_coeffs[r])
                               .frobenius(qdeg*order)
                               for i in range(r)] + [1], format='top'))
            M[0, r-1] += poly_K.gen() / drin_coeffs[r].frobenius(qdeg*order)
            return M

        companion_initial = prod([companion(i) for i in range(nrem, 0, -1)])
        companion_step = prod([companion(i)
                              for i in range(nstar + nrem, nrem, -1)])
        reduced_companions = []
        for k in range(nquo - 1, 0, -1):
            M = Matrix(poly_K, r, r)
            modulus = poly_K([c.frobenius(qdeg*(-k*nstar % n))
                                                for c in mu_coeffs])
            for i, row in enumerate(companion_step):
                for j, entry in enumerate(row):
                    reduction = entry % modulus
                    M[i, j] = poly_K([c.frobenius(qdeg*(k*nstar))
                                     for c in reduction
                                              .coefficients(sparse=False)])
            reduced_companions.append(M)
        return (prod(reduced_companions) * companion_step * companion_initial)

    def frobenius_endomorphism(self):
        r"""
        Return the Frobenius endomorphism of the Drinfeld module as a
        morphism object.

        Let `q` be the order of the base field of the function ring. The
        *Frobenius endomorphism* is defined as the endomorphism whose
        defining Ore polynomial is `t^q`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.frobenius_endomorphism()
            Endomorphism of Drinfeld module defined by T |--> z6*t^2 + 1
              Defn: t^2

        TESTS::

            sage: from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
            sage: isinstance(phi.frobenius_endomorphism(), DrinfeldModuleMorphism)
            True
        """
        t = self.ore_polring().gen()
        deg = self.base_over_constants_field().degree_over()
        return self._Hom_(self, category=self.category())(t**deg)

    def frobenius_charpoly(self, var='X', algorithm='crystalline'):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism.

        Let `\mathbb{F}_q` be the base field of the function ring. The
        *characteristic polynomial* `\chi` *of the Frobenius endomorphism*
        is defined in [Gek1991]_. An important feature of this
        polynomial is that it is monic, univariate, and has coefficients
        in the function ring. As in our case the function
        ring is a univariate polynomial ring, it is customary to see the
        characteristic polynomial of the Frobenius endomorphism as a
        bivariate polynomial.

        Let `\chi = X^r + \sum_{i=0}^{r-1} A_{i}(T)X^{i}` be the
        characteristic polynomial of the Frobenius endomorphism, and
        let `t^n` be the Ore polynomial that defines the Frobenius
        endomorphism of `\phi`; by definition, `n` is the degree of `K`
        over the base field `\mathbb{F}_q`. Then we have

        .. MATH::

            \chi(t^n)(\phi(T))
            = t^{nr} + \sum_{i=1}^{r} \phi_{A_{i}}t^{n(i)}
            = 0,

        with `\deg(A_i) \leq \frac{n(r-i)}{r}`.

        Note that the *Frobenius trace* is defined as `A_{r-1}(T)` and the
        *Frobenius norm* is defined as `A_0(T)`.

        INPUT:

        - ``var`` -- (default: ``'X'``) the name of the second variable
        - ``algorithm`` -- (default: ``'crystalline'``) the algorithm
          used to compute the characteristic polynomial

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.frobenius_charpoly()
            X^2 + ((4*z2 + 4)*T^3 + (z2 + 3)*T^2 + 3*T + 2*z2 + 3)*X + 3*z2*T^6 + (4*z2 + 3)*T^5 + (4*z2 + 4)*T^4 + 2*T^3 + (3*z2 + 3)*T^2 + (z2 + 2)*T + 4*z2

        ::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: chi = phi.frobenius_charpoly()
            sage: chi
            X^2 + ((3*z3^2 + z3 + 4)*T + 4*z3^2 + 6*z3 + 3)*X + (5*z3^2 + 2*z3)*T^2 + (4*z3^2 + 3*z3)*T + 5*z3^2 + 2*z3

        ::

            sage: frob_pol = phi.frobenius_endomorphism().ore_polynomial()
            sage: chi(frob_pol, phi(T))
            0

        ::

            sage: phi.frobenius_charpoly(algorithm="NotImplemented")
            Traceback (most recent call last):
            ...
            NotImplementedError: algorithm "NotImplemented" not implemented

        ALGORITHM:

        By default, this method uses the so-called *crystalline*
        algorithm which computes the characteristic polynomial of the
        Frobenius acting on the crystalline cohomology of the Drinfeld
        module. For further details, see [Ang1997]_.

        The available options for 'algorithm' are:

        - ``'crystalline'`` -- computes the characteristic polynomial of the
          Frobenius endomorphism on the crystalline cohomology of a Drinfeld
          module.

        - ``'motive'`` -- based on computing the characteristic polynomial of
          the Frobenius endomorphism on the motive of a Drinfeld module. This
          instantiates the Frobenius as a morphism object and calls its
          ``'characteristic_polynomial'`` method.
        """
        # Throw an error if the user asks for an unimplemented algorithm
        # even if the char poly has already been computed
        method_name = f'_frobenius_charpoly_{algorithm}'
        if hasattr(self, method_name):
            if self._frobenius_charpoly is not None:
                return self._frobenius_charpoly
            self._frobenius_charpoly = getattr(self, method_name)(var)
            return self._frobenius_charpoly
        raise NotImplementedError(f'algorithm \"{algorithm}\" not implemented')

    def _frobenius_charpoly_crystalline(self, var):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism using Crystalline cohomology.

        The algorithm works for Drinfeld `\mathbb{F}_q[T]`-modules of
        any rank.

        This method is private and should not be directly called.
        Instead, use :meth:`frobenius_charpoly` with the option
        `algorithm='crystalline'`.

        INPUT:

        - ``var`` -- the name of the second variable

        OUTPUT: a univariate polynomial with coefficients in the
                function ring

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.frobenius_charpoly(algorithm='crystalline') # indirect doctest
            X^2 + ((4*z2 + 4)*T^3 + (z2 + 3)*T^2 + 3*T + 2*z2 + 3)*X + 3*z2*T^6 + (4*z2 + 3)*T^5 + (4*z2 + 4)*T^4 + 2*T^3 + (3*z2 + 3)*T^2 + (z2 + 2)*T + 4*z2

        ::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(8)
            sage: phi = DrinfeldModule(A, [z, 4, 1, z, z+1, 2, z+2, 1, 1, 3, 1])
            sage: phi.frobenius_charpoly(algorithm='crystalline') # indirect doctest
            X^10 + X^9 + (3*T + z2 + 1)*X^8 + (4*T^2 + z2*T + 2*z2 + 1)*X^7 + ... + (4*z2 + 4)*T^4 + 4*z2*T^2 + (z2 + 2)*T + z2

        ::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(10)
            sage: phi = DrinfeldModule(A, [z, z^2 + z, 2, 1, z, z+1, 2, z+2, 0, 1, 1, z^2 + z])
            sage: phi.frobenius_charpoly(algorithm='crystalline') # indirect doctest
            X^11 + (z3^2 + 2*z3)*X^10 + ((z3 + 1)*T + z3)*X^9 + ((2*z3^2 + z3 + 2)*T^2 + ... + (2*z3^2 + 2*z3 + 2)*T + z3^2

        ALGORITHM:

        Compute the characteristic polynomial of the Frobenius endomorphism
        acting on the crystalline cohomology of a Drinfeld module, which
        is equal to that of the Frobenius endomorphism on the Drinfeld
        module. A recurrence on elements of the cohomology allows us to
        compute a matrix representation of the Frobenius endomorphism
        efficiently using a companion matrix method. Based on the algorithm
        of section 6.3 in [MS2023]_.
        """
        A = self.function_ring()
        K = self.base_over_constants_field()
        charpoly_K = self._frobenius_crystalline_matrix() \
                         .charpoly(var).coefficients(sparse=False)

        # The above line obtains the char poly with coefficients in K[T]
        # This maps them into A = Fq[T]

        coeffs_A = [A([x.in_base() for x in coeff]) for coeff in charpoly_K]
        return PolynomialRing(A, name=var)(coeffs_A)

    def _frobenius_charpoly_motive(self, var):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism using Motivic cohomology.

        The algorithm works for Drinfeld `\mathbb{F}_q[T]`-modules of
        any rank.

        This method is private and should not be directly called.
        Instead, use :meth:`frobenius_charpoly` with the option
        `algorithm='motive'`.

        INPUT:

        - ``var`` -- the name of the second variable

        OUTPUT: a univariate polynomial with coefficients in the
                function ring

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.frobenius_charpoly(algorithm='motive') # indirect doctest
            X^2 + ((4*z2 + 4)*T^3 + (z2 + 3)*T^2 + 3*T + 2*z2 + 3)*X + 3*z2*T^6 + (4*z2 + 3)*T^5 + (4*z2 + 4)*T^4 + 2*T^3 + (3*z2 + 3)*T^2 + (z2 + 2)*T + 4*z2

        ::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(8)
            sage: phi = DrinfeldModule(A, [z, 4, 1, z, z+1, 2, z+2, 1, 1, 3, 1])
            sage: phi.frobenius_charpoly(algorithm='motive') # indirect doctest
            X^10 + X^9 + (3*T + z2 + 1)*X^8 + (4*T^2 + z2*T + 2*z2 + 1)*X^7 + ... + (4*z2 + 4)*T^4 + 4*z2*T^2 + (z2 + 2)*T + z2

        ::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(10)
            sage: phi = DrinfeldModule(A, [z, z^2 + z, 2, 1, z, z+1, 2, z+2, 0, 1, 1, z^2 + z])
            sage: phi.frobenius_charpoly(algorithm='motive') # indirect doctest
            X^11 + (z3^2 + 2*z3)*X^10 + ((z3 + 1)*T + z3)*X^9 + ((2*z3^2 + z3 + 2)*T^2 + ... + (2*z3^2 + 2*z3 + 2)*T + z3^2
        """
        return self.frobenius_endomorphism().characteristic_polynomial(var)

    def frobenius_norm(self):
        r"""
        Return the Frobenius norm of the Drinfeld module.

        Let `C(X) = \sum_{i=0}^r a_iX^{i}` denote the characteristic
        polynomial of the Frobenius endomorphism. The Frobenius norm
        is `(-1)^r a_{0}`. This is an element of the regular function ring
        and if `n` is the degree of the base field over `\mathbb{F}_q`,
        then the Frobenius norm has degree `n`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: B = phi.frobenius_norm()
            sage: B
            (5*z3^2 + 2*z3)*T^2 + (4*z3^2 + 3*z3)*T + 5*z3^2 + 2*z3

        ::

            sage: n = 2  # Degree of the base field over Fq
            sage: B.degree() == n
            True

        ::

            sage: B == phi.frobenius_charpoly()[0]
            True

        ALGORITHM:

            The Frobenius norm is computed using the formula, by
            Gekeler, given in [MS2019]_, Section 4, Proposition 3.
        """
        if self._frobenius_norm is not None:
            return self._frobenius_norm
        K = self.base_over_constants_field()
        n = K.degree(self._Fq)
        char = self.characteristic()
        norm = K(self.coefficients()[-1]).norm()
        self._frobenius_norm = ((-1)**n)*(char**(n/char.degree())) / norm
        return self._frobenius_norm

    def frobenius_trace(self):
        r"""
        Return the Frobenius trace of the Drinfeld module.

        Let `C(X) = \sum_{i=0}^r a_iX^{i}` denote the characteristic
        polynomial of the Frobenius endomorphism. The Frobenius trace
        is `-a_{r-1}`. This is an element of the regular function ring
        and if `n` is the degree of the base field over `\mathbb{F}_q`,
        then the Frobenius trace has degree at most `\frac{n}{r}`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: A = phi.frobenius_trace()
            sage: A
            (4*z3^2 + 6*z3 + 3)*T + 3*z3^2 + z3 + 4

        ::

            sage: n = 2  # Degree over Fq of the base codomain
            sage: A.degree() <= n/2
            True

        ::

            sage: A == -phi.frobenius_charpoly()[1]
            True

        ALGORITHM:

            We extract the coefficient of `X^{r-1}` from the characteristic
            polynomial if it has been previously computed, otherwise we compute
            the trace of the matrix of the Frobenius acting on the crystalline
            cohomology.
        """
        K = self.base_over_constants_field()
        A = self.function_ring()
        if self._frobenius_trace is not None:
            return self._frobenius_trace
        if self._frobenius_charpoly is not None:
            self._frobenius_trace = -self._frobenius_charpoly \
                                    .coefficients(sparse=False)[-2]
        self._frobenius_trace = self._frobenius_crystalline_matrix().trace()
        self._frobenius_trace = A([x.in_base() for x in self._frobenius_trace])
        return self._frobenius_trace

    def invert(self, ore_pol):
        r"""
        Return the preimage of the input under the Drinfeld module, if it
        exists.

        INPUT:

        - ``ore_pol`` -- the Ore polynomial whose preimage we want to
          compute

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: a = A.random_element()
            sage: phi.invert(phi(a)) == a
            True
            sage: phi.invert(phi(T)) == T
            True
            sage: phi.invert(phi(Fq.gen())) == Fq.gen()
            True

        When the input is not in the image of the Drinfeld module, an
        exception is raised::

            sage: t = phi.ore_polring().gen()
            sage: phi.invert(t + 1)
            Traceback (most recent call last):
            ...
            ValueError: input must be in the image of the Drinfeld module

        ::

            sage: phi.invert(t^4 + t^2 + 1)
            Traceback (most recent call last):
            ...
            ValueError: input must be in the image of the Drinfeld module

        ALGORITHM:

            The algorithm relies on the inversion of a linear algebra
            system. See [MS2019]_, 3.2.5 for details.

        TESTS::

            sage: a = A.random_element()
            sage: cat = phi.category()
            sage: phi_r1 = cat.random_object(1)
            sage: phi_r1.invert(phi_r1(a)) == a
            True
            sage: phi_r2 = cat.random_object(2)
            sage: phi_r2.invert(phi_r2(a)) == a
            True
            sage: phi_r3 = cat.random_object(3)
            sage: phi_r3.invert(phi_r3(a)) == a
            True
            sage: phi_r4 = cat.random_object(4)
            sage: phi_r4.invert(phi_r4(a)) == a
            True
            sage: phi_r5 = cat.random_object(5)
            sage: phi_r5.invert(phi_r5(a)) == a
            True

        ::

            sage: B.<X> = Fq[]
            sage: phi_r5.invert(X)
            Traceback (most recent call last):
            ...
            TypeError: input must be an Ore polynomial
        """
        deg = ore_pol.degree()
        r = self.rank()
        E = self.base()
        K = self.base_over_constants_field()
        if ore_pol in self._ore_polring:
            ore_pol = self._ore_polring(ore_pol)
        else:
            raise TypeError('input must be an Ore polynomial')

        # Trivial cases
        if deg <= 0:
            return K(ore_pol[0]).in_base()
        if deg % r != 0:
            raise ValueError('input must be in the image of the Drinfeld module')
        # Write the system and solve it
        k = deg // r
        A = self._function_ring
        T = A.gen()
        mat_rows = [[E.zero() for _ in range(k+1)] for _ in range(k+1)]
        mat_rows[0][0] = E.one()
        phiT = self.gen()
        phiTi = self.ore_polring().one()
        for i in range(1, k+1):
            phiTi *= phiT
            for j in range(i+1):
                mat_rows[j][i] = phiTi[r*j]
        mat = Matrix(mat_rows)
        vec = vector([ore_pol[r*j] for j in range(k+1)])
        coeffs_K = list(mat.inverse() * vec)
        # Cast the coefficients to Fq
        try:
            coeffs_Fq = [K(c).in_base() for c in coeffs_K]
            pre_image = A(coeffs_Fq)
            if self(pre_image) == ore_pol:
                return pre_image
        except ValueError:
            pass
        raise ValueError('input must be in the image of the Drinfeld module')

    def is_isogenous(self, psi):
        r"""
        Return ``True`` when ``self`` is isogenous to the other
        Drinfeld module.

        If the Drinfeld modules do not belong to the same category, an
        exception is raised.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: phi.is_isogenous(psi)
            True

        ::

            sage: chi = DrinfeldModule(A, [z, z + 1, z^2 + z])
            sage: phi.is_isogenous(chi)
            False

        ::

            sage: mu = DrinfeldModule(A, [z + 1, z^2 + z + 1, z^2 + z])
            sage: phi.is_isogenous(mu)
            Traceback (most recent call last):
            ...
            TypeError: Drinfeld modules are not in the same category

        ::

            sage: mu = 1
            sage: phi.is_isogenous(mu)
            Traceback (most recent call last):
            ...
            TypeError: input must be a Drinfeld module

        ALGORITHM:

        Two Drinfeld A-modules of equal characteristic are isogenous
        if and only if:

        - they have the same rank
        - the characteristic polynomial of the Frobenius endomorphism
          for both Drinfeld modules are equal.
        """
        if not isinstance(psi, DrinfeldModule):
            raise TypeError("input must be a Drinfeld module")
        if self.category() != psi.category():
            raise TypeError("Drinfeld modules are not in the same category")
        return self.rank() == psi.rank() \
               and self.frobenius_charpoly() == psi.frobenius_charpoly()

    def is_supersingular(self):
        r"""
        Return ``True`` if this Drinfeld module is supersingular.

        A Drinfeld module is supersingular if and only if its
        height equals its rank.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.is_supersingular()
            True
            sage: phi(phi.characteristic())   # Purely inseparable
            z6*t^2

        In rank two, a Drinfeld module is either ordinary or
        supersinguler. In higher ranks, it could be neither of
        the two::

            sage: psi = DrinfeldModule(A, [1, 0, z6, z6])
            sage: psi.is_ordinary()
            False
            sage: psi.is_supersingular()
            False
        """
        return self.height() == self.rank()

    def is_ordinary(self):
        r"""
        Return ``True`` if this Drinfeld module is ordinary.

        A Drinfeld module is ordinary if and only if its
        height is one.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.is_ordinary()
            False

        ::

            sage: phi = DrinfeldModule(A, [1, z6, 0, z6])
            sage: phi.is_ordinary()
            True
        """
        return self.height() == 1
