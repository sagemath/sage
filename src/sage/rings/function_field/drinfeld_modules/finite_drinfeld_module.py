r"""
Finite Drinfeld modules

This module provides the class
:class:`sage.rings.function_fields.drinfeld_module.finite_drinfeld_module.FiniteDrinfeldModule`,
which inherits
:class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

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

from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule


class FiniteDrinfeldModule(DrinfeldModule):
    r"""
    This class implements finite Drinfeld `\mathbb{F}_q[T]`-modules.

    A *finite Drinfeld module* is a Drinfeld module whose base field is
    finite. In this case, the function field characteristic is a prime
    ideal.

    For general definitions and help on Drinfeld modules, see class
    :class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

    .. RUBRIC:: Construction:

    The user does not ever need to directly call
    ``FiniteDrinfeldModule`` --- the metaclass ``DrinfeldModule`` is
    responsible for instantiating ``DrinfeldModule`` or
    ``FiniteDrinfeldModule`` depending on the input::

        sage: Fq = GF(343)
        sage: A.<T> = Fq[]
        sage: K.<z6> = Fq.extension(2)
        sage: phi = DrinfeldModule(A, [z6, 0, 5])
        sage: phi
        Drinfeld module defined by T |--> 5*t^2 + z6

    ::

        sage: isinstance(phi, DrinfeldModule)
        True
        sage: from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
        sage: isinstance(phi, FiniteDrinfeldModule)
        True

    The user should never use ``FiniteDrinfeldModule`` to test if a
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
        Initialize `self`.

        Validity of the input is checked in `__classcall_private__`. The
        `__init__` just saves attributes.

        INPUT:

        - ``function_ring`` -- a univariate polynomial ring whose base
          is a finite field

        - ``gen`` -- the generator of the Drinfeld module as a list of
          coefficients or an Ore polynomial

        - ``name`` (default: `'t'`) -- the name of the Ore polynomial
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
        # added one to ensure that FiniteDrinfeldModule would always
        # have _frobenius_norm and _frobenius_trace attributes.
        super().__init__(gen, category)
        self._frobenius_norm = None
        self._frobenius_trace = None

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
        Fq = self._function_ring.base()
        deg = self._base.over(Fq).degree(Fq)
        return self._Hom_(self, category=self.category())(t**deg)

    def frobenius_charpoly(self, var='X'):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism if the rank is two. Raise a NotImplementedError
        otherwise.

        Let `\mathbb{F}_q` be the base field of the function ring. The
        *characteristic polynomial `\chi` of the Frobenius endomorphism*
        is defined in [Gek1991]_. An important feature of this
        polynomial is that it is a monic univariate polynomial with
        coefficients in the function ring. As in our case the function
        ring is a univariate polynomial ring, it is customary to see the
        characteristic polynomial of the Frobenius endomorphism as a
        bivariate polynomial.

        Let `\chi = X^2 - A(T)X + B(T)` be the characteristic polynomial
        of the Frobenius endomorphism, and let `t^n` be the Ore polynomial
        that defines the Frobenius endomorphism of `\phi`; by
        definition, `n` is the degree over of the base field over
        `\mathbb{F}_q`. We have `\chi(t^n)(\phi(T)) = t^{2n} - \phi_A
        t^n + \phi_B = 0`, with `\deg(A) \leq \frac{n}{2}` and `\deg(B)
        = n`.

        Note that the *Frobenius trace* is defined as `A(T)` and the
        *Frobenius norm* is defined as `B(T)`.

        INPUT:

        - ``var`` (default: ``'X'``) -- the name of the second variable

        EXAMPLES::

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

            sage: trace = phi.frobenius_trace()
            sage: trace
            (4*z3^2 + 6*z3 + 3)*T + 3*z3^2 + z3 + 4
            sage: norm = phi.frobenius_norm()
            sage: norm
            (5*z3^2 + 2*z3)*T^2 + (4*z3^2 + 3*z3)*T + 5*z3^2 + 2*z3

        ::

            sage: n = 2  # Degree of the base field over Fq
            sage: trace.degree() <= n/2
            True
            sage: norm.degree() == n
            True

        ALGORITHM:

            We compute the Frobenius norm, and with it the Frobenius
            trace. This gives the Frobenius characteristic polynomial.
            See [MS2019]_, Section 4.

            See docstrings of methods :meth:`frobenius_norm` and
            :meth:`frobenius_trace` for further details on the
            computation of the norm and of the trace.
        """
        self._check_rank_two()
        A = self._function_ring  # Fq[T]
        S = PolynomialRing(A, name=var)  # Fq[T][X]
        # Does not work when Fq is not a prime field...
        # chi = self._gen.reduced_charpoly()
        # return -chi(A.gen(), S.gen())
        return S([self.frobenius_norm(), -self.frobenius_trace(), 1])

    def frobenius_norm(self):
        r"""
        Return Frobenius norm of the Drinfeld module, if the rank is
        two, raise a NotImplementedError otherwise.

        Let `\mathbb{F}_q[T]` be the function ring, write `\chi = X^2 -
        A(T)X + B(T) \in \mathbb{F}_q[T][X]` for the characteristic
        polynomial of the Frobenius endomorphism. The *Frobenius norm*
        is defined as the polynomial `B(T) \in \mathbb{F}_q[T]`.

        Let `n` be the degree of the base field over `\mathbb{F}_q` Then the
        Frobenius norm has degree `n`.

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
        self._check_rank_two()
        L = self._base.over(self._Fq)
        # Notations from Schost-Musleh:
        if self._frobenius_norm is None:
            n = L.degree_over(self._Fq)
            d = self.characteristic().degree()
            m = n // d
            delta = self._gen[2]
            norm = L(delta).norm()
            char = self.characteristic()
            self._frobenius_norm = ((-1)**n) * (char**m) / norm
        return self._frobenius_norm

    def frobenius_trace(self):
        r"""
        Return Frobenius norm of the Drinfeld module, if the rank is
        two; raise a NotImplementedError otherwise.

        Let `\mathbb{F}_q[T]` be the function ring, write `\chi = T^2 -
        A(X)T + B(X) \in \mathbb{F}_q[T][X]` for the characteristic
        polynomial of the Frobenius endomorphism. The *Frobenius trace*
        is defined as the polynomial `A(T) \in \mathbb{F}_q[T]`.

        Let `n` be the degree over `\mathbb{F}_q` of the base codomain.
        Then the Frobenius trace has degree at most `\frac{n}{2}`.

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

            Let `A(T)` denote the Frobenius trace and `B(T)` denote the
            Frobenius norm. We begin by computing `B(T)`, see docstring
            of method :meth:`frobenius_norm` for details. The
            characteristic polynomial of the Frobenius yields `t^{2n} -
            \phi_A t^n + \phi_B = 0`, where `t^n` is the Frobenius
            endomorphism. As `\phi_B` is now known, we can compute
            `\phi_A = (t^{2n} + \phi_B) / t^n`. We get `A(T)` by
            inverting this quantity, using the method
            :meth:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule.invert`,
            see its docstring for details.
        """
        self._check_rank_two()
        # Notations from Schost-Musleh:
        if self._frobenius_trace is None:
            n = self._base.over(self._Fq).degree_over(self._Fq)
            B = self.frobenius_norm()
            t = self.ore_polring().gen()
            self._frobenius_trace = self.invert(t**n + (self(B) // t**n))
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

            sage: phi.invert(t^3 + t^2 + 1)
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
        """
        deg = ore_pol.degree()
        r = self.rank()
        base_over_Fq = self.base_over_constants_field()
        if ore_pol not in self._ore_polring:
            raise TypeError('input must be an Ore polynomial')
        if ore_pol.degree() == 0:
            coord = base_over_Fq(self._base(ore_pol)).vector()
            if coord.nonzero_positions == [0]:
                return self._Fq(coord[0])
        if ore_pol == 0:
            return self._Fq.zero()
        if deg % r != 0:
            raise ValueError('input must be in the image of the Drinfeld '
                             'module')

        k = deg // r
        T = self._function_ring.gen()
        mat_lines = [[0 for _ in range(k+1)] for _ in range(k+1)]
        for i in range(k+1):
            phi_T_i = self(T**i)
            for j in range(i+1):
                mat_lines[j][i] = phi_T_i[r*j]
        mat = Matrix(mat_lines)
        vec = vector([list(ore_pol)[r*j] for j in range(k+1)])
        coeffs_K = list((mat**(-1)) * vec)
        coeffs_Fq = list(map(lambda x: base_over_Fq(x).vector()[0], coeffs_K))
        pre_image = self._function_ring(coeffs_Fq)

        if self(pre_image) == ore_pol:
            return pre_image
        else:
            raise ValueError('input must be in the image of the Drinfeld '
                             'module')

    def is_ordinary(self):
        r"""
        Return ``True`` whether the Drinfeld module is ordinary; raise a
        NotImplementedError if the rank is not two.

        A rank two finite Drinfeld module is *ordinary* if and only if
        the function ring-characteristic does not divide the Frobenius
        trace. A *supersingular* rank two finite Drinfeld module is a
        Drinfeld module that is not ordinary.

        A rank two Drinfeld module is *ordinary* if and only if it is
        not supersingular; see :meth:`is_supersingular`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.is_ordinary()
            False
            sage: phi_p = phi(phi.characteristic())
            sage: phi_p  # Purely inseparable
            z6*t^2

        ALGORITHM:

            Compute the Frobenius trace and test if the
            `\mathbb{F}_q[T]` characteristic divides it.

            We could also test if the image of the
            `\mathbb{F}_q[T]`-characteristic under the Drinfeld module
            is purely inseparable; see [Gek1991]_, Proposition 4.1.
        """
        self._check_rank_two()
        return not self.is_supersingular()

    def is_supersingular(self):
        r"""
        Return ``True`` whether the Drinfeld module is supersingular; raise a
        NotImplementedError if the rank is not two.

        A rank two finite Drinfeld module is *supersingular* if and only
        if the function ring-characteristic divides the Frobenius
        trace. An *ordinary* rank two finite Drinfeld module is a
        Drinfeld module that is not supersingular.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.is_supersingular()
            True
            sage: phi(phi.characteristic()) # Purely inseparable
            z6*t^2

        ALGORITHM:

            Compute the Frobenius trace and test if the function
            ring-characteristic divides it.

            We could also test if the image of the function
            ring-characteristic under the Drinfeld module is purely
            inseparable; see [Gek1991]_, Proposition 4.1.
        """
        self._check_rank_two()
        return self.characteristic().divides(self.frobenius_trace())
