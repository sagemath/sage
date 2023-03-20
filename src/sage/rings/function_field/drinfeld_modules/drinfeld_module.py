r"""
Drinfeld modules

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.drinfeld_module.DrinfeldModule`.

For finite Drinfeld modules and their theory of complex multiplication, see
class
:class:`sage.rings.function_field.drinfeld_module.finite_drinfeld_module.DrinfeldModule`.

AUTHORS:

- Antoine Leudière (2022-04)
- Xavier Caruso (2022-06)
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

from sage.categories.drinfeld_modules import DrinfeldModules
from sage.categories.homset import Hom
from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name
from sage.misc.lazy_string import _LazyString
from sage.rings.integer import Integer
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.ring_extension import RingExtension_generic
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence
from sage.structure.unique_representation import UniqueRepresentation


class DrinfeldModule(Parent, UniqueRepresentation):
    r"""
    This class implements Drinfeld `\mathbb{F}_q[T]`-modules.

    Let `\mathbb{F}_q[T]` be a polynomial ring with coefficients in a
    finite field `\mathbb{F}_q` and let `K` be a field. Fix a ring
    morphism `\gamma: \mathbb{F}_q[T] \to K`; we say that `K` is an
    `\mathbb{F}_q[T]`-*field*. Let `K\{\tau\}` be the ring of Ore
    polynomials with coefficients in `K`, whose multiplication is given
    by the rule `\tau \lambda = \lambda^q \tau` for any `\lambda \in K`.

    A Drinfeld `\mathbb{F}_q[T]`-module over the base
    `\mathbb{F}_q[T]`-field `K` is an `\mathbb{F}_q`-algebra morphism
    `\phi: \mathbb{F}_q[T] \to K\{\tau\}` such that `\mathrm{Im}(\phi)
    \not\subset K` and `\phi` agrees with `\gamma` on `\mathbb{F}_q`.

    For `a` in `\mathbb{F}_q[T]`, `\phi(a)` is denoted `\phi_a`.

    The Drinfeld `\mathbb{F}_q[T]`-module `\phi` is uniquely determined
    by the image `\phi_T` of `T`; this serves as input of the class.

    .. NOTE::

        See also :class:`sage.categories.drinfeld_modules`.

    The *base morphism* is the morphism `\gamma: \mathbb{F}_q[T] \to K`.
    The monic polynomial that generates the kernel of `\gamma` is called
    the `\mathbb{F}_q[T]`-*characteristic*, or *function-field
    characteristic*, of the base field. We say that `\mathbb{F}_q[T]` is
    the *function ring* of `\phi`; `K\{\tau\}` is the *Ore polynomial
    ring*. Further, the *generator* is `\phi_T` and the *constant
    coefficient* is the constant coefficient of `\phi_T`.

    A Drinfeld module is said to be *finite* if the field `K` is.
    Despite an emphasis on this case, the base field can be any
    extension of `\mathbb{F}_q`::

        sage: Fq = GF(25)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(6)
        sage: phi = DrinfeldModule(A, [z, 4, 1])
        sage: phi
        Drinfeld module defined by T |--> t^2 + 4*t + z

    ::

        sage: Fq = GF(49)
        sage: A.<T> = Fq[]
        sage: K = Frac(A)
        sage: psi = DrinfeldModule(A, [K(T), T+1])
        sage: psi
        Drinfeld module defined by T |--> (T + 1)*t + T

    .. NOTE::

        Finite Drinfeld modules are implemented in the class
        :class:`sage.rings.function_field.drinfeld_modules.finite_drinfeld_module`.

    Classical references on Drinfeld modules include [Gos1998]_,
    [Rosen2002]_, [VS06]_ and [Gek1991]_.

    .. NOTE::

        Drinfeld modules are defined in a larger setting, in which the
        polynomial ring `\mathbb{F}_q[T]` is replaced by a more general
        function ring: the ring of functions in `k` that are regular
        outside `\infty`, where `k` is a function field over
        `\mathbb{F}_q` with transcendence degree `1` and `\infty` is a
        fixed place of `k`. This is out of the scope of this
        implementation.

    INPUT:

    - ``function_ring`` -- a univariate polynomial ring whose base field
      is a finite field

    - ``gen`` -- the generator of the Drinfeld module; as a list of
      coefficients or an Ore polynomial

    - ``name`` (default: ``'t'``) -- the name of the Ore polynomial ring
      generator

    .. RUBRIC:: Construction

    A Drinfeld module object is constructed by giving the function ring
    and the generator::

        sage: Fq.<z2> = GF(3^2)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(6)
        sage: phi = DrinfeldModule(A, [z, 1, 1])
        sage: phi
        Drinfeld module defined by T |--> t^2 + t + z

    .. NOTE::

        Note that the definition of the base field is implicit; it is
        automatically defined as the compositum of all the parents of
        the coefficients.

    The above Drinfeld module is finite; it can also be infinite::

        sage: L = Frac(A)
        sage: psi = DrinfeldModule(A, [L(T), 1, T^3 + T + 1])
        sage: psi
        Drinfeld module defined by T |--> (T^3 + T + 1)*t^2 + t + T

    ::

        sage: phi.is_finite()
        True
        sage: psi.is_finite()
        False

    In those examples, we used a list of coefficients (``[z, 1, 1]``) to
    represent the generator `\phi_T = z + t + t^2`. One can also use
    regular Ore polynomials::

        sage: ore_polring = phi.ore_polring()
        sage: t = ore_polring.gen()
        sage: rho_T = z + t^3
        sage: rho = DrinfeldModule(A, rho_T)
        sage: rho
        Drinfeld module defined by T |--> t^3 + z
        sage: rho(T) == rho_T
        True

    Images under the Drinfeld module are computed by calling the object::

        sage: phi(T)  # phi_T, the generator of the Drinfeld module
        t^2 + t + z
        sage: phi(T^3 + T + 1)  # phi_(T^3 + T + 1)
        t^6 + (z^11 + z^9 + 2*z^6 + 2*z^4 + 2*z + 1)*t^4 + (2*z^11 + 2*z^10 + z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^3)*t^3 + (2*z^11 + z^10 + z^9 + 2*z^7 + 2*z^6 + z^5 + z^4 + 2*z^3 + 2*z + 2)*t^2 + (2*z^11 + 2*z^8 + 2*z^6 + z^5 + z^4 + 2*z^2)*t + z^3 + z + 1
        sage: phi(1)  # phi_1
        1

    .. RUBRIC:: The category of Drinfeld modules

    Drinfeld modules have their own category (see class
    :class:`sage.categories.drinfeld_modules.DrinfeldModules`)::

        sage: phi.category()
        Category of Drinfeld modules over Finite Field in z of size 3^12 over its base
        sage: phi.category() is psi.category()
        False
        sage: phi.category() is rho.category()
        True

    One can use the category to directly create new objects::

        sage: cat = phi.category()
        sage: cat.object([z, 0, 0, 1])
        Drinfeld module defined by T |--> t^3 + z

    .. RUBRIC:: The base field of a Drinfeld module

    The base field of the Drinfeld module is retrieved using :meth:`base`::

        sage: phi.base()
        Finite Field in z of size 3^12 over its base

    The base morphism is retrieved using :meth:`base_morphism`::

        sage: phi.base_morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in T over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12 over its base
          Defn: T |--> z

    Note that the base field is *not* the field `K`. Rather, it is a ring
    extension (see :class:`sage.rings.ring_extension.RingExtension`) whose
    underlying ring is `K` and whose base is the base morphism::

        sage: phi.base() is K
        False

    .. RUBRIC:: Getters

    One can retrieve basic properties::

        sage: phi.base_morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in T over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12 over its base
          Defn: T |--> z

    ::

        sage: phi.ore_polring()  # K{t}
        Ore Polynomial Ring in t over Finite Field in z of size 3^12 over its base twisted by Frob^2

    ::

        sage: phi.function_ring()  # Fq[T]
        Univariate Polynomial Ring in T over Finite Field in z2 of size 3^2

    ::

        sage: phi.gen()  # phi_T
        t^2 + t + z
        sage: phi.gen() == phi(T)
        True

    ::

        sage: phi.constant_coefficient()  # Constant coefficient of phi_T
        z

    ::

        sage: phi.morphism()  # The Drinfeld module as a morphism
        Ring morphism:
          From: Univariate Polynomial Ring in T over Finite Field in z2 of size 3^2
          To:   Ore Polynomial Ring in t over Finite Field in z of size 3^12 over its base twisted by Frob^2
          Defn: T |--> t^2 + t + z

    One can compute the rank and height::

        sage: phi.rank()
        2
        sage: phi.height()
        1

    As well as the j-invariant if the rank is two::

        sage: phi.j_invariant()  # j-invariant
        1

    A Drinfeld `\mathbb{F}_q[T]`-module can be seen as an Ore polynomial
    with positive degree and constant coefficient `\gamma(T)`, where
    `\gamma` is the base morphism. This analogy is the motivation for
    the following methods::

        sage: phi.coefficients()
        [z, 1, 1]

    ::

        sage: phi.coefficient(1)
        1


    .. RUBRIC:: Morphisms and isogenies

    A *morphism* of Drinfeld modules `\phi \to \psi` is an Ore
    polynomial `f \in K\{\tau\}` such that `f \phi_a = \psi_a f` for
    every `a` in the function ring. In our case, this is equivalent to
    `f \phi_T = \psi_T f`. An *isogeny* is a nonzero morphism.

    Use the ``in`` syntax to test if an Ore polynomial defines a
    morphism::

        sage: phi(T) in Hom(phi, phi)
        True
        sage: t^6 in Hom(phi, phi)
        True
        sage: t^5 + 2*t^3 + 1 in Hom(phi, phi)
        False
        sage: 1 in Hom(phi, rho)
        False
        sage: 1 in Hom(phi, phi)
        True
        sage: 0 in Hom(phi, rho)
        True

    To create a SageMath object representing the morphism, call the
    homset (``hom``)::

        sage: hom = Hom(phi, phi)
        sage: frobenius_endomorphism = hom(t^6)
        sage: identity_morphism = hom(1)
        sage: zero_morphism = hom(0)
        sage: frobenius_endomorphism
        Endomorphism of Drinfeld module defined by T |--> t^2 + t + z
          Defn: t^6
        sage: identity_morphism
        Identity morphism of Drinfeld module defined by T |--> t^2 + t + z
        sage: zero_morphism
        Endomorphism of Drinfeld module defined by T |--> t^2 + t + z
          Defn: 0

    The underlying Ore polynomial is retrieved with the method
    :meth:`ore_polynomial`::

        sage: frobenius_endomorphism.ore_polynomial()
        t^6
        sage: identity_morphism.ore_polynomial()
        1

    One checks if a morphism is an isogeny, endomorphism or
    isomorphism::

        sage: frobenius_endomorphism.is_isogeny()
        True
        sage: identity_morphism.is_isogeny()
        True
        sage: zero_morphism.is_isogeny()
        False
        sage: frobenius_endomorphism.is_isomorphism()
        False
        sage: identity_morphism.is_isomorphism()
        True
        sage: zero_morphism.is_isomorphism()
        False

    .. RUBRIC:: The Vélu formula

    Let ``P`` be a nonzero Ore polynomial. We can decide if ``P``
    defines an isogeny with a given domain and, if it does, find
    the codomain::

        sage: P = (2*z^6 + z^3 + 2*z^2 + z + 2)*t + z^11 + 2*z^10 + 2*z^9 + 2*z^8 + z^7 + 2*z^6 + z^5 + z^3 + z^2 + z
        sage: psi = phi.velu(P)
        sage: psi
        Drinfeld module defined by T |--> (2*z^11 + 2*z^9 + z^6 + 2*z^5 + 2*z^4 + 2*z^2 + 1)*t^2 + (2*z^11 + 2*z^10 + 2*z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z)*t + z
        sage: P in Hom(phi, psi)
        True
        sage: P * phi(T) == psi(T) * P
        True

    If the input does not define an isogeny, an exception is raised::

        sage: phi.velu(0)
        Traceback (most recent call last):
        ...
        ValueError: the input does not define an isogeny
        sage: phi.velu(t)
        Traceback (most recent call last):
        ...
        ValueError: the input does not define an isogeny

    .. RUBRIC:: The action of a Drinfeld module

    The `\mathbb{F}_q[T]`-Drinfeld module `\phi` induces a special left
    `\mathbb{F}_q[T]`-module structure on any field extension `L/K`. Let
    `x \in L` and `a` be in the function ring; the action is defined as
    `(a, x) \mapsto \phi_a(x)`. The method :meth:`action` returns a
    :class:`sage.rings.function_field.drinfeld_modules.action.Action`
    object representing the Drinfeld module action.

    .. NOTE::

        In this implementation, `L` is `K`::

            sage: action = phi.action()
            sage: action
            Action on Finite Field in z of size 3^12 over its base induced by Drinfeld module defined by T |--> t^2 + t + z

    The action on elements is computed by calling the action object::

        sage: P = T + 1
        sage: a = z
        sage: action(P, a)
        ...
        z^9 + 2*z^8 + 2*z^7 + 2*z^6 + 2*z^3 + z^2
        sage: action(0, K.random_element())
        0
        sage: action(A.random_element(), 0)
        0

    .. WARNING::

        The class ``DrinfeldModuleAction`` may be replaced later on. See
        issues #34833 and #34834.

    TESTS:

    The generator must have positive degree::

        sage: Fq = GF(2)
        sage: K.<z> = Fq.extension(2)
        sage: A.<T> = Fq[]
        sage: DrinfeldModule(A, [K(1)])
        Traceback (most recent call last):
        ...
        ValueError: generator must have positive degree

    The constant coefficient must be nonzero::

        sage: Fq = GF(2)
        sage: K.<z> = Fq.extension(2)
        sage: A.<T> = Fq[]
        sage: DrinfeldModule(A, [K(0), K(1)])
        Traceback (most recent call last):
        ...
        ValueError: constant coefficient must be nonzero

    The coefficients of the generator must lie in an
    `\mathbb{F}_q[T]`-field, where `\mathbb{F}_q[T]` is the function
    ring of the Drinfeld module::

        sage: Fq = GF(2)
        sage: K.<z> = Fq.extension(2)
        sage: A.<T> = Fq[]
        sage: DrinfeldModule(A, [z, QQ(1)])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce into base field

    ::

        sage: Fq = GF(2)
        sage: K.<z> = Fq.extension(2)
        sage: A.<T> = Fq[]
        sage: DrinfeldModule(A, [1, QQ(1)])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce into base field

    The function ring must be an univariate polynomial ring whose
    base is a finite field::

        sage: Fq = GF(2)
        sage: K.<z> = Fq.extension(2)
        sage: A.<T> = Fq[]
        sage: DrinfeldModule(K, [z, 1, 1])
        Traceback (most recent call last):
        ...
        NotImplementedError: function ring must be a polynomial ring

    ::

        sage: Fq = GF(2)
        sage: K.<z> = Fq.extension(2)
        sage: A.<T> = Fq[]
        sage: AY.<Y> = A[]
        sage: DrinfeldModule(AY, [z, 1, 1])
        Traceback (most recent call last):
        ...
        TypeError: function ring base must be a finite field

    If you already defined a category of Drinfeld modules, and you
    create a Drinfeld module through this category, you must ensure that
    the constant coefficient is that of the category::

        sage: Fq = GF(2)
        sage: K.<z> = Fq.extension(2)
        sage: A.<T> = Fq[]
        sage: phi = DrinfeldModule(A, [z, 1])
        sage: phi.category().object([1, 1, K(1)])
        Traceback (most recent call last):
        ...
        ValueError: constant coefficient must equal that of the category

    ::

        sage: Fq = K = GF(2)
        sage: A.<T> = Fq[]
        sage: phi = DrinfeldModule(A, [1, 1])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce into base field

    ::

        sage: Fq = K = GF(2)
        sage: A.<T> = Fq[]
        sage: phi = DrinfeldModule(A, [K(1), 1])
        sage: isinstance(phi.ore_polring(), OrePolynomialRing)
        True
    """

    @staticmethod
    def __classcall_private__(cls, function_ring, gen, name='t'):
        """
        Check input validity and return a ``DrinfeldModule`` or
        ``FiniteDrinfeldModule`` object accordingly.

        INPUT:

        - ``function_ring`` -- a univariate polynomial ring whose base
          is a finite field

        - ``gen`` -- the generator of the Drinfeld module; as a list of
          coefficients or an Ore polynomial

        - ``name`` (default: ``'t'``) -- the name of the Ore polynomial
          ring gen

        OUTPUT:

        A DrinfeldModule or FiniteDrinfeldModule.

        TESTS::

            sage: from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: isinstance(phi, FiniteDrinfeldModule)
            True

        ::

            sage: K = Frac(A)
            sage: phi = DrinfeldModule(A, [K(T), 1])
            sage: isinstance(psi, FiniteDrinfeldModule)
            False
        """

        # FIXME: function_ring must be checked before calling base_ring
        # on it. But then it is checked twice: firstly here, secondly in
        # the category. Another problem is that those lines are
        # duplicate. As a general comment, there are sanity checks both
        # here and in the category constructor, which is not ideal.
        # Check domain is Fq[T]
        if not isinstance(function_ring, PolynomialRing_general):
            raise NotImplementedError('function ring must be a polynomial '
                                      'ring')
        function_ring_base = function_ring.base_ring()
        if not function_ring_base.is_field() \
                or not function_ring_base.is_finite():
            raise TypeError('function ring base must be a finite field')

        # Check all possible input types for gen
        # `gen` is an Ore polynomial:
        if isinstance(gen, OrePolynomial):
            ore_polring = gen.parent()
            # Base ring without morphism structure:
            base_field_noext = ore_polring.base()
            name = ore_polring.variable_name()
        # `gen` is a list of coefficients (function_ring = Fq[T]):
        elif isinstance(gen, (list, tuple)):
            ore_polring = None
            # Base ring without morphism structure:
            base_field_noext = Sequence(gen).universe()
        else:
            raise TypeError('generator must be list of coefficients or Ore '
                            'polynomial')
        # Constant coefficient must be nonzero:
        if gen[0].is_zero():
            raise ValueError('constant coefficient must be nonzero')
        # The coefficients are in a base field that has coercion from Fq:
        if not (hasattr(base_field_noext, 'has_coerce_map_from') and
                base_field_noext.has_coerce_map_from(function_ring.base_ring())):
            raise ValueError('function ring base must coerce into base field')

        # Build the category
        T = function_ring.gen()
        if isinstance(base_field_noext, RingExtension_generic):
            base_field = base_field_noext
        elif base_field_noext.has_coerce_map_from(function_ring) \
                and T == gen[0]:
            base_morphism = base_field_noext.coerce_map_from(function_ring)
            base_field = base_field_noext.over(base_morphism)
        else:
            base_morphism = Hom(function_ring, base_field_noext)(gen[0])
            base_field = base_field_noext.over(base_morphism)

        # This test is also done in the category. We put it here also
        # to have a friendlier error message
        if not base_field.is_field():
            raise ValueError('generator coefficients must live in a field')

        category = DrinfeldModules(base_field, name=name)

        # Check gen as Ore polynomial
        ore_polring = category.ore_polring()  # Sanity cast
        gen = ore_polring(gen)
        if gen.degree() <= 0:
            raise ValueError('generator must have positive degree')

        # Instantiate the appropriate class
        if base_field.is_finite():
            from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
            return FiniteDrinfeldModule(gen, category)
        return cls.__classcall__(cls, gen, category)

    def __init__(self, gen, category):
        """
        Initialize ``self``.

        Validity of the input is checked in meth:`__classcall_private__`.
        The meth:`__init__` just saves attributes.

        INPUT:

        - ``function_ring`` -- a univariate polynomial ring whose base
          is a finite field

        - ``gen`` -- the generator of the Drinfeld module; as a list of
          coefficients or an Ore polynomial

        - ``name`` (default: ``'t'``) -- the name of the Ore polynomial
          ring gen

        TESTS::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: gen = [p_root, z12^3, z12^5]
            sage: phi = DrinfeldModule(A, gen)
            sage: ore_polring = phi.ore_polring()
            sage: phi._base == phi.category().base()
            True
            sage: phi._function_ring == A
            True
            sage: phi._gen == ore_polring(gen)
            True
            sage: phi._ore_polring == ore_polring
            True
            sage: phi._morphism == Hom(A, ore_polring)(phi._gen)
            True

        ::

            sage: TestSuite(phi).run()
        """
        self._base = category.base()
        self._function_ring = category.function_ring()
        self._gen = gen
        self._morphism = category._function_ring.hom([gen])
        self._ore_polring = gen.parent()
        self._Fq = self._function_ring.base_ring()  # Must be last
        super().__init__(base=self._base, category=category)

    def __call__(self, a):
        r"""
        Return the image of input ``a`` by the morphism that defines the
        Drinfeld module; i.e. `\phi_a` if the Drinfeld module is denoted
        `phi`.

        INPUT:

        - ``a`` -- a function ring element

        OUTPUT: an element in the base codomain

        TESTS::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])

        ::
            sage: a = T^3 + 4*T + 2
            sage: phi(a) == phi(T)^3 + 4*phi(T) + 2
            True
            sage: phi(a)[0] == p_root^3 + 4*p_root + 2
            True

        ::

            sage: phi(0)
            0
            sage: phi(1)
            1
            sage: phi(T) == phi._gen
            True

        ::

            sage: a = A.random_element(5)
            sage: phi(a)[0] == phi.category().base()(a)
            True
        """
        return self._morphism(a)

    def _Hom_(self, other, category):
        r"""
        Return ``DrinfeldModuleHomset(self, other, category)``.

        Validity of the input is checked at the instantiation of
        ``DrinfeldModuleHomset``. ``self`` and ``other`` only need be in
        the same category.

        INPUT:

        - ``other`` -- the codomain of the homset

        - ``category`` -- the category in which we consider the
          morphisms, usually ``self.category()``

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: t = phi.ore_polring().gen()
            sage: isog = t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4
            sage: psi = phi.velu(isog)
            sage: hom = phi._Hom_(psi, category=phi.category())
            sage: hom is Hom(phi, psi)  # known bug
            True
            sage: from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
            sage: isinstance(hom, DrinfeldModuleHomset)
            True
        """
        from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        return DrinfeldModuleHomset(self, other, category)

    def _check_rank_two(self):
        r"""
        Raise ``NotImplementedError`` if the rank is not two.

        TESTS::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi._check_rank_two()
            sage: phi = DrinfeldModule(A, [p_root, 1])
            sage: phi._check_rank_two()
            Traceback (most recent call last):
            ...
            NotImplementedError: rank must be 2
        """
        if self.rank() != 2:
            raise NotImplementedError('rank must be 2')

    def _latex_(self):
        r"""
        Return a LaTeX representation of the Drinfeld module.

        If a representation name is given with meth:`rename`, it is
        taken into account for LaTeX representation.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: latex(phi)
            \phi: T \mapsto z_{12}^{5} t^{2} + z_{12}^{3} t + 2 z_{12}^{11} + 2 z_{12}^{10} + z_{12}^{9} + 3 z_{12}^{8} + z_{12}^{7} + 2 z_{12}^{5} + 2 z_{12}^{4} + 3 z_{12}^{3} + z_{12}^{2} + 2 z_{12}

        ::

            sage: phi.rename('phi')
            sage: latex(phi)
            \phi
            sage: phi.reset_name()
        """
        if hasattr(self, '__custom_name'):
            return latex_variable_name(getattr(self, '__custom_name'))
        else:
            return f'\\phi: {latex(self._function_ring.gen())} \\mapsto ' \
                   f'{latex(self._gen)}'

    def _repr_(self):
        r"""
        Return a string representation of this Drinfeld module.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi
            Drinfeld module defined by T |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
        """
        return f'Drinfeld module defined by {self._function_ring.gen()} ' \
               f'|--> {self._gen}'

    def _test_category(self, **options):
        """
        Run generic tests on the method :meth:`.category`.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi._test_category()

        .. NOTE::

            We reimplemented this method because Drinfeld modules are
            parents, and
            meth:`sage.structure.parent.Parent._test_category` requires
            parents' categories to be subcategories of ``Sets()``.
        """
        tester = self._tester(**options)
        SageObject._test_category(self, tester=tester)
        category = self.category()
        # Tests that self inherits methods from the categories
        tester.assertTrue(isinstance(self, category.parent_class),
                _LazyString("category of %s improperly initialized", (self,), {}))

    def __hash__(self):
        r"""
        Return a hash of ``self``.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: hash(phi)  # random
            -6894299335185957188
        """
        return hash((self.base(), self._gen))

    def action(self):
        r"""
        Return the action object
        (:class:`sage.rings.function_field.drinfeld_modules.action.Action`)
        that represents the module action, on the base codomain, that is
        induced by the Drinfeld module.

        OUTPUT: a Drinfeld module action object

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: action = phi.action()
            sage: action
            Action on Finite Field in z12 of size 5^12 over its base induced by Drinfeld module defined by T |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12

        The action on elements is computed as follows::

            sage: P = T^2 + T + 1
            sage: a = z12 + 1
            sage: action(P, a)
            3*z12^11 + 2*z12^10 + 3*z12^9 + 3*z12^7 + 4*z12^5 + z12^4 + z12^3 + 2*z12 + 1
            sage: action(0, a)
            0
            sage: action(P, 0)
            0
        """
        from sage.rings.function_field.drinfeld_modules.action import DrinfeldModuleAction
        return DrinfeldModuleAction(self)

    def coefficient(self, n):
        r"""
        Return the `n`-th coefficient of the generator.

        INPUT:

        - ``n`` -- a nonnegative integer

        OUTPUT: an element in the base codomain

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.coefficient(0)
            2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi.coefficient(0) == p_root
            True
            sage: phi.coefficient(1)
            z12^3
            sage: phi.coefficient(2)
            z12^5
            sage: phi.coefficient(5)
            Traceback (most recent call last):
            ...
            ValueError: input must be >= 0 and <= rank
        """
        if not isinstance(n, Integer) and not isinstance(n, int):
            raise TypeError('input must be an integer')
        if not 0 <= n <= self.rank():
            raise ValueError('input must be >= 0 and <= rank')
        return self.coefficients(sparse=False)[n]

    def coefficients(self, sparse=True):
        r"""
        Return the coefficients of the generator, as a list.

        If the flag ``sparse`` is ``True`` (default), only return the
        nonzero coefficients; otherwise, return all of them.

        INPUT:

        - ``sparse`` -- a boolean

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.coefficients()
            [2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12,
             z12^3,
             z12^5]

        Careful, the method only returns the nonzero coefficients,
        unless otherwise specified::

            sage: rho = DrinfeldModule(A, [p_root, 0, 0, 0, 1])
            sage: rho.coefficients()
            [2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12,
             1]
            sage: rho.coefficients(sparse=False)
            [2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12,
             0,
             0,
             0,
             1]
        """
        return self._gen.coefficients(sparse=sparse)

    def gen(self):
        r"""
        Return the generator of the Drinfeld module.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.gen() == phi(T)
            True
        """
        return self._gen

    def height(self):
        r"""
        Return the height of the Drinfeld module if the function field
        characteristic is a prime ideal; raise ValueError otherwise.

        The height of a Drinfeld module is defined when the function
        field characteristic is a prime ideal. In our case, this ideal
        is even generated by a monic polynomial `\mathfrak{p}` in the
        function field. Write `\phi_\mathfrak{p} = a_s \tau^s + \dots +
        \tau^{r*\deg(\mathfrak{p})}`. The height of the Drinfeld module
        is the well-defined positive integer `h =
        \frac{s}{\deg(\mathfrak{p})}`.

        .. NOTE::

            See [Gos1998]_, Definition 4.5.8 for the general definition.

        A rank two Drinfeld module is supersingular if and only if its
        height equals its rank.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.height() == 1
            True
            sage: phi.is_ordinary()
            True

        ::

            sage: B.<Y> = Fq[]
            sage: L = Frac(B)
            sage: phi = DrinfeldModule(A, [L(2), L(1)])
            sage: phi.height()
            Traceback (most recent call last):
            ...
            NotImplementedError: height not implemented in this case

        ::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.height()
            2
            sage: phi.is_supersingular()
            True

        """
        try:
            if self.characteristic().is_zero():
                raise ValueError('height is defined for prime '
                                 'function field characteristic')
            else:
                p = self.characteristic()
                return Integer(self(p).valuation() // p.degree())
        except NotImplementedError:
            raise NotImplementedError('height not implemented in this case')

    def is_finite(self) -> bool:
        r"""
        Return ``True`` if this Drinfeld module is finite,
        ``False`` otherwise.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.is_finite()
            True
            sage: B.<Y> = Fq[]
            sage: L = Frac(B)
            sage: psi = DrinfeldModule(A, [L(2), L(1)])
            sage: psi.is_finite()
            False
        """
        from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
        return isinstance(self, FiniteDrinfeldModule)

    def j_invariant(self):
        r"""
        Return the j-invariant of the Drinfeld module if the rank is
        two; raise a NotImplementedError otherwise.

        Assume the rank is two. Write the generator `\phi_T = \omega +
        g\tau + \Delta\tau^2`. The j-invariant is defined by
        `\frac{g^{q+1}}{\Delta}`, `q` being the order of the base field
        of the function ring. In our case, this field is always finite.

        OUTPUT: an element in the base codomain

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.j_invariant()
            z12^10 + 4*z12^9 + 3*z12^8 + 2*z12^7 + 3*z12^6 + z12^5 + z12^3 + 4*z12^2 + z12 + 2
            sage: psi = DrinfeldModule(A, [p_root, 1, 1])
            sage: psi.j_invariant()
            1
            sage: rho = DrinfeldModule(A, [p_root, 0, 1])
            sage: rho.j_invariant()
            0

        The rank must be two::

            sage: sigma = DrinfeldModule(A, [p_root, 1, 0])
            sage: sigma.j_invariant()
            Traceback (most recent call last):
            ...
            NotImplementedError: rank must be 2
        """
        self._check_rank_two()
        g = self.coefficient(1)
        delta = self.coefficient(2)
        q = self._Fq.order()
        return (g**(q+1)) / delta

    def morphism(self):
        r"""
        Return the morphism object that defines the Drinfeld module.

        OUTPUT: a ring morphism from the function ring to the Ore
        polynomial ring

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.morphism()
            Ring morphism:
              From: Univariate Polynomial Ring in T over Finite Field in z2 of size 5^2
              To:   Ore Polynomial Ring in t over Finite Field in z12 of size 5^12 over its base twisted by Frob^2
              Defn: T |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: from sage.rings.morphism import RingHomomorphism
            sage: isinstance(phi.morphism(), RingHomomorphism)
            True

        Actually, the ``DrinfeldModule`` method :meth:`__call__` simply
        class the ``__call__`` method of this morphism::

            sage: phi.morphism()(T) == phi(T)
            True
            sage: a = A.random_element()
            sage: phi.morphism()(a) == phi(a)
            True

        And many methods of the Drinfeld module have a counterpart in
        the morphism object::

            sage: m = phi.morphism()
            sage: m.domain() is phi.function_ring()
            True
            sage: m.codomain() is phi.ore_polring()
            True
            sage: m.im_gens()
            [z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12]
            sage: phi(T) == m.im_gens()[0]
            True
        """
        return self._morphism

    def rank(self):
        r"""
        Return the rank of the Drinfeld module.

        In our case, the rank is the degree of the generator.

        OUTPUT: an integer

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.rank()
            2
            sage: psi = DrinfeldModule(A, [p_root, 2])
            sage: psi.rank()
            1
            sage: rho = DrinfeldModule(A, [p_root, 0, 0, 0, 1])
            sage: rho.rank()
            4
        """
        return self._gen.degree()

    def velu(self, isog):
        r"""
        Return a new Drinfeld module such that input is an
        isogeny to this module with domain ``self``; if no such isogeny
        exists, raise an exception.

        INPUT:

        - ``isog`` -- the Ore polynomial that defines the isogeny

        OUTPUT: a Drinfeld module

        ALGORITHM:

            The input defines an isogeny if only if:

            1. The degree of the characteristic divides the height of
            the input. (The height of an Ore polynomial `P(\tau)` is the
            maximum `n` such that `\tau^n` right-divides `P(\tau)`.)

            2. The input right-divides the generator, which can
            be tested with Euclidean division.

            We test if the input is an isogeny, and, if it is, we
            return the quotient of the Euclidean division.

            Height and Euclidean division of Ore polynomials are
            implemented as methods of class
            :class:`sage.rings.polynomial.ore_polynomial_element.OrePolynomial`.

            Another possible algorithm is to recursively solve a system,
            see :arxiv:`2203.06970`, Eq. 1.1.

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: t = phi.ore_polring().gen()
            sage: isog = t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4
            sage: psi = phi.velu(isog)
            sage: psi
            Drinfeld module defined by T |--> (z12^11 + 3*z12^10 + z12^9 + z12^7 + z12^5 + 4*z12^4 + 4*z12^3 + z12^2 + 1)*t^2 + (2*z12^11 + 4*z12^10 + 2*z12^8 + z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + z12^2 + z12 + 4)*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: isog in Hom(phi, psi)
            True

        This method works for endomorphisms as well::

            sage: phi.velu(phi(T)) is phi
            True
            sage: phi.velu(t^6) is phi
            True

        The following inputs do not define isogenies, and the method
        returns ``None``::

            sage: phi.velu(0)
            Traceback (most recent call last):
            ...
            ValueError: the input does not define an isogeny
            sage: phi.velu(t)
            Traceback (most recent call last):
            ...
            ValueError: the input does not define an isogeny
            sage: phi.velu(t^3 + t + 2)
            Traceback (most recent call last):
            ...
            ValueError: the input does not define an isogeny
        """
        if isog not in self.ore_polring():
            raise TypeError('input must be an Ore polynomial')
        e = ValueError('the input does not define an isogeny')
        if isog == 0:
            raise e
        quo, rem = (isog * self.gen()).right_quo_rem(isog)
        char_deg = self.characteristic().degree()
        if not char_deg.divides(isog.valuation()) \
                or rem != 0:
            raise e
        else:
            return self.category().object(quo)
