# sage.doctest: optional - sage.rings.finite_rings
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
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name
from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_string import _LazyString
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.ring_extension import RingExtension_generic
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')


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
        t^6 + (z^11 + z^9 + 2*z^6 + 2*z^4 + 2*z + 1)*t^4
        + (2*z^11 + 2*z^10 + z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^3)*t^3
        + (2*z^11 + z^10 + z^9 + 2*z^7 + 2*z^6 + z^5 + z^4 + 2*z^3 + 2*z + 2)*t^2
        + (2*z^11 + 2*z^8 + 2*z^6 + z^5 + z^4 + 2*z^2)*t + z^3 + z + 1
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
        Ore Polynomial Ring in t over Finite Field in z of size 3^12 over its base
         twisted by Frob^2

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
          To:   Ore Polynomial Ring in t
                over Finite Field in z of size 3^12 over its base
                twisted by Frob^2
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
        Drinfeld module defined by T |--> (2*z^11 + 2*z^9 + z^6 + 2*z^5 + 2*z^4 + 2*z^2 + 1)*t^2
         + (2*z^11 + 2*z^10 + 2*z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z)*t + z
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
            Action on Finite Field in z of size 3^12 over its base
             induced by Drinfeld module defined by T |--> t^2 + t + z

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
        Return the set of morphisms from ``self`` to ``other``.

        Validity of the input is checked at the instantiation of
        ``DrinfeldModuleHomset``; ``self`` and ``other`` only need be in
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
            Action on Finite Field in z12 of size 5^12 over its base
             induced by Drinfeld module defined by T |--> z12^5*t^2 + z12^3*t + 2*z12^11
              + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12

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
            2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5
            + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
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
        n = Integer(n)
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
            [2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7
               + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12,
             z12^3,
             z12^5]

        Careful, the method only returns the nonzero coefficients,
        unless otherwise specified::

            sage: rho = DrinfeldModule(A, [p_root, 0, 0, 0, 1])
            sage: rho.coefficients()
            [2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7
               + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12,
             1]
            sage: rho.coefficients(sparse=False)
            [2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7
               + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12,
             0,
             0,
             0,
             1]
        """
        return self._gen.coefficients(sparse=sparse)

    @cached_method
    def _compute_coefficient_exp(self, k):
        r"""
        Return the `q^k`-th coefficient of the exponential of this Drinfeld module.

        INPUT:

        - ``k`` (integer) -- the index of the coefficient

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: phi._compute_coefficient_exp(0)
            1
            sage: phi._compute_coefficient_exp(1)
            1/(T^2 + T)
            sage: phi._compute_coefficient_exp(2)
            1/(T^8 + T^6 + T^5 + T^3)
            sage: phi._compute_coefficient_exp(3)
            1/(T^24 + T^20 + T^18 + T^17 + T^14 + T^13 + T^11 + T^7)
        """
        k = ZZ(k)
        if k.is_zero():
            return self._base.one()
        q = self._Fq.cardinality()
        c = self._base.zero()
        for i in range(k):
            j = k - i
            c += self._compute_coefficient_exp(i)*self._compute_coefficient_log(j)**(q**i)
        return -c

    def exponential(self, name='z'):
        r"""
        Return the exponential of this Drinfeld module.

        Note that the exponential is only defined when the
        `\mathbb{F}_q[T]`-characteristic is zero.

        INPUT:

        - ``name`` (string, default: ``'z'``) -- the name of the
          generator of the lazy power series ring.

        OUTPUT:

        A lazy power series over the base field.

        EXAMPLES::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: exp = phi.exponential(); exp
            z + ((1/(T^2+T))*z^2) + ((1/(T^8+T^6+T^5+T^3))*z^4) + O(z^8)

        The exponential is returned as a lazy power series, meaning that
        any of its coefficients can be computed on demands::

            sage: exp[2^4]
            1/(T^64 + T^56 + T^52 + ... + T^27 + T^23 + T^15)
            sage: exp[2^5]
            1/(T^160 + T^144 + T^136 + ... + T^55 + T^47 + T^31)

        Example in higher rank::

            sage: A = GF(5)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T^2, T + T^2 + T^4, 1])
            sage: exp = phi.exponential(); exp
            z + ((T/(T^4+4))*z^5) + O(z^8)

        The exponential is the compositional inverse of the logarithm
        (see :meth:`logarithm`)::

            sage: log = phi.logarithm(); log
            z + ((4*T/(T^4+4))*z^5) + O(z^8)
            sage: exp.compose(log)
            z + O(z^8)
            sage: log.compose(exp)
            z + O(z^8)

        ::

            sage: Fq.<w> = GF(3)
            sage: A = Fq['T']
            sage: phi = DrinfeldModule(A, [w, 1])
            sage: phi.exponential()
            Traceback (most recent call last):
            ...
            ValueError: characteristic must be zero (=T + 2)

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: exp = phi.exponential()
            sage: exp[2] == 1/(T**q - T)  # expected value
            True
            sage: exp[2^2] == 1/((T**(q**2) - T)*(T**q - T)**q)  # expected value
            True
            sage: exp[2^3] == 1/((T**(q**3) - T)*(T**(q**2) - T)**q*(T**q - T)**(q**2))  # expected value
            True

        REFERENCE:

        See section 4.6 of [Gos1998]_ for the definition of the
        exponential.
        """
        if self.category()._characteristic:
            raise ValueError(f"characteristic must be zero (={self.characteristic()})")
        L = LazyPowerSeriesRing(self._base, name)
        zero = self._base.zero()
        q = self._Fq.cardinality()

        def coeff_exp(k):
            # Return the k-th coefficient of the exponential.
            k = ZZ(k)
            if k.is_power_of(q):
                return self._compute_coefficient_exp(k.log(q))
            else:
                return zero
        return L(coeff_exp, valuation=1)

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

    @cached_method
    def _compute_goss_polynomial(self, n, q, poly_ring, X):
        r"""
        Utility function for computing the n-th Goss polynomial.

        The user should not call this method directly, but
        :meth:`goss_polynomial` instead.

        TESTS::

            sage: A = GF(2^2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T+1, T^2, 1])
            sage: poly_ring = phi.base()['X']
            sage: X = poly_ring.gen()
            sage: phi._compute_goss_polynomial(0, 2^2, poly_ring, X)
            0
            sage: phi._compute_goss_polynomial(3, 2^2, poly_ring, X)
            X^3
            sage: phi._compute_goss_polynomial(4*3, 2^2, poly_ring, X)
            X^12
            sage: phi._compute_goss_polynomial(9, 2^2, poly_ring, X)
            X^9 + (1/(T^3 + T^2 + T))*X^6 + (1/(T^6 + T^4 + T^2))*X^3

        """
        # Trivial cases
        if n.is_zero():
            return poly_ring.zero()
        if n <= q - 1:
            return X**n
        if n%q == 0:
            return self.goss_polynomial(ZZ(n/q))**q
        # General case
        pol = sum(self._compute_coefficient_exp(i+1)
                  *self._compute_goss_polynomial(n - q**(i+1), q, poly_ring, X)
                  for i in range(0, (n.log(q).n()).floor()))
        return X*(self._compute_goss_polynomial(n - 1, q, poly_ring, X) + pol)

    def goss_polynomial(self, n, var='X'):
        r"""
        Return the `n`-th Goss polynomial of the Drinfeld module.

        Note that Goss polynomials are only defined for Drinfeld modules
        of characteristic zero.

        INPUT:

        - ``n`` (integer) -- the index of the Goss polynomial

        - ``name`` (str, default: ``'X'``) -- the name of polynomial
          variable.

        OUTPUT:

        - a univariate polynomial in ``name`` over the base `A`-field.

        EXAMPLES::

            sage: from drinfeld_modular_forms import goss_polynomial
            sage: A = GF(3)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, K.one()])  # Carlitz module
            sage: phi.goss_polynomial(1)
            X
            sage: phi.goss_polynomial(2)
            X^2
            sage: phi.goss_polynomial(4)
            X^4 + (1/(T^3 + 2*T))*X^2
            sage: phi.goss_polynomial(5)
            X^5 + (2/(T^3 + 2*T))*X^3
            sage: phi.goss_polynomial(10)
            X^10 + (1/(T^3 + 2*T))*X^8 + (1/(T^6 + T^4 + T^2))*X^6 + (1/(T^9 + 2*T^3))*X^4 + (1/(T^18 + 2*T^12 + 2*T^10 + T^4))*X^2

        REFERENCE:

        Section 3 of [Gek1988]_ provides an exposition of Goss
        polynomials.
        """
        if self.category()._characteristic:
            raise ValueError(f"characteristic must be zero (={self.characteristic()})")
        n = ZZ(n)
        K = self.base()
        poly_ring = K[var]
        X = poly_ring.gen()
        q = self._Fq.cardinality()
        return self._compute_goss_polynomial(n, q, poly_ring, X)

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

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.height()
            2
            sage: phi.is_supersingular()
            True

        In characteristic zero, height is not defined::

            sage: L = A.fraction_field()
            sage: phi = DrinfeldModule(A, [L(T), L(1)])
            sage: phi.height()
            Traceback (most recent call last):
            ...
            ValueError: height is only defined for prime function field characteristic

        TESTS:

        In the following case, sage is unable to determine the
        characteristic; that is why an error is raised::

            sage: B.<Y> = Fq[]
            sage: L = Frac(B)
            sage: phi = DrinfeldModule(A, [L(2), L(1)])
            sage: phi.height()
            Traceback (most recent call last):
            ...
            NotImplementedError: height not implemented in this case

        """
        try:
            if self.characteristic().is_zero():
                raise ValueError('height is only defined for prime '
                                 'function field characteristic')
            else:
                p = self.characteristic()
                return Integer(self(p).valuation() // p.degree())
        except NotImplementedError:
            raise NotImplementedError('height not implemented in this case')

    def is_isomorphic(self, other, absolutely=False):
        r"""
        Return ``True`` if this Drinfeld module is isomorphic to ``other``;
        return ``False`` otherwise.

        INPUT:

        - ``absolutely`` -- a boolean (default: ``False``); if ``True``,
          check the existence of an isomorphism defined on the base
          field; if ``False``, check over an algebraic closure.

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: t = phi.ore_variable()

        We create a second Drinfeld module, which is isomorphic to `\phi`
        and then check that they are indeed isomorphic::

            sage: psi = phi.velu(z)
            sage: phi.is_isomorphic(psi)
            True

        In the example below, `\phi` and `\psi` are isogenous but not
        isomorphic::

            sage: psi = phi.velu(t + 1)
            sage: phi.is_isomorphic(psi)
            False

        Here is an example of two Drinfeld modules which are isomorphic
        on an algebraic closure but not on the base field::

            sage: phi = DrinfeldModule(A, [z, 1])
            sage: psi = DrinfeldModule(A, [z, z])
            sage: phi.is_isomorphic(psi)
            False
            sage: phi.is_isomorphic(psi, absolutely=True)
            True

        On certain fields, testing isomorphisms over the base field may
        fail::

            sage: L = A.fraction_field()
            sage: T = L.gen()
            sage: phi = DrinfeldModule(A, [T, 0, 1])
            sage: psi = DrinfeldModule(A, [T, 0, T])
            sage: psi.is_isomorphic(phi)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot solve the equation u^24 == T

        However, it never fails over the algebraic closure::

            sage: psi.is_isomorphic(phi, absolutely=True)
            True

        Note finally that when the constant coefficients of `\phi_T` and
        `\psi_T` differ, `\phi` and `\psi` do not belong to the same category
        and checking whether they are isomorphic does not make sense; in this
        case, an error is raised::

            sage: phi = DrinfeldModule(A, [z, 0, 1])
            sage: psi = DrinfeldModule(A, [z^2, 0, 1])
            sage: phi.is_isomorphic(psi)
            Traceback (most recent call last):
            ...
            ValueError: Drinfeld modules are not in the same category

        TESTS:

        A Drinfeld module is always isomorphic to itself::

            sage: phi = DrinfeldModule(A, [z] + [K.random_element() for _ in range(3)] + [1])
            sage: phi.is_isomorphic(phi)
            True

        Two Drinfeld modules of different ranks are never isomorphic::

            sage: psi = DrinfeldModule(A, [z] + [K.random_element() for _ in range(5)] + [1])
            sage: phi.is_isomorphic(psi)
            False

        Two Drinfeld modules which are not patterned-alike are also not
        isomorphic::

            sage: phi = DrinfeldModule(A, [T, 1, 0, 1])
            sage: psi = DrinfeldModule(A, [T, 1, 1, 1])
            sage: phi.is_isomorphic(psi)
            False
        """
        # Trivial cases
        if self.category() is not other.category():
            raise ValueError("Drinfeld modules are not in the same category")
        if self is other:
            return True
        r = self.rank()
        if other.rank() != r:
            return False
        # Check if there exists u such that u*self_X = other_X*u:
        # if self_X = a_0 + a_1*t + ... + a_r*t^r
        # and other_x = b_0 + b_1*t + ... + b_r*t^r
        # this reduces to find a solution of the system:
        #   b_i * u^(q^i - 1) = a_i  (0 <= i <= r)
        # which, using gcds, can be reduced to a unique equation
        # of the form u^e = ue.
        q = self._Fq.cardinality()
        A = self.gen()
        B = other.gen()
        e = Integer(0)
        ue = self._base(1)
        for i in range(1, r+1):
            ai = A[i]
            bi = B[i]
            if ai == 0 and bi == 0:
                continue
            if ai == 0 or bi == 0:
                return False
            if e != q - 1:
                # u^e = ue
                # u^(q^i - 1) = ai/bi
                e, s, t = e.xgcd(q**i - 1)
                ue = ue**s * (ai/bi)**t
        for i in range(1, r+1):
            if A[i]:
                f = (q**i - 1) // e
                if A[i] != B[i] * ue**f:
                    return False
        # Solve the equation u^e = ue
        # - when absolutely=True, over the algebraic closure (then a
        #   solution always exists)
        # - when absolutely=False, over the ground field.
        if absolutely:
            return True
        else:
            ue = ue.backend(force=True)
            try:
                _ = ue.nth_root(e)
            except ValueError:
                return False
            except (AttributeError, NotImplementedError):
                raise NotImplementedError(f"cannot solve the equation u^{e} == {ue}")
            return True

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

    @cached_method
    def _compute_coefficient_log(self, k):
        r"""
        Return the `q^k`-th coefficient of the logarithm of this Drinfeld module.

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: phi._compute_coefficient_log(0)
            1
            sage: phi._compute_coefficient_log(1)
            1/(T^2 + T)
            sage: phi._compute_coefficient_log(2)
            1/(T^6 + T^5 + T^3 + T^2)
            sage: phi._compute_coefficient_log(3)
            1/(T^14 + T^13 + T^11 + T^10 + T^7 + T^6 + T^4 + T^3)
        """
        k = ZZ(k)
        if k.is_zero():
            return self._base.one()
        r = self._gen.degree()
        T = self._gen[0]
        q = self._Fq.cardinality()
        c = self._base.zero()
        for i in range(k):
            j = k - i
            if j < r + 1:
                c += self._compute_coefficient_log(i)*self._gen[j]**(q**i)
        return c/(T - T**(q**k))

    def logarithm(self, name='z'):
        r"""
        Return the logarithm of the given Drinfeld module.

        By definition, the logarithm is the compositional inverse of the
        exponential (see :meth:`exponential`). Note that the logarithm
        is only defined when the `\mathbb{F}_q[T]`-characteristic is
        zero.

        INPUT:

        - ``name`` (string, default: ``'z'``) -- the name of the
          generator of the lazy power series ring.

        OUTPUT:

        A lazy power series over the base field.

        EXAMPLES::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: log = phi.logarithm(); log
            z + ((1/(T^2+T))*z^2) + ((1/(T^6+T^5+T^3+T^2))*z^4) + O(z^8)

        The logarithm is returned as a lazy power series, meaning that
        any of its coefficients can be computed on demands::

            sage: log[2^4]
            1/(T^30 + T^29 + T^27 + ... + T^7 + T^5 + T^4)
            sage: log[2^5]
            1/(T^62 + T^61 + T^59 + ... + T^8 + T^6 + T^5)

        Example in higher rank::

            sage: A = GF(5)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T^2, T + T^2 + T^4, 1])
            sage: phi.logarithm()
            z + ((4*T/(T^4+4))*z^5) + O(z^8)

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = 2
            sage: log[2] == -1/((T**q - T))  # expected value
            True
            sage: log[2**2] == 1/((T**q - T)*(T**(q**2) - T))  # expected value
            True
            sage: log[2**3] == -1/((T**q - T)*(T**(q**2) - T)*(T**(q**3) - T))  # expected value
            True

        ::

            sage: Fq.<w> = GF(3)
            sage: A = Fq['T']
            sage: phi = DrinfeldModule(A, [w, 1])
            sage: phi.logarithm()
            Traceback (most recent call last):
            ...
            ValueError: characteristic must be zero (=T + 2)
        """
        if self.category()._characteristic:
            raise ValueError(f"characteristic must be zero (={self.characteristic()})")
        L = LazyPowerSeriesRing(self._base, name)
        zero = self._base.zero()
        q = self._Fq.cardinality()

        def coeff_log(k):
            # Return the k-th coefficient of the logarithm
            k = ZZ(k)
            if k.is_power_of(q):
                return self._compute_coefficient_log(k.log(q))
            else:
                return self._base.zero()
        return L(coeff_log, valuation=1)

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
              To:   Ore Polynomial Ring in t over Finite Field in z12 of size 5^12
                    over its base twisted by Frob^2
              Defn: T |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8
                           + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
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
            [z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8
             + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12]
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

        TESTS:

        The rank must be an ``Integer`` (see PR #35519)::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: isinstance(phi.rank(), Integer)
            True
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
            Drinfeld module defined by T |-->
             (z12^11 + 3*z12^10 + z12^9 + z12^7 + z12^5 + 4*z12^4 + 4*z12^3 + z12^2 + 1)*t^2
             + (2*z12^11 + 4*z12^10 + 2*z12^8 + z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + z12^2 + z12 + 4)*t
             + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
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
        if rem.is_zero() and quo[0] == self.gen()[0]:
            return self.category().object(quo)
        else:
            raise e

    def hom(self, x, codomain=None):
        r"""
        Return the homomorphism defined by ``x`` having this Drinfeld
        module as domain.

        We recall that a homomorphism `f : \phi \to \psi` between
        two Drinfeld modules is defined by an Ore polynomial `u`,
        which is subject to the relation `phi_T u = u \psi_T`.

        INPUT:

        - ``x`` -- an element of the ring of functions, or an
          Ore polynomial

        - ``codomain`` -- a Drinfeld module or ``None`` (default:
          ``None``)

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: phi
            Drinfeld module defined by T |--> z*t^3 + t^2 + z

        An important class of endomorphisms of a Drinfeld module
        `\phi` is given by scalar multiplications, that are endomorphisms
        corresponding to the Ore polynomials `\phi_a` with `a` in the function
        ring `A`. We construct them as follows::

            sage: phi.hom(T)
            Endomorphism of Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: z*t^3 + t^2 + z

        ::

            sage: phi.hom(T^2 + 1)
            Endomorphism of Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: z^2*t^6 + (3*z^2 + z + 1)*t^5 + t^4 + 2*z^2*t^3 + (3*z^2 + z + 1)*t^2 + z^2 + 1

        We can also define a morphism by passing in the Ore polynomial
        defining it.
        For example, below, we construct the Frobenius endomorphism
        of `\phi`::

            sage: t = phi.ore_variable()
            sage: phi.hom(t^3)
            Endomorphism of Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: t^3

        If the input Ore polynomial defines a morphism to another
        Drinfeld module, the latter is determined automatically::

            sage: phi.hom(t + 1)
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> z*t^3 + t^2 + z
              To:   Drinfeld module defined by T |--> (2*z^2 + 4*z + 4)*t^3 + (3*z^2 + 2*z + 2)*t^2 + (2*z^2 + 3*z + 4)*t + z
              Defn: t + 1

        TESTS::

            sage: phi.hom(t)
            Traceback (most recent call last):
            ...
            ValueError: the input does not define an isogeny

        ::

            sage: phi.hom(T + z)
            Traceback (most recent call last):
            ...
            ValueError: the input does not define an isogeny

        ::

            sage: phi.hom(t + 1, codomain=phi)
            Traceback (most recent call last):
            ...
            ValueError: Ore polynomial does not define a morphism

        """
        # When `x` is in the function ring (or something that coerces to it):
        if self.function_ring().has_coerce_map_from(x.parent()):
            return self.Hom(self)(x)
        if codomain is None:
            try:
                codomain = self.velu(x)
            except TypeError:
                raise ValueError("the input does not define an isogeny")
        H = self.Hom(codomain)
        return H(x)

        # TODO: implement the method `moh`, the analogue of `hom`
        # with fixed codomain.
        # It's not straightforward because `left_divides` does not
        # work currently because Sage does not know how to invert the
        # Frobenius endomorphism; this is due to the RingExtension stuff.

    def scalar_multiplication(self, x):
        r"""
        Return the endomorphism of this Drinfeld module, which is
        the multiplication by `x`, i.e. the isogeny defined by the
        Ore polynomial `\phi_x`.

        INPUT:

        - ``x`` -- an element in the ring of functions

        EXAMPLES::

            sage: Fq = GF(5)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: phi = DrinfeldModule(A, [z, 0, 1, z])
            sage: phi
            Drinfeld module defined by T |--> z*t^3 + t^2 + z
            sage: phi.hom(T)
            Endomorphism of Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: z*t^3 + t^2 + z

        ::

            sage: phi.hom(T^2 + 1)
            Endomorphism of Drinfeld module defined by T |--> z*t^3 + t^2 + z
              Defn: z^2*t^6 + (3*z^2 + z + 1)*t^5 + t^4 + 2*z^2*t^3 + (3*z^2 + z + 1)*t^2 + z^2 + 1

        """
        if not self.function_ring().has_coerce_map_from(x.parent()):
            raise ValueError("%s is not element of the function ring" % x)
        return self.Hom(self)(x)
