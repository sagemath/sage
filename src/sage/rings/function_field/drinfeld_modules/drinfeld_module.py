r"""
Drinfeld modules

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.drinfeld_module.DrinfeldModule`.
For *finite* Drinfeld modules and their theory of complex multiplication, see
class
:class:`sage.rings.function_field.drinfeld_module.finite_drinfeld_module.DrinfeldModule`.

AUTHORS:

- Antoine Leudière (2022-04)
- Xavier Caruso (2022-06)
"""

#*****************************************************************************
#       Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.category_object import CategoryObject

from sage.categories.drinfeld_modules import DrinfeldModules
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general

from sage.misc.latex import latex
from sage.structure.sequence import Sequence
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector

class DrinfeldModule(UniqueRepresentation, CategoryObject):
    r"""
    This class represents a Drinfeld module.

    Let `\Fq` be a finite field with order `q`. Let `K` be a field
    equiped a ring morphism `\gamma: \Fq[X] \to K` --- the field `K` is
    said to be an *`\Fq[X]`-field*, and the monic polynomial that
    generates `\Ker(\gamma)` is called the *`\Fq[X]`-characteristic of
    the `\Fq[X]`-field `K`* (this characteristic plays the role of the
    characteristic of a function field). Let `K\{\tau\}` be the ring of
    Ore polynomials with coefficients in `K` and Frobenius
    variable `\tau: x \mapsto x^q`. A *Drinfeld `\Fq[X]`-module over the
    `\Fq[X]`-field `K`* is a ring morphism `\phi: \Fq[X] \to K\{\tau\}`
    such that:

        1. The image of `\phi` contains non-constant Ore polynomials.
        2. For every `a \in \Fq[X]`, the constant coefficient `\phi(a)`
           is `\gamma(a)`.

    For `a \in \Fq[X]`, `\phi(a)` is denoted `\phi_a`.

    We say that `\Fq[X]` is the *function ring of `\phi`*; `K` is the
    *base ring of `\phi`*, or simply its base or base field; `\Fq[X]` is
    the *function ring of `\phi`*; *K\{\tau\}* is the *Ore polynomial
    ring of `\phi`*; `t` is the *Ore variable of `\phi`*. The *generator of `\phi`* is
    `\phi_X`, its *constant coefficient* is the constant coefficient of
    `\phi_X`.

    The Drinfeld module `\phi` is uniquely determined by the image
    `\phi_X` of `X`. This Ore polynomial is an input of the class
    constructor.

    Classical references on Drinfeld modules include [Gos1998]_,
    [Rosen2002]_, [VS06]_ and [Gek1998]_.

    .. NOTE::

        Drinfeld modules are defined in a larger setting, in which
        `\Fq[X]` is replaced by a more general ring: the ring of
        functions in `k` that are regular outside `\infty`, where `k` is
        a function field over `\Fq` with transcendance degree `1` and
        `\infty` is a fixed place of `k`. This is out of the scope of
        this implementation.

    INPUT:

    - ``function_ring`` -- the polynomial ring `\Fq[X]`
    - ``gen`` -- the generator `\phi_X`, as a list of coefficients or an
      Ore polynomial
    - ``name`` (optional) the name of the Ore variable

    .. RUBRIC:: Construction and input

    A Drinfeld module object is constructed as follows::

        sage: Fq.<z2> = GF(3^2)
        sage: FqX.<X> = Fq[]
        sage: K.<z> = Fq.extension(6)
        sage: phi = DrinfeldModule(FqX, [z, 1, 1])
        sage: phi
        Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12

    In this example, we used a list of coefficients (``[z, 1, 1]``) to
    represent the generator `\phi_X = z + t + t^2`, `K = \Fq(z)`. One can
    also use regular Ore polynomials::

        sage: ore_polring = phi.ore_polring()
        sage: t = phi.ore_variable()  # same as ore_polring.gen()
        sage: psi_X = z + t^3
        sage: psi = DrinfeldModule(FqX, psi_X)
        sage: psi
        Drinfeld module defined by X |--> t^3 + z over Finite Field in z of size 3^12
        sage: psi(X) == psi_X
        True

    The generator must have positive degree::

        sage: DrinfeldModule(FqX, [z])
        Traceback (most recent call last):
        ...
        ValueError: generator must have positive degree

    The coefficients of the generator must live in some field `K` that
    is the codomain of a morphism `\gamma` with domain `\Fq[X]`::

        sage: DrinfeldModule(FqX, [z, QQ(1)])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce to base ring

        sage: DrinfeldModule(FqX, [1, QQ(1)])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce to base ring

    If the coefficients are regular integers, an exception is raised.
    One needs to manually cast them to the field of their choice::

        sage: DrinfeldModule(FqX, [1, 1, 1])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce to base ring

        sage: DrinfeldModule(FqX, [K(1), 1, 1])
        Drinfeld module defined by X |--> t^2 + t + 1 over Finite Field in z of size 3^12

        sage: DrinfeldModule(FqX, [Fq(1), 1, 1])
        Drinfeld module defined by X |--> t^2 + t + 1 over Finite Field in z2 of size 3^2

    The function ring must be an `\Fq[X]`::

        sage: DrinfeldModule(K, [z, 1, 1])
        Traceback (most recent call last):
        ...
        NotImplementedError: function ring must be a polynomial ring

        sage: FqXY.<Y> = FqX[]
        sage: DrinfeldModule(FqXY, [z, 1, 1])
        Traceback (most recent call last):
        ...
        TypeError: function ring base must be a finite field

    If you already defined a category of Drinfeld modules, you must
    ensure that the constant coefficient is a root of the
    `\Fq[X]-characteristic of the category base::

        sage: cat = phi.category()
        sage: cat([1, 1, K(1)])
        Traceback (most recent call last):
        ...
        ValueError: constant coefficient must be a root of the characteristic

    .. NOTE::

        The reader may think that it is odd to build a Drinfeld module
        without specifying the `\Fq[X]`-characteristic of the base (or,
        more generally, its parent category). Indeed, the
        `\Fq[X]`-characteristic plays the role of the characteristic of
        a function field, and thus preexists Drinfeld modules. The base
        field `K` should rather be seen as an `\Fq[X]`-field, i.e. the
        field `K` equiped with a morphism `\gamma: \Fq[X] \to K`,
        instead of just a field.

        However, as the characteristic may be deduced from the constant
        coefficient of the Drinfeld module, we chose to ommit the
        characteristic in the input of the class in order to have
        concise definitions.

    .. RUBRIC:: Possible base rings

    The morphism `\gamma` does not need be surjective like in the above
    examples. This is equivalent to say that the constant coefficient of
    the Drinfeld module may be different to the generator of `K` over
    `\Fq`. In the following example, `K` is still a degree six extension
    of `\Fq`, but `\gamma` is a projection over a degree two extension
    with modulus `X^3 + (z_2 + 2)X^2 + (6*z_2 + 1)X + 3z_2 + 5`::

        sage: p = X^2 + z2 + 2
        sage: p_root = z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z
        sage: rho = DrinfeldModule(FqX, [p_root, 1, 1])
        sage: rho
        Drinfeld module defined by X |--> t^2 + t + z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z over Finite Field in z of size 3^12

    Then one can check that the morphisms `\gamma` are not the same for
    ``phi`` and ``rho``, and that the `\gamma` associated to `\phi` is
    surjective, while the other one is not::

        sage: rho.category().morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z
        sage: phi.category().morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z

    Note that ``phi`` and ``psi`` are *finite* Drinfeld modules, in the
    sense that `K` is finite. But `K` can be infinite::

        sage: sigma = DrinfeldModule(FqX, [Frac(FqX).gen(), 1, 1])
        sage: sigma
        Drinfeld module defined by X |--> t^2 + t + X over Fraction Field of Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
        sage: sigma.is_finite()
        False
        sage: phi.is_finite()
        True

    It is possible to change the base ring::

        sage: L = K.extension(2)
        sage: phi_rebased = phi.change_ring(L)
        sage: Ltau = phi_rebased.ore_polring()
        sage: Ltau(phi(X)) == phi_rebased(X)
        True

    .. RUBRIC:: The category of Drinfeld modules

    Drinfeld modules have their own category (see class
    :class:`sage.categories.drinfeld_modules.DrinfeldModules`)::

        sage: phi.category()
        Category of Drinfeld modules defined by Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z
        sage: phi.category() is psi.category()
        True
        sage: phi.category() is sigma.category()
        False

    This category holds crucial information, like the
    `\Fq[X]`-characteristic of `K`::

        sage: char = phi.category().characteristic()

    As the output of
    :meth:`sage.rings.function_field.drinfeld_module.finite_drinfeld_module.DrinfeldModule.category`
    suggests, the morphism `\gamma` uniquely determines the category of a Drinfeld
    module.

    .. RUBRIC:: Basic methods

    For a polynomial `a \in \Fq[X]`, compute `\phi_a` by calling `phi`::

        sage: phi(X)  # phi_X
        t^2 + t + z
        sage: phi(X^3 + X + 1)  # phi_X^3 +X + 1
        t^6 + (z^11 + z^9 + 2*z^6 + 2*z^4 + 2*z + 1)*t^4 + (2*z^11 + 2*z^10 + z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^3)*t^3 + (2*z^11 + z^10 + z^9 + 2*z^7 + 2*z^6 + z^5 + z^4 + 2*z^3 + 2*z + 2)*t^2 + (2*z^11 + 2*z^8 + 2*z^6 + z^5 + z^4 + 2*z^2)*t + z^3 + z + 1
        sage: phi(1)  # phi_1
        1

    One can retrieve basic properties::

        sage: phi.base_ring()  # K
        Finite Field in z of size 3^12

        sage: phi.ore_polring()  # K{t}
        Ore Polynomial Ring in t over Finite Field in z of size 3^12 twisted by z |--> z^(3^2)

        sage: phi.ore_variable()  # t
        t

        sage: phi.function_ring()  # Fq[X]
        Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2

        sage: phi.gen()  # phi_X
        t^2 + t + z
        sage: phi.gen() == phi(X)
        True

        sage: phi.constant_coefficient()  # Constant coefficient of phi_X
        z

        sage: phi.morphism()  # The Drinfeld module as a morphism
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Ore Polynomial Ring in t over Finite Field in z of size 3^12 twisted by z |--> z^(3^2)
          Defn: X |--> t^2 + t + z

    One can compute the rank and height::

        sage: phi.rank()
        2
        sage: phi.height()
        1

    As well as the j-invariant if the rank is two::

        sage: phi.j_invariant()  # j-invariant
        1

    .. RUBRIC:: Morphisms, isogenies

    A *morphism of Drinfeld modules `\phi \to \psi`* is an Ore
    polynomial `f \in K\{\tau\}` such that `f \phi_a = \psi_a f` for
    every `a \in \Fq[X]`. In our case, this is equivalent to verifying
    `f \phi_X = \psi_X f`. An *isogeny* is a non-zero morphism.

    Use the ``in`` syntax to test if an Ore polynomial defines a
    morphism::

        sage: phi(X) in Hom(phi, phi)
        True
        sage: t^6 in Hom(phi, phi)
        True
        sage: t^5 + 2*t^3 + 1 in Hom(phi, phi)
        False
        sage: 1 in Hom(phi, psi)
        False
        sage: 1 in Hom(phi, phi)
        True
        sage: 0 in Hom(phi, psi)
        True

    To create a SageMath object representing the morphism, call the
    homset (``hom`` in the next example)::

        sage: hom = Hom(phi, phi)
        sage: frobenius_endomorphism = hom(t^6)
        sage: identity_morphism = hom(1)
        sage: zero_morphism = hom(0)
        sage: frobenius_endomorphism
        Drinfeld Module morphism:
          From: Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12
          To:   Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12
          Defn: t^6
        sage: identity_morphism
        Drinfeld Module morphism:
          From: Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12
          To:   Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12
          Defn: 1
        sage: zero_morphism
        Drinfeld Module morphism:
          From: Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12
          To:   Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12
          Defn: 0

    One can retrieve the underlying Ore polynomial with the method
    :meth:`sage.rings.function_field.drinfeld_module.finite_drinfeld_module.DrinfeldModule.ore_polynomial`::

        sage: frobenius_endomorphism.ore_polynomial()
        t^6

    And one can easily check if a morphism defines an isogeny or an
    isomorphism (i.e. an isogeny whose underlying Ore polynomial has
    degree `0`)::

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

    Let ``ore_pol`` be a non-zero Ore polynomial. For Drinfeld module,
    it is easy to decide if there exists a Drinfeld module ``psi`` such
    that ``ore_pol`` is an isogeny from ``self`` to ``psi``. If so, we
    find ``psi``::

        sage: ore_pol = (2*z^6 + z^3 + 2*z^2 + z + 2)*t + z^11 + 2*z^10 + 2*z^9 + 2*z^8 + z^7 + 2*z^6 + z^5 + z^3 + z^2 + z
        sage: psi = phi.velu(ore_pol)
        sage: psi
        Drinfeld module defined by X |--> (2*z^11 + 2*z^9 + z^6 + 2*z^5 + 2*z^4 + 2*z^2 + 1)*t^2 + (2*z^11 + 2*z^10 + 2*z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z)*t + z over Finite Field in z of size 3^12
        sage: ore_pol in Hom(phi, psi)
        True
        sage: ore_pol * phi(X) == psi(X) * ore_pol
        True

    If the input does not define an isogeny, an exception is raised:

        sage: phi.velu(0)
        Traceback (most recent call last):
        ...
        ValueError: the input does not define an isogeny
        sage: phi.velu(t)
        Traceback (most recent call last):
        ...
        ValueError: the input does not define an isogeny

    .. RUBRIC:: The action of a Drinfeld module

    An `\Fq[X]`-Drinfeld module `\phi` notoriously makes any field
    extension `L/K` a left `\Fq[X]`-module. Let `x \in L`, let `P \in
    \Fq[X]`; the action is defined as `(P, a) \mapsto \phi_P(a)`, where
    `\phi_P(a)`. The method :meth:`action` returns an ``Action`` object
    representing the Drinfeld module action; in this implementation, `K
    = L`.

        sage: action = phi.action()
        sage: action
        Action on Finite Field in z of size 3^12 induced by Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12

    The action on elements is computed as follows::

        sage: P = X + 1
        sage: a = z
        sage: action(P, a)
        ...
        z^9 + 2*z^8 + 2*z^7 + 2*z^6 + 2*z^3 + z^2
        sage: action(0, K.random_element())
        0
        sage: action(FqX.random_element(), 0)
        0

    To act on a field larger than `K`, one can change the ring of the
    Drinfeld module, then create the action::

        sage: extended_action = phi.change_ring(K.extension(5)).action()
        sage: extended_action
        Action on Finite Field in z60 of size 3^60 induced by Drinfeld module defined by X |--> t^2 + t + 2*z60^59 + z60^56 + 2*z60^55 + 2*z60^54 + 2*z60^53 + z60^49 + z60^48 + z60^47 + 2*z60^45 + z60^44 + 2*z60^41 + 2*z60^40 + 2*z60^39 + 2*z60^37 + 2*z60^36 + z60^34 + z60^33 + z60^32 + 2*z60^31 + 2*z60^30 + 2*z60^27 + 2*z60^25 + z60^23 + z60^22 + z60^21 + 2*z60^20 + z60^19 + z60^18 + z60^17 + z60^16 + z60^15 + 2*z60^14 + z60^12 + 2*z60^11 + 2*z60^10 + z60^8 + z60^6 + 2*z60^5 + z60^4 + z60^3 + z60 + 1 over Finite Field in z60 of size 3^60

    .. RUBRIC:: Inversion of the Drinfeld module

    Given an Ore polynomial that equals `\phi_a` for some `a \in
    \Fq[X]`, one can retrieve `a` (as a morphism, a Drinfeld
    module is injective, see [Gos1998]_, cor. 4.5.2.)::

        sage: a = FqX.random_element()
        sage: phi.invert(phi(a)) == a
        True
    """

    @staticmethod
    def __classcall_private__(cls, function_ring, gen, name='t'):

        # FIXME: function_ring must be checked before calling base_ring
        # on it. But then it is checked twice: firstly here, secondly in
        # the category. Another problem is that those lines are
        # duplicate. As a general comment, there are sanity checks both
        # here and in the category constructor, which is not ideal.
        # Check domain is Fq[X]
        if not isinstance(function_ring, PolynomialRing_general):
            raise NotImplementedError('function ring must be a polynomial ring')
        function_ring_base = function_ring.base_ring()
        if not function_ring_base.is_field() or not function_ring_base.is_finite() :
            raise TypeError('function ring base must be a finite field')
        Fq = function_ring_base
        FqX = function_ring
        X = FqX.gen()

        # Check all possible input types for gen
        # `gen` is an Ore polynomial:
        if isinstance(gen, OrePolynomial):
            ore_polring = gen.parent()
            ore_polring_base = ore_polring.base_ring()
            name = ore_polring.variable_name()
        # `gen` is a list of coefficients (function_ring = Fq[X]):
        elif isinstance(gen, (list, tuple)):
            ore_polring = None
            ore_polring_base = Sequence(gen).universe()
        else:
            raise TypeError('generator must be list of coefficients or Ore ' \
                    'polynomial')
        # The coefficients are in a base ring that has coercion from Fq:
        if not (hasattr(ore_polring_base, 'has_coerce_map_from') and \
                ore_polring_base.has_coerce_map_from(function_ring.base_ring())):
            raise ValueError('function ring base must coerce to base ring')

        # Build the morphism that defines the category
        gamma = function_ring.hom([ore_polring_base(gen[0])])

        # Other checks in the category definition
        category = DrinfeldModules(gamma, name=name)

        # Check gen as Ore polynomial
        if ore_polring is not None and \
                ore_polring != category.ore_polring():
            raise ValueError(f'generator must lie in {category.ore_polring()}')
        ore_polring = category.ore_polring()  # Sanity cast
        gen = ore_polring(gen)
        if gen.degree() <= 0:
            raise ValueError('generator must have positive degree')

        # Instantiate the appropriate class
        if ore_polring_base.is_finite():
            from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
            return FiniteDrinfeldModule(gen, category)
        return cls.__classcall__(cls, gen, category)

    def __init__(self, gen, category):
        super().__init__(category=category)
        self._base_ring = category.base()
        self._function_ring = category.function_ring()
        self._gen = gen
        self._morphism = category._function_ring.hom([gen])
        self._ore_polring = gen.parent()
        self._Fq = self._function_ring.base_ring()  # Must be last

    def __call__(self, a):
        r"""
        Return the image of ``a`` by the morphism that defines the
        Drinfeld module, i.e. `\phi_a` if the Drinfeld module is denoted
        `phi`.

        INPUT:

        - ``a`` -- an element in the function ring

        OUTPUT: an element of the base ring
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
          morphisms, usually `self.category()`

        OUTPUT: an homset

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: t = phi.ore_variable()
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
         """
         if self.rank() != 2:
             raise NotImplementedError('rank must be 2')

    def _latex_(self):
        r"""
        Return a LaTeX representation of the Drinfeld module.

        OUTPUT: a string

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: latex(phi)
            \text{Drinfeld{ }module{ }defined{ }by{ }} X \mapsto z_{12}^{5} t^{2} + z_{12}^{3} t + 2 z_{12}^{11} + 2 z_{12}^{10} + z_{12}^{9} + 3 z_{12}^{8} + z_{12}^{7} + 2 z_{12}^{5} + 2 z_{12}^{4} + 3 z_{12}^{3} + z_{12}^{2} + 2 z_{12}\text{{ }over{ }}\Bold{F}_{5^{12}}
        """
        return f'\\text{{Drinfeld{{ }}module{{ }}defined{{ }}by{{ }}}} ' \
                f'{latex(self._function_ring.gen())} '\
                f'\\mapsto {latex(self._gen)}' \
                f'\\text{{{{ }}over{{ }}}}{latex(self._base_ring)}'

    def _repr_(self):
        r"""
        Return a string representation of the Drinfeld module.

        OUTPUT: a string

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi
            Drinfeld module defined by X |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over Finite Field in z12 of size 5^12
        """
        return f'Drinfeld module defined by {self._function_ring.gen()} ' \
                f'|--> {self._gen} over {self._base_ring}'

    def action(self):
        r"""
        Return the action object that represents the action on the base that is
        induced by the Drinfeld module.

        OUTPUT: a Drinfeld module action object

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: action = phi.action()
            sage: action
            Action on Finite Field in z12 of size 5^12 induced by Drinfeld module defined by X |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over Finite Field in z12 of size 5^12

    The action on elements is computed as follows::

            sage: P = X^2 + X + 1
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

    def base_ring(self):
        r"""
        Return the base ring of the Drinfeld module.
    
        This is the Ore polynomial ring base. In particular, the base
        ring is always a field.

        OUTPUT: a field

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.base_ring()
            Finite Field in z12 of size 5^12

        This is always true::

            sage: phi.base_ring() is phi.ore_polring().base_ring()
            True
            sage: phi.base_ring() is K
            True

        Note that in the above example, the constant coefficient
        generates a strict sub-extension of `K/\Fq`. In fact, the base
        ring may also be the same as ``Fq``::

            sage: psi = DrinfeldModule(FqX, [Fq(1), Fq.gen()])
            sage: psi.base_ring()
            Finite Field in z2 of size 5^2
            sage: psi.base_ring() is Fq
            True

        In this case the Ore polynomial ring is isomorphic to a regular
        polynomial ring::

            sage: psi.ore_polring()
            Ore Polynomial Ring in t over Finite Field in z2 of size 5^2 twisted by Identity
            sage: psi.ore_polring().twisting_morphism()
            Identity endomorphism of Finite Field in z2 of size 5^2

        TESTS::

            sage: psi.ore_polring().twisting_morphism().is_identity()
            True

            sage: psi.base_ring() is psi.function_ring().base_ring()
            True

        """
        return self._base_ring

    def change_ring(self, new_field, name=None):
        r"""
        Return a Drinfeld module defined like the current one, but with
        base ring ``new_field``.

        The new base can either be a field extension of the base ring,
        or field that has a coercion map from the field of definitions
        of the coefficients of the generator.

        INPUT:

        - ``new_field`` -- the field extension of the base ring that
          serves as base ring for the new Drinfeld module

        OUTPUT: a Drinfeld module

        EXAMPLES:

        The new ring can be an extension of the base ring::

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, p_root^3, 2])
            sage: K2 = K.extension(2)
            sage: phi_2 = phi.change_ring(K2)
            sage: phi_2
            Drinfeld module defined by X |--> 2*t^2 + (3*z24^23 + 2*z24^22 + 2*z24^20 + z24^19 + 4*z24^18 + 3*z24^17 + 4*z24^15 + 2*z24^13 + 4*z24^12 + 4*z24^11 + 3*z24^10 + 3*z24^9 + 4*z24^8 + 4*z24^6 + 3*z24^5 + 4*z24^4 + 4*z24^3 + 2*z24)*t + 2*z24^23 + 2*z24^22 + z24^21 + 2*z24^20 + z24^19 + 2*z24^18 + 3*z24^17 + 2*z24^16 + 4*z24^12 + 3*z24^11 + 4*z24^10 + z24^9 + z24^8 + 3*z24^7 + 2*z24^6 + z24^4 + 4*z24^3 + 3*z24^2 + 3*z24 + 2 over Finite Field in z24 of size 5^24

        And one can check various things::

            sage: phi.change_ring(K2).change_ring(K) is phi
            True
            sage: phi_2.base_ring() is K2
            True

        Naturally, the category has changed::

            sage: phi_2.category()
            Category of Drinfeld modules defined by Ring morphism:
              From: Univariate Polynomial Ring in X over Finite Field in z2 of size 5^2
              To:   Finite Field in z24 of size 5^24
              Defn: X |--> 2*z24^23 + 2*z24^22 + z24^21 + 2*z24^20 + z24^19 + 2*z24^18 + 3*z24^17 + 2*z24^16 + 4*z24^12 + 3*z24^11 + 4*z24^10 + z24^9 + z24^8 + 3*z24^7 + 2*z24^6 + z24^4 + 4*z24^3 + 3*z24^2 + 3*z24 + 2

        One can also change the base ring to a subfield, even though some things
        do not work as expected::

            sage: K0 = Fq.extension(2)
            sage: phi_0 = phi.change_ring(K0)
            sage: phi_0.base_ring() is K0
            True
            sage: phi.change_ring(K0).change_ring(K)  # known bug
            Traceback (most recent call last)
            ...
            TypeError: no coercion defined

        Furthermore::

            sage: phi.change_ring(K) is phi
            True
        """
        coeffs = self._gen.coefficients()
        new_coeffs = list(map(new_field, coeffs))
        if name == None:
            name = self._ore_polring.variable_name()
        return DrinfeldModule(self._function_ring, new_coeffs, name=name)

    def constant_coefficient(self):
        r"""
        Return the constant coefficient of the generator (`\phi_X`).

        From the definition of a Drinfeld module, the constant coefficient equals
        `\gamma(X)`. Hence, it is a root of the `\Fq[X]`-characteristic
        of the base ring.
    
        OUTPUT: an element in the base ring

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.constant_coefficient() == p_root
            True

        The constant coefficient equals `\gamma(X)`::

            sage: cat = phi.category()
            sage: gamma = cat.morphism()
            sage: gamma(X) == phi.constant_coefficient()
            True

        Naturally, two Drinfeld modules in the same category have the
        same constant coefficient::

            sage: t = phi.ore_variable()
            sage: psi = cat(phi.constant_coefficient() + t^3)
            sage: psi
            Drinfeld module defined by X |--> t^3 + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over Finite Field in z12 of size 5^12

        Reciprocally, it is impossible to create two Drinfeld modules in
        this category if they do not share the same constant
        coefficient::

            sage: rho = cat(phi.constant_coefficient() + 1 + t^3)
            Traceback (most recent call last):
            ...
            ValueError: constant coefficient must be a root of the characteristic

        One can also retrieve the constant coefficient using
        ``phi(X)[0]`::

            sage: phi.constant_coefficient() == phi(X)[0]
            True
        """
        return self.gen()[0]

    def function_ring(self):
        r"""
        Return the function ring of the Drinfeld module.

        In our case, the function ring is an `\Fq[X]`.

        OUTPUT: a polynomial ring

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.function_ring() is FqX
            True
        """
        return self._function_ring

    def gen(self):
        r"""
        Return the generator of the Drinfeld module, i.e. `\phi_X`.

        This method makes sense because, in our case, the function ring
        `\Fq[X]`, which is `\Fq`-generated by a single element, whose
        image characterizes the Drinfeld module.
            
        OUTPUT: an Ore polynomial

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.gen() == phi(X)
            True
        """
        return self._gen

    def height(self):
        r"""
        Return the height of the Drinfeld module.

        In our case, the height is always 1.

        OUTPUT: an integer

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.height() == 1
            True
        """
        return Integer(1)

    def invert(self, ore_pol):
        r"""
        Return the preimage of the input under the Drinfeld module;
        raise an exception if the input is not in the image of the
        Drinfeld module.

        INPUT:

        - ``ore_pol`` -- the Ore polynomial whose preimage we want to
          compute

        OUTPUT: a polynomial

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: a = FqX.random_element()
            sage: phi.invert(phi(a)) == a
            True
            sage: phi.invert(phi(X)) == X
            True
            sage: phi.invert(phi(Fq.gen())) == Fq.gen()
            True

        When the input is not in the image of the Drinfeld module, an
        exception is raised::

            sage: t = phi.ore_variable()
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

            sage: a = FqX.random_element()
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
        if not ore_pol in self._ore_polring:
            raise TypeError('input must be an Ore polynomial')
        if ore_pol in self._base_ring:
            return self._Fq(ore_pol)
        if deg % r != 0:
            raise ValueError('input must be in the image of the Drinfeld module')

        k = deg // r
        X = self._function_ring.gen()
        mat_lines = [[0 for _ in range(k+1)] for _ in range(k+1)]
        for i in range(k+1):
            phi_X_i = self(X**i)
            for j in range(i+1):
                mat_lines[j][i] = phi_X_i[r*j]
        mat = Matrix(mat_lines)
        vec = vector([list(ore_pol)[r*j] for j in range(k+1)])
        pre_image = self._function_ring(list((mat**(-1)) * vec))

        if self(pre_image) == ore_pol:
            return pre_image
        else:
            return None

    def is_finite(self):
        r"""
        Return ``True`` whether the Drinfeld module is finite.

        OUTPUT: a boolean

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.is_finite()
            True
            sage: L = Frac(FqX)
            sage: psi = DrinfeldModule(FqX, [L(2), L(1)])
            sage: psi.is_finite()
            False
        """
        from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
        return isinstance(self, FiniteDrinfeldModule)

    def j_invariant(self):
        r"""
        Return the j-invariant of the Drinfeld module; only the rank two
        case has been implemented, a NotImplementedError is raised if
        the rank is not two.

        Assume the rank is two. Write the generator `\phi_X = \gamma(X)
        + g\tau + \Delta\tau^2`. The j-invariant is defined by
        `\frac{g^{q+1}}{\Delta}`, `q` being the order of the base field
        of the polynomial ring. In our case, this base field is always
        finite, as we force the function ring to be of the form
        `\Fq[X]`.

        OUTPUT: an element in the base ring if the rank is two; an
        exception is raised otherwise

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.j_invariant()
            z12^10 + 4*z12^9 + 3*z12^8 + 2*z12^7 + 3*z12^6 + z12^5 + z12^3 + 4*z12^2 + z12 + 2
            sage: psi = DrinfeldModule(FqX, [p_root, 1, 1])
            sage: psi.j_invariant()
            1
            sage: rho = DrinfeldModule(FqX, [p_root, 0, 1])
            sage: rho.j_invariant()
            0

        The rank must be two::

            sage: theta = DrinfeldModule(FqX, [p_root, 1, 0])
            sage: theta.j_invariant()
            Traceback (most recent call last):
            ...
            NotImplementedError: rank must be 2
        """
        self._check_rank_two()
        g = self._gen[1]
        delta = self._gen[2]
        q = self._Fq.order()
        return (g**(q+1)) / delta

    def morphism(self):
        r"""
        Return the morphism object that defines the Drinfeld module.

        OUTPUT: a ring morphism, from the function ring to the Ore
        polynomial ring

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.morphism()
            Ring morphism:
              From: Univariate Polynomial Ring in X over Finite Field in z2 of size 5^2
              To:   Ore Polynomial Ring in t over Finite Field in z12 of size 5^12 twisted by z12 |--> z12^(5^2)
              Defn: X |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: from sage.rings.morphism import RingHomomorphism
            sage: isinstance(phi.morphism(), RingHomomorphism)
            True

        Actually, the ``DrinfeldModule`` method ``__call__`` simply
        class the ``__call__`` method of this morphism::

            sage: phi.morphism()(X) == phi(X)
            True
            sage: a = FqX.random_element()
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
            sage: phi(X) == m.im_gens()[0]
            True
            """
        return self._morphism

    def ore_polring(self):
        r"""
        Return the Ore polynomial ring of the Drinfeld module.

        OUTPUT: an Ore polynomial ring

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: ore_polring = phi.ore_polring()
            sage: ore_polring
            Ore Polynomial Ring in t over Finite Field in z12 of size 5^12 twisted by z12 |--> z12^(5^2)

        The Ore polynomial ring can also be retrieved from the category
        of the Drinfeld module::

            sage: ore_polring is phi.category().ore_polring()
            True

        The generator of the Drinfeld module is in the Ore polynomial
        ring::

            sage: phi(X) in ore_polring
            True

        The Ore variable is just the generator of the Ore polynomial
        ring::

            sage: ore_polring.gen()
            t
            sage: phi.ore_variable() is ore_polring.gen()
            True
        """
        return self._ore_polring

    def ore_variable(self):
        r"""
        Return the Ore variable.

        The Ore variable is defined as the generator of the Ore
        polynomial ring.

        OUTPUT: an Ore polynomial

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.ore_variable()
            t
            sage: phi.ore_variable() is phi.ore_polring().gen()
            True

        One can use the Ore variable to instanciate new Drinfeld
        modules...::

            sage: t = phi.ore_variable()
            sage: psi_X = phi.constant_coefficient() + 3*t + 2*t^4
            sage: psi = DrinfeldModule(FqX, psi_X)

        ...or morphisms and isogenies::

            sage: t^6 in End(phi)  # Frobenius endomorphism
            True
        """
        return self._ore_polring.gen()

    def rank(self):
        r"""
        Return the rank of the Drinfeld module.

        When the function ring is a polynomial ring, the rank is the
        degree of the generator.

        OUTPUT: an integer

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.rank()
            2
            sage: psi = DrinfeldModule(FqX, [p_root, 2])
            sage: psi.rank()
            1
            sage: rho = DrinfeldModule(FqX, [p_root, 0, 0, 0, 1])
            sage: rho.rank()
            4
        """
        return self.gen().degree()

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
                1. The degree of the characteristic divides the height
                of the input. (The height of an Ore polynomial
                `P(t)` is the maximum `n` such that `t^n` right-divides
                `P(t)`.)
                2. The input right-divides the generator, which can
                be tested with Euclidean division.

            We test if the input is an isogeny, and, if it is, we
            return the quotient of the Euclidean division.

            Height and Euclidean division of Ore polynomials are
            implemented as methods of class
            :class:`sage.rings.polynomial.ore_polynomial_element.OrePolynomial`.

            Another possible algorithm is to recursively solve a system,
            see :arxiv:`2203.06970`, eq. 1.1.

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: t = phi.ore_variable()
            sage: isog = t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4
            sage: psi = phi.velu(isog)
            sage: psi
            Drinfeld module defined by X |--> (z12^11 + 3*z12^10 + z12^9 + z12^7 + z12^5 + 4*z12^4 + 4*z12^3 + z12^2 + 1)*t^2 + (2*z12^11 + 4*z12^10 + 2*z12^8 + z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + z12^2 + z12 + 4)*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over Finite Field in z12 of size 5^12
            sage: isog in Hom(phi, psi)
            True

        This method works for endomorphisms as well::

            sage: phi.velu(phi(X)) is phi
            True
            sage: phi.velu(t^6) is phi
            True

        The following inputs do not define isogenies, and the method
        returns None::

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
        if not isog in self.ore_polring():
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
            return self.category()(quo)
