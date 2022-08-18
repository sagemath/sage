r"""
Drinfeld modules

This module provides the class
:class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.
For *finite* Drinfeld modules and their theory of complex multiplication, see
class
:class:`sage.rings.function_fields.drinfeld_module.finite_drinfeld_module.DrinfeldModule`.

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
    This class handles Drinfeld modules.

    Let `q` be the order of a finite field `\Fq`. Let `K` be a field
    equiped a ring morphism `\gamma: \Fq[X] \to K` --- the field `K` is
    said to be an *`\Fq[X]`-field*, and the monic polynomial that
    generates `\Ker(\gamma)` is called the *`\Fq[X]`-characteristic of
    the `\Fq[X]`-field `K`* (this `\Fq[X]`-characteristic plays the role in
    `\Fq[X]` of the standard characteristic, in `\ZZ`, of a finite
    field). Let `K\{\tau\}` be the ring of Ore polynomials with
    coefficients in `K` and Frobenius variable `\tau: x \mapsto x^q`. A
    *Drinfeld `\Fq[X]`-module over the `\Fq[X]`-field `K`* is a ring
    morphism `\phi: \Fq[X] \to K\{\tau\}` such that:

        1. The image of `\phi` has non-constant Ore polynomials.
        2. For every `a \in \Fq[X]`, the constant coefficient of the
        Ore polynomial `\phi(a)` is `\gamma(a)`.

    For `a \in \Fq[X]`, `\phi(a)` is denoted `\phi_a`.

    We say that `K` is the *base ring of `\phi`*, `\Fq[X]` is the
    *function ring of*, *K\{\tau\}* is the *Ore polynomial ring*, 
    `t` is the *Ore variable*. The *generator of `\phi`* is `\phi_X`,
    its *constant coefficient* is the constant coefficient of `\phi_X`.

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

    .. RUBRIC:: Construction

    A Drinfeld module object (class
    :class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`)
    is constructed as follows::

        sage: Fq.<z2> = GF(3^2)
        sage: FqX.<X> = Fq[]
        sage: K.<z> = Fq.extension(6)
        sage: phi = DrinfeldModule(FqX, [K.gen(), 1, 1])
        sage: phi
        Drinfeld module defined by X |--> t^2 + t + z over Finite Field in z of size 3^12

    In this example, we used a list of coefficients (``[K.gen(), 1,
    1]``) to represent the Ore polynomial `\phi_X = z + t + t^2`, `K
    = \Fq(z)`. We can also use regular Ore polynomials::

        sage: ore_polring = phi.ore_polring()
        sage: t = phi.ore_variable()  # Equals ore_polring.gen()
        sage: psi_X = K.gen() + t^3
        sage: psi = DrinfeldModule(FqX, psi_X)
        sage: psi
        Drinfeld module defined by X |--> t^3 + z over Finite Field in z of size 3^12

    .. NOTE::

        If the Ore polynomial has coefficients in the integers, the
        constructor does not try to guess if the user wants to see the
        coefficients as elements of Fq, or an extension like `K`.

    Note that ``phi`` and ``psi`` are *finite* Drinfeld modules, in the
    sense that `K` is finite. This is not mandatory::

        sage: K_infinite = Frac(FqX)
        sage: phi_infinite = DrinfeldModule(FqX, [K_infinite.gen(), 1, 1])
        sage: phi_infinite
        Drinfeld module defined by X |--> t^2 + t + X over Fraction Field of Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
        sage: phi_infinite.is_finite()
        False
        sage: phi.is_finite()
        True

    Drinfeld modules have their own category (see class
    :class:`sage.categories.drinfeld_modules.DrinfeldModules`)::

        sage: phi.category()
        Category of Drinfeld modules defined by Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z
        sage: phi.category() is psi.category()
        True
        sage: phi.category() is phi_infinite.category()
        False

    This category holds crucial information, like the
    `\Fq[X]`-characteristic of `K`::

        sage: char = phi.category().characteristic()

    .. NOTE::   

        As the output of
        :meth:`sage.rings.function_fields.drinfeld_module.finite_drinfeld_module.DrinfeldModule.category`
        suggests, the morphism `\gamma` uniquely determines the category of a Drinfeld
        module.

    .. RUBRIC:: More general `K`

    The field `K` does not need be generated by the constant coefficient
    (i.e. generated, as an extension of `\Fq`, by the image
    `\gamma(X)`). In the following example, `K` is still a
    degree six extension of `\Fq`, but `\gamma` is a projection over
    `\Fq[X]/p(X)`, with `p(X) = X^3 + (z_2 + 2)X^2 + (6*z_2 + 1)X + 3z_2
    + 5`::

        sage: p = X^2 + z2 + 2  # Prime polynomial
        sage: p_root = z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z  # Root of p
        sage: phi_inter = DrinfeldModule(FqX, [p_root, 1, 1])
        sage: phi_inter
        Drinfeld module defined by X |--> t^2 + t + z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z over Finite Field in z of size 3^12

    We can check that the morphisms `\gamma` are not the same for
    ``phi`` and ``phi_inter``, and that the `\gamma` associated to
    `\phi` is surjective, while the other one is not::

        sage: phi_inter.category().morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z
        sage: phi.category().morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z

    .. RUBRIC:: Basic methods

    For a polynomial `a \in \Fq[X]`, compute `\phi_a` by calling `phi`::

        sage: phi(X)  # phi_X
        t^2 + t + z
        sage: phi(X^3 + X + 1)  # phi_X^3 +X + 1
        t^6 + (z^11 + z^9 + 2*z^6 + 2*z^4 + 2*z + 1)*t^4 + (2*z^11 + 2*z^10 + z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^3)*t^3 + (2*z^11 + z^10 + z^9 + 2*z^7 + 2*z^6 + z^5 + z^4 + 2*z^3 + 2*z + 2)*t^2 + (2*z^11 + 2*z^8 + 2*z^6 + z^5 + z^4 + 2*z^2)*t + z^3 + z + 1
        sage: phi(1)  # phi_1
        1

    We can retrieve basic properties::

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

    We can compute the rank and height::

        sage: phi.rank()
        2
        sage: phi.height()
        1

    As well as the j-invariant if the rank is two::

        sage: phi.j()  # j-invariant
        1

    .. RUBRIC:: Morphisms, isogenies

    A *morphism of Drinfeld modules `\phi \to \psi`* is an Ore
    polynomial `f \in K\{\tau\}` such that `f \phi_a = \psi_a f` for
    every `a \in \Fq[X]`. In our case, this is equivalent to verifying
    `f \phi_X = \psi_X f`. An *isogeny* is a non-zero morphism.

    Use the ``in`` syntax to test if an Ore polynomial defines an
    isogeny::

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
    homset ``hom``::

        sage: hom = Hom(phi, phi)
        sage: frob = hom(t^6)
        sage: identity_morphism = hom(1)
        sage: zero_morphism = hom(0)
        sage: frob
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

    We can retrieve the underlying Ore polynomial with the method
    :meth:`sage.rings.function_fields.drinfeld_module.finite_drinfeld_module.DrinfeldModule.ore_polynomial`::

        sage: frob.ore_polynomial()
        t^6

    And we can easily check if a morphism defines an isogeny or an
    isomorphism (i.e. an isogeny whose underlying Ore polynomial has
    degree `0`)::

        sage: frob.is_isogeny()
        True
        sage: identity_morphism.is_isogeny()
        True
        sage: zero_morphism.is_isogeny()
        False
        sage: frob.is_isomorphism()
        False
        sage: identity_morphism.is_isomorphism()
        True
        sage: zero_morphism.is_isomorphism()
        False

    .. RUBRIC:: The Vélu formula

    Let ``ore_pol`` be a non-zero Ore polynomial ``ore_pol``. For
    Drinfeld module, it is easy to decide --- and find as the case may
    be --- if there exists a Drinfeld module ``psi``, such that
    ``ore_pol`` is an isogeny from ``self`` to ``psi``. If this
    Drinfeld module exists, it is unique.

        sage: ore_pol = (2*z^6 + z^3 + 2*z^2 + z + 2)*t + z^11 + 2*z^10 + 2*z^9 + 2*z^8 + z^7 + 2*z^6 + z^5 + z^3 + z^2 + z
        sage: psi = phi.velu(ore_pol)
        sage: psi
        Drinfeld module defined by X |--> (2*z^11 + 2*z^9 + z^6 + 2*z^5 + 2*z^4 + 2*z^2 + 1)*t^2 + (2*z^11 + 2*z^10 + 2*z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z)*t + z over Finite Field in z of size 3^12
        sage: ore_pol in Hom(phi, psi)
        True
        sage: ore_pol * phi(X) == psi(X) * ore_pol
        True

    .. RUBRIC:: Other methods

    It is possible to change the base ring::

        sage: L = K.extension(2)
        sage: phi_rebased = phi.change_ring(L)
        sage: Ltau = phi_rebased.ore_polring()
        sage: Ltau(phi(X)) == phi_rebased(X)
        True

    Given an Ore polynomial that equals `\phi_a` for some `a \in
    \Fq[X]`, we can retrieve `a` (as a morphism, a Drinfeld
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
            raise NotImplementedError('domain must be a polynomial ring')
        function_ring_base = function_ring.base_ring()
        if not function_ring_base.is_field() or not function_ring_base.is_finite() :
            raise TypeError('the base ring of the domain must be a finite field')
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
            raise TypeError('generator must be a list of coefficients '\
                    'or an Ore polynomial')

        # Build the morphism that defines the category
        if not ore_polring_base.has_coerce_map_from(function_ring.base_ring()):
            raise TypeError('base ring of function ring must coerce to base ' \
                    'ring of Ore polynomial ring')
        gamma = function_ring.hom([ore_polring_base(gen[0])])

        # Other checks in the category definition
        category = DrinfeldModules(gamma, name=name)

        # Check gen as Ore polynomial
        if ore_polring not in (None, category.ore_polring()):
            raise ValueError(f'generator must lie in {category.ore_polring()}')
        ore_polring = category.ore_polring()  # Sanity cast
        gen = ore_polring(gen)
        if gen.degree() <= 0:
            raise ValueError('generator must have positive degree')

        # Instantiate the appropriate class
        if ore_polring_base.is_finite():
            from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
            return FiniteDrinfeldModule(gen, category)
        else:
            return cls.__classcall__(cls, gen, category)

    def __init__(self, gen, category):
        CategoryObject.__init__(self, category=category)
        self._base_ring = category.base()
        self._function_ring = category.function_ring()
        self._gen = gen
        self._morphism = category._function_ring.hom([gen])
        self._ore_polring = gen.parent()
        self._Fq = self._function_ring.base_ring()  # Must be last

    ##########################
    # Special Sage functions #
    ##########################

    def _get_action_(self):
        from sage.rings.function_fields.drinfeld_modules.action import DrinfeldModuleAction
        return DrinfeldModuleAction(self)

    def _latex_(self):
        return f'\\text{{Drinfeld{{ }}module{{ }}defined{{ }}by{{ }}}} ' \
                f'{latex(self._function_ring.gen())} '\
                f'\\mapsto {latex(self._gen)}' \
                f'\\text{{{{ }}over{{ }}}}{latex(self._base_ring)}'

    def _repr_(self):
        return f'Drinfeld module defined by {self._function_ring.gen()} ' \
                f'|--> {self._gen} over {self._base_ring}'

    ###########
    # Getters #
    ###########

    def base_ring(self):
        r"""
        Return the base ring of the Drinfeld module.
    
        This is the base field of Ore polynomial ring. In particular, the base
        ring is always a field.

        A Drinfeld module is said to be finite if the base ring is
        finite.

        OUTPUT:
        
        - a field

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

        Note that in the above example, the base ring does is not
        generated by the constant coefficient (i.e. generated, as an
        extension of `\Fq`, by the image `\gamma(X)`). The base
        ring may also be the same as ``Fq``::

            sage: psi = DrinfeldModule(FqX, [Fq(1), Fq.gen()])
            sage: psi.base_ring()
            Finite Field in z2 of size 5^2
            sage: psi.base_ring() is Fq
            True

        In which case the Ore polynomial ring is isomorphic to a regular
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

    def constant_coefficient(self):
        r"""
        Return the constant coefficient of the generator (`\phi_X`).

        The `A`-characteristic of the base field (see
        :meth:`sage.categories.drinfeld_modules.DrinfeldModules.characteristic`)
        is the minimal polynomial of this constant term, over the base
        ring of the function ring. Equivalently, the constant term is
        the image, by the morphism (`\gamma`) that defines the category,
        of the generator (`X`) of the polynomial ring.
    
        OUTPUT:

        - an element in the base ring

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.constant_coefficient() == p_root
            True

        The constant coefficient is the image of ``X`` by the
        morphism that defines the category of ``phi``::

            sage: cat = phi.category()
            sage: gamma = cat.morphism()
            sage: gamma(X) == phi.constant_coefficient()
            True

        Two Drinfeld modules in the same category have the same constant
        coefficient::

            sage: t = phi.ore_variable()
            sage: psi = cat(phi.constant_coefficient() + t^3)
            sage: psi
            Drinfeld module defined by X |--> t^3 + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over Finite Field in z12 of size 5^12

        Reciprocally, it is impossible to create two Drinfeld modules in
        this category if they don't share the same constant
        coefficient::

            sage: rho = cat(phi.constant_coefficient() + 1 + t^3)
            Traceback (most recent call last):
            ...
            ValueError: constant coefficient is not a root of the characteristic

        One cal also retrieve the constant coefficient using
        ``phi(X)[0]`::

            sage: phi.constant_coefficient() == phi(X)[0]
            True
        """
        return self.gen()[0]

    def gen(self):
        r"""
        Return the generator (`\phi_X`) of the Drinfeld module.

        This method makes sense because, in our case, the function ring is
        a polynomial ring; it is generated by one element, whose image
        characterizes the morphism that defines the Drinfeld module.
            
        OUTPUT:

        - an Ore polynomial

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

    def morphism(self):
        r"""
        Return the morphism object that defines the Drinfeld module.

        OUTPUT:

        - a ring morphism, from the function ring to the Ore polynomial
          ring

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

        If the Drinfeld module is defined by a morphism `A \to
        K\{\tau\}`, this is the codomain `K\{\tau\}`.

        OUTPUT:

        - an Ore polynomial ring

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

        This is generator of the Ore polynomial ring.

        OUTPUT:

        - an Ore polynomial

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

        We can use the Ore variable to instanciate new Drinfeld
        modules...::

            sage: t = phi.ore_variable()
            sage: psi_X = phi.constant_coefficient() + 3*t + 2*t^4
            sage: psi = DrinfeldModule(FqX, psi_X)

        ...or morphisms and isogenies::

            sage: t^6 in End(phi)  # Frobenius endomorphism
            True
        """
        return self._ore_polring.gen()

    def function_ring(self):
        r"""
        Return the function ring of the Drinfeld module.

        If the Drinfeld module is defined by a morphism `A \to
        K\{\tau\}`, this is the domain `A`.

        In our case, the function ring is a polynomial ring.

        OUTPUT:

        - a polynomial ring

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

    ###########
    # Methods #
    ###########

    def __call__(self, a):
        r"""
        Return the image of ``a`` by the morphism that defines the
        Drinfeld module, i.e. `\phi_a` if the Drinfeld module is denoted
        `phi`.

        INPUT:

        - ``a`` -- an element in the function ring

        OUTPUT:

        - an element of the base ring
        """

        return self._morphism(a)

    def change_ring(self, new_field, name=None):
        r"""
        Return a Drinfeld module defined like ``self``, but with base
        ring ``new_field``.

        The new base can either be a field extension of the base ring,
        or field that has a coercion map from the field of definitions
        of the coefficients of the generator.

        INPUT:

        - ``new_field`` -- the field extension of the base ring that
          serves as base ring for the new Drinfeld module

        OUTPUT:

        - a Drinfeld module

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

        And we can check various things::

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

        We can also change the base ring to a subfield, even though some things
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

    def height(self):
        r"""
        Return the height of the Drinfeld module.

        In our case, the height is always 1.

        OUTPUT:

        - an integer

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
        Return the inverse ``ore_pol`` by the morphism that defines the
        Drinfeld module; if ``ore_pol`` is not in the image of the
        morphism, return ``None``.

        Said otherwise, return `a` if ``ore_pol`` is `phi_a`, otherwise
        return ``None``.

        INPUT:

        - ``ore_pol`` -- the Ore polynomial whose preimage we want to
          compute

        OUTPUT:

        - a polynomial

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

        When the input is not in the image of the Drinfeld module, the
        method returns None::

            sage: t = phi.ore_variable()
            sage: phi.invert(t + 1) is None
            True
            sage: phi.invert(t^3 + t^2 + 1) is None
            True

        ALGORITHM:

            The algorithm relies on the inversion of a linear algebra
            system. See [MS2019]_, 3.2.5 for details.

        TESTS::

            sage: a = FqX.random_element()
            sage: cat = phi.category()
            sage: phi_r1 = cat.random_element(1)
            sage: phi_r1.invert(phi_r1(a)) == a
            True
            sage: phi_r2 = cat.random_element(2)
            sage: phi_r2.invert(phi_r2(a)) == a
            True
            sage: phi_r3 = cat.random_element(3)
            sage: phi_r3.invert(phi_r3(a)) == a
            True
            sage: phi_r4 = cat.random_element(4)
            sage: phi_r4.invert(phi_r4(a)) == a
            True
            sage: phi_r5 = cat.random_element(5)
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
            return None

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
        Return ``True`` if the Drinfeld module is finite; return
        ``False`` otherwise.

        OUTPUT:

        - ``True`` or ``False``

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

    def j(self):
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

        OUTPUT:

        - an element in the base ring if the rank is two; an
          exception is raised otherwise

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.j()
            z12^10 + 4*z12^9 + 3*z12^8 + 2*z12^7 + 3*z12^6 + z12^5 + z12^3 + 4*z12^2 + z12 + 2
            sage: psi = DrinfeldModule(FqX, [p_root, 1, 1])
            sage: psi.j()
            1
            sage: rho = DrinfeldModule(FqX, [p_root, 0, 1])
            sage: rho.j()
            0

        The rank must be two::

            sage: theta = DrinfeldModule(FqX, [p_root, 1, 0])
            sage: theta.j()
            Traceback (most recent call last):
            ...
            NotImplementedError: the rank must be 2
        """
        self._check_rank_two()
        g = self._gen[1]
        delta = self._gen[2]
        q = self._Fq.order()
        return (g**(q+1)) / delta

    def rank(self):
        r"""
        Return the rank of the Drinfeld module.

        When the function ring is a polynomial ring, the rank is the
        degree of the generator.

        OUTPUT:

        - an integer

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
        Return a new Drinfeld module such that ``isog`` is an
        isogeny to this module with domain ``self``; if no such isogeny
        exists, return ``None``.

        If the input is zero, return ``None``, as an isogeny is
        required to be non zero.

        INPUT:

        - ``isog`` -- the Ore polynomial that defines the isogeny

        OUTPUT:

        - a Drinfeld module

        ALGORITHM:
        
            The input defines an isogeny if only if:
                1. The degree of the characteristic devides the height
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
            sage: phi.velu(t)
            sage: phi.velu(t^3 + t + 2)
        """
        if not isog in self.ore_polring():
            raise TypeError('input must be an Ore polynomial')
        if isog == 0:
            return None
        if not self.characteristic().degree().divides(isog.valuation()):
            return None
        quo, rem = (isog * self.gen()).right_quo_rem(isog)
        return None if rem != 0 else DrinfeldModule(self._function_ring, quo)

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

        OUTPUT:

        - an homset

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
             raise NotImplementedError('the rank must be 2')
