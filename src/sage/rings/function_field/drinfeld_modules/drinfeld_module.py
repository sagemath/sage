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
from sage.matrix.constructor import Matrix
from sage.misc.latex import latex
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.structure.parent import Parent
from sage.structure.sequence import Sequence
from sage.structure.unique_representation import UniqueRepresentation


class DrinfeldModule(Parent, UniqueRepresentation):
    r"""
    This class represents a Drinfeld module.

    Let `\Fq[X]` be a polynomial ring with coefficients in a finite
    field `\Fq` and let `K` be a field. We fix a ring morphism `\gamma:
    \Fq[X] \to K`, which we call the *base* of the Drinfeld module.
    Please note that the base is not a ring; in particular, it is
    not the field `K`. We also call `K` an *`\Fq[X]`-field*.

    .. NOTE::

        The base of the Drinfeld module is the base of the category of
        the Drinfeld module. 

    The monic polynomial that generates the kernel of the base is called
    the *`\Fq[X]`-characteristic of the `\Fq[X]`-field `K`*.

    Let `K\{\tau\}` be the ring of Ore polynomials with coefficients in
    `K` and Frobenius variable `\tau: x \mapsto x^q`. A *Drinfeld
    `\Fq[X]`-module over the base `\gamma`* is an `\Fq`-algebra
    morphism `\phi: \Fq[X] \to K\{\tau\}` such that:

        1. The image of `\phi` contains non-constant Ore polynomials.
        2. For every element `a` in the function ring, the constant
           coefficient `\phi(a)` is `\gamma(a)`.

    For `a` in the function ring, `\phi(a)` is denoted `\phi_a`.

    The Drinfeld module `\phi` is uniquely determined by the image
    `\phi_X` of `X`, which is an input of the class.

    We say that `\Fq[X]` is the *function ring of `\phi`*; *K\{\tau\}*
    is the *Ore polynomial ring of `\phi`*. Further, the *generator of
    `\phi`* is `\phi_X` and its *constant coefficient* is the constant
    coefficient of `\phi_X`. The `\Fq[X]`-characteristic of the
    `\Fq[X]`-field `K` can also be referred to as its *function
    ring-characteristic*. Finally, `K` is just refered to as the
    codomain base.

    Classical references on Drinfeld modules include [Gos1998]_,
    [Rosen2002]_, [VS06]_ and [Gek1998]_.

    .. NOTE::

        Drinfeld modules are defined in a larger setting, in which the
        polynomial ring `\Fq[X]` is replaced by a more general function
        ring: the ring of functions in `k` that are regular outside
        `\infty`, where `k` is a function field over `\Fq` with
        transcendence degree `1` and `\infty` is a fixed place of `k`.
        This is out of the scope of this implementation.

    INPUT:

    - ``function_ring`` -- a univariate polynomial ring whose base is a
      finite field
    - ``gen`` -- the generator of the Drinfeld module; as a list of
      coefficients or an Ore polynomial
    - ``name`` (optional) -- the name of the Ore polynomial ring gen

    .. RUBRIC:: Construction

    A Drinfeld module object is constructed as follows::

        sage: Fq.<z2> = GF(3^2)
        sage: FqX.<X> = Fq[]
        sage: K.<z> = Fq.extension(6)
        sage: phi = DrinfeldModule(FqX, [z, 1, 1])
        sage: phi
        Drinfeld module defined by X |--> t^2 + t + z over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z

    In this example, we used a list of coefficients (``[z, 1, 1]``) to
    represent the generator `\phi_X = z + t + t^2`, `K = \Fq(z)`. One can
    also use regular Ore polynomials::

        sage: ore_polring = phi.ore_polring()
        sage: t = phi.ore_polring().gen()
        sage: psi_X = z + t^3
        sage: psi = DrinfeldModule(FqX, psi_X)
        sage: psi
        Drinfeld module defined by X |--> t^3 + z over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z
        sage: psi(X) == psi_X
        True

    The generator must have positive degree::

        sage: DrinfeldModule(FqX, [z])
        Traceback (most recent call last):
        ...
        ValueError: generator must have positive degree

    The constant coefficient must be non zero::

        sage: DrinfeldModule(FqX, [K(0), K(1)])
        Traceback (most recent call last):
        ...
        ValueError: base must be a non zero morphism

    The coefficients of the generator must lie in an `\Fq[X]`-field,
    where `\Fq[X]` is the function ring of the Drinfeld module::

        sage: DrinfeldModule(FqX, [z, QQ(1)])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce into base codomain

        sage: DrinfeldModule(FqX, [1, QQ(1)])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce into base codomain

    The function ring must be an univariate polynomial ring whose
    base is a finite field::

        sage: DrinfeldModule(K, [z, 1, 1])
        Traceback (most recent call last):
        ...
        NotImplementedError: function ring must be a polynomial ring

        sage: FqXY.<Y> = FqX[]
        sage: DrinfeldModule(FqXY, [z, 1, 1])
        Traceback (most recent call last):
        ...
        TypeError: function ring base must be a finite field

    If you already defined a category of Drinfeld modules, and you
    create a Drinfeld module through this category, you must
    ensure that the constant coefficient is the generator of the algebra
    morphism that defines the category::

        sage: cat = phi.category()
        sage: cat.object([1, 1, K(1)])
        Traceback (most recent call last):
        ...
        ValueError: constant coefficient must be the generator of the morphism that defines the category

    .. NOTE::

        The reader may think that it is odd to build a Drinfeld module
        without explicitly specifying the base. However, the base can be
        deduced from the generator, and we omit the base in the input
        of the class for conciseness.

    .. RUBRIC:: Possible bases

    The base does not need be surjective like in the above examples. In
    the following example, the base codomain is still a degree six
    extension of `\Fq`, but the base is a projection over a degree two
    extension with modulus `X^3 + (z_2 + 2)X^2 + (6*z_2 + 1)X + 3z_2 +
    5`::

        sage: p = X^2 + z2 + 2
        sage: p_root = z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z
        sage: rho = DrinfeldModule(FqX, [p_root, 1, 1])
        sage: rho
        Drinfeld module defined by X |--> t^2 + t + z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z

    Drinfeld modules `\phi` and `\rho` have different based. That of
    `\phi` is surjective while that of `\rho` is note::

        sage: rho.category().base()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z^10 + 2*z^9 + z^8 + z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z
        sage: phi.category().base()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z

    Note that ``phi`` and ``psi`` are *finite* Drinfeld modules, in the
    sense that `K` is finite. But `K` can be infinite::

        sage: sigma = DrinfeldModule(FqX, [Frac(FqX).gen(), 1, 1])
        sage: sigma
        Drinfeld module defined by X |--> t^2 + t + X over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Fraction Field of Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          Defn: X |--> X
        sage: sigma.is_finite()
        False
        sage: phi.is_finite()
        True

    .. RUBRIC:: The category of Drinfeld modules

    Drinfeld modules have their own category (see class
    :class:`sage.categories.drinfeld_modules.DrinfeldModules`)::

        sage: phi.category()
        Category of Drinfeld modules defined over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z
        sage: phi.category() is psi.category()
        True
        sage: phi.category() is sigma.category()
        False

    This category holds crucial information, like the
    function ring-characteristic of the base::

        sage: char = phi.category().characteristic()

    As the output of :meth:`category` suggests, the base uniquely
    determines the category.

    .. RUBRIC:: Basics

    Images under the Drinfeld module are computed by calling the object::

        sage: phi(X)  # phi_X
        t^2 + t + z
        sage: phi(X^3 + X + 1)  # phi_X^3 +X + 1
        t^6 + (z^11 + z^9 + 2*z^6 + 2*z^4 + 2*z + 1)*t^4 + (2*z^11 + 2*z^10 + z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^3)*t^3 + (2*z^11 + z^10 + z^9 + 2*z^7 + 2*z^6 + z^5 + z^4 + 2*z^3 + 2*z + 2)*t^2 + (2*z^11 + 2*z^8 + 2*z^6 + z^5 + z^4 + 2*z^2)*t + z^3 + z + 1
        sage: phi(1)  # phi_1
        1

    This is useful to quickly retrieve the generator of the Drinfeld
    module. Furthermore, a Drinfeld `\Fq[X]`-module can be seen as an
    Ore polynomial with positive degree and constant coefficient
    `\gamma(X)`, where `\gamma` is the base. This analogy is the
    motivation for the following methods::

        sage: phi.coefficients()
        [z, 1, 1]

        sage: phi.coefficient(1)
        1
        sage: phi.coefficient(1)
        1

    One can retrieve basic properties::

        sage: phi.base()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z

        sage: phi.ore_polring()  # K{t}
        Ore Polynomial Ring in t over Finite Field in z of size 3^12 twisted by z |--> z^(3^2)

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
    every `a` in the function ring. In our case, this is equivalent to
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
    homset (``hom``)::

        sage: hom = Hom(phi, phi)
        sage: frobenius_endomorphism = hom(t^6)
        sage: identity_morphism = hom(1)
        sage: zero_morphism = hom(0)
        sage: frobenius_endomorphism
        Drinfeld Module morphism:
          From (gen): t^2 + t + z
          To (gen):   t^2 + t + z
          Defn:       t^6
        sage: identity_morphism
        Drinfeld Module morphism:
          From (gen): t^2 + t + z
          To (gen):   t^2 + t + z
          Defn:       1
        sage: zero_morphism
        Drinfeld Module morphism:
          From (gen): t^2 + t + z
          To (gen):   t^2 + t + z
          Defn:       0

    The underlying Ore polynomial is retrieved with the method
    :meth:`ore_polynomial`::

        sage: frobenius_endomorphism.ore_polynomial()
        t^6
        sage: identity_morphism.ore_polynomial()
        1

    It is easy to check if a morphism is an isogeny, endomorphism or
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

    Let ``ore_pol`` be a non-zero Ore polynomial. We can decide if there
    exists a Drinfeld module ``psi`` such that ``ore_pol`` is an isogeny
    from ``self`` to ``psi``. If so, we find ``psi``::

        sage: ore_pol = (2*z^6 + z^3 + 2*z^2 + z + 2)*t + z^11 + 2*z^10 + 2*z^9 + 2*z^8 + z^7 + 2*z^6 + z^5 + z^3 + z^2 + z
        sage: psi = phi.velu(ore_pol)
        sage: psi
        Drinfeld module defined by X |--> (2*z^11 + 2*z^9 + z^6 + 2*z^5 + 2*z^4 + 2*z^2 + 1)*t^2 + (2*z^11 + 2*z^10 + 2*z^9 + z^8 + 2*z^7 + 2*z^6 + z^5 + 2*z^4 + 2*z^2 + 2*z)*t + z over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z
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

    The `\Fq[X]`-Drinfeld module `\phi` induces a special left
    `\Fq[X]`-module structure on any field extension `L/K`. Let `x \in
    L` and `a` be in the function ring; the action is defined as `(a,
    x) \mapsto \phi_a(x)`. The method :meth:`action` returns an
    ``Action`` object representing the Drinfeld module action.

    .. NOTE::

        In this implementation, `L` is `L`.

        sage: action = phi.action()
        sage: action
        Action on Finite Field in z of size 3^12 induced by Drinfeld module defined by X |--> t^2 + t + z over base Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field in z2 of size 3^2
          To:   Finite Field in z of size 3^12
          Defn: X |--> z

    The action on elements is computed by calling the action object::

        sage: P = X + 1
        sage: a = z
        sage: action(P, a)
        ...
        z^9 + 2*z^8 + 2*z^7 + 2*z^6 + 2*z^3 + z^2
        sage: action(0, K.random_element())
        0
        sage: action(FqX.random_element(), 0)
        0

    .. RUBRIC:: Inverting the Drinfeld module

    Given an Ore polynomial that equals `\phi_a` for some function ring
    elelement `a`, one can retrieve `a` (as a morphism, a Drinfeld
    module is injective, see [Gos1998]_, cor. 4.5.2.)::

        sage: a = FqX.random_element()
        sage: phi.invert(phi(a)) == a
        True

    TESTS:

        sage: Fq = K = GF(2)
        sage: FqX.<X> = Fq[]
        sage: phi = DrinfeldModule(FqX, [1, 1])
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce into base codomain

        sage: Fq = K = GF(2)
        sage: FqX.<X> = Fq[]
        sage: phi = DrinfeldModule(FqX, [K(1), 1])
        sage: isinstance(phi.ore_polring(), OrePolynomialRing)
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
            ore_polring_base = ore_polring.base_ring()
            name = ore_polring.variable_name()
        # `gen` is a list of coefficients (function_ring = Fq[X]):
        elif isinstance(gen, (list, tuple)):
            ore_polring = None
            ore_polring_base = Sequence(gen).universe()
        else:
            raise TypeError('generator must be list of coefficients or Ore '
                            'polynomial')
        # The coefficients are in a base ring that has coercion from Fq:
        if not (hasattr(ore_polring_base, 'has_coerce_map_from')
                and ore_polring_base.has_coerce_map_from(
                        function_ring.base_ring())):
            raise ValueError('function ring base must coerce into base codomain')

        # Build the morphism that defines the category
        base = function_ring.hom([ore_polring_base(gen[0])])

        # Other checks in the category definition
        category = DrinfeldModules(base, name=name)

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
        self._base = category.base()
        self._function_ring = category.function_ring()
        self._gen = gen
        self._morphism = category._function_ring.hom([gen])
        self._ore_polring = gen.parent()
        self._Fq = self._function_ring.base_ring()  # Must be last
        super().__init__(base=self._base, category=category)

    def __call__(self, a):
        r"""
        Return the image of ``a`` by the morphism that defines the
        Drinfeld module; i.e. `\phi_a` if the Drinfeld module is denoted
        `phi`.

        INPUT:

        - ``a`` -- a function ring element

        OUTPUT: an element in the base codomain

        TESTS:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])

            sage: a = X^3 + 4*X + 2
            sage: phi(a) == phi(X)^3 + 4*phi(X) + 2
            True
            sage: phi(a)[0] == p_root^3 + 4*p_root + 2
            True

            sage: phi(0)
            0
            sage: phi(1)
            1
            sage: phi(X) == phi._gen
            True

            sage: a = FqX.random_element(5)
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
          morphisms, usually `self.category()`

        OUTPUT: an homset

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
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
            \text{Drinfeld{ }module{ }defined{ }by{ }} X \mapsto z_{12}^{5} t^{2} + z_{12}^{3} t + 2 z_{12}^{11} + 2 z_{12}^{10} + z_{12}^{9} + 3 z_{12}^{8} + z_{12}^{7} + 2 z_{12}^{5} + 2 z_{12}^{4} + 3 z_{12}^{3} + z_{12}^{2} + 2 z_{12}\text{{ }over{ }base{ }}\begin{array}{l}
            \text{\texttt{Ring{ }morphism:}}\\
            \text{\texttt{{ }{ }From:{ }Univariate{ }Polynomial{ }Ring{ }in{ }X{ }over{ }Finite{ }Field{ }in{ }z2{ }of{ }size{ }5{\char`\^}2}}\\
            \text{\texttt{{ }{ }To:{ }{ }{ }Finite{ }Field{ }in{ }z12{ }of{ }size{ }5{\char`\^}12}}\\
            \text{\texttt{{ }{ }Defn:{ }X{ }|{-}{-}>{ }2*z12{\char`\^}11{ }+{ }2*z12{\char`\^}10{ }+{ }z12{\char`\^}9{ }+{ }3*z12{\char`\^}8{ }+{ }z12{\char`\^}7{ }+{ }2*z12{\char`\^}5{ }+{ }2*z12{\char`\^}4{ }+{ }3*z12{\char`\^}3{ }+{ }z12{\char`\^}2{ }+{ }2*z12}}
            \end{array}
        """
        return f'\\text{{Drinfeld{{ }}module{{ }}defined{{ }}by{{ }}}} ' \
               f'{latex(self._function_ring.gen())} '\
               f'\\mapsto {latex(self._gen)}' \
               f'\\text{{{{ }}over{{ }}base{{ }}}}{latex(self._base)}'

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
            Drinfeld module defined by X |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over base Ring morphism:
              From: Univariate Polynomial Ring in X over Finite Field in z2 of size 5^2
              To:   Finite Field in z12 of size 5^12
              Defn: X |--> 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
        """
        return f'Drinfeld module defined by {self._function_ring.gen()} ' \
               f'|--> {self._gen} over base {self._base}'

    def action(self):
        r"""
        Return the action object
        (:class:`sage.rings.function_field.drinfeld_modules.action.Action`)
        that represents the module action, on the base codomain, that is
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
            Action on Finite Field in z12 of size 5^12 induced by Drinfeld module defined by X |--> z12^5*t^2 + z12^3*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over base Ring morphism:
              From: Univariate Polynomial Ring in X over Finite Field in z2 of size 5^2
              To:   Finite Field in z12 of size 5^12
              Defn: X |--> 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12

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

    def coefficient(self, n):
        r"""
        Return the n-th coefficient of the generator.

        INPUT:

        - ``n`` -- a non-negative integer

        OUTPUT: an element in the base codomain

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
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
        non-zero coefficients; otherwise, return all of them.

        INPUT:

        - ``sparse`` -- a boolean

        OUTPUT: a list of elements in the base codomain

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.coefficients()
            [2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12,
             z12^3,
             z12^5]

        Careful, the method only returns the non-zero coefficients,
        unless otherwise specified::

            sage: rho = DrinfeldModule(FqX, [p_root, 0, 0, 0, 1])
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
        Return the height of the Drinfeld module if the function field
        characteristic is a prime ideal; raise ValueError otherwise.

        The height of a Drinfeld module is defined when the function
        field characteristic is a prime ideal. In our case, this ideal
        is even generated by a monic polynomial `\mathfrak{p}` in the
        function field. Write `\phi_\mathfrak{p} = a_s \tau^s + \dots +
        \tau^{r*\deg(\mathfrak{p})}`. The *height* of the Drinfeld
        module is the well-defined positive integer `h =
        \frac{s}{\deg(\mathfrak{p})}`.

        .. NOTE::

            See [Gos1998]_, Definition 4.5.8 for the general definition.

        A rank two Drinfeld module is supersingular if and only if its
        height equals its rank.

        OUTPUT: an integer

        EXAMPLES:

            sage: Fq = GF(25)
            sage: FqX.<X> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(FqX, [p_root, z12^3, z12^5])
            sage: phi.height() == 1
            True
            sage: phi.is_ordinary()
            True

            sage: L = Frac(FqX)
            sage: phi = DrinfeldModule(FqX, [L(2), L(1)])
            sage: phi.height()
            Traceback (most recent call last):
            ...
            ValueError: height is defined for prime function field characteristic

            sage: Fq = GF(343)
            sage: FqX.<X> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(FqX, [1, 0, z6])
            sage: phi.height()
            2
            sage: phi.is_supersingular()
            True

        """
        try:
            if self.characteristic().is_zero():
                raise ValueError('height is defined for prime ' \
                                 'function field characteristic')
            else:
                p = self.characteristic()
                return Integer((self(p).valuation()) // (p.degree()))
        except NotImplementedError:
            raise NotImplementedError('height not implemented in this case')

    def invert(self, ore_pol):
        r"""
        Return the preimage of the input under the Drinfeld module;
        raise an exception if the input is not in the image of the
        Drinfeld module.

        INPUT:

        - ``ore_pol`` -- the Ore polynomial whose preimage we want to
          compute

        OUTPUT: a function ring element

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
        if ore_pol not in self._ore_polring:
            raise TypeError('input must be an Ore polynomial')
        if ore_pol in self._base.codomain():
            return self._Fq(ore_pol)
        if deg % r != 0:
            raise ValueError('input must be in the image of the Drinfeld '
                             'module')

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
        Return the j-invariant of the Drinfeld module if the rank is
        two; raise a NotImplementedError otherwise.

        Assume the rank is two. Write the generator `\phi_X = \omega +
        g\tau + \Delta\tau^2`. The j-invariant is defined by
        `\frac{g^{q+1}}{\Delta}`, `q` being the order of the base ring
        of the function ring. In our case, this field is always finite.

        OUTPUT: an element in the base codomain

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
        g = self.coefficient(1)
        delta = self.coefficient(2)
        q = self._Fq.order()
        return (g**(q+1)) / delta

    def morphism(self):
        r"""
        Return the morphism object that defines the Drinfeld module.

        OUTPUT: a ring morphism from the function ring to the Ore
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

    def rank(self):
        r"""
        Return the rank of the Drinfeld module.

        In our case, the rank is the degree of the generator.

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
            sage: t = phi.ore_polring().gen()
            sage: isog = t + 2*z12^11 + 4*z12^9 + 2*z12^8 + 2*z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + 4*z12^2 + 4*z12 + 4
            sage: psi = phi.velu(isog)
            sage: psi
            Drinfeld module defined by X |--> (z12^11 + 3*z12^10 + z12^9 + z12^7 + z12^5 + 4*z12^4 + 4*z12^3 + z12^2 + 1)*t^2 + (2*z12^11 + 4*z12^10 + 2*z12^8 + z12^6 + 3*z12^5 + z12^4 + 2*z12^3 + z12^2 + z12 + 4)*t + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over base Ring morphism:
              From: Univariate Polynomial Ring in X over Finite Field in z2 of size 5^2
              To:   Finite Field in z12 of size 5^12
              Defn: X |--> 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
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
