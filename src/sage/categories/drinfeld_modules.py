r"""
Category of Drinfeld modules

This module provides the class
:class:`sage.category.drinfeld_modules.DrinfeldModules`.

AUTHORS:

- Antoine Leudière (2022-04)
- Xavier Caruso (2022-06)
"""

#*****************************************************************************
#  Copyright (C) 2022      Xavier Caruso <xavier.caruso@normalesup.org>
#                          Antoine Leudière <antoine.leudiere@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import CategoryWithParameters
from sage.categories.homsets import Homsets
from sage.misc.functional import log
from sage.misc.latex import latex
from sage.rings.integer import Integer
from sage.rings.morphism import RingHomomorphism
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general


class DrinfeldModules(CategoryWithParameters):
    r"""
    This class represents the category of Drinfeld modules on a given
    `\Fq[X]`-field `K`.

    The `\Fq[X]`-field structure on `K` is given by a ring morphism
    `\gamma: \Fq[X] \to K`.

    We say that `\Fq[X]` is the *function ring of the category*; `K` is the
    *base of the category*, or simply its base ring or base field; `\Fq[X]` is
    the *function ring of the category*; *K\{\tau\}* is the *Ore
    polynomial ring of the category*;
    `t` is the *Ore variable of the category*. The *constant coefficient
    of the category* is `\gamma(X)`.

    INPUT: a ring morphism `\Fq[X] \to K`

    EXAMPLES:

    Generally, Drinfeld modules objects are created before their
    category, and the category is retrieved as an attribute of the
    Drinfeld module::

        sage: Fq = GF(11)
        sage: FqX.<X> = Fq[]
        sage: K.<z> = Fq.extension(4)
        sage: p_root = z^3 + 7*z^2 + 6*z + 10
        sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
        sage: cat = phi.category()
        sage: cat
        Category of Drinfeld modules defined by Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field of size 11
          To:   Finite Field in z of size 11^4
          Defn: X |--> z^3 + 7*z^2 + 6*z + 10

    The output tells the user that the category is only defined by the
    ring morphism `\gamma`.

    .. RUBRIC:: Properties of the category

    The defining morphism is retrieved using the method
    :meth:`morphism`::

        sage: cat.morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in X over Finite Field of size 11
          To:   Finite Field in z of size 11^4
          Defn: X |--> z^3 + 7*z^2 + 6*z + 10

    The so-called *constant coefficient* --- which is the same for all
    Drinfeld modules in the category --- is simply the image of `X` by
    this morphism:

        sage: cat.constant_coefficient()
        z^3 + 7*z^2 + 6*z + 10
        sage: cat.morphism()(X) == cat.constant_coefficient()
        True

    Similarly, the *`\Fq[X]`-characteristic* of the category is either
    `0` or the unique monic polynomial in `\Fq[X]` that generates
    `\mathrm{Ker}(\gamma)`::

        sage: cat.characteristic()
        X^2 + 7*X + 2
        sage: cat.morphism()(cat.characteristic())
        0

    Like for its Drinfeld modules, the *base* of the category is the
    field `K`::

        sage: cat.base()
        Finite Field in z of size 11^4
        sage: cat.base() is K
        True

    And the *function ring* is the polynomial ring `\Fq[X]`::

        sage: cat.function_ring()
        Univariate Polynomial Ring in X over Finite Field of size 11
        sage: cat.function_ring() is FqX
        True

    And as expected, the *Ore polynomial ring* is that of
    the Drinfeld modules in the category:

        sage: cat.ore_polring()
        Ore Polynomial Ring in t over Finite Field in z of size 11^4 twisted by z |--> z^11
        sage: cat.ore_polring() is phi.ore_polring()
        True

    .. RUBRIC:: Creating Drinfeld module objects from the category

    Calling the category with an Ore polynomial creates a Drinfeld
    module object in the category whose generator is the input::

        sage: psi = cat([p_root, 1])
        sage: psi
        Drinfeld module defined by X |--> t + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4
        sage: psi.category() is cat
        True

    Of course, the constant coefficient of the input must be the same as
    the category'::

        sage: cat([z, 1])
        Traceback (most recent call last):
        ...
        ValueError: constant coefficient must be a root of the characteristic

    It is also possible to create a random object in the category, with
    a given rank::

        sage: rho = cat.random_object(2)
        sage: rho  # random
        Drinfeld module defined by X |--> (7*z^3 + 7*z^2 + 10*z + 2)*t^2 + (9*z^3 + 5*z^2 + 2*z + 7)*t + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4
        sage: rho.rank() == 2
        True
        sage: rho.category() is cat
        True
    """

    def __init__(self, morphism, name='t'):
        gamma = morphism
        # Check input is a ring Morphism
        if not isinstance(gamma, RingHomomorphism):
            raise TypeError('category input must be a Ring morphism')
        self._morphism = morphism
        self._function_ring = gamma.domain()
        # Check domain is Fq[X]
        function_ring = self._function_ring
        if not isinstance(function_ring, PolynomialRing_general):
            raise NotImplementedError('function ring must be a polynomial ring')
        function_ring_base = function_ring.base_ring()
        if not function_ring_base.is_field() or not function_ring_base.is_finite() :
            raise TypeError('function ring base must be a finite field')
        Fq = function_ring_base
        FqX = function_ring
        X = FqX.gen()
        # Check codomain of gamma is field
        K = gamma.codomain()
        self._base = K
        if not K.is_field():
            raise TypeError('base must be a field')
        # Build K{t}
        d = log(Fq.cardinality(), Fq.characteristic())
        tau = K.frobenius_endomorphism(d)
        self._ore_polring = OrePolynomialRing(K, tau, names=name,
                polcast=False)
        # Create constant coefficient
        self._constant_coefficient = gamma(X)
        # Create characteristic
        self._characteristic = None
        if K.is_finite():
            f = gamma * FqX.coerce_map_from(Fq)  # Fq -> K
            E = K.over(f)
            self._characteristic = FqX(E(gamma(X)).minpoly())
        elif FqX.is_subring(K):
            self._characteristic = Integer(0)

    def _call_(self, gen):
        r"""
        Return a Drinfeld module object, in the category, whose
        generator is the input.

        INPUT: the generator of the Drinfeld module, given as an Ore
        polynomial or a list of coefficients

        OUTPUT: a Drinfeld module in the category

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: psi = cat([p_root, 0, 1])
            sage: psi
            Drinfeld module defined by X |--> t^2 + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4
            sage: t = phi.ore_variable()
            sage: cat(t^3 + z^3 + 7*z^2 + 6*z + 10) is phi
            True
        """
        from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule
        # If gen is not in the Ore polring, an exception is raised
        gen = self._ore_polring(gen)
        if self.characteristic()(gen[0]) != 0:
            raise ValueError('constant coefficient must be a root of the characteristic')
        return DrinfeldModule(self._function_ring, gen)

    # Somehow required for the class definition
    def _make_named_class_key(self, name):
        return self._function_ring.category()

    def _latex_(self):
        r"""
        Return a latex representation of the category

        OUTPUT: a string

        EXAMPLE:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: latex(cat)
            sage: latex(cat)
            \text{Category{ }of{ }Drinfeld{ }modules{ }defined{ }by\begin{array}{l}
            \text{\texttt{Ring{ }morphism:}}\\
            \text{\texttt{{ }{ }From:{ }Univariate{ }Polynomial{ }Ring{ }in{ }X{ }over{ }Finite{ }Field{ }of{ }size{ }11}}\\
            \text{\texttt{{ }{ }To:{ }{ }{ }Finite{ }Field{ }in{ }z{ }of{ }size{ }11{\char`\^}4}}\\
            \text{\texttt{{ }{ }Defn:{ }X{ }|{-}{-}>{ }z{\char`\^}3{ }+{ }7*z{\char`\^}2{ }+{ }6*z{ }+{ }10}}
            \end{array}
        """
        return f'\\text{{Category{{ }}of{{ }}Drinfeld{{ }}modules{{ }}' \
                f'defined{{ }}by{latex(self._morphism)}'

    def _repr_(self):
        r"""
        Return a string representation of the category

        OUTPUT: a string

        EXAMPLE:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat
            Category of Drinfeld modules defined by Ring morphism:
              From: Univariate Polynomial Ring in X over Finite Field of size 11
              To:   Finite Field in z of size 11^4
              Defn: X |--> z^3 + 7*z^2 + 6*z + 10
        """
        return f'Category of Drinfeld modules defined by {self._morphism}'

    # Somehow required for the class definition
    def Homsets(self):
        return Homsets()

    # Somehow required for the class definition
    def Endsets(self):
        return Homsets()

    def base(self):
        r"""
        Return the base of the category.

        OUTPUT: a ring

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.base()
            Finite Field in z of size 11^4
            sage: cat.base() is K
            True
        """
        return self._base

    def characteristic(self):
        r"""
        Return the `\Fq[X]`-characteristic of the category.

        OUTPUT: `0` or a monic prime polynomial in `\Fq[X]`

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.characteristic()
            X^2 + 7*X + 2

            sage: L = Frac(FqX)
            sage: psi = DrinfeldModule(FqX, [L.gen(), 1])
            sage: psi
            Drinfeld module defined by X |--> t + X over Fraction Field of Univariate Polynomial Ring in X over Finite Field of size 11
            sage: fox = psi.category()
            sage: fox.characteristic()
            0
        """
        if self._characteristic is None:
            raise NotImplementedError
        return self._characteristic 

    def constant_coefficient(self):
        r"""
        Return the constant coefficient of the category.

        OUTPUT: an element in the base

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.constant_coefficient()
            z^3 + 7*z^2 + 6*z + 10
            sage: cat.constant_coefficient() == cat.morphism()(X)
            True
        """
        return self._constant_coefficient

    def function_ring(self):
        r"""
        Return the function ring of the category.

        OUTPUT: the ring `\Fq[X]`

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.function_ring()
            Univariate Polynomial Ring in X over Finite Field of size 11
            sage: cat.function_ring() is FqX
            True
        """
        return self._function_ring

    def morphism(self):
        r"""
        Return the morphism that defines the category.

        OUTPUT: a ring morphism

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: gamma = cat.morphism()
            sage: gamma
            Ring morphism:
              From: Univariate Polynomial Ring in X over Finite Field of size 11
              To:   Finite Field in z of size 11^4
              Defn: X |--> z^3 + 7*z^2 + 6*z + 10
            sage: gamma(X) == cat.constant_coefficient()
            True
        """
        return self._morphism

    def ore_polring(self):
        r"""
        Return the Ore polynomial ring of the category.

        OUTPUT: the Ore polynomial ring `K\{\tau\}`

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.ore_polring()
            Ore Polynomial Ring in t over Finite Field in z of size 11^4 twisted by z |--> z^11
            sage: cat.ore_polring() is phi.ore_polring()
            True
        """
        return self._ore_polring

    def ore_variable(self):
        r"""
        Return the Ore variable of the category.

        OUTPUT: the generator of the Ore polynomial ring

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.ore_variable()
            t
            sage: cat.ore_variable() is phi.ore_variable()
            True
        """
        return self._ore_polring.gen()

    def random_object(self, rank):
        r"""
        Return a random Drinfeld module in the category, whose rank is
        the input.

        INPUT: an integer, the rank of the Drinfeld module

        OUTPUT: a Drinfeld module in the category

        EXAMPLES:

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(FqX, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: psi = cat.random_object(3) # random
            Drinfeld module defined by X |--> (6*z^3 + 4*z^2 + 10*z + 9)*t^3 + (4*z^3 + 8*z^2 + 8*z)*t^2 + (10*z^3 + 3*z^2 + 6*z)*t + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4
            sage: psi.rank() == 3
            True
        """
        if not isinstance(rank, Integer):
            raise TypeError('rank must be a positive integer')
        if rank <= 0:
            raise ValueError('rank must be a positive integer')

        K = self._base
        coeffs = [self._constant_coefficient]
        for _ in range(rank-1):
            coeffs.append(K.random_element())
        dom_coeff = 0
        while dom_coeff == 0:
            dom_coeff = K.random_element()
        coeffs.append(dom_coeff)

        return self(coeffs)

    # Somehow required for the class definition
    def super_categories(self):
        return []

    # Somehow required for the class definition
    class ParentMethods:

        def characteristic(self):
            return self.category().characteristic()

    # Somehow required for the class definition
    class ElementMethods:
        pass
