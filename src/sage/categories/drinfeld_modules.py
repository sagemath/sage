r"""
Drinfeld modules over a base

This module provides the class
:class:`sage.category.drinfeld_modules.DrinfeldModules`.

AUTHORS:

- Antoine Leudière (2022-04)
- Xavier Caruso (2022-06)
"""

# *****************************************************************************
#   Copyright (C) 2022      Xavier Caruso <xavier.caruso@normalesup.org>
#                           Antoine Leudière <antoine.leudiere@inria.fr>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#                   http://www.gnu.org/licenses/
# ******************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.homsets import Homsets
from sage.misc.functional import log
from sage.misc.latex import latex
from sage.rings.integer import Integer
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.ring_extension import RingExtension_generic


class DrinfeldModules(Category_over_base_ring):
    r"""
    This class implements the category of Drinfeld
    `\mathbb{F}_q[T]`-modules on a given base field.

    Let `\mathbb{F}_q[T]` be a polynomial ring with coefficients in a
    finite field `\mathbb{F}_q` and let `K` be a field. Fix a ring
    morphism `\gamma: \mathbb{F}_q[T] \to K`; we say that `K` is an
    `\mathbb{F}_q[T]`*-field*. Let `K\{\tau\}` be the ring of Ore
    polynomials with coefficients in `K`, whose multiplication is given
    by the rule `\tau \lambda = \lambda^q \tau` for any `\lambda \in K`.

    We call `K` the *base field* of the category, and `\gamma` its *base
    morphism*.

    .. NOTE::

        Equivalently, the base of the category could be defined as the
        base morphism `\gamma: \mathbb{F}_q[T] \to K`.

    The monic polynomial that generates the kernel of `\gamma` is called
    the `\mathbb{F}_q[T]`-*characteristic*, or *function-field
    characteristic*, of the base field. We say that `\mathbb{F}_q[T]` is
    the *function ring* of the category; `K\{\tau\}` is the *Ore
    polynomial ring*. The constant coefficient of the category is the
    image of `T` under the base morphism.

    INPUT: the base field

    .. RUBRIC:: Construction

    Generally, Drinfeld modules objects are created before their
    category, and the category is retrieved as an attribute of the
    Drinfeld module::

        sage: Fq = GF(11)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(4)
        sage: p_root = z^3 + 7*z^2 + 6*z + 10
        sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
        sage: cat = phi.category()
        sage: cat
        Category of Drinfeld modules defined over Finite Field in z of size 11^4 over its base

    The output tells the user that the category is only defined by its
    base.

    .. RUBRIC:: Properties of the category

    The base field is retrieved using the method :meth:`base`.

        sage: cat.base()
        Finite Field in z of size 11^4 over its base

    Equivalently, one can use :meth:`base_morphism` to retrieve the base
    morphism::

        sage: cat.base_morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in T over Finite Field of size 11
          To:   Finite Field in z of size 11^4 over its base
          Defn: T |--> z^3 + 7*z^2 + 6*z + 10

    The so-called constant coefficient --- which is the same for all
    Drinfeld modules in the category --- is simply the image of `T` by
    the base morphism::

        sage: cat.constant_coefficient()
        z^3 + 7*z^2 + 6*z + 10
        sage: cat.base_morphism()(T) == cat.constant_coefficient()
        True

    Similarly, the function ring-characteristic of the category is
    either `0` or the unique monic polynomial in `\mathbb{F}_q[T]` that
    generates the kernel of the base::

        sage: cat.characteristic()
        T^2 + 7*T + 2
        sage: cat.base_morphism()(cat.characteristic())
        0

    The base field, base morphism, function ring and Ore polynomial ring
    are the same for the category and its objects::

        sage: cat.base() is phi.base()
        True
        sage: cat.base_morphism() is phi.base_morphism()
        True
        sage: cat.function_ring() is phi.function_ring()
        True
        sage: cat.function_ring()
        Univariate Polynomial Ring in T over Finite Field of size 11
        sage: cat.ore_polring() is phi.ore_polring()
        True
        sage: cat.ore_polring()
        Ore Polynomial Ring in t over Finite Field in z of size 11^4 over its base twisted by Frob

    .. RUBRIC:: Creating Drinfeld module objects from the category

    Calling :meth:`object` with an Ore polynomial creates a Drinfeld module
    object in the category whose generator is the input::

        sage: psi = cat.object([p_root, 1])
        sage: psi
        Drinfeld module defined by T |--> t + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4 over its base
        sage: psi.category() is cat
        True

    Of course, the constant coefficient of the input must be the same as
    the category::

        sage: cat.object([z, 1])
        Traceback (most recent call last):
        ...
        ValueError: constant coefficient must equal that of the category

    It is also possible to create a random object in the category. The
    input is the desired rank::

        sage: rho = cat.random_object(2)
        sage: rho  # random
        Drinfeld module defined by T |--> (7*z^3 + 7*z^2 + 10*z + 2)*t^2 + (9*z^3 + 5*z^2 + 2*z + 7)*t + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4
        sage: rho.rank() == 2
        True
        sage: rho.category() is cat
        True

    TESTS::

        sage: Fq = GF(11)
        sage: A.<T> = Fq[]
        sage: K.<z> = Fq.extension(4)
        sage: from sage.categories.drinfeld_modules import DrinfeldModules
        sage: base = Hom(A, K)(0)
        sage: cat = DrinfeldModules(base)
        Traceback (most recent call last):
        ...
        TypeError: base field must be a ring extension
    
    ::

        sage: cat.base().defining_morphism() == cat.base_morphism()
        True

    ::

        sage: base = Hom(A, A)(1)
        sage: cat = DrinfeldModules(base)
        Traceback (most recent call last):
        ...
        TypeError: base field must be a ring extension

    ::

        sage: base = 'I hate Rostropovitch'
        sage: cat = DrinfeldModules(base)  # known bug (blankline)
        <BLANKLINE>
        Traceback (most recent call last):
        ...
        TypeError: input must be a ring morphism

    ::

        sage: ZZT.<T> = ZZ[]
        sage: base = Hom(ZZT, K)(1)
        sage: cat = DrinfeldModules(base)  # known bug (blankline)
        <BLANKLINE>
        Traceback (most recent call last):
        ...
        TypeError: function ring base must be a finite field
    """

    def __init__(self, base_field, name='t'):
        r"""
        Initialize `self`.

        INPUT:

        - ``base_field`` -- the base field, which is a ring extension
          over a base
        - ``name`` (default: `'t'`) -- the name of the Ore polynomial
          variable

        TESTS::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: ore_polring.<t> = OrePolynomialRing(phi.base(), phi.base().frobenius_endomorphism())
            sage: cat._ore_polring is ore_polring
            True
            sage: i = phi.base().coerce_map_from(K)
            sage: base_morphism = Hom(A, K)(p_root)
            sage: cat.base() == K.over(base_morphism)
            True
            sage: cat._base_morphism == i * base_morphism
            True
            sage: cat._function_ring is A
            True
            sage: cat._constant_coefficient == base_morphism(T)
            True
            sage: cat._characteristic(cat._constant_coefficient)
            0
        """
        # Check input is a ring extension
        if not isinstance(base_field, RingExtension_generic):
            raise TypeError('base field must be a ring extension')
        base_morphism = base_field.defining_morphism()
        self._base_morphism = base_morphism
        # Check input is a field
        if not base_field.is_field():
            raise TypeError('input must be a field')
        self._base_field = base_field
        self._function_ring = base_morphism.domain()
        # Check domain of base morphism is Fq[T]
        function_ring = self._function_ring
        if not isinstance(function_ring, PolynomialRing_general):
            raise NotImplementedError('function ring must be a polynomial '
                                      'ring')
        function_ring_base = function_ring.base_ring()
        if not function_ring_base.is_field() \
                or not function_ring_base.is_finite():
            raise TypeError('function ring base must be a finite field')
        # Shortcuts
        Fq = function_ring_base
        A = function_ring
        T = A.gen()
        K = base_field  # A ring extension
        # Build K{t}
        d = log(Fq.cardinality(), Fq.characteristic())
        tau = K.frobenius_endomorphism(d)
        self._ore_polring = OrePolynomialRing(K, tau, names=name,
                                              polcast=False)
        # Create constant coefficient
        self._constant_coefficient = base_morphism(T)
        # Create characteristic
        self._characteristic = None
        if K.is_finite():
            self._characteristic = A(K.over(Fq)(base_morphism(T)).minpoly())
        try:
            if A.is_subring(K):
                self._characteristic = Integer(0)
        except NotImplementedError:
            pass
        super().__init__(base=base_field)

    def _latex_(self):
        r"""
        Return a latex representation of the category.

        OUTPUT: a string

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: latex(cat)
            \text{Category{ }of{ }Drinfeld{ }modules{ }defined{ }over{ }\Bold{F}_{11^{4}}
        """
        return f'\\text{{Category{{ }}of{{ }}Drinfeld{{ }}modules{{ }}' \
               f'defined{{ }}over{{ }}{latex(self._base_field)}'

    def _repr_(self):
        r"""
        Return a string representation of the category.

        OUTPUT: a string

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat
            Category of Drinfeld modules defined over Finite Field in z of size 11^4 over its base
        """
        return f'Category of Drinfeld modules defined over {self._base_field}'

    def Homsets(self):
        r"""
        Return the category of homsets.

        OUTPUT: the category of homsets

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: from sage.categories.homsets import Homsets
            sage: cat.Homsets() is Homsets()
            True
        """
        return Homsets()

    def Endsets(self):
        r"""
        Return the category of endsets.

        OUTPUT: the category of endsets

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: from sage.categories.homsets import Homsets
            sage: cat.Endsets() is Homsets().Endsets()
            True
        """
        return Homsets().Endsets()

    def base_morphism(self):
        r"""
        Return the base morphism of the category.

        OUTPUT: a ring morphism

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.base_morphism()
            Ring morphism:
              From: Univariate Polynomial Ring in T over Finite Field of size 11
              To:   Finite Field in z of size 11^4 over its base
              Defn: T |--> z^3 + 7*z^2 + 6*z + 10
            sage: cat.constant_coefficient() == cat.base_morphism()(T)
            True
        """
        return self._base_morphism

    def characteristic(self):
        r"""
        Return the function ring-characteristic.

        OUTPUT: `0` or a monic prime polynomial in the function ring

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.characteristic()
            T^2 + 7*T + 2

        ::

            sage: psi = DrinfeldModule(A, [Frac(A).gen(), 1])  # todo: not tested
            sage: psi.category().characteristic()  # todo: not tested
            0
        """
        if self._characteristic is None:
            raise NotImplementedError('function ring characteristic not' \
                                      'implemented in this case')
        return self._characteristic

    def constant_coefficient(self):
        r"""
        Return the constant coefficient of the category.

        OUTPUT: an element in the base field

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.constant_coefficient()
            z^3 + 7*z^2 + 6*z + 10
            sage: cat.constant_coefficient() == cat.base()(T)
            True
        """
        return self._constant_coefficient

    def function_ring(self):
        r"""
        Return the function ring of the category.

        OUTPUT: a univariate polynomial ring

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.function_ring()
            Univariate Polynomial Ring in T over Finite Field of size 11
            sage: cat.function_ring() is A
            True
        """
        return self._function_ring

    def object(self, gen):
        r"""
        Return a Drinfeld module object in the category whose generator
        is the input.

        INPUT: the generator of the Drinfeld module, given as an Ore
        polynomial or a list of coefficients

        OUTPUT: a Drinfeld module in the category

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: psi = cat.object([p_root, 0, 1])
            sage: psi
            Drinfeld module defined by T |--> t^2 + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4 over its base
            sage: t = phi.ore_polring().gen()
            sage: cat.object(t^3 + z^3 + 7*z^2 + 6*z + 10) is phi
            True
        """
        from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule
        # If gen is not in the Ore polring, an exception is raised
        gen = self._ore_polring(gen)
        T = self._function_ring.gen()
        if gen[0] != self._base_morphism(T):
            raise ValueError('constant coefficient must equal that of the ' \
                             'category')
        return DrinfeldModule(self._function_ring, gen)

    def ore_polring(self):
        r"""
        Return the Ore polynomial ring of the category.

        OUTPUT: an Ore polynomial ring

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.ore_polring()
            Ore Polynomial Ring in t over Finite Field in z of size 11^4 over its base twisted by Frob
            sage: cat.ore_polring() is phi.ore_polring()
            True
        """
        return self._ore_polring

    def random_object(self, rank):
        r"""
        Return a random Drinfeld module in the category, whose rank is
        the input.

        INPUT: an integer, the rank of the Drinfeld module

        OUTPUT: a Drinfeld module in the category

        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: psi = cat.random_object(3) # random
            Drinfeld module defined by T |--> (6*z^3 + 4*z^2 + 10*z + 9)*t^3 + (4*z^3 + 8*z^2 + 8*z)*t^2 + (10*z^3 + 3*z^2 + 6*z)*t + z^3 + 7*z^2 + 6*z + 10 over Finite Field in z of size 11^4
            sage: psi.rank() == 3
            True
        """
        if not isinstance(rank, Integer):
            raise TypeError('rank must be a positive integer')
        if rank <= 0:
            raise ValueError('rank must be a positive integer')

        K = self._base_field
        coeffs = [self._constant_coefficient]
        for _ in range(rank-1):
            coeffs.append(K.random_element())
        dom_coeff = 0
        while dom_coeff == 0:
            dom_coeff = K.random_element()
        coeffs.append(dom_coeff)

        return self.object(coeffs)

    def super_categories(self):
        """
        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: cat = phi.category()
            sage: cat.super_categories()
            []
        """
        return []

    class ParentMethods:

        def base(self):
            r"""
            Return the base field of the Drinfeld module.

            OUTPUT: a field, which is a ring extension over a base

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z12> = Fq.extension(6)
                sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
                sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
                sage: phi.base()
                Finite Field in z12 of size 5^12 over its base

            The base can be infinite::

                sage: sigma = DrinfeldModule(A, [Frac(A).gen(), 1])  # todo: not tested
                sage: sigma.base()  # todo: not tested
                Fraction Field of Univariate Polynomial Ring in T over Finite Field in z2 of size 5^2 over its base
            """
            return self.category().base()

        def base_morphism(self):
            r"""
            Return the base morphism of the Drinfeld module.

            OUTPUT: a ring morphism

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z12> = Fq.extension(6)
                sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
                sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
                sage: phi.base_morphism()
                Ring morphism:
                  From: Univariate Polynomial Ring in T over Finite Field in z2 of size 5^2
                  To:   Finite Field in z12 of size 5^12 over its base
                  Defn: T |--> 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12

            The base field can be infinite::

                sage: sigma = DrinfeldModule(A, [Frac(A).gen(), 1])  # todo: not tested
                sage: sigma.base_morphism()  # todo: not tested
            """
            return self.category().base_morphism()

        def characteristic(self):
            r"""
            Return the function ring-characteristic.

            OUTPUT: a univariate polynomial ring

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z12> = Fq.extension(6)
                sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
                sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
                sage: phi.characteristic()  # todo: not tested
                T^2 + (4*z2 + 2)*T + 2
                sage: phi.base_morphism()(phi.characteristic())
                0

            ::

                sage: L = Frac(A)  # todo: not tested
                sage: psi = DrinfeldModule(A, [L(1), 0, 0, L(1)])  # todo: not tested
                sage: psi.characteristic()  # todo: not tested
                0
            """
            return self.category().characteristic()

        def function_ring(self):
            r"""
            Return the function ring of the Drinfeld module.

            OUTPUT: a univariate polynomial ring

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z12> = Fq.extension(6)
                sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
                sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
                sage: phi.function_ring() is A
                True
            """
            return self.category().function_ring()

        def constant_coefficient(self):
            r"""
            Return the constant coefficient of the generator.

            OUTPUT: an element in the base field

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z12> = Fq.extension(6)
                sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
                sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
                sage: phi.constant_coefficient() == p_root
                True

            Let `\mathbb{F}_q[T]` be the function ring, and let `\gamma`
            the base of the Drinfeld module. The constant coefficient
            equals `\gamma(T)`::

                sage: cat = phi.category()
                sage: base = cat.base()
                sage: base(T) == phi.constant_coefficient()
                True

            Naturally, two Drinfeld modules in the same category have the
            same constant coefficient::

                sage: t = phi.ore_polring().gen()
                sage: psi = cat.object(phi.constant_coefficient() + t^3)
                sage: psi
                Drinfeld module defined by T |--> t^3 + 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12 over Finite Field in z12 of size 5^12 over its base

            Reciprocally, it is impossible to create two Drinfeld modules in
            this category if they do not share the same constant
            coefficient::

                sage: rho = cat.object(phi.constant_coefficient() + 1 + t^3)
                Traceback (most recent call last):
                ...
                ValueError: constant coefficient must equal that of the category
            """
            return self.category().constant_coefficient()

        def ore_polring(self):
            r"""
            Return the Ore polynomial ring of the Drinfeld module.

            OUTPUT: an Ore polynomial ring

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z12> = Fq.extension(6)
                sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
                sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
                sage: ore_polring = phi.ore_polring()
                sage: ore_polring
                Ore Polynomial Ring in t over Finite Field in z12 of size 5^12 over its base twisted by Frob^2

            The Ore polynomial ring can also be retrieved from the category
            of the Drinfeld module::

                sage: ore_polring is phi.category().ore_polring()
                True

            The generator of the Drinfeld module is in the Ore polynomial
            ring::

                sage: phi(T) in ore_polring
                True
            """
            return self.category().ore_polring()
