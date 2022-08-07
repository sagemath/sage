r"""
Drinfeld modules

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

from sage.misc.latex import latex
from sage.structure.sequence import Sequence
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector

class DrinfeldModule(UniqueRepresentation, CategoryObject):
    @staticmethod
    def __classcall_private__(cls, function_ring, gen, name='t'):
        # Check all possible input types
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

            raise TypeError('generator must be list of coefficients or an ' \
                    'Ore polynomial')

        # Build the morphism that defines the category
        if not ore_polring_base.has_coerce_map_from(function_ring.base_ring()):
            raise TypeError('base ring of function ring must coerce to base ' \
                    'ring of Ore polynomial ring')
        gamma = function_ring.hom([ore_polring_base(gen[0])])

        # Mathematical integrity of the data is delegated to the category
        category = DrinfeldModules(gamma, name=name)
        # Check gen as Ore polynomial
        if ore_polring is not None and ore_polring is not category.codomain():
            raise ValueError(f'generator must lie in {category.codomain()}')
        # Sanity cast
        ore_polring = category.codomain()
        # Be sure to have a generator that is an Ore polynomial
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
        self._function_ring = category.domain()
        self._gen = gen
        self._morphism = self._function_ring.hom([gen])
        self._ore_polring = gen.parent()
        self._Fq = self._function_ring.base_ring()

    #################
    # Private utils #
    #################

    def _test_rank_two(self):
        r"""
        Raise ``NotImplementedError`` if the rank is not two.
        """
        if self.rank() != 2:
            raise NotImplementedError('the rank must be 2')

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
        """
        return self._base_ring

    def constant_term(self):
        r"""
        Return the constant term of the generator (`\phi_X`).

        The `A`-characteristic of the base field (see
        :meth:`sage.categories.drinfeld_modules.DrinfeldModules.characteristic`)
        is the minimal polynomial of this constant term, over the base
        ring of the function ring. Equivalently, the constant term is
        the image, by the morphism (`\gamma`) that defines the category,
        of the generator (`X`) of the polynomial ring.
    
        OUTPUT:

        - an element in the base ring
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
        """
        return self._gen

    def morphism(self):
        r"""
        Return the morphism object that defines the Drinfeld module.

        OUTPUT:

        - a ring morphism, from the function ring to the Ore polynomial
          ring
        """
        return self._morphism

    def ore_polring(self):
        r"""
        Return the Ore polynomial ring of the Drinfeld module.

        If the Drinfeld module is defined by a morphism `A \to
        K\{\tau\}`, this is the codomain `K\{\tau\}`.

        OUTPUT:

        - an Ore polynomial ring
        """
        return self._ore_polring

    def ore_variable(self):
        r"""
        Return the Ore variable.

        This is generator of the Ore polynomial ring.

        OUTPUT:

        - an Ore polynomial
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

    def change_ring(self, new_field):
        r"""
        If ``new_field`` is a field extension of the base ring, return a
        new Drinfeld module that extends ``self`` to the base ring
        ``new_field``.

        Let `f` be the morphism that defines ``self``, let `i` be the
        inclusion of ``self.ore_polring()`` into the Ore pol. ring whose
        base is ``new_field``. The morphism that defines the new
        Drinfeld module is the composition `i \circ f`.

        INPUT:

        - ``new_field`` -- the field extension of the base ring that
          serves as base ring for the new Drinfeld module

        OUTPUT:

        - a Drinfeld module
        """
        R = new_field
        if not (R.is_field() and R.is_finite() and self._Fq.is_subring(R)) \
                and self.ore_polring().base_ring().is_subring(R):
            raise ValueError('new base field must be a finite extension of ' \
                    'the base ring')
        frobenius = self.ore_polring().twisting_morphism()
        new_frobenius = R.frobenius_endomorphism(frobenius.power())
        new_ore_polring = OrePolynomialRing(R, new_frobenius,
                names=self.ore_polring().variable_names())
        return DrinfeldModule(self.function_ring(),
                new_ore_polring(self.gen()), self.characteristic())

    def height(self):
        r"""
        Return the height of the Drinfeld module.

        When the function ring is a polynomial ring, the height is 1.

        OUTPUT:

        - an integer
        """
        return Integer(1)

    def invert(self, ore_pol):
        r"""
        Find the inverse of ``ore_pol`` by the morphism that defines the
        Drinfeld module. If ``ore_pol`` is not in the image of the
        morphism, return ``None``.

        Said otherwise, return `a` if ``ore_pol`` is `phi_a`, otherwise
        return ``None``.

        INPUT:

        - ``ore_pol`` -- the Ore polynomial whose preimage we want to
          compute

        OUTPUT:

        - a polynomial
        
        ALGORITHM:
            TODO
            TODO
            TODO
            TODO
            TODO
            TODO
            TODO
            TODO
            TODO
        """
        if not ore_pol in self.ore_polring():
            raise TypeError('input must be an Ore polynomial')
        if ore_pol in self._base_ring:
            return self._Fq(ore_pol)
        r = self.rank()
        X = self.function_ring().gen()
        k = ore_pol.degree() // r
        m_lines = [[0 for _ in range(k+1)] for _ in range(k+1)]
        for i in range(k+1):
            phi_X_i = self(X**i)
            for j in range(i+1):
                m_lines[j][i] = phi_X_i[r*j]
        m = Matrix(m_lines)
        v = vector([list(ore_pol)[r*j] for j in range(k+1)])
        pre_image = self.function_ring()(list((m**(-1)) * v))
        if self(pre_image) == ore_pol:
            return pre_image
        else:
            return None

    def is_finite(self):
        r"""
        Return ``True`` if the Drinfeld module is finite, return
        ``False`` otherwise.

        OUTPUT:

        - ``True`` or ``False``
        """
        from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
        return isinstance(self, FiniteDrinfeldModule)

    def rank(self):
        r"""
        Return the rank of the Drinfeld module.

        When the function ring is a polynomial ring, the rank is the
        degree of the generator.

        OUTPUT:

        - an integer
        """
        return self.gen().degree()

    def velu(self, candidate):
        r"""
        Return a new Drinfeld module such that ``candidate`` is an
        isogeny to this module with domain ``self``. If no such isogeny
        exists, return ``None``.

        If the candidate is zero, return ``None``, as an isogeny is
        required to be non zero.

        INPUT:

        - ``candidate`` -- an Ore polynomial that defines the isogeny
          with domain ``self`` and codomain the output of the method

        OUTPUT:

        - a Drinfeld module

        ALGORITHM:
        
            We write the Ore Euclidean division `\phi_X =
            \mathrm{candidate}*q + r`, and return 
            The candidate is an isogeny if only if:

                1. The degree of the characteristic devides the height
                of the candidate. (The height of an Ore polynomial
                `P(t)` is the maximum `n` such that `t^n` right-divides
                `P(t)`.)

                2. The candidate right-divides the generator, which can
                be tested with Euclidean division.

            We test if the candidate is an isogeny, and, if it is, we
            return the quotient of the Euclidean division.

            Height and Euclidean division of Ore polynomials are
            implemented as methods of class
            :class:`sage.rings.polynomial.ore_polynomial_element.OrePolynomial`.

            Another possible algorithm is to recursively solve a system,
            see :arxiv:`2203.06970`, eq. 1.1.
        """
        if not candidate in self.ore_polring():
            raise TypeError('input must be an Ore polynomial')
        if candidate == 0:
            return None
        if not self.characteristic().degree().divides(candidate.valuation()):
            return None
        quo, rem = (candidate * self.gen()).right_quo_rem(candidate)
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
        """
        from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        return DrinfeldModuleHomset(self, other, category)

    # Rank two methods

    def delta(self):
        r"""
        If the rank is two, return `\Delta` such that the generator is
        `phi_X = \gamma(X) + g\tau + \Delta\tau^2`; if the rank is not
        two, raise an exception.

        OUTPUT:

        - an element in the base ring if the rank is two; an
          exception is raised otherwise
        """
        self._test_rank_two()
        return self.gen()[2]

    def g(self):
        r"""
        If the rank is two, return `g` such that the generator is `phi_X
        = \gamma(X) + g\tau + \Delta\tau^2`; if the rank is not two,
        raise an exception.

        OUTPUT:

        - an element in the base ring if the rank is two; an
          exception is raised otherwise
        """
        self._test_rank_two()
        return self.gen()[1]

    def j(self):
        r"""
        If the rank is two, return the j-invariant of the Drinfeld
        module; if the rank is not two, raise an exception.

        Write the generator `\phi_X = \gamma(X) + g\tau + \Delta\tau^2`.
        The j-invariant is defined by `\frac{g^{q+1}}{\Delta}`, `q`
        being the order of the base field of the polynomial ring. In our
        case, this base field is always finite, as we force the function
        ring to be of the form `\Fq[X]`.

        OUTPUT:

        - an element in the base ring if the rank is two; an
          exception is raised otherwise
        """
        self._test_rank_two()
        return (self.g()**(self._Fq.order()+1)) / self.delta()
