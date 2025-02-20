# -*- coding: utf-8 -*-
r"""
Hecke characters

A Hecke character is a character of the idele class group of a number field.

REFERENCE:

- :wikipedia:`Hecke_character`

AUTHORS:

- Robert Harron (2016-08-15): initial version

"""
# ****************************************************************************
#       Copyright (C) 2016 Robert Harron <robert.harron@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from cypari2.handle_error import PariError

from sage.arith.functions import lcm
from sage.groups.abelian_gps.dual_abelian_group import DualAbelianGroup_class
from sage.groups.abelian_gps.dual_abelian_group_element import DualAbelianGroupElement
from sage.lfunctions.dokchitser import Dokchitser
from sage.lfunctions.pari import lfun_generic, LFunction
from sage.libs.pari import pari
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import CyclotomicField
from sage.structure.unique_representation import UniqueRepresentation


class HeckeCharacter(DualAbelianGroupElement):
    def __call__(self, g):
        """
        Evaluate this Hecke character on the input ``g``.

        INPUT:

        - ``g`` -- either an element of the ray class group on which this
          character is defined, or something that can be turned into an ideal

        OUTPUT:

        The value of this character at ``g``

        EXAMPLES:

        Evaluating on an element of the ray class group. ::

            sage: F.<a> = QuadraticField(5)
            sage: mf = F.modulus(F.ideal(16), [0,1])
            sage: H = HeckeCharacterGroup(mf)
            sage: chi = H.gen(0)
            sage: chi(H.group().gen(0))
            zeta4

        Evaluating at an element of the base field. ::

            sage: chi(a / 2 + 33 / 2)
            1

        Evaluating at ideals of the base field. ::

            sage: chi(F.ideal(6))
            0
            sage: chi(F.ideal(a))
            zeta4
        """
        R = self.parent().group()
        if g.parent() is R:
            return DualAbelianGroupElement.__call__(self, g)
        K = R.number_field()
        g = K.ideal(g)
        if self.parent().modulus().finite_part().is_coprime(g):
            return DualAbelianGroupElement.__call__(self, R(g))
        return self.parent().base_ring().zero()

    def __mul__(self, right):
        """
        Multiply the two characters, if necessary by extending them to
        a common modulus.

        EXAMPLES:

        Multiplying two characters with the same parent. ::

            sage: F.<a> = QuadraticField(11)
            sage: H = HeckeCharacterGroup(F.modulus(5, [0, 1]))
            sage: prod(H.gens())
            χ0*χ1

        Multiplying two characters with different moduli. ::

            sage: chi1 = HeckeCharacterGroup(F.modulus(3, [0])).gen()^4
            sage: chi2 = HeckeCharacterGroup(F.modulus(3, [1])).gen()^4
            sage: chi = chi1 * chi2; chi
            1
            sage: chi.parent()
            Group of finite order Hecke characters modulo (Fractional ideal (3)) * ∞_0 * ∞_1
        """
        if self.parent() is right.parent():
            return DualAbelianGroupElement._mul_(self, right)
        Lmod = self.modulus()
        Rmod = right.modulus()
        if Rmod.divides(Lmod):
            return self * right.extend(self.parent())
        if Lmod.divides(Rmod):
            return self.extend(right.parent()) * right
        m = self.modulus().lcm(right.modulus())
        return self.extend(m) * right.extend(m)

    def _log_values_on_gens(self) -> tuple:
        r"""
        Return a tuple of integers `(a_j)` such that the value of this
        character on the j-th generator of the ray class group is
        `\exp(2\pi i a_j/d_j)`, where `d_j` is the order of the jth
        generator.

        This tuple is simply the exponents of the character
        with respect the generators of its Hecke character group.

        EXAMPLES::

            sage: F = QuadraticField(5)
            sage: H = HeckeCharacterGroup(F.ideal(16).modulus([0,1]))
            sage: prod(H.gens())._log_values_on_gens()
            (1, 1, 1)
        """
        return self.exponents()

    def modulus(self):
        """
        Return the modulus modulo which this character is defined.

        EXAMPLES::

            sage: F = QuadraticField(2)
            sage: H = HeckeCharacterGroup(F.modulus(F.ideal(8), [0,1]))
            sage: chi = H.gens()[0]
            sage: chi.modulus()
            (Fractional ideal (8)) * ∞_0 * ∞_1
            sage: chi.conductor()
            (Fractional ideal (2)) * ∞_0 * ∞_1
        """
        return self.parent().modulus()

    level = modulus

    def conductor(self):
        """
        Return the conductor of this character.

        This uses :pari:`bnrconductor` to compute the conductor of
        this character, i.e. the smallest modulus for which this
        character is defined.

        EXAMPLES::

            sage: F.<a> = NumberField(x^3 - 3*x -1)
            sage: H = HeckeCharacterGroup(F.modulus(3, [0,1,2]))
            sage: H.gen().conductor()
            (Fractional ideal (a^2 - 1)) * ∞_0 * ∞_1 * ∞_2
            sage: H.gen().modulus()
            (Fractional ideal (3)) * ∞_0 * ∞_1 * ∞_2
        """
        R = self.parent().group()
        K = R.number_field()
        bnr = R.pari_bnr()
        modulus = bnr.bnrconductor(self._log_values_on_gens())
        m1 = modulus[1]
        conversion = self.parent().number_field()._pari_real_places_to_sage()
        infinite = [conversion[i] for i, m1i in enumerate(m1) if m1i != 0]
        infinite.sort()
        return K.modulus(K.ideal(modulus[0]), infinite)

    def is_primitive(self) -> bool:
        """
        Determine whether this character is primitive.

        A character is primitive if its conductor is its modulus.

        EXAMPLES::

            sage: F.<a> = NumberField(x^4-5)
            sage: H = HeckeCharacterGroup(F.modulus(2, [0, 1]))
            sage: bools = [chi.is_primitive() for chi in H]; bools.sort()
            sage: bools
            [False, False, False, True]
        """
        # When a newer version of pari is part of Sage, call
        # bnrisconductor instead.
        # R = self.parent().group()
        # bnr = R.pari_bnr()
        # return bnr.bnrisconductor(self.modulus())  # not correct
        return self.conductor() == self.modulus()

    def primitive_character(self):
        """
        Return the associated primitive character.

        EXAMPLES::

            sage: F.<a> = QuadraticField(13)
            sage: H = HeckeCharacterGroup(F.ideal(-1/2*a + 7/2).modulus([0, 1]))
            sage: chi = H.gen()
            sage: chi.conductor()
            (Fractional ideal (-1/2*a + 1/2)) * ∞_0
            sage: chi0 = chi.primitive_character()
            sage: chi0.conductor()
            (Fractional ideal (-1/2*a + 1/2)) * ∞_0
            sage: chi0.parent()
            Group of finite order Hecke characters modulo (Fractional ideal (-1/2*a + 1/2)) * ∞_0
        """
        cond = self.conductor()
        if cond == self.modulus():
            return self
        H = HeckeCharacterGroup(cond)
        Rgens = [cond.equivalent_ideal_coprime_to_other(I, self.modulus())
                 for I in H.ray_class_gens()]
        return H.element_from_values_on_gens([self(I) for I in Rgens])

    def extend(self, m, check=True):
        r"""
        Return the extension of this character to the larger modulus ``m``.

        INPUT:

        - ``m`` -- a modulus that is a multiple of this Hecke
          character's modulus.

        - ``check`` -- (default: ``True``) if ``True``, ensure that this
          Hecke character's modulus divides ``m``; if it does not,
          a :exc:`ValueError` is raised.

        OUTPUT:

        The extension of this Hecke character to the one of modulus ``m``.

        EXAMPLES::

            sage: F.<a> = QuadraticField(13)
            sage: H = HeckeCharacterGroup(F.ideal(-1/2*a + 1/2).modulus([0]))
            sage: m_big = F.ideal(-1/2*a + 7/2).modulus([0, 1])
            sage: chi = H.gen().extend(m_big)
            sage: chi.parent()
            Group of finite order Hecke characters modulo (Fractional ideal (-1/2*a + 7/2)) * ∞_0 * ∞_1
            sage: chi.exponents()
            (1,)
        """
        if isinstance(m, HeckeCharacterGroup):
            if check and not self.modulus().divides(m.modulus()):
                raise ValueError("Hecke character can only be extended to "
                                 "a modulus that is a multiple of its own modulus")
            return m.element_from_values_on_gens([self(I) for I in m.ray_class_gens()])
        if check and not self.modulus().divides(m):
            raise ValueError("Hecke character can only be extended to "
                             "a modulus that is a multiple of its own modulus")
        if self.modulus() == m:
            return self
        # if not isinstance(m, HeckeCharacterGroup):
        m = HeckeCharacterGroup(m)
        return m.element_from_values_on_gens([self(I)
                                              for I in m.ray_class_gens()])

    def analytic_conductor(self):
        r"""
        Return the analytic conductor of this Hecke character.

        This is the conductor that appears in the functional equation of this
        character's `L`-function.

        The analytic conductor of a Hecke character is merely the norm
        of the finite part of its conductor times the absolute value
        of the discriminant of the number field which it is over.

        EXAMPLES:

        Over a real quadratic field::

            sage: F.<a> = QuadraticField(5)
            sage: H = HeckeCharacterGroup(F.modulus(F.ideal(a), [0,1]))
            sage: H.gens()[0].analytic_conductor()
            25

        Over an imaginary quadratic field::

            sage: F.<a> = QuadraticField(-13)
            sage: H = HeckeCharacterGroup(F.ideal(7).modulus())
            sage: H.gens()[0].analytic_conductor()
            364
        """
        N = self.conductor().finite_part().norm()
        return N * self.parent().group().number_field().disc().abs()

    def root_number(self):
        r"""
        Return the Artin root number of this Hecke character.

        EXAMPLES:

        Some root numbers over a pure cubic field::

            sage: F.<a> = NumberField(x^3-3)
            sage: H = HeckeCharacterGroup(F.ideal(5*a).modulus([0]))
            sage: chis = [chi for chi in H if chi.order() == 2]
            sage: rns = [chi.root_number() for chi in chis]; rns.sort()
            sage: rns
            [1, 1, 1]

        A quadratic character whose root number is not 1::

            sage: F.<a> = QuadraticField(-11)
            sage: H = HeckeCharacterGroup(F.modulus(a + 2))
            sage: (H.gen()^2).root_number()
            1
        """
        rn = self.parent().group().pari_bnr().bnrrootnumber(self._log_values_on_gens())
        return rn.sage()

    def _pari_init_(self):
        """
        Return this Hecke character as a Pari object.

        cf https://pari.math.u-bordeaux.fr/Events/PARI2024/talks/gchar.pdf

        and https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html

        EXAMPLES::

            sage: F.<a> = NumberField(x^3 - 3*x -1)
            sage: H = HeckeCharacterGroup(F.modulus(3, [0,1,2]))
            sage: chi = H.gen(0)
            sage: p_chi = pari(chi)

        TESTS::

            sage: K.<z> = CyclotomicField(5)
            sage: p = K.prime_above(5)
            sage: H = HeckeCharacterGroup(K.modulus(p**2))
            sage: chi = [1,1,-1,0]

        see https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html#se:bnrinit
        """
        return pari([self.parent(), self.exponents()])

    def dirichlet_series_coefficients(self, max_n):
        """
        Return some coefficients of the associated Dirichlet series.

        EXAMPLES::

            sage: F.<a> = NumberField(x^3-3)
            sage: H = HeckeCharacterGroup(F.ideal(5*a).modulus([0]))
            sage: chi = next(chi for chi in H if chi.order() == 2)
            sage: chi.dirichlet_series_coefficients(50)
            [1, 1, 0, 2, 0, 0, 0, 2, 0, 0, -1, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0,
             -1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, -1,
             0, 0, -2, 0, 1, 1, 0, 0, 0]
        """
        Idict = self.parent().number_field().ideals_of_bdd_norm(max_n)
        ans = [ZZ.one()] + [ZZ.zero()] * (max_n - 1)
        for n in range(2, max_n + 1):
            Is = Idict[n]
            if len(Is) == 0:
                continue
            ans[n - 1] = sum(self(I) for I in Is)
        return ans

    def lfunction(self, prec=53):
        r"""
        Return this character's L-function.

        This is handled by PARI.

        INPUT:

        - ``prec`` -- (default: 53) the number of bits of real precision.

        OUTPUT:

        A L-function object used to compute values of the L-function
        of this Hecke character.

        EXAMPLES:

        A totally odd character of a real quadratic field::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: mf = F.modulus(F.ideal(4), [0, 1])
            sage: H = HeckeCharacterGroup(mf)
            sage: chi = H.gens()[1]
            sage: L = chi.lfunction(prec=100); L
            Hecke L-function of χ1
            sage: [L(-n) for n in range(3)]
            [0.99999999999999999999999999999,
             0.00000000000000000000000000000,
             15.000000000000000000000000000]
        """
        from sage.lfunctions.pari import lfun_hecke
        L = LFunction(lfun_hecke(self), prec=prec)
        L.rename(f'Hecke L-function of {self}')
        return L


class HeckeCharacterGroup(DualAbelianGroup_class, UniqueRepresentation):
    r"""
    Hecke character group of a given modulus.

    EXAMPLES::

        sage: F.<a> = NumberField(x^2 - 5)
        sage: mf = F.modulus(F.prime_above(5) * F.prime_above(29), [0,1])
        sage: H = HeckeCharacterGroup(mf); H
        Group of finite order Hecke characters modulo (Fractional ideal (-11/2*a - 5/2)) * ∞_0 * ∞_1
        sage: [[chi(F.ideal(31)), chi(F.ideal(-12672))] for chi in H.gens()]
        [[zeta4, 1], [1, -1]]

    TESTS::

        sage: F.<a> = QuadraticField(5)
        sage: H = HeckeCharacterGroup(F.modulus(F.ideal(16), [0,1]))
        sage: H.gens()
        (χ0, χ1, χ2)
    """
    Element = HeckeCharacter

    @staticmethod
    def __classcall_private__(cls, modulus, base_ring=None, names=None):
        """
        Prepare input for unique representation.

        EXAMPLES::

            sage: F.<a> = NumberField(x^2 - 7)
            sage: mf = F.modulus(F.prime_above(7), [0,1])
            sage: H = HeckeCharacterGroup(mf); H  # indirect doctest
            Group of finite order Hecke characters modulo (Fractional ideal (a)) * ∞_0 * ∞_1
        """
        ray_class_group = modulus.number_field().ray_class_group(modulus)
        if base_ring is None:
            base_ring = CyclotomicField(lcm(ray_class_group.gens_orders()))
        if names is None:
            names = 'χ'
        return super().__classcall__(cls, ray_class_group, base_ring, names)

    def __init__(self, ray_class_group, base_ring, names) -> None:
        """
        Create the Hecke character group.

        EXAMPLES::

            sage: F.<a> = NumberField(x^2 - 61)
            sage: mf = F.modulus(F.prime_above(5), [0,1])
            sage: H = HeckeCharacterGroup(mf); H  # indirect doctest
            Group of finite order Hecke characters modulo (Fractional ideal (1/2*a + 9/2)) * ∞_0 * ∞_1
        """
        DualAbelianGroup_class.__init__(self, ray_class_group, names, base_ring)

    def _repr_(self) -> str:
        """
        Return the string representation.

        EXAMPLES::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: mf = F.modulus(F.prime_above(5) * F.prime_above(29), [0,1])
            sage: H = HeckeCharacterGroup(mf); H
            Group of finite order Hecke characters modulo (Fractional ideal (-11/2*a - 5/2)) * ∞_0 * ∞_1
        """
        return f'Group of finite order Hecke characters modulo {self.modulus()}'

    def modulus(self):
        """
        Return the modulus.

        EXAMPLES::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: mf = F.modulus(F.prime_above(5) * F.prime_above(29), [0,1])
            sage: H = HeckeCharacterGroup(mf); H.modulus()
            (Fractional ideal (-11/2*a - 5/2)) * ∞_0 * ∞_1
        """
        return self.group().modulus()

    level = modulus

    def number_field(self):
        """
        Return the number field.

        EXAMPLES::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: mf = F.modulus(F.prime_above(5) * F.prime_above(29), [0,1])
            sage: H = HeckeCharacterGroup(mf); H.number_field()
            Number Field in a with defining polynomial x^2 - 5
        """
        return self.group().number_field()

    def _pari_init_(self):
        """
        Return this Hecke character group as a Pari object.

        cf https://pari.math.u-bordeaux.fr/Events/PARI2024/talks/gchar.pdf

        and https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html

        EXAMPLES::

            sage: F.<a> = NumberField(x^3 - 3*x -1)
            sage: H = HeckeCharacterGroup(F.modulus(3, [0,1,2]))
            sage: pH = pari(H)
        """
        field = self.number_field().pari_nf().bnfinit()
        return pari.bnrinit(field, self.modulus())

    def ray_class_gens(self) -> tuple:
        """
        Return the ray class generators.

        EXAMPLES::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: mf = F.modulus(F.prime_above(5) * F.prime_above(29), [0,1])
            sage: H = HeckeCharacterGroup(mf); H.ray_class_gens()
            (Fractional ideal (31), Fractional ideal (88))
        """
        return self.group().gens_ideals()

    def element_from_values_on_gens(self, vals):
        """
        missing documentation
        """
        gens_orders = self.gens_orders()
        if len(vals) != len(gens_orders):
            raise ValueError("Incorrect number of values specified. %s specified, but needed %s" % (len(vals), len(gens_orders)))
        exponents = [gens_orders[i].divide_knowing_divisible_by(vals[i].multiplicative_order())
                     if vals[i] != 1 else ZZ.zero()
                     for i in range(len(gens_orders))]
        return self.element_class(self, exponents)

    # def _element_constructor_(self, *args, **kwds):
    #    if isinstance(args[0], (str,  bytes)):
    #        raise TypeError("wrong type to coerce into HeckeCharacterGroup")
    #    try:
    #        n = len(args[0])
    #    except TypeError:
    #        return DualAbelianGroup_class._element_constructor_(self, *args, **kwds)
    #    gens_orders = self.gens_orders()
    #    if n != len(gens_orders):
    #        return DualAbelianGroup_class._element_constructor_(self, *args, **kwds)
    #    exponents = [gens_orders[i].divide_knowing_divisible_by(args[0][i].multiplicative_order()) for i in range(n)]
    #    return self.element_class(self, exponents)
