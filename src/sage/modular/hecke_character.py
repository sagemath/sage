# -*- coding: utf-8 -*-
r"""
Hecke characters

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
from sage.structure.factory import UniqueFactory
from sage.arith.all import lcm
from sage.rings.integer_ring import ZZ
from sage.groups.abelian_gps.dual_abelian_group_element import DualAbelianGroupElement
from sage.groups.abelian_gps.dual_abelian_group import DualAbelianGroup_class
from sage.lfunctions.dokchitser import Dokchitser

from cypari2.handle_error import PariError


class HeckeCharacter(DualAbelianGroupElement):
    def __call__(self, g):
        """
        Evaluate this Hecke character on the input ``g``.

        INPUT:

        - ``g`` -- either an element of the ray class group on which this character
        is defined, or something that can be turned into an ideal.

        OUTPUT:

        The value of this character at ``g``.

        EXAMPLES:

        Evaluating on an element of the ray class group. ::

            sage: F.<a> = QuadraticField(5)
            sage: H = HeckeCharacterGroup(F.modulus(F.ideal(16), [0,1]))
            sage: chi = H.gens()[0]; chi(H.group().gens()[0])
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
            chi0*chi1

        Multiplying two characters with different moduli. ::

            sage: chi1 = HeckeCharacterGroup(F.modulus(3, [0])).gen()^4
            sage: chi2 = HeckeCharacterGroup(F.modulus(3, [1])).gen()^4
            sage: chi = chi1 * chi2; chi
            1
            sage: chi.parent()
            Group of finite order Hecke characters modulo (Fractional ideal (3)) * infinity_0 * infinity_1
        """
        if self.parent() is right.parent():
            return DualAbelianGroupElement._mul_(self, right)
        Lmod = self.modulus()
        Rmod = right.modulus()
        if Rmod.divides(Lmod):
            return self * right.extend(self.parent())
        elif Lmod.divides(Rmod):
            return self.extend(right.parent()) * right
        m = self.modulus().lcm(right.modulus())
        return self.extend(m) * right.extend(m)

    def _log_values_on_gens(self):
        r"""
        Return a tuple of integers `(a_j)` such that the value of this
        character on the jth generator of the ray class group is
        `\exp(2\pi ia_j/d_j)`, where `d_j` is the order of the jth
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
            (Fractional ideal (8)) * infinity_0 * infinity_1
            sage: chi.conductor()
            (Fractional ideal (2)) * infinity_0 * infinity_1
        """
        return self.parent().modulus()

    def level(self):
        """
        Return the modulus modulo which this character is defined.

        An alias for :meth:`modulus`.

        EXAMPLES::

            sage: F = QuadraticField(5)
            sage: H = HeckeCharacterGroup(F.ideal(F.gen()).modulus([0,1]))
            sage: H.gens()[0].level()
            (Fractional ideal (a)) * infinity_0 * infinity_1
        """
        return self.modulus()

    def conductor(self):
        """
        Return the conductor of this character.

        Uses pari to compute the conductor of this character, i.e. the smallest
        modulus for which this character is defined.

        EXAMPLES::

            sage: F.<a> = NumberField(x^3 - 3*x -1)
            sage: H = HeckeCharacterGroup(F.modulus(3, [0,1,2]))
            sage: H.gen().conductor()
            (Fractional ideal (-a^2 + 1)) * infinity_0 * infinity_1 * infinity_2
            sage: H.gen().modulus()
            (Fractional ideal (3)) * infinity_0 * infinity_1 * infinity_2
        """
        R = self.parent().group()
        K = R.number_field()
        bnr = R.pari_bnr()
        modulus = bnr.bnrconductor(self._log_values_on_gens())
        # in a newer version of pari, this is replaced with bnrconductor
        infinite = []
        m1 = modulus[1]
        conversion = self.parent().number_field()._pari_real_places_to_sage()
        for i in range(len(m1)):
            if m1[i] != 0:
                infinite.append(conversion[i])
        infinite.sort()
        return K.modulus(K.ideal(modulus[0]), infinite)

    def is_primitive(self):
        """
        Determine whether this character is primitive.

        This means it checks whether its conductor is its modulus.

        EXAMPLES::

            sage: F.<a> = NumberField(x^4-5)
            sage: H = HeckeCharacterGroup(F.modulus(2, [0, 1]))
            sage: bools = [chi.is_primitive() for chi in H]; bools.sort()
            sage: bools
            [False, False, False, True]
        """
        # When a newer version of pari is part of Sage, call
        # bnrisconductor instead.
        return self.conductor() == self.modulus()

    def primitive_character(self):
        """
        Return the associated primitive character.

        EXAMPLES::

            sage: F.<a> = QuadraticField(13)
            sage: H = HeckeCharacterGroup(F.ideal(-1/2*a + 7/2).modulus([0, 1]))
            sage: chi = H.gen()
            sage: chi.conductor()
            (Fractional ideal (-1/2*a + 1/2)) * infinity_0
            sage: chi0 = chi.primitive_character()
            sage: chi0.conductor()
            (Fractional ideal (-1/2*a + 1/2)) * infinity_0
            sage: chi0.parent()
            Group of finite order Hecke characters modulo (Fractional ideal (-1/2*a + 1/2)) * infinity_0
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

        - ``m`` -- a modulus that is a multiple of this Hecke character's modulus.
        - ``check`` -- (default: ``True``) if ``True``, ensure that this Hecke character's
        modulus divides ``m``; if it does not, a ``ValueError`` is raised.

        OUTPUT:

        The extension of this Hecke character to the one of modulus ``m``.

        EXAMPLES::

            sage: F.<a> = QuadraticField(13)
            sage: H = HeckeCharacterGroup(F.ideal(-1/2*a + 1/2).modulus([0]))
            sage: m_big = F.ideal(-1/2*a + 7/2).modulus([0, 1])
            sage: chi = H.gen().extend(m_big)
            sage: chi.parent()
            Group of finite order Hecke characters modulo (Fractional ideal (-1/2*a + 7/2)) * infinity_0 * infinity_1
            sage: chi.exponents()
            (1,)
        """
        if isinstance(m, HeckeCharacterGroup_class):
            if check and not self.modulus().divides(m.modulus()):
                raise ValueError("Hecke character can only be extended to a modulus that is a multiple of its own modulus.")
            return m.element_from_values_on_gens([self(I) for I in m.ray_class_gens()])
        if check and not self.modulus().divides(m):
            raise ValueError("Hecke character can only be extended to a modulus that is a multiple of its own modulus.")
        if self.modulus() == m:
            return self
        #if not isinstance(m, HeckeCharacterGroup_class):
        m = HeckeCharacterGroup(m)
        return m.element_from_values_on_gens([self(I) for I in m.ray_class_gens()])

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
        #try-except only necessary because of http://pari.math.u-bordeaux.fr/cgi-bin/bugreport.cgi?bug=1848
        #According to Belabas this bug seems to occur because of a problem dealing with
        #imprimitive characters, thus the code in the except clause below. Once we incorporate a version
        #a pari where this bug is fixed, the command in the try itself should be sufficient.
        #The first example above (with Q(cubert(3))) fails without the except clause and so can
        #be used to test if the pari bug has been fixed.
        try:
            rn = self.parent().group().pari_bnr().bnrrootnumber(self._log_values_on_gens())
        except PariError:
            chi = self.primitive_character()
            rn = chi.parent().group().pari_bnr().bnrrootnumber(chi._log_values_on_gens())
        return rn.sage()

    def dirichlet_series_coefficients(self, max_n):
        """
        documentation missing here
        """
        Idict = self.parent().number_field().ideals_of_bdd_norm(max_n)
        ans = [ZZ(1)] + [ZZ(0)] * (max_n - 1)
        for n in range(2, max_n + 1):
            Is = Idict[n]
            if len(Is) == 0:
                continue
            ans[n - 1] = sum(self(I) for I in Is)
        return ans

    def Lfunction(self, prec=53):
        r"""
        Return this character's L-function as a :class:`sage.lfunctions.dokchitser.Dokchitser` object.

        INPUT:

        - ``prec`` -- (default: 53) the number of bits of real precision.

        OUTPUT:

        A Dokchitser L-function object used to compute values of the L-function of this Hecke character.

        EXAMPLES:

        A totally odd character of a real quadratic field::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: mf = F.modulus(F.ideal(4), [0, 1])
            sage: H = HeckeCharacterGroup(mf)
            sage: chi = H.gens()[1]
            sage: L = chi.Lfunction(); L
            Hecke L-function of chi1
            sage: [L(-n) for n in range(3)]
            [1.00000000000000, 0.000000000000000, 15.0000000000000]
        """
        # Figure out Gamma factors for more general characters
        gamma_factors = [0] * self.parent().number_field().degree()
        for i in self.conductor().infinite_part():
            gamma_factors[i] = 1
        ana_cond = self.analytic_conductor()
        rn = self.root_number()
        it_worked = False
        number_of_allocs = 0
        while not it_worked:
            #Sage automatically ups the memory allocation, but this messes with L.num_coeffs for some reason I have yet to figure out, so L needs to be remade completely every time memory is auto-allocated. If we can control the auto-allocation process/rewrite the Dokchitser code, then we can remove the extra stuff here.
            L = Dokchitser(ana_cond, gamma_factors, 1, rn, prec=prec)
            for n in range(number_of_allocs):
                L.gp().eval('allocatemem()')
            try:
                L.init_coeffs(self.dirichlet_series_coefficients(L.num_coeffs()))
                it_worked = True
            except RuntimeError:
                number_of_allocs += 1
        L.rename('Hecke L-function of %s' % self)
        return L


class HeckeCharacterGroup_class(DualAbelianGroup_class):
    r"""
    EXAMPLES::

        sage: F.<a> = NumberField(x^2 - 5)
        sage: mf = F.modulus(F.prime_above(5) * F.prime_above(29), [0,1])
        sage: H = HeckeCharacterGroup(mf); H
        Group of finite order Hecke characters modulo (Fractional ideal (-11/2*a - 5/2)) * infinity_0 * infinity_1
        sage: [[chi(F.ideal(31)), chi(F.ideal(-12672))] for chi in H.gens()]
        [[zeta4, 1], [1, -1]]
    """
    Element = HeckeCharacter

    def __init__(self, ray_class_group, base_ring=None, names=None):
        if base_ring is None:
            from sage.rings.number_field.number_field import CyclotomicField
            base_ring = CyclotomicField(lcm(ray_class_group.gens_orders()))
        if names is None:
            names = 'chi'
        DualAbelianGroup_class.__init__(self, ray_class_group, names, base_ring)

    def _repr_(self):
        return 'Group of finite order Hecke characters modulo ' + str(self.modulus())

    def modulus(self):
        return self.group().modulus()

    def level(self):
        return self.modulus()

    def number_field(self):
        return self.group().number_field()

    def ray_class_gens(self):
        return self.group().gens_ideals()

    def element_from_values_on_gens(self, vals):
        gens_orders = self.gens_orders()
        if len(vals) != len(gens_orders):
            raise ValueError("Incorrect number of values specified. %s specified, but needed %s" % (len(vals), len(gens_orders)))
        exponents = [gens_orders[i].divide_knowing_divisible_by(vals[i].multiplicative_order()) if vals[i] != 1 else ZZ.zero()
                     for i in range(len(gens_orders))]
        return self.element_class(self, exponents)

    #def _element_constructor_(self, *args, **kwds):
    #    if isinstance(args[0], (str,  bytes)):
    #        raise TypeError("Wrong type to coerce into HeckeCharacterGroup.")
    #    try:
    #        n = len(args[0])
    #    except TypeError:
    #        return DualAbelianGroup_class._element_constructor_(self, *args, **kwds)
    #    gens_orders = self.gens_orders()
    #    if n != len(gens_orders):
    #        return DualAbelianGroup_class._element_constructor_(self, *args, **kwds)
    #    exponents = [gens_orders[i].divide_knowing_divisible_by(args[0][i].multiplicative_order()) for i in range(n)]
    #    return self.element_class(self, exponents)

class HeckeCharacterGroupFactory(UniqueFactory):
    def create_key(self, modulus, base_ring=None, names=None):#(m, base_ring=None, names=None):
        #from sage.structure.category_object import normalize_names
        #names = normalize_names()
        if names is None:
            names = 'chi'
        return (base_ring, modulus, names)

    def create_object(self, version, key, **extra_args):
        base_ring, modulus, names = key
        return HeckeCharacterGroup_class(modulus.number_field().ray_class_group(modulus), base_ring, names)


HeckeCharacterGroup = HeckeCharacterGroupFactory("HeckeCharacterGroup")
