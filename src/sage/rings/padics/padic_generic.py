"""
p-Adic Generic

A generic superclass for all p-adic parents.

AUTHORS:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
- Julian Rueth (2013-03-16): test methods for basic arithmetic

"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#                               Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from sage.misc.prandom import sample
from sage.misc.misc import some_tuples
from copy import copy

from sage.structure.richcmp import richcmp
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.morphism import Morphism
from sage.categories.fields import Fields
from sage.rings.infinity import infinity
from .local_generic import LocalGeneric
from sage.rings.ring import PrincipalIdealDomain
from sage.rings.integer import Integer
from sage.rings.padics.padic_printing import pAdicPrinter
from sage.rings.padics.precision_error import PrecisionError
from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import richcmp_not_equal


class pAdicGeneric(PrincipalIdealDomain, LocalGeneric):
    def __init__(self, base, p, prec, print_mode, names, element_class, category=None):
        """
        Initialization.

        INPUT:

            - base -- Base ring.
            - p -- prime
            - print_mode -- dictionary of print options
            - names -- how to print the uniformizer
            - element_class -- the class for elements of this ring

        EXAMPLES::

            sage: R = Zp(17) #indirect doctest
        """
        if category is None:
            if self.is_field():
                category = Fields()
            else:
                category = PrincipalIdealDomains()
        category = category.Metric().Complete()
        LocalGeneric.__init__(self, base, prec, names, element_class, category)
        self._printer = pAdicPrinter(self, print_mode)

    def some_elements(self):
        r"""
        Returns a list of elements in this ring.

        This is typically used for running generic tests (see :class:`TestSuite`).

        EXAMPLES::

            sage: Zp(2,4).some_elements()
            [0, 1 + O(2^4), 2 + O(2^5), 1 + 2^2 + 2^3 + O(2^4), 2 + 2^2 + 2^3 + 2^4 + O(2^5)]
        """
        p = self(self.prime())
        a = self.gen()
        one = self.one()
        L = [self.zero(), one, p, (one+p+p).inverse_of_unit(), p-p**2]
        if a != p:
            L.extend([a, (one + a + p).inverse_of_unit()])
        if self.is_field():
            L.extend([~(p-p-a),p**(-20)])
        return L

    def _modified_print_mode(self, print_mode):
        """
        Returns a dictionary of print options, starting with self's
        print options but modified by the options in the dictionary
        print_mode.

        INPUT:

            - print_mode -- dictionary with keys in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']

        EXAMPLES::

            sage: R = Zp(5)
            sage: R._modified_print_mode({'mode': 'bars'})['ram_name']
            '5'
        """
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet', 'show_prec']:
            if option not in print_mode:
                print_mode[option] = self._printer.dict()[option]
        return print_mode

    def ngens(self):
        """
        Returns the number of generators of self.

        We conventionally define this as 1: for base rings, we take a
        uniformizer as the generator; for extension rings, we take a
        root of the minimal polynomial defining the extension.

        EXAMPLES::

            sage: Zp(5).ngens()
            1
            sage: Zq(25,names='a').ngens()
            1
        """
        return 1

    def gens(self):
        """
        Returns a list of generators.

        EXAMPLES::

            sage: R = Zp(5); R.gens()
            [5 + O(5^21)]
            sage: Zq(25,names='a').gens()
            [a + O(5^20)]
            sage: S.<x> = ZZ[]; f = x^5 + 25*x -5; W.<w> = R.ext(f); W.gens()
            [w + O(w^101)]
        """
        return [self.gen()]

    def __richcmp__(self, other, op):
        """
        Return 0 if self == other, and 1 or -1 otherwise.

        We consider two p-adic rings or fields to be equal if they are
        equal mathematically, and also have the same precision cap and
        printing parameters.

        EXAMPLES::

            sage: R = Qp(7)
            sage: S = Qp(7,print_mode='val-unit')
            sage: R == S
            False
            sage: S = Qp(7,type='capped-rel')
            sage: R == S
            True
            sage: R is S
            True
        """
        if not isinstance(other, pAdicGeneric):
            return NotImplemented

        lx = self.prime()
        rx = other.prime()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.precision_cap()
        rx = other.precision_cap()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        return self._printer.richcmp_modes(other._printer, op)

    #def ngens(self):
    #    return 1

    #def gen(self, n = 0):
    #    if n != 0:
    #        raise IndexError, "only one generator"
    #    return self(self.prime())

    def print_mode(self):
        r"""
        Returns the current print mode as a string.

        INPUT:

            self -- a p-adic field

        OUTPUT:

            string -- self's print mode

        EXAMPLES::

            sage: R = Qp(7,5, 'capped-rel')
            sage: R.print_mode()
            'series'
        """
        return self._printer._print_mode()

    #def element_class(self):
    #    return self._element_class

    def characteristic(self):
        r"""
        Returns the characteristic of self, which is always 0.

        INPUT:

            self -- a p-adic parent

        OUTPUT:

            integer -- self's characteristic, i.e., 0

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.characteristic()
            0
        """
        return Integer(0)

    def prime(self):
        """
        Returns the prime, ie the characteristic of the residue field.

        INPUT:

            self -- a p-adic parent

        OUTPUT:

            integer -- the characteristic of the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: R.prime()
            3
        """
        return self.prime_pow._prime()

    def uniformizer_pow(self, n):
        """
        Returns p^n, as an element of self.

        If n is infinity, returns 0.

        EXAMPLES::

            sage: R = Zp(3, 5, 'fixed-mod')
            sage: R.uniformizer_pow(3)
            3^3 + O(3^5)
            sage: R.uniformizer_pow(infinity)
            O(3^5)
        """
        if n is infinity:
            return self(0)
        return self(self.prime_pow.pow_Integer_Integer(n))

    def _unram_print(self):
        """
        For printing.  Will be None if the unramified subextension of self is of degree 1 over Z_p or Q_p.

        EXAMPLES::

            sage: Zp(5)._unram_print()
        """
        return None

    def residue_characteristic(self):
        """
        Return the prime, i.e., the characteristic of the residue field.

        OUTPUT:

        integer -- the characteristic of the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: R.residue_characteristic()
            3
        """
        return self.prime()

    def residue_class_field(self):
        """
        Returns the residue class field.

        INPUT:

            self -- a p-adic ring

        OUTPUT:

            the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: k = R.residue_class_field()
            sage: k
            Finite Field of size 3
        """
        from sage.rings.finite_rings.finite_field_constructor import GF
        return GF(self.prime())

    def residue_field(self):
        """
        Returns the residue class field.

        INPUT:

            self -- a p-adic ring

        OUTPUT:

            the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: k = R.residue_field()
            sage: k
            Finite Field of size 3
        """
        return self.residue_class_field()

    def residue_ring(self, n):
        """
        Returns the quotient of the ring of integers by the nth power of the maximal ideal.

        EXAMPLES::

            sage: R = Zp(11)
            sage: R.residue_ring(3)
            Ring of integers modulo 1331
        """
        from sage.rings.finite_rings.integer_mod_ring import Zmod
        return Zmod(self.prime()**n)

    def residue_system(self):
        """
        Returns a list of elements representing all the residue classes.

        INPUT:

            self -- a p-adic ring

        OUTPUT:

            list of elements -- a list of elements representing all the residue classes

        EXAMPLES::

            sage: R = Zp(3, 5,'fixed-mod')
            sage: R.residue_system()
            [O(3^5), 1 + O(3^5), 2 + O(3^5)]
        """
        return [self(i) for i in self.residue_class_field()]

    def fraction_field(self, print_mode=None):
        r"""
        Returns the fraction field of this ring or field.

        For `\ZZ_p`, this is the `p`-adic field with the same options,
        and for extensions, it is just the extension of the fraction
        field of the base determined by the same polynomial.

        The fraction field of a capped absolute ring is capped relative,
        and that of a fixed modulus ring is floating point.

        INPUT:

        - ``print_mode`` -- a dictionary containing print options.
          Defaults to the same options as this ring.

        OUTPUT:

        - the fraction field of this ring.

        EXAMPLES::

            sage: R = Zp(5, print_mode='digits')
            sage: K = R.fraction_field(); repr(K(1/3))[3:]
            '31313131313131313132'
            sage: L = R.fraction_field({'max_ram_terms':4}); repr(L(1/3))[3:]
            '3132'
            sage: U.<a> = Zq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)
            sage: U.fraction_field()
            Unramified Extension in a defined by x^4 + 7*x^2 + 10*x + 3 with capped relative precision 6 over 17-adic Field
            sage: U.fraction_field({"pos":False}) == U.fraction_field()
            False
        """
        if self.is_field() and print_mode is None:
            return self
        if print_mode is None:
            return self.change(field=True)
        else:
            return self.change(field=True, **print_mode)

    def integer_ring(self, print_mode=None):
        r"""
        Returns the ring of integers of this ring or field.

        For `\QQ_p`, this is the `p`-adic ring with the same options,
        and for extensions, it is just the extension of the ring
        of integers of the base determined by the same polynomial.

        INPUT:

        - ``print_mode`` -- a dictionary containing print options.
          Defaults to the same options as this ring.

        OUTPUT:

        - the ring of elements of this field with nonnegative valuation.

        EXAMPLES::

            sage: K = Qp(5, print_mode='digits')
            sage: R = K.integer_ring(); repr(R(1/3))[3:]
            '31313131313131313132'
            sage: S = K.integer_ring({'max_ram_terms':4}); repr(S(1/3))[3:]
            '3132'
            sage: U.<a> = Qq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)
            sage: U.integer_ring()
            Unramified Extension in a defined by x^4 + 7*x^2 + 10*x + 3 with capped relative precision 6 over 17-adic Ring
            sage: U.fraction_field({"pos":False}) == U.fraction_field()
            False
        """
        #Currently does not support fields with non integral defining polynomials.  This should change when the padic_general_extension framework gets worked out.
        if not self.is_field() and print_mode is None:
            return self
        if print_mode is None:
            return self.change(field=False)
        else:
            return self.change(field=False, **print_mode)

    def teichmuller(self, x, prec = None):
        r"""
        Returns the teichmuller representative of x.

        INPUT:

            - self -- a p-adic ring
            - x -- something that can be cast into self

        OUTPUT:

            - element -- the teichmuller lift of x

        EXAMPLES::

            sage: R = Zp(5, 10, 'capped-rel', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Qp(5, 10,'capped-rel','series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5, 10, 'capped-abs', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5, 10, 'fixed-mod', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W.teichmuller(3); y
            3 + 3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + 3*w^15 + 2*w^16 + 3*w^17 + w^18 + 3*w^19 + 3*w^20 + 2*w^21 + 2*w^22 + 3*w^23 + 4*w^24 + O(w^25)
            sage: y^5 == y
            True
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = A.teichmuller(1 + 2*a - a^2); b
            (4*a^2 + 2*a + 1) + 2*a*5 + (3*a^2 + 1)*5^2 + (a + 4)*5^3 + (a^2 + a + 1)*5^4 + O(5^5)
            sage: b^125 == b
            True

        We check that :trac:`23736` is resolved::

            sage: R.teichmuller(GF(5)(2))
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + O(5^5)

        AUTHORS:

        - Initial version: David Roe
        - Quadratic time version: Kiran Kedlaya <kedlaya@math.mit.edu> (3/27/07)
        """
        ans = self(x) if prec is None else self(x, prec)
        # Since teichmuller representatives are defined at infinite precision,
        # we can lift to precision prec, as long as the absolute precision of ans is positive.
        if ans.precision_absolute() <= 0:
            raise ValueError("Not enough precision to determine Teichmuller representative")
        if ans.valuation() > 0:
            return self(0) if prec is None else self(0, prec)
        ans = ans.lift_to_precision(prec)
        if ans is x:
            ans = copy(ans)
        ans._teichmuller_set_unsafe()
        return ans

    def teichmuller_system(self):
        r"""
        Returns a set of teichmuller representatives for the invertible elements of `\ZZ / p\ZZ`.

        INPUT:

        - self -- a p-adic ring

        OUTPUT:

        - list of elements -- a list of teichmuller representatives for the invertible elements of `\ZZ / p\ZZ`

        EXAMPLES::

            sage: R = Zp(3, 5,'fixed-mod', 'terse')
            sage: R.teichmuller_system()
            [1 + O(3^5), 242 + O(3^5)]

        Check that :trac:`20457` is fixed::

            sage: F.<a> = Qq(5^2,6)
            sage: F.teichmuller_system()[3]
            (2*a + 2) + (4*a + 1)*5 + 4*5^2 + (2*a + 1)*5^3 + (4*a + 1)*5^4 + (2*a + 3)*5^5 + O(5^6)

        NOTES:

        Should this return 0 as well?
        """
        R = self.residue_class_field()
        prec = self.precision_cap()
        return [self.teichmuller(self(i).lift_to_precision(prec)) for i in R if i != 0]

#     def different(self):
#         raise NotImplementedError

#     def automorphisms(self):
#         r"""
#         Returns the group of automorphisms of `\ZZ_p`, i.e. the trivial group.
#         """
#         raise NotImplementedError

#     def galois_group(self):
#         r"""
#         Returns the Galois group of `\ZZ_p`, i.e. the trivial group.
#         """
#         raise NotImplementedError

#     def hasGNB(self):
#         r"""
#         Returns whether or not `\ZZ_p` has a Gauss Normal Basis.
#         """
#         raise NotImplementedError

    def extension(self, modulus, prec = None, names = None, print_mode = None, implementation='FLINT', **kwds):
        """
        Create an extension of this p-adic ring.

        EXAMPLES::

            sage: k = Qp(5)
            sage: R.<x> = k[]
            sage: l.<w> = k.extension(x^2-5); l
            Eisenstein Extension in w defined by x^2 - 5 with capped relative precision 40 over 5-adic Field

            sage: F = list(Qp(19)['x'](cyclotomic_polynomial(5)).factor())[0][0]
            sage: L = Qp(19).extension(F, names='a')
            sage: L
            Unramified Extension in a defined by x^2 + 8751674996211859573806383*x + 1 with capped relative precision 20 over 19-adic Field
        """
        from sage.rings.padics.factory import ExtensionFactory
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'print_mode': print_mode}
        else:
            if not isinstance(print_mode, dict):
                print_mode = dict(print_mode)
            for option in ['mode', 'pos', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
                if option in print_mode:
                    print_mode["print_" + option] = print_mode[option]
                    del print_mode[option]
                elif "print_" + option not in print_mode:
                    if "print_" + option in kwds:
                        print_mode["print_" + option] = kwds["print_" + option]
                    else:
                        print_mode["print_" + option] = self._printer.dict()[option]
            for option in ['ram_name', 'unram_name', 'var_name']:
                if option not in print_mode:
                    if option in kwds:
                        print_mode[option] = kwds[option]
                    else:
                        print_mode[option] = self._printer.dict()[option]
        return ExtensionFactory(base=self, modulus=modulus, prec=prec, names=names, check = True, implementation=implementation, **print_mode)

    def _test_add(self, **options):
        """
        Test addition of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_add()

        .. SEEALSO::

            :class:`TestSuite`

        """
        tester = self._tester(**options)
        elements = tester.some_elements()

        for x in elements:
            y = x + self.zero()
            tester.assertEqual(y,x)
            tester.assertEqual(y.precision_absolute(),x.precision_absolute())
            tester.assertEqual(y.precision_relative(),x.precision_relative())

        for x,y in some_tuples(elements, 2, tester._max_runs):
            z = x + y
            tester.assertIs(z.parent(), self)
            zprec = min(x.precision_absolute(), y.precision_absolute())
            if not self.is_floating_point():
                tester.assertEqual(z.precision_absolute(), zprec)
            tester.assertGreaterEqual(z.valuation(), min(x.valuation(),y.valuation()))
            if x.valuation() != y.valuation():
                tester.assertEqual(z.valuation(), min(x.valuation(),y.valuation()))
            tester.assert_(y.is_equal_to(z-x,zprec))
            tester.assert_(x.is_equal_to(z-y,zprec))

    def _test_sub(self, **options):
        """
        Test subtraction on elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_sub()

        .. SEEALSO::

            :class:`TestSuite`

        """
        tester = self._tester(**options)

        elements = list(tester.some_elements())
        for x in elements:
            y = x - self.zero()
            tester.assertEqual(y, x)
            tester.assertEqual(y.precision_absolute(), x.precision_absolute())
            tester.assertEqual(y.precision_relative(), x.precision_relative())

        for x,y in some_tuples(elements, 2, tester._max_runs):
            z = x - y
            tester.assertIs(z.parent(), self)
            zprec = min(x.precision_absolute(), y.precision_absolute())
            if not self.is_floating_point():
                tester.assertEqual(z.precision_absolute(), zprec)
            tester.assertGreaterEqual(z.valuation(), min(x.valuation(),y.valuation()))
            if x.valuation() != y.valuation():
                tester.assertEqual(z.valuation(), min(x.valuation(),y.valuation()))
            tester.assert_((-y).is_equal_to(z - x,zprec))
            tester.assert_(x.is_equal_to(z + y,zprec))

    def _test_invert(self, **options):
        """
        Test multiplicative inversion of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_invert()

        .. SEEALSO::

            :class:`TestSuite`

        """
        tester = self._tester(**options)

        elements = tester.some_elements()
        for x in elements:
            try:
                y = ~x
            except (ZeroDivisionError, PrecisionError, ValueError):
                tester.assertFalse(x.is_unit())
                if not self.is_fixed_mod(): tester.assertTrue(x.is_zero())
            else:
                try:
                    e = y * x
                except ZeroDivisionError:
                    tester.assertTrue(self.is_floating_point() and (x.is_zero() or y.is_zero()))
                else:
                    tester.assertFalse(x.is_zero())
                    tester.assertIs(y.parent(), self if self.is_fixed_mod() else self.fraction_field())
                    tester.assertTrue(e.is_one())
                    tester.assertEqual(e.precision_relative(), x.precision_relative())
                    tester.assertEqual(y.valuation(), -x.valuation())

    def _test_mul(self, **options):
        """
        Test multiplication of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_mul()

        .. SEEALSO::

            :class:`TestSuite`

        """
        tester = self._tester(**options)

        elements = list(tester.some_elements())
        for x,y in some_tuples(elements, 2, tester._max_runs):
            z = x * y
            tester.assertIs(z.parent(), self)
            if self.is_capped_relative() or self.is_floating_point():
                tester.assertEqual(z.precision_relative(), min(x.precision_relative(), y.precision_relative()))
            else:
                tester.assertLessEqual(z.precision_relative(), min(x.precision_relative(), y.precision_relative()))
            if not z.is_zero():
                tester.assertEqual(z.valuation(), x.valuation() + y.valuation())

    def _test_div(self, **options):
        """
        Test division of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_div()

        .. SEEALSO::

            :class:`TestSuite`

        """
        tester = self._tester(**options)

        elements = list(tester.some_elements())
        for x,y in some_tuples(elements, 2, tester._max_runs):
            try:
                z = x / y
            except (ZeroDivisionError, PrecisionError, ValueError):
                if self.is_fixed_mod(): tester.assertFalse(y.is_unit())
                else: tester.assertTrue(y.is_zero())
            else:
                try:
                    xx = z*y
                except ZeroDivisionError:
                    tester.assertTrue(self.is_floating_point() and (z.is_zero() or y.is_zero()))
                else:
                    tester.assertFalse(y.is_zero())
                    tester.assertIs(z.parent(), self if self.is_fixed_mod() else self.fraction_field())
                    tester.assertEqual(z.precision_relative(), min(x.precision_relative(), y.precision_relative()))
                    tester.assertEqual(z.valuation(), x.valuation() - y.valuation())
                    tester.assertEqual(xx, x)

    def _test_neg(self, **options):
        """
        Test the negation operator on elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_neg()

        .. SEEALSO::

            :class:`TestSuite`
        """
        tester = self._tester(**options)
        for x in tester.some_elements():
            y = -x
            tester.assertIs(y.parent(), self)
            tester.assertTrue((x+y).is_zero())
            tester.assertEqual(y.valuation(),x.valuation())
            tester.assertEqual(x.precision_absolute(),y.precision_absolute())
            tester.assertEqual(x.precision_relative(),y.precision_relative())
            tester.assertEqual(x.is_zero(),y.is_zero())
            tester.assertEqual(x.is_unit(),y.is_unit())

    def _test_log(self, **options):
        """
        Test the log operator on elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_log()

        .. SEEALSO::

            :class:`TestSuite`
        """
        tester = self._tester(**options)
        for x in tester.some_elements():
            if x.is_zero(): continue
            l = x.log(p_branch=0)
            tester.assertIs(l.parent(), self)
            tester.assertGreater(l.valuation(), 0)
            if self.is_capped_absolute() or self.is_capped_relative():
                tester.assertEqual(x.precision_relative(), l.precision_absolute())

        if self.is_capped_absolute() or self.is_capped_relative():
            # In the fixed modulus setting, rounding errors may occur
            elements = list(tester.some_elements())
            for x, y, b in some_tuples(elements, 3, tester._max_runs):
                if x.is_zero() or y.is_zero(): continue
                r1 = x.log(pi_branch=b) + y.log(pi_branch=b)
                r2 = (x*y).log(pi_branch=b)
                tester.assertEqual(r1, r2)

            p = self.prime()
            for x in tester.some_elements():
                if x.is_zero(): continue
                if p == 2:
                    a = 4 * x.unit_part()
                else:
                    a = p * x.unit_part()
                b = a.exp().log()
                c = (1+a).log().exp()
                tester.assertEqual(a, b)
                tester.assertEqual(1+a, c)

    def _test_teichmuller(self, **options):
        """
        Test Teichmuller lifts.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_teichmuller()

        .. SEEALSO::

            :class:`TestSuite`
        """
        tester = self._tester(**options)

        for x in tester.some_elements():
            try:
                y = self.teichmuller(x)
            except ValueError:
                tester.assertTrue(x.valuation() < 0 or x.precision_absolute()==0)
            else:
                try:
                    tester.assertEqual(x.residue(), y.residue())
                except (NotImplementedError, AttributeError):
                    pass
                tester.assertEqual(y**self.residue_field().order(), y)

    def _test_convert_residue_field(self, **options):
        r"""
        Test that conversion of residue field elements back to this ring works.

        INPUT:

         - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

        EXAMPLES::

            sage: Zp(3)._test_convert_residue_field()

        .. SEEALSO::

            :class:`TestSuite`
        """
        tester = self._tester(**options)

        for x in tester.some_elements():
            if x.valuation() < 0:
                continue
            if x.precision_absolute() <= 0:
                continue
            y = x.residue()
            z = self(y)
            tester.assertEqual(z.residue(), y)

    @cached_method
    def _log_unit_part_p(self):
        """
        Compute the logarithm of the unit-part of `p`.

        If `\pi` is the uniformizer in this ring, then we can uniquely write
        `p=\pi^e u` where `u` is a `\pi`-adic unit. This method computes the
        logarithm of `u`.

        This is a helper method for
        :meth:`sage.rings.padics.padic_generic_element.pAdicGenericElement.log`.

        TESTS::

            sage: R = Qp(3,5)
            sage: R._log_unit_part_p()
            O(3^5)

            sage: S.<x> = ZZ[]
            sage: W.<pi> = R.extension(x^3-3)
            sage: W._log_unit_part_p()
            O(pi^15)

            sage: W.<pi> = R.extension(x^3-3*x-3)
            sage: W._log_unit_part_p()
            2 + pi + 2*pi^2 + pi^4 + pi^5 + 2*pi^7 + 2*pi^8 + pi^9 + 2*pi^10 + pi^11 + pi^12 + 2*pi^14 + O(pi^15)

        """
        return self(self.prime()).unit_part().log()

    def frobenius_endomorphism(self, n=1):
        """
        INPUT:

        -  ``n`` -- an integer (default: 1)

        OUTPUT:

        The `n`-th power of the absolute arithmetic Frobenius
        endomorphism on this field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism on Unramified Extension ... lifting a |--> a^3 on the residue field
            sage: Frob(a) == a.frobenius()
            True

        We can specify a power::

            sage: K.frobenius_endomorphism(2)
            Frobenius endomorphism on Unramified Extension ... lifting a |--> a^(3^2) on the residue field

        The result is simplified if possible::

            sage: K.frobenius_endomorphism(6)
            Frobenius endomorphism on Unramified Extension ... lifting a |--> a^3 on the residue field
            sage: K.frobenius_endomorphism(5)
            Identity endomorphism of Unramified Extension ...

        Comparisons work::

            sage: K.frobenius_endomorphism(6) == Frob
            True
        """
        from .morphism import FrobeniusEndomorphism_padics
        return FrobeniusEndomorphism_padics(self, n)

    def _test_elements_eq_transitive(self, **options):
        """
        The operator ``==`` is not transitive for `p`-adic numbers. We disable
        the check of the category framework by overriding this method.

        EXAMPLES:

            sage: R = Zp(3)
            sage: R(3) == R(0,1)
            True
            sage: R(0,1) == R(6)
            True
            sage: R(3) == R(6)
            False
            sage: R._test_elements_eq_transitive()

        """
        pass

class ResidueReductionMap(Morphism):
    """
    Reduction map from a p-adic ring or field to its residue field or ring.

    These maps must be created using the :meth:`_create_` method in order
    to support categories correctly.

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic import ResidueReductionMap
        sage: R.<a> = Zq(125); k = R.residue_field()
        sage: f = ResidueReductionMap._create_(R, k); f
        Reduction morphism:
          From: Unramified Extension in a defined by x^3 + 3*x + 3 with capped relative precision 20 over 5-adic Ring
          To:   Finite Field in a0 of size 5^3
    """
    @staticmethod
    def _create_(R, k):
        """
        Initialization.  We have to implement this as a static method
        in order to call ``__make_element_class__``.

        INPUT:

        - ``R`` -- a `p`-adic ring or field.
        - ``k`` -- the residue field of ``R``, or a residue ring of ``R``.

        EXAMPLES::

            sage: f = Zmod(49).convert_map_from(Zp(7))
            sage: TestSuite(f).run()
            sage: K.<a> = Qq(125); k = K.residue_field(); f = k.convert_map_from(K)
            sage: TestSuite(f).run()
        """
        if R.is_field():
            from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
            cat = SetsWithPartialMaps()
        else:
            from sage.categories.rings import Rings
            cat = Rings()
        from sage.categories.homset import Hom
        kfield = R.residue_field()
        N = k.cardinality()
        q = kfield.cardinality()
        n = N.exact_log(q)
        if N != q**n:
            raise RuntimeError("N must be a power of q")
        H = Hom(R, k, cat)
        f = H.__make_element_class__(ResidueReductionMap)(H)
        f._n = n
        if kfield is k:
            f._field = True
        else:
            f._field = False
        return f

    def is_surjective(self):
        """
        The reduction map is surjective.

        EXAMPLES::

            sage: GF(7).convert_map_from(Qp(7)).is_surjective()
            True
        """
        return True

    def is_injective(self):
        """
        The reduction map is far from injective.

        EXAMPLES::

            sage: GF(5).convert_map_from(ZpCA(5)).is_injective()
            False
        """
        return False

    def _call_(self, x):
        """
        Evaluate this morphism.

        EXAMPLES::

            sage: R.<a> = Zq(125); k = R.residue_field()
            sage: f = k.convert_map_from(R)
            sage: f(15)
            0
            sage: f(1/(1+a))
            a0^2 + 4*a0 + 4

            sage: Zmod(121).convert_map_from(Qp(11))(3/11)
            Traceback (most recent call last):
            ...
            ValueError: element must have non-negative valuation in order to compute residue.
        """
        return x.residue(self._n, field=self._field, check_prec=self._field)

    def section(self):
        """
        Returns the section from the residue ring or field
        back to the p-adic ring or field.

        EXAMPLES::

            sage: GF(3).convert_map_from(Zp(3)).section()
            Lifting morphism:
              From: Finite Field of size 3
              To:   3-adic Ring with capped relative precision 20
        """
        return ResidueLiftingMap._create_(self.codomain(), self.domain())

    def _repr_type(self):
        """
        Type of morphism, for printing.

        EXAMPLES::

            sage: GF(3).convert_map_from(Zp(3))._repr_type()
            'Reduction'
        """
        return "Reduction"

    def _richcmp_(self, other, op):
        r"""
        Compare this element to ``other`` with respect to ``op``.

        EXAMPLES::

            sage: from sage.rings.padics.padic_generic import ResidueReductionMap
            sage: f = ResidueReductionMap._create_(Zp(3), GF(3))
            sage: g = ResidueReductionMap._create_(Zp(3), GF(3))
            sage: f is g
            False
            sage: f == g
            True
        """
        if type(self) != type(other):
            return NotImplemented
        return richcmp((self.domain(), self.codomain()), (other.domain(), other.codomain()), op)

# A class for the Teichmuller lift would also be reasonable....

class ResidueLiftingMap(Morphism):
    """
    Lifting map to a p-adic ring or field from its residue field or ring.

    These maps must be created using the :meth:`_create_` method in order
    to support categories correctly.

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic import ResidueLiftingMap
        sage: R.<a> = Zq(125); k = R.residue_field()
        sage: f = ResidueLiftingMap._create_(k, R); f
        Lifting morphism:
          From: Finite Field in a0 of size 5^3
          To:   Unramified Extension in a defined by x^3 + 3*x + 3 with capped relative precision 20 over 5-adic Ring
    """
    @staticmethod
    def _create_(k, R):
        """
        Initialization.  We have to implement this as a static method
        in order to call ``__make_element_class__``.

        INPUT:

        - ``k`` -- the residue field of ``R``, or a residue ring of ``R``.
        - ``R`` -- a `p`-adic ring or field.

        EXAMPLES::

            sage: f = Zp(3).convert_map_from(Zmod(81))
            sage: TestSuite(f).run()
        """
        from sage.categories.sets_cat import Sets
        from sage.categories.homset import Hom
        kfield = R.residue_field()
        N = k.cardinality()
        q = kfield.cardinality()
        n = N.exact_log(q)
        if N != q**n:
            raise RuntimeError("N must be a power of q")
        H = Hom(k, R, Sets())
        f = H.__make_element_class__(ResidueLiftingMap)(H)
        f._n = n
        return f

    def _call_(self, x):
        """
        Evaluate this morphism.

        EXAMPLES::

            sage: R.<a> = Zq(27); k = R.residue_field(); a0 = k.gen()
            sage: f = R.convert_map_from(k); f
            Lifting morphism:
              From: Finite Field in a0 of size 3^3
              To:   Unramified Extension in a defined by x^3 + 2*x + 1 with capped relative precision 20 over 3-adic Ring
            sage: f(a0 + 1)
            (a + 1) + O(3)

            sage: Zp(3)(Zmod(81)(0))
            O(3^4)
        """
        R = self.codomain()
        if R.degree() == 1:
            return R.element_class(R, x, self._n)
        elif R.f() == 1:
            return R([x], self._n)
        elif R.e() == 1:
            return R(x.polynomial().list(), self._n)
        else:
            raise NotImplementedError

    def _call_with_args(self, x, args=(), kwds={}):
        """
        Evaluate this morphism with extra arguments.

        EXAMPLES::

            sage: f = Zp(2).convert_map_from(Zmod(128))
            sage: f(7, 5) # indirect doctest
            1 + 2 + 2^2 + O(2^5)
        """
        R = self.codomain()
        if args:
            args = (min(args[0], self._n),) + args[1:]
        else:
            kwds['absprec'] = min(kwds.get('absprec', self._n), self._n)
        if R.degree() == 1:
            return R.element_class(R, x, *args, **kwds)
        elif R.f() == 1:
            return R([x], *args, **kwds)
        elif R.e() == 1:
            return R(x.polynomial().list(), *args, **kwds)
        else:
            raise NotImplementedError

    def _repr_type(self):
        """
        Type of morphism, for printing.

        EXAMPLES::

            sage: Zp(3).convert_map_from(GF(3))._repr_type()
            'Lifting'
        """
        return "Lifting"

    def _richcmp_(self, other, op):
        r"""
        Compare this element to ``other`` with respect to ``op``.

        EXAMPLES::

            sage: from sage.rings.padics.padic_generic import ResidueLiftingMap
            sage: f = ResidueLiftingMap._create_(GF(3), Zp(3))
            sage: g = ResidueLiftingMap._create_(GF(3), Zp(3))
            sage: f is g
            False
            sage: f == g
            True
        """
        if type(self) != type(other):
            return NotImplemented
        return richcmp((self.domain(), self.codomain()), (other.domain(), other.codomain()), op)

def local_print_mode(obj, print_options, pos = None, ram_name = None):
    r"""
    Context manager for safely temporarily changing the print_mode
    of a p-adic ring/field.

    EXAMPLES::

        sage: R = Zp(5)
        sage: R(45)
        4*5 + 5^2 + O(5^21)
        sage: with local_print_mode(R, 'val-unit'):
        ....:     print(R(45))
        5 * 9 + O(5^21)

    .. NOTE::

        For more documentation see ``localvars`` in parent_gens.pyx
    """
    if isinstance(print_options, str):
        print_options = {'mode': print_options}
    elif not isinstance(print_options, dict):
        raise TypeError("print_options must be a dictionary or a string")
    if pos is not None:
        print_options['pos'] = pos
    if ram_name is not None:
        print_options['ram_name'] = ram_name
    for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
        if option not in print_options:
            print_options[option] = obj._printer.dict()[option]
    return pAdicPrinter(obj, print_options)
