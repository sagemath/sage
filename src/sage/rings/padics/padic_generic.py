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

from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.fields import Fields
from sage.rings.infinity import infinity
from .local_generic import LocalGeneric
from sage.rings.ring import PrincipalIdealDomain
from sage.rings.integer import Integer
from sage.rings.padics.padic_printing import pAdicPrinter
from sage.rings.padics.precision_error import PrecisionError
from sage.misc.cachefunc import cached_method


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
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
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

    def __cmp__(self, other):
        """
        Returns 0 if self == other, and 1 or -1 otherwise.

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
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        if self.prime() < other.prime():
            return -1
        elif self.prime() > other.prime():
            return 1
        try:
            if self.halting_parameter() < other.halting_parameter():
                return -1
            elif self.halting_parameter() > other.halting_parameter():
                return 1
        except AttributeError:
            pass
        if self.precision_cap() < other.precision_cap():
            return -1
        elif self.precision_cap() > other.precision_cap():
            return 1
        return self._printer.cmp_modes(other._printer)

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

        AUTHORS:

        - Initial version: David Roe
        - Quadratic time version: Kiran Kedlaya <kedlaya@math.mit.edu> (3/27/07)
        """
        if prec is None:
            prec = self.precision_cap()
        else:
            prec = min(Integer(prec), self.precision_cap())
        ans = self(x, prec)
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

    def extension(self, modulus, prec = None, names = None, print_mode = None, halt = None, implementation='FLINT', **kwds):
        """
        Create an extension of this p-adic ring.

        EXAMPLES::

            sage: k = Qp(5)
            sage: R.<x> = k[]
            sage: l.<w> = k.extension(x^2-5); l
            Eisenstein Extension of 5-adic Field with capped relative precision 20 in w defined by (1 + O(5^20))*x^2 + (O(5^21))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + 4*5^20 + O(5^21))

            sage: F = list(Qp(19)['x'](cyclotomic_polynomial(5)).factor())[0][0]
            sage: L = Qp(19).extension(F, names='a')
            sage: L
            Unramified Extension of 19-adic Field with capped relative precision 20 in a defined by (1 + O(19^20))*x^2 + (5 + 2*19 + 10*19^2 + 14*19^3 + 7*19^4 + 13*19^5 + 5*19^6 + 12*19^7 + 8*19^8 + 4*19^9 + 14*19^10 + 6*19^11 + 5*19^12 + 13*19^13 + 16*19^14 + 4*19^15 + 17*19^16 + 8*19^18 + 4*19^19 + O(19^20))*x + (1 + O(19^20))
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
        return ExtensionFactory(base=self, premodulus=modulus, prec=prec, halt=halt, names=names, check = True, implementation=implementation, **print_mode)

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
            tester.assertEqual(z.precision_absolute(), min(x.precision_absolute(), y.precision_absolute()))
            tester.assertGreaterEqual(z.valuation(), min(x.valuation(),y.valuation()))
            if x.valuation() != y.valuation():
                tester.assertEqual(z.valuation(), min(x.valuation(),y.valuation()))
            tester.assertEqual(z - x, y)
            tester.assertEqual(z - y, x)

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
            tester.assertEqual(z.precision_absolute(), min(x.precision_absolute(), y.precision_absolute()))
            tester.assertGreaterEqual(z.valuation(), min(x.valuation(),y.valuation()))
            if x.valuation() != y.valuation():
                tester.assertEqual(z.valuation(), min(x.valuation(),y.valuation()))
            tester.assertEqual(z - x, -y)
            tester.assertEqual(z + y, x)

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
                e = y * x

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
                tester.assertFalse(y.is_zero())
                tester.assertIs(z.parent(), self if self.is_fixed_mod() else self.fraction_field())
                tester.assertEqual(z.precision_relative(), min(x.precision_relative(), y.precision_relative()))
                tester.assertEqual(z.valuation(), x.valuation() - y.valuation())

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

    @cached_method
    def _exp_p(self):
        """
        Compute the exponential of `p`.

        This is a helper method for
        :meth:`sage.rings.padics.padic_generic_element.pAdicGenericElement.exp`.

        TESTS::

            sage: R = Qp(3, 5)
            sage: R._exp_p()
            1 + 3 + 3^2 + 2*3^3 + 2*3^4 + O(3^5)

            sage: S.<x> = ZZ[]
            sage: W.<pi> = R.extension(x^3-3)
            sage: W._exp_p()
            1 + pi^3 + pi^6 + 2*pi^9 + 2*pi^12 + O(pi^15)
            sage: R._exp_p() == W._exp_p()
            True

            sage: W.<pi> = R.extension(x^3-3*x-3)
            sage: W._exp_p()
            1 + pi^3 + 2*pi^4 + pi^5 + pi^7 + pi^9 + pi^10 + 2*pi^11 + pi^12 + pi^13 + 2*pi^14 + O(pi^15)
            sage: R._exp_p() == W._exp_p()
            True

        """
        p = self.prime()
        if p == 2:
            # the exponential of 2 does not exist, so we compute the
            # exponential of 4 instead.
            p = 4
        return self(p)._exp(self.precision_cap())

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
            Frobenius endomorphism on Unramified Extension of 3-adic Field ... lifting a |--> a^3 on the residue field
            sage: Frob(a) == a.frobenius()
            True

        We can specify a power:: 

            sage: K.frobenius_endomorphism(2)
            Frobenius endomorphism on Unramified Extension of 3-adic Field ... lifting a |--> a^(3^2) on the residue field

        The result is simplified if possible::

            sage: K.frobenius_endomorphism(6)
            Frobenius endomorphism on Unramified Extension of 3-adic Field ... lifting a |--> a^3 on the residue field
            sage: K.frobenius_endomorphism(5)
            Identity endomorphism of Unramified Extension of 3-adic Field ...

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

    NOTES::

        For more documentation see localvars in parent_gens.pyx
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
