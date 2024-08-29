# sage.doctest: needs sage.rings.padics
r"""
`p`-adic Generic

A generic superclass for all `p`-adic parents.

AUTHORS:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
- Julian Rueth (2013-03-16): test methods for basic arithmetic
"""

# ****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#                               Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.misc import some_tuples
from copy import copy

from sage.structure.richcmp import richcmp
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.morphism import Morphism
from sage.categories.fields import Fields
from sage.rings.infinity import infinity
from .local_generic import LocalGeneric
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.rings.padics.precision_error import PrecisionError
from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import richcmp_not_equal


class pAdicGeneric(LocalGeneric):
    def __init__(self, base, p, prec, print_mode, names, element_class, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base`` -- base ring
        - ``p`` -- prime
        - ``print_mode`` -- dictionary of print options
        - ``names`` -- how to print the uniformizer
        - ``element_class`` -- the class for elements of this ring

        EXAMPLES::

            sage: R = Zp(17)  # indirect doctest
        """
        from sage.rings.padics.padic_printing import pAdicPrinter

        if category is None:
            if self.is_field():
                category = Fields()
            else:
                category = PrincipalIdealDomains()
            category = category.Metric().Complete()
        LocalGeneric.__init__(self, base, prec, names, element_class, category)
        self._printer = pAdicPrinter(self, print_mode)
        self._qth_roots_of_unity = [(1, Infinity)]

    def some_elements(self):
        r"""
        Return a list of elements in this ring.

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
        r"""
        Return a dictionary of print options, starting with ``self``'s
        print options but modified by the options in the dictionary
        ``print_mode``.

        INPUT:

        - ``print_mode`` -- dictionary with keys in

          * ``mode``
          * ``pos``
          * ``ram_name``
          * ``unram_name``
          * ``var_name``
          * ``max_ram_terms``
          * ``max_unram_terms``
          * ``max_terse_terms``
          * ``sep``
          * ``alphabet``

        EXAMPLES::

            sage: R = Zp(5)
            sage: R._modified_print_mode({'mode': 'bars'})['ram_name']
            '5'
        """
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name',
                       'max_ram_terms', 'max_unram_terms', 'max_terse_terms',
                       'sep', 'alphabet', 'show_prec']:
            if option not in print_mode:
                print_mode[option] = self._printer.dict()[option]
        return print_mode

    def ngens(self):
        r"""
        Return the number of generators of ``self``.

        We conventionally define this as 1: for base rings, we take a
        uniformizer as the generator; for extension rings, we take a
        root of the minimal polynomial defining the extension.

        EXAMPLES::

            sage: Zp(5).ngens()
            1
            sage: Zq(25,names='a').ngens()                                              # needs sage.libs.ntl
            1
        """
        return 1

    def gens(self):
        r"""
        Return a list of generators.

        EXAMPLES::

            sage: R = Zp(5); R.gens()
            [5 + O(5^21)]
            sage: Zq(25,names='a').gens()                                               # needs sage.libs.ntl
            [a + O(5^20)]
            sage: S.<x> = ZZ[]; f = x^5 + 25*x -5; W.<w> = R.ext(f); W.gens()           # needs sage.libs.ntl
            [w + O(w^101)]
        """
        return [self.gen()]

    def __richcmp__(self, other, op):
        r"""
        Rich comparison of ``self`` with ``other``.

        We consider two `p`-adic rings or fields to be equal if they are
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

    # def ngens(self):
    #     return 1

    # def gen(self, n = 0):
    #     if n != 0:
    #         raise IndexError, "only one generator"
    #     return self(self.prime())

    def print_mode(self):
        r"""
        Return the current print mode as a string.

        EXAMPLES::

            sage: R = Qp(7,5, 'capped-rel')
            sage: R.print_mode()
            'series'
        """
        return self._printer._print_mode()

    def characteristic(self):
        r"""
        Return the characteristic of ``self``, which is always 0.

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.characteristic()
            0
        """
        return Integer(0)

    def prime(self):
        r"""
        Return the prime, ie the characteristic of the residue field.

        OUTPUT: the characteristic of the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: R.prime()
            3
        """
        return self.prime_pow._prime()

    def uniformizer_pow(self, n):
        r"""
        Return `p^n`, as an element of ``self``.

        If ``n`` is infinity, returns 0.

        EXAMPLES::

            sage: R = Zp(3, 5, 'fixed-mod')
            sage: R.uniformizer_pow(3)
            3^3
            sage: R.uniformizer_pow(infinity)
            0
        """
        if n is infinity:
            return self(0)
        return self(self.prime_pow.pow_Integer_Integer(n))

    def _unram_print(self):
        r"""
        For printing.  Will be ``None`` if the unramified subextension
        of ``self`` is of degree 1 over `\ZZ_p` or `\QQ_p`.

        EXAMPLES::

            sage: Zp(5)._unram_print()
        """
        return None

    def residue_characteristic(self):
        r"""
        Return the prime, i.e., the characteristic of the residue field.

        OUTPUT: the characteristic of the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: R.residue_characteristic()
            3
        """
        return self.prime()

    def residue_class_field(self):
        r"""
        Return the residue class field.

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: k = R.residue_class_field()
            sage: k
            Finite Field of size 3
        """
        from sage.rings.finite_rings.finite_field_constructor import GF
        return GF(self.prime())

    def residue_field(self):
        r"""
        Return the residue class field.

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: k = R.residue_field()
            sage: k
            Finite Field of size 3
        """
        return self.residue_class_field()

    def residue_ring(self, n):
        r"""
        Return the quotient of the ring of integers by the ``n``-th
        power of the maximal ideal.

        EXAMPLES::

            sage: R = Zp(11)
            sage: R.residue_ring(3)
            Ring of integers modulo 1331
        """
        from sage.rings.finite_rings.integer_mod_ring import Zmod
        return Zmod(self.prime()**n)

    def residue_system(self):
        r"""
        Return a list of elements representing all the residue classes.

        EXAMPLES::

            sage: R = Zp(3, 5,'fixed-mod')
            sage: R.residue_system()
            [0, 1, 2]
        """
        return [self(i) for i in self.residue_class_field()]

    def _fraction_field_key(self, print_mode=None):
        r"""
        Change ``print_mode`` from a dictionary to a tuple, raising
        a deprecation warning if it is present.

        EXAMPLES::

            sage: Zp(5)._fraction_field_key()
            sage: Zp(5)._fraction_field_key({"pos":False})
            doctest:warning
            ...
            DeprecationWarning: Use the change method if you want to change print options in fraction_field()
            See https://github.com/sagemath/sage/issues/23227 for details.
            (('pos', False),)
        """
        if print_mode is not None:
            from sage.misc.superseded import deprecation
            deprecation(23227, "Use the change method if you want to change print options in fraction_field()")
            return tuple(sorted(print_mode.items()))

    @cached_method(key=_fraction_field_key)
    def fraction_field(self, print_mode=None):
        r"""
        Return the fraction field of this ring or field.

        For `\ZZ_p`, this is the `p`-adic field with the same options,
        and for extensions, it is just the extension of the fraction
        field of the base determined by the same polynomial.

        The fraction field of a capped absolute ring is capped relative,
        and that of a fixed modulus ring is floating point.

        INPUT:

        - ``print_mode`` -- (optional) a dictionary containing print options;
          defaults to the same options as this ring

        OUTPUT: the fraction field of this ring

        EXAMPLES::

            sage: R = Zp(5, print_mode='digits', show_prec=False)
            sage: K = R.fraction_field(); K(1/3)
            31313131313131313132
            sage: L = R.fraction_field({'max_ram_terms':4}); L(1/3)
            doctest:warning
            ...
            DeprecationWarning: Use the change method if you want to change print options in fraction_field()
            See https://github.com/sagemath/sage/issues/23227 for details.
            3132
            sage: U.<a> = Zq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)   # needs sage.libs.ntl
            sage: U.fraction_field()                                                    # needs sage.libs.ntl
            17-adic Unramified Extension Field in a defined by x^4 + 7*x^2 + 10*x + 3
            sage: U.fraction_field({"pos":False}) == U.fraction_field()                 # needs sage.libs.ntl
            False

        TESTS::

            sage: R = ZpLC(2); R
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/23505 for details.
            2-adic Ring with lattice-cap precision
            sage: K = R.fraction_field(); K
            2-adic Field with lattice-cap precision

            sage: K = QpLC(2); K2 = K.fraction_field({'mode':'terse'})
            sage: K2 is K
            False
            sage: K = QpLC(2, label='test'); K
            2-adic Field with lattice-cap precision (label: test)
            sage: K.fraction_field()
            2-adic Field with lattice-cap precision (label: test)
            sage: K.fraction_field({'mode':'series'}) is K
            True
        """
        if self.is_field() and print_mode is None:
            return self
        if print_mode is None:
            return self.change(field=True)
        else:
            return self.change(field=True, **print_mode)

    def integer_ring(self, print_mode=None):
        r"""
        Return the ring of integers of this ring or field.

        For `\QQ_p`, this is the `p`-adic ring with the same options,
        and for extensions, it is just the extension of the ring
        of integers of the base determined by the same polynomial.

        INPUT:

        - ``print_mode`` -- (optional) a dictionary containing print options;
          defaults to the same options as this ring

        OUTPUT: the ring of elements of this field with nonnegative valuation

        EXAMPLES::

            sage: K = Qp(5, print_mode='digits', show_prec=False)
            sage: R = K.integer_ring(); R(1/3)
            31313131313131313132
            sage: S = K.integer_ring({'max_ram_terms':4}); S(1/3)
            doctest:warning
            ...
            DeprecationWarning: Use the change method if you want to change print options in integer_ring()
            See https://github.com/sagemath/sage/issues/23227 for details.
            3132
            sage: U.<a> = Qq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)   # needs sage.libs.ntl
            sage: U.integer_ring()                                                      # needs sage.libs.ntl
            17-adic Unramified Extension Ring in a defined by x^4 + 7*x^2 + 10*x + 3
            sage: U.fraction_field({"print_mode":"terse"}) == U.fraction_field()        # needs sage.libs.ntl
            doctest:warning
            ...
            DeprecationWarning: Use the change method if you want to change print options in fraction_field()
            See https://github.com/sagemath/sage/issues/23227 for details.
            False

        TESTS::

            sage: K = QpLC(2); K
            2-adic Field with lattice-cap precision
            sage: R = K.integer_ring(); R
            2-adic Ring with lattice-cap precision

            sage: R = ZpLC(2); R2 = R.integer_ring({'mode':'terse'})
            sage: R2 is R
            False
            sage: R = ZpLC(2, label='test'); R
            2-adic Ring with lattice-cap precision (label: test)
            sage: R.integer_ring()
            2-adic Ring with lattice-cap precision (label: test)
            sage: R.integer_ring({'mode':'series'}) is R
            True

        The `secure` attribute for relaxed type is preserved::

            sage: K = QpER(5, secure=True)                                              # needs sage.libs.flint
            sage: K.integer_ring().is_secure()
            True
        """
        # Currently does not support fields with non integral defining
        # polynomials.  This should change when the padic_general_extension
        # framework gets worked out.
        if not self.is_field() and print_mode is None:
            return self
        if print_mode is None:
            return self.change(field=False, check=False)
        else:
            from sage.misc.superseded import deprecation
            deprecation(23227, "Use the change method if you want to change print options in integer_ring()")
            return self.change(field=False, **print_mode)

    def teichmuller(self, x, prec=None):
        r"""
        Return the Teichmüller representative of ``x``.

        - ``x`` -- something that can be cast into ``self``

        OUTPUT: the Teichmüller lift of ``x``

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
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9

            sage: # needs sage.libs.ntl
            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W.teichmuller(3); y
            3 + 3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + 3*w^15
             + 2*w^16 + 3*w^17 + w^18 + 3*w^19 + 3*w^20 + 2*w^21 + 2*w^22
             + 3*w^23 + 4*w^24 + O(w^25)
            sage: y^5 == y
            True
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = A.teichmuller(1 + 2*a - a^2); b
            (4*a^2 + 2*a + 1) + 2*a*5 + (3*a^2 + 1)*5^2 + (a + 4)*5^3
             + (a^2 + a + 1)*5^4 + O(5^5)
            sage: b^125 == b
            True

        We check that :issue:`23736` is resolved::

            sage: # needs sage.libs.ntl
            sage: R.teichmuller(GF(5)(2))
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + O(5^5)

        AUTHORS:

        - Initial version: David Roe
        - Quadratic time version: Kiran Kedlaya <kedlaya@math.mit.edu> (2007-03-27)
        """
        ans = self(x) if prec is None else self(x, prec)
        # Since Teichmüller representatives are defined at infinite precision,
        # we can lift to precision prec, as long as the absolute precision of ans is positive.
        if ans.precision_absolute() <= 0:
            raise ValueError("not enough precision to determine Teichmuller representative")
        if ans.valuation() > 0:
            return self(0) if prec is None else self(0, prec)
        ans = ans.lift_to_precision(prec)
        if ans is x:
            ans = copy(ans)
        ans._teichmuller_set_unsafe()
        return ans

    def teichmuller_system(self):
        r"""
        Return a set of Teichmüller representatives for the invertible
        elements of `\ZZ / p\ZZ`.

        OUTPUT:

        A list of Teichmüller representatives for the invertible elements
        of `\ZZ / p\ZZ`.

        EXAMPLES::

            sage: R = Zp(3, 5,'fixed-mod', 'terse')
            sage: R.teichmuller_system()
            [1, 242]

        Check that :issue:`20457` is fixed::

            sage: F.<a> = Qq(5^2,6)                                                     # needs sage.libs.ntl
            sage: F.teichmuller_system()[3]                                             # needs sage.libs.ntl
            (2*a + 2) + (4*a + 1)*5 + 4*5^2 + (2*a + 1)*5^3 + (4*a + 1)*5^4 + (2*a + 3)*5^5 + O(5^6)

        .. NOTE::

            Should this return 0 as well?
        """
        R = self.residue_class_field()
        prec = self.precision_cap()
        return [self.teichmuller(self(i).lift_to_precision(prec))
                for i in R if i != 0]

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

    def extension(self, modulus, prec=None, names=None, print_mode=None, implementation='FLINT', **kwds):
        r"""
        Create an extension of this `p`-adic ring.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: k = Qp(5)
            sage: R.<x> = k[]
            sage: l.<w> = k.extension(x^2 - 5); l
            5-adic Eisenstein Extension Field in w defined by x^2 - 5
            sage: F = list(Qp(19)['x'](cyclotomic_polynomial(5)).factor())[0][0]
            sage: L = Qp(19).extension(F, names='a'); L
            19-adic Unramified Extension Field in a defined by x^2 + 8751674996211859573806383*x + 1
        """
        if isinstance(modulus, list):
            if len(modulus) == 0:
                return self
            else:
                return self.extension(modulus[-1], prec=prec[-1],
                                      names=names[-1],
                                      implementation=implementation[-1],
                                      print_mode=print_mode, **kwds).extension(
                                          modulus[:-1], prec=prec[:-1],
                                          names=names[:-1],
                                          implementation=implementation[:-1],
                                          print_mode=print_mode, **kwds)
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
        return ExtensionFactory(base=self, modulus=modulus, prec=prec, names=names, check=True, implementation=implementation, **print_mode)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        r"""
        Check whether the given images and map on the base ring determine a
        valid homomorphism to the codomain.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<x> = ZZ[]
            sage: K.<a> = Qq(25, modulus=x^2-2)
            sage: L.<b> = Qq(625, modulus=x^4-2)
            sage: K._is_valid_homomorphism_(L, [b^2])
            True
            sage: L._is_valid_homomorphism_(L, [b^3])
            False
            sage: z = L(-1).sqrt()
            sage: L._is_valid_homomorphism_(L, [z*b])
            True
            sage: L._is_valid_homomorphism_(L, [-b])
            True
            sage: W.<w> = K.extension(x^2 - 5)
            sage: cc = K.hom([-a])
            sage: W._is_valid_homomorphism_(W, [w], base_map=cc)
            True
            sage: W._is_valid_homomorphism_(W, [-w], base_map=cc)
            True
            sage: W._is_valid_homomorphism_(W, [w+1])
            False
        """
        K = self.base_ring()
        if base_map is None and not codomain.has_coerce_map_from(K):
            return False
        if len(im_gens) != 1:
            raise ValueError("Wrong number of generators")
        if self is K:
            # Qp or Zp, so either base_map is not None or there's a coercion to the codomain
            # We check that the im_gens has the right length and value
            return im_gens[0] == codomain(self.prime())
        # Now we're an extension.  We check that the defining polynomial maps to zero
        f = self.modulus()
        if base_map is not None:
            f = f.change_ring(base_map)
        return f(im_gens[0]) == 0

    def _test_add(self, **options):
        r"""
        Test addition of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
            if self.is_lattice_prec():
                tester.assertGreaterEqual(z.precision_absolute(), zprec)
            elif not self.is_floating_point():
                tester.assertEqual(z.precision_absolute(), zprec)
            tester.assertGreaterEqual(z.valuation(), min(x.valuation(),y.valuation()))
            if x.valuation() != y.valuation():
                tester.assertEqual(z.valuation(), min(x.valuation(),y.valuation()))
            tester.assertTrue(y.is_equal_to(z-x,zprec))
            tester.assertTrue(x.is_equal_to(z-y,zprec))

    def _test_sub(self, **options):
        r"""
        Test subtraction on elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
            if self.is_lattice_prec():
                tester.assertGreaterEqual(z.precision_absolute(), zprec)
            elif not self.is_floating_point():
                tester.assertEqual(z.precision_absolute(), zprec)
            tester.assertGreaterEqual(z.valuation(), min(x.valuation(),y.valuation()))
            if x.valuation() != y.valuation():
                tester.assertEqual(z.valuation(), min(x.valuation(),y.valuation()))
            tester.assertTrue((-y).is_equal_to(z - x,zprec))
            tester.assertTrue(x.is_equal_to(z + y,zprec))

    def _test_invert(self, **options):
        """
        Test multiplicative inversion of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
                if not self.is_fixed_mod():
                    tester.assertTrue(x.is_zero())
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
        r"""
        Test multiplication of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
        r"""
        Test division of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
                if self.is_fixed_mod():
                    tester.assertFalse(y.is_unit())
                else:
                    tester.assertTrue(y.is_zero())
            else:
                try:
                    xx = z*y
                except ZeroDivisionError:
                    tester.assertTrue(self.is_floating_point() and (z.is_zero() or y.is_zero()))
                else:
                    tester.assertFalse(y.is_zero())
                    tester.assertIs(z.parent(), self if self.is_fixed_mod() else self.fraction_field())
                    # The following might be false if there is an absolute cap
                    # tester.assertEqual(z.precision_relative(), min(x.precision_relative(), y.precision_relative()))
                    if not x.is_zero():
                        tester.assertEqual(z.valuation(), x.valuation() - y.valuation())
                    tester.assertEqual(xx, x)

    def _test_neg(self, **options):
        r"""
        Test the negation operator on elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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

    def _test_shift(self, **options):
        """
        Test the shift operator on elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

        EXAMPLES::

            sage: Zp(3)._test_shift()

        .. SEEALSO::

            :class:`TestSuite`
        """
        tester = self._tester(**options)
        if self.is_relaxed():
            cap = self.default_prec()
        else:
            cap = self.precision_cap()
        k = self.residue_field()
        for v in range(min(cap,10)):
            if self.is_capped_absolute() or self.is_fixed_mod():
                prec = cap - v
            else:
                prec = cap
            b = self.uniformizer_pow(v)
            for x in tester.some_elements():
                y = (x << v) >> v
                if x._is_exact_zero() or self.is_field():
                    tester.assertEqual(x, y)
                else:
                    tester.assertTrue(x.is_equal_to(y, prec))
                y = (x >> v) << v
                if x._is_exact_zero() or self.is_field():
                    tester.assertEqual(x, y)
                else:
                    for i in range(min(v,prec)):
                        tester.assertEqual(k(y.expansion(i)), 0)
                    for i in range(v,prec):
                        tester.assertEqual(y.expansion(i), x.expansion(i))
                    xx = y + (x % b)
                    tester.assertTrue(xx.is_equal_to(x,prec))

    def _test_log(self, **options):
        r"""
        Test the log operator on elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

        EXAMPLES::

            sage: Zp(3)._test_log()

        .. SEEALSO::

            :class:`TestSuite`
        """
        tester = self._tester(**options)
        for x in tester.some_elements():
            if x.is_zero():
                continue
            try:
                l = x.log(p_branch=0)
                tester.assertIs(l.parent(), self)
            except ValueError:
                l = x.log(p_branch=0, change_frac=True)
            if self.is_capped_absolute() or self.is_capped_relative():
                if self.absolute_e() == 1:
                    tester.assertEqual(l.precision_absolute(), x.precision_relative())
                else:
                    tester.assertLessEqual(l.precision_absolute(), x.precision_relative())

        if self.is_capped_absolute() or self.is_capped_relative():
            # In the fixed modulus setting, rounding errors may occur
            for x, y, b in tester.some_elements(repeat=3):
                if (x*y).is_zero():
                    continue
                r1 = x.log(pi_branch=b) + y.log(pi_branch=b)
                r2 = (x*y).log(pi_branch=b)
                tester.assertEqual(r1, r2)

            p = self.prime()
            for x in tester.some_elements():
                if x.is_zero():
                    continue
                if p == 2:
                    a = 4 * x.unit_part()
                else:
                    a = p * x.unit_part()
                b = a.exp().log()
                c = (1+a).log().exp()
                tester.assertEqual(a, b)
                tester.assertEqual(1+a, c)

    def _test_teichmuller(self, **options):
        r"""
        Test Teichmüller lifts.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
                tester.assertTrue(x.valuation() < 0 or x.precision_absolute() == 0)
            else:
                try:
                    tester.assertEqual(x.residue(), y.residue())
                except (NotImplementedError, AttributeError):
                    pass
                if self.is_relaxed():
                    tester.assertTrue(y.is_equal_at_precision(y**self.residue_field().order(), self.default_prec()))
                else:
                    tester.assertEqual(y**self.residue_field().order(), y)

    def _test_convert_residue_field(self, **options):
        r"""
        Test that conversion of residue field elements back to this ring works.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
        r"""
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
            sage: W.<pi> = R.extension(x^3-3)                                           # needs sage.libs.ntl
            sage: W._log_unit_part_p()                                                  # needs sage.libs.ntl
            O(pi^15)

            sage: W.<pi> = R.extension(x^3-3*x-3)                                       # needs sage.libs.ntl
            sage: W._log_unit_part_p()                                                  # needs sage.libs.ntl
            2 + pi + 2*pi^2 + pi^4 + pi^5 + 2*pi^7 + 2*pi^8 + pi^9 + 2*pi^10 + pi^11 + pi^12 + 2*pi^14 + O(pi^15)
        """
        return self(self.prime()).unit_part().log()

    def frobenius_endomorphism(self, n=1):
        r"""
        Return the `n`-th power of the absolute arithmetic Frobeninus
        endomorphism on this field.

        INPUT:

        - ``n`` -- integer (default: 1)

        EXAMPLES::

            sage: K.<a> = Qq(3^5)                                                       # needs sage.libs.ntl
            sage: Frob = K.frobenius_endomorphism(); Frob                               # needs sage.libs.ntl
            Frobenius endomorphism on 3-adic Unramified Extension
            ... lifting a |--> a^3 on the residue field
            sage: Frob(a) == a.frobenius()                                              # needs sage.libs.ntl
            True

        We can specify a power::

            sage: K.frobenius_endomorphism(2)                                           # needs sage.libs.ntl
            Frobenius endomorphism on 3-adic Unramified Extension
            ... lifting a |--> a^(3^2) on the residue field

        The result is simplified if possible::

            sage: K.frobenius_endomorphism(6)                                           # needs sage.libs.ntl
            Frobenius endomorphism on 3-adic Unramified Extension
            ... lifting a |--> a^3 on the residue field
            sage: K.frobenius_endomorphism(5)                                           # needs sage.libs.ntl
            Identity endomorphism of 3-adic Unramified Extension ...

        Comparisons work::

            sage: K.frobenius_endomorphism(6) == Frob                                   # needs sage.libs.ntl
            True
        """
        from .morphism import FrobeniusEndomorphism_padics
        return FrobeniusEndomorphism_padics(self, n)

    def _test_elements_eq_transitive(self, **options):
        r"""
        The operator ``==`` is not transitive for `p`-adic numbers.

        We disable the check of the category framework by overriding
        this method.

        EXAMPLES::

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

    def valuation(self):
        r"""
        Return the `p`-adic valuation on this ring.

        OUTPUT:

        A valuation that is normalized such that the rational prime `p` has
        valuation 1.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^3 - 3)
            sage: v = L.valuation(); v
            3-adic valuation
            sage: v(3)
            1
            sage: L(3).valuation()
            3

        The normalization is chosen such that the valuation restricts to the
        valuation on the base ring::

            sage: v(3) == K.valuation()(3)                                              # needs sage.libs.ntl
            True
            sage: v.restriction(K) == K.valuation()                                     # needs sage.libs.ntl
            True

        .. SEEALSO::

            :meth:`NumberField_generic.valuation() <sage.rings.number_field.number_field.NumberField_generic.valuation>`,
            :meth:`Order.valuation() <sage.rings.number_field.order.Order.valuation>`
        """
        from sage.rings.padics.padic_valuation import pAdicValuation
        return pAdicValuation(self)

    def _primitive_qth_root_of_unity(self, exponent):
        r"""
        Compute the ``p^exponent``-th roots of unity in this ring.

        INPUT:

        - ``exponent`` -- integer or ``Infinity``

        OUTPUT:

        A triple ``(zeta, n, nextzeta)`` where

        - ``zeta`` is a generator of the group of ``p^exponent``-th
          roots of unity in this ring, and
        - ``p^n`` is the order of ``zeta``.
        - ``nextzeta`` is the result of ``zeta._inverse_pth_root()``
          if ``n`` is positive and ``None`` otherwise

        TESTS::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(2^3, 5)
            sage: S.<x> = K[]
            sage: L.<pi> = K.extension(x^2 + 2*x + 2)
            sage: zeta = L.primitive_root_of_unity(); zeta  # indirect doctest
            a + a*pi + pi^2 + a*pi^4 + a*pi^5 + a^2*pi^8 + a^2*pi^9 + O(pi^10)
            sage: zeta.parent() is L
            True
        """
        n = len(self._qth_roots_of_unity)

        # We check if the result is cached
        if exponent < n-1:
            return self._qth_roots_of_unity[exponent][0], exponent, self._qth_roots_of_unity[exponent+1]
        zeta, accuracy = self._qth_roots_of_unity[-1]
        if accuracy is not Infinity:
            return self._qth_roots_of_unity[-2][0], n-2, (zeta, accuracy)

        # It is not, so we compute it
        while accuracy is Infinity and n <= exponent + 1:
            self._qth_roots_of_unity[-1] = (self(zeta), Infinity)  # to avoid multiple conversions
            if n == 1:  # case of pth root of unity
                p = self.prime()
                e = self.absolute_e()
                k = self.residue_field()
                if e % (p-1) != 0:
                    # No pth root of unity in this ring
                    zeta = accuracy = None
                else:
                    rho = -k(self(p).expansion(e))
                    try:
                        r = rho.nth_root(p-1)
                    except ValueError:
                        # No pth root of unity in this ring
                        zeta = accuracy = None
                    else:
                        # We compute a primitive pth root of unity
                        m = e // (p-1)
                        prec = self.precision_cap() + e * (1 + m.valuation(p))
                        ring = self.change(prec=prec)
                        zeta = 1 + (ring(r).lift_to_precision() << m)
                        curprec = m*p + 1
                        while curprec < prec:
                            curprec -= e
                            curprec = min(2*curprec + e, p*curprec)
                            zeta = zeta.lift_to_precision(min(prec,curprec))
                            zeta += zeta * (1 - zeta**p) // p
            else:
                zeta, accuracy = zeta._inverse_pth_root()
                assert accuracy is not None
            self._qth_roots_of_unity.append((zeta, accuracy))
            n += 1
        return self._qth_roots_of_unity[-2][0], n-2, self._qth_roots_of_unity[-1]

    def primitive_root_of_unity(self, n=None, order=False):
        r"""
        Return a generator of the group of ``n``-th roots of unity
        in this ring.

        INPUT:

        - ``n`` -- integer or ``None`` (default: ``None``)

        - ``order`` -- boolean (default: ``False``)

        OUTPUT:

        A generator of the group of ``n``-th roots of unity.
        If ``n`` is ``None``, a generator of the full group of roots
        of unity is returned.

        If ``order`` is ``True``, the order of the above group is
        returned as well.

        EXAMPLES::

            sage: R = Zp(5, 10)
            sage: zeta = R.primitive_root_of_unity(); zeta
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: zeta == R.teichmuller(2)
            True

        Now we consider an example with non trivial ``p``-th roots of unity::

            sage: # needs sage.libs.ntl
            sage: W = Zp(3, 2)
            sage: S.<x> = W[]
            sage: R.<pi> = W.extension((x+1)^6 + (x+1)^3 + 1)
            sage: zeta, order = R.primitive_root_of_unity(order=True)
            sage: zeta
            2 + 2*pi + 2*pi^3 + 2*pi^7 + 2*pi^8 + 2*pi^9 + pi^11 + O(pi^12)
            sage: order
            18
            sage: zeta.multiplicative_order()
            18
            sage: zeta, order = R.primitive_root_of_unity(24, order=True)
            sage: zeta
            2 + pi^3 + 2*pi^7 + 2*pi^8 + 2*pi^10 + 2*pi^11 + O(pi^12)
            sage: order   # equal to gcd(18,24)
            6
            sage: zeta.multiplicative_order()
            6
        """
        p = self.prime()
        k = self.residue_field()
        prec = self.precision_cap()
        c = k.cardinality()

        # We compute a primitive qth root of unity
        # where q is the highest power of p dividing exponent
        if n is None:
            qthzeta, s, _ = self._primitive_qth_root_of_unity(Infinity)
            m = c - 1
        else:
            qthzeta, s, _ = self._primitive_qth_root_of_unity(n.valuation(p))
            m = n.gcd(c - 1)
        qthzeta = self(qthzeta)

        # We now compute a primitive mth root of qthzeta
        if m == 1:
            zeta = qthzeta
        else:
            zeta = self(k.multiplicative_generator() ** ((c-1) // m))
            invm = self(1/m)
            curprec = 1
            while curprec < prec:
                curprec *= 2
                zeta = zeta.lift_to_precision(min(prec,curprec))
                zeta += invm * zeta * (1 - qthzeta*zeta**m)

        if order:
            return zeta, m * p**s
        else:
            return zeta

    def roots_of_unity(self, n=None):
        r"""
        Return all the ``n``-th roots of unity in this ring.

        INPUT:

        - ``n`` -- integer or ``None`` (default: ``None``); if
          ``None``, the full group of roots of unity is returned

        EXAMPLES::

            sage: R = Zp(5, 10)
            sage: roots = R.roots_of_unity(); roots
            [1 + O(5^10),
             2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10),
             4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10),
             3 + 3*5 + 2*5^2 + 3*5^3 + 5^4 + 2*5^6 + 5^7 + 4*5^8 + 5^9 + O(5^10)]

            sage: R.roots_of_unity(10)
            [1 + O(5^10),
             4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)]

        In this case, the roots of unity are the Teichmüller representatives::

            sage: R.teichmuller_system()
            [1 + O(5^10),
             2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10),
             3 + 3*5 + 2*5^2 + 3*5^3 + 5^4 + 2*5^6 + 5^7 + 4*5^8 + 5^9 + O(5^10),
             4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)]

        In general, there might be more roots of unity (it happens when the ring has non
        trivial ``p``-th roots of unity)::

            sage: # needs sage.libs.ntl
            sage: W.<a> = Zq(3^2, 2)
            sage: S.<x> = W[]
            sage: R.<pi> = W.extension((x+1)^2 + (x+1) + 1)
            sage: roots = R.roots_of_unity(); roots
            [1 + O(pi^4),
             a + 2*a*pi + 2*a*pi^2 + a*pi^3 + O(pi^4),
             ...
             1 + pi + O(pi^4),
             a + a*pi^2 + 2*a*pi^3 + O(pi^4),
             ...
             1 + 2*pi + pi^2 + O(pi^4),
             a + a*pi + a*pi^2 + O(pi^4),
             ...]
            sage: len(roots)
            24

        We check that the logarithm of each root of unity vanishes::

            sage: # needs sage.libs.ntl
            sage: for root in roots:
            ....:     if root.log() != 0:
            ....:         raise ValueError
        """
        zeta, order = self.primitive_root_of_unity(n, order=True)
        return [ zeta**i for i in range(order) ]

    def _roots_univariate_polynomial(self, P, ring, multiplicities, algorithm, secure=False):
        r"""
        Return the roots of ``P`` in the ring ``ring``.

        INPUT:

        - ``P`` -- a polynomial defined over this ring

        - ``ring`` -- a ring into which this ring coerces

        - ``multiplicities`` -- boolean (default: ``True``);
          whether we have to return the multiplicities of each
          root or not

        - ``algorithm`` -- ``'pari'``, ``'sage'`` or ``None`` (default:
          ``None``); Sage provides an implementation for any extension of
          `\QQ_p` whereas only roots of polynomials over `\QQ_p` is implemented
          in Pari; the default is ``'pari'`` if ``ring`` is `\ZZ_p` or `\QQ_p`,
          ``'sage'`` otherwise.

        - ``secure`` -- boolean (default: ``False``)

        .. NOTE::

            When ``secure`` is ``True``, this method raises an error when
            the precision on the input polynomial is not enough to determine
            the number of roots in the ground field. This happens when two
            roots cannot be separated.
            A typical example is the polynomial

            .. MATH::

                 (1 + O(p^10))*X^2 + O(p^10)*X + O(p^10).

            Indeed its discriminant might be any `p`-adic integer divisible
            by `p^{10}` (resp. `p^{11}` when `p=2`) and so can be as well
            zero, a square and a non-square.
            In the first case, the polynomial has one double root; in the
            second case, it has two roots; in the third case, it has no
            root in `\QQ_p`.

            When ``secure`` is ``False``, this method assumes that two
            inseparable roots actually collapse. In the above example,
            it then answers that the given polynomial has a double root
            `O(p^5)`.

        This keyword is ignored when using the ``pari`` algorithm.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: A = Zp(3, prec=10, print_mode='terse')
            sage: S.<x> = A[]
            sage: P = x^2 - 7
            sage: P.roots()
            [(30793 + O(3^10), 1), (28256 + O(3^10), 1)]
            sage: P.roots(multiplicities=False)
            [30793 + O(3^10), 28256 + O(3^10)]

        We compare with the result given by the method
        :meth:`sage.rings.padics.padic_generic_element.square_root`::

            sage: A(7).square_root(all=True)                                            # needs sage.libs.ntl
            [30793 + O(3^10), 28256 + O(3^10)]

        Here is another example::

            sage: P = x * (x-1) * (x-2) * (x-3) * (x-4)                                 # needs sage.libs.ntl
            sage: P.roots(multiplicities=False)                                         # needs sage.libs.ntl
            [39370 + O(3^10),
             19684 + O(3^10),
             2 + O(3^10),
             3 + O(3^10),
             O(3^10)]

        The result is not quite what we expected.
        In fact, the roots are correct but the precision is not::

            sage: [ root.add_bigoh(9) for root in P.roots(multiplicities=False) ]       # needs sage.libs.ntl
            [4 + O(3^9),
             1 + O(3^9),
             2 + O(3^9),
             3 + O(3^9),
             O(3^9)]

        This is due to the fact that we are using ``'pari'`` which does not
        track precision (it can only compute `p`-adic roots of exact polynomials).
        If we are switching to ``'sage'`` then the precision on the result
        becomes correct (but the computation is much slower)::

            sage: P.roots(multiplicities=False, algorithm='sage')                       # needs sage.geometry.polyhedron sage.libs.ntl
            [0,
             3 + O(3^11),
             1 + O(3^9),
             4 + O(3^9),
             2 + O(3^9)]

        We check that the keyword ``secure`` works as explained above::

            sage: P = x^2 + O(3^10)*x + O(3^10)                                         # needs sage.libs.ntl
            sage: P.roots(algorithm='sage')                                             # needs sage.geometry.polyhedron sage.libs.ntl
            [(O(3^5), 2)]
            sage: P.roots(algorithm='sage', secure=True)                                # needs sage.libs.ntl
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to determine the number of roots

        An example over an extension::

            sage: B.<b> = Zq(3^3, prec=10, print_mode='terse')                          # needs sage.libs.ntl
            sage: P = B.modulus()                                                       # needs sage.libs.ntl

        We check that `P` has no root in `A`::

            sage: P.roots()                                                             # needs sage.libs.ntl
            []

        but that it has roots in `B`::

            sage: P.roots(B)                                                            # needs sage.geometry.polyhedron sage.libs.ntl
            [(35149 + 57730*b + 41124*b^2 + O(3^10), 1),
             (23900 + 1318*b + 17925*b^2 + O(3^10), 1),
             (b + O(3^10), 1)]

        We check further that the other roots are the conjugates
        of ``b`` under Frobenius::

            sage: b.frobenius()                                                         # needs sage.libs.ntl
            23900 + 1318*b + 17925*b^2 + O(3^10)
            sage: b.frobenius().frobenius()                                             # needs sage.libs.ntl
            35149 + 57730*b + 41124*b^2 + O(3^10)

        Root finding works over ramified extensions also::

            sage: # needs sage.libs.ntl
            sage: E = x^3 - 3*x + 3*b
            sage: C.<pi> = B.extension(E)
            sage: E.roots(C)                                                            # needs sage.geometry.polyhedron
            [(pi + O(pi^30), 1)]
            sage: S.<x> = C[]
            sage: P = prod(x - (pi+i) for i in range(5))
            sage: P.roots()                                                             # needs sage.geometry.polyhedron
            [(pi + O(pi^29), 1),
             (3 + pi + O(pi^29), 1),
             (1 + pi + O(pi^27), 1),
             (4 + pi + O(pi^27), 1),
             (2 + pi + O(pi^30), 1)]

        TESTS::

            sage: S(0).roots()                                                          # needs sage.libs.ntl
            Traceback (most recent call last):
            ...
            ArithmeticError: factorization of 0 is not defined
        """
        if P.is_zero():
            raise ArithmeticError("factorization of 0 is not defined")
        if ring is None:
            ring = self
        if algorithm is None:
            try:
                return self._roots_univariate_polynomial(P, ring, multiplicities, "pari", secure)
            except (NotImplementedError, PrecisionError):
                return self._roots_univariate_polynomial(P, ring, multiplicities, "sage", secure)
        elif algorithm == "pari":
            P = P.change_ring(ring)
            try:
                # note that P.factor() calls pari
                return P._roots_from_factorization(P.factor(), multiplicities)
            except (AttributeError, TypeError):
                raise NotImplementedError("root finding for this polynomial is not implemented in pari")
        elif algorithm == "sage":
            if ring.is_field():
                roots = P.change_ring(ring)._roots(secure, -Infinity, None)
            else:
                K = ring.fraction_field()
                roots = P.change_ring(K)._roots(secure, 0, None)
            if multiplicities:
                return [ (ring(root), m) for (root, m) in roots ]
            else:
                return [ ring(root) for (root, m) in roots ]


class ResidueReductionMap(Morphism):
    r"""
    Reduction map from a `p`-adic ring or field to its residue field or ring.

    These maps must be created using the :meth:`_create_` method in order
    to support categories correctly.

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic import ResidueReductionMap
        sage: R.<a> = Zq(125); k = R.residue_field()                                    # needs sage.libs.ntl
        sage: f = ResidueReductionMap._create_(R, k); f                                 # needs sage.libs.ntl
        Reduction morphism:
          From: 5-adic Unramified Extension Ring in a defined by x^3 + 3*x + 3
          To:   Finite Field in a0 of size 5^3
    """
    @staticmethod
    def _create_(R, k):
        r"""
        Initialization.  We have to implement this as a static method
        in order to call ``__make_element_class__``.

        INPUT:

        - ``R`` -- a `p`-adic ring or field
        - ``k`` -- the residue field of ``R``, or a residue ring of ``R``

        EXAMPLES::

            sage: f = Zmod(49).convert_map_from(Zp(7))                                  # needs sage.rings.finite_rings
            sage: TestSuite(f).run()                                                    # needs sage.rings.finite_rings
            sage: K.<a> = Qq(125); k = K.residue_field(); f = k.convert_map_from(K)     # needs sage.libs.ntl
            sage: TestSuite(f).run()                                                    # needs sage.rings.finite_rings
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
        r"""
        The reduction map is surjective.

        EXAMPLES::

            sage: GF(7).convert_map_from(Qp(7)).is_surjective()                         # needs sage.rings.finite_rings
            True
        """
        return True

    def is_injective(self):
        r"""
        The reduction map is far from injective.

        EXAMPLES::

            sage: GF(5).convert_map_from(ZpCA(5)).is_injective()                        # needs sage.rings.finite_rings
            False
        """
        return False

    def _call_(self, x):
        r"""
        Evaluate this morphism.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<a> = Zq(125); k = R.residue_field()
            sage: f = k.convert_map_from(R)
            sage: f(15)
            0
            sage: f(1/(1+a))
            a0^2 + 4*a0 + 4

            sage: Zmod(121).convert_map_from(Qp(11))(3/11)                              # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue
        """
        return x.residue(self._n, field=self._field, check_prec=self._field)

    def section(self):
        r"""
        Return the section from the residue ring or field
        back to the `p`-adic ring or field.

        EXAMPLES::

            sage: GF(3).convert_map_from(Zp(3)).section()                               # needs sage.rings.finite_rings
            Lifting morphism:
              From: Finite Field of size 3
              To:   3-adic Ring with capped relative precision 20
        """
        return ResidueLiftingMap._create_(self.codomain(), self.domain())

    def _repr_type(self):
        r"""
        Type of morphism, for printing.

        EXAMPLES::

            sage: GF(3).convert_map_from(Zp(3))._repr_type()                            # needs sage.rings.finite_rings
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
        if type(self) is not type(other):
            return NotImplemented
        return richcmp((self.domain(), self.codomain()), (other.domain(), other.codomain()), op)

# A class for the Teichmüller lift would also be reasonable....

class ResidueLiftingMap(Morphism):
    r"""
    Lifting map to a `p`-adic ring or field from its residue field or ring.

    These maps must be created using the :meth:`_create_` method in order
    to support categories correctly.

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic import ResidueLiftingMap
        sage: R.<a> = Zq(125); k = R.residue_field()                                    # needs sage.libs.ntl
        sage: f = ResidueLiftingMap._create_(k, R); f                                   # needs sage.libs.ntl
        Lifting morphism:
          From: Finite Field in a0 of size 5^3
          To:   5-adic Unramified Extension Ring in a defined by x^3 + 3*x + 3
    """
    @staticmethod
    def _create_(k, R):
        r"""
        Initialization.  We have to implement this as a static method
        in order to call ``__make_element_class__``.

        INPUT:

        - ``k`` -- the residue field of ``R``, or a residue ring of ``R``
        - ``R`` -- a `p`-adic ring or field

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
        r"""
        Evaluate this morphism.

        EXAMPLES::

            sage: R.<a> = Zq(27); k = R.residue_field(); a0 = k.gen()                   # needs sage.libs.ntl
            sage: f = R.convert_map_from(k); f                                          # needs sage.libs.ntl
            Lifting morphism:
              From: Finite Field in a0 of size 3^3
              To:   3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
            sage: f(a0 + 1)                                                             # needs sage.libs.ntl
            (a + 1) + O(3)

            sage: Zp(3)(Zmod(81)(0))
            O(3^4)
        """
        R = self.codomain()
        K = R.maximal_unramified_subextension()
        if self._n == 1 or K is R:
            unram_n = self._n
            if K.absolute_degree() == 1:
                lift = K._element_constructor_(x, unram_n)
            else:
                lift = K(x.polynomial().list(), unram_n)
            return R(lift, self._n)
        else:
            #unram_n = (self._n - 1) // R.absolute_e() + 1
            raise NotImplementedError

    def _call_with_args(self, x, args=(), kwds={}):
        r"""
        Evaluate this morphism with extra arguments.

        EXAMPLES::

            sage: f = Zp(2).convert_map_from(Zmod(128))
            sage: f(7, 5)  # indirect doctest
            1 + 2 + 2^2 + O(2^5)
        """
        R = self.codomain()
        kwds = dict(kwds) # we're changing it
        if args:
            args = (min(args[0], self._n),) + args[1:]
            absprec = args[0]
        else:
            absprec = kwds['absprec'] = min(kwds.get('absprec', self._n), self._n)
        K = R.maximal_unramified_subextension()
        if absprec == 1 or K is R:
            if K.absolute_degree() == 1:
                lift = K._element_constructor_(x, *args, **kwds)
            else:
                lift = K(x.polynomial().list(), *args, **kwds)
            return R(lift, *args, **kwds)
        else:
            raise NotImplementedError

    def _repr_type(self):
        r"""
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

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.padics.padic_generic import ResidueLiftingMap
            sage: f = ResidueLiftingMap._create_(GF(3), Zp(3))
            sage: g = ResidueLiftingMap._create_(GF(3), Zp(3))
            sage: f is g
            False
            sage: f == g
            True
        """
        if type(self) is not type(other):
            return NotImplemented
        return richcmp((self.domain(), self.codomain()), (other.domain(), other.codomain()), op)

def local_print_mode(obj, print_options, pos=None, ram_name=None):
    r"""
    Context manager for safely temporarily changing the print_mode
    of a `p`-adic ring/field.

    EXAMPLES::

        sage: R = Zp(5)
        sage: R(45)
        4*5 + 5^2 + O(5^21)
        sage: with local_print_mode(R, 'val-unit'):
        ....:     print(R(45))
        5 * 9 + O(5^21)

    .. NOTE::

        For more documentation see :class:`sage.structure.parent_gens.localvars`.
    """
    from sage.rings.padics.padic_printing import pAdicPrinter

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
