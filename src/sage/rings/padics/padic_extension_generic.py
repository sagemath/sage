"""
p-Adic Extension Generic

A common superclass for all extensions of Qp and Zp.

AUTHORS:

- David Roe
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .padic_generic import pAdicGeneric
from .padic_base_generic import pAdicBaseGeneric
from functools import reduce

class pAdicExtensionGeneric(pAdicGeneric):
    def __init__(self, poly, prec, print_mode, names, element_class):
        """
        Initialization

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f) #indirect doctest
        """
        #type checking done in factory
        self._given_poly = poly
        R = poly.base_ring()
        # We'll deal with the different names better later.
        # Using a tuple here is mostly needed for more general extensions
        # (ie not eisenstein or unramified)
        print_mode['unram_name'] = names[2]
        print_mode['ram_name'] = names[3]
        print_mode['var_name'] = names[0]
        names = names[0]
        pAdicGeneric.__init__(self, R, R.prime(), prec, print_mode, names, element_class)
        self._populate_coercion_lists_(coerce_list=[R], element_constructor=element_class)

#     def __reduce__(self):
#         """
#         For pickling.

#         This function is provided because prime_pow needs to be set before _printer, so the standard unpickling fails.
#         """
#         from sage.rings.padics.factory import ExtensionFactory
#         return ExtensionFactory, (self.base_ring(), self._pre_poly, self.precision_cap(), self.print_mode(), None, self.variable_name())

    def _coerce_map_from_(self, R):
        """
        Returns True if there is a coercion map from R to self.

        EXAMPLES::

            sage: R = Zp(5); S.<x> = ZZ[]; f = x^5 + 25*x - 5; W.<w> = R.ext(f)
            sage: L = W.fraction_field()
            sage: w + L(w) #indirect doctest
            2*w + O(w^101)
            sage: w + R(5,2)
            w + w^5 + O(w^10)
        """
        # Far more functionality needs to be added here later.
        if isinstance(R, pAdicExtensionGeneric) and R.fraction_field() is self:
            if self._implementation == 'NTL':
                return True
            elif R._prec_type() == 'capped-abs':
                from sage.rings.padics.qadic_flint_CA import pAdicCoercion_CA_frac_field as coerce_map
            elif R._prec_type() == 'capped-rel':
                from sage.rings.padics.qadic_flint_CR import pAdicCoercion_CR_frac_field as coerce_map
            return coerce_map(R, self)

    def __cmp__(self, other):
        """
        Returns 0 if self == other, and 1 or -1 otherwise.

        We consider two p-adic rings or fields to be equal if they are equal mathematically, and also have the same precision cap and printing parameters.

        EXAMPLES::

            sage: R.<a> = Qq(27)
            sage: S.<a> = Qq(27,print_mode='val-unit')
            sage: R == S
            False
            sage: S.<a> = Qq(27,type='capped-rel')
            sage: R == S
            True
            sage: R is S
            True
        """
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        groundcmp = self.ground_ring().__cmp__(other.ground_ring())
        if groundcmp != 0:
            return groundcmp
        c = cmp(self.defining_polynomial(), other.defining_polynomial())
        if c != 0:
            return c
        c = cmp(self.precision_cap(), other.precision_cap())
        if c != 0:
            return c
        return self._printer.cmp_modes(other._printer)

    #def absolute_discriminant(self):
    #    raise NotImplementedError

    #def discriminant(self):
    #    raise NotImplementedError

    #def is_abelian(self):
    #    raise NotImplementedError

    #def is_normal(self):
    #    raise NotImplementedError

    def degree(self):
        """
        Returns the degree of this extension.

        EXAMPLES::

            sage: R.<a> = Zq(125); R.degree()
            3
            sage: R = Zp(5); S.<x> = ZZ[]; f = x^5 - 25*x^3 + 5; W.<w> = R.ext(f)
            sage: W.degree()
            5
        """
        return self._given_poly.degree()

    def defining_polynomial(self):
        """
        Returns the polynomial defining this extension.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W.defining_polynomial()
            (1 + O(5^5))*x^5 + (O(5^6))*x^4 + (3*5^2 + O(5^6))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^6))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
        """
        return self._given_poly

    def modulus(self):
        """
        Returns the polynomial defining this extension.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W.modulus()
            (1 + O(5^5))*x^5 + (O(5^6))*x^4 + (3*5^2 + O(5^6))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^6))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
        """
        return self._given_poly

    def ground_ring(self):
        """
        Returns the ring of which this ring is an extension.

        EXAMPLE::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W.ground_ring()
            5-adic Ring with capped relative precision 5
        """
        return self._given_poly.base_ring()

    def ground_ring_of_tower(self):
        """
        Returns the p-adic base ring of which this is ultimately an
        extension.

        Currently this function is identical to ground_ring(), since
        relative extensions have not yet been implemented.

        EXAMPLES::

            sage: Qq(27,30,names='a').ground_ring_of_tower()
            3-adic Field with capped relative precision 30
        """
        if isinstance(self.ground_ring(), pAdicBaseGeneric):
            return self.ground_ring()
        else:
            return self.ground_ring().ground_ring_of_tower()

    #def is_isomorphic(self, ring):
    #    raise NotImplementedError

    def polynomial_ring(self):
        """
        Returns the polynomial ring of which this is a quotient.

        EXAMPLES::

            sage: Qq(27,30,names='a').polynomial_ring()
            Univariate Polynomial Ring in x over 3-adic Field with capped relative precision 30
        """
        return self._given_poly.parent()

    #def teichmuller(self, x, prec = None):
    #    if prec is None:
    #        prec = self.precision_cap()
    #    x = self(x, prec)
    #    if x.valuation() > 0:
    #        return self(0)
    #    q = self.residue_class_field().order()
    #    u = 1 / self(1 - q, prec)
    #    delta = u * (1 - x ** (q - 1))
    #    xnew = x - x*delta*(1 - q * delta)
    #    while x != xnew:
    #        x = xnew
    #        delta = u*(1-x**(q-1))
    #        xnew = x - x*delta*(1-q*delta)
    #    return x

    def fraction_field(self, print_mode=None):
        r"""
        Returns the fraction field of this extension, which is just
        the extension of base.fraction_field() determined by the
        same polynomial.

        INPUT:

        - print_mode -- a dictionary containing print options.
          Defaults to the same options as this ring.

        OUTPUT:

        - the fraction field of self.

        EXAMPLES::

            sage: U.<a> = Zq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)
            sage: U.fraction_field()
            Unramified Extension of 17-adic Field with capped relative precision 6 in a defined by (1 + O(17^6))*x^4 + (O(17^6))*x^3 + (7 + O(17^6))*x^2 + (10 + O(17^6))*x + (3 + O(17^6))
            sage: U.fraction_field({"pos":False}) == U.fraction_field()
            False
        """
        if self.is_field() and print_mode is None:
            return self
        print_mode = self._modified_print_mode(print_mode)
        ground_mode = print_mode.copy()
        # We don't want to confuse the ground ring with different names.
        ground_mode['ram_name'] = None
        ground_mode['unram_name'] = None
        K = self.ground_ring().fraction_field(ground_mode)
        #we don't want to set the print options due to the ground ring since
        #different extension fields (with different options) can share the same ground ring.
        if self.is_lazy():
            return K.extension(self._pre_poly, prec = self.precision_cap(), halt = self.halting_parameter(), res_name = self.residue_field().variable_name(), print_mode=print_mode, implementation=self._implementation)
        else:
            return K.extension(self._pre_poly, prec = self.precision_cap(), res_name = self.residue_field().variable_name(), print_mode=print_mode, implementation=self._implementation)

    def integer_ring(self, print_mode=None):
        r"""
        Returns the ring of integers of self, which is just the
        extension of base.integer_ring() determined by the same
        polynomial.

        INPUT:

            - print_mode -- a dictionary containing print options.
              Defaults to the same options as this ring.

        OUTPUT:

            - the ring of elements of self with nonnegative valuation.

        EXAMPLES::

            sage: U.<a> = Qq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)
            sage: U.integer_ring()
            Unramified Extension of 17-adic Ring with capped relative precision 6 in a defined by (1 + O(17^6))*x^4 + (O(17^6))*x^3 + (7 + O(17^6))*x^2 + (10 + O(17^6))*x + (3 + O(17^6))
            sage: U.fraction_field({"pos":False}) == U.fraction_field()
            False
        """
        #Currently does not support fields with non integral defining polynomials.  This should change when the padic_general_extension framework gets worked out.
        if not self.is_field() and print_mode is None:
            return self
        print_mode = self._modified_print_mode(print_mode)
        ground_mode = print_mode.copy()
        # We don't want to confuse the ground ring with different names.
        ground_mode['ram_name'] = None
        ground_mode['unram_name'] = None
        K = self.ground_ring().integer_ring(ground_mode)
        #we don't want to set the print options due to the ground ring since
        #different extension fields (with different options) can share the same ground ring.
        if self.is_lazy():
            return K.extension(self._pre_poly, prec = self.precision_cap(), halt = self.halting_parameter(), res_name = self.residue_field().variable_name(), print_mode=print_mode)
        else:
            return K.extension(self._pre_poly, prec = self.precision_cap(), res_name = self.residue_field().variable_name(), print_mode=print_mode)

    #def hasGNB(self):
    #    raise NotImplementedError

    def random_element(self):
        """
        Returns a random element of self.

        This is done by picking a random element of the ground ring
        self.degree() times, then treating those elements as
        coefficients of a polynomial in self.gen().

        EXAMPLES::

            sage: R.<a> = Zq(125, 5); R.random_element()
            3*a + (2*a + 1)*5 + (4*a^2 + 3*a + 4)*5^2 + (a^2 + 2*a)*5^3 + (a + 2)*5^4 + O(5^5)
            sage: R = Zp(5,3); S.<x> = ZZ[]; f = x^5 + 25*x^2 - 5; W.<w> = R.ext(f)
            sage: W.random_element()
            3 + 4*w + 3*w^2 + w^3 + 4*w^4 + w^5 + w^6 + 3*w^7 + w^8 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + 4*w^14 + O(w^15)
        """
        return reduce(lambda x,y: x+y,
                      [self.ground_ring().random_element() * self.gen()**i for i in
                           range(self.modulus().degree())],
                      0)

    #def unit_group(self):
    #    raise NotImplementedError

    #def unit_group_gens(self):
    #    raise NotImplementedError

    #def principal_unit_group(self):
    #    raise NotImplementedError

    #def zeta(self, n = None):
    #    raise NotImplementedError

    #def zeta_order(self):
    #    raise NotImplementedError

