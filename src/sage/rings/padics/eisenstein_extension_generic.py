"""
Eisenstein Extension Generic

This file implements the shared functionality for Eisenstein extensions.

AUTHORS:

- David Roe
"""

# ****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .padic_extension_generic import pAdicExtensionGeneric
from sage.rings.infinity import infinity


class EisensteinExtensionGeneric(pAdicExtensionGeneric):
    def __init__(self, poly, prec, print_mode, names, element_class):
        """
        Initializes self.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7) #indirect doctest
        """
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        #self._precompute()

    def _extension_type(self):
        """
        Return the type (``Unramified``, ``Eisenstein``) of this 
        extension as a string, if any.

        Used for printing.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: K._extension_type()
            'Unramified'

            sage: L.<pi> = Qp(5).extension(x^2 - 5)
            sage: L._extension_type()
            'Eisenstein'
        """
        return "Eisenstein"

    def absolute_e(self):
        """
        Return the absolute ramification index of this ring or field

        EXAMPLES::

            sage: K.<a> = Qq(3^5)
            sage: K.absolute_e()
            1

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.absolute_e()
            2
        """
        return self.modulus().degree() * self.base_ring().absolute_e()

    def inertia_subring(self):
        """
        Returns the inertia subring.

        Since an Eisenstein extension is totally ramified, this is
        just the ground field.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.inertia_subring()
            7-adic Ring with capped relative precision 10
        """
        return self.ground_ring()

    def residue_class_field(self):
        """
        Returns the residue class field.

        INPUT:

        - self -- a p-adic ring

        OUTPUT:

        - the residue field

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.residue_class_field()
            Finite Field of size 7
        """
        return self.ground_ring().residue_class_field()

    def residue_ring(self, n):
        """
        Return the quotient of the ring of integers by the nth power of its maximal ideal.

        EXAMPLES::

            sage: S.<x> = ZZ[]
            sage: W.<w> = Zp(5).extension(x^2 - 5)
            sage: W.residue_ring(1)
            Ring of integers modulo 5

        The following requires implementing more general Artinian rings::

            sage: W.residue_ring(2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if n == 1:
            return self.ground_ring().residue_ring(1)
        else:
            raise NotImplementedError

    #def discriminant(self, K=None):
    #    if K is self:
    #        return 1
    #    else:
    #        raise NotImplementedError

    #def automorphisms(self):
    #    raise NotImplementedError

    #def galois_group(self):
    #    r"""
    #    Returns the Galois group of self's fraction field over Qp.
    #    """
    #    ##
    #    ## If K is a number field, then K.galois_group() can return
    #    ## other variants, i.e. via Pari or KASH. We could consider
    #    ## doing this.
    #    ##
    #    raise NotImplementedError

    #def is_abelian(self):
    #    raise NotImplementedError

    #def is_normal(self):
    #    raise NotImplementedError

    def gen(self, n=0):
        """
        Returns a generator for self as an extension of its ground ring.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.gen()
            t + O(t^21)
        """
        if n != 0:
            raise IndexError("only one generator")
        return self([0,1])

    def uniformizer_pow(self, n):
        """
        Returns the nth power of the uniformizer of self (as an
        element of self).

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.uniformizer_pow(5)
            t^5 + O(t^25)
        """
        if n is infinity:
            return self(0)
        else:
            return self(1) << n

    def uniformizer(self):
        """
        Returns the uniformizer of self, ie a generator for the unique
        maximal ideal.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.uniformizer()
            t + O(t^21)
        """
        return self.gen()

    def _uniformizer_print(self):
        """
        Returns a string representation of how the uniformizer of self
        prints.  Mainly for internal use.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B._uniformizer_print()
            't'
        """
        return self.variable_name()

#     def has_pth_root(self):
#         raise NotImplementedError

#     def has_root_of_unity(self, n):
#         raise NotImplementedError
