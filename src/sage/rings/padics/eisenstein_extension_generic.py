# sage_setup: distribution = sagemath-pari
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
        Initialize ``self``.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]                                                           # needs sage.libs.ntl
            sage: B.<t> = A.ext(x^2+7)  # indirect doctest                              # needs sage.libs.ntl sage.rings.padics
        """
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        #self._precompute()

    def _extension_type(self):
        """
        Return the type (``Unramified``, ``Eisenstein``) of this
        extension as a string, if any.

        Used for printing.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)                                                       # needs sage.libs.ntl
            sage: K._extension_type()                                                   # needs sage.libs.ntl
            'Unramified'

            sage: x = polygen(ZZ, 'x')
            sage: L.<pi> = Qp(5).extension(x^2 - 5)                                     # needs sage.libs.ntl
            sage: L._extension_type()                                                   # needs sage.libs.ntl
            'Eisenstein'
        """
        return "Eisenstein"

    def absolute_e(self):
        """
        Return the absolute ramification index of this ring or field.

        EXAMPLES::

            sage: K.<a> = Qq(3^5)                                                       # needs sage.libs.ntl
            sage: K.absolute_e()                                                        # needs sage.libs.ntl
            1

            sage: x = polygen(ZZ, 'x')
            sage: L.<pi> = Qp(3).extension(x^2 - 3)                                     # needs sage.libs.ntl
            sage: L.absolute_e()                                                        # needs sage.libs.ntl
            2
        """
        return self.modulus().degree() * self.base_ring().absolute_e()

    def inertia_subring(self):
        """
        Return the inertia subring.

        Since an Eisenstein extension is totally ramified, this is
        just the ground field.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]                                                           # needs sage.libs.ntl
            sage: B.<t> = A.ext(x^2 + 7)                                                # needs sage.libs.ntl
            sage: B.inertia_subring()                                                   # needs sage.libs.ntl
            7-adic Ring with capped relative precision 10
        """
        return self.ground_ring()

    def residue_class_field(self):
        """
        Return the residue class field.

        INPUT:

        - ``self`` -- a `p`-adic ring

        OUTPUT: the residue field

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]                                                           # needs sage.libs.ntl
            sage: B.<t> = A.ext(x^2 + 7)                                                # needs sage.libs.ntl
            sage: B.residue_class_field()                                               # needs sage.libs.ntl
            Finite Field of size 7
        """
        return self.ground_ring().residue_class_field()

    def residue_ring(self, n):
        """
        Return the quotient of the ring of integers by the `n`-th power of its maximal ideal.

        EXAMPLES::

            sage: S.<x> = ZZ[]
            sage: W.<w> = Zp(5).extension(x^2 - 5)                                      # needs sage.libs.ntl
            sage: W.residue_ring(1)                                                     # needs sage.libs.ntl
            Ring of integers modulo 5

        The following requires implementing more general Artinian rings::

            sage: W.residue_ring(2)                                                     # needs sage.libs.ntl
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
    #    Returns the Galois group of ``self``'s fraction field over Qp.
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
        Return a generator for ``self`` as an extension of its ground ring.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]                                                           # needs sage.libs.ntl
            sage: B.<t> = A.ext(x^2 + 7)                                                # needs sage.libs.ntl
            sage: B.gen()                                                               # needs sage.libs.ntl
            t + O(t^21)
        """
        if n != 0:
            raise IndexError("only one generator")
        return self([0,1])

    def uniformizer_pow(self, n):
        """
        Return the `n`-th power of the uniformizer of ``self`` (as an
        element of ``self``).

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]                                                           # needs sage.libs.ntl
            sage: B.<t> = A.ext(x^2 + 7)                                                # needs sage.libs.ntl
            sage: B.uniformizer_pow(5)                                                  # needs sage.libs.ntl
            t^5 + O(t^25)
        """
        if n is infinity:
            return self(0)
        else:
            return self(1) << n

    def uniformizer(self):
        """
        Return the uniformizer of ``self``, i.e., a generator for the unique
        maximal ideal.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]                                                           # needs sage.libs.ntl
            sage: B.<t> = A.ext(x^2 + 7)                                                # needs sage.libs.ntl
            sage: B.uniformizer()                                                       # needs sage.libs.ntl
            t + O(t^21)
        """
        return self.gen()

    def _uniformizer_print(self):
        """
        Return a string representation of how the uniformizer of self
        prints.  Mainly for internal use.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]                                                           # needs sage.libs.ntl
            sage: B.<t> = A.ext(x^2+7)                                                  # needs sage.libs.ntl
            sage: B._uniformizer_print()                                                # needs sage.libs.ntl
            't'
        """
        return self.variable_name()

#     def has_pth_root(self):
#         raise NotImplementedError

#     def has_root_of_unity(self, n):
#         raise NotImplementedError
