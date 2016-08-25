"""
Eisenstein Extension Generic

This file implements the shared functionality for Eisenstein extensions.

AUTHORS:

- David Roe
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .padic_extension_generic import pAdicExtensionGeneric
from sage.rings.infinity import infinity
from sage.misc.latex import latex
from sage.rings.integer import Integer

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

    def _repr_(self, do_latex = False):
        """
        Returns a print representation of this extension.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B #indirect doctest
            Eisenstein Extension of 7-adic Ring with capped relative precision 10 in t defined by (1 + O(7^10))*x^2 + (O(7^11))*x + (7 + O(7^11))
        """
        if do_latex:
            return "Eisenstein Extension of %s in %s defined by %s"%(latex(self.ground_ring()), self.latex_name(), latex(self.modulus()))
        else:
            return "Eisenstein Extension of %s in %s defined by %s"%(self.ground_ring(), self.variable_name(), self.modulus())

    def ramification_index(self, K = None):
        """
        Returns the ramification index of self over K, or over the
        ground ring if K is None.

        The ramification index is the index of the image of the
        valuation map on K in the image of the valuation map on self
        (both normalized so that the valuation of p is 1).

        INPUT:

        - self -- an Eisenstein extension
        - K -- a subring of self (default None -> self.ground_ring())

        OUTPUT:

        - The ramification index of the extension self/K

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.ramification_index()
            2
        """
        if K is None or K is self.ground_ring():
            return self.modulus().degree()
        elif K is self:
            return 1
        else:
            raise NotImplementedError

    def inertia_degree(self, K = None):
        """
        Returns the inertia degree of self over K, or the ground ring
        if K is None.

        The inertia degree is the degree of the extension of residue
        fields induced by this extensions.  Since Eisenstein
        extensions are totally ramified, this will be 1 for K=None.

        INPUT:

        - self -- an Eisenstein extension
        - K -- a subring of self (default None -> self.ground_ring())

        OUTPUT:

        - The degree of the induced extensions of residue fields.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.inertia_degree()
            1
        """
        if K is None or K is self.ground_ring():
            return Integer(1)
        elif K is self:
            return Integer(1)
        else:
            raise NotImplementedError

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
