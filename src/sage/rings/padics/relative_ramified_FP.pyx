# sage.doctest: needs sage.libs.flint
include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "FP_template.pxi"

cdef class RelativeRamifiedFloatingPointElement(FPElement):
    def _poly_rep(self):
        """
        Return the underlying polynomial representation of the element
        (which is used for computations).

        For debugging and printing purpose.

        EXAMPLES::

            sage: K.<a> = QqFP(125)
            sage: S.<x> = PolynomialRing(K)
            sage: W.<w> = K.extension(x^3 - 25*x^2 - 5*a*x + 5)
            sage: w._poly_rep()
            x
            sage: W(5)._poly_rep()
            5

        The coefficients of P are floating point `p`-adics::

            sage: P = W.random_element()._poly_rep()
            sage: ring = P.parent().base_ring()
            sage: ring
            5-adic Unramified Extension Field in a defined by x^3 + 3*x + 3
            sage: ring._prec_type()
            'floating-point'
        """
        return self.unit.parent()(ccoefficients(self.unit, self.ordp, self.prime_pow.ram_prec_cap, self.prime_pow))
