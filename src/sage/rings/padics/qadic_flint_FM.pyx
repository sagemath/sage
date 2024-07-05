include "sage/libs/linkages/padics/fmpz_poly_unram.pxi"
include "sage/libs/linkages/padics/unram_shared.pxi"
include "FM_template.pxi"

cdef class PowComputer_(PowComputer_flint_unram):
    """
    A PowComputer for a fixed-modulus unramified ring.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        """
        Initialization.

        EXAMPLES::

            sage: R.<a> = ZqFM(125)
            sage: type(R.prime_pow)
            <class 'sage.rings.padics.qadic_flint_FM.PowComputer_'>
            sage: R.prime_pow._prec_type
            'fixed-mod'
        """
        self._prec_type = 'fixed-mod'
        PowComputer_flint_unram.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

cdef class qAdicFixedModElement(FMElement):
    frobenius = frobenius_unram
    trace = trace_unram
    norm = norm_unram

    def matrix_mod_pn(self):
        r"""
        Return the matrix of right multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for this
        extension field.  Thus the *rows* of this matrix give the
        images of each of the `x^i`.  The entries of the matrices are
        IntegerMod elements, defined modulo ``p^(self.absprec() / e)``.

        EXAMPLES::

            sage: R.<a> = ZqFM(5^5,5)
            sage: b = (5 + 15*a)^3
            sage: b.matrix_mod_pn()
            [ 125 1125  250  250    0]
            [   0  125 1125  250  250]
            [2375 2125  125 1125  250]
            [2375 1375 2125  125 1125]
            [2875 1000 1375 2125  125]

            sage: M = R(0,3).matrix_mod_pn(); M == 0
            True
            sage: M.base_ring()
            Ring of integers modulo 3125
        """
        return cmatrix_mod_pn(self.value, self.prime_pow.prec_cap, 0, self.prime_pow)

    def _flint_rep(self, var='x'):
        """
        Replacement for _ntl_rep for use in printing and debugging.

        EXAMPLES::

            sage: R.<a> = ZqFM(27, 4)
            sage: (1+a).inverse_of_unit()._flint_rep()
            41*x^2 + 40*x + 42
            sage: (1+a)*(41*a^2+40*a+42)
            1
        """
        return self.prime_pow._new_fmpz_poly(self.value, var)

    def _flint_rep_abs(self, var='x'):
        """
        Replacement for _ntl_rep_abs for use in printing and debugging.

        EXAMPLES::

            sage: R.<a> = ZqFM(27, 4)
            sage: (3+3*a)._flint_rep_abs()
            (3*x + 3, 0)
        """
        return self._flint_rep(var), Integer(0)

    def _modp_rep(self, use_smallest_mode=False, return_list=True):
        r"""
        Return the element with the same reduction mod p that can be expressed
        with coefficients between 0 and p-1.  The absolute precision will be maximal.

        This method is used in printing and computing `p`-adic expansions.

        INPUT:

        - ``use_smallest_mode`` -- if ``True``, use reps between -p/2 and p/2 instead
        - ``return_list`` -- if ``True``, return a list of coefficients (as integers);
          for use in printing

        EXAMPLES::

            sage: R.<a> = ZqFM(27,4)
            sage: b = a^2 + 5*a - 3
            sage: b._modp_rep()
            ((a^2 + 2*a), [0, 2, 1])
            sage: b._modp_rep(use_smallest_mode=True)[1]
            [0, -1, 1]
        """
        cdef FMElement rep = self._new_c()
        L = cmodp_rep(rep.value, self.value, smallest_mode if use_smallest_mode else simple_mode, return_list, self.prime_pow)
        if return_list:
            return rep, L
        else:
            return rep
