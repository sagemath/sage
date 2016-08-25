"Cremona modular symbols"
from __future__ import print_function

include "cysignals/signals.pxi"

from ..eclib cimport mat
from .mat cimport MatrixFactory

cdef MatrixFactory MF = MatrixFactory()

cdef class ModularSymbols:
    """
    Class of Cremona Modular Symbols of given level and sign (and weight 2).

    EXAMPLES::

        sage: M = CremonaModularSymbols(225)
        sage: type(M)
        <type 'sage.libs.eclib.homspace.ModularSymbols'>
    """
    def __init__(self, long level, int sign=0, bint cuspidal=False, int verbose=0):
        """
        Called when creating a space of Cremona modular symbols.

        INPUT:

        - ``level`` (int) -- the level: an integer, at least 2.
        - ``sign`` (int, default 0) -- the sign: 0, +1 or -1
        - ``cuspidal`` (boolean, default False) -- True for cuspidal homology
        - ``verbose`` (int, default 0) -- verbosity level

        EXAMPLES::

            sage: CremonaModularSymbols(123, sign=1, cuspidal=True)
            Cremona Cuspidal Modular Symbols space of dimension 13 for Gamma_0(123) of weight 2 with sign 1
            sage: CremonaModularSymbols(123, sign=-1, cuspidal=True)
            Cremona Cuspidal Modular Symbols space of dimension 13 for Gamma_0(123) of weight 2 with sign -1
            sage: CremonaModularSymbols(123, sign=0, cuspidal=True)
            Cremona Cuspidal Modular Symbols space of dimension 26 for Gamma_0(123) of weight 2 with sign 0
            sage: CremonaModularSymbols(123, sign=0, cuspidal=False)
            Cremona Modular Symbols space of dimension 29 for Gamma_0(123) of weight 2 with sign 0
        """
        if not (sign == 0 or sign == 1 or sign == -1):
            raise ValueError("sign (= %s) is not supported; use 0, +1 or -1" % sign)
        if level <= 1:
            raise ValueError("the level (= %s) must be at least 2" % level)
        sig_on()
        self.H = new homspace(level, sign, cuspidal, verbose)
        sig_off()

    def __dealloc__(self):
        del self.H

    def __repr__(self):
        """
        String representation of space of Cremona modular symbols.

        EXAMPLES:

        We test various types of spaces that impact printing::

            sage: M = CremonaModularSymbols(37, sign=1)
            sage: M.__repr__()
            'Cremona Modular Symbols space of dimension 3 for Gamma_0(37) of weight 2 with sign 1'

            sage: CremonaModularSymbols(37, cuspidal=True).__repr__()
            'Cremona Cuspidal Modular Symbols space of dimension 4 for Gamma_0(37) of weight 2 with sign 0'
        """
        return "Cremona %sModular Symbols space of dimension %s for Gamma_0(%s) of weight 2 with sign %s"%(
            'Cuspidal ' if self.is_cuspidal() else '',
            self.dimension(), self.level(), self.sign())

    #cpdef long level(self):
    def level(self):
        """
        Return the level of this modular symbols space.

        EXAMPLES::

            sage: M = CremonaModularSymbols(1234, sign=1)
            sage: M.level()
            1234
        """
        return self.H.modulus

    #cpdef int dimension(self):
    def dimension(self):
        """
        Return the dimension of this modular symbols space.

        EXAMPLES::

            sage: M = CremonaModularSymbols(1234, sign=1)
            sage: M.dimension()
            156
        """
        if self.is_cuspidal():
           return self.H.h1cuspdim()
        else:
           return self.H.h1dim()

    def number_of_cusps(self):
        r"""
        Return the number of cusps for $\Gamma_0(N)$, where $N$ is the
        level.

        EXAMPLES::

            sage: M = CremonaModularSymbols(225)
            sage: M.number_of_cusps()
            24
        """
        return self.H.h1ncusps()

    #cpdef int sign(self):
    def sign(self):
        """
        Return the sign of this Cremona modular symbols space.  The sign is either 0, +1 or -1.

        EXAMPLES::

            sage: M = CremonaModularSymbols(1122, sign=1); M
            Cremona Modular Symbols space of dimension 224 for Gamma_0(1122) of weight 2 with sign 1
            sage: M.sign()
            1
            sage: M = CremonaModularSymbols(1122); M
            Cremona Modular Symbols space of dimension 433 for Gamma_0(1122) of weight 2 with sign 0
            sage: M.sign()
            0
            sage: M = CremonaModularSymbols(1122, sign=-1); M
            Cremona Modular Symbols space of dimension 209 for Gamma_0(1122) of weight 2 with sign -1
            sage: M.sign()
            -1
        """
        return self.H.plusflag

    #cpdef bint is_cuspidal(self):
    def is_cuspidal(self):
        """
        Return whether or not this space is cuspidal.

        EXAMPLES::

            sage: M = CremonaModularSymbols(1122); M.is_cuspidal()
            0
            sage: M = CremonaModularSymbols(1122, cuspidal=True); M.is_cuspidal()
            1
        """
        return self.H.cuspidal

    def hecke_matrix(self, long p, dual=False, verbose=False):
        """
        Return the matrix of the ``p``-th Hecke operator acting on
        this space of modular symbols.

        The result of this command is not cached.

        INPUT:

        - ``p`` -- a prime number

        - ``dual`` -- (default: False) whether to compute the Hecke
                    operator acting on the dual space, i.e., the
                    transpose of the Hecke operator

        - ``verbose`` -- (default: False) print verbose output

        OUTPUT:

        (matrix) If ``p`` divides the level, the matrix of the
        Atkin-Lehner involution `W_p` at ``p``; otherwise the matrix of the
        Hecke operator `T_p`,

        EXAMPLES::

            sage: M = CremonaModularSymbols(37)
            sage: t = M.hecke_matrix(2); t
            5 x 5 Cremona matrix over Rational Field
            sage: print(t.str())
            [ 3  0  0  0  0]
            [-1 -1  1  1  0]
            [ 0  0 -1  0  1]
            [-1  1  0 -1 -1]
            [ 0  0  1  0 -1]
            sage: t.charpoly().factor()
            (x - 3) * x^2 * (x + 2)^2
            sage: print(M.hecke_matrix(2, dual=True).str())
            [ 3 -1  0 -1  0]
            [ 0 -1  0  1  0]
            [ 0  1 -1  0  1]
            [ 0  1  0 -1  0]
            [ 0  0  1 -1 -1]
            sage: w = M.hecke_matrix(37); w
            5 x 5 Cremona matrix over Rational Field
            sage: w.charpoly().factor()
            (x - 1)^2 * (x + 1)^3
            sage: sw = w.sage_matrix_over_ZZ()
            sage: st = t.sage_matrix_over_ZZ()
            sage: sw^2 == sw.parent()(1)
            True
            sage: st*sw == sw*st
            True
        """
        sig_on()
        cdef mat M = self.H.heckeop(p, dual, verbose)
        sig_off()
        return MF.new_matrix(M)
