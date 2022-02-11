"""
Local and Global Genus Symbols
"""
#############################################################
#                                                           #
#   Wrappers for the Genus/Genus Symbol Code in ../genera/  #
#                                                           #
#############################################################

from sage.quadratic_forms.genera.genus import Genus, LocalGenusSymbol
from sage.rings.integer_ring import ZZ
from sage.arith.all import is_prime, prime_divisors


def global_genus_symbol(self):
    r"""
    Return the genus of two times a quadratic form over `\ZZ`.

    These are defined by a collection of local genus symbols (a la
    Chapter 15 of Conway-Sloane [CS1999]_), and a signature.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3,4])
        sage: Q.global_genus_symbol()
        Genus of
        [2 0 0 0]
        [0 4 0 0]
        [0 0 6 0]
        [0 0 0 8]
        Signature:  (4, 0)
        Genus symbol at 2:    [2^-2 4^1 8^1]_6
        Genus symbol at 3:     1^3 3^-1

    ::

        sage: Q = QuadraticForm(ZZ, 4, range(10))
        sage: Q.global_genus_symbol()
        Genus of
        [ 0  1  2  3]
        [ 1  8  5  6]
        [ 2  5 14  8]
        [ 3  6  8 18]
        Signature:  (3, 1)
        Genus symbol at 2:    1^-4
        Genus symbol at 563:     1^3 563^-1
    """
    if self.base_ring() is not ZZ:
        raise TypeError("the quadratic form is not defined over the integers")
    return Genus(self.Hessian_matrix())


def local_genus_symbol(self, p):
    r"""
    Return the Conway-Sloane genus symbol of 2 times a quadratic form
    defined over `\ZZ` at a prime number `p`.

    This is defined (in the
    Genus_Symbol_p_adic_ring() class in the quadratic_forms/genera
    subfolder) to be a list of tuples (one for each Jordan component
    p^m*A at p, where A is a unimodular symmetric matrix with
    coefficients the p-adic integers) of the following form:

        1. If p>2 then return triples of the form [`m`, `n`, `d`] where

            `m` = valuation of the component

            `n` = rank of A

            `d` = det(A) in {1,u} for normalized quadratic non-residue u.

        2. If p=2 then return quintuples of the form [`m`,`n`,`s`, `d`, `o`] where

            `m` = valuation of the component

            `n` = rank of A

            `d` = det(A) in {1,3,5,7}

            `s` = 0 (or 1) if A is even (or odd)

            `o` = oddity of A (= 0 if s = 0) in Z/8Z
              = the trace of the diagonalization of A

    .. NOTE::

        The Conway-Sloane convention for describing the prime 'p = -1'
        is not supported here, and neither is the convention for
        including the 'prime' Infinity.  See note on p370 of Conway-Sloane
        (3rd ed) [CS1999]_ for a discussion of this convention.

    INPUT:

    - `p` -- a prime number > 0

    OUTPUT:

    a Conway-Sloane genus symbol at `p`, which is an
    instance of the Genus_Symbol_p_adic_ring class.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3,4])
        sage: Q.local_genus_symbol(2)
        Genus symbol at 2:    [2^-2 4^1 8^1]_6
        sage: Q.local_genus_symbol(3)
        Genus symbol at 3:     1^3 3^-1
        sage: Q.local_genus_symbol(5)
        Genus symbol at 5:     1^4
    """
    if not is_prime(p):
        raise TypeError("the number " + str(p) + " is not prime")
    if self.base_ring() is not ZZ:
        raise TypeError("the quadratic form is not defined over the integers")
    return LocalGenusSymbol(self.Hessian_matrix(), p)


def CS_genus_symbol_list(self, force_recomputation=False):
    """
    Return the list of Conway-Sloane genus symbols in increasing order of primes dividing 2*det.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3,4])
        sage: Q.CS_genus_symbol_list()
        [Genus symbol at 2:    [2^-2 4^1 8^1]_6, Genus symbol at 3:     1^3 3^-1]
    """
    # Try to use the cached list
    if not force_recomputation:
        try:
            return self.__CS_genus_symbol_list
        except AttributeError:
            pass

    # Otherwise recompute and cache the list
    list_of_CS_genus_symbols = []

    for p in prime_divisors(2 * self.det()):
        list_of_CS_genus_symbols.append(self.local_genus_symbol(p))

    self.__CS_genus_symbol_list = list_of_CS_genus_symbols
    return list_of_CS_genus_symbols
