# sage_setup: distribution = sagemath-categories
"""
The set of prime numbers

AUTHORS:

 - William Stein (2005): original version
 - Florent Hivert (2009-11): adapted to the category framework.
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.integer_ring import ZZ
from .set import Set_generic
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.arith.misc import nth_prime
from sage.structure.unique_representation import UniqueRepresentation


class Primes(Set_generic, UniqueRepresentation):
    """
    The set of prime numbers.

    EXAMPLES::

        sage: P = Primes(); P
        Set of all prime numbers: 2, 3, 5, 7, ...

    We show various operations on the set of prime numbers::

        sage: P.cardinality()
        +Infinity
        sage: R = Primes()
        sage: P == R
        True
        sage: 5 in P
        True
        sage: 100 in P
        False

        sage: len(P)
        Traceback (most recent call last):
        ...
        NotImplementedError: infinite set
    """
    @staticmethod
    def __classcall__(cls, proof=True):
        """
        TESTS::

            sage: Primes(proof=True) is Primes()
            True
            sage: Primes(proof=False) is Primes()
            False
        """
        return super().__classcall__(cls, proof)

    def __init__(self, proof):
        """
        EXAMPLES::

            sage: P = Primes(); P
            Set of all prime numbers: 2, 3, 5, 7, ...

            sage: Q = Primes(proof=False); Q
            Set of all prime numbers: 2, 3, 5, 7, ...

        TESTS::

            sage: P.category()
            Category of facade infinite enumerated sets
            sage: TestSuite(P).run()                                                    # needs sage.libs.pari

            sage: Q.category()
            Category of facade infinite enumerated sets
            sage: TestSuite(Q).run()                                                    # needs sage.libs.pari

        The set of primes can be compared to various things,
        but is only equal to itself::

            sage: P = Primes()
            sage: R = Primes()
            sage: P == R
            True
            sage: P != R
            False
            sage: Q = [1,2,3]
            sage: Q != P        # indirect doctest
            True
            sage: R.<x> = ZZ[]
            sage: P != x^2+x
            True
        """
        super().__init__(facade=ZZ,
                         category=InfiniteEnumeratedSets())
        self.__proof = proof

    def _repr_(self):
        """
        Representation of the set of primes.

        EXAMPLES::

            sage: P = Primes(); P
            Set of all prime numbers: 2, 3, 5, 7, ...
        """
        return "Set of all prime numbers: 2, 3, 5, 7, ..."

    def __contains__(self, x):
        """
        Check whether an object is a prime number.

        EXAMPLES::

            sage: P = Primes()
            sage: 5 in P
            True
            sage: 100 in P
            False
            sage: 1.5 in P
            False
            sage: e in P                                                                # needs sage.symbolic
            False
        """
        try:
            if x not in ZZ:
                return False
            return ZZ(x).is_prime()
        except TypeError:
            return False

    def _an_element_(self):
        """
        Return a typical prime number.

        EXAMPLES::

            sage: P = Primes()
            sage: P._an_element_()
            43
        """
        return ZZ(43)

    def first(self):
        """
        Return the first prime number.

        EXAMPLES::

            sage: P = Primes()
            sage: P.first()
            2
        """
        return ZZ(2)

    def next(self, pr):
        """
        Return the next prime number.

        EXAMPLES::

            sage: P = Primes()
            sage: P.next(5)                                                             # needs sage.libs.pari
            7
        """
        pr = pr.next_prime(self.__proof)
        return pr

    def unrank(self, n):
        """
        Return the n-th prime number.

        EXAMPLES::

            sage: P = Primes()
            sage: P.unrank(0)                                                           # needs sage.libs.pari
            2
            sage: P.unrank(5)                                                           # needs sage.libs.pari
            13
            sage: P.unrank(42)                                                          # needs sage.libs.pari
            191
        """
        return nth_prime(n + 1)
