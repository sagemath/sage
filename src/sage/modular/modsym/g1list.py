r"""
List of coset representatives for `\Gamma_1(N)` in `\SL_2(\ZZ)`
"""

# ****************************************************************************
#       Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import GCD
from sage.structure.richcmp import richcmp_method, richcmp
from sage.structure.sage_object import SageObject
from sage.misc.persist import register_unpickle_override


@richcmp_method
class G1list(SageObject):
    r"""
    A class representing a list of coset representatives for `\Gamma_1(N)` in
    `\SL_2(\ZZ)`. What we actually calculate is a list of elements of
    `(\ZZ/N\ZZ)^2` of exact order `N`.

    TESTS::

        sage: L = sage.modular.modsym.g1list.G1list(18)
        sage: loads(dumps(L)) == L
        True
    """
    def __init__(self, N):
        """
        EXAMPLES::

            sage: L = sage.modular.modsym.g1list.G1list(6); L # indirect doctest
            List of coset representatives for Gamma_1(6) in SL_2(Z)
        """
        self.__N = N
        self.__list = [(u, v) for u in range(N) for v in range(N)
                       if GCD(GCD(u, v), N) == 1]

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: L1 = sage.modular.modsym.g1list.G1list(6)
            sage: L2 = sage.modular.modsym.g1list.G1list(7)
            sage: L1 < L2
            True
            sage: L1 == QQ
            False
        """
        if not isinstance(other, G1list):
            return NotImplemented
        else:
            return richcmp(self.__N, other.__N, op)

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: L = sage.modular.modsym.g1list.G1list(19); L[100] # indirect doctest
            (5, 6)
        """
        return self.__list[i]

    def __len__(self):
        """
        Return the length of the underlying list.

        EXAMPLES::

            sage: L = sage.modular.modsym.g1list.G1list(24); len(L) # indirect doctest
            384
        """
        return len(self.__list)

    def __repr__(self):
        """
        String representation of ``self``.

        EXAMPLES::

            sage: L = sage.modular.modsym.g1list.G1list(3); L.__repr__()
            'List of coset representatives for Gamma_1(3) in SL_2(Z)'
        """
        return f"List of coset representatives for Gamma_1({self.__N}) in SL_2(Z)"

    def list(self):
        r"""
        Return a list of vectors representing the cosets.

        Do not change the returned list!

        EXAMPLES::

            sage: L = sage.modular.modsym.g1list.G1list(4); L.list()
            [(0, 1), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (2, 1), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3)]
        """
        return self.__list

    def normalize(self, u, v):
        r"""
        Given a pair `(u,v)` of integers, return the unique pair `(u', v')`
        such that the pair `(u', v')` appears in ``self.list()`` and `(u, v)`
        is equivalent to `(u', v')`. This is rather trivial, but is here for
        consistency with the ``P1List`` class which is the equivalent for
        `\Gamma_0` (where the problem is rather harder).

        This will only make sense if `{\rm gcd}(u, v, N) = 1`; otherwise the
        output will not be an element of ``self``.

        EXAMPLES::

            sage: L = sage.modular.modsym.g1list.G1list(4); L.normalize(6, 1)
            (2, 1)
            sage: L = sage.modular.modsym.g1list.G1list(4); L.normalize(6, 2) # nonsense!
            (2, 2)
        """
        return u % self.__N, v % self.__N


class _G1list_old_pickle(G1list):
    """
    This class exists only for dealing with old pickles.

    This needs to handle both old-style class pickles, where there is
    no input to the class on the initial ``__init__`` call, and the
    new class pickles, we need to have ``__setstate__`` handle it.
    """
    def __init__(self):
        """
        For unpickling old pickles.

        TESTS::

            sage: from sage.modular.modsym.g1list import _G1list_old_pickle
            sage: L = _G1list_old_pickle()
            sage: type(L) == G1list
            True
        """
        self.__class__ = G1list

    def __setstate__(self, state):
        """
        For unpickling new pickles.

        TESTS::

            sage: from sage.modular.modsym.g1list import G1list
            sage: L = G1list(6)
            sage: Lp = loads(dumps(L))
            sage: L == Lp
            True
            sage: type(Lp) == G1list
            True
        """
        # We don't really want this class, but we want to handle new
        #   pickles without creating a new class
        self.__class__ = G1list
        self.__dict__ = state  # Default pickling is ``state = self.__dict__``


register_unpickle_override('sage.modular.modsym.g1list', 'G1list',
                           _G1list_old_pickle)
