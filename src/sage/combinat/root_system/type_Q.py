# sage_setup: distribution = sagemath-modules
"""
Root system data for type Q
"""
#*****************************************************************************
#       Copyright (C) 2018 Wencin Poh
#       Copyright (C) 2018 Anne Schilling
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from .cartan_type import CartanType_standard_finite
from sage.combinat.root_system.root_system import RootSystem


class CartanType(CartanType_standard_finite):
    """
    Cartan Type `Q_n`.

    .. SEEALSO:: :func:`~sage.combinat.root_systems.cartan_type.CartanType`
    """

    def __init__(self, m):
        """
        EXAMPLES::

            sage: ct = CartanType(['Q',4])
            sage: ct
            ['Q', 4]
            sage: ct._repr_(compact = True)
            'Q4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_simply_laced()
            True
            sage: ct.dual()
            ['Q', 4]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert m >= 2
        CartanType_standard_finite.__init__(self, "Q", m-1)

    def _repr_(self, compact=False):
        """
        TESTS::

            sage: ct = CartanType(['Q',4])
            sage: repr(ct)
            "['Q', 4]"
            sage: ct._repr_(compact=True)
            'Q4'
        """
        format = '%s%s' if compact else "['%s', %s]"
        return format % (self.letter, self.n+1)

    def __reduce__(self):
        """
        TESTS::

            sage: T = CartanType(['Q', 4])
            sage: T.__reduce__()
            (CartanType, ('Q', 4))
            sage: T == loads(dumps(T))
            True
        """
        from .cartan_type import CartanType
        return (CartanType, (self.letter, self.n+1))

    def index_set(self):
        r"""
        Return the index set for Cartan type Q.

        The index set for type Q is of the form
        `\{-n, \ldots, -1, 1, \ldots, n\}`.

        EXAMPLES::

            sage: CartanType(['Q', 3]).index_set()
            (1, 2, -2, -1)
        """
        return tuple(list(range(1, self.n + 1)) + list(range(-self.n, 0)))

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['Q',4]))
            Q_{4}
        """
        return "Q_{%s}" % (self.n + 1)

    def root_system(self):
        """
        Return the root system of ``self``.

        EXAMPLES::

            sage: Q = CartanType(['Q',3])
            sage: Q.root_system()
            Root system of type ['A', 2]
        """
        return RootSystem(['A',self.n])

    def is_irreducible(self):
        """
        Return whether this Cartan type is irreducible.

        EXAMPLES::

            sage: Q = CartanType(['Q',3])
            sage: Q.is_irreducible()
            True
        """
        return True

    def is_simply_laced(self):
        """
        Return whether this Cartan type is simply-laced.

        EXAMPLES::

            sage: Q = CartanType(['Q',3])
            sage: Q.is_simply_laced()
            True
        """
        return True

    def dual(self):
        """
        Return dual of ``self``.

        EXAMPLES::

            sage: Q = CartanType(['Q',3])
            sage: Q.dual()
            ['Q', 3]
        """
        return self
