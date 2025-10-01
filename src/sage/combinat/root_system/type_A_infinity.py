"""
Root system data for type A infinity
"""
# ***************************************************************************
# Copyright (C) 2016 Andrew Mathas <Andrew dot Mathas at Sydney dot edu dot au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from .cartan_type import CartanType_standard, CartanType_simple
from sage.rings.integer_ring import ZZ


class CartanType(CartanType_standard, CartanType_simple):
    r"""
    The Cartan type `A_{\infty}`.

    We use ``NN`` and ``ZZ`` to explicitly differentiate between the
    `A_{+\infty}` and `A_{\infty}` root systems, respectively.
    While ``oo`` is the same as ``+Infinity`` in Sage, it is used as
    an alias for ``ZZ``.
    """
    # We do not inherit from CartanType_crystallographic because it provides
    # methods that are not implemented for A_oo.

    def __init__(self, index_set):
        """
        EXAMPLES::

            sage: CartanType(['A',oo]) is CartanType(['A', ZZ])
            True
            sage: CartanType(['A',oo]) is CartanType(['A', NN])
            False
            sage: ct=CartanType(['A',ZZ])
            sage: ct
            ['A', ZZ]
            sage: ct._repr_(compact = True)
            'A_ZZ'
            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            False
            sage: ct.is_affine()
            False
            sage: ct.is_untwisted_affine()
            False
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.dual()
            ['A', ZZ]

        TESTS::

            sage: TestSuite(ct).run()
        """
        super().__init__()
        self.letter = 'A'
        self.n = index_set

    def _repr_(self, compact=False):
        """
        Return a representation of ``self``.

        TESTS::

            sage: CartanType(['A',ZZ])
            ['A', ZZ]
            sage: CartanType(['A',NN])._repr_(compact=True)
            'A_NN'
        """
        ret = '%s_%s' if compact else "['%s', %s]"
        return ret % (self.letter, 'ZZ' if self.n == ZZ else 'NN')

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex( CartanType(['A',NN]) )
            A_{\Bold{N}}
            sage: latex( CartanType(['A',ZZ]) )
            A_{\Bold{Z}}
        """
        return 'A_{{{}}}'.format(self.n._latex_())

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A', ZZ]).ascii_art())
            ..---O---O---O---O---O---O---O---..
                -3  -2  -1   0   1   2   3
            sage: print(CartanType(['A', NN]).ascii_art())
            O---O---O---O---O---O---O---..
            0   1   2   3   4   5   6
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node

        if self.n == ZZ:
            ret = '..---' + '---'.join(node(label(i)) for i in range(7)) + '---..\n'
            ret += '  ' + ''.join("{:4}".format(label(i)) for i in range(-3, 4))
        else:
            ret = '---'.join(node(label(i)) for i in range(7)) + '---..\n'
            ret += '0' + ''.join("{:4}".format(label(i)) for i in range(1, 7))

        return ret

    def dual(self):
        """
        Simply laced Cartan types are self-dual, so return ``self``.

        EXAMPLES::

            sage: CartanType(["A", NN]).dual()
            ['A', NN]
            sage: CartanType(["A", ZZ]).dual()
            ['A', ZZ]
        """
        return self

    def is_simply_laced(self):
        """
        Return ``True`` because ``self`` is simply laced.

        EXAMPLES::

            sage: CartanType(['A', NN]).is_simply_laced()
            True
            sage: CartanType(['A', ZZ]).is_simply_laced()
            True
        """
        return True

    def is_crystallographic(self):
        """
        Return ``False`` because ``self`` is not crystallographic.

        EXAMPLES::

            sage: CartanType(['A', NN]).is_crystallographic()
            True
            sage: CartanType(['A', ZZ]).is_crystallographic()
            True
        """
        return True

    def is_finite(self):
        """
        Return ``True`` because ``self`` is not finite.

        EXAMPLES::

            sage: CartanType(['A', NN]).is_finite()
            False
            sage: CartanType(['A', ZZ]).is_finite()
            False
        """
        return False

    def is_affine(self):
        """
        Return ``False`` because ``self`` is not (untwisted) affine.

        EXAMPLES::

            sage: CartanType(['A', NN]).is_affine()
            False
            sage: CartanType(['A', ZZ]).is_affine()
            False
        """
        return False

    def is_untwisted_affine(self):
        """
        Return ``False`` because ``self`` is not (untwisted) affine.

        EXAMPLES::

            sage: CartanType(['A', NN]).is_untwisted_affine()
            False
            sage: CartanType(['A', ZZ]).is_untwisted_affine()
            False
        """
        return False

    def rank(self):
        """
        Return the rank of ``self`` which for type `X_n` is `n`.

        EXAMPLES::

            sage: CartanType(['A', NN]).rank()
            +Infinity
            sage: CartanType(['A', ZZ]).rank()
            +Infinity

        As this example shows, the rank is slightly ambiguous because the root
        systems of type `['A',NN]` and type `['A',ZZ]` have the same rank.
        Instead, it is better ot use :meth:`index_set` to differentiate between
        these two root systems.
        """
        return self.n.cardinality()

    def type(self):
        """
        Return the type of ``self``.

        EXAMPLES::

            sage: CartanType(['A', NN]).type()
            'A'
            sage: CartanType(['A', ZZ]).type()
            'A'
        """
        return self.letter

    def index_set(self):
        r"""
        Return the index set for the Cartan type ``self``.

        The index set for all standard finite Cartan types is of the form
        `\{1, \ldots, n\}`. (See :mod:`~sage.combinat.root_system.type_I`
        for a slight abuse of this).

        EXAMPLES::

            sage: CartanType(['A', NN]).index_set()
            Non negative integer semiring
            sage: CartanType(['A', ZZ]).index_set()
            Integer Ring
        """
        return self.n
