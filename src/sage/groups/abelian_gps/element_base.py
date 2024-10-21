# sage_setup: distribution = sagemath-modules
"""
Base class for abelian group elements

This is the base class for both
:mod:`~sage.groups.abelian_gps.abelian_group_element` and
:mod:`~sage.groups.abelian_gps.dual_abelian_group_element`.

As always, elements are immutable once constructed.
"""


###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner  <wdjoyner@gmail.com>
#  Copyright (C) 2012 Volker Braun  <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
###########################################################################

from sage.structure.element import MultiplicativeGroupElement
from sage.misc.cachefunc import cached_method
from sage.arith.misc import GCD
from sage.arith.functions import lcm as LCM
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.structure.richcmp import richcmp


class AbelianGroupElementBase(MultiplicativeGroupElement):
    """
    Base class for abelian group elements.

    The group element is defined by a tuple whose ``i``-th entry is an
    integer in the range from 0 (inclusively) to ``G.gen(i).order()``
    (exclusively) if the `i`-th generator is of finite order, and an
    arbitrary integer if the `i`-th generator is of infinite order.

    INPUT:

    - ``exponents`` -- ``1`` or a list/tuple/iterable of integers; the
      exponent vector (with respect to the parent generators) defining
      the group element

    - ``parent`` -- abelian group; the parent of the group element

    EXAMPLES::

        sage: F = AbelianGroup(3,[7,8,9])
        sage: Fd = F.dual_group(names='ABC')                                            # needs sage.rings.number_field
        sage: A,B,C = Fd.gens()                                                         # needs sage.rings.number_field
        sage: A*B^-1 in Fd                                                              # needs sage.rings.number_field
        True
    """

    def __init__(self, parent, exponents):
        """
        Create an element.

        EXAMPLES::

            sage: F = AbelianGroup(3,[7,8,9])
            sage: Fd = F.dual_group(names='ABC')                                        # needs sage.rings.number_field
            sage: A,B,C = Fd.gens()                                                     # needs sage.rings.number_field
            sage: A*B^-1 in Fd                                                          # needs sage.rings.number_field
            True

        Check that :issue:`35216` is fixed::

            sage: M = AbelianGroup([3])
            sage: M([5]) == M([2])
            True
            sage: M([3]).is_trivial()
            True
        """
        MultiplicativeGroupElement.__init__(self, parent)
        n = parent.ngens()
        if exponents == 1:
            self._exponents = tuple(ZZ.zero() for i in range(n))
        else:
            if len(exponents) != n:
                raise IndexError('argument length (= %s) must be %s'
                                 % (len(exponents), n))
            self._exponents = tuple(ZZ(e % o if o else e) for e, o in
                                    zip(exponents, parent.gens_orders()))

    def __hash__(self):
        r"""
        TESTS::

            sage: F = AbelianGroup(3,[7,8,9])
            sage: hash(F.an_element()) # random
            1024
        """
        return hash(self.parent()) ^ hash(self._exponents)

    def exponents(self):
        """
        The exponents of the generators defining the group element.

        OUTPUT:

        A tuple of integers for an abelian group element. The integer
        can be arbitrary if the corresponding generator has infinite
        order. If the generator is of finite order, the integer is in
        the range from 0 (inclusive) to the order (exclusive).

        EXAMPLES::

            sage: F.<a,b,c,f> = AbelianGroup([7,8,9,0])
            sage: (a^3*b^2*c).exponents()
            (3, 2, 1, 0)
            sage: F([3, 2, 1, 0])
            a^3*b^2*c
            sage: (c^42).exponents()
            (0, 0, 6, 0)
            sage: (f^42).exponents()
            (0, 0, 0, 42)
        """
        return self._exponents

    def _libgap_(self):
        r"""
        TESTS::

            sage: F.<a,b,c> = AbelianGroup([7,8,9])
            sage: libgap(a**2 * c) * libgap(b * c**2)                                   # needs sage.libs.gap
            f1^2*f2*f6
        """
        from sage.misc.misc_c import prod
        from sage.libs.gap.libgap import libgap
        G = libgap(self.parent())
        return prod(g**i for g,i in zip(G.GeneratorsOfGroup(), self._exponents))

    def list(self):
        """
        Return a copy of the exponent vector.

        Use :meth:`exponents` instead.

        OUTPUT:

        The underlying coordinates used to represent this element.  If
        this is a word in an abelian group on `n` generators, then
        this is a list of nonnegative integers of length `n`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: F = AbelianGroup(5,[2, 3, 5, 7, 8], names='abcde')
            sage: a,b,c,d,e = F.gens()
            sage: Ad = F.dual_group(names='ABCDE')
            sage: A,B,C,D,E = Ad.gens()
            sage: (A*B*C^2*D^20*E^65).exponents()
            (1, 1, 2, 6, 1)
            sage: X = A*B*C^2*D^2*E^-6
            sage: X.exponents()
            (1, 1, 2, 2, 2)
        """
        # to be deprecated (really, return a list??). Use exponents() instead.
        return list(self._exponents)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = AbelianGroup([2])
            sage: G.gen(0)._repr_()
            'f'
            sage: G.one()._repr_()
            '1'
        """
        s = ""
        G = self.parent()
        for v_i, x_i in zip(self.exponents(), G.variable_names()):
            if v_i == 0:
                continue
            if len(s) > 0:
                s += '*'
            if v_i == 1:
                s += str(x_i)
            else:
                s += str(x_i) + '^' + str(v_i)
        if s:
            return s
        else:
            return '1'

    def _richcmp_(self, other, op):
        """
        Compare ``self`` and ``other``.

        The comparison is based on the exponents.

        OUTPUT: boolean

        EXAMPLES::

            sage: G.<a,b> = AbelianGroup([2,3])
            sage: a > b
            True

            sage: Gd.<A,B> = G.dual_group()                                             # needs sage.rings.number_field
            sage: A > B                                                                 # needs sage.rings.number_field
            True
        """
        return richcmp(self._exponents, other._exponents, op)

    @cached_method
    def order(self):
        """
        Return the order of this element.

        OUTPUT: integer or ``infinity``

        EXAMPLES::

            sage: F = AbelianGroup(3,[7,8,9])
            sage: Fd = F.dual_group()                                                   # needs sage.rings.number_field
            sage: A,B,C = Fd.gens()                                                     # needs sage.rings.number_field
            sage: (B*C).order()                                                         # needs sage.rings.number_field
            72

            sage: F = AbelianGroup(3,[7,8,9]); F
            Multiplicative Abelian group isomorphic to C7 x C8 x C9
            sage: F.gens()[2].order()
            9
            sage: a,b,c = F.gens()
            sage: (b*c).order()
            72
            sage: G = AbelianGroup(3,[7,8,9])
            sage: type((G.0 * G.1).order())==Integer
            True
        """
        M = self.parent()
        order = M.gens_orders()
        L = self.exponents()
        N = LCM([order[i]/GCD(order[i],L[i]) for i in range(len(order)) if L[i] != 0])
        if N == 0:
            return infinity
        else:
            return ZZ(N)

    multiplicative_order = order

    def _div_(left, right):
        """
        Divide ``left`` and ``right``.

        TESTS::

            sage: G.<a,b> = AbelianGroup(2)
            sage: a/b
            a*b^-1
            sage: a._div_(b)
            a*b^-1
        """
        G = left.parent()
        assert G is right.parent()
        exponents = [x - y for x, y in
                     zip(left._exponents, right._exponents)]
        return G.element_class(G, exponents)

    def _mul_(left, right):
        """
        Multiply ``left`` and ``right``.

        TESTS::

            sage: G.<a,b> = AbelianGroup(2)
            sage: a*b
            a*b
            sage: a._mul_(b)
            a*b
        """
        G = left.parent()
        assert G is right.parent()
        exponents = [x + y for x, y in
                     zip(left._exponents, right._exponents)]
        return G.element_class(G, exponents)

    def __pow__(self, n):
        """
        Exponentiate ``self``.

        TESTS::

            sage: G.<a,b> = AbelianGroup(2)
            sage: a^3
            a^3
        """
        m = Integer(n)
        if n != m:
            raise TypeError('argument n (= '+str(n)+') must be an integer.')
        G = self.parent()
        exponents = [m * e for e in self._exponents]
        return G.element_class(G, exponents)

    def __invert__(self):
        """
        Return the inverse element.

        EXAMPLES::

            sage: G.<a,b> = AbelianGroup([0,5])
            sage: a.inverse()  # indirect doctest
            a^-1
            sage: a.__invert__()
            a^-1
            sage: a^-1
            a^-1
            sage: ~a
            a^-1
            sage: (a*b).exponents()
            (1, 1)
            sage: (a*b).inverse().exponents()
            (-1, 4)
        """
        G = self.parent()
        exponents = [-e for e in self._exponents]
        return G.element_class(G, exponents)

    def is_trivial(self):
        """
        Test whether ``self`` is the trivial group element ``1``.

        OUTPUT: boolean

        EXAMPLES::

            sage: G.<a,b> = AbelianGroup([0,5])
            sage: (a^5).is_trivial()
            False
            sage: (b^5).is_trivial()
            True
        """
        return all(e == 0 for e in self._exponents)
