"""
Formal sums

AUTHORS:

 - David Harvey (2006-09-20): changed FormalSum not to derive from
   "list" anymore, because that breaks new Element interface

 - Nick Alexander (2006-12-06): added test cases.

 - William Stein (2006, 2009): wrote the first version in 2006, documented it in 2009.

 - Volker Braun (2010-07-19): new-style coercions, documentation
   added. FormalSums now derives from UniqueRepresentation.

FUNCTIONS:
    - ``FormalSums(ring)`` -- create the module of formal finite sums with
                          coefficients in the given ring.

    - ``FormalSum(list of pairs (coeff, number))`` -- create a formal sum

EXAMPLES::

    sage: A = FormalSum([(1, 2/3)]); A
    2/3
    sage: B = FormalSum([(3, 1/5)]); B
    3*1/5
    sage: -B
    -3*1/5
    sage: A + B
    2/3 + 3*1/5
    sage: A - B
    2/3 - 3*1/5
    sage: B*3
    9*1/5
    sage: 2*A
    2*2/3
    sage: list(2*A + A)
    [(3, 2/3)]

TESTS::

    sage: R = FormalSums(QQ)
    sage: loads(dumps(R)) == R
    True
    sage: a = R(2/3) + R(-5/7); a
    2/3 + -5/7
    sage: loads(dumps(a)) == a
    True
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.repr import repr_lincomb
import operator
from collections import OrderedDict

from sage.modules.module import Module
from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp
from sage.rings.integer_ring import ZZ
from sage.structure.parent import Parent
from sage.structure.coerce import LeftModuleAction, RightModuleAction
from sage.categories.action import PrecomposedAction
from sage.structure.unique_representation import UniqueRepresentation


class FormalSum(ModuleElement):
    """
    A formal sum over a ring.
    """
    def __init__(self, x, parent=None, check=True, reduce=True):
        """
        INPUT:
            - ``x`` -- object
            - ``parent`` -- FormalSums(R) module (default: FormalSums(ZZ))
            - ``check`` -- bool (default: True) if False, might not coerce
                           coefficients into base ring, which can speed
                           up constructing a formal sum.
            - ``reduce`` -- reduce (default: True) if False, do not
                            combine common terms

        EXAMPLES::

            sage: FormalSum([(1,2/3), (3,2/3), (-5, 7)])
            4*2/3 - 5*7
            sage: a = FormalSum([(1,2/3), (3,2/3), (-5, 7)], reduce=False); a
            2/3 + 3*2/3 - 5*7
            sage: a.reduce()
            sage: a
            4*2/3 - 5*7
            sage: FormalSum([(1,2/3), (3,2/3), (-5, 7)], parent=FormalSums(GF(5)))
            4*2/3

        Notice below that the coefficient 5 doesn't get reduced modulo 5::

            sage: FormalSum([(1,2/3), (3,2/3), (-5, 7)], parent=FormalSums(GF(5)), check=False)
            4*2/3 - 5*7

        Make sure we first reduce before checking coefficient types::

            sage: x,y = var('x, y')
            sage: FormalSum([(1/2,x), (2,y)], FormalSums(QQ))
            1/2*x + 2*y
            sage: FormalSum([(1/2,x), (2,y)], FormalSums(ZZ))
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: FormalSum([(1/2,x), (1/2,x), (2,y)], FormalSums(ZZ))
            x + 2*y
        """
        if x == 0:
            x = []
        self._data = x
        if parent is None:
            parent = formal_sums
        ModuleElement.__init__(self, parent)
        assert isinstance(parent, parent.category().parent_class)
        if reduce:  # first reduce
            self.reduce()
        if check:   # then check
            k = parent.base_ring()
            try:
                self._data = [(k(t[0]), t[1]) for t in self._data]
            except (IndexError, KeyError) as msg:
                raise TypeError("%s\nInvalid formal sum"%msg)

    def __iter__(self):
        """
        EXAMPLES::

            sage: for z in FormalSum([(1, 2), (-3, 7), (5, 1000)]):
            ....:     print(z)
            (1, 2)
            (-3, 7)
            (5, 1000)
        """
        return iter(self._data)

    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: v = FormalSum([(1, 2), (-3, 7), (5, 1000)]); v
            2 - 3*7 + 5*1000
            sage: v[0]
            (1, 2)
            sage: v[1]
            (-3, 7)
            sage: v[2]
            (5, 1000)
            sage: list(v)
            [(1, 2), (-3, 7), (5, 1000)]
        """
        return self._data[n]

    def __len__(self):
        """
        EXAMPLES::

            sage: v = FormalSum([(1, 2), (5, 1000), (-3, 7)])
            sage: len(v)
            3
        """
        return len(self._data)

    def _repr_(self):
        """
        EXAMPLES::

            sage: a = FormalSum([(1,2/3), (-3,4/5), (7,Mod(2,3))])
            sage: a   # random, since comparing Mod(2,3) and rationals ill-defined
            sage: a._repr_()    # random
            '2/3 - 3*4/5 + 7*2'
        """
        return repr_lincomb([t, c] for c, t in self)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(FormalSum([(1,2), (5, 8/9), (-3, 7)]))
            2 + 5\cdot \frac{8}{9} - 3\cdot 7
        """
        from sage.misc.latex import repr_lincomb
        symbols = [z[1] for z in self]
        coeffs = [z[0] for z in self]
        return repr_lincomb(symbols, coeffs)
        # TODO: finish merging sage.misc.latex.repr_lincomb and
        # sage.misc.misc.repr_lincomb and use instead:
        # return repr_lincomb([[t,c] for c,t in self], is_latex=True)

    def _richcmp_(self, other, op):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- a :class:`FormalSum` with the same parent

        - ``op`` -- a comparison operator

        EXAMPLES::

            sage: a = FormalSum([(1,3),(2,5)]); a
            3 + 2*5
            sage: b = FormalSum([(1,3),(2,7)]); b
            3 + 2*7
            sage: a != b
            True
            sage: a_QQ = FormalSum([(1,3),(2,5)],parent=FormalSums(QQ))
            sage: a == a_QQ       # a is coerced into FormalSums(QQ)
            True
            sage: a == 0          # 0 is coerced into a.parent()(0)
            False
        """
        return richcmp(self._data, other._data, op)

    def _neg_(self):
        """
        EXAMPLES::

            sage: -FormalSum([(1,3),(2,5)])
            -3 - 2*5
        """
        return self.__class__([(-c, s) for (c, s) in self._data], check=False, parent=self.parent())

    def _add_(self, other):
        """
        EXAMPLES::

            sage: FormalSum([(1,3/7),(2,5/8)]) + FormalSum([(1,3/7),(-2,5)])  # indirect doctest
            2*3/7 + 2*5/8 - 2*5
        """
        return self.__class__(self._data + other._data, check=False, parent=self.parent())

    def _lmul_(self, s):
        """
        EXAMPLES::

            sage: FormalSum([(1,3/7),(-2,5)])*(-3)
            -3*3/7 + 6*5
        """
        return self.__class__([(c*s, x) for (c, x) in self], check=False, parent=self.parent())

    def _rmul_(self, s):
        """
        EXAMPLES::

            sage: -3*FormalSum([(1,3/7),(-2,5)])
            -3*3/7 + 6*5
        """
        return self.__class__([(s*c, x) for (c, x) in self], check=False, parent=self.parent())

    def __bool__(self):
        """
        EXAMPLES::

            sage: bool(FormalSum([(1,3/7),(-2,5)]))
            True
            sage: bool(FormalSums(QQ)(0))
            False
            sage: bool(FormalSums(QQ)(1))
            True
        """
        if not len(self._data):
            return False
        for c, _ in self._data:
            if not c.is_zero():
                return True
        return False

    __nonzero__ = __bool__

    def reduce(self):
        """
        EXAMPLES::

            sage: a = FormalSum([(-2,3), (2,3)], reduce=False); a
            -2*3 + 2*3
            sage: a.reduce()
            sage: a
            0
        """
        new = OrderedDict()
        for coeff, x in self:
            try:
                coeff += new[x]
            except KeyError:
                pass
            new[x] = coeff
        self._data = [(c, x) for (x, c) in new.items() if c]


class FormalSums(UniqueRepresentation, Module):
    """
    The R-module of finite formal sums with coefficients in some ring R.

    EXAMPLES::

        sage: FormalSums()
        Abelian Group of all Formal Finite Sums over Integer Ring
        sage: FormalSums(ZZ)
        Abelian Group of all Formal Finite Sums over Integer Ring
        sage: FormalSums(GF(7))
        Abelian Group of all Formal Finite Sums over Finite Field of size 7
        sage: FormalSums(ZZ[sqrt(2)])
        Abelian Group of all Formal Finite Sums over Order in Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?
        sage: FormalSums(GF(9,'a'))
        Abelian Group of all Formal Finite Sums over Finite Field in a of size 3^2

    TESTS::

        sage: TestSuite(FormalSums(QQ)).run()

    """
    Element = FormalSum
    @staticmethod
    def __classcall__(cls, base_ring = ZZ):
        """
        Set the default value for the base ring.

        EXAMPLES::

            sage: FormalSums(ZZ) == FormalSums()    # indirect test
            True
        """
        return UniqueRepresentation.__classcall__(cls, base_ring)

    def _repr_(self):
        """
        EXAMPLES::

            sage: FormalSums(GF(7))
            Abelian Group of all Formal Finite Sums over Finite Field of size 7
            sage: FormalSums(GF(7))._repr_()
            'Abelian Group of all Formal Finite Sums over Finite Field of size 7'
        """
        return "Abelian Group of all Formal Finite Sums over %s"%self.base_ring()

    def _element_constructor_(self, x, check=True, reduce=True):
        """
        Make a formal sum in self from x.

        INPUT:

        - ``x`` -- formal sum, list or number

        - ``check`` -- bool (default: True)

        - ``reduce`` -- bool (default: True); whether to combine terms

        EXAMPLES::

            sage: P = FormalSum([(1,2/3)]).parent()
            sage: P([(1,2/3), (5,-2/9)])  # indirect test
            2/3 + 5*-2/9
        """
        if isinstance(x, FormalSum):
            P = x.parent()
            if P is self:
                return x
            else:
                x = x._data
        if isinstance(x, list):
            return self.element_class(x, check=check,reduce=reduce,parent=self)
        if x == 0:
            return self.element_class([], check=False, reduce=False, parent=self)
        else:
            return self.element_class([(self.base_ring()(1), x)], check=False, reduce=False, parent=self)

    def _coerce_map_from_(self, X):
        r"""
        Return whether there is a coercion from ``X``

        EXAMPLES::

            sage: FormalSums(QQ).has_coerce_map_from( FormalSums(ZZ) )   # indirect test
            True

            sage: FormalSums(ZZ).get_action(QQ)   # indirect test
            Right scalar multiplication by Rational Field on Abelian Group of all Formal Finite Sums over Rational Field
            with precomposition on left by Coercion map:
              From: Abelian Group of all Formal Finite Sums over Integer Ring
              To:   Abelian Group of all Formal Finite Sums over Rational Field
        """
        if isinstance(X,FormalSums):
            if self.base_ring().has_coerce_map_from(X.base_ring()):
                return True
        return False

    def base_extend(self, R):
        """
        EXAMPLES::

            sage: F7 = FormalSums(ZZ).base_extend(GF(7)); F7
            Abelian Group of all Formal Finite Sums over Finite Field of size 7

        The following tests against a bug that was fixed at :trac:`18795`::

            sage: isinstance(F7, F7.category().parent_class)
            True
        """
        if self.base_ring().has_coerce_map_from(R):
            return self
        elif R.has_coerce_map_from(self.base_ring()):
            return FormalSums(R)

    def _get_action_(self, other, op, self_is_left):
        """
        EXAMPLES::

            sage: A = FormalSums(RR);  A.get_action(RR)     # indirect doctest
            Right scalar multiplication by Real Field with 53 bits of precision on Abelian Group of all Formal Finite Sums over Real Field with 53 bits of precision

            sage: A = FormalSums(ZZ);  A.get_action(QQ)
            Right scalar multiplication by Rational Field on Abelian Group of all Formal Finite Sums over Rational Field
            with precomposition on left by Coercion map:
              From: Abelian Group of all Formal Finite Sums over Integer Ring
              To:   Abelian Group of all Formal Finite Sums over Rational Field
            sage: A = FormalSums(QQ);  A.get_action(ZZ)
            Right scalar multiplication by Integer Ring on Abelian Group of all Formal Finite Sums over Rational Field
        """
        if op is operator.mul and isinstance(other, Parent):
            extended = self.base_extend(other)
            if self_is_left:
                action = RightModuleAction(other, extended)
                if extended is not self:
                    action = PrecomposedAction(action, extended._internal_coerce_map_from(self), None)
            else:
                action = LeftModuleAction(other, extended)
                if extended is not self:
                    action = PrecomposedAction(action, None, extended._internal_coerce_map_from(self))
            return action

    def _an_element_(self, check=False, reduce=False):
        """
        EXAMPLES::

            sage: FormalSums(ZZ).an_element()     # indirect test
            1
            sage: FormalSums(QQ).an_element()
            1/2*1
            sage: QQ.an_element()
            1/2
        """
        return self.element_class([(self.base_ring().an_element(), 1)],
                         check=check, reduce=reduce, parent=self)


formal_sums = FormalSums()

# Formal sums now derives from UniqueRepresentation, which makes the
# factory function unnecessary. This is why the name was changed from
# class FormalSums_generic to class FormalSums.
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.structure.formal_sum', 'FormalSums_generic', FormalSums)
