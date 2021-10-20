r"""
Elements of finitely generated modules over a PID

AUTHOR:
    - William Stein, 2009
"""

#*****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp

# This adds extra maybe-not-necessary checks in the code, but could
# slow things down.  It can impact what happens in more than just this
# file.
DEBUG=True


class FGP_Element(ModuleElement):
    """
    An element of a finitely generated module over a PID.

    INPUT:

    - ``parent`` -- parent module M

    - ``x`` -- element of M.V()

    EXAMPLES::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
        sage: Q = V/W
        sage: x = Q(V.0-V.1); x #indirect doctest
        (0, 9)
        sage: isinstance(x, sage.modules.fg_pid.fgp_element.FGP_Element)
        True
        sage: type(x)
        <class 'sage.modules.fg_pid.fgp_module.FGP_Module_class_with_category.element_class'>
        sage: x is Q(x)
        True
        sage: x.parent() is Q
        True

    TESTS::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
        sage: loads(dumps(Q.0)) == Q.0
        True
    """
    def __init__(self, parent, x, check=DEBUG):
        """
        INPUT:

        - ``parent`` -- parent module M

        - ``x`` -- element of M.V()

        - ``check`` -- (default: True) if True, verify that x in M.V()

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: x = Q(V.0-V.1); type(x)
            <class 'sage.modules.fg_pid.fgp_module.FGP_Module_class_with_category.element_class'>
            sage: isinstance(x,sage.modules.fg_pid.fgp_element.FGP_Element)
            True

        For full documentation, see :class:`FGP_Element`.
        """
        if check:
            assert x in parent.V(), 'The argument x='+str(x)+' is not in the covering module!'
        ModuleElement.__init__(self, parent)
        self._x = x

    def lift(self):
        """
        Lift self to an element of V, where the parent of self is the quotient module V/W.

        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.0
            (1, 0)
            sage: Q.1
            (0, 1)
            sage: Q.0.lift()
            (0, 6, 1)
            sage: Q.1.lift()
            (0, -2, 0)
            sage: x = Q(V.0); x
            (0, 8)
            sage: x.lift()
            (1/2, 0, 0)
            sage: x == 8*Q.1
            True
            sage: x.lift().parent() == V
            True

        A silly version of the integers modulo 100::

            sage: A = (ZZ^1)/span([[100]], ZZ); A
            Finitely generated module V/W over Integer Ring with invariants (100)
            sage: x = A([5]); x
            (5)
            sage: v = x.lift(); v
            (5)
            sage: v.parent()
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
        """
        return self._x


    def __neg__(self):
        """
        EXAMPLES::

            sage: V1 = ZZ^2; W1 = V1.span([[1,2],[3,4]]); A1 = V1/W1; A1
            Finitely generated module V/W over Integer Ring with invariants (2)
            sage: -A1.0
            (1)
            sage: -A1.0 == A1.0           # order 2
            True
        """
        P = self.parent()
        return P.element_class(P, -self._x)


    def _add_(self, other):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0; x
            (1, 0)
            sage: y = Q.1; y
            (0, 1)
            sage: x + y                        # indirect doctest
            (1, 1)
            sage: x + x + x + x
            (0, 0)
            sage: x + 0
            (1, 0)
            sage: 0 + x
            (1, 0)

        We test canonical coercion from V and W.

            sage: Q.0 + V.0
            (1, 8)
            sage: V.0 + Q.0
            (1, 8)
            sage: W.0 + Q.0
            (1, 0)
            sage: W.0 + Q.0 == Q.0
            True
        """
        P = self.parent()
        return P.element_class(P, self._x + other._x)


    def _sub_(self, other):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0; x
            (1, 0)
            sage: y = Q.1; y
            (0, 1)
            sage: x - y                                # indirect doctest
            (1, 11)
            sage: x - x
            (0, 0)
        """
        P = self.parent()
        return P.element_class(P, self._x - other._x)


    def _rmul_(self, c):
        """
        Multiplication by a scalar from the left (``self`` is on the right).

        INPUT:

        - ``c`` -- an element of ``self.parent().base_ring()``.

        OUTPUT:

        The product ``c * self`` as a new instance of a module
        element.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0; x
            (1, 0)
            sage: 2 * x                            # indirect doctest
            (2, 0)
            sage: x._rmul_(4)
            (0, 0)
            sage: V = V.base_extend(QQ); W = V.span([2*V.0+4*V.1])
            sage: Q = V/W; Q
            Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 0]
            sage: x = Q.0; x
            (1, 0)
            sage: (1/2) * x                            # indirect doctest
            (1/2, 0)
            sage: x._rmul_(1/4)
            (1/4, 0)
        """
        # print "_rmul_"
        P = self.parent()
        return P.element_class(P, self._x._rmul_(c))

    def _lmul_(self, s):
        """
        Multiplication by a scalar from the right (``self`` is on the left).

        INPUT:

        - ``c`` -- an element of ``self.parent().base_ring()``.

        OUTPUT:

        The product ``self * c`` as a new instance of a module
        element.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0; x
            (1, 0)
            sage: x * 2                          # indirect doctest
            (2, 0)
            sage: x._lmul_(4)
            (0, 0)
            sage: V = V.base_extend(QQ); W = V.span([2*V.0+4*V.1])
            sage: Q = V/W; Q
            Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 0]
            sage: x = Q.0; x
            (1, 0)
            sage: x * (1/2)                            # indirect doctest
            (1/2, 0)
            sage: x._lmul_(1/4)
            (1/4, 0)
        """
        # print '_lmul_'
        P = self.parent()
        return P.element_class(P, self._x._lmul_(s))


    def _repr_(self):
        """

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q(V.1)._repr_()
            '(0, 11)'
        """
        return repr(self.vector())


    def __getitem__(self, *args):
        """
        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0 + 3*Q.1; x
            (1, 3)
            sage: x[0]
            1
            sage: x[1]
            3
            sage: x[-1]
            3
        """
        return self.vector().__getitem__(*args)

    def vector(self):
        """
        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0 + 3*Q.1; x
            (1, 3)
            sage: x.vector()
            (1, 3)
            sage: tuple(x)
            (1, 3)
            sage: list(x)
            [1, 3]
            sage: x.vector().parent()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
        """
        try:
            return self.__vector
        except AttributeError:
            self.__vector = self.parent().coordinate_vector(self, reduce=True)
            self.__vector.set_immutable()
            return self.__vector

    def __hash__(self):
        r"""
        TESTS::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ)
            sage: W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: x = Q.0 + 3*Q.1
            sage: hash(x) == hash((1,3))
            True

            sage: A = AdditiveAbelianGroup([3])
            sage: hash(A.an_element()) == hash((1,))
            True
        """
        return hash(self.vector())

    def _vector_(self, base_ring=None):
        """
        Support for conversion to vectors.

        INPUT:

        - ``base_ring`` -- the desired base ring of the vector.

        OUTPUT:
        
        A vector over the base ring.

        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0 + 3*Q.1
            sage: vector(x)
            (1, 3)
            sage: vector(CDF, x)
            (1.0, 3.0)

        TESTS::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ)
            sage: W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: x = Q.0 + 3*Q.1
            sage: vector(x).is_mutable()
            True
            sage: vector(CDF,x).is_mutable()
            True
        """
        v = self.vector()
        if base_ring is None or v.base_ring() is base_ring:
            return v.__copy__()
        else:
            return v.change_ring(base_ring)

    def _richcmp_(self, right, op):
        """
        Compare self and right.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: x = Q.0; x
            (1, 0)
            sage: y = Q.1; y
            (0, 1)
            sage: x == y
            False
            sage: x == x
            True
            sage: x + x == 2*x
            True
        """
        return richcmp(self.vector(), right.vector(), op)

    def additive_order(self):
        """
        Return the additive order of this element.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.0.additive_order()
            4
            sage: Q.1.additive_order()
            12
            sage: (Q.0+Q.1).additive_order()
            12
            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (12, 0)
            sage: Q.0.additive_order()
            12
            sage: type(Q.0.additive_order())
            <class 'sage.rings.integer.Integer'>
            sage: Q.1.additive_order()
            +Infinity
        """
        Q = self.parent()
        I = Q.invariants()
        v = self.vector()

        from sage.rings.infinity import infinity
        from sage.rings.finite_rings.integer_mod import Mod
        from sage.rings.integer import Integer
        from sage.arith.all import lcm
        n = Integer(1)
        for i, a in enumerate(I):
            if a == 0:
                if v[i] != 0:
                    return infinity
            else:
                n = lcm(n, Mod(v[i],a).additive_order())
        return n
