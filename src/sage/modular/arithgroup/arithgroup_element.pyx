# sage_setup: distribution = sagemath-schemes
"""
Elements of arithmetic subgroups
"""
# ##############################################################################
#
#       Copyright (C) 2009, The Sage Group -- https://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
#
# ##############################################################################

from sage.structure.element cimport MultiplicativeGroupElement
from sage.structure.richcmp cimport richcmp
from sage.rings.integer_ring import ZZ

from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense as MatrixClass
# FIXME: This was Matrix_integer_dense; changed for modularization to avoid flint dep

M2Z = MatrixSpace(ZZ, 2, implementation='generic')


cdef class ArithmeticSubgroupElement(MultiplicativeGroupElement):
    r"""
    An element of the group `\SL_2(\ZZ)`, i.e. a 2x2 integer matrix of
    determinant 1.
    """

    cdef MatrixClass __x

    def __init__(self, parent, x, check=True):
        """
        Create an element of an arithmetic subgroup.

        INPUT:

        - ``parent`` -- an arithmetic subgroup

        - ``x`` -- data defining a 2x2 matrix over ZZ
          which lives in ``parent``

        - ``check`` -- if ``True``, check that parent is an arithmetic
          subgroup, and that `x` defines a matrix of determinant `1`

        We tend not to create elements of arithmetic subgroups that are not
        SL2Z, in order to avoid coercion issues (that is, the other arithmetic
        subgroups are "facade parents").

        EXAMPLES::

            sage: G = Gamma0(27)
            sage: sage.modular.arithgroup.arithgroup_element.ArithmeticSubgroupElement(G, [4,1,27,7])
            [ 4  1]
            [27  7]
            sage: sage.modular.arithgroup.arithgroup_element.ArithmeticSubgroupElement(Integers(3), [4,1,27,7])
            Traceback (most recent call last):
            ...
            TypeError: parent (= Ring of integers modulo 3) must be an arithmetic subgroup
            sage: sage.modular.arithgroup.arithgroup_element.ArithmeticSubgroupElement(G, [2,0,0,2])
            Traceback (most recent call last):
            ...
            TypeError: matrix must have determinant 1
            sage: sage.modular.arithgroup.arithgroup_element.ArithmeticSubgroupElement(G, [2,0,0,2], check=False)
            [2 0]
            [0 2]
            sage: x = Gamma0(11)([2,1,11,6])
            sage: TestSuite(x).run()

            sage: x = Gamma0(5).0
            sage: SL2Z(x)
            [1 1]
            [0 1]
            sage: x in SL2Z
            True
        """
        if check:
            from .arithgroup_generic import ArithmeticSubgroup
            if not isinstance(parent, ArithmeticSubgroup):
                raise TypeError("parent (= %s) must be an arithmetic subgroup" % parent)

            x = M2Z(x, copy=True, coerce=True)
            if x.determinant() != 1:
                raise TypeError("matrix must have determinant 1")
        else:
            x = M2Z(x, copy=True, coerce=False)
            # Getting rid of this would result in a small speed gain for
            # arithmetic operations, but it would have the disadvantage of
            # causing strange and opaque errors when inappropriate data types
            # are used with "check=False".

        x.set_immutable()
        MultiplicativeGroupElement.__init__(self, parent)
        self.__x = x

    def __setstate__(self, state):
        r"""
        For unpickling objects pickled with the old ArithmeticSubgroupElement class.

        EXAMPLES::

            sage: si = unpickle_newobj(sage.modular.arithgroup.arithgroup_element.ArithmeticSubgroupElement, ())
            sage: x = matrix(ZZ,2,[1,1,0,1])
            sage: unpickle_build(si, (Gamma0(13), {'_ArithmeticSubgroupElement__x': x}))
        """
        from .congroup_sl2z import SL2Z
        _, kwdict = state
        self._parent = SL2Z
        if '_ArithmeticSubgroupElement__x' in kwdict:
            self.__x = M2Z(kwdict['_ArithmeticSubgroupElement__x'])
        elif '_CongruenceSubgroupElement__x' in kwdict:
            self.__x = M2Z(kwdict['_CongruenceSubgroupElement__x'])
        else:
            raise ValueError("Don't know how to unpickle %s" % repr(state))

    def __iter__(self):
        """
        EXAMPLES::

            sage: Gamma0(2).0
            [1 1]
            [0 1]
            sage: list(Gamma0(2).0)
            [1, 1, 0, 1]

        Warning: this is different from the iteration on the matrix::

            sage: list(Gamma0(2).0.matrix())
            [(1, 1), (0, 1)]
        """
        yield self.__x[0, 0]
        yield self.__x[0, 1]
        yield self.__x[1, 0]
        yield self.__x[1, 1]

    def __repr__(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: Gamma1(5)([6,1,5,1]).__repr__()
            '[6 1]\n[5 1]'
        """
        return "%s" % self.__x

    def _latex_(self):
        r"""
        Return latex representation of ``self``.

        EXAMPLES::

            sage: Gamma1(5)([6,1,5,1])._latex_()
            '\\left(\\begin{array}{rr}\n6 & 1 \\\\\n5 & 1\n\\end{array}\\right)'
        """
        return '%s' % self.__x._latex_()

    cpdef _richcmp_(self, right_r, int op):
        """
        Compare ``self`` to ``right``, where ``right`` is guaranteed to have
        the same parent as ``self``.

        EXAMPLES::

            sage: x = Gamma0(18)([19,1,18,1])
            sage: x == 3
            False
            sage: x == x
            True

            sage: x = Gamma0(5)([1,1,0,1])
            sage: y = Gamma0(5)([1,4,0,1])
            sage: x == 0
            False
            sage: x == y
            False
            sage: x != y
            True

        This once caused a segfault (see :issue:`5443`)::

            sage: r,s,t,u,v = map(SL2Z, [[1, 1, 0, 1], [-1, 0, 0, -1], [1, -1, 0, 1], [1, -1, 2, -1], [-1, 1, -2, 1]])
            sage: v == s*u
            True
            sage: s*u == v
            True
        """
        cdef ArithmeticSubgroupElement right = <ArithmeticSubgroupElement>right_r
        return richcmp(self.__x, right.__x, op)

    def __bool__(self):
        r"""
        Return ``True``, since the ``self`` lives in SL(2,\Z), which does not
        contain the zero matrix.

        EXAMPLES::

            sage: x = Gamma0(5)([1,1,0,1])
            sage: bool(x)
            True
        """
        return True

    cpdef _mul_(self, right):
        """
        Return ``self * right``.

        EXAMPLES::

            sage: x = Gamma0(7)([1,0,7,1]) * Gamma0(7)([3,2,7,5]) ; x # indirect doctest
            [ 3  2]
            [28 19]
            sage: x.parent()
            Modular Group SL(2,Z)

        We check that :issue:`5048` is fixed::

            sage: a = Gamma0(10).1 * Gamma0(5).2; a # random
            sage: a.parent()
            Modular Group SL(2,Z)
        """
        return self.__class__(self.parent(), self.__x * (<ArithmeticSubgroupElement> right).__x, check=False)

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: Gamma0(11)([1,-1,0,1]).__invert__()
            [1 1]
            [0 1]
        """
        return self._parent(
                [self.__x.get_unsafe(1, 1), -self.__x.get_unsafe(0, 1),
                 -self.__x.get_unsafe(1, 0), self.__x.get_unsafe(0, 0)],
                check=False)

    def matrix(self):
        """
        Return the matrix corresponding to ``self``.

        EXAMPLES::

            sage: x = Gamma1(3)([4,5,3,4]) ; x
            [4 5]
            [3 4]
            sage: x.matrix()
            [4 5]
            [3 4]
            sage: type(x.matrix())
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
        """
        return self.__x

    def determinant(self):
        """
        Return the determinant of ``self``, which is always 1.

        EXAMPLES::

            sage: Gamma0(691)([1,0,691,1]).determinant()
            1
        """
        return ZZ.one()

    def det(self):
        """
        Return the determinant of ``self``, which is always 1.

        EXAMPLES::

            sage: Gamma1(11)([12,11,-11,-10]).det()
            1
        """
        return self.determinant()

    def a(self):
        """
        Return the upper left entry of ``self``.

        EXAMPLES::

            sage: Gamma0(13)([7,1,13,2]).a()
            7
        """
        return self.__x.get_unsafe(0, 0)

    def b(self):
        """
        Return the upper right entry of ``self``.

        EXAMPLES::

            sage: Gamma0(13)([7,1,13,2]).b()
            1
        """
        return self.__x.get_unsafe(0, 1)

    def c(self):
        """
        Return the lower left entry of ``self``.

        EXAMPLES::

            sage: Gamma0(13)([7,1,13,2]).c()
            13
        """
        return self.__x.get_unsafe(1, 0)

    def d(self):
        """
        Return the lower right entry of ``self``.

        EXAMPLES::

            sage: Gamma0(13)([7,1,13,2]).d()
            2
        """
        return self.__x.get_unsafe(1, 1)

    def acton(self, z):
        """
        Return the result of the action of ``self`` on z as a fractional linear
        transformation.

        EXAMPLES::

            sage: G = Gamma0(15)
            sage: g = G([1, 2, 15, 31])

        An example of g acting on a symbolic variable::

            sage: z = var('z')                                                          # needs sage.symbolic
            sage: g.acton(z)                                                            # needs sage.symbolic
            (z + 2)/(15*z + 31)

        An example involving the Gaussian numbers::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: g.acton(i)
            1/1186*i + 77/1186

        An example with complex numbers::

            sage: C.<i> = ComplexField()
            sage: g.acton(i)
            0.0649241146711636 + 0.000843170320404721*I

        An example with the cusp infinity::

            sage: g.acton(infinity)
            1/15

        An example which maps a finite cusp to infinity::

            sage: g.acton(-31/15)
            +Infinity

        Note that when acting on instances of cusps the return value
        is still a rational number or infinity (Note the presence of
        '+', which does not show up for cusp instances)::

            sage: g.acton(Cusp(-31/15))
            +Infinity

        TESTS:

        We cover the remaining case, i.e., infinity mapped to infinity::

            sage: G([1, 4, 0, 1]).acton(infinity)
            +Infinity
        """
        from sage.rings.infinity import InfinityElement, infinity
        if isinstance(z, InfinityElement):
            if self.c() != 0:
                return self.a() / self.c()
            else:
                return infinity
        if hasattr(z, 'denominator') and hasattr(z, 'numerator'):
            p = z.numerator()
            q = z.denominator()
            P = self.a() * p + self.b() * q
            Q = self.c() * p + self.d() * q
            if not Q and P:
                return infinity
            else:
                return P / Q
        return (self.a() * z + self.b()) / (self.c() * z + self.d())

    def __getitem__(self, q):
        r"""
        Fetch entries by direct indexing.

        EXAMPLES::

            sage: SL2Z([3,2,1,1])[0,0]
            3
        """
        return self.__x[q]

    def __hash__(self):
        r"""
        Return a hash value.

        EXAMPLES::

            sage: hash(SL2Z.0)
            -8192788425652673914  # 64-bit
            -1995808122           # 32-bit
        """
        return hash(self.__x)

    def __reduce__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: (SL2Z.1).__reduce__()
            (Modular Group SL(2,Z), (
            [1 1]
            [0 1]
            ))
        """
        from .congroup_sl2z import SL2Z
        return SL2Z, (self.__x,)

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of this element.

        EXAMPLES::

            sage: SL2Z.one().multiplicative_order()
            1
            sage: SL2Z([-1,0,0,-1]).multiplicative_order()
            2
            sage: s,t = SL2Z.gens()
            sage: ((t^3*s*t^2) * s * ~(t^3*s*t^2)).multiplicative_order()
            4
            sage: (t^3 * s * t * t^-3).multiplicative_order()
            6
            sage: (t^3 * s * t * s * t^-2).multiplicative_order()
            3
            sage: SL2Z([2,1,1,1]).multiplicative_order()
            +Infinity
            sage: SL2Z([-2,1,1,-1]).multiplicative_order()
            +Infinity
        """
        m = self.matrix()

        if m.is_one():
            return ZZ(1)
        elif (-m).is_one():
            return ZZ(2)

        t = m.trace()
        if t <= -2 or t >= 2:
            from sage.rings.infinity import infinity
            return infinity
        elif t == 0:
            return ZZ(4)
        elif t == 1:
            return ZZ(6)
        elif t == -1:
            return ZZ(3)
