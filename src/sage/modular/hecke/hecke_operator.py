# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Hecke operators
"""
# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.element import AlgebraElement
from sage.structure.richcmp import richcmp, rich_to_bool
from sage.categories.homset import End
import sage.arith.all as arith
from sage.rings.integer import Integer

from . import algebra
from . import morphism


def is_HeckeOperator(x):
    r"""
    Return ``True`` if x is of type HeckeOperator.

    EXAMPLES::

        sage: from sage.modular.hecke.hecke_operator import is_HeckeOperator
        sage: M = ModularSymbols(Gamma0(7), 4)
        sage: is_HeckeOperator(M.T(3))
        doctest:warning...
        DeprecationWarning: the function is_HeckeOperator is deprecated;
        use 'isinstance(..., HeckeOperator)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        True
        sage: is_HeckeOperator(M.T(3) + M.T(5))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_HeckeOperator is deprecated; use 'isinstance(..., HeckeOperator)' instead")
    return isinstance(x, HeckeOperator)


def is_HeckeAlgebraElement(x):
    r"""
    Return ``True`` if x is of type HeckeAlgebraElement.

    EXAMPLES::

        sage: from sage.modular.hecke.hecke_operator import is_HeckeAlgebraElement
        sage: M = ModularSymbols(Gamma0(7), 4)
        sage: is_HeckeAlgebraElement(M.T(3))
        doctest:warning...
        DeprecationWarning: the function is_HeckeAlgebraElement is deprecated;
        use 'isinstance(..., HeckeAlgebraElement)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        True
        sage: is_HeckeAlgebraElement(M.T(3) + M.T(5))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_HeckeAlgebraElement is deprecated; use 'isinstance(..., HeckeAlgebraElement)' instead")
    return isinstance(x, HeckeAlgebraElement)


class HeckeAlgebraElement(AlgebraElement):
    r"""
    Base class for elements of Hecke algebras.
    """
    def __init__(self, parent):
        r"""
        Create an element of a Hecke algebra.

        EXAMPLES::

            sage: R = ModularForms(Gamma0(7), 4).hecke_algebra()
            sage: sage.modular.hecke.hecke_operator.HeckeAlgebraElement(R)  # please don't do this!
            Generic element of a structure
        """
        if not isinstance(parent, algebra.HeckeAlgebra_base):
            raise TypeError("parent (=%s) must be a Hecke algebra" % parent)
        AlgebraElement.__init__(self, parent)

    def domain(self):
        r"""
        The domain of this operator. This is the Hecke module associated to the
        parent Hecke algebra.

        EXAMPLES::

            sage: R = ModularForms(Gamma0(7), 4).hecke_algebra()
            sage: sage.modular.hecke.hecke_operator.HeckeAlgebraElement(R).domain()
            Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7)
             of weight 4 over Rational Field
        """
        return self.parent().module()

    def codomain(self):
        r"""
        The codomain of this operator. This is the Hecke module associated to the
        parent Hecke algebra.

        EXAMPLES::

            sage: R = ModularForms(Gamma0(7), 4).hecke_algebra()
            sage: sage.modular.hecke.hecke_operator.HeckeAlgebraElement(R).codomain()
            Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7)
             of weight 4 over Rational Field
        """
        return self.parent().module()

    def hecke_module_morphism(self):
        """
        Return the endomorphism of Hecke modules defined by the matrix
        attached to this Hecke operator.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(13))
            sage: t = M.hecke_operator(2)
            sage: t
            Hecke operator T_2 on Modular Symbols space of dimension 15 for Gamma_1(13)
             of weight 2 with sign 0 over Rational Field
            sage: t.hecke_module_morphism()
            Hecke module morphism T_2 defined by the matrix
            [ 2  0  0  0  0  0  0  1  0  0  1  0  0  0  0]
            [ 0  2  0  1  0  1  0  0 -1  0  0  0  0  0  1]
            [ 0  1  2  0  0  0  0  0  0  0  0 -1  1  0  0]
            [ 1  0  0  2  0 -1  1  0  1  0 -1  1 -1  0  0]
            [ 0  0  1  0  2  0 -1  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  1 -2  2 -1]
            [ 0  0  0  0  0  2 -1  0 -1  0  0  0  0  1  0]
            [ 0  0  0  0  1  0  0  2  0  0  0  0  0  0 -1]
            [ 0  0  0  0  0  1  0  0 -1  0  2 -1  0  2 -1]
            [ 0  0  0  0  0  1  1  0  0 -1  0  1 -1  2  0]
            [ 0  0  0  0  0  2  0  0 -1 -1  1 -1  0  1  0]
            [ 0  0  0  0  0  1  1  0  1  0  0  0 -1  1  0]
            [ 0  0  0  0  0  1  1  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  1  0  0  1 -1  2  0  0  0 -1]
            [ 0  0  0  0  0  0  0  0  0  1  0 -1  2  0 -1]
            Domain:   Modular Symbols space of dimension 15 for Gamma_1(13) of weight ...
            Codomain: Modular Symbols space of dimension 15 for Gamma_1(13) of weight ...
        """
        try:
            return self.__hecke_module_morphism
        except AttributeError:
            T = self.matrix()
            M = self.domain()
            H = End(M)
            if isinstance(self, HeckeOperator):
                name = "T_%s" % self.index()
            else:
                name = ""
            self.__hecke_module_morphism = morphism.HeckeModuleMorphism_matrix(H, T, name)
            return self.__hecke_module_morphism

    def _add_(self, other):
        """
        Add ``self`` to ``other``.

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: t = M.hecke_operator(2)
            sage: t
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: t + t # indirect doctest
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [ 6  0 -2]
            [ 0 -4  0]
            [ 0  0 -4]

        We can also add Hecke operators with different indexes::

            sage: M = ModularSymbols(Gamma1(6),4)
            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t2 + t3
            Hecke operator on Modular Symbols space of dimension 6 for Gamma_1(6) of weight 4 with sign 0 over Rational Field defined by:
            [   35     0     0   8/5   8/5 -16/5]
            [    4    28     0 -19/5 -19/5  38/5]
            [   18     0     9    -6     8    -2]
            [    0    18     4 -23/5 -13/5  46/5]
            [    0    18     4   2/5 -38/5  46/5]
            [    0    18     4   2/5 -13/5  21/5]
            sage: (t2 - t3).charpoly('x')
            x^6 + 36*x^5 + 104*x^4 - 3778*x^3 + 7095*x^2 - 3458*x
        """
        return self.parent()(self.matrix() + other.matrix(), check=False)

    def __call__(self, x):
        """
        Apply this Hecke operator to `x`.

        EXAMPLES::

            sage: M = ModularSymbols(11); t2 = M.hecke_operator(2)
            sage: t2(M.gen(0))
            3*(1,0) - (1,9)

        ::

            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t3(t2(M.gen(0)))
            12*(1,0) - 2*(1,9)
            sage: (t3*t2)(M.gen(0))
            12*(1,0) - 2*(1,9)
        """
        T = self.hecke_module_morphism()
        return T(x)

    def __rmul__(self, left):
        """
        EXAMPLES::

            sage: M = ModularSymbols(11); t2 = M.hecke_operator(2)
            sage: 2*t2
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [ 6  0 -2]
            [ 0 -4  0]
            [ 0  0 -4]
        """
        return self.parent()(left * self.matrix())

    def _sub_(self, other):
        """
        Compute the difference of ``self`` and ``other``, where ``other`` has
        already been coerced into the parent of ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(6),4)
            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t2 - t3 # indirect doctest
            Hecke operator on Modular Symbols space of dimension 6 for Gamma_1(6) of weight 4 with sign 0 over Rational Field defined by:
            [  -19     0     0  -4/5  -4/5   8/5]
            [    4   -26     0  17/5  17/5 -34/5]
            [  -18     0     7 -18/5  12/5   6/5]
            [    0   -18     4   3/5  23/5 -26/5]
            [    0   -18     4  -2/5  28/5 -26/5]
            [    0   -18     4  -2/5  23/5 -21/5]
        """
        return self.parent()(self.matrix() - other.matrix(), check=False)

    def apply_sparse(self, x):
        """
        Apply this Hecke operator to x, where we avoid computing the matrix
        of x if possible.

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: T = M.hecke_operator(23)
            sage: T.apply_sparse(M.gen(0))
            24*(1,0) - 5*(1,9)
        """
        if x not in self.domain():
            raise TypeError("x (=%s) must be in %s" % (x, self.domain()))
        # Generic implementation which doesn't actually do anything
        # special regarding sparseness.  Override this for speed.
        T = self.hecke_module_morphism()
        return T(x)

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of this Hecke operator.

        INPUT:

        - ``var`` -- string (default: ``'x'``)

        OUTPUT: a monic polynomial in the given variable

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(6),4)
            sage: M.hecke_operator(2).charpoly('x')
            x^6 - 14*x^5 + 29*x^4 + 172*x^3 - 124*x^2 - 320*x + 256
        """
        return self.matrix().charpoly(var)

    def decomposition(self):
        """
        Decompose the Hecke module under the action of this Hecke
        operator.

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: t2 = M.hecke_operator(2)
            sage: t2.decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            ]

        ::

            sage: M = ModularSymbols(33, sign=1).new_submodule()
            sage: T = M.hecke_operator(2)
            sage: T.decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field
            ]
        """
        try:
            return self.__decomposition
        except AttributeError:
            pass
        if isinstance(self, HeckeOperator) and \
           arith.gcd(self.index(), self.domain().level()) == 1:
            D = self.hecke_module_morphism().decomposition(is_diagonalizable=True)
        else:
            # TODO: There are other weaker hypotheses that imply diagonalizability.
            D = self.hecke_module_morphism().decomposition()
        D.sort()
        D.set_immutable()
        self.__decomposition = D
        return D

    def det(self):
        """
        Return the determinant of this Hecke operator.

        EXAMPLES::

            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.det()
            100
        """
        return self.hecke_module_morphism().det()

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of this
        Hecke operator.

        EXAMPLES::

            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.fcp('x')
            (x - 4) * (x^2 - 5)^2
        """
        return self.hecke_module_morphism().fcp(var)

    def image(self):
        """
        Return the image of this Hecke operator.

        EXAMPLES::

            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.fcp('x')
            (x - 4) * (x^2 - 5)^2
            sage: T.image()
            Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T-4).image()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T**2-5).image()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
        """
        return self.hecke_module_morphism().image()

    def kernel(self):
        """
        Return the kernel of this Hecke operator.

        EXAMPLES::

            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.fcp('x')
            (x - 4) * (x^2 - 5)^2
            sage: T.kernel()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T-4).kernel()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T**2-5).kernel()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
        """
        return self.hecke_module_morphism().kernel()

    def trace(self):
        """
        Return the trace of this Hecke operator.

        ::

            sage: M = ModularSymbols(1,12)
            sage: T = M.hecke_operator(2)
            sage: T.trace()
            2001
        """
        return self.hecke_module_morphism().trace()

    def __getitem__(self, ij):
        """
        EXAMPLES::

            sage: M = ModularSymbols(1,12)
            sage: T = M.hecke_operator(2).matrix_form()
            sage: T[0,0]
            -24
        """
        return self.matrix()[ij]


class HeckeAlgebraElement_matrix(HeckeAlgebraElement):
    r"""
    An element of the Hecke algebra represented by a matrix.
    """
    def __init__(self, parent, A):
        r"""
        Initialise an element from a matrix. This *must* be over the base ring
        of ``self`` and have the right size.

        This is a bit overkill as similar checks will be performed by the call
        and coerce methods of the parent of self, but it can't hurt to be
        paranoid. Any fancy coercion / base_extension / etc happens there, not
        here.

        TESTS::

            sage: T = ModularForms(Gamma0(7), 4).hecke_algebra()
            sage: M = sage.modular.hecke.hecke_operator.HeckeAlgebraElement_matrix(T, matrix(QQ,3,[2,3,0,1,2,3,7,8,9])); M
            Hecke operator on Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) of weight 4 over Rational Field defined by:
            [2 3 0]
            [1 2 3]
            [7 8 9]
            sage: loads(dumps(M)) == M
            True
            sage: sage.modular.hecke.hecke_operator.HeckeAlgebraElement_matrix(T, matrix(Integers(2),3,[2,3,0,1,2,3,7,8,9]))
            Traceback (most recent call last):
            ...
            TypeError: base ring of matrix (Ring of integers modulo 2) does not match base ring of space (Rational Field)
            sage: sage.modular.hecke.hecke_operator.HeckeAlgebraElement_matrix(T, matrix(QQ,2,[2,3,0,1]))
            Traceback (most recent call last):
            ...
            TypeError: A must be a square matrix of rank 3
        """
        HeckeAlgebraElement.__init__(self, parent)
        from sage.structure.element import Matrix
        if not isinstance(A, Matrix):
            raise TypeError("A must be a matrix")
        if not A.base_ring() == self.parent().base_ring():
            raise TypeError("base ring of matrix (%s) does not match base ring of space (%s)" % (A.base_ring(), self.parent().base_ring()))
        if not A.nrows() == A.ncols() == self.parent().module().rank():
            raise TypeError("A must be a square matrix of rank %s" % self.parent().module().rank())
        self.__matrix = A

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` to ``other``, where the coercion model has already ensured
        that ``other`` has the same parent as ``self``.

        EXAMPLES::

            sage: T = ModularForms(SL2Z, 12).hecke_algebra()
            sage: m = T(matrix(QQ, 2, [1,2,0,1]), check=False); n = T.hecke_operator(14)
            sage: m == n
            False
            sage: m == n.matrix_form()
            False
            sage: n.matrix_form() == T(matrix(QQ, 2, [401856,0,0,4051542498456]), check=False)
            True
        """
        if not isinstance(other, HeckeAlgebraElement_matrix):
            if isinstance(other, HeckeOperator):
                return richcmp(self, other.matrix_form(), op)
            else:
                raise RuntimeError("Bug in coercion code")  # can't get here

        return richcmp(self.__matrix, other.__matrix, op)

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(1,12)
            sage: M.hecke_operator(2).matrix_form()._repr_()
            'Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field defined by:\n[ -24    0    0]\n[   0  -24    0]\n[4860    0 2049]'
            sage: ModularForms(Gamma0(100)).hecke_operator(4).matrix_form()._repr_()
            'Hecke operator on Modular Forms space of dimension 24 for Congruence Subgroup Gamma0(100) of weight 2 over Rational Field defined by:\n24 x 24 dense matrix over Rational Field'
        """
        return "Hecke operator on %s defined by:\n%r" % (self.parent().module(), self.__matrix)

    def _latex_(self):
        r"""
        Latex representation of ``self`` (just prints the matrix).

        EXAMPLES::

            sage: M = ModularSymbols(1,12)
            sage: M.hecke_operator(2).matrix_form()._latex_()
            '\\left(\\begin{array}{rrr}\n-24 & 0 & 0 \\\\\n0 & -24 & 0 \\\\\n4860 & 0 & 2049\n\\end{array}\\right)'
        """
        return self.__matrix._latex_()

    def matrix(self):
        """
        Return the matrix that defines this Hecke algebra element.

        EXAMPLES::

            sage: M = ModularSymbols(1,12)
            sage: T = M.hecke_operator(2).matrix_form()
            sage: T.matrix()
            [ -24    0    0]
            [   0  -24    0]
            [4860    0 2049]
        """
        return self.__matrix

    def _mul_(self, other):
        r"""
        Multiply ``self`` by ``other`` (which has already been coerced into an element
        of the parent of ``self``).

        EXAMPLES::

            sage: M = ModularSymbols(1,12)
            sage: T = M.hecke_operator(2).matrix_form()
            sage: T * T # indirect doctest
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field defined by:
            [    576       0       0]
            [      0     576       0]
            [9841500       0 4198401]
        """
        return self.parent()(other.matrix() * self.matrix(), check=False)


class DiamondBracketOperator(HeckeAlgebraElement_matrix):
    r"""
    The diamond bracket operator `\langle d \rangle` for some `d \in \ZZ /
    N\ZZ` (which need not be a unit, although if it is not, the operator will
    be zero).
    """
    def __init__(self, parent, d):
        r"""
        Standard init function.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(5),6)
            sage: d = M.diamond_bracket_operator(2); d # indirect doctest
            Diamond bracket operator <2> on Modular Symbols space of dimension 10 for Gamma_1(5) of weight 6 with sign 0 over Rational Field
            sage: type(d)
            <class 'sage.modular.hecke.hecke_operator.DiamondBracketOperator'>
            sage: d.matrix()
            [    0     1     0     0     0     0     0     0     0     0]
            [    1     0     0     0     0     0     0     0     0     0]
            [    0     0     0     0     0     0     1     0     0     0]
            [    0     0     0     0     0     0     0     0     0     1]
            [    0     0     0     0     0     0     0     1     0     0]
            [    0     0 17/16 11/16  -3/4    -1 17/16  -3/4     0 11/16]
            [    0     0     1     0     0     0     0     0     0     0]
            [    0     0     0     0     1     0     0     0     0     0]
            [    0     0  -1/2   1/2     1     0  -1/2     1    -1   1/2]
            [    0     0     0     1     0     0     0     0     0     0]
            sage: d**4 == 1
            True
        """
        self.__d = d
        A = parent.diamond_bracket_matrix(d)
        HeckeAlgebraElement_matrix.__init__(self, parent, A)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: ModularSymbols(Gamma1(5), 6).diamond_bracket_operator(2)._repr_()
            'Diamond bracket operator <2> on Modular Symbols space of dimension 10 for Gamma_1(5) of weight 6 with sign 0 over Rational Field'
        """
        return "Diamond bracket operator <%s> on %s" % (self.__d, self.domain())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(ModularSymbols(Gamma1(5), 12).diamond_bracket_operator(2)) # indirect doctest
            \langle 2 \rangle
        """
        return r"\langle %s \rangle" % self.__d


class HeckeOperator(HeckeAlgebraElement):
    r"""
    The Hecke operator `T_n` for some `n` (which need not be coprime to the
    level). The matrix is not computed until it is needed.
    """
    def __init__(self, parent, n):
        """
        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: H = M.hecke_operator(2005); H
            Hecke operator T_2005 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: H == loads(dumps(H))
            True

        We create a Hecke operator of large index (greater than 32 bits)::

            sage: M1 =  ModularSymbols(21,2)
            sage: M1.hecke_operator(13^9)
            Hecke operator T_10604499373 on Modular Symbols space of dimension 5 for Gamma_0(21) of weight 2 with sign 0 over Rational Field
        """
        HeckeAlgebraElement.__init__(self, parent)
        if not isinstance(n, (int, Integer)):
            raise TypeError("n must be an int")
        self.__n = int(n)

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` (where the coercion model has already ensured
        that ``self`` and ``other`` have the same parent). Hecke operators on the same
        space compare as equal if and only if their matrices are equal, so we
        check if the indices are the same and if not we compute the matrices
        (which is potentially expensive).

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(7), 4)
            sage: m = M.hecke_operator(3)
            sage: m == m
            True
            sage: m == 2*m
            False
            sage: m == M.hecke_operator(5)
            False

        These last two tests involve a coercion::

            sage: m == m.matrix_form()
            True
            sage: m == m.matrix()
            False
        """
        if not isinstance(other, HeckeOperator):
            if isinstance(other, HeckeAlgebraElement_matrix):
                return richcmp(self.matrix_form(), other, op)
            else:
                raise RuntimeError("Bug in coercion code")  # can't get here

        if self.__n == other.__n:
            return rich_to_bool(op, 0)
        return richcmp(self.matrix(), other.matrix(), op)

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(7), 4).hecke_operator(6)._repr_()
            'Hecke operator T_6 on Modular Symbols space of dimension 4 for Gamma_0(7) of weight 4 with sign 0 over Rational Field'
        """
        return "Hecke operator T_%s on %s" % (self.__n, self.domain())

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(7), 4).hecke_operator(6)._latex_()
            'T_{6}'
        """
        return "T_{%s}" % self.__n

    def _mul_(self, other):
        r"""
        Multiply this Hecke operator by another element of the same algebra.

        If the other element is of the form `T_m` for some m, we check
        whether the product is equal to `T_{mn}` and return that; if
        the product is not (easily seen to be) of the form `T_{mn}`,
        then we calculate the product of the two matrices and return a
        Hecke algebra element defined by that.

        EXAMPLES: We create the space of modular symbols of level
        `11` and weight `2`, then compute `T_2`
        and `T_3` on it, along with their composition.

        ::

            sage: M = ModularSymbols(11)
            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t2*t3 # indirect doctest
            Hecke operator T_6 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: t3.matrix() * t2.matrix()
            [12  0 -2]
            [ 0  2  0]
            [ 0  0  2]
            sage: (t2*t3).matrix()
            [12  0 -2]
            [ 0  2  0]
            [ 0  0  2]

        When we compute `T_2^5` the result is not (easily seen to
        be) a Hecke operator of the form `T_n`, so it is returned
        as a Hecke module homomorphism defined as a matrix::

            sage: t2**5
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [243   0 -55]
            [  0 -32   0]
            [  0   0 -32]
        """
        if isinstance(other, HeckeOperator) and other.parent() == self.parent():
            n = None
            if arith.gcd(self.__n, other.__n) == 1:
                n = self.__n * other.__n
            else:
                P = set(arith.prime_divisors(self.domain().level()))
                if P.issubset(set(arith.prime_divisors(self.__n))) and \
                   P.issubset(set(arith.prime_divisors(other.__n))):
                    n = self.__n * other.__n
            if n:
                return HeckeOperator(self.parent(), n)
        # otherwise
        return self.matrix_form() * other

    def index(self):
        """
        Return the index of this Hecke operator, i.e., if this Hecke
        operator is `T_n`, return the int `n`.

        EXAMPLES::

            sage: T = ModularSymbols(11).hecke_operator(17)
            sage: T.index()
            17
        """
        return self.__n

    def matrix(self, *args, **kwds):
        """
        Return the matrix underlying this Hecke operator.

        EXAMPLES::

            sage: T = ModularSymbols(11).hecke_operator(17)
            sage: T.matrix()
            [18  0 -4]
            [ 0 -2  0]
            [ 0  0 -2]
        """
        try:
            return self.__matrix
        except AttributeError:
            self.__matrix = self.parent().hecke_matrix(self.__n, *args, **kwds)
            return self.__matrix

    def matrix_form(self):
        """
        Return the matrix form of this element of a Hecke algebra.

        ::

            sage: T = ModularSymbols(11).hecke_operator(17)
            sage: T.matrix_form()
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [18  0 -4]
            [ 0 -2  0]
            [ 0  0 -2]
        """
        try:
            return self.__matrix_form
        except AttributeError:
            self.__matrix_form = self.parent()(self.matrix(), check=False)
            return self.__matrix_form
