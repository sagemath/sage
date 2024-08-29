# sage.doctest: needs sage.combinat sage.modules
"""
Free algebra quotient elements

AUTHORS:
    - William Stein (2011-11-19): improved doctest coverage to 100%
    - David Kohel (2005-09): initial version
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.repr import repr_lincomb
from sage.structure.element import RingElement, AlgebraElement
from sage.structure.parent_gens import localvars
from sage.structure.richcmp import richcmp
from sage.rings.integer import Integer
from sage.modules.free_module_element import FreeModuleElement
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.algebras.free_algebra_element import FreeAlgebraElement


def is_FreeAlgebraQuotientElement(x):
    """
    EXAMPLES::

        sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
        sage: sage.algebras.free_algebra_quotient_element.is_FreeAlgebraQuotientElement(i)
        doctest:warning...
        DeprecationWarning: The function is_FreeAlgebraQuotientElement is deprecated;
        use 'isinstance(..., FreeAlgebraQuotientElement)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        True

    Of course this is testing the data type::

        sage: sage.algebras.free_algebra_quotient_element.is_FreeAlgebraQuotientElement(1)
        False
        sage: sage.algebras.free_algebra_quotient_element.is_FreeAlgebraQuotientElement(H(1))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_FreeAlgebraQuotientElement is deprecated; "
                "use 'isinstance(..., FreeAlgebraQuotientElement)' instead.")
    return isinstance(x, FreeAlgebraQuotientElement)


class FreeAlgebraQuotientElement(AlgebraElement):
    def __init__(self, A, x):
        """
        Create the element x of the FreeAlgebraQuotient A.

        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(ZZ)
            sage: sage.algebras.free_algebra_quotient.FreeAlgebraQuotientElement(H, i)
            i
            sage: a = sage.algebras.free_algebra_quotient.FreeAlgebraQuotientElement(H, 1); a
            1
            sage: a in H
            True

        TESTS::

            sage: TestSuite(i).run()
        """
        AlgebraElement.__init__(self, A)
        Q = self.parent()

        if isinstance(x, FreeAlgebraQuotientElement) and x.parent() == Q:
            self.__vector = Q.module()(x.vector())
            return
        if isinstance(x, (Integer, int)):
            self.__vector = Q.module().gen(0) * x
            return
        elif isinstance(x, FreeModuleElement) and x.parent() is Q.module():
            self.__vector = x
            return
        elif isinstance(x, FreeModuleElement) and x.parent() == A.module():
            self.__vector = x
            return
        R = A.base_ring()
        M = A.module()
        F = A.monoid()
        B = A.monomial_basis()

        if isinstance(x, (Integer, int)):
            self.__vector = x*M.gen(0)
        elif isinstance(x, RingElement) and not isinstance(x, AlgebraElement) and x in R:
            self.__vector = x * M.gen(0)
        elif isinstance(x, FreeMonoidElement) and x.parent() is F:
            if x in B:
                self.__vector = M.gen(B.index(x))
            else:
                raise AttributeError("argument x (= %s) is not in monomial basis" % x)
        elif isinstance(x, list) and len(x) == A.dimension():
            try:
                self.__vector = M(x)
            except TypeError:
                raise TypeError("argument x (= %s) is of the wrong type" % x)
        elif isinstance(x, FreeAlgebraElement) and x.parent() is A.free_algebra():
            # Need to do more work here to include monomials not
            # represented in the monomial basis.
            self.__vector = M(0)
            for m, c in x._FreeAlgebraElement__monomial_coefficients.items():
                self.__vector += c*M.gen(B.index(m))
        elif isinstance(x, dict):
            self.__vector = M(0)
            for m, c in x.items():
                self.__vector += c*M.gen(B.index(m))
        elif isinstance(x, AlgebraElement) and x.parent().ambient_algebra() is A:
            self.__vector = x.ambient_algebra_element().vector()
        else:
            raise TypeError("argument x (= %s) is of the wrong type" % x)

    def _repr_(self):
        """
        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(ZZ)
            sage: i._repr_()
            'i'
        """
        Q = self.parent()
        M = Q.monoid()
        with localvars(M, Q.variable_names()):
            cffs = list(self.__vector)
            mons = Q.monomial_basis()
            return repr_lincomb(zip(mons, cffs), strip_one=True)

    def _latex_(self):
        """
        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: ((2/3)*i - j)._latex_()
            '\\frac{2}{3} i - j'
        """
        Q = self.parent()
        M = Q.monoid()
        with localvars(M, Q.variable_names()):
            cffs = tuple(self.__vector)
            mons = Q.monomial_basis()
            return repr_lincomb(zip(mons, cffs), is_latex=True, strip_one=True)

    def vector(self):
        """
        Return underlying vector representation of this element.

        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: ((2/3)*i - j).vector()
            (0, 2/3, -1, 0)
        """
        return self.__vector

    def _richcmp_(self, right, op):
        """
        Compare two quotient algebra elements; done by comparing the
        underlying vector representatives.

        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: i > j
            True
            sage: i == i
            True
            sage: i == 1
            False
            sage: i + j == j + i
            True
        """
        return richcmp(self.vector(), right.vector(), op)

    def __neg__(self):
        """
        Return negative of ``self``.

        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: -i
            -i
            sage: -(2/3*i - 3/7*j + k)
            -2/3*i + 3/7*j - k
        """
        y = self.parent()(0)
        y.__vector = -self.__vector
        return y

    def _add_(self, y):
        """
        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: 2/3*i + 4*j + k
            2/3*i + 4*j + k
        """
        A = self.parent()
        z = A(0)
        z.__vector = self.__vector + y.__vector
        return z

    def _sub_(self, y):
        """
        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: 2/3*i - 4*j
            2/3*i - 4*j
            sage: a = 2/3*i - 4*j; a
            2/3*i - 4*j
            sage: a - a
            0
        """
        A = self.parent()
        z = A(0)
        z.__vector = self.__vector - y.__vector
        return z

    def _mul_(self, y):
        """
        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: a = (5 + 2*i - 3/5*j + 17*k); a*(a+10)
            -5459/25 + 40*i - 12*j + 340*k

        Double check that the above is actually right::

            sage: R.<i,j,k> = QuaternionAlgebra(QQ,-1,-1)
            sage: a = (5 + 2*i - 3/5*j + 17*k); a*(a+10)
            -5459/25 + 40*i - 12*j + 340*k
        """
        A = self.parent()

        def monomial_product(X, w, m):
            mats = X._FreeAlgebraQuotient__matrix_action
            for j, k in m._element_list:
                M = mats[int(j)]
                for _ in range(k):
                    w *= M
            return w
        u = self.__vector.__copy__()
        v = y.__vector
        z = A(0)
        B = A.monomial_basis()
        for i in range(A.dimension()):
            c = v[i]
            if c != 0:
                z.__vector += monomial_product(A,c*u,B[i])
        return z

    def _rmul_(self, c):
        """
        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: 3 * (-1+i-2*j+k)
            -3 + 3*i - 6*j + 3*k
            sage: (-1+i-2*j+k)._rmul_(3)
            -3 + 3*i - 6*j + 3*k
        """
        return self.parent([c*a for a in self.__vector])

    def _lmul_(self, c):
        """
        EXAMPLES::

            sage: H, (i,j,k) = sage.algebras.free_algebra_quotient.hamilton_quatalg(QQ)
            sage: (-1+i-2*j+k) * 3
            -3 + 3*i - 6*j + 3*k
            sage: (-1+i-2*j+k)._lmul_(3)
            -3 + 3*i - 6*j + 3*k
        """
        return self.parent([a*c for a in self.__vector])
