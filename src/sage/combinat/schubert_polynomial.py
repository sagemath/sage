r"""
Schubert Polynomials
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from __future__ import absolute_import

from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra
from sage.categories.all import GradedAlgebrasWithBasis
from sage.rings.all import Integer, PolynomialRing
from sage.rings.polynomial.multi_polynomial import is_MPolynomial
from . import permutation
import sage.libs.symmetrica.all as symmetrica

from sage.combinat.permutation import Permutations

def SchubertPolynomialRing(R):
    """
    Returns the Schubert polynomial ring over R on the X basis.

    EXAMPLES::

        sage: X = SchubertPolynomialRing(ZZ); X
        Schubert polynomial ring with X basis over Integer Ring
        sage: X(1)
        X[1]
        sage: X([1,2,3])*X([2,1,3])
        X[2, 1]
        sage: X([2,1,3])*X([2,1,3])
        X[3, 1, 2]
        sage: X([2,1,3])+X([3,1,2,4])
        X[2, 1] + X[3, 1, 2]
        sage: a = X([2,1,3])+X([3,1,2,4])
        sage: a^2
        X[3, 1, 2] + 2*X[4, 1, 2, 3] + X[5, 1, 2, 3, 4]
    """
    return SchubertPolynomialRing_xbasis(R)

def is_SchubertPolynomial(x):
    """
    Returns True if x is a Schubert polynomial and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.schubert_polynomial import is_SchubertPolynomial
        sage: X = SchubertPolynomialRing(ZZ)
        sage: a = 1
        sage: is_SchubertPolynomial(a)
        False
        sage: b = X(1)
        sage: is_SchubertPolynomial(b)
        True
        sage: c = X([2,1,3])
        sage: is_SchubertPolynomial(c)
        True
    """
    return isinstance(x, SchubertPolynomial_class)

class SchubertPolynomial_class(CombinatorialFreeModule.Element):
    def expand(self):
        """
        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: X([2,1,3]).expand()
            x0
            sage: map(lambda x: x.expand(), [X(p) for p in Permutations(3)])
            [1, x0 + x1, x0, x0*x1, x0^2, x0^2*x1]

        TESTS:

        Calling .expand() should always return an element of an
        MPolynomialRing

        ::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: f = X([1]); f
            X[1]
            sage: type(f.expand())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f.expand()
            1
            sage: f = X([1,2])
            sage: type(f.expand())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f = X([1,3,2,4])
            sage: type(f.expand())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
        """
        p = symmetrica.t_SCHUBERT_POLYNOM(self)
        if not is_MPolynomial(p):
            R = PolynomialRing(self.parent().base_ring(), 1, 'x')
            p = R(p)
        return p

    def divided_difference(self, i):
        """
        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: a = X([3,2,1])
            sage: a.divided_difference(1)
            X[2, 3, 1]
            sage: a.divided_difference([3,2,1])
            X[1]
        """
        if isinstance(i, Integer):
            return symmetrica.divdiff_schubert(i, self)
        elif i in permutation.Permutations():
            return symmetrica.divdiff_perm_schubert(i, self)
        else:
            raise TypeError("i must either be an integer or permutation")

    def scalar_product(self, x):
        """
        Returns the standard scalar product of self and x.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: a = X([3,2,4,1])
            sage: a.scalar_product(a)
            0
            sage: b = X([4,3,2,1])
            sage: b.scalar_product(a)
            X[1, 3, 4, 6, 2, 5]
            sage: Permutation([1, 3, 4, 6, 2, 5, 7]).to_lehmer_code()
            [0, 1, 1, 2, 0, 0, 0]
            sage: s = SymmetricFunctions(ZZ).schur()
            sage: c = s([2,1,1])
            sage: b.scalar_product(a).expand()
            x0^2*x1*x2 + x0*x1^2*x2 + x0*x1*x2^2 + x0^2*x1*x3 + x0*x1^2*x3 + x0^2*x2*x3 + 3*x0*x1*x2*x3 + x1^2*x2*x3 + x0*x2^2*x3 + x1*x2^2*x3 + x0*x1*x3^2 + x0*x2*x3^2 + x1*x2*x3^2
            sage: c.expand(4)
            x0^2*x1*x2 + x0*x1^2*x2 + x0*x1*x2^2 + x0^2*x1*x3 + x0*x1^2*x3 + x0^2*x2*x3 + 3*x0*x1*x2*x3 + x1^2*x2*x3 + x0*x2^2*x3 + x1*x2^2*x3 + x0*x1*x3^2 + x0*x2*x3^2 + x1*x2*x3^2
        """
        if is_SchubertPolynomial(x):
            return symmetrica.scalarproduct_schubert(self, x)
        else:
            raise TypeError("x must be a Schubert polynomial")

    def multiply_variable(self, i):
        """
        Returns the Schubert polynomial obtained by multiplying self by the
        variable `x_i`.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(ZZ)
            sage: a = X([3,2,4,1])
            sage: a.multiply_variable(0)
            X[4, 2, 3, 1]
            sage: a.multiply_variable(1)
            X[3, 4, 2, 1]
            sage: a.multiply_variable(2)
            X[3, 2, 5, 1, 4] - X[3, 4, 2, 1] - X[4, 2, 3, 1]
            sage: a.multiply_variable(3)
            X[3, 2, 4, 5, 1]
        """
        if isinstance(i, Integer):
            return symmetrica.mult_schubert_variable(self, i)
        else:
            raise TypeError("i must be an integer")

# FIXME: inherit from CombinatorialFreeModule once the
# coercion from ground ring is implemented in the category
class SchubertPolynomialRing_xbasis(CombinatorialAlgebra):

    Element = SchubertPolynomial_class

    def __init__(self, R):
        """
        EXAMPLES::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X == loads(dumps(X))
            True
        """
        self._name = "Schubert polynomial ring with X basis"
        self._repr_option_bracket = False
        self._one = permutation.Permutation([1])
        CombinatorialAlgebra.__init__(self, R, cc = permutation.Permutations(), category = GradedAlgebrasWithBasis(R))
        self.print_options(prefix='X')

    def _element_constructor_(self, x):
        """
        Coerce x into self.

        EXAMPLES::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X._element_constructor_([2,1,3])
            X[2, 1]
            sage: X._element_constructor_(Permutation([2,1,3]))
            X[2, 1]

            sage: R.<x1, x2, x3> = QQ[]
            sage: X(x1^2*x2)
            X[3, 2, 1]

        TESTS:

        We check that :trac:`12924` is fixed::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X._element_constructor_([1,2,1])
            Traceback (most recent call last):
            ...
            ValueError: The input [1, 2, 1] is not a valid permutation
        """
        if isinstance(x, list):
            #checking the input to avoid symmetrica crashing Sage, see trac 12924
            if not x in Permutations():
                raise ValueError("The input %s is not a valid permutation"%(x))
            perm = permutation.Permutation(x).remove_extra_fixed_points()
            return self._from_dict({ perm: self.base_ring()(1) })
        elif isinstance(x, permutation.Permutation):
            if not list(x) in Permutations():
                raise ValueError("The input %s is not a valid permutation"%(x))
            perm = x.remove_extra_fixed_points()
            return self._from_dict({ perm: self.base_ring()(1) })
        elif is_MPolynomial(x):
            return symmetrica.t_POLYNOM_SCHUBERT(x)
        else:
            raise TypeError

    def _multiply_basis(self, left, right):
        """
        EXAMPLES::

            sage: p1 = Permutation([3,2,1])
            sage: p2 = Permutation([2,1,3])
            sage: X = SchubertPolynomialRing(QQ)
            sage: X._multiply_basis(p1,p2)
            {[4, 2, 1, 3]: 1}
        """
        return symmetrica.mult_schubert_schubert(left, right).monomial_coefficients()
