# sage.doctest: needs sage.combinat sage.modules
"""
Free algebra elements

AUTHORS:

- David Kohel (2005-09)

TESTS::

    sage: R.<x,y> = FreeAlgebra(QQ,2)
    sage: x == loads(dumps(x))
    True
    sage: x*y
    x*y
    sage: (x*y)^0
    1
    sage: (x*y)^3
    x*y*x*y*x*y
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
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement
from sage.structure.element import AlgebraElement


class FreeAlgebraElement(IndexedFreeModuleElement, AlgebraElement):
    """
    A free algebra element.

    TESTS:

    The ordering is inherited from ``IndexedFreeModuleElement``::

        sage: R.<x,y> = FreeAlgebra(QQ,2)
        sage: x < y
        True
        sage: x * y < y * x
        True
        sage: y * x < x * y
        False
    """
    def __init__(self, A, x):
        """
        Create the element ``x`` of the FreeAlgebra ``A``.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: elt = x^3 * y - z^2*x
            sage: TestSuite(elt).run()
        """
        if isinstance(x, FreeAlgebraElement):
            # We should have an input for when we know we don't need to
            #   convert the keys/values
            x = x._monomial_coefficients
        R = A.base_ring()
        if isinstance(x, AlgebraElement):  # and x.parent() == A.base_ring():
            x = {A.monoid()(1): R(x)}
        elif isinstance(x, FreeMonoidElement):
            x = {x: R(1)}
        elif True:
            x = {A.monoid()(e1): R(e2) for e1, e2 in x.items()}
        else:
            raise TypeError("argument x (= {}) is of the wrong type".format(x))

        IndexedFreeModuleElement.__init__(self, A, x)

    def _repr_(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(ZZ,3)
            sage: repr(-x+3*y*z)    # indirect doctest
            '-x + 3*y*z'

        Github issue :issue:`11068` enables the use of local variable names::

            sage: from sage.structure.parent_gens import localvars
            sage: with localvars(A, ['a','b','c']):
            ....:    print(-x+3*y*z)
            -a + 3*b*c
        """
        v = sorted(self._monomial_coefficients.items())
        P = self.parent()
        M = P.monoid()
        from sage.structure.parent_gens import localvars
        with localvars(M, P.variable_names(), normalize=False):
            x = repr_lincomb(v, strip_one=True)
        return x

    def _latex_(self):
        r"""
        Return latex representation of ``self``.

        EXAMPLES::

            sage: A.<x,y,z>=FreeAlgebra(ZZ,3)
            sage: latex(-x+3*y^20*z)   # indirect doctest
            -x + 3 y^{20}z
            sage: alpha,beta,gamma=FreeAlgebra(ZZ,3,'alpha,beta,gamma').gens()
            sage: latex(alpha-beta)
            \alpha - \beta
        """
        v = sorted(self._monomial_coefficients.items())
        return repr_lincomb(v, strip_one=True, is_latex=True)

    def __call__(self, *x, **kwds):
        """
        EXAMPLES::

            sage: A.<x,y,z>=FreeAlgebra(ZZ,3)
            sage: (x+3*y).subs(x=1,y=2,z=14)
            7
            sage: (2*x+y).subs({x:1,y:z})
            2 + z
            sage: f=x+3*y+z
            sage: f(1,2,1/2)
            15/2
            sage: f(1,2)
            Traceback (most recent call last):
            ...
            ValueError: must specify as many values as generators in parent

        AUTHORS:

        - Joel B. Mohler (2007-10-27)
        """
        if kwds and x:
            raise ValueError("must not specify both a keyword and positional argument")

        if kwds:
            p = self.parent()

            def extract_from(kwds, g):
                for x in g:
                    try:
                        return kwds[x]
                    except KeyError:
                        pass
                return None

            x = [extract_from(kwds,(p.gen(i),p.variable_name(i))) for i in range(p.ngens())]
        elif isinstance(x[0], tuple):
            x = x[0]

        if len(x) != self.parent().ngens():
            raise ValueError("must specify as many values as generators in parent")

        # I don't start with 0, because I don't want to preclude evaluation with
        # arbitrary objects (e.g. matrices) because of funny coercion.
        result = None
        for m, c in self._monomial_coefficients.items():
            if result is None:
                result = c*m(x)
            else:
                result += c*m(x)

        if result is None:
            return self.parent().zero()
        return result

    def _mul_(self, y):
        """
        Return the product of ``self`` and ``y`` (another free algebra
        element with the same parent).

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(ZZ,3)
            sage: (x+y+x*y)*(x+y+1)
            x + y + x^2 + 2*x*y + y*x + y^2 + x*y*x + x*y^2
        """
        A = self.parent()
        z_elt = {}
        for mx, cx in self:
            for my, cy in y:
                key = mx*my
                if key in z_elt:
                    z_elt[key] += cx*cy
                else:
                    z_elt[key] = cx*cy
                if not z_elt[key]:
                    del z_elt[key]
        return A._from_dict(z_elt)

    def is_unit(self):
        r"""
        Return ``True`` if ``self`` is invertible.

        EXAMPLES::

            sage: A.<x, y, z> = FreeAlgebra(ZZ)
            sage: A(-1).is_unit()
            True
            sage: A(2).is_unit()
            False
            sage: A(1 + x).is_unit()
            False
            sage: A.<x, y> = FreeAlgebra(QQ, degrees=(1,-1))
            sage: A(x * y).is_unit()
            False
            sage: A(2).is_unit()
            True
        """
        mc = self._monomial_coefficients
        if not mc or len(mc) > 1:
            return False
        m, c = next(iter(mc.items()))
        return m.is_one() and c.is_unit()

    def __invert__(self):
        """
        EXAMPLES::

            sage: A.<x, y, z> = FreeAlgebra(QQ)
            sage: ~A(1)
            1

        TESTS::

            sage: ~A(0)
            Traceback (most recent call last):
            ...
            ArithmeticError: element is not invertible

            sage: ~A(1 + x)
            Traceback (most recent call last):
            ...
            ArithmeticError: element is not invertible
        """
        if self.is_unit():
            m, c = next(iter(self._monomial_coefficients.items()))
            return type(self)(self.parent(), {m: c.inverse_of_unit()})
        raise ArithmeticError("element is not invertible")

    def _acted_upon_(self, scalar, self_on_left=False):
        """
        Return the action of a scalar on ``self``.

        EXAMPLES::

            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: f = Factorization([(x,2),(y,3)]); f
            x^2 * y^3
            sage: x * f
            x^3 * y^3
            sage: f * x
            x^2 * y^3 * x
        """
        from sage.structure.factorization import Factorization
        # FIXME: Make factorization work properly in the coercion framework
        # Keep factorization since we want to "coerce" into a factorization
        if isinstance(scalar, Factorization):
            if self_on_left:
                return Factorization([(self, 1)]) * scalar
            return scalar * Factorization([(self, 1)])
        return super()._acted_upon_(scalar, self_on_left)

    # For backward compatibility
    # _lmul_ = _acted_upon_
    # _rmul_ = _acted_upon_

    def _im_gens_(self, codomain, im_gens, base_map):
        """
        Apply a morphism defined by its values on the generators.

        EXAMPLES::

            sage: ring = algebras.Free(QQ, ['a', 'b'])
            sage: a, b = ring.gens()
            sage: A = matrix(QQ, 2, 2, [2, 3, 4, 1])
            sage: B = matrix(QQ, 2, 2, [1, 7, 7, 1])
            sage: f = ring.hom([A, B])
            sage: f(a*b+1)
            [24 17]
            [11 30]
        """
        n = self.parent().ngens()
        if n == 0:
            cf = next(iter(self._monomial_coefficients.values()))
            return codomain.coerce(cf)

        if base_map is None:
            base_map = codomain

        return codomain.sum(base_map(c) * m(*im_gens)
                            for m, c in self._monomial_coefficients.items())

    def variables(self):
        """
        Return the variables used in ``self``.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(ZZ,3)
            sage: elt = x + x*y + x^3*y
            sage: elt.variables()
            [x, y]
            sage: elt = x + x^2 - x^4
            sage: elt.variables()
            [x]
            sage: elt = x + z*y + z*x
            sage: elt.variables()
            [x, y, z]
        """
        v = set()
        for s in self._monomial_coefficients:  # Only gets the keys
            for var, _ in s:
                v.add(var)
        A = self.parent()
        return sorted(map(A, v))

    def to_pbw_basis(self):
        """
        Return ``self`` in the Poincaré-Birkhoff-Witt (PBW) basis.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(ZZ, 3)
            sage: p = x^2*y + 3*y*x + 2
            sage: p.to_pbw_basis()
            2*PBW[1] + 3*PBW[y]*PBW[x] + PBW[x^2*y]
             + 2*PBW[x*y]*PBW[x] + PBW[y]*PBW[x]^2
        """
        return self.parent().pbw_element(self)
