# sage.doctest: needs sage.combinat sage.modules
r"""
Weyl Algebras

AUTHORS:

- Travis Scrimshaw (2013-09-06): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex, LatexExpr
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.action import Action
from sage.categories.rings import Rings
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.sets.family import Family
import sage.data_structures.blas_dict as blas
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.global_options import GlobalOptions
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement


def repr_from_monomials(monomials, term_repr, use_latex=False) -> str:
    r"""
    Return a string representation of an element of a free module
    from the dictionary ``monomials``.

    INPUT:

    - ``monomials`` -- list of pairs ``[m, c]`` where ``m`` is the index
      and ``c`` is the coefficient
    - ``term_repr`` -- a function which returns a string given an index
      (can be ``repr`` or ``latex``, for example)
    - ``use_latex`` -- boolean (default: ``False``); if ``True`` then the
      output is in latex format

    EXAMPLES::

        sage: from sage.algebras.weyl_algebra import repr_from_monomials
        sage: R.<x,y,z> = QQ[]
        sage: d = [(z, 4/7), (y, sqrt(2)), (x, -5)]                                     # needs sage.symbolic
        sage: repr_from_monomials(d, lambda m: repr(m))                                 # needs sage.symbolic
        '4/7*z + sqrt(2)*y - 5*x'
        sage: a = repr_from_monomials(d, lambda m: latex(m), True); a                   # needs sage.symbolic
        \frac{4}{7} z + \sqrt{2} y - 5 x
        sage: type(a)                                                                   # needs sage.symbolic
        <class 'sage.misc.latex.LatexExpr'>

    The zero element::

        sage: repr_from_monomials([], lambda m: repr(m))
        '0'
        sage: a = repr_from_monomials([], lambda m: latex(m), True); a
        0
        sage: type(a)
        <class 'sage.misc.latex.LatexExpr'>

    A "unity" element::

        sage: repr_from_monomials([(1, 1)], lambda m: repr(m))
        '1'
        sage: a = repr_from_monomials([(1, 1)], lambda m: latex(m), True); a
        1
        sage: type(a)
        <class 'sage.misc.latex.LatexExpr'>

    ::

        sage: repr_from_monomials([(1, -1)], lambda m: repr(m))
        '-1'
        sage: a = repr_from_monomials([(1, -1)], lambda m: latex(m), True); a
        -1
        sage: type(a)
        <class 'sage.misc.latex.LatexExpr'>

    Leading minus signs are dealt with appropriately::

        sage: # needs sage.symbolic
        sage: d = [(z, -4/7), (y, -sqrt(2)), (x, -5)]
        sage: repr_from_monomials(d, lambda m: repr(m))
        '-4/7*z - sqrt(2)*y - 5*x'
        sage: a = repr_from_monomials(d, lambda m: latex(m), True); a
        -\frac{4}{7} z - \sqrt{2} y - 5 x
        sage: type(a)
        <class 'sage.misc.latex.LatexExpr'>

    Indirect doctests using a class that uses this function::

        sage: R.<x,y> = QQ[]
        sage: A = CliffordAlgebra(QuadraticForm(R, 3, [x,0,-1,3,-4,5]))
        sage: a,b,c = A.gens()
        sage: a*b*c
        e0*e1*e2
        sage: b*c
        e1*e2
        sage: (a*a + 2)
        x + 2
        sage: c*(a*a + 2)*b
        (-x - 2)*e1*e2 - 4*x - 8
        sage: latex(c*(a*a + 2)*b)
        \left( -x - 2 \right)  e_{1} e_{2} - 4 x - 8
    """
    if not monomials:
        if use_latex:
            return latex(0)
        else:
            return '0'

    ret = ''
    for m, c in monomials:
        # Get the monomial portion
        term = term_repr(m)

        # Determine what to do with the coefficient
        if use_latex:
            coeff = latex(c)
        else:
            coeff = repr(c)

        if not term or term == '1':
            term = coeff
        elif coeff == '-1':
            term = '-' + term
        elif coeff != '1':
            atomic_repr = c.parent()._repr_option('element_is_atomic')
            if not atomic_repr and (coeff.find("+") != -1 or coeff.rfind("-") > 0):
                if use_latex:
                    term = '\\left(' + coeff + '\\right) ' + term
                elif coeff not in ['', '-']:
                    term = '(' + coeff + ')*' + term
            else:
                if use_latex:
                    term = coeff + ' ' + term
                else:
                    term = coeff + '*' + term

        # Append this term with the correct sign
        if ret:
            if term[0] == '-':
                ret += ' - ' + term[1:]
            else:
                ret += ' + ' + term
        else:
            ret = term
    return ret


def repr_factored(w, latex_output=False) -> str:
    r"""
    Return a string representation of ``w`` with the `dx_i` generators
    factored on the right.

    EXAMPLES::

        sage: from sage.algebras.weyl_algebra import repr_factored
        sage: R.<t> = QQ[]
        sage: D = DifferentialWeylAlgebra(R)
        sage: t, dt = D.gens()
        sage: x = dt^3*t^3 + dt^2*t^4
        sage: x
        t^3*dt^3 + t^4*dt^2 + 9*t^2*dt^2 + 8*t^3*dt + 18*t*dt + 12*t^2 + 6
        sage: print(repr_factored(x))
        (12*t^2 + 6) + (8*t^3 + 18*t)*dt + (t^4 + 9*t^2)*dt^2 + (t^3)*dt^3
        sage: repr_factored(x, True)
        (12 t^{2} + 6) + (8 t^{3} + 18 t) \frac{\partial}{\partial t}
         + (t^{4} + 9 t^{2}) \frac{\partial^{2}}{\partial t^{2}}
         + (t^{3}) \frac{\partial^{3}}{\partial t^{3}}
        sage: repr_factored(D.zero())
        '0'

    With multiple variables::

        sage: R.<x,y,z> = QQ[]
        sage: D = DifferentialWeylAlgebra(R)
        sage: x, y, z, dx, dy, dz = D.gens()
        sage: elt = dx^3*x^3 + (y^3-z*x)*dx^3 + dy^3*x^3 + dx*dy*dz*x*y*z
        sage: elt
        x^3*dy^3 + x*y*z*dx*dy*dz + y^3*dx^3 + x^3*dx^3 - x*z*dx^3 + y*z*dy*dz
         + x*z*dx*dz + x*y*dx*dy + 9*x^2*dx^2 + z*dz + y*dy + 19*x*dx + 7
        sage: print(repr_factored(elt))
        (7) + (z)*dz + (y)*dy + (y*z)*dy*dz + (x^3)*dy^3 + (19*x)*dx
         + (x*z)*dx*dz + (x*y)*dx*dy + (x*y*z)*dx*dy*dz
         + (9*x^2)*dx^2 + (x^3 + y^3 - x*z)*dx^3
        sage: repr_factored(D.zero(), True)
        0
    """
    f = w.factor_differentials()
    gens = w.parent().polynomial_ring().gens()

    if latex_output:
        def exp(e):
            return '^{{{}}}'.format(e) if e > 1 else ''

        def repr_dx(k):
            total = sum(k)
            if total == 0:
                return ''
            denom = ' '.join('\\partial {}{}'.format(latex(g), exp(e))
                             for e, g in zip(k, gens) if e != 0)
            return ''.join(' \\frac{{\\partial{}}}{{{}}}'.format(exp(total), denom))
        repr_x = latex
    else:
        def exp(e):
            return '^{}'.format(e) if e > 1 else ''

        def repr_dx(k):
            return ''.join('*d{}{}'.format(g, exp(e)) for e, g in zip(k, gens) if e != 0)
        repr_x = repr
    ret = " + ".join("({}){}".format(repr_x(f[k]), repr_dx(k))
                     for k in sorted(f))
    if not ret:
        ret = '0'
    if latex_output:
        return LatexExpr(ret)
    return ret


class DifferentialWeylAlgebraElement(IndexedFreeModuleElement):
    """
    An element in a differential Weyl algebra.

    TESTS:

    Some computations in the Weyl algebra, using the implicit
    coercion from the polynomial ring into the Weyl algebra::

        sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
        sage: dx,dy,dz = W.differentials()
        sage: elt = ((x^3-z)*dx + dy)^2
        sage: TestSuite(elt).run()

        sage: R.<x,y,z> =  QQ[]
        sage: W = DifferentialWeylAlgebra(R)
        sage: dx,dy,dz = W.differentials()
        sage: dy*(x^3-y*z)*dx == -z*dx + x^3*dx*dy - y*z*dx*dy
        True
        sage: W.zero() == 0
        True
        sage: W.one() == 1
        True
        sage: x == 1
        False
        sage: x + 1 == 1
        False
        sage: W(x^3 - y*z) == x^3 - y*z
        True
        sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
        sage: dx,dy,dz = W.differentials()
        sage: dx != dy
        True
        sage: W.one() != 1
        False
        sage: dy - (3*x - z)*dx
        dy + z*dx - 3*x*dx
        sage: (dx*dy) + dz + x^3 - 2
        dx*dy + dz + x^3 - 2
        sage: elt = (dy - (3*x - z)*dx)
        sage: sorted(elt.monomial_coefficients().items())
        [(((0, 0, 0), (0, 1, 0)), 1),
         (((0, 0, 1), (1, 0, 0)), 1),
         (((1, 0, 0), (1, 0, 0)), -3)]
        sage: elt = dy - (3*x - z)*dx + 1
        sage: sorted(elt.support())
        [((0, 0, 0), (0, 0, 0)),
         ((0, 0, 0), (0, 1, 0)),
         ((0, 0, 1), (1, 0, 0)),
         ((1, 0, 0), (1, 0, 0))]

    Hashing works::

        sage: hash(dx) == hash(dx) # hashing works
        True

    Comparison works, even if mathematically meaningless
    (but useful e.g. for sorting)::

        sage: dx < dy or dy < dx
        True
    """
    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: dx,dy,dz = W.differentials()
            sage: ((x^3-z)*dx + dy)^2
            dy^2 + 2*x^3*dx*dy - 2*z*dx*dy + x^6*dx^2 - 2*x^3*z*dx^2
             + z^2*dx^2 + 3*x^5*dx - 3*x^2*z*dx
        """
        if self.parent().options.factor_representation:
            return repr_factored(self, False)

        def term(m):
            ret = ''
            for i, power in enumerate(m[0] + m[1]):
                if power == 0:
                    continue
                name = self.parent().variable_names()[i]
                if ret:
                    ret += '*'
                if power == 1:
                    ret += '{}'.format(name)
                else:
                    ret += '{}^{}'.format(name, power)
            return ret
        return repr_from_monomials(self.list(), term)

    def _latex_(self) -> str:
        r"""
        Return a `\LaTeX` representation of ``self``.

        TESTS::

            sage: R = PolynomialRing(QQ, 'x', 3)
            sage: W = DifferentialWeylAlgebra(R)
            sage: x0,x1,x2,dx0,dx1,dx2 = W.gens()
            sage: latex( ((x0^3-x2)*dx0 + dx1)^2 )
            \frac{\partial^{2}}{\partial x_{1}^{2}}
             + 2 x_{0}^{3} \frac{\partial^{2}}{\partial x_{0} \partial x_{1}}
             - 2 x_{2} \frac{\partial^{2}}{\partial x_{0} \partial x_{1}}
             + x_{0}^{6} \frac{\partial^{2}}{\partial x_{0}^{2}}
             - 2 x_{0}^{3} x_{2} \frac{\partial^{2}}{\partial x_{0}^{2}}
             + x_{2}^{2} \frac{\partial^{2}}{\partial x_{0}^{2}}
             + 3 x_{0}^{5} \frac{\partial}{\partial x_{0}}
             - 3 x_{0}^{2} x_{2} \frac{\partial}{\partial x_{0}}
        """
        if self.parent().options.factor_representation:
            return repr_factored(self, True)

        def exp(e):
            return '^{{{}}}'.format(e) if e > 1 else ''

        def term(m):
            R = self.parent()._poly_ring

            def half_term(mon, polynomial):
                total = sum(mon)
                if total == 0:
                    return '1'
                ret = ' '.join('{}{}'.format(latex(R.gen(i)), exp(power)) if polynomial
                               else '\\partial {}{}'.format(latex(R.gen(i)), exp(power))
                               for i, power in enumerate(mon) if power > 0)
                if not polynomial:
                    return '\\frac{{\\partial{}}}{{{}}}'.format(exp(total), ret)
                return ret
            p = half_term(m[0], True)
            d = half_term(m[1], False)
            if p == '1':  # No polynomial part
                return d
            elif d == '1':  # No differential part
                return p
            else:
                return p + ' ' + d
        return repr_from_monomials(self.list(), term, True)

    def _mul_(self, other):
        """
        Return ``self`` multiplied by ``other``.

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: dx,dy,dz = W.differentials()
            sage: dx*(x*y + z)
            x*y*dx + z*dx + y
            sage: ((x^3-z)*dx + dy) * (dx*dz^2 - 10*x)
            dx*dy*dz^2 + x^3*dx^2*dz^2 - z*dx^2*dz^2 - 10*x*dy - 10*x^4*dx
             + 10*x*z*dx - 10*x^3 + 10*z
        """
        def add_tuples(x, y):
            return tuple(a + y[i] for i, a in enumerate(x))

        d = {}
        n = self.parent()._n
        t = tuple([0] * n)
        zero = self.parent().base_ring().zero()
        for ml in self._monomial_coefficients:
            cl = self._monomial_coefficients[ml]
            for mr in other._monomial_coefficients:
                cr = other._monomial_coefficients[mr]
                cur = [((mr[0], t), cl * cr)]
                for i, p in enumerate(ml[1]):
                    for _ in range(p):
                        next = []
                        for m, c in cur:  # Distribute and apply the derivative
                            diff = list(m[1])
                            diff[i] += 1
                            next.append(((m[0], tuple(diff)), c))
                            if m[0][i] != 0:
                                poly = list(m[0])
                                c *= poly[i]
                                poly[i] -= 1
                                next.append(((tuple(poly), m[1]), c))
                        cur = next

                for m, c in cur:
                    # multiply the resulting term by the other term
                    m = (add_tuples(ml[0], m[0]), add_tuples(mr[1], m[1]))
                    d[m] = d.get(m, zero) + c
                    if d[m] == zero:
                        del d[m]
        return self.__class__(self.parent(), d)

    def _rmul_(self, other):
        """
        Multiply ``self`` on the right side of ``other``.

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: dx,dy,dz = W.differentials()
            sage: a = (x*y + z) * dx
            sage: 3/2 * a
            3/2*x*y*dx + 3/2*z*dx
        """
        if other == 0:
            return self.parent().zero()
        M = self._monomial_coefficients
        return self.__class__(self.parent(), {t: other * M[t] for t in M})

    def _lmul_(self, other):
        """
        Multiply ``self`` on the left side of ``other``.

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: dx,dy,dz = W.differentials()
            sage: a = (x*y + z) * dx
            sage: a * 3/2
            3/2*x*y*dx + 3/2*z*dx
        """
        if other == 0:
            return self.parent().zero()
        M = self._monomial_coefficients
        return self.__class__(self.parent(), {t: M[t] * other for t in M})

    def __iter__(self):
        """
        Return an iterator of ``self``.

        This is the iterator of ``self.list()``.

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: dx,dy,dz = W.differentials()
            sage: list(dy - (3*x - z)*dx)
            [(((0, 0, 0), (0, 1, 0)), 1),
             (((0, 0, 1), (1, 0, 0)), 1),
             (((1, 0, 0), (1, 0, 0)), -3)]
        """
        return iter(self.list())

    def list(self) -> list:
        """
        Return ``self`` as a list.

        This list consists of pairs `(m, c)`, where `m` is a pair of
        tuples indexing a basis element of ``self``, and `c` is the
        coordinate of ``self`` corresponding to this basis element.
        (Only nonzero coordinates are shown.)

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: dx,dy,dz = W.differentials()
            sage: elt = dy - (3*x - z)*dx
            sage: elt.list()
            [(((0, 0, 0), (0, 1, 0)), 1),
             (((0, 0, 1), (1, 0, 0)), 1),
             (((1, 0, 0), (1, 0, 0)), -3)]
        """
        return sorted(self._monomial_coefficients.items(),
                      key=lambda x: (-sum(x[0][1]), x[0][1], -sum(x[0][0]), x[0][0]))

    # This is essentially copied from
    #   sage.combinat.free_module.CombinatorialFreeModuleElement
    def __truediv__(self, x):
        """
        Division by coefficients.

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: x / 2
            1/2*x
            sage: W.<x,y,z> = DifferentialWeylAlgebra(ZZ)
            sage: a = 2*x + 4*y*z
            sage: a / 2
            2*y*z + x
        """
        F = self.parent()
        D = self._monomial_coefficients
        if F.base_ring().is_field():
            x = F.base_ring()(x)
            x_inv = x**-1
            D = blas.linear_combination([(D, x_inv)])

            return self.__class__(F, D)

        return self.__class__(F, {t: D[t]._divide_if_possible(x) for t in D})

    def factor_differentials(self) -> dict:
        """
        Return a dict representing ``self`` with the differentials
        factored out.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: D = DifferentialWeylAlgebra(R)
            sage: t, dt = D.gens()
            sage: x = dt^3*t^3 + dt^2*t^4
            sage: x
            t^3*dt^3 + t^4*dt^2 + 9*t^2*dt^2 + 8*t^3*dt + 18*t*dt + 12*t^2 + 6
            sage: x.factor_differentials()
            {(0,): 12*t^2 + 6, (1,): 8*t^3 + 18*t, (2,): t^4 + 9*t^2, (3,): t^3}
            sage: D.zero().factor_differentials()
            {}

            sage: R.<x,y,z> = QQ[]
            sage: D = DifferentialWeylAlgebra(R)
            sage: x, y, z, dx, dy, dz = D.gens()
            sage: elt = dx^3*x^3 + (y^3-z*x)*dx^3 + dy^3*x^3 + dx*dy*dz*x*y*z
            sage: elt
            x^3*dy^3 + x*y*z*dx*dy*dz + y^3*dx^3 + x^3*dx^3 - x*z*dx^3 + y*z*dy*dz
             + x*z*dx*dz + x*y*dx*dy + 9*x^2*dx^2 + z*dz + y*dy + 19*x*dx + 7
            sage: elt.factor_differentials()
            {(0, 0, 0): 7,
             (0, 0, 1): z,
             (0, 1, 0): y,
             (0, 1, 1): y*z,
             (0, 3, 0): x^3,
             (1, 0, 0): 19*x,
             (1, 0, 1): x*z,
             (1, 1, 0): x*y,
             (1, 1, 1): x*y*z,
             (2, 0, 0): 9*x^2,
             (3, 0, 0): x^3 + y^3 - x*z}
        """
        ret = {}
        DW = self.parent()
        P = DW.polynomial_ring()
        gens = P.gens()
        for m, c in self:
            x, dx = m
            if dx not in ret:
                ret[dx] = P.zero()
            ret[dx] += c * prod(g**e for e, g in zip(x, gens))
        return ret

    def diff(self, p):
        """
        Apply this differential operator to a polynomial.

        INPUT:

        - ``p`` -- polynomial of the underlying polynomial ring

        OUTPUT:

        The result of the left action of the Weyl algebra on the polynomial
        ring via differentiation.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: W = R.weyl_algebra()
            sage: dx, dy = W.differentials()
            sage: dx.diff(x^3)
            3*x^2
            sage: (dx*dy).diff(W(x^3*y^3))
            9*x^2*y^2
            sage: (x*dx + dy + 1).diff(x^4*y^4 + 1)
            5*x^4*y^4 + 4*x^4*y^3 + 1
        """
        return self.parent().diff_action(self, p)


class DifferentialWeylAlgebra(UniqueRepresentation, Parent):
    r"""
    The differential Weyl algebra of a polynomial ring.

    Let `R` be a commutative ring. The (differential) Weyl algebra `W` is
    the algebra generated by `x_1, x_2, \ldots x_n, \partial_{x_1},
    \partial_{x_2}, \ldots, \partial_{x_n}` subject to the relations:
    `[x_i, x_j] = 0`, `[\partial_{x_i}, \partial_{x_j}] = 0`, and
    `\partial_{x_i} x_j = x_j \partial_{x_i} + \delta_{ij}`. Therefore
    `\partial_{x_i}` is acting as the partial differential operator on `x_i`.

    The Weyl algebra can also be constructed as an iterated Ore extension
    of the polynomial ring `R[x_1, x_2, \ldots, x_n]` by adding `x_i` at
    each step. It can also be seen as a quantization of the symmetric algebra
    `Sym(V)`, where `V` is a finite dimensional vector space over a field
    of characteristic zero, by using a modified Groenewold-Moyal
    product in the symmetric algebra.

    The Weyl algebra (even for `n = 1`) over a field of characteristic 0
    has many interesting properties.

    - It's a non-commutative domain.
    - It's a simple ring (but not in positive characteristic) that is not
      a matrix ring over a division ring.
    - It has no finite-dimensional representations.
    - It's a quotient of the universal enveloping algebra of the
      Heisenberg algebra `\mathfrak{h}_n`.

    REFERENCES:

    - :wikipedia:`Weyl_algebra`

    INPUT:

    - ``R`` -- a (polynomial) ring
    - ``names`` -- (default: ``None``) if ``None`` and ``R`` is a
      polynomial ring, then the variable names correspond to
      those of ``R``; otherwise if ``names`` is specified, then ``R``
      is the base ring

    EXAMPLES:

    There are two ways to create a Weyl algebra, the first is from
    a polynomial ring::

        sage: R.<x,y,z> = QQ[]
        sage: W = DifferentialWeylAlgebra(R); W
        Differential Weyl algebra of polynomials in x, y, z over Rational Field

    We can call ``W.inject_variables()`` to give the polynomial ring
    variables, now as elements of ``W``, and the differentials::

        sage: W.inject_variables()
        Defining x, y, z, dx, dy, dz
        sage: (dx * dy * dz) * (x^2 * y * z + x * z * dy + 1)
        x*z*dx*dy^2*dz + z*dy^2*dz + x^2*y*z*dx*dy*dz + dx*dy*dz
         + x*dx*dy^2 + 2*x*y*z*dy*dz + dy^2 + x^2*z*dx*dz + x^2*y*dx*dy
         + 2*x*z*dz + 2*x*y*dy + x^2*dx + 2*x

    Or directly by specifying a base ring and variable names::

        sage: W.<a,b> = DifferentialWeylAlgebra(QQ); W
        Differential Weyl algebra of polynomials in a, b over Rational Field

    .. TODO::

        Implement the :meth:`graded_algebra` as a polynomial ring once
        they are considered to be graded rings (algebras).
    """
    @staticmethod
    def __classcall__(cls, R, names=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: W1.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: W2 = DifferentialWeylAlgebra(QQ['x,y,z'])
            sage: W1 is W2
            True
        """
        if isinstance(R, (PolynomialRing_general, MPolynomialRing_base)):
            if names is None:
                names = R.variable_names()
                R = R.base_ring()
        elif names is None:
            raise ValueError("the names must be specified")
        elif R not in Rings().Commutative():
            raise TypeError("argument R must be a commutative ring")
        return super().__classcall__(cls, R, names)

    def __init__(self, R, names=None) -> None:
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: TestSuite(W).run()
        """
        self._n = len(names)
        self._poly_ring = PolynomialRing(R, names)
        names = names + tuple('d' + n for n in names)
        if len(names) != self._n * 2:
            raise ValueError("variable names cannot differ by a leading 'd'")
        # TODO: Make this into a filtered algebra under the natural grading of
        #   x_i and dx_i have degree 1
        # Filtered is not included because it is a supercategory of super
        if R.is_field():
            cat = AlgebrasWithBasis(R).NoZeroDivisors().Super()
        else:
            cat = AlgebrasWithBasis(R).Super()
        Parent.__init__(self, base=R, names=names, category=cat)

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: DifferentialWeylAlgebra(R)
            Differential Weyl algebra of polynomials in x, y, z over Rational Field
        """
        poly_gens = ', '.join(repr(x) for x in self.variables())
        return "Differential Weyl algebra of polynomials in {} over {}".format(
            poly_gens, self.base_ring())

    # add options to class
    class options(GlobalOptions):
        r"""
        Set the global options for elements of the differential Weyl
        algebra class. The default is to have the factored
        representations turned off.

        @OPTIONS@

        If no parameters are set, then the function returns a copy of the
        options dictionary.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: D = DifferentialWeylAlgebra(R)
            sage: t,dt = D.gens()
            sage: x = dt^3*t^3 + dt^2*t^4
            sage: x
            t^3*dt^3 + t^4*dt^2 + 9*t^2*dt^2 + 8*t^3*dt + 18*t*dt + 12*t^2 + 6

            sage: D.options.factor_representation = True
            sage: x
            (12*t^2 + 6) + (8*t^3 + 18*t)*dt + (t^4 + 9*t^2)*dt^2 + (t^3)*dt^3

            sage: D.options._reset()
        """
        NAME = 'DifferentialWeylAlgebra'
        module = 'sage.algebras.weyl_algebra'
        factor_representation = {'default': False,
                                 'description': 'Controls whether to factor the differentials out or not in the output representations',
                                 'checker': lambda x: x in [True, False]}

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: a = W(2); a
            2
            sage: a.parent() is W
            True
            sage: W(x^2 - y*z)
            -y*z + x^2
        """
        t = tuple([0] * (self._n))
        if x in self.base_ring():
            if x == self.base_ring().zero():
                return self.zero()
            return self.element_class(self, {(t, t): x})
        if isinstance(x, DifferentialWeylAlgebraElement):
            R = self.base_ring()
            if x.parent().base_ring() is R:
                return self.element_class(self, dict(x))
            zero = R.zero()
            return self.element_class(self, {i: R(c) for i, c in x if R(c) != zero})
        x = self._poly_ring(x)
        return self.element_class(self, {(tuple(m), t): c
                                         for m, c in x.dict().items()})

    def _coerce_map_from_(self, R):
        """
        Return data which determines if there is a coercion map
        from ``R`` to ``self``.

        If such a map exists, the output could be a map, callable,
        or ``True``, which constructs a generic map. Otherwise the output
        must be ``False`` or ``None``.

        EXAMPLES::

            sage: R.<x,y,z> =  QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: W._coerce_map_from_(R)
            True
            sage: W._coerce_map_from_(QQ)
            True
            sage: W._coerce_map_from_(ZZ['x'])
            True

        Order of the names matter::

            sage: Wp = DifferentialWeylAlgebra(QQ['x,z,y'])
            sage: W.has_coerce_map_from(Wp)
            False
            sage: Wp.has_coerce_map_from(W)
            False

        Zero coordinates are handled appropriately::

            sage: R.<x,y,z> = ZZ[]
            sage: W3 = DifferentialWeylAlgebra(GF(3)['x,y,z'])
            sage: W3.has_coerce_map_from(R)
            True

            sage: W.<x,y,z> = DifferentialWeylAlgebra(ZZ)
            sage: W3.has_coerce_map_from(W)
            True
            sage: W3(3*x + y)
            y
        """
        if self._poly_ring.has_coerce_map_from(R):
            return True
        if isinstance(R, DifferentialWeylAlgebra):
            return (R.variable_names() == self.variable_names()
                    and self.base_ring().has_coerce_map_from(R.base_ring()))
        return super()._coerce_map_from_(R)

    def degree_on_basis(self, i):
        """
        Return the degree of the basis element indexed by ``i``.

        EXAMPLES::

            sage: W.<a,b> = DifferentialWeylAlgebra(QQ)
            sage: W.degree_on_basis( ((1, 3, 2), (0, 1, 3)) )
            10

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: dx,dy,dz = W.differentials()
            sage: elt = y*dy - (3*x - z)*dx
            sage: elt.degree()
            2
        """
        return sum(i[0]) + sum(i[1])

    def polynomial_ring(self):
        """
        Return the associated polynomial ring of ``self``.

        EXAMPLES::

            sage: W.<a,b> = DifferentialWeylAlgebra(QQ)
            sage: W.polynomial_ring()
            Multivariate Polynomial Ring in a, b over Rational Field

        ::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: W.polynomial_ring() == R
            True
        """
        return self._poly_ring

    @cached_method
    def basis(self):
        """
        Return a basis of ``self``.

        EXAMPLES::

            sage: W.<x,y> = DifferentialWeylAlgebra(QQ)
            sage: B = W.basis()
            sage: it = iter(B)
            sage: [next(it) for i in range(20)]
            [1, x, y, dx, dy, x^2, x*y, x*dx, x*dy, y^2, y*dx, y*dy,
             dx^2, dx*dy, dy^2, x^3, x^2*y, x^2*dx, x^2*dy, x*y^2]
            sage: dx, dy = W.differentials()
            sage: sorted((dx*x).monomials(), key=str)
            [1, x*dx]
            sage: B[(x*y).support()[0]]
            x*y
            sage: sorted((dx*x).monomial_coefficients().items())
            [(((0, 0), (0, 0)), 1), (((1, 0), (1, 0)), 1)]
        """
        n = self._n
        from sage.combinat.integer_lists.nn import IntegerListsNN

        def elt_map(u):
            return (tuple(u[:n]), tuple(u[n:]))

        I = IntegerListsNN(length=2 * n, element_constructor=elt_map)
        one = self.base_ring().one()

        def f(x):
            return self.element_class(self, {(x[0], x[1]): one})

        return Family(I, f, name="basis map")

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        .. SEEALSO::

            :meth:`variables`, :meth:`differentials`

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: W.algebra_generators()
            Finite family {'x': x, 'y': y, 'z': z, 'dx': dx, 'dy': dy, 'dz': dz}
        """
        d = {x: self.gen(i) for i, x in enumerate(self.variable_names())}
        return Family(self.variable_names(), lambda x: d[x])

    gens = algebra_generators

    @cached_method
    def variables(self):
        """
        Return the variables of ``self``.

        .. SEEALSO::

            :meth:`algebra_generators`, :meth:`differentials`

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: W.variables()
            Finite family {'x': x, 'y': y, 'z': z}
        """
        N = self.variable_names()[:self._n]
        d = {x: self.gen(i) for i, x in enumerate(N)}
        return Family(N, lambda x: d[x])

    @cached_method
    def differentials(self):
        """
        Return the differentials of ``self``.

        .. SEEALSO::

            :meth:`algebra_generators`, :meth:`variables`

        EXAMPLES::

            sage: W.<x,y,z> = DifferentialWeylAlgebra(QQ)
            sage: W.differentials()
            Finite family {'dx': dx, 'dy': dy, 'dz': dz}
        """
        N = self.variable_names()[self._n:]
        d = {x: self.gen(self._n + i) for i, x in enumerate(N)}
        return Family(N, lambda x: d[x])

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        .. SEEALSO::

            :meth:`algebra_generators`

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: [W.gen(i) for i in range(6)]
            [x, y, z, dx, dy, dz]
        """
        P = [0] * self._n
        D = [0] * self._n
        if i < self._n:
            P[i] = 1
        else:
            D[i - self._n] = 1
        return self.element_class(self, {(tuple(P), tuple(D)): self.base_ring().one()})

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: W.ngens()
            6
        """
        return self._n * 2

    @cached_method
    def one(self):
        """
        Return the multiplicative identity element `1`.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: W.one()
            1
        """
        t = tuple([0] * self._n)
        return self.element_class(self, {(t, t): self.base_ring().one()})

    @cached_method
    def zero(self):
        """
        Return the additive identity element `0`.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: W = DifferentialWeylAlgebra(R)
            sage: W.zero()
            0
        """
        return self.element_class(self, {})

    @lazy_attribute
    def diff_action(self):
        """
        Left action of this Weyl algebra on the underlying polynomial ring by
        differentiation.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: W = R.weyl_algebra()
            sage: dx, dy = W.differentials()
            sage: W.diff_action
            Left action by Differential Weyl algebra of polynomials in x, y
            over Rational Field on Multivariate Polynomial Ring in x, y over
            Rational Field
            sage: W.diff_action(dx^2 + dy + 1, x^3*y^3)
            x^3*y^3 + 3*x^3*y^2 + 6*x*y^3
        """
        return DifferentialWeylAlgebraAction(self)

    Element = DifferentialWeylAlgebraElement


class DifferentialWeylAlgebraAction(Action):
    """
    Left action of a Weyl algebra on its underlying polynomial ring by
    differentiation.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: W = R.weyl_algebra()
        sage: dx, dy = W.differentials()
        sage: W.diff_action
        Left action by Differential Weyl algebra of polynomials in x, y
        over Rational Field on Multivariate Polynomial Ring in x, y over
        Rational Field

    ::

        sage: g = dx^2 + x*dy
        sage: p = x^5 + x^3 + y^2*x^2 + 1
        sage: W.diff_action(g, p)
        2*x^3*y + 20*x^3 + 2*y^2 + 6*x

    The action is a left action::

        sage: h = dx*x + x*y
        sage: W.diff_action(h, W.diff_action(g, p)) == W.diff_action(h*g, p)
        True

    The action endomorphism of a differential operator::

        sage: dg = W.diff_action(g); dg
        Action of dx^2 + x*dy on Multivariate Polynomial Ring in x, y over
        Rational Field under Left action by Differential Weyl algebra...
        sage: dg(p) == W.diff_action(g, p) == g.diff(p)
        True
    """

    def __init__(self, G) -> None:
        """
        INPUT:

        - ``G`` -- Weyl algebra

        EXAMPLES::

            sage: from sage.algebras.weyl_algebra import DifferentialWeylAlgebraAction
            sage: W.<x,y> = DifferentialWeylAlgebra(QQ)
            sage: DifferentialWeylAlgebraAction(W)
            Left action by Differential Weyl algebra of polynomials in x, y
            over Rational Field on Multivariate Polynomial Ring in x, y over
            Rational Field
        """
        super().__init__(G, G.polynomial_ring(), is_left=True)

    def _act_(self, g, x):
        """
        Apply a differential operator to a polynomial.

        EXAMPLES::

            sage: W.<x,y> = DifferentialWeylAlgebra(QQ)
            sage: dx, dy = W.differentials()
            sage: W.diff_action(dx^3 + dx, x^3*y^3 + x*y)
            3*x^2*y^3 + 6*y^3 + y
        """
        f = g * x
        D = {y: c for (y, dy), c in f.monomial_coefficients(copy=False).items()
             if all(dyi == 0 for dyi in dy)}
        return self.right_domain()(D)
