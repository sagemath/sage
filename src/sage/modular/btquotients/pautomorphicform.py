# sage.doctest: needs sage.libs.pari
#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#########################################################################
r"""
Spaces of `p`-adic automorphic forms

Compute with harmonic cocycles and `p`-adic automorphic forms, including
overconvergent `p`-adic automorphic forms.

For a discussion of nearly rigid analytic modular forms and
the rigid analytic Shimura-Maass operator, see [Fra2011]_. It is worth also
looking at [FM2014]_ for information on how these are implemented in this code.

EXAMPLES:

Create a quotient of the Bruhat-Tits tree::

    sage: X = BruhatTitsQuotient(13, 11)

Declare the corresponding space of harmonic cocycles::

    sage: H = X.harmonic_cocycles(2, prec=5)

And the space of `p`-adic automorphic forms::

    sage: A = X.padic_automorphic_forms(2, prec=5, overconvergent=True)                 # needs sage.rings.padics

Harmonic cocycles, unlike `p`-adic automorphic forms, can be used to compute a basis::

    sage: a = H.gen(0)                                                                  # needs sage.rings.padics

This can then be lifted to an overconvergent `p`-adic modular form::

    sage: A.lift(a)  # long time                                                        # needs sage.rings.padics
    p-adic automorphic form of cohomological weight 0
"""

from copy import copy

import sage.modular.hecke.hecke_operator

from sage.matrix.constructor import Matrix, zero_matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.lazy_import import lazy_import
from sage.misc.verbose import verbose
from sage.modular.btquotients.btquotient import DoubleCosetReduction
from sage.modular.hecke.all import AmbientHeckeModule, HeckeModuleElement
from sage.modular.pollack_stevens.distributions import OverconvergentDistributions, Symk
from sage.modular.pollack_stevens.sigma0 import Sigma0ActionAdjuster
from sage.modules.module import Module
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.structure.element import ModuleElement
from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.rings.padics.factory', 'Qp')


# Need this to be pickleable


class _btquot_adjuster(Sigma0ActionAdjuster):
    """
    Callable object that turns matrices into 4-tuples.

    Since the modular symbol and harmonic cocycle code use different
    conventions for group actions, this function is used to make sure
    that actions are correct for harmonic cocycle computations.

    EXAMPLES::

        sage: from sage.modular.btquotients.pautomorphicform import _btquot_adjuster
        sage: adj = _btquot_adjuster()
        sage: adj(matrix(ZZ,2,2,[1..4]))
        (4, 2, 3, 1)
    """

    def __call__(self, g):
        """
        Turn matrices into 4-tuples.

        INPUT:

        - ``g`` -- a 2x2 matrix

        OUTPUT: a 4-tuple encoding the entries of ``g``

        EXAMPLES::

            sage: from sage.modular.btquotients.pautomorphicform import _btquot_adjuster
            sage: adj = _btquot_adjuster()
            sage: adj(matrix(ZZ,2,2,[0, 1, 2, 3]))
            (3, 1, 2, 0)
        """
        a, b, c, d = g.list()
        return (d, b, c, a)


def eval_dist_at_powseries(phi, f):
    """
    Evaluate a distribution on a powerseries.

    A distribution is an element in the dual of the Tate ring. The
    elements of coefficient modules of overconvergent modular symbols
    and overconvergent `p`-adic automorphic forms give examples of
    distributions in Sage.

    INPUT:

    - ``phi`` -- a distribution

    - ``f`` -- a power series over a ring coercible into a `p`-adic field

    OUTPUT:

    The value of ``phi`` evaluated at ``f``, which will be an element in the
    ring of definition of ``f``

    EXAMPLES::

        sage: from sage.modular.btquotients.pautomorphicform import eval_dist_at_powseries
        sage: R.<X> = PowerSeriesRing(ZZ,10)
        sage: f = (1 - 7*X)^(-1)

        sage: D = OverconvergentDistributions(0,7,10)                                   # needs sage.rings.padics
        sage: phi = D(list(range(1,11)))                                                # needs sage.rings.padics
        sage: eval_dist_at_powseries(phi,f)                                             # needs sage.rings.padics
        1 + 2*7 + 3*7^2 + 4*7^3 + 5*7^4 + 6*7^5 + 2*7^7 + 3*7^8 + 4*7^9 + O(7^10)
    """
    nmoments = phi.parent().precision_cap()
    K = f.parent().base_ring()
    if K.is_exact():
        K = phi.parent().base_ring()
    return sum(a * K(phi.moment(i))
               for a, i in zip(f.coefficients(), f.exponents())
               if i >= 0 and i < nmoments)


class BruhatTitsHarmonicCocycleElement(HeckeModuleElement):
    r"""
    `\Gamma`-invariant harmonic cocycles on the Bruhat-Tits
    tree. `\Gamma`-invariance is necessary so that the cocycle can be
    stored in terms of a finite amount of data.

    More precisely, given a ``BruhatTitsQuotient`` `T`, harmonic cocycles are stored as
    a list of values in some coefficient module (e.g. for weight 2 forms
    can take `\CC_p`) indexed by edges of a fundamental domain for `T` in the
    Bruhat-Tits tree. Evaluate the cocycle at other edges using Gamma
    invariance (although the values may not be equal over an orbit of
    edges as the coefficient module action may be nontrivial).

    EXAMPLES:

    Harmonic cocycles form a vector space, so they can be added and/or
    subtracted from each other::

        sage: X = BruhatTitsQuotient(5,23)
        sage: H = X.harmonic_cocycles(2,prec=10)
        sage: v1 = H.basis()[0]; v2 = H.basis()[1] # indirect doctest
        sage: v3 = v1+v2
        sage: v1 == v3-v2
        True

    and rescaled::

        sage: v4 = 2*v1
        sage: v1 == v4 - v1
        True

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self, _parent, vec):
        """
        Create a harmonic cocycle element.

        INPUT:

        - ``_parent`` : the parent space of harmonic cocycles.
        - ``vec`` : a list of elements in the coefficient module.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(31,7)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: v = H.basis()[0] # indirect doctest
            sage: TestSuite(v).run()
        """
        HeckeModuleElement.__init__(self, _parent, None)
        self._parent = _parent
        assert type(vec) is list
        assert all(v.parent() is _parent._U for v in vec)
        self._R = _parent._U.base_ring()
        self._wt = _parent._k
        self._nE = len(_parent._E)
        self._F = copy(vec)

    def _add_(self, g):
        r"""
        Add two cocycles componentwise.

        INPUT:

        - ``g`` -- a harmonic cocycle

        OUTPUT: a harmonic cocycle

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: v1 = H.basis()[0]; v2 = H.basis()[1]
            sage: v3 = v1+v2 # indirect doctest
            sage: v1 == v3-v2
            True
        """
        return self.parent()(self.element() + g.element())

    def _sub_(self, g):
        r"""
        Compute the difference of two cocycles.

        INPUT:

        - ``g`` -- a harmonic cocycle

        OUTPUT: a harmonic cocycle

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,11)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: v1 = H.basis()[0]; v2 = H.basis()[1]
            sage: v3 = v1-v2 # indirect doctest
            sage: v1 == v3+v2
            True
        """
        # Should ensure that self and g are modular forms of the same
        # weight and on the same curve
        return self.parent()(self.element() - g.element())

    def _lmul_(self, a):
        r"""
        Multiply a cocycle by a scalar.

        INPUT:

        - ``a`` -- a ring element

        OUTPUT: a harmonic cocycle

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: v1 = H.basis()[0]
            sage: v2 = 2*v1 # indirect doctest
            sage: v1 == v2-v1
            True
        """
        # Should ensure that 'a' is a scalar
        return self.parent()(a * self.element())

    def _richcmp_(self, other, op):
        r"""
        General comparison method for ``HarmonicCocycles``.

        INPUT:

        - ``other`` -- another harmonic cocycle

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: v1 = H.basis()[0]
            sage: v2 = 3*v1 # indirect doctest
            sage: 2*v1 == v2-v1
            True
        """
        if op not in [op_EQ, op_NE]:
            return NotImplemented

        b = all(self._F[e] == other._F[e] for e in range(self._nE))
        if op == op_EQ:
            return b
        return not b

    def _repr_(self):
        r"""
        Return a string describing the cocycle.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,13)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H.basis()[0] # indirect doctest
            Harmonic cocycle with values in Sym^0 Q_5^2
        """
        return 'Harmonic cocycle with values in %s' % self.parent()._U

    def monomial_coefficients(self, copy=True):
        r"""
        Return a dictionary whose keys are indices of basis elements
        in the support of ``self`` and whose values are the
        corresponding coefficients.

        EXAMPLES::

            sage: M = BruhatTitsQuotient(3,5).harmonic_cocycles(2, prec=10)
            sage: M.monomial_coefficients()
            {}
        """
        return self.element().monomial_coefficients(copy=copy)

    def print_values(self):
        r"""
        Print the values of the cocycle on all of the edges.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H.basis()[0].print_values()
            0   |1 + O(5^10)
            1   |0
            2   |0
            3   |4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
            4   |0
            5   |0
            6   |0
            7   |0
            8   |0
            9   |0
            10  |0
            11  |0
        """
        tmp = ''
        for e in range(self._nE):
            tmp += str(e) + '\t|' + str(self._F[e]) + '\n'
        print(tmp[:-1])

    def valuation(self):
        r"""
        Return the valuation of the cocycle, defined as the
        minimum of the values it takes on a set of representatives.

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: b1 = H.basis()[0]
            sage: b2 = 3*b1
            sage: b1.valuation()
            0
            sage: b2.valuation()
            1
            sage: H(0).valuation()
            +Infinity
        """
        if self == 0:
            return Infinity
        else:
            return min(self._F[e].valuation() for e in range(self._nE))

    def _compute_element(self):
        r"""
        Express a harmonic cocycle in a coordinate vector.

        OUTPUT:

        A coordinate vector encoding ``self`` in terms of the ambient
        basis in ``self.parent``

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H.basis()[0]._compute_element()
            (1 + O(3^10), 0, 0)
            sage: H.basis()[1]._compute_element()
            (0, 1 + O(3^10), 0)
            sage: H.basis()[2]._compute_element()
            (0, 0, 1 + O(3^10))
        """
        R = self._R
        A = self.parent().basis_matrix().transpose()
        B = Matrix(R, self._nE * (self.parent()._k - 1), 1,
                   [self._F[e].moment(ii) for e in range(self._nE)
                    for ii in range(self.parent()._k - 1)])
        try:
            res = (A.solve_right(B)).transpose()
        except ValueError:
            rest = (A.transpose() * A).solve_right(A.transpose() * B)
            err = A * rest - B
            if err != 0:
                try:
                    if hasattr(err.parent().base_ring().an_element(),
                               'valuation'):
                        minval = min([o.valuation() for o in err.list()
                                      if o != 0])
                    else:
                        minval = sum([RR(o.norm() ** 2) for o in err.list()])
                    verbose('Error = %s' % minval)
                except AttributeError:
                    verbose('Warning: something did not work in the '
                            'computation')
            res = rest.transpose()
        return self.parent().free_module()(res.row(0))

    # In BruhatTitsHarmonicCocycle
    def evaluate(self, e1):
        r"""
        Evaluate a harmonic cocycle on an edge of the Bruhat-Tits tree.

        INPUT:

        - ``e1`` -- a matrix corresponding to an edge of the
          Bruhat-Tits tree

        OUTPUT:

        - An element of the coefficient module of the cocycle which
          describes the value of the cocycle on ``e1``

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,17)
            sage: e0 = X.get_edge_list()[0]
            sage: e1 = X.get_edge_list()[1]
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: b = H.basis()[0]
            sage: b.evaluate(e0.rep)
            1 + O(5^10)
            sage: b.evaluate(e1.rep)
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
        """
        X = self.parent()._X
        p = X._p
        u = DoubleCosetReduction(X, e1)
        if u.label < self._nE:
            val = self._F[u.label]
        else:
            val = -self._F[u.label - self._nE]

        return u.igamma(self.parent().embed_quaternion, scale=p ** (-u.power)) * val

    # In BruhatTitsHarmonicCocycle
    def riemann_sum(self, f, center=1, level=0, E=None):
        r"""
        Evaluate the integral of the function ``f`` with respect
        to the measure determined by ``self`` over `\mathbf{P}^1(\QQ_p)`.

        INPUT:

        - ``f`` -- a function on `\mathbf{P}^1(\QQ_p)`

        - ``center`` -- integer (default: 1); Center of integration

        - ``level`` -- integer (default: 0); Determines the size of
          the covering when computing the Riemann sum. Runtime is
          exponential in the level.

        - ``E`` -- list of edges (default: ``None``); They should describe
          a covering of `\mathbf{P}^1(\QQ_p)`

        OUTPUT: a `p`-adic number

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: b = H.basis()[0]
            sage: R.<z> = PolynomialRing(QQ,1)
            sage: f = z^2

        Note that `f` has a pole at infinity, so that the result will
        be meaningless::

            sage: b.riemann_sum(f,level=0)
            1 + 5 + 2*5^3 + 4*5^4 + 2*5^5 + 3*5^6 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
        """
        R1 = LaurentSeriesRing(f.base_ring(), 'r1')
        if E is None:
            E = self.parent()._X._BT.get_balls(center, level)
        else:
            E = self.parent()._X._BT.subdivide(E, level)
        value = 0
        for e in E:
            expansion = ((R1([e[1, 1], e[1, 0]]) ** (self.parent()._k - 2) * e.determinant() ** (-(self.parent()._k - 2) / 2)) * f(R1([e[0, 1], e[0, 0]]) / R1([e[1, 1], e[1, 0]]))).truncate(self.parent()._k - 1)
            dist = self.parent()._Sigma0(e.inverse(), check=False) * self.evaluate(e)
            value += eval_dist_at_powseries(dist, expansion)
        return value

    def modular_form(self, z=None, level=0):
        r"""
        Integrate Teitelbaum's `p`-adic Poisson kernel against
        the measure corresponding to ``self`` to evaluate the associated
        modular form at ``z``.

        If ``z`` = None, a function is returned that encodes the modular form.

        .. NOTE::

            This function uses the integration method of Riemann
            summation and is incredibly slow! It should only be used for
            testing and bug-finding. Overconvergent methods are quicker.

        INPUT:

        - ``z`` -- an element in the quadratic unramified extension of
          `\QQ_p` that is not contained in `\QQ_p` (default: ``None``)

        - ``level`` -- integer; how fine of a mesh should the Riemann
          sum use

        OUTPUT: an element of the quadratic unramified extension of `\QQ_p`

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,23)
            sage: H = X.harmonic_cocycles(2,prec = 8)
            sage: b = H.basis()[0]
            sage: R.<a> = Qq(9,prec=10)
            sage: x1 = b.modular_form(a,level = 0); x1
            a + (2*a + 1)*3 + (a + 1)*3^2 + (a + 1)*3^3 + 3^4 + (a + 2)*3^5 + a*3^7 + O(3^8)
            sage: x2 = b.modular_form(a,level = 1); x2
            a + (a + 2)*3 + (2*a + 1)*3^3 + (2*a + 1)*3^4 + 3^5 + (a + 2)*3^6 + a*3^7 + O(3^8)
            sage: x3 = b.modular_form(a,level = 2); x3
            a + (a + 2)*3 + (2*a + 2)*3^2 + 2*a*3^4 + (a + 1)*3^5 + 3^6 + O(3^8)
            sage: x4 = b.modular_form(a,level = 3);x4
            a + (a + 2)*3 + (2*a + 2)*3^2 + (2*a + 2)*3^3 + 2*a*3^5 + a*3^6 + (a + 2)*3^7 + O(3^8)
            sage: (x4-x3).valuation()
            3

        TESTS:

        Check that :issue:`22634` is fixed::

            sage: X = BruhatTitsQuotient(7,2)
            sage: H = X.harmonic_cocycles(4,20)
            sage: f0, g0 = H.basis()
            sage: A = X.padic_automorphic_forms(4,20,overconvergent=True)
            sage: f = A.lift(f0).modular_form(method='moments')
            sage: T.<x> = Qq(7^2,20)
            sage: a,b,c,d = X.embed_quaternion(X.get_units_of_order()[1]).change_ring(Qp(7,20)).list()
            sage: (c*x + d)^4 * f(x) == f((a*x + b)/(c*x + d))
            True
            sage: g = A.lift(g0).modular_form(method='moments')
            sage: (c*x + d)^4 * f(x) == f((a*x + b)/(c*x + d))
            True
        """
        return self.derivative(z, level, order=0)

    # In BruhatTitsHarmonicCocycle
    def derivative(self, z=None, level=0, order=1):
        r"""
        Integrate Teitelbaum's `p`-adic Poisson kernel against
        the measure corresponding to ``self`` to evaluate the rigid
        analytic Shimura-Maass derivatives of the associated modular
        form at `z`.

        If ``z = None``, a function is returned that encodes the
        derivative of the modular form.

        .. NOTE::

            This function uses the integration method of Riemann
            summation and is incredibly slow! It should only be used for
            testing and bug-finding. Overconvergent methods are quicker.

        INPUT:

        - ``z`` -- an element in the quadratic unramified extension of
          `\QQ_p` that is not contained in `\QQ_p` (default: ``None``); if
          ``z = None`` then a function encoding the derivative is returned.

        - ``level`` -- integer; how fine of a mesh should the Riemann
          sum use

        - ``order`` -- integer; how many derivatives to take

        OUTPUT:

        An element of the quadratic unramified extension of `\QQ_p`, or
        a function encoding the derivative.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,23)
            sage: H = X.harmonic_cocycles(2,prec=5)
            sage: b = H.basis()[0]
            sage: R.<a> = Qq(9,prec=10)
            sage: b.modular_form(a,level=0) == b.derivative(a,level=0,order=0)
            True
            sage: b.derivative(a,level=1,order=1)
            (2*a + 2)*3 + (a + 2)*3^2 + 2*a*3^3 + 2*3^4 + O(3^5)
            sage: b.derivative(a,level=2,order=1)
            (2*a + 2)*3 + 2*a*3^2 + 3^3 + a*3^4 + O(3^5)
        """
        def F(z):
            R = PolynomialRing(z.parent(), 'x,y').fraction_field()
            Rx = PolynomialRing(z.parent(), 'x1').fraction_field()
            x1 = Rx.gen()
            subst = R.hom([x1, z], codomain=Rx)
            x, y = R.gens()
            center = self.parent()._X._BT.find_containing_affinoid(z)
            zbar = z.trace() - z
            f = R(1) / (x - y)
            k = self.parent()._k
            V = [f]
            for ii in range(order):
                V = [v.derivative(y) for v in V] + [k / (y - zbar) * v
                                                    for v in V]
                k += 2
            return sum([self.riemann_sum(subst(v), center, level) for v in V])
        if z is None:
            return F
        else:
            return F(z)


class BruhatTitsHarmonicCocycles(AmbientHeckeModule, UniqueRepresentation):
    r"""
    Ensure unique representation.

    EXAMPLES::

        sage: X = BruhatTitsQuotient(3,5)
        sage: M1 = X.harmonic_cocycles( 2, prec = 10)
        sage: M2 = X.harmonic_cocycles( 2, 10)
        sage: M1 is M2
        True
    """
    Element = BruhatTitsHarmonicCocycleElement

    @staticmethod
    def __classcall__(cls, X, k, prec=None, basis_matrix=None, base_field=None):
        r"""
        Represent a space of Gamma invariant harmonic
        cocycles valued in a coefficient module.

        INPUT:

        - ``X`` -- a BruhatTitsQuotient object

        - ``k`` -- integer; the weight. It must be even.

        - ``prec`` -- integer (default: ``None``); if specified, the
          precision for the coefficient module

        - ``basis_matrix`` -- a matrix (default: ``None``)

        - ``base_field`` -- a ring (default: ``None``)

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,23)
            sage: H = X.harmonic_cocycles(2,prec = 5)
            sage: H.dimension()
            3
            sage: X.genus()
            3

        Higher even weights are implemented::

            sage: H = X.harmonic_cocycles(8, prec = 10)
            sage: H.dimension()
            26

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu
        """
        return super().__classcall__(cls, X, k, prec,
                                     basis_matrix,
                                     base_field)

    def __init__(self, X, k, prec=None, basis_matrix=None, base_field=None):
        """
        Compute the space of harmonic cocycles.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,37)
            sage: H = X.harmonic_cocycles(4,prec=10)
            sage: TestSuite(H).run()
        """
        self._k = k
        self._X = X
        self._E = self._X.get_edge_list()
        self._V = self._X.get_vertex_list()

        if base_field is not None and not base_field.is_exact():
            prec = base_field.precision_cap()

        if prec is None:
            self._prec = None  # Be careful!
            if base_field is None:
                try:
                    self._R = X.get_splitting_field()
                except AttributeError:
                    raise ValueError("It looks like you are not using Magma as"
                                     " backend...and still we don't know how "
                                     "to compute splittings in that case!")
            else:
                pol = X.get_splitting_field().defining_polynomial().factor()[0][0]
                self._R = base_field.extension(pol, pol.variable_name()).absolute_field(name='r')
        else:
            self._prec = prec
            if base_field is None:
                self._R = Qp(self._X._p, prec=prec)
            else:
                self._R = base_field

        self._U = Symk(self._k - 2, base=self._R, act_on_left=True,
                       adjuster=_btquot_adjuster(),
                       dettwist=-ZZ((self._k - 2) // 2), act_padic=True)

        if basis_matrix is None:
            self.__rank = self._X.dimension_harmonic_cocycles(self._k)
        else:
            self.__rank = basis_matrix.nrows()
        if basis_matrix is not None:
            self.__matrix = basis_matrix
            self.__matrix.set_immutable()
            assert self.__rank == self.__matrix.nrows()

        self._Sigma0 = self._U._act._Sigma0

        AmbientHeckeModule.__init__(self, self._R, self.__rank,
                                    self._X.prime() * self._X.Nplus() * self._X.Nminus(), weight=self._k)
        self._populate_coercion_lists_()

    def monomial_coefficients(self):
        r"""
        Void method to comply with pickling.

        EXAMPLES::

            sage: M = BruhatTitsQuotient(3,5).harmonic_cocycles(2,prec=10)
            sage: M.monomial_coefficients()
            {}
        """
        return {}

    def base_extend(self, base_ring):
        r"""
        Extend the base ring of the coefficient module.

        INPUT:

        - ``base_ring`` -- a ring that has a coerce map from the
          current base ring

        OUTPUT: a new space of HarmonicCocycles with the base extended

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,19)
            sage: H = X.harmonic_cocycles(2,10)
            sage: H.base_ring()
            3-adic Field with capped relative precision 10
            sage: H1 = H.base_extend(Qp(3,prec=15))
            sage: H1.base_ring()
            3-adic Field with capped relative precision 15
        """
        if not base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError("No coercion defined")
        else:
            return self.change_ring(base_ring)

    def change_ring(self, new_base_ring):
        r"""
        Change the base ring of the coefficient module.

        INPUT:

        - ``new_base_ring`` -- a ring that has a coerce map from the
          current base ring

        OUTPUT: new space of HarmonicCocycles with different base ring

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,17)
            sage: H = X.harmonic_cocycles(2,10)
            sage: H.base_ring()
            5-adic Field with capped relative precision 10
            sage: H1 = H.base_extend(Qp(5,prec=15)) # indirect doctest
            sage: H1.base_ring()
            5-adic Field with capped relative precision 15
        """
        if not new_base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError("No coercion defined")

        basis_matrix = self.basis_matrix().change_ring(new_base_ring)
        basis_matrix.set_immutable()
        return self.__class__(self._X, self._k, prec=None,
                              basis_matrix=basis_matrix,
                              base_field=new_base_ring)

    def rank(self):
        r"""
        Return the rank (dimension) of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,11)
            sage: H = X.harmonic_cocycles(2,prec = 10)
            sage: X.genus() == H.rank()
            True
            sage: H1 = X.harmonic_cocycles(4,prec = 10)
            sage: H1.rank()
            16
        """
        return self.__rank

    def submodule(self, v, check=False):
        r"""
        Return the submodule of ``self`` spanned by ``v``.

        INPUT:

        - ``v`` -- submodule of ``self.free_module()``

        - ``check`` -- boolean (default: ``False``)

        OUTPUT: subspace of harmonic cocycles

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H.rank()
            3
            sage: v = H.gen(0)
            sage: N = H.free_module().span([v.element()])
            sage: H1 = H.submodule(N)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # return BruhatTitsHarmonicCocyclesSubmodule(self, v)
        raise NotImplementedError

    def is_simple(self):
        r"""
        Whether ``self`` is irreducible.

        OUTPUT: boolean; ``True`` if and only if ``self`` is irreducible

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3, 29)
            sage: H = X.harmonic_cocycles(4, prec=10)
            sage: H.rank()
            14
            sage: H.is_simple()
            False
            sage: X = BruhatTitsQuotient(7, 2)
            sage: H = X.harmonic_cocycles(2, prec=10)
            sage: H.rank()
            1
            sage: H.is_simple()
            True
        """
        return self.rank() == 1

    def _repr_(self):
        r"""
        This returns the representation of ``self`` as a string.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H
            Space of harmonic cocycles of weight 2 on Quotient of the Bruhat
            Tits tree of GL_2(QQ_5) with discriminant 23 and level 1
        """
        return 'Space of harmonic cocycles of weight %s on %s' % (self._k,
                                                                  self._X)

    def _latex_(self):
        r"""
        A LaTeX representation of ``self``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: latex(H) # indirect doctest
            \text{Space of harmonic cocycles of weight } 2 \text{ on } X(5 \cdot 23,1)\otimes_{\Bold{Z}} \Bold{F}_{5}
        """
        s = '\\text{Space of harmonic cocycles of weight } '
        s += (self._k)._latex_() + ' \\text{ on } ' + self._X._latex_()
        return s

    def _an_element_(self):
        r"""
        Return an element of the ambient space.

        OUTPUT: a harmonic cocycle in self

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H.an_element() # indirect doctest
            Harmonic cocycle with values in Sym^0 Q_5^2
        """
        return self.basis()[0]

    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other BruhatTitsHarmonicCocycles or from
        pAdicAutomorphicForms, also from 0.

        OUTPUT: boolean; ``True`` if and only if ``self`` is a space of
        BruhatTitsHarmonicCocycles or pAdicAutomorphicForms

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: A(H.basis()[0]) # indirect doctest
            p-adic automorphic form of cohomological weight 0
        """
        if isinstance(S, (BruhatTitsHarmonicCocycles, pAdicAutomorphicForms)):
            if S._k != self._k:
                return False
            if S._X != self._X:
                return False
            return True
        return False

    def __eq__(self, other):
        r"""
        Test whether two BruhatTitsHarmonicCocycle spaces are equal.

        INPUT:

        - ``other`` -- a BruhatTitsHarmonicCocycles class

        OUTPUT: boolean

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: H1 = X.harmonic_cocycles(2,prec=10)
            sage: H2 = X.harmonic_cocycles(2,prec=10)
            sage: H1 == H2
            True
        """
        if not isinstance(other, BruhatTitsHarmonicCocycles):
            return False

        return (self.base_ring() == other.base_ring() and
                self._X == other._X and
                self._k == other._k)

    def __ne__(self, other):
        r"""
        Test whether two BruhatTitsHarmonicCocycle spaces are not equal.

        INPUT:

        - ``other`` -- a BruhatTitsHarmonicCocycles class

        OUTPUT: boolean

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: H1 = X.harmonic_cocycles(2,prec=10)
            sage: H2 = X.harmonic_cocycles(2,prec=10)
            sage: H1 != H2
            False
        """
        return not self.__eq__(other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: H1 = X.harmonic_cocycles(2,prec=10)
            sage: H2 = X.harmonic_cocycles(2,prec=10)
            sage: hash(H1) == hash(H2)
            True
        """
        return hash((self.base_ring(), self._X, self._k))

    def _element_constructor_(self, x):
        r"""
        Constructor for harmonic cocycles.

        INPUT:

        - ``x`` -- an object coercible into a harmonic cocycle

        OUTPUT: a harmonic cocycle

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H(H.an_element()) # indirect doctest
            Harmonic cocycle with values in Sym^0 Q_3^2
            sage: H(0)
            Harmonic cocycle with values in Sym^0 Q_3^2
        """
        if type(x) is sage.modules.free_module_element.FreeModuleElement_generic_dense:
            vmat = MatrixSpace(self._R, 1, self.dimension())(x)
            tmp = (vmat * self.ambient_module().basis_matrix()).row(0)
            vec = [self._U(tmp[e * (self._k - 1):(e + 1) * (self._k - 1)])
                   for e in range(len(self._E))]
            return self.element_class(self, vec)

        if type(x) is list:
            return self.element_class(self, [self._U(o) for o in x])

        if hasattr(x, 'parent'):
            parent = x.parent()
            if isinstance(parent, BruhatTitsHarmonicCocycles):
                return self.element_class(self, [self._U(o) for o in x._F])
            elif isinstance(parent, pAdicAutomorphicForms):
                tmp = [self._E[ii].rep * self._U(x._F[ii]) for ii in range(self._nE)]
                return self.element_class(self, tmp)
        if x == 0:
            tmp = [self._U([0] * (self.weight() - 1))] * self._X._num_edges
            return self.element_class(self, tmp)
        else:
            raise TypeError

    def free_module(self):
        r"""
        Return the underlying free module.

        OUTPUT: a free module

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H.free_module()
            Vector space of dimension 1 over 3-adic Field with
            capped relative precision 10
        """
        try:
            return self.__free_module
        except AttributeError:
            pass
        V = self.base_ring() ** self.dimension()
        self.__free_module = V
        return V

    def character(self):
        r"""
        The trivial character.

        OUTPUT: the identity map

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: H = X.harmonic_cocycles(2,prec = 10)
            sage: f = H.character()
            sage: f(1)
            1
            sage: f(2)
            2
        """
        return lambda x: x

    def embed_quaternion(self, g, scale=1, exact=None):
        r"""
        Embed the quaternion element ``g`` into the matrix algebra.

        INPUT:

        - ``g`` -- a quaternion, expressed as a 4x1 matrix

        OUTPUT: a 2x2 matrix with `p`-adic entries

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,2)
            sage: q = X.get_stabilizers()[0][1][0]
            sage: H = X.harmonic_cocycles(2,prec = 5)
            sage: Hmat = H.embed_quaternion(q)
            sage: Hmat.matrix().trace() == X._conv(q).reduced_trace() and Hmat.matrix().determinant() == 1
            True
        """
        if exact is None:
            exact = self._R.is_exact()
        return self._Sigma0(scale * self._X.embed_quaternion(g, exact=exact,
                                                             prec=self._prec),
                            check=False)

    def basis_matrix(self):
        r"""
        Return a basis of ``self`` in matrix form.

        If the coefficient module `M` is of finite rank then the space
        of Gamma invariant `M` valued harmonic cocycles can be
        represented as a subspace of the finite rank space of all
        functions from the finitely many edges in the corresponding
        BruhatTitsQuotient into `M`. This function computes this
        representation of the space of cocycles.

        OUTPUT:

        - A basis matrix describing the cocycles in the spaced of all
          `M` valued Gamma invariant functions on the tree.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,3)
            sage: M = X.harmonic_cocycles(4,prec = 20)
            sage: B = M.basis() # indirect doctest
            sage: len(B) == X.dimension_harmonic_cocycles(4)
            True

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        try:
            return self.__matrix
        except AttributeError:
            pass
        nV = len(self._V)
        nE = len(self._E)
        stab_conds = []
        S = self._X.get_edge_stabilizers()
        p = self._X._p
        d = self._k - 1
        for e in self._E:
            try:
                g = next(g for g in S[e.label] if g[2])
                C = self._U.acting_matrix(self._Sigma0(self.embed_quaternion(g[0])), d).transpose()  # Warning - Need to allow the check = True
                C -= self._U.acting_matrix(self._Sigma0(Matrix(QQ, 2, 2, p ** g[1])), d).transpose()  # Warning - Need to allow the check = True
                stab_conds.append([e.label, C])
            except StopIteration:
                pass

        n_stab_conds = len(stab_conds)
        self._M = Matrix(self._R, (nV + n_stab_conds) * d, nE * d, 0,
                         sparse=True)
        for v in self._V:
            for e in v.leaving_edges:
                if e.parity:
                    continue
                C = sum([self._U.acting_matrix(self.embed_quaternion(x[0]), d)
                         for x in e.links],
                        Matrix(self._R, d, d, 0)).transpose()
                self._M.set_block(v.label * d, e.label * d, C)
            for e in v.entering_edges:
                if e.parity:
                    continue
                C = sum([self._U.acting_matrix(self.embed_quaternion(x[0]), d)
                         for x in e.opposite.links],
                        Matrix(self._R, d, d, 0)).transpose()
                self._M.set_block(v.label * d, e.opposite.label * d, C)

        for kk in range(n_stab_conds):
            v = stab_conds[kk]
            self._M.set_block((nV + kk) * d, v[0] * d, v[1])

        x1 = self._M.right_kernel().matrix()

        if x1.nrows() != self.rank():
            raise RuntimeError('The computed dimension does not agree with '
                               'the expectation. Consider increasing '
                               'precision!')

        K = [c.list() for c in x1.rows()]

        if not self._R.is_exact():
            for ii in range(len(K)):
                s = min([t.valuation() for t in K[ii]])
                for jj in range(len(K[ii])):
                    K[ii][jj] = (p ** (-s)) * K[ii][jj]

        self.__matrix = Matrix(self._R, len(K), nE * d, K)
        self.__matrix.set_immutable()
        return self.__matrix

    def __apply_atkin_lehner(self, q, f):
        r"""
        Apply an Atkin-Lehner involution to a harmonic cocycle.

        INPUT:

        - ``q`` -- integer dividing the full level ``p*Nminus*Nplus``

        - ``f`` -- a harmonic cocycle

        OUTPUT:

        The harmonic cocycle obtained by hitting ``f`` with the
        Atkin-Lehner at ``q``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,17)
            sage: H = X.harmonic_cocycles(2,prec = 10)
            sage: A = H.atkin_lehner_operator(5).matrix() # indirect doctest
            sage: A**2 == 1
            True
        """
        Data = self._X._get_atkin_lehner_data(q)
        p = self._X._p
        tmp = [self._U(0) for jj in range(len(self._E))]
        d1 = Data[1]
        mga = self.embed_quaternion(Data[0])
        nE = len(self._E)
        for jj in range(nE):
            t = d1[jj]
            if t.label < nE:
                tmp[jj] += mga * t.igamma(self.embed_quaternion, scale=p ** -t.power) * f._F[t.label]
            else:
                tmp[jj] += mga * t.igamma(self.embed_quaternion, scale=p ** -t.power) * (-f._F[t.label - nE])

        return self(tmp)

    def __apply_hecke_operator(self, l, f):
        r"""
        This function applies a Hecke operator to a harmonic cocycle.

        INPUT:

        - ``l`` -- integer

        - ``f`` -- a harmonic cocycle

        OUTPUT:

        A harmonic cocycle which is the result of applying the `l`-th
        Hecke operator to ``f``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,17)
            sage: H = X.harmonic_cocycles(2,prec=50)
            sage: A = H.hecke_operator(7).matrix() # indirect doctest
            sage: [o.rational_reconstruction() for o in A.charpoly().coefficients()]
            [-8, -12, 12, 20, 8, 1]
        """
        HeckeData, alpha = self._X._get_hecke_data(l)
        if self.level() % l == 0:
            factor = QQ(l ** (Integer((self._k - 2) // 2)) / (l + 1))
        else:
            factor = QQ(l ** (Integer((self._k - 2) // 2)))
        p = self._X._p
        alphamat = self.embed_quaternion(alpha)
        tmp = [self._U(0) for jj in range(len(self._E))]
        for d0, d1 in HeckeData:
            mga = self.embed_quaternion(d0) * alphamat
            nE = len(self._E)
            for jj in range(nE):
                t = d1[jj]
                if t.label < nE:
                    tmp[jj] += mga * t.igamma(self.embed_quaternion, scale=p ** -t.power) * f._F[t.label]
                else:
                    tmp[jj] += mga * t.igamma(self.embed_quaternion, scale=p ** -t.power) * (-f._F[t.label - nE])
        return self([factor * x for x in tmp])

    def _compute_atkin_lehner_matrix(self, d):
        r"""
        When the underlying coefficient module is finite, this
        function computes the matrix of an Atkin-Lehner involution in
        the basis provided by the function basis_matrix

        INPUT:

        - ``d`` -- integer dividing ``p*Nminus*Nplus``, where these
          quantities are associated to the BruhatTitsQuotient ``self._X``

        OUTPUT:

        The matrix of the Atkin-Lehner involution at ``d`` in the basis given by
        ``self.basis_matrix``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,13)
            sage: H = X.harmonic_cocycles(2,prec=5)
            sage: A = H.atkin_lehner_operator(5).matrix() # indirect doctest
            sage: A**2 == 1
            True
        """
        return self.__compute_operator_matrix(lambda f: self.__apply_atkin_lehner(d, f))

    def _compute_hecke_matrix_prime(self, l):
        r"""
        When the underlying coefficient module is finite, this
        function computes the matrix of a (prime) Hecke operator in
        the basis provided by the function basis_matrix

        INPUT:

        - ``l`` -- prime integer

        OUTPUT:

        The matrix of `T_l` acting on the cocycles in the basis given by
        ``self.basis_matrix``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,11)
            sage: H = X.harmonic_cocycles(4,prec=60)
            sage: A = H.hecke_operator(7).matrix() # long time, indirect doctest
            sage: [o.rational_reconstruction() for o in A.charpoly().coefficients()] # long time
            [6496256, 1497856, -109040, -33600, -904, 32, 1]
        """
        return self.__compute_operator_matrix(lambda f: self.__apply_hecke_operator(l, f))

    def __compute_operator_matrix(self, T):
        r"""
        Compute the matrix of the operator `T`.

        Used primarily to compute matrices of Hecke operators
        in a streamlined way.

        INPUT:

        - ``T`` -- a linear function on the space of harmonic cocycles

        OUTPUT: the matrix of ``T`` acting on the space of harmonic cocycles

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: A = H.hecke_operator(11).matrix() # indirect doctest
            sage: [o.rational_reconstruction() for o in A.charpoly().coefficients()]
            [-12, -1, 4, 1]
        """
        R = self._R
        A = self.basis_matrix().transpose()
        basis = self.basis()
        B = zero_matrix(R, len(self._E) * (self._k - 1), self.dimension())
        for rr in range(len(basis)):
            g = T(basis[rr])
            B.set_block(0, rr, Matrix(R, len(self._E) * (self._k - 1), 1, [g._F[e].moment(ii) for e in range(len(self._E)) for ii in range(self._k - 1)]))
        try:
            res = (A.solve_right(B)).transpose()
        except ValueError:
            rest = (A.transpose() * A).solve_right(A.transpose() * B)
            err = A * rest - B
            if err != 0:
                try:
                    if hasattr(err.parent().base_ring().an_element(),
                               'valuation'):
                        minval = min([o.valuation() for o in err.list()
                                      if o != 0])
                    else:
                        minval = sum([RR(o.norm() ** 2) for o in err.list()])
                    verbose('Error = %s' % minval)
                except AttributeError:
                    verbose('Warning: something did not work in the computation')
            res = rest.transpose()
        res.set_immutable()
        return res

# class BruhatTitsHarmonicCocyclesSubmodule(BruhatTitsHarmonicCocycles,sage.modular.hecke.submodule.HeckeSubmodule):
#     r"""
#     Submodule of a space of BruhatTitsHarmonicCocycles.
#
#     INPUT:
#
#     - ``x`` -- integer (default: 1)
#
#     - ``y`` -- integer (default: 2)
#
#     EXAMPLES::
#
#         sage: X = BruhatTitsQuotient(3,17)
#         sage: H = X.harmonic_cocycles(2,prec=10)
#         sage: N = H.free_module().span([H.an_element().element()])
#         sage: H1 = H.submodule(N) # indirect doctest
#         sage: H1
#         Subspace of Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_3) with discriminant 17 and level 1 of dimension 1
#
#     AUTHOR:
#
#     - Marc Masdeu (2012-02-20)
#     """
#     def __init__(self, ambient_module, submodule, check):
#         """
#         Submodule of harmonic cocycles.
#
#         INPUT:
#
#         - ``ambient_module`` -- BruhatTitsHarmonicCocycles
#
#         - ``submodule`` -- submodule of the ambient space
#
#         - ``check`` -- boolean (default: ``False``); whether to check that the
#           submodule is Hecke equivariant
#
#         EXAMPLES::
#
#             sage: X = BruhatTitsQuotient(3,17)
#             sage: H = X.harmonic_cocycles(2,prec=10)
#             sage: N = H.free_module().span([H.an_element().element()])
#             sage: H1 = H.submodule(N)
#             sage: TestSuite(H1).run()
#         """
#         A = ambient_module
#         self.__rank = submodule.dimension()
#         basis_matrix = submodule.basis_matrix()*A.basis_matrix()
#         basis_matrix.set_immutable()
#         BruhatTitsHarmonicCocycles.__init__(self,A._X,A._k,A._prec,basis_matrix,A.base_ring())
#
#     def rank(self):
#         r"""
#         Returns the rank (dimension) of the submodule.
#
#         OUTPUT:
#
#         Integer - The rank of ``self``.
#
#         EXAMPLES::
#
#             sage: X = BruhatTitsQuotient(3,17)
#             sage: H = X.harmonic_cocycles(2,prec=10)
#             sage: N = H.free_module().span([H.an_element().element()])
#             sage: H1 = H.submodule(basis = [H.an_element()])
#             sage: H1.rank()
#             1
#         """
#         return self.__rank
#
#     def _repr_(self):
#         r"""
#         Returns the representation of self as a string.
#
#         OUTPUT:
#
#         String representation of self.
#
#         EXAMPLES::
#
#             sage: X = BruhatTitsQuotient(3,17)
#             sage: H = X.harmonic_cocycles(2,prec=10)
#             sage: N = H.free_module().span([H.an_element().element()])
#             sage: H1=H.submodule(N)
#             sage: H1
#             Subspace of Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_3) with discriminant 17 and level 1 of dimension 1
#         """
#         return "Subspace of %s of dimension %s"%(self.ambient(),self.dimension())


class pAdicAutomorphicFormElement(ModuleElement):
    r"""
    Rudimentary implementation of a class for a `p`-adic
    automorphic form on a definite quaternion algebra over `\QQ`. These
    are required in order to compute moments of measures associated to
    harmonic cocycles on the Bruhat-Tits tree using the overconvergent modules
    of Darmon-Pollack and Matt Greenberg. See Greenberg's thesis [Gr2006]_ for
    more details.

    INPUT:

    - ``vec`` -- a preformatted list of data

    EXAMPLES::

        sage: X = BruhatTitsQuotient(17,3)
        sage: H = X.harmonic_cocycles(2,prec=10)
        sage: h = H.an_element()
        sage: HH = X.padic_automorphic_forms(2,10)
        sage: a = HH(h)
        sage: a
        p-adic automorphic form of cohomological weight 0

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self, parent, vec):
        """
        Create a pAdicAutomorphicFormElement.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(17,3)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: TestSuite(A.an_element()).run()
        """
        self._num_generators = len(parent._list)
        self._cached_values = {}
        self._R = Qp(parent.prime(), prec=parent._prec)
        self._value = [parent._U(v) for v in vec]
        ModuleElement.__init__(self, parent)

    def _add_(self, g):
        r"""
        This function adds two `p`-adic automorphic forms.

        INPUT:

        - ``g`` -- a `p`-adic automorphic form

        OUTPUT: the result of adding ``g`` to ``self``

        EXAMPLES::

            sage: X = BruhatTitsQuotient(17,3)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: a = A.an_element()
            sage: b = a + a # indirect doctest
        """
        # Should ensure that self and g are of the same weight and on
        # the same curve
        vec = [self._value[e] + g._value[e]
               for e in range(self._num_generators)]
        return self.parent()(vec)

    def _sub_(self, g):
        r"""
        This function subtracts a `p`-adic automorphic form from another.

        INPUT:

        - ``g`` -- a `p`-adic automorphic form

        OUTPUT: the result of subtracting ``g`` from ``self``

        EXAMPLES::

            sage: X = BruhatTitsQuotient(17,3)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: a = A.an_element()
            sage: b = a - a # indirect doctest
            sage: b == 0
            True
        """
        # Should ensure that self and g are of the same weight and on
        # the same curve
        vec = [self._value[e] - g._value[e]
               for e in range(self._num_generators)]
        return self.parent()(vec)

    def _richcmp_(self, other, op):
        r"""
        Test for equality of pAdicAutomorphicForm elements.

        INPUT:

        - ``other`` -- another `p`-automorphic form

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,23)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: v1 = A(H.basis()[0])
            sage: v2 = 3*v1
            sage: 2*v1 == v2-v1 # indirect doctest
            True
        """
        if op not in [op_EQ, op_NE]:
            return NotImplemented

        b = all(self._value[e] == other._value[e]
                for e in range(self._num_generators))
        if op == op_EQ:
            return b
        return not b

    def __bool__(self):
        """
        Tell whether the form is zero or not.

        OUTPUT: boolean; whether ``self`` is zero

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,23)
            sage: H = X.harmonic_cocycles(4,prec = 20)
            sage: A = X.padic_automorphic_forms(4,prec = 20)
            sage: v1 = A(H.basis()[1])
            sage: bool(v1)
            True
            sage: v2 = v1-v1
            sage: bool(v2)
            False
        """
        return any(not o.is_zero() for o in self._value)

    def __getitem__(self, e1):
        r"""
        Evaluate a `p`-adic automorphic form on a matrix in `GL_2(\QQ_p)`.

        INPUT:

        - ``e1`` -- a matrix in `GL_2(\QQ_p)`

        OUTPUT: the value of ``self`` evaluated on ``e1``

        EXAMPLES::

            sage: X = BruhatTitsQuotient(17,3)
            sage: M = X.harmonic_cocycles(2,prec=5)
            sage: A = X.padic_automorphic_forms(2,prec=5)
            sage: a = A(M.gen(0))
            sage: a[Matrix(ZZ,2,2,[1,2,3,4])]
            8 + 8*17 + 8*17^2 + 8*17^3 + 8*17^4 + O(17^5)
        """
        return self.evaluate(e1)

    def evaluate(self, e1):
        r"""
        Evaluate a `p`-adic automorphic form on a matrix in `GL_2(\QQ_p)`.

        INPUT:

        - ``e1`` -- a matrix in `GL_2(\QQ_p)`

        OUTPUT: the value of ``self`` evaluated on ``e1``

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,5)
            sage: M = X.harmonic_cocycles(2,prec=5)
            sage: A = X.padic_automorphic_forms(2,prec=5)
            sage: a = A(M.basis()[0])
            sage: a.evaluate(Matrix(ZZ,2,2,[1,2,3,1]))
            4 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + O(7^5)
            sage: a.evaluate(Matrix(ZZ,2,2,[17,0,0,1]))
            1 + O(7^5)
        """
        X = self.parent()._source
        p = self.parent().prime()
        u = DoubleCosetReduction(X, e1)
        tmp = ((u.t(self.parent()._U.base_ring().precision_cap())) * p ** (u.power)).adjugate()
        S0 = self.parent()._Sigma0
        return S0(tmp, check=False) * self._value[u.label]
        # Warning! Should remove check=False...

    def _lmul_(self, a):
        r"""
        Multiply the automorphic form by a scalar.

        INPUT:

        - ``a`` -- a scalar

        EXAMPLES::

            sage: X = BruhatTitsQuotient(17,3)
            sage: M = X.harmonic_cocycles(2,prec=5)
            sage: A = X.padic_automorphic_forms(2,prec=5)
            sage: a = A(M.basis()[0])
            sage: a.evaluate(Matrix(ZZ,2,2,[1,2,3,4]))
            8 + 8*17 + 8*17^2 + 8*17^3 + 8*17^4 + O(17^5)
            sage: b = 2*a # indirect doctest
            sage: b.evaluate(Matrix(ZZ,2,2,[1,2,3,4]))
            16 + 16*17 + 16*17^2 + 16*17^3 + 16*17^4 + O(17^5)
        """
        # Should ensure that 'a' is a scalar
        return self.parent()([a * self._value[e]
                              for e in range(self._num_generators)])

    def _repr_(self):
        r"""
        This returns the representation of ``self`` as a string.

        If ``self`` corresponds to a modular form of weight `k`, then the
        cohomological weight is `k-2`.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(17,3)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: a = A.an_element()
            sage: a # indirect doctest
            p-adic automorphic form of cohomological weight 0
        """
        return 'p-adic automorphic form of cohomological weight %s' % self.parent()._U.weight()

    def valuation(self):
        r"""
        The valuation of ``self``, defined as the minimum of the
        valuations of the values that it takes on a set of edge
        representatives.

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(17,3)
            sage: M = X.harmonic_cocycles(2,prec=10)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: a = A(M.gen(0))
            sage: a.valuation()
            0
            sage: (17*a).valuation()
            1
        """
        return min(self._value[e].valuation()
                   for e in range(self._num_generators))

    def _improve(self, hc):
        r"""
        Repeatedly apply the `U_p` operator to a `p`-adic
        automorphic form. This is used to compute moments of a measure
        associated to a rigid modular form in the following way: lift
        a rigid modular form to an overconvergent `p`-adic
        automorphic form in any way, and then repeatedly apply `U_p`
        to project to the ordinary part.  The resulting form encodes
        the moments of the measure of the original rigid modular form
        (assuming it is ordinary).

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,2)
            sage: H = X.harmonic_cocycles(2,prec = 10)
            sage: h = H.gen(0)
            sage: A = X.padic_automorphic_forms(2,prec = 10,overconvergent=True)
            sage: a = A.lift(h) # indirect doctest

        REFERENCES:

        For details see [Gr2006]_. Alternatively, one can look at
        [DP]_ for the analogous algorithm in the case of modular symbols.

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu
        """
        MMM = self.parent()
        U = MMM._U

        h1 = MMM([o.lift(M=MMM.precision_cap()) for o in self._value])
        h2 = MMM._apply_Up_operator(h1, True)
        verbose("Applied Up once")
        ii = 0
        current_val = 0
        init_val = self.valuation()
        old_val = init_val - 1
        while current_val > old_val:
            old_val = current_val
            ii += 1
            h1._value = [U(c) for c in h2._value]
            h2 = MMM._apply_Up_operator(h1, True)
            current_val = (h2 - h1).valuation() - init_val
            verbose('val  = %s' % current_val)
            if current_val is Infinity:
                break
            verbose('Applied Up %s times' % (ii + 1))
        return h2

    def integrate(self, f, center=1, level=0, method='moments'):
        r"""
        Calculate

        .. MATH::

            \int_{\mathbf{P}^1(\QQ_p)} f(x)d\mu(x)

        were `\mu` is the measure associated to ``self``.

        INPUT:

        - ``f`` -- an analytic function

        - ``center`` -- 2x2 matrix over `\QQ_p` (default: 1)

        - ``level`` -- integer (default: 0)

        - ``method`` -- string (default: ``'moments'``); which method of
          integration to use. Either ``'moments'`` or ``'riemann_sum'``.

        EXAMPLES:

        Integrating the Poisson kernel against a measure yields a
        value of the associated modular form. Such values can be
        computed efficiently using the overconvergent method, as long
        as one starts with an ordinary form::

            sage: X = BruhatTitsQuotient(7,2)
            sage: X.genus()
            1

        Since the genus is 1, the space of weight 2 forms is 1
        dimensional. Hence any nonzero form will be a `U_7`
        eigenvector. By Jacquet-Langlands and Cerednik-Drinfeld, in
        this case the Hecke eigenvalues correspond to that of any
        nonzero form on `\Gamma_0(14)` of weight `2`. Such a form is
        ordinary at `7`, and so we can apply the overconvergent method
        directly to this form without `p`-stabilizing::

            sage: H = X.harmonic_cocycles(2,prec = 5)
            sage: h = H.gen(0)
            sage: A = X.padic_automorphic_forms(2,prec = 5,overconvergent=True)
            sage: a = A.lift(h)
            sage: a._value[0].moment(2)
            2 + 6*7 + 4*7^2 + 4*7^3 + 6*7^4 + O(7^5)

        Now that we've lifted our harmonic cocycle to an
        overconvergent automorphic form we simply need to define the
        Teitelbaum-Poisson Kernel, and then integrate::

            sage: Kp.<x> = Qq(49,prec = 5)
            sage: z = Kp['z'].gen()
            sage: f = 1/(z-x)
            sage: a.integrate(f)
            (5*x + 5) + (4*x + 4)*7 + (5*x + 5)*7^2 + (5*x + 6)*7^3 + O(7^5)

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        E = self.parent()._source._BT.get_balls(center, level)
        R1 = LaurentSeriesRing(f.base_ring(), 'r1', default_prec=self.parent()._U.base_ring().precision_cap() + 1)
        R2 = PolynomialRing(f.base_ring(), 'x')
        x = R2.gen()
        value = 0
        if method == 'riemann_sum':
            for e in E:
                exp = ((R1([e[1, 1], e[1, 0]])) ** (self.parent()._U.weight()) * e.determinant() ** (-(self.parent()._U.weight()) / 2)) * f(R1([e[0, 1], e[0, 0]]) / R1([e[1, 1], e[1, 0]]))
                new = eval_dist_at_powseries(self.evaluate(e), exp.truncate(self.parent()._U.weight() + 1))
                value += new
        elif method == 'moments':
            n = self.parent()._U.weight()
            for e in E:
                a, b, c, d = e.list()
                delta = e.determinant()
                verbose('%s' % (R2([e[0, 1], e[0, 0]])
                                / R2([e[1, 1], e[1, 0]])))
                tmp = ((c * x + d) ** n * delta ** -ZZ(n // 2)) * f((a * x + b) / (c * x + d))
                exp = R1(tmp.numerator()) / R1(tmp.denominator())
                new = eval_dist_at_powseries(self.evaluate(e), exp)

                value += new
        else:
            print('The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should never be used.')
            return False
        return value

    def modular_form(self, z=None, level=0, method='moments'):
        r"""
        Return the modular form corresponding to ``self``.

        INPUT:

        - ``z`` -- (default: ``None``) if specified, returns the value of
          the form at the point ``z`` in the `p`-adic upper half
          plane

        - ``level`` -- integer (default: 0); if ``method`` is
          'riemann_sum', will use a covering of `P^1(\QQ_p)` with
          balls of size `p^-\mbox{level}`

        - ``method`` -- string (default: ``'moments'``); it must be
          either ``'moments'`` or ``'riemann_sum'``

        OUTPUT:

        A function from the `p`-adic upper half plane to `\CC_p`. If
        an argument ``z`` was passed, returns instead the value at
        that point.

        EXAMPLES:

        Integrating the Poisson kernel against a measure yields a
        value of the associated modular form. Such values can be
        computed efficiently using the overconvergent method, as long
        as one starts with an ordinary form::

            sage: X = BruhatTitsQuotient(7, 2)
            sage: X.genus()
            1

        Since the genus is 1, the space of weight 2 forms is 1
        dimensional. Hence any nonzero form will be a `U_7`
        eigenvector. By Jacquet-Langlands and Cerednik-Drinfeld, in
        this case the Hecke eigenvalues correspond to that of any
        nonzero form on `\Gamma_0(14)` of weight `2`. Such a form is
        ordinary at `7`, and so we can apply the overconvergent method
        directly to this form without `p`-stabilizing::

            sage: H = X.harmonic_cocycles(2,prec = 5)
            sage: A = X.padic_automorphic_forms(2,prec = 5,overconvergent=True)
            sage: f0 = A.lift(H.basis()[0])

        Now that we've lifted our harmonic cocycle to an
        overconvergent automorphic form, we extract the associated
        modular form as a function and test the modular property::

            sage: T.<x> = Qq(7^2,prec = 5)
            sage: f = f0.modular_form(method = 'moments')
            sage: a,b,c,d = X.embed_quaternion(X.get_units_of_order()[1]).change_ring(T.base_ring()).list()
            sage: ((c*x + d)^2*f(x)-f((a*x + b)/(c*x + d))).valuation()
            5
        """
        return self.derivative(z, level, method, order=0)

    def derivative(self, z=None, level=0, method='moments', order=1):
        r"""
        Return the derivative of the modular form corresponding to
        ``self``.

        INPUT:

        - ``z`` -- (default: ``None``) if specified, evaluates the derivative
          at the point ``z`` in the `p`-adic upper half plane

        - ``level`` -- integer (default: 0); if ``method`` is
          'riemann_sum', will use a covering of `P^1(\QQ_p)` with
          balls of size `p^-\mbox{level}`.

        - ``method`` -- string (default: ``'moments'``); it must be
          either ``'moments'`` or ``'riemann_sum'``

        - ``order`` -- integer (default: 1); the order of the
          derivative to be computed

        OUTPUT:

        A function from the `p`-adic upper half plane to `\CC_p`. If
        an argument ``z`` was passed, returns instead the value of
        the derivative at that point.

        EXAMPLES:

        Integrating the Poisson kernel against a measure yields a
        value of the associated modular form. Such values can be
        computed efficiently using the overconvergent method, as long
        as one starts with an ordinary form::

            sage: X = BruhatTitsQuotient(7, 2)
            sage: X.genus()
            1

        Since the genus is 1, the space of weight 2 forms is 1
        dimensional. Hence any nonzero form will be a `U_7`
        eigenvector. By Jacquet-Langlands and Cerednik-Drinfeld, in
        this case the Hecke eigenvalues correspond to that of any
        nonzero form on `\Gamma_0(14)` of weight `2`. Such a form is
        ordinary at `7`, and so we can apply the overconvergent method
        directly to this form without `p`-stabilizing::

            sage: H = X.harmonic_cocycles(2,prec=5)
            sage: h = H.gen(0)
            sage: A = X.padic_automorphic_forms(2,prec=5,overconvergent=True)
            sage: f0 = A.lift(h)

        Now that we've lifted our harmonic cocycle to an
        overconvergent automorphic form, we extract the associated
        modular form as a function and test the modular property::

            sage: T.<x> = Qq(49,prec=10)
            sage: f = f0.modular_form()
            sage: g = X.get_embedding_matrix()*X.get_units_of_order()[1]
            sage: a,b,c,d = g.change_ring(T).list()
            sage: (c*x +d)^2*f(x)-f((a*x + b)/(c*x + d))
            O(7^5)

        We can also compute the Shimura-Maass derivative, which is a
        nearly rigid analytic modular forms of weight 4::

            sage: f = f0.derivative()
            sage: (c*x + d)^4*f(x)-f((a*x + b)/(c*x + d))
            O(7^5)
        """
        def F(z, level=level, method=method):
            R = PolynomialRing(z.parent(), 'x,y').fraction_field()
            Rx = PolynomialRing(z.parent(), 'x1').fraction_field()
            x1 = Rx.gen()
            subst = R.hom([x1, z], codomain=Rx)
            x, y = R.gens()
            center = self.parent()._source._BT.find_containing_affinoid(z)
            zbar = z.trace() - z
            f = R(1) / (x - y)
            k = self.parent()._n + 2
            V = [f]
            for ii in range(order):
                V = [v.derivative(y) for v in V] + [k / (y - zbar) * v
                                                    for v in V]
                k += 2
            return sum(self.integrate(subst(v), center, level, method)
                       for v in V)
        if z is None:
            return F

        return F(z, level, method)

    # So far we cannot break it into two integrals because of the pole
    # at infinity.
    def coleman(self, t1, t2, E=None, method='moments', mult=False):
        r"""
        If ``self`` is a `p`-adic automorphic form that
        corresponds to a rigid modular form, then this computes the
        Coleman integral of this form between two points on the
        boundary `P^1(\QQ_p)` of the `p`-adic upper half plane.

        INPUT:

        - ``t1``, ``t2`` -- elements of `P^1(\QQ_p)` (the endpoints
          of integration)

        - ``E`` -- (default: ``None``) if specified, will not compute the
          covering adapted to ``t1`` and ``t2`` and instead use the
          given one. In that case, ``E`` should be a list of matrices
          corresponding to edges describing the open balls to be
          considered.

        - ``method`` -- string (default: ``'moments'``); tells which
          algorithm to use (alternative is ``'riemann_sum'``, which is
          unsuitable for computations requiring high precision)

        - ``mult`` -- boolean (default: ``False``); whether to compute the
          multiplicative version

        OUTPUT: the result of the Coleman integral

        EXAMPLES::

            sage: p = 7
            sage: lev = 2
            sage: prec = 10
            sage: X = BruhatTitsQuotient(p, lev)
            sage: k = 2
            sage: M = X.harmonic_cocycles(k, prec)
            sage: B = M.basis()
            sage: f = 3*B[0]
            sage: MM = X.padic_automorphic_forms(k, prec, overconvergent=True)
            sage: D = -11
            sage: X.is_admissible(D)
            True
            sage: K.<a> = QuadraticField(D)
            sage: Kp.<g> = Qq(p**2, prec)
            sage: P = Kp.gen()
            sage: Q = 2 + Kp.gen() + p*(Kp.gen()+1)
            sage: F = MM.lift(f)  # long time
            sage: J0 = F.coleman(P, Q, mult=True)  # long time

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        p = self.parent().prime()
        K = t1.parent()
        R = PolynomialRing(K, 'x')
        x = R.gen()
        R1 = LaurentSeriesRing(K, 'r1', default_prec=self.parent()._U.base_ring().precision_cap())
        if E is None:
            E = self.parent()._source._BT.find_covering(t1, t2)
        value = 0
        value_exp = K(1)
        if method == 'riemann_sum':
            for e in E:
                b = e[0, 1]
                d = e[1, 1]
                y = (b - d * t1) / (b - d * t2)
                poly = R1(y.log())
                c_e = self.evaluate(e)
                new = eval_dist_at_powseries(c_e, poly)
                value += new
                if mult:
                    value_exp *= K.teichmuller(y) ** Integer(c_e.moment(0).rational_reconstruction())

        elif method == 'moments':
            for e in E:
                f = (x - t1) / (x - t2)
                a, b, c, d = e.list()
                y0 = f(R1([b, a]) / R1([d, c]))  # f( (ax+b)/(cx+d) )
                y0 = p ** (-y0(ZZ(0)).valuation()) * y0
                mu = K.teichmuller(y0(ZZ(0)))
                y = y0 / mu - 1
                poly = R1(0)
                ypow = y
                for jj in range(1, R1.default_prec() + 10):
                    poly += (-1) ** (jj + 1) * ypow / jj
                    ypow *= y
                c_e = self.evaluate(e)
                new = eval_dist_at_powseries(c_e, poly)
                if hasattr(new, 'degree'):
                    assert 0
                value += new
                if mult:
                    value_exp *= K.teichmuller((b - d * t1) / (b - d * t2)) ** Integer(c_e.moment(0).rational_reconstruction())

        else:
            print('The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should not be used in practice.')
            return False
        if mult:
            return K.teichmuller(value_exp) * value.exp()
        return value


class pAdicAutomorphicForms(Module, UniqueRepresentation):
    Element = pAdicAutomorphicFormElement

    @staticmethod
    def __classcall__(cls, domain, U, prec=None, t=None, R=None,
                      overconvergent=False):
        r"""
        The module of (quaternionic) `p`-adic automorphic forms.

        INPUT:

        - ``domain`` -- a BruhatTitsQuotient

        - ``U`` -- a distributions module or an integer. If ``U`` is a
          distributions module then this creates the relevant space of
          automorphic forms. If ``U`` is an integer then the coefficients
          are the (`U-2`)nd power of the symmetric representation of
          `GL_2(\QQ_p)`.

        - ``prec`` -- a precision (default: ``None``); if not ``None`` should
          be a positive integer

        - ``t`` -- (default: ``None``) the number of additional moments to
          store. If ``None``, determine it automatically from ``prec``, ``U``
          and the ``overconvergent`` flag.

        - ``R`` -- (default: ``None``) if specified, coefficient field of the automorphic forms.
          If not specified it defaults to the base ring of the distributions ``U``, or to `Q_p`
          with the working precision ``prec``.

        - ``overconvergent`` -- boolean (default: ``False``); if ``True``, will construct overconvergent
          `p`-adic automorphic forms. Otherwise it constructs the finite dimensional space of
          `p`-adic automorphic forms which is isomorphic to the space of harmonic cocycles.

        EXAMPLES:

        The space of weight 2 p-automorphic forms is isomorphic with
        the space of scalar valued invariant harmonic cocycles::

            sage: X = BruhatTitsQuotient(11,5)
            sage: H0 = X.padic_automorphic_forms(2,10)
            sage: H1 = X.padic_automorphic_forms(2,prec = 10)
            sage: H0 == H1
            True

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        return super().__classcall__(cls, domain, U,
                                     prec, t, R,
                                     overconvergent)

    def __init__(self, domain, U, prec=None, t=None, R=None,
                 overconvergent=False):
        """
        Create a space of `p`-automorphic forms.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,5)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: TestSuite(A).run()
        """
        if R is None:
            if not isinstance(U, Integer):
                self._R = U.base_ring()
            else:
                if prec is None:
                    prec = 100
                self._R = Qp(domain._p, prec)
        else:
            self._R = R
        # U is a CoefficientModuleSpace
        if isinstance(U, Integer):
            if t is None:
                if overconvergent:
                    t = prec - U + 1
                else:
                    t = 0
            if overconvergent:
                self._U = OverconvergentDistributions(U - 2, base=self._R,
                                                      prec_cap=U - 1 + t,
                                                      act_on_left=True,
                                                      adjuster=_btquot_adjuster(),
                                                      dettwist=-ZZ((U - 2) // 2),
                                                      act_padic=True)
            else:
                self._U = Symk(U - 2, base=self._R, act_on_left=True,
                               adjuster=_btquot_adjuster(),
                               dettwist=-ZZ((U - 2) // 2),
                               act_padic=True)
        else:
            self._U = U
        self._source = domain
        self._list = self._source.get_list()  # Contains also the opposite edges
        self._prec = self._R.precision_cap()
        self._n = self._U.weight()
        self._p = self._source._p

        self._Sigma0 = self._U._act._Sigma0

        Module.__init__(self, base=self._R)
        self._populate_coercion_lists_()

    def prime(self):
        """
        Return the underlying prime.

        OUTPUT: ``p`` -- prime integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,5)
            sage: H = X.harmonic_cocycles(2,prec = 10)
            sage: A = X.padic_automorphic_forms(2,prec = 10)
            sage: A.prime()
            11
        """
        return self._p

    def zero(self):
        r"""
        Return the zero element of ``self``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5, 7)
            sage: H1 = X.padic_automorphic_forms( 2, prec=10)
            sage: H1.zero() == 0
            True
        """
        return self.element_class(self, [self._U(0) for o in self._list])

    def __eq__(self, other):
        r"""
        Test whether two pAdicAutomorphicForm spaces are equal.

        INPUT:

        - ``other`` -- another space of `p`-automorphic forms

        OUTPUT: boolean

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: H1 = X.padic_automorphic_forms(2,prec = 10)
            sage: H2 = X.padic_automorphic_forms(2,prec = 10)
            sage: H1 == H2
            True
        """
        if not isinstance(other, pAdicAutomorphicForms):
            return False

        return (self.base_ring() == other.base_ring() and
                self._source == other._source and
                self._U == other._U)

    def __ne__(self, other):
        r"""
        Test whether two pAdicAutomorphicForm spaces are not equal.

        INPUT:

        - ``other`` -- another space of `p`-automorphic forms

        OUTPUT: boolean

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: H1 = X.padic_automorphic_forms(2,prec = 10)
            sage: H2 = X.padic_automorphic_forms(2,prec = 10)
            sage: H1 == H2
            True
        """
        return not self.__eq__(other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: H1 = X.padic_automorphic_forms(2,prec = 10)
            sage: H2 = X.padic_automorphic_forms(2,prec = 10)
            sage: hash(H1) == hash(H2)
            True
        """
        return hash((self.base_ring(), self._source, self._U))

    def _repr_(self):
        r"""
        Return the representation of ``self`` as a string.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: A = X.padic_automorphic_forms(2,prec = 10)
            sage: A # indirect doctest
            Space of automorphic forms on Quotient of the Bruhat Tits tree of GL_2(QQ_3) with discriminant 7 and level 1 with values in Sym^0 Q_3^2
        """
        s = 'Space of automorphic forms on '
        s += str(self._source)
        s += ' with values in ' + str(self._U)
        return s

    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other BruhatTitsHarmonicCocycles or from pAdicAutomorphicForms.

        INPUT:

        - ``S`` -- a BruhatTitsHarmonicCocycle or pAdicAutomorphicForm

        OUTPUT: boolean; ``True`` if and only if `S` is coercible into ``self``

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: A._coerce_map_from_(H)
            True
        """
        if isinstance(S, BruhatTitsHarmonicCocycles):
            if S.weight() - 2 != self._n:
                return False
            if S._X != self._source:
                return False
            return True
        if isinstance(S, pAdicAutomorphicForms):
            if S._n != self._n:
                return False
            if S._source != self._source:
                return False
            return True
        return False

    def _element_constructor_(self, data):
        r"""
        Construct a `p`-automorphic form.

        INPUT:

        - ``data`` -- defining data. Can be either a harmonic cocycle, or a `p`-adic automorphic form,
          or a list of elements coercible into the module of coefficients of ``self``

        OUTPUT: a `p`-adic automorphic form

        EXAMPLES::

            sage: X = BruhatTitsQuotient(13,5)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: h = H.an_element()  # indirect doctest
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: A(h)
            p-adic automorphic form of cohomological weight 0
        """
        # Code how to coerce x into the space
        # Admissible values of x?
        if isinstance(data, list):
            return self.element_class(self, [self._U(o, normalize=False) for o in data])

        if isinstance(data, pAdicAutomorphicFormElement):
            vals = [self._U(o, normalize=False) for o in data._value]
            return self.element_class(self, vals)

        if isinstance(data, BruhatTitsHarmonicCocycleElement):
            E = self._list
            tmp = []
            F = []
            Uold = data.parent()._U
            for ii in range(len(data._F)):
                newtmp = data.parent()._Sigma0(E[ii].rep.inverse(), check=False) * Uold(data._F[ii],
                                                                                        normalize=False)
                tmp.append(newtmp)
                F.append(newtmp)
            A = data.parent()._Sigma0(Matrix(QQ, 2, 2, [0, ~self.prime(), 1, 0]), check=False)
            F.extend(-(A * tmp[ii]) for ii in range(len(data._F)))
            vals = self._make_invariant([self._U(o, normalize=False) for o in F])
            return self.element_class(self, vals)
        if data == 0:
            return self.zero()

    def _an_element_(self):
        r"""
        Return an element of the module.

        OUTPUT: a harmonic cocycle

        EXAMPLES::

            sage: X = BruhatTitsQuotient(13,5)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: A.an_element() # indirect doctest
            p-adic automorphic form of cohomological weight 0
        """
        return self(0)

    def precision_cap(self):
        """
        Return the precision of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(13,11)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: A.precision_cap()
            10
        """
        return self._prec

    def lift(self, f):
        r"""
        Lift the harmonic cocycle ``f`` to a p-automorphic form.

        If one is using overconvergent coefficients, then this will
        compute all of the moments of the measure associated to ``f``.

        INPUT:

        - ``f`` -- a harmonic cocycle

        OUTPUT: a `p`-adic automorphic form

        EXAMPLES:

        If one does not work with an overconvergent form then lift
        does nothing::

            sage: X = BruhatTitsQuotient(13,5)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: h = H.gen(0)
            sage: A = X.padic_automorphic_forms(2,prec=10)
            sage: A.lift(h) # long time
            p-adic automorphic form of cohomological weight 0

        With overconvergent forms, the input is lifted naively and its
        moments are computed::

            sage: X = BruhatTitsQuotient(13,11)
            sage: H = X.harmonic_cocycles(2,prec=5)
            sage: A2 = X.padic_automorphic_forms(2,prec=5,overconvergent=True)
            sage: a = H.gen(0)
            sage: A2.lift(a) # long time
            p-adic automorphic form of cohomological weight 0
        """
        return self(f)._improve(f)

    def _make_invariant(self, F):
        r"""
        Naively lift a ``classical`` automorphic form to an
        overconvergent form.

        INPUT:

        - ``F`` -- a classical (nonoverconvergent) pAdicAutomorphicForm or
          BruhatTitsHarmonicCocycle

        OUTPUT: an overconvergent pAdicAutomorphicForm

        EXAMPLES::

            sage: X = BruhatTitsQuotient(13,11)
            sage: H = X.harmonic_cocycles(2,prec = 5)
            sage: A = X.padic_automorphic_forms(2,prec = 5)
            sage: h = H.basis()[0]
            sage: A.lift(h) # indirect doctest long time
            p-adic automorphic form of cohomological weight 0
        """
        S = self._source.get_stabilizers()
        M = [e.rep for e in self._list]
        newF = []
        for ii in range(len(S)):
            Si = S[ii]
            x = self._U(F[ii], normalize=False)

            if any(v[2] for v in Si):
                newFi = self._U(0)
                s = QQ(0)
                m = M[ii]
                for v in Si:
                    s += 1
                    g = self._Sigma0(m.adjugate() * self._source.embed_quaternion(v[0], prec=self._prec).adjugate() * m,
                                     check=False)
                    newFi += g * x
                newF.append((QQ(1) / s) * newFi)
            else:
                newF.append(self._U(x, normalize=False))
        return newF

    def _apply_Up_operator(self, f, scale=False, original_moments=None):
        r"""
        Apply the Up operator to ``f``.

        INPUT:

        - ``f`` -- a `p`-adic automorphic form
        - ``scale`` -- boolean (default: ``True``); whether to scale by the appropriate power of `p`
          at each iteration

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,11)
            sage: M = X.harmonic_cocycles(4,10)
            sage: A = X.padic_automorphic_forms(4,10, overconvergent = True)
            sage: F = A.lift(M.basis()[0]); F # indirect doctest
            p-adic automorphic form of cohomological weight 2
        """
        HeckeData = self._source._get_Up_data()
        S0 = f._value[0].parent()._act._Sigma0
        prec_cap = self._U.base_ring().precision_cap()

        if not scale:
            factor = self._p ** (self._U.weight() // 2)
        else:
            factor = 1

        # Save original moments
        if original_moments is None:
            original_moments = [[fval._moments[ii] for ii in range(self._n + 1)]
                                for fval in f._value]

        Tf = []
        for jj in range(len(self._list)):
            tmp = self._U(0, normalize=False)
            for gg, edge_list in HeckeData:
                u = edge_list[jj]
                tprec = 2 * (prec_cap + u.power) + 1
                r = S0(self._p ** -u.power * (u.t(tprec) * gg).adjugate(),
                       check=False)
                tmp += r * f._value[u.label]
            tmp *= factor
            for ii in range(self._n + 1):
                tmp._moments[ii] = original_moments[jj][ii]
            Tf.append(tmp)
        return self(Tf)
