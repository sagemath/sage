r"""
Nilpotent Lie groups

AUTHORS:

- Eero Hakavuori (2018-09-25): initial version of nilpotent Lie groups
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.categories.lie_groups import LieGroups
from sage.categories.lie_algebras import LieAlgebras
from sage.groups.group import Group
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.structure import(DifferentialStructure,
                                     RealDifferentialStructure)
from sage.misc.cachefunc import cached_method
from sage.misc.repr import repr_lincomb
from sage.modules.free_module_element import vector
import sage.rings.abc
from sage.structure.element import MultiplicativeGroupElement
from sage.symbolic.ring import SR, var


def _symbolic_lie_algebra_copy(L):
    r"""
    Create a copy of the Lie algebra ``L`` admitting symbolic coefficients.

    This is used internally to compute a symbolic expression for the group law.

    INPUT:

    - ``L`` -- a finite dimensional Lie algebra with basis

    EXAMPLES::

        sage: from sage.groups.lie_gps.nilpotent_lie_group import _symbolic_lie_algebra_copy
        sage: L = LieAlgebra(QQ, 2, step=2)
        sage: L_SR = _symbolic_lie_algebra_copy(L)
        sage: L.structure_coefficients()
        Finite family {((1,), (2,)): X_12}
        sage: L_SR.structure_coefficients()
        Finite family {((1,), (2,)): L[(1, 2)]}

    TESTS:

    Verify that copying works with something that is not an instance of
    :class:`LieAlgebraWithStructureCoefficients`::

        sage: from sage.groups.lie_gps.nilpotent_lie_group import _symbolic_lie_algebra_copy
        sage: L = lie_algebras.Heisenberg(QQ, 1)
        sage: hasattr(L, 'change_ring')
        False
        sage: L_SR = _symbolic_lie_algebra_copy(L)
        sage: L_SR.structure_coefficients()
        Finite family {('p1', 'q1'): z}
    """
    try:
        return L.change_ring(SR)
    except AttributeError:
        s_coeff = L.structure_coefficients()
        index_set = L.basis().keys()
        names = L.variable_names()
        return LieAlgebraWithStructureCoefficients(SR, s_coeff, names=names,
                                                   index_set=index_set)


class NilpotentLieGroup(Group, DifferentiableManifold):
    r"""
    A nilpotent Lie group.

    INPUT:

    - ``L`` -- the Lie algebra of the Lie group; must be a finite
      dimensional nilpotent Lie algebra with basis over a topological
      field, e.g. `\QQ` or `\RR`
    - ``name`` -- a string; name (symbol) given to the Lie group

    Two types of exponential coordinates are defined on any
    nilpotent Lie group using the basis of the Lie algebra,
    see :meth:`chart_exp1` and :meth:`chart_exp2`.

    EXAMPLES:

    Creation of a nilpotent Lie group::

        sage: L = lie_algebras.Heisenberg(QQ, 1)
        sage: G = L.lie_group(); G
        Lie group G of Heisenberg algebra of rank 1 over Rational Field

    Giving a different name to the group::

        sage: L.lie_group('H')
        Lie group H of Heisenberg algebra of rank 1 over Rational Field

    Elements can be created using the exponential map::

        sage: p,q,z = L.basis()
        sage: g = G.exp(p); g
        exp(p1)
        sage: h = G.exp(q); h
        exp(q1)

    Lie group multiplication has the usual product syntax::

        sage: k = g*h; k
        exp(p1 + q1 + 1/2*z)

    The identity element is given by :meth:`one`::

        sage: e = G.one(); e
        exp(0)
        sage: e*k == k and k*e == k
        True

    The default coordinate system is exponential coordinates of the first kind::

        sage: G.default_chart() == G.chart_exp1()
        True
        sage: G.chart_exp1()
        Chart (G, (x_0, x_1, x_2))

    Changing the default coordinates to exponential coordinates of the second
    kind will change how elements are printed::

        sage: G.set_default_chart(G.chart_exp2())
        sage: k
        exp(z)exp(q1)exp(p1)
        sage: G.set_default_chart(G.chart_exp1())
        sage: k
        exp(p1 + q1 + 1/2*z)

    The frames of left- or right-invariant vector fields are created using
    :meth:`left_invariant_frame` and :meth:`right_invariant_frame`::

        sage: X = G.left_invariant_frame(); X
        Vector frame (G, (X_0,X_1,X_2))
        sage: X[0]
        Vector field X_0 on the Lie group G of Heisenberg algebra of rank 1 over Rational Field

    A vector field can be displayed with respect to a coordinate frame::

        sage: exp1_frame = G.chart_exp1().frame()
        sage: exp2_frame = G.chart_exp2().frame()
        sage: X[0].display(exp1_frame)
        X_0 = ∂/∂x_0 - 1/2*x_1 ∂/∂x_2
        sage: X[0].display(exp2_frame)
        X_0 = ∂/∂y_0
        sage: X[1].display(exp1_frame)
        X_1 = ∂/∂x_1 + 1/2*x_0 ∂/∂x_2
        sage: X[1].display(exp2_frame)
        X_1 = ∂/∂y_1 + x_0 ∂/∂y_2

    Defining a left translation by a generic point::

        sage: g = G.point([var('a'), var('b'), var('c')]); g
        exp(a*p1 + b*q1 + c*z)
        sage: L_g = G.left_translation(g); L_g
        Diffeomorphism of the Lie group G of Heisenberg algebra of rank 1 over Rational Field
        sage: L_g.display()
        G → G
           (x_0, x_1, x_2) ↦ (a + x_0, b + x_1, -1/2*b*x_0 + 1/2*a*x_1 + c + x_2)
           (x_0, x_1, x_2) ↦ (y_0, y_1, y_2) = (a + x_0, b + x_1,
                                                1/2*a*b + 1/2*(2*a + x_0)*x_1 + c + x_2)
           (y_0, y_1, y_2) ↦ (x_0, x_1, x_2) = (a + y_0, b + y_1,
                                                -1/2*b*y_0 + 1/2*(a - y_0)*y_1 + c + y_2)
           (y_0, y_1, y_2) ↦ (a + y_0, b + y_1, 1/2*a*b + a*y_1 + c + y_2)

    Verifying the left-invariance of the left-invariant frame::

        sage: x = G(G.chart_exp1()[:])
        sage: L_g.differential(x)(X[0].at(x)) == X[0].at(L_g(x))
        True
        sage: L_g.differential(x)(X[1].at(x)) == X[1].at(L_g(x))
        True
        sage: L_g.differential(x)(X[2].at(x)) == X[2].at(L_g(x))
        True

    An element of the Lie algebra can be extended to a left or right invariant
    vector field::

        sage: X_L = G.left_invariant_extension(p + 3*q); X_L
        Vector field p1 + 3*q1 on the Lie group G of Heisenberg algebra of rank 1 over Rational Field
        sage: X_L.display(exp1_frame)
        p1 + 3*q1 = ∂/∂x_0 + 3 ∂/∂x_1 + (3/2*x_0 - 1/2*x_1) ∂/∂x_2
        sage: X_R = G.right_invariant_extension(p + 3*q)
        sage: X_R.display(exp1_frame)
        p1 + 3*q1 = ∂/∂x_0 + 3 ∂/∂x_1 + (-3/2*x_0 + 1/2*x_1) ∂/∂x_2

    The nilpotency step of the Lie group is the nilpotency step of its algebra.
    Nilpotency for Lie groups means that group commutators that are longer than
    the nilpotency step vanish::

        sage: G.step()
        2
        sage: g = G.exp(p); h = G.exp(q)
        sage: c = g*h*~g*~h; c
        exp(z)
        sage: g*c*~g*~c
        exp(0)
    """

    def __init__(self, L, name, **kwds):
        r"""
        Initialize ``self``.

        TESTS::

            sage: L = lie_algebras.Heisenberg(QQ, 2)
            sage: G = L.lie_group()
            sage: TestSuite(G).run()
        """
        required_cat = LieAlgebras(L.base_ring()).FiniteDimensional()
        required_cat = required_cat.WithBasis().Nilpotent()
        if L not in required_cat:
            raise TypeError("L needs to be a finite dimensional nilpotent "
                            "Lie algebra with basis")
        self._lie_algebra = L

        R = L.base_ring()
        category = kwds.pop('category', None)
        category = LieGroups(R).or_subcategory(category)
        if isinstance(R, sage.rings.abc.RealField):
            structure = RealDifferentialStructure()
        else:
            structure = DifferentialStructure()

        DifferentiableManifold.__init__(self, L.dimension(), name, R,
                                        structure, category=category)

        # initialize exponential coordinates of the first kind
        basis_strs = [str(X) for X in L.basis()]
        split = list(zip(*[s.split('_') for s in basis_strs]))
        if len(split) == 2 and all(sk == split[0][0] for sk in split[0]):
            self._var_indexing = split[1]
        else:
            self._var_indexing = [str(k) for k in range(L.dimension())]
        variables = ' '.join('x_%s' % k for k in self._var_indexing)
        self._Exp1 = self.chart(variables)

        # compute a symbolic formula for the group law
        L_SR = _symbolic_lie_algebra_copy(L)
        n = L.dimension()
        a, b = (tuple(var('%s_%d' % (s, j)) for j in range(n))
                for s in ['a', 'b'])
        self._group_law_vars = (a, b)
        bch = L_SR.bch(L_SR.from_vector(a), L_SR.from_vector(b), L.step())
        self._group_law = vector(SR, (zk.expand() for zk in bch.to_vector()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(RR, 1)
            sage: L.lie_group()
            Lie group G of Heisenberg algebra of rank 1 over
             Real Field with 53 bits of precision
        """
        return "Lie group %s of %s" % (self._name, self.lie_algebra())

    def _dLx(self):
        r"""
        Return the matrix of the differential at the identity of a left
        translation by a generic point in the default coordinate system.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: G._dLx()
            [       1       0  0]
            [       0       1  0]
            [-1/2*x_2 1/2*x_1  1]
        """
        a = self._group_law_vars[0]
        x = self.default_chart()[:]
        asubs = dict(zip(a, x))
        L_a = self.left_translation(self.point(a))
        return L_a.differential(self.one()).matrix().subs(asubs)

    def _dRx(self):
        r"""
        Return the matrix of the differential at the identity of a right
        translation by a generic point in the default coordinate system.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: G._dRx()
            [       1        0        0]
            [       0        1        0]
            [ 1/2*x_2 -1/2*x_1        1]
        """
        a = self._group_law_vars[0]
        x = self.default_chart()[:]
        asubs = dict(zip(a, x))
        R_a = self.right_translation(self.point(a))
        return R_a.differential(self.one()).matrix().subs(asubs)

    @cached_method
    def gens(self):
        r"""
        Return a tuple of elements whose one-parameter subgroups generate
        the Lie group.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1)
            sage: G = L.lie_group()
            sage: G.gens()
            (exp(p1), exp(q1), exp(z))
        """
        return tuple(self.exp(X) for X in self.lie_algebra().basis())

    def lie_algebra(self):
        r"""
        Return the Lie algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: G.lie_algebra() == L
            True
        """
        return self._lie_algebra

    def step(self):
        r"""
        Return the nilpotency step of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=4)
            sage: G = L.lie_group()
            sage: G.step()
            4
        """
        return self._lie_algebra.step()

    def chart_exp1(self):
        r"""
        Return the chart of exponential coordinates of the first kind.

        Exponential coordinates of the first kind are

        .. MATH ::

            \exp(x_1X_1 + \cdots + x_nX_n) \mapsto (x_1, \ldots, x_n).

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: G.chart_exp1()
            Chart (G, (x_1, x_2, x_12))
        """
        return self._Exp1

    @cached_method
    def chart_exp2(self):
        r"""
        Return the chart of exponential coordinates of the second kind.

        Exponential coordinates of the second kind are

        .. MATH ::

            \exp(x_nX_n) \cdots \exp(x_1X_1) \mapsto (x_1, \ldots, x_n).

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: G.chart_exp2()
            Chart (G, (y_1, y_2, y_12))
        """
        variables = ' '.join('y_%s' % k for k in self._var_indexing)
        ret = self.chart(variables)

        # compute transitions between exponential coordinates
        # compute exp-2 to exp-1
        n = self.dimension()
        Z = self.one()
        for k, yk in enumerate(ret[:]):
            v = [0] * n
            v[k] = yk
            Z = self.point(v, chart=self._Exp1) * Z
        f = [zk.expand() for zk in Z.coordinates(chart=self._Exp1)]
        ret.transition_map(self._Exp1, f)

        # compute exp-1 to exp-2 by inverting the previous map
        inv_subs = {}
        for xk, yk, fk in zip(self._Exp1[:], ret[:], f):
            inv_subs[yk] = xk - (fk - yk).subs(inv_subs)
        f_inv = [inv_subs[yk].expand() for yk in ret[:]]
        self._Exp1.transition_map(ret, f_inv)
        return ret

    def exp(self, X):
        r"""
        Return the group element `exp(X)`.

        INPUT:

        - ``X`` -- an element of the Lie algebra of ``self``

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: G.exp(X)
            exp(X)
            sage: G.exp(Y)
            exp(Y)
            sage: G.exp(X + Y)
            exp(X + Y)
        """
        return self.point(X.to_vector(), chart=self._Exp1)

    def log(self, x):
        r"""
        Return the logarithm of the element ``x`` of ``self``.

        INPUT:

        - ``x`` -- an element of ``self``

        The logarithm is by definition the inverse of :meth:`exp`.

        If the Lie algebra of ``self`` does not admit symbolic coefficients,
        the logarithm is not defined for abstract, i.e. symbolic, points.

        EXAMPLES:

        The logarithm is the inverse of the exponential::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: G.log(G.exp(X)) == X
            True
            sage: G.log(G.exp(X)*G.exp(Y))
            X + Y + 1/2*Z

        The logarithm is not defined for abstract (symbolic) points::

            sage: g = G.point([var('a'), 1, 2]); g
            exp(a*X + Y + 2*Z)
            sage: G.log(g)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert a to a rational
        """
        xvec = x.coordinates(chart=self._Exp1)
        return self.lie_algebra().from_vector(xvec)

    def one(self):
        r"""
        Return the identity element of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 2, step=4)
            sage: G = L.lie_group()
            sage: G.one()
            exp(0)
        """
        return self.exp(self._lie_algebra.zero())

    def left_translation(self, g):
        r"""
        Return the left translation by ``g`` as an automorphism of ``self``.

        The left translation by `g` on a Lie group `G` is the map

        .. MATH::

            G \to G, \qquad
            h \mapsto gh.

        INPUT:

        - ``g`` -- an element of ``self``

        EXAMPLES:

        A left translation in the Heisenberg group::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: p,q,z = H.basis()
            sage: G = H.lie_group()
            sage: g = G.exp(p)
            sage: L_g = G.left_translation(g); L_g
            Diffeomorphism of the Lie group G of Heisenberg algebra of rank 1 over Rational Field
            sage: L_g.display(chart1=G.chart_exp1(), chart2=G.chart_exp1())
            G → G
                (x_0, x_1, x_2) ↦ (x_0 + 1, x_1, 1/2*x_1 + x_2)

        Left translation by a generic element::

            sage: h = G.point([var('a'), var('b'), var('c')])
            sage: L_h = G.left_translation(h)
            sage: L_h.display(chart1=G.chart_exp1(), chart2=G.chart_exp1())
            G → G
               (x_0, x_1, x_2) ↦ (a + x_0, b + x_1, -1/2*b*x_0 + 1/2*a*x_1 + c + x_2)
        """
        chart = self.default_chart()
        x = self.point(chart[:])
        L_g_expr = (g * x).coordinates()
        return self.diffeomorphism(self, coord_functions=L_g_expr)

    def left_invariant_frame(self, **kwds):
        r"""
        Return the frame of left-invariant vector fields of ``self``.

        The labeling of the frame and the dual frame can be customized using
        keyword parameters as described in
        :meth:`sage.manifolds.differentiable.manifold.DifferentiableManifold.vector_frame`.

        EXAMPLES:

        The default left-invariant frame::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: livf = G.left_invariant_frame(); livf
            Vector frame (G, (X_1,X_2,X_12))
            sage: coord_frame = G.chart_exp1().frame()
            sage: livf[0].display(coord_frame)
            X_1 = ∂/∂x_1 - 1/2*x_2 ∂/∂x_12
            sage: livf[1].display(coord_frame)
            X_2 = ∂/∂x_2 + 1/2*x_1 ∂/∂x_12
            sage: livf[2].display(coord_frame)
            X_12 = ∂/∂x_12

        Examples of custom labeling for the frame::

            sage: G.left_invariant_frame(symbol='Y')
            Vector frame (G, (Y_1,Y_2,Y_12))
            sage: G.left_invariant_frame(symbol='Z', indices=None)
            Vector frame (G, (Z_0,Z_1,Z_2))
            sage: G.left_invariant_frame(symbol='W', indices=('a','b','c'))
            Vector frame (G, (W_a,W_b,W_c))
        """
        dLx_field = self.automorphism_field()
        dLx_field[:] = self._dLx()
        coord_frame = self._Exp1.frame()
        symbol = kwds.pop('symbol', 'X')
        indices = kwds.pop('indices', self._var_indexing)
        return coord_frame.new_frame(dLx_field, symbol=symbol,
                                     indices=indices, **kwds)

    livf = left_invariant_frame

    def left_invariant_extension(self, X, name=None):
        r"""
        Return the left-invariant vector field that has the value ``X``
        at the identity.

        INPUT:

        - ``X`` -- an element of the Lie algebra of ``self``
        - ``name`` -- (optional) a string to use as a name for the vector field;
          if nothing is given, the name of the vector ``X`` is used

        EXAMPLES:

        A left-invariant extension in the Heisenberg group::

            sage: L = lie_algebras.Heisenberg(QQ, 1)
            sage: p, q, z = L.basis()
            sage: H = L.lie_group('H')
            sage: X = H.left_invariant_extension(p); X
            Vector field p1 on the Lie group H of Heisenberg algebra of rank 1 over Rational Field
            sage: X.display(H.chart_exp1().frame())
            p1 = ∂/∂x_0 - 1/2*x_1 ∂/∂x_2

        Default vs. custom naming for the invariant vector field::

            sage: Y = H.left_invariant_extension(p + q); Y
            Vector field p1 + q1 on the Lie group H of Heisenberg algebra of rank 1 over Rational Field
            sage: Z = H.left_invariant_extension(p + q, 'Z'); Z
            Vector field Z on the Lie group H of Heisenberg algebra of rank 1 over Rational Field
        """
        if name is None:
            name = str(X)
        X_vf = self.vector_field(name)
        frame = self._Exp1.frame()
        X_vf[frame, :] = self._dLx() * X.to_vector()
        return X_vf

    def right_translation(self, g):
        r"""
        Return the right translation by ``g`` as an automorphism of ``self``.

        The right translation by `g` on a Lie group `G` is the map

        .. MATH::

            G \to G, \qquad
            h\mapsto hg.

        INPUT:

        - ``g`` -- an element of ``self``

        EXAMPLES:

        A right translation in the Heisenberg group::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: p,q,z = H.basis()
            sage: G = H.lie_group()
            sage: g = G.exp(p)
            sage: R_g = G.right_translation(g); R_g
            Diffeomorphism of the Lie group G of Heisenberg algebra of rank 1 over Rational Field
            sage: R_g.display(chart1=G.chart_exp1(), chart2=G.chart_exp1())
            G → G
               (x_0, x_1, x_2) ↦ (x_0 + 1, x_1, -1/2*x_1 + x_2)

        Right translation by a generic element::

            sage: h = G.point([var('a'), var('b'), var('c')])
            sage: R_h = G.right_translation(h)
            sage: R_h.display(chart1=G.chart_exp1(), chart2=G.chart_exp1())
            G → G
               (x_0, x_1, x_2) ↦ (a + x_0, b + x_1, 1/2*b*x_0 - 1/2*a*x_1 + c + x_2)
        """
        chart = self.default_chart()
        x = self.point(chart[:])
        R_g_expr = (x * g).coordinates()
        return self.diffeomorphism(self, coord_functions=R_g_expr)

    def right_invariant_frame(self, **kwds):
        r"""
        Return the frame of right-invariant vector fields of ``self``.

        The labeling of the frame and the dual frame can be customized using
        keyword parameters as described in
        :meth:`sage.manifolds.differentiable.manifold.DifferentiableManifold.vector_frame`.

        EXAMPLES:

        The default right-invariant frame::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: rivf = G.right_invariant_frame(); rivf
            Vector frame (G, (XR_1,XR_2,XR_12))
            sage: coord_frame = G.chart_exp1().frame()
            sage: rivf[0].display(coord_frame)
            XR_1 = ∂/∂x_1 + 1/2*x_2 ∂/∂x_12
            sage: rivf[1].display(coord_frame)
            XR_2 = ∂/∂x_2 - 1/2*x_1 ∂/∂x_12
            sage: rivf[2].display(coord_frame)
            XR_12 = ∂/∂x_12

        Examples of custom labeling for the frame::

            sage: G.right_invariant_frame(symbol='Y')
            Vector frame (G, (Y_1,Y_2,Y_12))
            sage: G.right_invariant_frame(symbol='Z', indices=None)
            Vector frame (G, (Z_0,Z_1,Z_2))
            sage: G.right_invariant_frame(symbol='W', indices=('a','b','c'))
            Vector frame (G, (W_a,W_b,W_c))
        """
        dRx_field = self.automorphism_field()
        dRx_field[:] = self._dRx()
        coord_frame = self._Exp1.frame()
        symbol = kwds.pop('symbol', 'XR')
        indices = kwds.pop('indices', self._var_indexing)
        return coord_frame.new_frame(dRx_field, symbol=symbol,
                                     indices=indices, **kwds)

    rivf = right_invariant_frame

    def right_invariant_extension(self, X, name=None):
        r"""
        Return the right-invariant vector field that has the value ``X``
        at the identity.

        INPUT:

        - ``X`` -- an element of the Lie algebra of ``self``
        - ``name`` -- (optional) a string to use as a name for the vector field;
          if nothing is given, the name of the vector ``X`` is used

        EXAMPLES:

        A right-invariant extension in the Heisenberg group::

            sage: L = lie_algebras.Heisenberg(QQ, 1)
            sage: p, q, z = L.basis()
            sage: H = L.lie_group('H')
            sage: X = H.right_invariant_extension(p); X
            Vector field p1 on the Lie group H of Heisenberg algebra of rank 1 over Rational Field
            sage: X.display(H.chart_exp1().frame())
            p1 = ∂/∂x_0 + 1/2*x_1 ∂/∂x_2

        Default vs. custom naming for the invariant vector field::

            sage: Y = H.right_invariant_extension(p + q); Y
            Vector field p1 + q1 on the Lie group H of Heisenberg algebra of rank 1 over Rational Field
            sage: Z = H.right_invariant_extension(p + q, 'Z'); Z
            Vector field Z on the Lie group H of Heisenberg algebra of rank 1 over Rational Field
        """
        if name is None:
            name = str(X)
        X_vf = self.vector_field(name)
        frame = self._Exp1.frame()
        X_vf[frame, :] = self._dRx() * X.to_vector()
        return X_vf

    def conjugation(self, g):
        r"""
        Return the conjugation by ``g`` as an automorphism of ``self``.

        The conjugation by `g` on a Lie group `G` is the map

        .. MATH::

            G \to G, \qquad
            h \mapsto ghg^{-1}.

        INPUT:

        - ``g`` -- an element of ``self``

        EXAMPLES:

        A generic conjugation in the Heisenberg group::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: p,q,z = H.basis()
            sage: G = H.lie_group()
            sage: g = G.point([var('a'), var('b'), var('c')])
            sage: C_g = G.conjugation(g); C_g
            Diffeomorphism of the Lie group G of Heisenberg algebra of rank 1 over Rational Field
            sage: C_g.display(chart1=G.chart_exp1(), chart2=G.chart_exp1())
            G → G
               (x_0, x_1, x_2) ↦ (x_0, x_1, -b*x_0 + a*x_1 + x_2)
        """
        chart = self.default_chart()
        x = self.point(chart[:])
        C_g_expr = (g * x * ~g).coordinates()
        return self.diffeomorphism(self, coord_functions=C_g_expr)

    def adjoint(self, g):
        r"""
        Return the adjoint map as an automorphism
        of the Lie algebra of ``self``.

        INPUT:

        - ``g`` -- an element of ``self``

        For a Lie group element `g`, the adjoint map `\operatorname{Ad}_g` is
        the map on the Lie algebra `\mathfrak{g}` given by the differential
        of the conjugation by `g` at the identity.

        If the Lie algebra of ``self`` does not admit symbolic coefficients,
        the adjoint is not in general defined for abstract points.

        EXAMPLES:

        An example of an adjoint map::

            sage: L = LieAlgebra(QQ, 2, step=3)
            sage: G = L.lie_group()
            sage: g = G.exp(L.basis().list()[0]); g
            exp(X_1)
            sage: Ad_g = G.adjoint(g); Ad_g
            Lie algebra endomorphism of Free Nilpotent Lie algebra on 5
            generators (X_1, X_2, X_12, X_112, X_122) over Rational Field
              Defn: X_1 |--> X_1
                    X_2 |--> X_2 + X_12 + 1/2*X_112
                    X_12 |--> X_12 + X_112
                    X_112 |--> X_112
                    X_122 |--> X_122

        Usually the adjoint map of a symbolic point is not defined::

            sage: L = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: g = G.point([var('a'), var('b'), var('c')]); g
            exp(a*X_1 + b*X_2 + c*X_12)
            sage: G.adjoint(g)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert -b to a rational

        However, if the adjoint map is independent from the symbolic terms,
        the map is still well defined::

            sage: g = G.point([0, 0, var('a')]); g
            exp(a*X_12)
            sage: G.adjoint(g)
            Lie algebra endomorphism of Free Nilpotent Lie algebra on 3 generators (X_1, X_2, X_12) over Rational Field
              Defn: X_1 |--> X_1
                    X_2 |--> X_2
                    X_12 |--> X_12
        """
        Adg_mat = self.conjugation(g).differential(self.one()).matrix()
        L = self.lie_algebra()
        basis_images = {X: L.from_vector(Adg_mat * X.to_vector())
                        for X in L.basis()}
        return L.morphism(basis_images, codomain=L)

    class Element(DifferentiableManifold.Element, MultiplicativeGroupElement):
        r"""
        A base class for an element of a Lie group.

        EXAMPLES:

        Elements of the  group are printed in the default
        exponential coordinates::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, 2, step=2)
            sage: G = L.lie_group()
            sage: g = G.exp(2*X + 3*Z); g
            exp(2*X + 3*Z)
            sage: h = G.point([ var('a'), var('b'), 0]); h
            exp(a*X + b*Y)
            sage: G.set_default_chart(G.chart_exp2())
            sage: g
            exp(3*Z)exp(2*X)
            sage: h
            exp(1/2*a*b*Z)exp(b*Y)exp(a*X)

        Multiplication of two elements uses the usual product syntax::

            sage: G.exp(Y)*G.exp(X)
            exp(Y)exp(X)
            sage: G.exp(X)*G.exp(Y)
            exp(Z)exp(Y)exp(X)
            sage: G.set_default_chart(G.chart_exp1())
            sage: G.exp(X)*G.exp(Y)
            exp(X + Y + 1/2*Z)
        """

        def __init__(self, parent, **kwds):
            r"""
            Initialize ``self``.

            TESTS::

                sage: L.<X,Y,Z> = LieAlgebra(QQ, 2, step=2)
                sage: G = L.lie_group()
                sage: g = G.exp(X)
                sage: TestSuite(g).run()
            """
            MultiplicativeGroupElement.__init__(self, parent)
            DifferentiableManifold.Element.__init__(self, parent, **kwds)

        def __invert__(self):
            r"""
            Return the inverse of ``self``.

            EXAMPLES::

                sage: L.<X,Y,Z> = LieAlgebra(QQ, 2, step=2)
                sage: G = L.lie_group('H')
                sage: g = G.point([var('a'), var('b'), var('c')]); g
                exp(a*X + b*Y + c*Z)
                sage: ~g
                exp((-a)*X + (-b)*Y + (-c)*Z)
                sage: g*~g
                exp(0)
            """
            G = self.parent()
            x = self.coordinates(chart=G._Exp1)
            return G.point(tuple(-xk for xk in x), chart=G._Exp1)

        def _mul_(self, other):
            r"""
            Return the product ``self`` * ``other``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 2, step=2)
                sage: G = L.lie_group('H')
                sage: g1 = G.point([2, 0, 0]); g1
                exp(2*X_1)
                sage: g2 = G.point([0, 1/3, 0]); g2
                exp(1/3*X_2)
                sage: g1*g2
                exp(2*X_1 + 1/3*X_2 + 1/3*X_12)
            """
            G = self.parent()
            a, b = G._group_law_vars
            self_c = list(zip(a, self.coordinates(chart=G._Exp1)))
            other_c = list(zip(b, other.coordinates(chart=G._Exp1)))
            sd = dict(self_c + other_c)
            return G.point([gk.expand() for gk in G._group_law.subs(sd)],
                           chart=G._Exp1)

        def _repr_(self):
            r"""
            Return a string representation of ``self``.

            Supports printing in exponential coordinates of the first and
            second kinds, depending on the default coordinate system.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 2, step=2)
                sage: G = L.lie_group('H')
                sage: g = G.point([1, 2, 3]); g
                exp(X_1 + 2*X_2 + 3*X_12)
                sage: G.set_default_chart(G.chart_exp2())
                sage: g
                exp(4*X_12)exp(2*X_2)exp(X_1)
            """
            G = self.parent()
            chart = G.default_chart()
            if chart != G._Exp1:
                if chart != G.chart_exp2():
                    chart = G._Exp1

            x = self.coordinates(chart=chart)
            B = G.lie_algebra().basis()
            nonzero_pairs = [(Xk, xk) for Xk, xk in zip(B, x) if xk]

            if chart == G._Exp1:
                s = repr_lincomb(nonzero_pairs)
            else:
                s = ")exp(".join(repr_lincomb([(Xk, xk)])
                                 for Xk, xk in reversed(nonzero_pairs))
                if not s:
                    s = "0"
            return "exp(%s)" % s
