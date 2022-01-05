r"""
Affine Connections

The class :class:`AffineConnection` implements affine connections on
smooth manifolds.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Marco Mancini (2015) : parallelization of some computations
- Florentin Jaffredo (2018) : series expansion with respect to a given
  parameter

REFERENCES:

- [Lee1997]_
- [KN1963]_
- [ONe1983]_

"""
# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.parallel.decorate import parallel
from sage.parallel.parallelism import Parallelism

class AffineConnection(SageObject):
    r"""
    Affine connection on a smooth manifold.

    Let `M` be a differentiable manifold of class `C^\infty` (smooth manifold)
    over a non-discrete topological field `K` (in most applications `K=\RR`
    or `K=\CC`), let `C^\infty(M)` be the algebra of smooth functions
    `M\rightarrow K` (cf.
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`)
    and let `\mathfrak{X}(M)` be the `C^\infty(M)`-module of vector fields on
    `M` (cf.
    :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`).
    An *affine connection* on `M` is an operator

    .. MATH::

        \begin{array}{cccc}
        \nabla: & \mathfrak{X}(M)\times \mathfrak{X}(M) & \longrightarrow &
                 \mathfrak{X}(M) \\
                & (u,v) & \longmapsto & \nabla_u v
        \end{array}

    that

    - is `K`-bilinear, i.e. is bilinear when considering `\mathfrak{X}(M)` as a
      vector space over `K`
    - is `C^\infty(M)`-linear w.r.t. the first argument:
      `\forall f\in C^\infty(M),\ \nabla_{fu} v = f\nabla_u v`
    - obeys Leibniz rule w.r.t. the second argument:
      `\forall f\in C^\infty(M),\ \nabla_u (f v) = \mathrm{d}f(u)\, v + f  \nabla_u v`

    The affine connection `\nabla` gives birth to the *covariant derivative
    operator* acting on tensor fields, denoted by the same symbol:

    .. MATH::

        \begin{array}{cccc}
        \nabla: &  T^{(k,l)}(M) & \longrightarrow & T^{(k,l+1)}(M)\\
                & t & \longmapsto & \nabla t
        \end{array}

    where `T^{(k,l)}(M)` stands for the `C^\infty(M)`-module of tensor fields
    of type `(k,l)` on `M` (cf.
    :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`),
    with the convention `T^{(0,0)}(M):=C^\infty(M)`.
    For a vector field `v`,  the covariant derivative `\nabla v` is a
    type-(1,1) tensor field such that

    .. MATH::

        \forall u \in\mathfrak{X}(M), \   \nabla_u v = \nabla v(., u)

    More generally for any tensor field `t\in T^{(k,l)}(M)`, we have

    .. MATH::

        \forall u \in\mathfrak{X}(M), \   \nabla_u t = \nabla t(\ldots, u)


    .. NOTE::

        The above convention means that, in terms of index notation,
        the "derivation index" in `\nabla t` is the *last* one:

        .. MATH::

            \nabla_c t^{a_1\ldots a_k}_{\quad\quad b_1\ldots b_l} =
                (\nabla t)^{a_1\ldots a_k}_{\quad\quad b_1\ldots b_l c}


    INPUT:

    - ``domain`` -- the manifold on which the connection is defined
      (must be an instance of class
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
    - ``name`` -- name given to the affine connection
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the affine
      connection; if ``None``, it is set to ``name``.

    EXAMPLES:

    Affine connection on a 3-dimensional manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: nab = M.affine_connection('nabla', r'\nabla') ; nab
        Affine connection nabla on the 3-dimensional differentiable manifold M

    A just-created connection has no connection coefficients::

        sage: nab._coefficients
        {}

    The connection coefficients relative to the manifold's default frame
    [here `(\partial/\partial x, \partial/\partial y, \partial/\partial z)`],
    are created by providing the relevant indices inside square brackets::

        sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
        sage: nab._coefficients
        {Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z)): 3-indices components w.r.t.
         Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))}

    If not the default one, the vector frame w.r.t. which the connection
    coefficients are defined can be specified as the first argument inside the
    square brackets; hence the above definition is equivalent to::

        sage: nab[c_xyz.frame(), 1,1,2], nab[c_xyz.frame(),3,2,3] = x^2, y*z
        sage: nab._coefficients
        {Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z)): 3-indices components w.r.t.
         Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))}

    Unset components are initialized to zero::

        sage: nab[:] # list of coefficients relative to the manifold's default vector frame
        [[[0, x^2, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, y*z], [0, 0, 0]]]

    The treatment of connection coefficients in a given vector frame is similar
    to that of tensor components; see therefore the class
    :class:`~sage.manifolds.differentiable.tensorfield.TensorField` for the
    documentation. In particular, the square brackets return the connection
    coefficients as instances of
    :class:`~sage.manifolds.chart_func.ChartFunction`,
    while the double square brackets return a scalar field::

        sage: nab[1,1,2]
        x^2
        sage: nab[1,1,2].display()
        (x, y, z) ↦ x^2
        sage: type(nab[1,1,2])
        <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
        sage: nab[[1,1,2]]
        Scalar field on the 3-dimensional differentiable manifold M
        sage: nab[[1,1,2]].display()
        M → ℝ
        (x, y, z) ↦ x^2
        sage: nab[[1,1,2]].coord_function() is nab[1,1,2]
        True

    Action on a scalar field::

        sage: f = M.scalar_field(x^2 - y^2, name='f')
        sage: Df = nab(f) ; Df
        1-form df on the 3-dimensional differentiable manifold M
        sage: Df[:]
        [2*x, -2*y, 0]

    The action of an affine connection on a scalar field must
    coincide with the differential::

        sage: Df == f.differential()
        True

    A generic affine connection has some torsion::

        sage: DDf = nab(Df) ; DDf
        Tensor field nabla(df) of type (0,2) on the 3-dimensional
         differentiable manifold M
        sage: DDf.antisymmetrize()[:] # nabla does not commute on scalar fields:
        [   0 -x^3    0]
        [ x^3    0    0]
        [   0    0    0]

    Let us check the standard formula

    .. MATH::

        \nabla_j \nabla_i \, f - \nabla_i \nabla_j \, f =
            T^k_{\ \, ij} \nabla_k \, f ,

    where the `T^k_{\ \, ij}`'s are the components of the connection's
    torsion tensor::

        sage: 2*DDf.antisymmetrize() == nab.torsion().contract(0,Df)
        True

    The connection acting on a vector field::

        sage: v = M.vector_field(y*z, x*z, x*y, name='v')
        sage: Dv = nab(v) ; Dv
        Tensor field nabla(v) of type (1,1) on the 3-dimensional differentiable
         manifold M
        sage: Dv[:]
        [            0 (x^2*y + 1)*z             y]
        [            z             0             x]
        [            y             x       x*y*z^2]

    Another example: connection on a non-parallelizable 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
        ....:                              restrictions1= x>0, restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: c_xyW = c_xy.restrict(W) ; c_uvW = c_uv.restrict(W)
        sage: eUW = c_xyW.frame() ; eVW = c_uvW.frame()
        sage: nab = M.affine_connection('nabla', r'\nabla')

    The connection is first defined on the open subset U by means of its
    coefficients w.r.t. the frame eU (the manifold's default frame)::

        sage: nab[0,0,0], nab[1,0,1] = x, x*y

    The coefficients w.r.t the frame eV are deduced by continuation of the
    coefficients w.r.t. the frame eVW on the open subset `W=U\cap V`::

        sage: for i in M.irange():
        ....:     for j in M.irange():
        ....:         for k in M.irange():
        ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()

    At this stage, the connection is fully defined on all the manifold::

        sage: nab.coef(eU)[:]
        [[[x, 0], [0, 0]], [[0, x*y], [0, 0]]]
        sage: nab.coef(eV)[:]
        [[[1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v, -1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v],
          [1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v, -1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v]],
         [[-1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v, 1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v],
          [-1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v, 1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v]]]

    We may let it act on a vector field defined globally on `M`::

        sage: a = M.vector_field({eU: [-y,x]}, name='a')
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = -y ∂/∂x + x ∂/∂y
        sage: a.display(eV)
        a = v ∂/∂u - u ∂/∂v
        sage: da = nab(a) ; da
        Tensor field nabla(a) of type (1,1) on the 2-dimensional differentiable
         manifold M
        sage: da.display(eU)
        nabla(a) = -x*y ∂/∂x⊗dx - ∂/∂x⊗dy + ∂/∂y⊗dx - x*y^2 ∂/∂y⊗dy
        sage: da.display(eV)
        nabla(a) = (-1/16*u^3 + 1/16*u^2*v + 1/16*(u + 2)*v^2 - 1/16*v^3 - 1/8*u^2) ∂/∂u⊗du
         + (1/16*u^3 - 1/16*u^2*v - 1/16*(u - 2)*v^2 + 1/16*v^3 - 1/8*u^2 + 1) ∂/∂u⊗dv
         + (1/16*u^3 - 1/16*u^2*v - 1/16*(u - 2)*v^2 + 1/16*v^3 - 1/8*u^2 - 1) ∂/∂v⊗du
         + (-1/16*u^3 + 1/16*u^2*v + 1/16*(u + 2)*v^2 - 1/16*v^3 - 1/8*u^2) ∂/∂v⊗dv

    A few tests::

        sage: nab(a.restrict(V)) == da.restrict(V)
        True
        sage: nab.restrict(V)(a) == da.restrict(V)
        True
        sage: nab.restrict(V)(a.restrict(U)) == da.restrict(W)
        True
        sage: nab.restrict(U)(a.restrict(V)) == da.restrict(W)
        True

    Same examples with SymPy as the engine for symbolic calculus::

        sage: M.set_calculus_method('sympy')
        sage: nab = M.affine_connection('nabla', r'\nabla')
        sage: nab[0,0,0], nab[1,0,1] = x, x*y
        sage: for i in M.irange():
        ....:     for j in M.irange():
        ....:         for k in M.irange():
        ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()

    At this stage, the connection is fully defined on all the manifold::

        sage: nab.coef(eU)[:]
        [[[x, 0], [0, 0]], [[0, x*y], [0, 0]]]
        sage: nab.coef(eV)[:]
        [[[u**2/16 + u/8 - v**2/16 + v/8, -u**2/16 + u/8 + v**2/16 + v/8],
         [u**2/16 + u/8 - v**2/16 + v/8, -u**2/16 + u/8 + v**2/16 + v/8]],
        [[-u**2/16 + u/8 + v**2/16 + v/8, u**2/16 + u/8 - v**2/16 + v/8],
         [-u**2/16 + u/8 + v**2/16 + v/8, u**2/16 + u/8 - v**2/16 + v/8]]]

    We may let it act on a vector field defined globally on `M`::

        sage: a = M.vector_field({eU: [-y,x]}, name='a')
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = -y ∂/∂x + x ∂/∂y
        sage: a.display(eV)
        a = v ∂/∂u - u ∂/∂v
        sage: da = nab(a) ; da
        Tensor field nabla(a) of type (1,1) on the 2-dimensional differentiable
         manifold M
        sage: da.display(eU)
        nabla(a) = -x*y ∂/∂x⊗dx - ∂/∂x⊗dy + ∂/∂y⊗dx - x*y**2 ∂/∂y⊗dy
        sage: da.display(eV)
        nabla(a) = (-u**3/16 + u**2*v/16 - u**2/8 + u*v**2/16 - v**3/16 + v**2/8) ∂/∂u⊗du
         + (u**3/16 - u**2*v/16 - u**2/8 - u*v**2/16 + v**3/16 + v**2/8 + 1) ∂/∂u⊗dv
         + (u**3/16 - u**2*v/16 - u**2/8 - u*v**2/16 + v**3/16 + v**2/8 - 1) ∂/∂v⊗du
         + (-u**3/16 + u**2*v/16 - u**2/8 + u*v**2/16 - v**3/16 + v**2/8) ∂/∂v⊗dv

    To make affine connections hashable, they have to be set immutable before::

        sage: nab.is_immutable()
        False
        sage: nab.set_immutable()
        sage: nab.is_immutable()
        True

    Immutable connections cannot be changed anymore::

        sage: nab.set_coef(eU)
        Traceback (most recent call last):
        ...
        ValueError: the coefficients of an immutable element cannot be
         changed

    However, they can now be used as keys for dictionaries::

        sage: {nab: 1}[nab]
        1

    The immutability process cannot be made undone. If a connection is
    needed to be changed again, a copy has to be created::

        sage: nab_copy = nab.copy('nablo'); nab_copy
        Affine connection nablo on the 2-dimensional differentiable manifold M
        sage: nab_copy is nab
        False
        sage: nab_copy == nab
        True
        sage: nab_copy.is_immutable()
        False

    """
    def __init__(self, domain, name, latex_name=None):
        r"""
        Construct an affine connection.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: from sage.manifolds.differentiable.affine_connection import \
                                                                AffineConnection
            sage: nab = AffineConnection(M, 'nabla', latex_name=r'\nabla')
            sage: nab
            Affine connection nabla on the 3-dimensional differentiable
             manifold M
            sage: X.<x,y,z> = M.chart()
            sage: nab[0,1,0] = x*y*z
            sage: TestSuite(nab).run()

        """
        if not isinstance(domain, DifferentiableManifold):
            raise TypeError("the first argument must be a differentiable " +
                            "manifold")
        self._is_immutable = False
        self._domain = domain
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._coefficients = {}  # dict. of connection coefficients, with the
                                 # vector frames as keys
        # Initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab._repr_()
            'Affine connection nabla on the 5-dimensional differentiable manifold M'
            sage: repr(nab)  # indirect doctest
            'Affine connection nabla on the 5-dimensional differentiable manifold M'

        """
        description = "Affine connection"
        if self._name is not None:
            description += " " + self._name
        description += " on the {}".format(self._domain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab._latex_()
            '\\nabla'
            sage: latex(nab)  # indirect doctest
            \nabla
            sage: nab = M.affine_connection('D')
            sage: nab._latex_()
            'D'
            sage: latex(nab)  # indirect doctest
            D

        """
        return self._latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(4, 'M')
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab._init_derived()

        """
        self._restrictions = {} # dict. of restrictions of ``self`` on some
                                # subdomains, with the subdomains as keys
        self._torsion = None
        self._riemann = None
        self._ricci = None
        self._connection_forms = {}  # dict. of dict. of connection 1-forms
                                     # (key: vector frame)
        self._torsion_forms = {}  # dict. of dict. of torsion 1-forms
                                  # (key: vector frame)
        self._curvature_forms = {}  # dict. of dict. of curvature 2-forms
                                    # (key: vector frame)

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(4, 'M')
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab._del_derived()

        """
        self._restrictions.clear()
        self._torsion = None
        self._riemann = None
        self._ricci = None
        self._connection_forms.clear()
        self._torsion_forms.clear()
        self._curvature_forms.clear()

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- an affine connection

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab[0,1,0], nab[0,1,1] = 1+x, x*y
            sage: nab.display()
            Gam^x_yx = x + 1
            Gam^x_yy = x*y
            sage: nab1 = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: (nab1 == nab) or (nab == nab1)
            False
            sage: nab1[0,1,0], nab1[0,1,1] = 2, 3-y
            sage: (nab1 == nab) or (nab == nab1)
            False
            sage: nab1[0,1,0], nab1[0,1,1] = 1+x, x*y
            sage: (nab1 == nab) and (nab == nab1)
            True
            sage: nab2 = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: a = M.automorphism_field()
            sage: a[:] = [[0,1], [1,0]]
            sage: e = X.frame().new_frame(a, 'e')
            sage: nab2.set_coef(e)[1,0,1] = 1+x
            sage: nab2.set_coef(e)[1,0,0] = x*y
            sage: (nab2 == nab) and (nab == nab2)
            True
            sage: f = M.vector_frame('f')
            sage: nab2.set_coef(f)[1,0,1] = x-y
            sage: (nab2 == nab) or (nab == nab2)
            False

        """
        if other is self:
            return True
        if not isinstance(other, AffineConnection):
            return False
        if other._domain != self._domain:
            return False
        if self._coefficients == {}:
            return False
        for frame, coef in self._coefficients.items():
            try:
                if other.coef(frame) != coef:
                    return False
            except ValueError:
                return False
        return True

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- an affine connection

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab[0,1,0], nab[0,1,1] = 1+x, x*y
            sage: nab1 = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: (nab1 != nab) and (nab != nab1)
            True
            sage: nab1[0,1,0], nab1[0,1,1] = 2, 3-y
            sage: (nab1 != nab) and (nab != nab1)
            True
            sage: nab1[0,1,0], nab1[0,1,1] = 1+x, x*y
            sage: (nab1 != nab) or (nab != nab1)
            False

        """
        return not (self == other)

    def domain(self):
        r"""
        Return the manifold subset on which the affine connection is defined.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          representing the manifold on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab.domain()
            3-dimensional differentiable manifold M
            sage: U = M.open_subset('U', coord_def={c_xyz: x>0})
            sage: nabU = U.affine_connection('D')
            sage: nabU.domain()
            Open subset U of the 3-dimensional differentiable manifold M

        """
        return self._domain

    def _new_coef(self, frame):
        r"""
        Create the connection coefficients w.r.t. the given frame.

        This method, to be called by :meth:`coef`, must be redefined by derived
        classes to adapt the output to the relevant subclass of
        :class:`~sage.tensor.modules.comp.Components`.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab._new_coef(X.frame())
            3-indices components w.r.t. Coordinate frame (M, (∂/∂x,∂/∂y))

        """
        from sage.tensor.modules.comp import Components
        from sage.manifolds.differentiable.scalarfield import DiffScalarField
        return Components(frame._domain.scalar_field_algebra(), frame, 3,
                          start_index=self._domain._sindex,
                          output_formatter=DiffScalarField.coord_function)

    def coef(self, frame=None):
        r"""
        Return the connection coefficients relative to the given frame.

        `n` being the manifold's dimension, the connection coefficients
        relative to the vector frame `(e_i)` are the `n^3` scalar fields
        `\Gamma^k_{\ \, ij}` defined by

        .. MATH::

            \nabla_{e_j} e_i = \Gamma^k_{\ \, ij} e_k


        If the connection coefficients are not known already, they are computed
        from the above formula.

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame relative to which the
          connection coefficients are required; if none is provided, the
          domain's default frame is assumed

        OUTPUT:

        - connection coefficients relative to the frame ``frame``, as an
          instance of the class :class:`~sage.tensor.modules.comp.Components`
          with 3 indices ordered as `(k,i,j)`

        EXAMPLES:

        Connection coefficient of an affine connection on a 3-dimensional
        manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: nab.coef()
            3-indices components w.r.t. Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))
            sage: type(nab.coef())
            <class 'sage.tensor.modules.comp.Components'>
            sage: M.default_frame()
            Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))
            sage: nab.coef() is nab.coef(c_xyz.frame())
            True
            sage: nab.coef()[:]  # full list of coefficients:
            [[[0, x^2, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, y*z], [0, 0, 0]]]

        """
        if frame is None:
            frame = self._domain.default_frame()
        if frame not in self._coefficients:
            # the coefficients must be computed
            #
            # Check whether frame is a subframe of a frame in which the
            # coefficients are already known:
            for oframe in self._coefficients:
                if frame in oframe._subframes:
                    self._coefficients[frame] = self._new_coef(frame)
                    comp_store = self._coefficients[frame]._comp
                    ocomp_store = self._coefficients[oframe]._comp
                    for ind, value in ocomp_store.items():
                        comp_store[ind] = value.restrict(frame._domain)
                    break
            else:
                # If not, the coefficients must be computed from scratch:
                manif = self._domain
                ev = frame        # the vector frame
                ef = ev.coframe() # the dual frame
                gam = self._new_coef(ev)
                for i in manif.irange():
                    nab_evi = self(ev[i])
                    for k in manif.irange():
                        for j in manif.irange():
                            gam[[k,i,j]] = nab_evi(ef[k],ev[j])
                self._coefficients[frame] = gam
        return self._coefficients[frame]

    def set_coef(self, frame=None):
        r"""
        Return the connection coefficients in a given frame for assignment.

        See method :meth:`coef` for details about the definition of the
        connection coefficients.

        The connection coefficients with respect to other frames are deleted,
        in order to avoid any inconsistency. To keep them, use the method
        :meth:`add_coef` instead.

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame in which the connection
          coefficients are defined; if ``None``, the default frame of the
          connection's domain is assumed.

        OUTPUT:

        - connection coefficients in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          connection coefficients did not exist previously, they are created.
          See method :meth:`coef` for the storage convention of the connection
          coefficients.

        EXAMPLES:

        Setting the coefficients of an affine connection w.r.t. some coordinate
        frame::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: eX = X.frame(); eX
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: nab.set_coef(eX)
            3-indices components w.r.t. Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: nab.set_coef(eX)[1,2,1] = x*y
            sage: nab.display(eX)
            Gam^x_yx = x*y

        Since ``eX`` is the manifold's default vector frame, its mention may
        be omitted::

            sage: nab.set_coef()[1,2,1] = x*y
            sage: nab.set_coef()
            3-indices components w.r.t. Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: nab.set_coef()[1,2,1] = x*y
            sage: nab.display()
            Gam^x_yx = x*y

        To set the coefficients in the default frame, one can even bypass the
        method ``set_coef()`` and call directly the operator ``[]`` on the
        connection object::

            sage: nab[1,2,1] = x*y
            sage: nab.display()
                Gam^x_yx = x*y

        Setting the connection coefficients w.r.t. to another vector frame::

            sage: e = M.vector_frame('e')
            sage: nab.set_coef(e)
            3-indices components w.r.t. Vector frame (M, (e_1,e_2))
            sage: nab.set_coef(e)[2,1,1] = x+y
            sage: nab.set_coef(e)[2,1,2] = x-y
            sage: nab.display(e)
            Gam^2_11 = x + y
            Gam^2_12 = x - y

        The coefficients w.r.t. the frame ``eX`` have been deleted::

            sage: nab.display(eX)
            Traceback (most recent call last):
            ...
            ValueError: no common frame found for the computation

        To keep them, use the method :meth:`add_coef` instead.

        """
        if self.is_immutable():
            raise ValueError("the coefficients of an immutable element "
                             "cannot be changed")
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._coefficients:
            if frame not in self._domain._frames:
                raise ValueError("the {} is not".format(frame) +
                                 " a frame on the {}".format(self._domain))
            self._coefficients[frame] = self._new_coef(frame)
        self._del_derived() # deletes the derived quantities
        self.del_other_coef(frame)
        return self._coefficients[frame]

    def add_coef(self, frame=None):
        r"""
        Return the connection coefficients in a given frame for assignment,
        keeping the coefficients in other frames.

        See method :meth:`coef` for details about the definition of the
        connection coefficients.

        To delete the connection coefficients in other frames, use the method
        :meth:`set_coef` instead.

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame in which the connection
          coefficients are defined; if ``None``, the default frame of the
          connection's domain is assumed.

        .. WARNING::

            If the connection has already coefficients in other frames, it
            is the user's responsibility to make sure that the coefficients
            to be added are consistent with them.

        OUTPUT:

        - connection coefficients in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          connection coefficients did not exist previously, they are created.
          See method :meth:`coef` for the storage convention of the connection
          coefficients.


        EXAMPLES:

        Setting the coefficients of an affine connection w.r.t. some coordinate
        frame::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: eX = X.frame(); eX
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: nab.add_coef(eX)
            3-indices components w.r.t. Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: nab.add_coef(eX)[1,2,1] = x*y
            sage: nab.display(eX)
            Gam^x_yx = x*y

        Since ``eX`` is the manifold's default vector frame, its mention may
        be omitted::

            sage: nab.add_coef()[1,2,1] = x*y
            sage: nab.add_coef()
            3-indices components w.r.t. Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: nab.add_coef()[1,2,1] = x*y
            sage: nab.display()
            Gam^x_yx = x*y

        Adding connection coefficients w.r.t. to another vector frame::

            sage: e = M.vector_frame('e')
            sage: nab.add_coef(e)
            3-indices components w.r.t. Vector frame (M, (e_1,e_2))
            sage: nab.add_coef(e)[2,1,1] = x+y
            sage: nab.add_coef(e)[2,1,2] = x-y
            sage: nab.display(e)
            Gam^2_11 = x + y
            Gam^2_12 = x - y

        The coefficients w.r.t. the frame ``eX`` have been kept::

            sage: nab.display(eX)
            Gam^x_yx = x*y

        To delete them, use the method :meth:`set_coef` instead.


        """
        if self.is_immutable():
            raise ValueError("the coefficients of an immutable element "
                             "cannot be changed")
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._coefficients:
            if frame not in self._domain._frames:
                raise ValueError("the {} is not".format(frame) +
                                 " a frame on the {}".format(self._domain))
            self._coefficients[frame] = self._new_coef(frame)
        self._del_derived() # deletes the derived quantities
        return self._coefficients[frame]

    def del_other_coef(self, frame=None):
        r"""
        Delete all the coefficients but those corresponding to ``frame``.

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame, the connection
          coefficients w.r.t. which are to be kept; if ``None``, the default
          frame of the connection's domain is assumed.

        EXAMPLES:

        We first create two sets of connection coefficients::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: eX = X.frame()
            sage: nab.set_coef(eX)[1,2,1] = x*y
            sage: e = M.vector_frame('e')
            sage: nab.add_coef(e)[2,1,1] = x+y
            sage: nab.display(eX)
            Gam^x_yx = x*y
            sage: nab.display(e)
            Gam^2_11 = x + y

        Let us delete the connection coefficients w.r.t. all frames except for
        frame ``eX``::

            sage: nab.del_other_coef(eX)
            sage: nab.display(eX)
            Gam^x_yx = x*y

        The connection coefficients w.r.t. frame ``e`` have indeed been
        deleted::

            sage: nab.display(e)
            Traceback (most recent call last):
            ...
            ValueError: no common frame found for the computation

        """
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._coefficients:
            raise ValueError("the coefficients w.r.t. {}".format(frame) +
                             " have not been defined")
        to_be_deleted = []
        for other_frame in self._coefficients:
            if other_frame != frame:
                to_be_deleted.append(other_frame)
        for other_frame in to_be_deleted:
            del self._coefficients[other_frame]

    def set_immutable(self):
        r"""
        Set ``self`` and all restrictions of ``self`` immutable.

        EXAMPLES:

        An affine connection can be set immutable::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1})
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: eX = X.frame()
            sage: nab.set_coef(eX)[1,2,1] = x*y
            sage: nab.is_immutable()
            False
            sage: nab.set_immutable()
            sage: nab.is_immutable()
            True

        The coefficients of immutable elements cannot be changed::

            sage: nab.add_coef(eX)[2,1,1] = x+y
            Traceback (most recent call last):
            ...
            ValueError: the coefficients of an immutable element cannot
             be changed

        The restriction are set immutable as well::

            sage: nabU = nab.restrict(U)
            sage: nabU.is_immutable()
            True

        """
        for rst in self._restrictions.values():
            rst.set_immutable()
        self._is_immutable = True

    def is_immutable(self):
        r"""
        Return ``True`` if this object is immutable, i.e. its coefficients
        cannot be chanced, and ``False`` if it is not.

        To set an affine connection immutable, use :meth:`set_immutable`.

        EXAMPLES::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab.is_immutable()
            False
            sage: nab.set_immutable()
            sage: nab.is_immutable()
            True

        """
        return self._is_immutable

    def is_mutable(self):
        r"""
        Return ``True`` if this object is mutable, i.e. its coefficients can
        be changed, and ``False`` if it is not.

        EXAMPLES::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab.is_mutable()
            True
            sage: nab.set_immutable()
            sage: nab.is_mutable()
            False

        """
        return not self._is_immutable

    def copy(self, name, latex_name=None):
        r"""
        Return an exact copy of ``self``.

        INPUT:

        - ``name`` -- name given to the copy
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          copy; if none is provided, the LaTeX symbol is set to ``name``

        .. NOTE::

            The name and the derived quantities are not copied.

        EXAMPLES::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: eX = X.frame()
            sage: nab.set_coef(eX)[1,2,1] = x*y
            sage: nab.set_coef(eX)[1,2,2] = x+y
            sage: nab.display()
            Gam^x_yx = x*y
            Gam^x_yy = x + y
            sage: nab_copy = nab.copy(name='nabla_1', latex_name=r'\nabla_1')
            sage: nab is nab_copy
            False
            sage: nab == nab_copy
            True
            sage: nab_copy.display()
            Gam^x_yx = x*y
            Gam^x_yy = x + y

        """
        copy = type(self)(self._domain, name, latex_name=latex_name)
        for dom, rst in self._restrictions.items():
            copy._restrictions[dom] = rst.copy(name, latex_name=latex_name)
        for frame, coef in self._coefficients.items():
            copy._coefficients[frame] = coef.copy()
        return copy

    def __getitem__(self, args):
        r"""
        Return the connection coefficient w.r.t. some frame corresponding to
        the given indices.

        INPUT:

        - ``args`` -- list of indices defining the coefficient; if ``[:]`` is
          provided, all the coefficients are returned. The frame can be passed
          as the first item of ``args``; if not, the default frame of the
          connection's domain is assumed

        OUTPUT:

        - the connection coefficient corresponding to the specified frame and
          indices, as an instance of
          :class:`~sage.manifolds.chart_func.ChartFunction`
          (or the list of all connection coefficients if ``args==[:]`` or
          ``args=[frame,:]``).

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab.set_coef(X.frame())[1,2,1] = x*y
            sage: nab.__getitem__((1,2,1))
            x*y
            sage: nab[1,2,1]  # equivalent to above
            x*y
            sage: type(nab.__getitem__((1,2,1)))
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: nab.__getitem__((X.frame(),1,2,1))
            x*y
            sage: nab[X.frame(),1,2,1]  # equivalent to above
            x*y

        Returning the full set of coefficients::

            sage: nab.__getitem__(slice(None))
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]
            sage: nab[:]  # equivalent to above
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]
            sage: nab.__getitem__((X.frame(), slice(None)))
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]
            sage: nab[X.frame(), :]  # equivalent to above
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]

        Returning a scalar field::

            sage: nab.__getitem__(([1,2,1]))
            Scalar field on the 2-dimensional differentiable manifold M
            sage: nab[[1,2,1]]  # equivalent to above
            Scalar field on the 2-dimensional differentiable manifold M
            sage: nab.__getitem__(([X.frame(),1,2,1])).coord_function() is nab[1,2,1]
            True

        """
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], (int, Integer, slice)):
                frame = self._domain._def_frame
            else:
                frame = args[0]
                args = args[1:]
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._domain._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
                if len(args) == 1:
                    args = args[0]  # to accommodate for [e,:] syntax
            else:
                frame = self._domain._def_frame
        return self.coef(frame)[args]

    def __setitem__(self, args, value):
        r"""
        Set the connection coefficient w.r.t. some frame corresponding to the
        given indices.

        INPUT:

        - ``args`` -- list of indices defining the coefficient; if ``[:]`` is
          provided, all the coefficients are set. The frame can be passed
          as the first item of ``args``; if not, the default frame of the
          connection's domain is assumed
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab.__setitem__((1,2,1), x*y)
            sage: nab[:]
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]
            sage: nab[1,2,1] = x*y  # equivalent to __setitem__ above
            sage: nab[:]
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]
            sage: nab.__setitem__((X.frame(),1,2,1), -x^2)
            sage: nab[1,2,1]
            -x^2
            sage: nab[X.frame(), 1,2,1] = -x^2  # equivalent to __setitem__ above
            sage: nab[1,2,1]
            -x^2

        Setting all the coefficients at once::

            sage: nab.__setitem__(slice(None),
            ....:                 [[[-x^2, 0], [x*y, 0]], [[0, 1+y], [0, 0]]])
            sage: nab[:]
            [[[-x^2, 0], [x*y, 0]], [[0, y + 1], [0, 0]]]
            sage: nab[:] = [[[-x^2, 0], [x*y, 0]], [[0, 1+y], [0, 0]]]  # equivalent to above
            sage: nab[:]
            [[[-x^2, 0], [x*y, 0]], [[0, y + 1], [0, 0]]]

        Providing a scalar field as value::

            sage: f = M.scalar_field({X: x*y})
            sage: nab.__setitem__((1,2,1), f)
            sage: nab[1,2,1]
            x*y

        """
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], (int, Integer, slice)):
                frame = self._domain._def_frame
            else:
                frame = args[0]
                args = args[1:]
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._domain._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
                if len(args) == 1:
                    args = args[0]  # to accommodate for [e,:] syntax
            else:
                frame = self._domain._def_frame
        self.set_coef(frame)[args] = value

    def display(self, frame=None, chart=None, symbol=None, latex_symbol=None,
                index_labels=None, index_latex_labels=None,
                coordinate_labels=True, only_nonzero=True,
                only_nonredundant=False):
        r"""
        Display all the connection coefficients w.r.t. to a given frame, one
        per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame relative to which the
          connection coefficients are defined; if ``None``, the
          default frame of the connection's domain is used
        - ``chart`` -- (default: ``None``) chart specifying the coordinate
          expression of the connection coefficients; if ``None``,
          the default chart of the domain of ``frame`` is used
        - ``symbol`` -- (default: ``None``) string specifying the
          symbol of the connection coefficients; if ``None``, 'Gam' is used
        - ``latex_symbol`` -- (default: ``None``) string specifying the LaTeX
          symbol for the components; if ``None``, '\\Gamma' is used
        - ``index_labels`` -- (default: ``None``) list of strings representing
          the labels of each index; if ``None``, integer labels are used,
          except if ``frame`` is a coordinate frame and ``coordinate_symbols``
          is set to ``True``, in which case the coordinate symbols are used
        - ``index_latex_labels`` -- (default: ``None``) list of strings
          representing the LaTeX labels of each index; if ``None``, integer
          labels are used, except if ``frame`` is a coordinate frame and
          ``coordinate_symbols`` is set to ``True``, in which case the
          coordinate LaTeX symbols are used
        - ``coordinate_labels`` -- (default: ``True``) boolean; if ``True``,
          coordinate symbols are used by default (instead of integers) as
          index labels whenever ``frame`` is a coordinate frame
        - ``only_nonzero`` -- (default: ``True``) boolean; if ``True``, only
          nonzero connection coefficients are displayed
        - ``only_nonredundant`` -- (default: ``False``) boolean; if ``True``,
          only nonredundant connection coefficients are displayed in case of
          symmetries

        EXAMPLES:

        Coefficients of a connection on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z

        By default, only the nonzero connection coefficients are displayed::

            sage: nab.display()
            Gam^x_xy = x^2
            Gam^z_yz = y*z
            sage: latex(nab.display())
            \begin{array}{lcl} \Gamma_{ \phantom{\, x} \, x \, y }^{ \, x \phantom{\, x} \phantom{\, y} }
            & = & x^{2} \\
            \Gamma_{ \phantom{\, z} \, y \, z }^{ \, z \phantom{\, y} \phantom{\, z} }
            & = & y z \end{array}

        By default, the displayed connection coefficients are those w.r.t.
        to the default frame of the connection's domain, so the above is
        equivalent to::

            sage: nab.display(frame=M.default_frame())
            Gam^x_xy = x^2
            Gam^z_yz = y*z

        Since the default frame is a coordinate frame, coordinate symbols are
        used to label the indices, but one may ask for integers instead::

            sage: M.default_frame() is c_xyz.frame()
            True
            sage: nab.display(coordinate_labels=False)
            Gam^1_12 = x^2
            Gam^3_23 = y*z

        The index labels can also be customized::

            sage: nab.display(index_labels=['(1)', '(2)', '(3)'])
            Gam^(1)_(1),(2) = x^2
            Gam^(3)_(2),(3) = y*z

        The symbol 'Gam' can be changed::

            sage: nab.display(symbol='C', latex_symbol='C')
            C^x_xy = x^2
            C^z_yz = y*z
            sage: latex(nab.display(symbol='C', latex_symbol='C'))
            \begin{array}{lcl} C_{ \phantom{\, x} \, x \, y }^{ \, x \phantom{\, x} \phantom{\, y} }
            & = & x^{2} \\
            C_{ \phantom{\, z} \, y \, z }^{ \, z \phantom{\, y} \phantom{\, z} }
            & = & y z \end{array}

        Display of Christoffel symbols, skipping the redundancy associated
        with the symmetry of the last two indices::

            sage: M = Manifold(3, 'R^3', start_index=1)
            sage: c_spher.<r,th,ph> = M.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2 , (r*sin(th))^2
            sage: g.display()
            g = dr⊗dr + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph
            sage: g.connection().display(only_nonredundant=True)
            Gam^r_th,th = -r
            Gam^r_ph,ph = -r*sin(th)^2
            Gam^th_r,th = 1/r
            Gam^th_ph,ph = -cos(th)*sin(th)
            Gam^ph_r,ph = 1/r
            Gam^ph_th,ph = cos(th)/sin(th)

        By default, the parameter ``only_nonredundant`` is set to ``False``::

            sage: g.connection().display()
            Gam^r_th,th = -r
            Gam^r_ph,ph = -r*sin(th)^2
            Gam^th_r,th = 1/r
            Gam^th_th,r = 1/r
            Gam^th_ph,ph = -cos(th)*sin(th)
            Gam^ph_r,ph = 1/r
            Gam^ph_th,ph = cos(th)/sin(th)
            Gam^ph_ph,r = 1/r
            Gam^ph_ph,th = cos(th)/sin(th)

        """
        from sage.misc.latex import latex
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        if frame is None:
            frame = self._domain.default_frame()
        if chart is None:
            chart = frame.domain().default_chart()
        if symbol is None:
            symbol = 'Gam'
        if latex_symbol is None:
            latex_symbol = r'\Gamma'
        if index_labels is None and isinstance(frame, CoordFrame) and \
          coordinate_labels:
            ch = frame.chart()
            index_labels = [str(z) for z in ch[:]]
            index_latex_labels = [latex(z) for z in ch[:]]
        return self.coef(frame=frame).display(symbol,
              latex_symbol=latex_symbol, index_positions='udd',
              index_labels=index_labels, index_latex_labels=index_latex_labels,
              format_spec=chart, only_nonzero=only_nonzero,
              only_nonredundant=only_nonredundant)

    def restrict(self, subdomain):
        r"""
        Return the restriction of the connection to some subdomain.

        If such restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of the connection's domain (must be
          an instance of
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)

        OUTPUT:

        - instance of :class:`AffineConnection` representing the restriction.

        EXAMPLES:

        Restriction of a connection on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[1,1,2], nab[2,1,1] = x^2, x+y
            sage: nab[:]
            [[[0, x^2], [0, 0]], [[x + y, 0], [0, 0]]]
            sage: U = M.open_subset('U', coord_def={c_xy: x>0})
            sage: nabU = nab.restrict(U) ; nabU
            Affine connection nabla on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: nabU.domain()
            Open subset U of the 2-dimensional differentiable manifold M
            sage: nabU[:]
            [[[0, x^2], [0, 0]], [[x + y, 0], [0, 0]]]

        The result is cached::

            sage: nab.restrict(U) is nabU
            True

        until the connection is modified::

            sage: nab[1,2,2] = -y
            sage: nab.restrict(U) is nabU
            False
            sage: nab.restrict(U)[:]
            [[[0, x^2], [0, -y]], [[x + y, 0], [0, 0]]]

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("The provided domains is not a subset of " +
                                 "the connection's domain.")
            resu = AffineConnection(subdomain, name=self._name,
                                    latex_name=self._latex_name)
            for frame in self._coefficients:
                for sframe in subdomain._top_frames:
                    if sframe in frame._subframes:
                        comp_store = self._coefficients[frame]._comp
                        scoef = resu._new_coef(sframe)
                        scomp_store = scoef._comp
                        # the coefficients of the restriction are evaluated
                        # index by index:
                        for ind, value in comp_store.items():
                            scomp_store[ind] = value.restrict(sframe._domain)
                        resu._coefficients[sframe] = scoef
            if self._torsion is not None:
                resu._torsion = self._torsion.restrict(subdomain)
            if self._riemann is not None:
                resu._riemann = self._riemann.restrict(subdomain)
            if self._ricci is not None:
                resu._ricci = self._ricci.restrict(subdomain)
            resu.set_immutable()  # restrictions must be immutable, too
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def _common_frame(self, other):
        r"""
        Find a common vector frame for the coefficients of ``self`` and
        the components of  ``other``.

        In case of multiple common frames, the default frame of ``self``'s
        domain is privileged.

        INPUT:

        - ``other`` -- a tensor field on parallelizable domain, as an
          instance of
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`

        OUTPUT:

        - common frame; if no common frame is found, None is returned.

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab[1,2,1] = x*y
            sage: v = M.vector_field()
            sage: v[:] = [-y, x]
            sage: nab._common_frame(v)
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: e = M.vector_frame('e')
            sage: u = M.vector_field()
            sage: u[e,:] = [-3, 2]
            sage: nab._common_frame(u)  # no common frame is found

        """
        # The domain of search is restricted to other._domain:
        dom = other._domain
        # 1/ Does each object have components on the domain's default frame ?
        def_frame = dom._def_frame
        if def_frame in self._coefficients and def_frame in other._components:
            return def_frame
        # 2/ Search for a common frame among the existing components, i.e.
        #    without performing any component transformation.
        #    -------------------------------------------------------------
        for frame in self._coefficients:
            if frame in other._components:
                return frame
        # 3/ Search for a common frame among the subframes of self's frames:
        #    --------------------------------------------------------------
        for frame in self._coefficients:
            for oframe in other._components:
                if oframe in frame._subframes:
                    self.coef(oframe) # update the coefficients of self in oframe
                    return oframe
        #
        # 4/ Search for a common frame via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is necessary to perform at least
        # one component transformation to get a common frame
        for frame in self._coefficients:
            for oframe in other._components:
                if (oframe, frame) in dom._frame_changes:
                    other.comp(frame, from_basis=oframe)
                    return frame
        # 5/ Search for a common frame via one component transformation to
        #    a subframe of self's frames:
        #    -------------------------------------------------------------
        for frame in self._coefficients:
            for oframe in other._components:
                for sframe in frame._subframes:
                    if (oframe, sframe) in dom._frame_changes:
                        self.coef(sframe)
                        other.comp(sframe, from_basis=oframe)
                        return sframe
        #
        # If this point is reached, no common frame could be found, even at
        # the price of a component transformation:
        return None

    def __call__(self, tensor):
        r"""
        Action of the connection on a tensor field.

        INPUT:

        - ``tensor`` -- a tensor field `T`, of type `(k,\ell)`

        OUTPUT:

        - tensor field `\nabla T`.

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab[1,2,1] = x*y
            sage: v = M.vector_field()
            sage: v[:] = [-y, x]
            sage: nab.__call__(v)
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M

        See documentation of
        :class:`~sage.manifolds.differentiable.affine_connection.AffineConnection`
        for more examples.

        """
        from sage.manifolds.differentiable.tensorfield_paral import \
                                                               TensorFieldParal
        from sage.tensor.modules.format_utilities import format_unop_latex
        dom_resu = self._domain.intersection(tensor._domain)
        tensor_r = tensor.restrict(dom_resu)
        if tensor_r._tensor_type == (0,0):  # scalar field case
            return tensor_r.differential()
        if isinstance(tensor_r, TensorFieldParal):
            return self._derive_paral(tensor_r)
        resu_rst = []
        for dom, rst in tensor_r._restrictions.items():
            # the computation is performed only if dom is not a subdomain
            # of another restriction:
            for odom in tensor_r._restrictions:
                if dom in odom._subsets and dom is not odom:
                    break
            else:
                # dom is a not a subdomain and the computation is performed:
                resu_rst.append(self.__call__(rst))
        tensor_type_resu = (tensor_r._tensor_type[0],
                            tensor_r._tensor_type[1]+1)
        if tensor_r._name is None:
            name_resu = None
        else:
            name_resu = self._name + '(' + tensor_r._name + ')'
        if tensor_r._latex_name is None:
            latex_name_resu = None
        else:
            latex_name_resu = format_unop_latex(self._latex_name + ' ',
                                                          tensor_r._latex_name)
        vmodule = dom_resu.vector_field_module()
        resu = vmodule.tensor(tensor_type_resu, name=name_resu,
                              latex_name=latex_name_resu,
                              sym=resu_rst[0]._sym,
                              antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def _derive_paral(self, tensor):
        r"""
        Action of the connection on a tensor field on a parallelizable domain.

        INPUT:

        - ``tensor`` -- a tensor field `T`, of type `(k,\ell)`

        OUTPUT:

        - tensor field `\nabla T`.

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: nab = M.affine_connection('nabla', latex_name=r'\nabla')
            sage: nab[1,2,1] = x*y
            sage: v = M.vector_field()
            sage: v[:] = [-y, x]
            sage: nab._derive_paral(v)
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M

        """
        from sage.manifolds.differentiable.scalarfield import DiffScalarField
        from sage.tensor.modules.comp import Components, CompWithSym
        from sage.tensor.modules.format_utilities import format_unop_latex
        manif = self._domain
        tdom = tensor._domain
        frame = self._common_frame(tensor)
        if frame is None:
            raise ValueError("no common frame found for the computation")
        # Component computation in the common frame:
        tc = tensor._components[frame]
        gam = self._coefficients[frame]
        if not tensor._sym and not tensor._antisym:
            resc = Components(tdom.scalar_field_algebra(), frame,
                              tensor._tensor_rank+1,
                              start_index=self._domain._sindex,
                              output_formatter=DiffScalarField.coord_function)
        else:
            resc = CompWithSym(tdom.scalar_field_algebra(), frame,
                              tensor._tensor_rank+1,
                              start_index=self._domain._sindex,
                              output_formatter=DiffScalarField.coord_function,
                              sym=tensor._sym, antisym=tensor._antisym)
        n_con = tensor._tensor_type[0]
        n_cov = tensor._tensor_type[1]

        if Parallelism().get('tensor') != 1:
            # parallel computation
            # !!!!! Seems to work only when a frame is chosen !!!!!!

            nproc = Parallelism().get('tensor')
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]

            ind_list = list(resc.non_redundant_index_generator())
            ind_step = max(1,int(len(ind_list)/nproc/2))
            local_list = lol(ind_list,ind_step)

            # definition of the list of input parameters
            listParalInput = []
            for ind_part in local_list:
                listParalInput.append((ind_part,tc,gam,frame,n_con,
                                       tensor._tensor_rank,manif))

            # definition of the parallel function
            @parallel(p_iter='multiprocessing',ncpus=nproc)
            def make_CovDerivative(ind_part,tc,gam,frame,n_con,rank,manif):
                partial = []
                for ind in ind_part:
                    p = ind[-1]  # derivation index
                    ind0 = ind[:-1]
                    rsum = frame[p](tc[[ind0]])
                    # loop on contravariant indices:
                    for k in range(n_con):
                        for i in manif.irange():
                            indk = list(ind0)
                            indk[k] = i
                            rsum += gam[[ind0[k], i, p]] * tc[[indk]]
                    # loop on covariant indices:
                    for k in range(n_con, rank):
                        for i in manif.irange():
                            indk = list(ind0)
                            indk[k] = i
                            rsum -= gam[[i, ind0[k], p]] * tc[[indk]]
                    partial.append([ind,rsum])
                return partial

            # Computation and Assignation of values
            for ii,val in make_CovDerivative(listParalInput):
                for jj in val:
                    resc[[jj[0]]] = jj[1]

        else:
            # sequential
            for ind in resc.non_redundant_index_generator():
                p = ind[-1]  # derivation index
                ind0 = ind[:-1]
                rsum = frame[p](tc[[ind0]])
                # loop on contravariant indices:
                for k in range(n_con):
                    for i in manif.irange():
                        indk = list(ind0)
                        indk[k] = i
                        rsum += gam[[ind0[k], i, p]] * tc[[indk]]
                # loop on covariant indices:
                for k in range(n_con, tensor._tensor_rank):
                    for i in manif.irange():
                        indk = list(ind0)
                        indk[k] = i
                        rsum -= gam[[i, ind0[k], p]] * tc[[indk]]
                resc[[ind]] = rsum

        # Resulting tensor field
        if tensor._name is None:
            name_resu = None
        else:
            name_resu = self._name + '(' + tensor._name + ')'
        if tensor._latex_name is None:
            latex_name_resu = None
        else:
            latex_name_resu = format_unop_latex(self._latex_name + ' ',
                                                            tensor._latex_name)
        return tdom.vector_field_module().tensor_from_comp((n_con, n_cov+1),
                              resc, name=name_resu, latex_name=latex_name_resu)

    def torsion(self):
        r"""
        Return the connection's torsion tensor.

        The torsion tensor is the tensor field `T` of type (1,2) defined by

        .. MATH::

            T(\omega, u, v) = \left\langle \omega, \nabla_u v - \nabla_v u
                - [u, v] \right\rangle

        for any 1-form  `\omega`  and any vector fields `u` and `v`.

        OUTPUT:

        - the torsion tensor `T`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`

        EXAMPLES:

        Torsion of an affine connection on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: t = nab.torsion() ; t
            Tensor field of type (1,2) on the 3-dimensional differentiable
             manifold M
            sage: t.symmetries()
            no symmetry;  antisymmetry: (1, 2)
            sage: t[:]
            [[[0, -x^2, 0], [x^2, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, -y*z], [0, y*z, 0]]]

        The torsion expresses the lack of commutativity of two successive
        derivatives of a scalar field::

            sage: f = M.scalar_field(x*z^2 + y^2 - z^2, name='f')
            sage: DDf = nab(nab(f)) ; DDf
            Tensor field nabla(df) of type (0,2) on the 3-dimensional
             differentiable manifold M
            sage: DDf.antisymmetrize()[:]  # two successive derivatives do not commute:
            [             0   -1/2*x^2*z^2              0]
            [   1/2*x^2*z^2              0 -(x - 1)*y*z^2]
            [             0  (x - 1)*y*z^2              0]
            sage: 2*DDf.antisymmetrize() == nab.torsion().contract(0,nab(f))
            True

        The above identity is the standard formula

        .. MATH::

            \nabla_j \nabla_i \, f - \nabla_i \nabla_j \, f = T^k_{\ \, ij} \nabla_k \, f ,

        where the `T^k_{\ \, ij}`'s are the components of the torsion tensor.

        The result is cached::

            sage: nab.torsion() is t
            True

        as long as the connection remains unchanged::

            sage: nab[2,1,3] = 1+x    # changing the connection
            sage: nab.torsion() is t  # a new computation of the torsion has been made
            False
            sage: (nab.torsion() - t).display()
            (-x - 1) ∂/∂y⊗dx⊗dz + (x + 1) ∂/∂y⊗dz⊗dx

        Another example: torsion of some connection on a non-parallelizable
        2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: c_xyW = c_xy.restrict(W) ; c_uvW = c_uv.restrict(W)
            sage: eUW = c_xyW.frame() ; eVW = c_uvW.frame()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[0,0,0], nab[0,1,0], nab[1,0,1] = x, x-y, x*y
            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         for k in M.irange():
            ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()
            sage: t = nab.torsion() ; t
            Tensor field of type (1,2) on the 2-dimensional differentiable
             manifold M
            sage: t.parent()
            Module T^(1,2)(M) of type-(1,2) tensors fields on the 2-dimensional
             differentiable manifold M
            sage: t[eU,:]
            [[[0, x - y], [-x + y, 0]], [[0, -x*y], [x*y, 0]]]
            sage: t[eV,:]
            [[[0, 1/8*u^2 - 1/8*v^2 - 1/2*v], [-1/8*u^2 + 1/8*v^2 + 1/2*v, 0]],
             [[0, -1/8*u^2 + 1/8*v^2 - 1/2*v], [1/8*u^2 - 1/8*v^2 + 1/2*v, 0]]]

        Check of the torsion formula::

            sage: f = M.scalar_field({c_xy: (x+y)^2, c_uv: u^2}, name='f')
            sage: DDf = nab(nab(f)) ; DDf
            Tensor field nabla(df) of type (0,2) on the 2-dimensional
             differentiable manifold M
            sage: DDf.antisymmetrize().display(eU)
            (-x^2*y - (x + 1)*y^2 + x^2) dx∧dy
            sage: DDf.antisymmetrize().display(eV)
            (1/8*u^3 - 1/8*u*v^2 - 1/2*u*v) du∧dv
            sage: 2*DDf.antisymmetrize() == nab(f).contract(nab.torsion())
            True

        """
        if self._torsion is None:
            manif = self._domain
            resu = self._domain.tensor_field(1, 2, antisym=(1,2))
            for frame, gam in self._coefficients.items():
                sc = frame.structure_coeff()
                res = resu.add_comp(frame)
                for k in manif.irange():
                    for i in manif.irange():
                         for j in manif.irange(start=i+1):
                             res[[k,i,j]] = gam[[k,j,i]] - gam[[k,i,j]] - \
                                            sc[[k,i,j]]
            self._torsion = resu
        return self._torsion

    def riemann(self):
        r"""
        Return the connection's Riemann curvature tensor.

        The *Riemann curvature tensor* is the tensor field `R` of type (1,3)
        defined by

        .. MATH::

            R(\omega, w, u, v) = \left\langle \omega, \nabla_u \nabla_v w
                - \nabla_v \nabla_u w - \nabla_{[u, v]} w \right\rangle

        for any 1-form  `\omega`  and any vector fields `u`, `v` and `w`.

        OUTPUT:

        - the Riemann curvature tensor `R`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`

        EXAMPLES:

        Curvature of an affine connection on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla') ; nab
            Affine connection nabla on the 3-dimensional differentiable
             manifold M
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: r = nab.riemann() ; r
            Tensor field of type (1,3) on the 3-dimensional differentiable
             manifold M
            sage: r.parent()
            Free module T^(1,3)(M) of type-(1,3) tensors fields on the
             3-dimensional differentiable manifold M

        By construction, the Riemann tensor is antisymmetric with respect to
        its last two arguments (denoted `u` and `v` in the definition above),
        which are at positions 2 and 3 (the first argument being at position
        0)::

            sage: r.symmetries()
            no symmetry;  antisymmetry: (2, 3)

        The components::

            sage: r[:]
            [[[[0, 2*x, 0], [-2*x, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, z], [0, -z, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]]]

        The result is cached (until the connection is modified via
        :meth:`set_coef` or :meth:`add_coef`)::

            sage: nab.riemann() is r
            True

        Another example: Riemann curvature tensor of some connection on a
        non-parallelizable 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: c_xyW = c_xy.restrict(W) ; c_uvW = c_uv.restrict(W)
            sage: eUW = c_xyW.frame() ; eVW = c_uvW.frame()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[0,0,0], nab[0,1,0], nab[1,0,1] = x, x-y, x*y
            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         for k in M.irange():
            ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()
            sage: r = nab.riemann() ; r
            Tensor field of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: r.parent()
            Module T^(1,3)(M) of type-(1,3) tensors fields on the 2-dimensional
             differentiable manifold M
            sage: r.display(eU)
            (x^2*y - x*y^2) ∂/∂x⊗dx⊗dx⊗dy + (-x^2*y + x*y^2) ∂/∂x⊗dx⊗dy⊗dx + ∂/∂x⊗dy⊗dx⊗dy
             - ∂/∂x⊗dy⊗dy⊗dx - (x^2 - 1)*y ∂/∂y⊗dx⊗dx⊗dy + (x^2 - 1)*y ∂/∂y⊗dx⊗dy⊗dx
             + (-x^2*y + x*y^2) ∂/∂y⊗dy⊗dx⊗dy + (x^2*y - x*y^2) ∂/∂y⊗dy⊗dy⊗dx
            sage: r.display(eV)
            (1/32*u^3 - 1/32*u*v^2 - 1/32*v^3 + 1/32*(u^2 + 4)*v - 1/8*u - 1/4) ∂/∂u⊗du⊗du⊗dv
             + (-1/32*u^3 + 1/32*u*v^2 + 1/32*v^3 - 1/32*(u^2 + 4)*v + 1/8*u + 1/4) ∂/∂u⊗du⊗dv⊗du
             + (1/32*u^3 - 1/32*u*v^2 + 3/32*v^3 - 1/32*(3*u^2 - 4)*v - 1/8*u + 1/4) ∂/∂u⊗dv⊗du⊗dv
             + (-1/32*u^3 + 1/32*u*v^2 - 3/32*v^3 + 1/32*(3*u^2 - 4)*v + 1/8*u - 1/4) ∂/∂u⊗dv⊗dv⊗du
             + (-1/32*u^3 + 1/32*u*v^2 + 5/32*v^3 - 1/32*(5*u^2 + 4)*v + 1/8*u - 1/4) ∂/∂v⊗du⊗du⊗dv
             + (1/32*u^3 - 1/32*u*v^2 - 5/32*v^3 + 1/32*(5*u^2 + 4)*v - 1/8*u + 1/4) ∂/∂v⊗du⊗dv⊗du
             + (-1/32*u^3 + 1/32*u*v^2 + 1/32*v^3 - 1/32*(u^2 + 4)*v + 1/8*u + 1/4) ∂/∂v⊗dv⊗du⊗dv
             + (1/32*u^3 - 1/32*u*v^2 - 1/32*v^3 + 1/32*(u^2 + 4)*v - 1/8*u - 1/4) ∂/∂v⊗dv⊗dv⊗du

        The same computation parallelized on 2 cores::

            sage: Parallelism().set(nproc=2)
            sage: r_backup = r
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[0,0,0], nab[0,1,0], nab[1,0,1] = x, x-y, x*y
            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         for k in M.irange():
            ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()
            sage: r = nab.riemann() ; r
            Tensor field of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: r.parent()
            Module T^(1,3)(M) of type-(1,3) tensors fields on the 2-dimensional
             differentiable manifold M
            sage: r == r_backup
            True
            sage: Parallelism().set(nproc=1)  # switch off parallelization

        """
        if self._riemann is None:
            manif = self._domain
            resu = self._domain.tensor_field(1, 3, antisym=(2,3))
            for frame, gam in self._coefficients.items():
                # The computation is performed only on the top frames:
                for oframe in self._coefficients:
                    if frame in oframe._subframes and frame is not oframe:
                        break
                else:
                    # frame in not a subframe and the computation is performed:
                    sc = frame.structure_coeff()
                    gam_gam = gam.contract(1, gam, 0)
                    gam_sc = gam.contract(2, sc, 0)
                    res = resu.add_comp(frame)
                    if Parallelism().get('tensor') != 1:
                        # parallel computation
                        nproc = Parallelism().get('tensor')
                        lol = lambda lst, sz: [lst[i:i+sz] for i in range(0,
                                                                 len(lst), sz)]
                        ind_list = []
                        for i in manif.irange():
                            for j in manif.irange():
                                ind_list.append((i,j))
                        ind_step = max(1,int(len(ind_list)/nproc/2))
                        local_list = lol(ind_list,ind_step)
                        # definition of the list of input parameters
                        listParalInput = []
                        for ind_part in local_list:
                            listParalInput.append((frame,gam,gam_gam,gam_sc,
                                                   manif.irange,ind_part))

                        # definition of the parallel function
                        @parallel(p_iter='multiprocessing',ncpus=nproc)
                        def make_Reim(frame,gam,gam_gam,gam_sc,indices,
                                      local_list_ij):
                            partial = []
                            for i,j in local_list_ij:
                                for k in indices():
                                    for l in indices(start=k+1):
                                        partial.append([i,j,k,l,
                                                frame[k](gam[[i,j,l]]) - \
                                                frame[l](gam[[i,j,k]]) + \
                                                gam_gam[[i,k,j,l]] -  \
                                                gam_gam[[i,l,j,k]] -  \
                                                gam_sc[[i,j,k,l]]]
                                            )
                            return partial
                        # Computation and assignation of values
                        for ii,val in make_Reim(listParalInput):
                            for jj in val:
                                res[jj[0],jj[1],jj[2],jj[3]] = jj[4]

                    else:
                        # sequential
                        for i in manif.irange():
                            for j in manif.irange():
                                for k in manif.irange():
                                    # antisymmetry of the Riemann tensor taken
                                    # into account by l>k:
                                    for l in manif.irange(start=k+1):
                                        res[i,j,k,l] = frame[k](gam[[i,j,l]]) - \
                                                       frame[l](gam[[i,j,k]]) + \
                                                       gam_gam[[i,k,j,l]] -  \
                                                       gam_gam[[i,l,j,k]] -  \
                                                       gam_sc[[i,j,k,l]]
            self._riemann = resu
        return self._riemann

    def ricci(self):
        r"""
        Return the connection's Ricci tensor.

        The *Ricci tensor* is the tensor field `Ric` of type (0,2)
        defined from the Riemann curvature tensor `R` by

        .. MATH::

            Ric(u, v) = R(e^i, u, e_i, v)

        for any vector fields `u` and `v`, `(e_i)` being any vector frame and
        `(e^i)` the dual coframe.

        OUTPUT:

        - the Ricci  tensor `Ric`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`

        EXAMPLES:

        Ricci tensor of an affine connection on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla') ; nab
            Affine connection nabla on the 3-dimensional differentiable
             manifold M
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: r = nab.ricci() ; r
            Tensor field of type (0,2) on the 3-dimensional differentiable
             manifold M
            sage: r[:]
            [  0 2*x   0]
            [  0  -z   0]
            [  0   0   0]

        The result is cached (until the connection is modified via
        :meth:`set_coef` or :meth:`add_coef`)::

            sage: nab.ricci() is r
            True

        """
        if self._ricci is None:
            self._ricci = self.riemann().trace(0,2)
        return self._ricci

    def connection_form(self, i, j, frame=None):
        r"""
        Return the connection 1-form corresponding to the given index and
        vector frame.

        The *connection 1-forms* with respect to the frame `(e_i)` are the
        `n^2` 1-forms `\omega^i_{\ \, j}` defined by

        .. MATH::

            \nabla_v e_j = \langle \omega^i_{\ \, j}, v \rangle
                \, e_i

        for any vector `v`.

        The components of `\omega^i_{\ \, j}` in the coframe `(e^i)` dual to
        `(e_i)` are nothing but the connection coefficients `\Gamma^i_{\ \, jk}`
        relative to the frame `(e_i)`:

        .. MATH::

            \omega^i_{\ \, j} = \Gamma^i_{\ \, jk} e^k


        INPUT:

        - ``i``, ``j`` -- indices identifying the 1-form `\omega^i_{\ \, j}`
        - ``frame`` -- (default: ``None``) vector frame relative to which the
          connection 1-forms are defined; if ``None``, the default frame of the
          connection's domain is assumed.

        OUTPUT:

        - the 1-form `\omega^i_{\ \, j}`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES:

        Connection 1-forms on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.connection_form(1,1)  # connection 1-form (i,j)=(1,1) w.r.t. M's default frame
            1-form nabla connection 1-form (1,1) on the 3-dimensional
             differentiable manifold M
            sage: nab.connection_form(1,1)[:]
            [x*y*z, x^2, -y*z]

        The result is cached (until the connection is modified via
        :meth:`set_coef` or :meth:`add_coef`)::

            sage: nab.connection_form(1,1) is nab.connection_form(1,1)
            True

        Connection 1-forms w.r.t. a non-holonomic frame::

            sage: ch_basis = M.automorphism_field()
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = M.default_frame().new_frame(ch_basis, 'e')
            sage: e[1][:], e[2][:], e[3][:]
            ([y, 0, 0], [0, z, 0], [0, 0, x])
            sage: nab.connection_form(1,1,e)
            1-form nabla connection 1-form (1,1) on the 3-dimensional
             differentiable manifold M
            sage: nab.connection_form(1,1,e).comp(e)[:]
            [x*y^2*z, (x^2*y + 1)*z/y, -x*y*z]

        Check of the formula `\omega^i_{\ \, j} = \Gamma^i_{\ \, jk} e^k`:

        First on the manifold's default frame (∂/∂x, ∂/∂y, d:dz)::

            sage: dx = M.default_frame().coframe() ; dx
            Coordinate coframe (M, (dx,dy,dz))
            sage: check = []
            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         check.append( nab.connection_form(i,j) == \
            ....:               sum( nab[[i,j,k]]*dx[k] for k in M.irange() ) )
            sage: check
            [True, True, True, True, True, True, True, True, True]

        Then on the frame e::

            sage: ef = e.coframe() ; ef
            Coframe (M, (e^1,e^2,e^3))
            sage: check = []
            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         s = nab.connection_form(i,j,e).comp(c_xyz.frame(), from_basis=e)
            ....:         check.append( nab.connection_form(i,j,e) == sum( nab.coef(e)[[i,j,k]]*ef[k] for k in M.irange() ) )
            sage: check
            [True, True, True, True, True, True, True, True, True]

        Check of the formula
        `\nabla_v e_j = \langle \omega^i_{\ \, j}, v \rangle e_i`::

            sage: v = M.vector_field()
            sage: v[:] = (x*y, z^2-3*x, z+2*y)
            sage: b = M.default_frame()
            sage: for j in M.irange():  # check on M's default frame
            ....:     nab(b[j]).contract(v) == \
            ....:      sum( nab.connection_form(i,j)(v)*b[i] for i in M.irange())
            True
            True
            True
            sage: for j in M.irange():  # check on frame e
            ....:     nab(e[j]).contract(v) == \
            ....:      sum( nab.connection_form(i,j,e)(v)*e[i] for i in M.irange())
            True
            True
            True

        """
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._connection_forms:
            forms = {}
            frame_dom = frame.domain()
            coef_frame = self.coef(frame)
            for i1 in self._domain.irange():
                for j1 in self._domain.irange():
                    name = self._name + " connection 1-form (" + str(i1) + \
                           "," + str(j1) + ")"
                    latex_name = r"\omega^" + str(i1) + r"_{\ \, " + \
                                 str(j1) + "}"
                    omega = frame_dom.one_form(name=name,
                                               latex_name=latex_name)
                    comega = omega.set_comp(frame)
                    for k in self._domain.irange():
                        comega[k] = coef_frame[[i1,j1,k]]
                    forms[(i1,j1)] = omega
            self._connection_forms[frame] = forms
        return  self._connection_forms[frame][(i,j)]

    def torsion_form(self, i, frame=None):
        r"""
        Return the torsion 2-form corresponding to the given index and
        vector frame.

        The *torsion 2-forms* with respect to the frame `(e_i)` are the
        `n` 2-forms `\theta^i` defined by

        .. MATH::

            \theta^i(u,v) = T(e^i, u, v)

        where `T` is the connection's torsion tensor (cf. :meth:`torsion`),
        `(e^i)` is the coframe dual to `(e_i)` and `(u,v)` is a generic pair of
        vectors.

        INPUT:

        - ``i`` -- index identifying the 2-form `\theta^i`
        - ``frame`` -- (default: ``None``) vector frame relative to which the
          torsion 2-forms are defined; if ``None``, the default frame of the
          connection's domain is assumed.

        OUTPUT:

        - the 2-form `\theta^i`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES:

        Torsion 2-forms on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.torsion_form(1)
            2-form torsion (1) of connection nabla w.r.t. Coordinate frame
             (M, (∂/∂x,∂/∂y,∂/∂z)) on the 3-dimensional differentiable manifold M
            sage: nab.torsion_form(1)[:]
            [               0             -x^2      (y^2 + y)*z]
            [             x^2                0  x^3 - x^2 + y^2]
            [    -(y^2 + y)*z -x^3 + x^2 - y^2                0]

        Torsion 2-forms w.r.t. a non-holonomic frame::

            sage: ch_basis = M.automorphism_field()
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = M.default_frame().new_frame(ch_basis, 'e')
            sage: e[1][:], e[2][:], e[3][:]
            ([y, 0, 0], [0, z, 0], [0, 0, x])
            sage: ef = e.coframe()
            sage: ef[1][:], ef[2][:], ef[3][:]
            ([1/y, 0, 0], [0, 1/z, 0], [0, 0, 1/x])
            sage: nab.torsion_form(1, e)
            2-form torsion (1) of connection nabla w.r.t. Vector frame
             (M, (e_1,e_2,e_3)) on the 3-dimensional differentiable manifold M
            sage: nab.torsion_form(1, e).comp(e)[:]
            [                       0                   -x^2*z          (x*y^2 + x*y)*z]
            [                   x^2*z                        0  (x^4 - x^3 + x*y^2)*z/y]
            [        -(x*y^2 + x*y)*z -(x^4 - x^3 + x*y^2)*z/y                        0]

        Cartan's first structure equation is

        .. MATH::

            \theta^i = \mathrm{d} e^i + \omega^i_{\ \, j} \wedge e^j

        where the `\omega^i_{\ \, j}`'s are the connection 1-forms (cf.
        :meth:`connection_form`). Let us check it on the frame e::

            sage: for i in M.irange():  # long time
            ....:     nab.torsion_form(i, e) == ef[i].exterior_derivative() + \
            ....:      sum(nab.connection_form(i,j,e).wedge(ef[j]) for j in M.irange())
            True
            True
            True
        """
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._torsion_forms:
            forms = {}
            frame_dom = frame.domain()
            torsion_comp = self.torsion().comp(frame)
            for i1 in self._domain.irange():
                name = "torsion ({}) of connection ".format(i1) + \
                       self._name + " w.r.t. {}".format(frame)
                latex_name = r"\theta^" + str(i1)
                theta = frame_dom.diff_form(2, name=name,
                                            latex_name=latex_name)
                ctheta = theta.set_comp(frame)
                for k in self._domain.irange():
                    for l in self._domain.irange(start=k+1):
                        ctheta[k,l] = torsion_comp[[i1,k,l]]
                forms[i1] = theta
            self._torsion_forms[frame] = forms
        return  self._torsion_forms[frame][i]

    def curvature_form(self, i, j, frame=None):
        r"""
        Return the curvature 2-form corresponding to the given index and
        vector frame.

        The *curvature 2-forms* with respect to the frame `(e_i)` are the
        `n^2` 2-forms `\Omega^i_{\ \, j}` defined by

        .. MATH::

            \Omega^i_{\ \, j}(u,v) = R(e^i, e_j, u, v)

        where `R` is the connection's Riemann curvature tensor (cf.
        :meth:`riemann`), `(e^i)` is the coframe dual to `(e_i)` and `(u,v)` is
        a generic pair of vectors.

        INPUT:

        - ``i``, ``j`` -- indices identifying the 2-form `\Omega^i_{\ \, j}`
        - ``frame`` -- (default: ``None``) vector frame relative to which the
          curvature 2-forms are defined; if ``None``, the default frame
          of the connection's domain is assumed.

        OUTPUT:

        - the 2-form `\Omega^i_{\ \, j}`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES:

        Curvature 2-forms on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.curvature_form(1,1)  # long time
            2-form curvature (1,1) of connection nabla w.r.t. Coordinate frame
             (M, (∂/∂x,∂/∂y,∂/∂z)) on the 3-dimensional differentiable manifold M
            sage: nab.curvature_form(1,1).display()  # long time (if above is skipped)
            curvature (1,1) of connection nabla w.r.t. Coordinate frame
             (M, (∂/∂x,∂/∂y,∂/∂z)) = (y^2*z^3 + (x*y^3 - x)*z + 2*x) dx∧dy
              + (x^3*z^2 - x*y) dx∧dz + (x^4*y*z^2 - z) dy∧dz

        Curvature 2-forms w.r.t. a non-holonomic frame::

            sage: ch_basis = M.automorphism_field()
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = M.default_frame().new_frame(ch_basis, 'e')
            sage: e[1].display(), e[2].display(), e[3].display()
            (e_1 = y ∂/∂x, e_2 = z ∂/∂y, e_3 = x ∂/∂z)
            sage: ef = e.coframe()
            sage: ef[1].display(), ef[2].display(), ef[3].display()
            (e^1 = 1/y dx, e^2 = 1/z dy, e^3 = 1/x dz)
            sage: nab.curvature_form(1,1,e)  # long time
            2-form curvature (1,1) of connection nabla w.r.t. Vector frame
             (M, (e_1,e_2,e_3)) on the 3-dimensional differentiable manifold M
            sage: nab.curvature_form(1,1,e).display(e)  # long time (if above is skipped)
             curvature (1,1) of connection nabla w.r.t. Vector frame
             (M, (e_1,e_2,e_3)) =
              (y^3*z^4 + 2*x*y*z + (x*y^4 - x*y)*z^2) e^1∧e^2
              + (x^4*y*z^2 - x^2*y^2) e^1∧e^3 + (x^5*y*z^3 - x*z^2) e^2∧e^3

        Cartan's second structure equation is

        .. MATH::

            \Omega^i_{\ \, j} = \mathrm{d} \omega^i_{\ \, j}
                                + \omega^i_{\ \, k} \wedge \omega^k_{\ \, j}

        where the `\omega^i_{\ \, j}`'s are the connection 1-forms (cf.
        :meth:`connection_form`). Let us check it on the frame e::

            sage: omega = nab.connection_form
            sage: check = []
            sage: for i in M.irange():  # long time
            ....:     for j in M.irange():
            ....:         check.append( nab.curvature_form(i,j,e) == \
            ....:                       omega(i,j,e).exterior_derivative() + \
            ....:         sum( omega(i,k,e).wedge(omega(k,j,e)) for k in M.irange()) )
            sage: check  # long time
            [True, True, True, True, True, True, True, True, True]

        """
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._curvature_forms:
            forms = {}
            frame_dom = frame.domain()
            riemann_comp = self.riemann().comp(frame)
            for i1 in self._domain.irange():
                for j1 in self._domain.irange():
                    name = "curvature ({},{}) of connection ".format(i1,j1) + \
                           self._name + " w.r.t. {}".format(frame)
                    latex_name = r"\Omega^" + str(i1) + r"_{\ \, " + \
                                str(j1) + "}"
                    omega = frame_dom.diff_form(2, name=name,
                                                latex_name=latex_name)
                    comega = omega.set_comp(frame)
                    for k in self._domain.irange():
                        for l in self._domain.irange(start=k+1):
                            comega[k,l] = riemann_comp[[i1,j1,k,l]]
                    forms[(i1,j1)] = omega
            self._curvature_forms[frame] = forms
        return  self._curvature_forms[frame][(i,j)]

    def set_calc_order(self, symbol, order, truncate=False):
        r"""
        Trigger a series expansion with respect to a small parameter in
        computations involving ``self``.

        This property is propagated by usual operations. The internal
        representation must be ``SR`` for this to take effect.

        INPUT:

        - ``symbol`` -- symbolic variable (the "small parameter" `\epsilon`)
          with respect to which the connection coefficients are expanded in
          power series
        - ``order`` -- integer; the order `n` of the expansion, defined as the
          degree of the polynomial representing the truncated power series in
          ``symbol``
        - ``truncate`` -- (default: ``False``) determines whether the
          connection coefficients are replaced by their expansions to the
          given order

        EXAMPLES::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: C.<t,x,y,z> = M.chart()
            sage: e = var('e')
            sage: g = M.metric()
            sage: h = M.tensor_field(0, 2, sym=(0,1))
            sage: g[0, 0], g[1, 1], g[2, 2], g[3, 3] = -1, 1, 1, 1
            sage: h[0, 1] = x
            sage: g.set(g + e*h)
            sage: g[:]
            [ -1 e*x   0   0]
            [e*x   1   0   0]
            [  0   0   1   0]
            [  0   0   0   1]
            sage: nab = g.connection()
            sage: nab[0, 1, 1]
            -e/(e^2*x^2 + 1)
            sage: nab.set_calc_order(e, 1, truncate=True)
            sage: nab[0, 1, 1]
            -e

        """
        for coef in self._coefficients.values():
            for ind in coef.non_redundant_index_generator():
                coef[ind]._expansion_symbol = symbol
                coef[ind]._order = order
                if truncate:
                    coef[ind].simplify()
        self._del_derived()

    @cached_method
    def __hash__(self):
        r"""
        Hash function.

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: eX = X.frame()
            sage: nab1 = M.affine_connection('nabla1', latex_name=r'\nabla_1')
            sage: nab1.set_coef(eX)[1,2,1] = x*y
            sage: nab2 = M.affine_connection('nabla2', latex_name=r'\nabla_2')
            sage: nab2.set_coef(eX)[1,2,1] = x*y
            sage: nab1.set_immutable(); nab2.set_immutable()
            sage: nab1 == nab2
            True
            sage: hash(nab1) == hash(nab2)
            True

        Let us check that affine connections can be used as dictionary keys::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: eX = X.frame()
            sage: nab1 = M.affine_connection('nabla1', latex_name=r'\nabla_1')
            sage: nab1.set_coef(eX)[1,2,1] = x*y
            sage: nab2 = M.affine_connection('nabla2', latex_name=r'\nabla_2')
            sage: nab2.set_coef(eX)[1,2,1] = x^2
            sage: nab1.set_immutable(); nab2.set_immutable()
            sage: d = {nab1: 1, nab2: 2}
            sage: d[nab1]
            1
            sage: d[nab2]
            2

        """
        if self.is_mutable():
            raise ValueError('element must be immutable in order to be '
                             'hashable')
        return hash((type(self).__name__, self._domain))
