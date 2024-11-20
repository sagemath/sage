r"""
Pseudo-Riemannian Metrics and Degenerate Metrics

The class :class:`PseudoRiemannianMetric` implements pseudo-Riemannian metrics
on differentiable manifolds over `\RR`. The derived class
:class:`PseudoRiemannianMetricParal` is devoted to metrics with values on a
parallelizable manifold.

The class :class:`DegenerateMetric` implements degenerate (or null or lightlike)
metrics on differentiable manifolds over `\RR`. The derived class
:class:`DegenerateMetricParal` is devoted to metrics with values on a
parallelizable manifold.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Pablo Angulo (2016) : Schouten, Cotton and Cotton-York tensors
- Florentin Jaffredo (2018) : series expansion for the inverse metric
- Hans Fotsing Tetsing (2019) : degenerate metrics
- Marius Gerbershagen (2022) : compute volume forms with contravariant indices
  only as needed

REFERENCES:

- [KN1963]_
- [Lee1997]_
- [ONe1983]_
- [DB1996]_
- [DS2010]_
"""
# *****************************************************************************
#  Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#  Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#  Copyright (C) 2016 Pablo Angulo <pang@cancamusa.net>
#  Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#  Copyright (C) 2019 Hans Fotsing Tetsing <hans.fotsing@aims-cameroon.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import annotations

from typing import TYPE_CHECKING, overload

from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
from sage.rings.integer import Integer

if TYPE_CHECKING:
    from sage.manifolds.differentiable.diff_form import DiffForm


class PseudoRiemannianMetric(TensorField):
    r"""
    Pseudo-Riemannian metric with values on an open subset of a
    differentiable manifold.

    An instance of this class is a field of nondegenerate symmetric bilinear
    forms (metric field) along a differentiable manifold `U` with
    values on a differentiable manifold `M` over `\RR`, via a differentiable
    mapping `\Phi: U \rightarrow M`.
    The standard case of a metric field *on* a manifold corresponds to `U=M`
    and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval
    of `\RR`).

    A *metric* `g` is a field on `U`, such that at each point `p\in U`, `g(p)`
    is a bilinear map of the type:

    .. MATH::

        g(p):\ T_q M\times T_q M  \longrightarrow \RR

    where `T_q M` stands for the tangent space to the
    manifold `M` at the point `q=\Phi(p)`, such that `g(p)` is symmetric:
    `\forall (u,v)\in  T_q M\times T_q M, \ g(p)(v,u) = g(p)(u,v)`
    and nondegenerate:
    `(\forall v\in T_q M,\ \ g(p)(u,v) = 0) \Longrightarrow u=0`.

    .. NOTE::

        If `M` is parallelizable, the class :class:`PseudoRiemannianMetricParal`
        should be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `\Phi(U)\subset M`
    - ``name`` -- name given to the metric
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      single integer: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number
      of positive terms (resp. number of negative terms) in any diagonal
      writing of the metric components; if ``signature`` is ``None``, `S` is
      set to the dimension of manifold `M` (Riemannian signature)
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the metric;
      if ``None``, it is formed from ``name``

    EXAMPLES:

    Let us construct the standard metric on the sphere `S^2`, described in
    terms of stereographic coordinates, from the North pole (open subset `U`)
    and from the South pole (open subset `V`)::

        sage: M = Manifold(2, 'S^2', start_index=1)
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coord
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
        ....:                 restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: g = M.metric('g') ; g
        Riemannian metric g on the 2-dimensional differentiable manifold S^2

    The metric is considered as a tensor field of type (0,2) on `S^2`::

        sage: g.parent()
        Module T^(0,2)(S^2) of type-(0,2) tensors fields on the 2-dimensional
         differentiable manifold S^2

    We define `g` by its components on domain `U`::

        sage: g[eU,1,1], g[eU,2,2] = 4/(1+x^2+y^2)^2, 4/(1+x^2+y^2)^2
        sage: g.display(eU)
        g = 4/(x^2 + y^2 + 1)^2 dx⊗dx + 4/(x^2 + y^2 + 1)^2 dy⊗dy

    A matrix view of the components::

        sage: g[eU,:]
        [4/(x^2 + y^2 + 1)^2                   0]
        [                  0 4/(x^2 + y^2 + 1)^2]

    The components of `g` on domain `V` expressed in terms of coordinates
    `(u,v)` are obtained by applying (i) the tensor transformation law on
    `W = U\cap V` and (ii) some analytical continuation::

        sage: W = U.intersection(V)
        sage: g.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: g.apply_map(factor, frame=eV, keep_other_components=True) # for a nicer display
        sage: g.display(eV)
        g = 4/(u^2 + v^2 + 1)^2 du⊗du + 4/(u^2 + v^2 + 1)^2 dv⊗dv

    At this stage, the metric is fully defined on the whole sphere. Its
    restriction to some subdomain is itself a metric (by default, it bears the
    same symbol)::

        sage: g.restrict(U)
        Riemannian metric g on the Open subset U of the 2-dimensional
         differentiable manifold S^2
        sage: g.restrict(U).parent()
        Free module T^(0,2)(U) of type-(0,2) tensors fields on the Open subset
         U of the 2-dimensional differentiable manifold S^2

    The parent of `g|_U` is a free module because is `U` is a parallelizable
    domain, contrary to `S^2`. Actually, `g` and `g|_U` have different Python
    type::

        sage: type(g)
        <class 'sage.manifolds.differentiable.metric.PseudoRiemannianMetric'>
        sage: type(g.restrict(U))
        <class 'sage.manifolds.differentiable.metric.PseudoRiemannianMetricParal'>

    As a field of bilinear forms, the metric acts on pairs of vector fields,
    yielding a scalar field::

        sage: a = M.vector_field({eU: [x, 2+y]}, name='a')
        sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: b = M.vector_field({eU: [-y, x]}, name='b')
        sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: s = g(a,b) ; s
        Scalar field g(a,b) on the 2-dimensional differentiable manifold S^2
        sage: s.display()
        g(a,b): S^2 → ℝ
        on U: (x, y) ↦ 8*x/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1)
        on V: (u, v) ↦ 8*(u^3 + u*v^2)/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1)

    The inverse metric is::

        sage: ginv = g.inverse() ; ginv
        Tensor field inv_g of type (2,0) on the 2-dimensional differentiable
         manifold S^2
        sage: ginv.parent()
        Module T^(2,0)(S^2) of type-(2,0) tensors fields on the 2-dimensional
         differentiable manifold S^2
        sage: latex(ginv)
        g^{-1}
        sage: ginv.display(eU)
        inv_g = (1/4*x^4 + 1/4*y^4 + 1/2*(x^2 + 1)*y^2 + 1/2*x^2 + 1/4) ∂/∂x⊗∂/∂x
         + (1/4*x^4 + 1/4*y^4 + 1/2*(x^2 + 1)*y^2 + 1/2*x^2 + 1/4) ∂/∂y⊗∂/∂y
        sage: ginv.display(eV)
        inv_g = (1/4*u^4 + 1/4*v^4 + 1/2*(u^2 + 1)*v^2 + 1/2*u^2 + 1/4) ∂/∂u⊗∂/∂u
         + (1/4*u^4 + 1/4*v^4 + 1/2*(u^2 + 1)*v^2 + 1/2*u^2 + 1/4) ∂/∂v⊗∂/∂v

    We have::

        sage: ginv.restrict(U) is g.restrict(U).inverse()
        True
        sage: ginv.restrict(V) is g.restrict(V).inverse()
        True
        sage: ginv.restrict(W) is g.restrict(W).inverse()
        True

    To get the volume form (Levi-Civita tensor) associated with `g`, we have
    first to define an orientation on `S^2`. The standard orientation is that
    in which ``eV`` is right-handed; indeed, once supplemented by the outward
    unit normal, ``eV`` give birth to a right-handed frame with respect to the
    standard orientation of the ambient Euclidean space `E^3`. With such an
    orientation, ``eU`` is then left-handed and in order to define an
    orientation on the whole of `S^2`, we introduce a vector frame
    on `U` by swapping ``eU``'s vectors::

        sage: f = U.vector_frame('f', (eU[2], eU[1]))
        sage: M.set_orientation([eV, f])

    We have then, factorizing the components for a nicer display::

        sage: eps = g.volume_form() ; eps
        2-form eps_g on the 2-dimensional differentiable manifold S^2
        sage: eps.apply_map(factor, frame=eU, keep_other_components=True)
        sage: eps.apply_map(factor, frame=eV, keep_other_components=True)
        sage: eps.display(eU)
        eps_g = -4/(x^2 + y^2 + 1)^2 dx∧dy
        sage: eps.display(eV)
        eps_g = 4/(u^2 + v^2 + 1)^2 du∧dv

    The unique non-trivial component of the volume form is, up to a sign
    depending of the chosen orientation, nothing but the square root of the
    determinant of `g` in the corresponding frame::

        sage: eps[[eU,1,2]] == -g.sqrt_abs_det(eU)
        True
        sage: eps[[eV,1,2]] == g.sqrt_abs_det(eV)
        True

    The Levi-Civita connection associated with the metric `g`::

        sage: nabla = g.connection() ; nabla
        Levi-Civita connection nabla_g associated with the Riemannian metric g
         on the 2-dimensional differentiable manifold S^2
        sage: latex(nabla)
        \nabla_{g}

    The Christoffel symbols `\Gamma^i_{\ \, jk}` associated with some
    coordinates::

        sage: g.christoffel_symbols(c_xy)
        3-indices components w.r.t. Coordinate frame (U, (∂/∂x,∂/∂y)), with
         symmetry on the index positions (1, 2)
        sage: g.christoffel_symbols(c_xy)[:]
        [[[-2*x/(x^2 + y^2 + 1), -2*y/(x^2 + y^2 + 1)],
          [-2*y/(x^2 + y^2 + 1), 2*x/(x^2 + y^2 + 1)]],
         [[2*y/(x^2 + y^2 + 1), -2*x/(x^2 + y^2 + 1)],
          [-2*x/(x^2 + y^2 + 1), -2*y/(x^2 + y^2 + 1)]]]
        sage: g.christoffel_symbols(c_uv)[:]
        [[[-2*u/(u^2 + v^2 + 1), -2*v/(u^2 + v^2 + 1)],
          [-2*v/(u^2 + v^2 + 1), 2*u/(u^2 + v^2 + 1)]],
         [[2*v/(u^2 + v^2 + 1), -2*u/(u^2 + v^2 + 1)],
          [-2*u/(u^2 + v^2 + 1), -2*v/(u^2 + v^2 + 1)]]]

    The Christoffel symbols are nothing but the connection coefficients w.r.t.
    the coordinate frame::

        sage: g.christoffel_symbols(c_xy) is nabla.coef(c_xy.frame())
        True
        sage: g.christoffel_symbols(c_uv) is nabla.coef(c_uv.frame())
        True

    Test that `\nabla` is the connection compatible with `g`::

        sage: t = nabla(g) ; t
        Tensor field nabla_g(g) of type (0,3) on the 2-dimensional
         differentiable manifold S^2
        sage: t.display(eU)
        nabla_g(g) = 0
        sage: t.display(eV)
        nabla_g(g) = 0
        sage: t == 0
        True

    The Riemann curvature tensor of `g`::

        sage: riem = g.riemann() ; riem
        Tensor field Riem(g) of type (1,3) on the 2-dimensional differentiable
         manifold S^2
        sage: riem.display(eU)
        Riem(g) = 4/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) ∂/∂x⊗dy⊗dx⊗dy
         - 4/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) ∂/∂x⊗dy⊗dy⊗dx
         - 4/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) ∂/∂y⊗dx⊗dx⊗dy
         + 4/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) ∂/∂y⊗dx⊗dy⊗dx
        sage: riem.display(eV)
        Riem(g) = 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) ∂/∂u⊗dv⊗du⊗dv
         - 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) ∂/∂u⊗dv⊗dv⊗du
         - 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) ∂/∂v⊗du⊗du⊗dv
         + 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) ∂/∂v⊗du⊗dv⊗du

    The Ricci tensor of `g`::

        sage: ric = g.ricci() ; ric
        Field of symmetric bilinear forms Ric(g) on the 2-dimensional
         differentiable manifold S^2
        sage: ric.display(eU)
        Ric(g) = 4/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) dx⊗dx
         + 4/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) dy⊗dy
        sage: ric.display(eV)
        Ric(g) = 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) du⊗du
         + 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) dv⊗dv
        sage: ric == g
        True

    The Ricci scalar of `g`::

        sage: r = g.ricci_scalar() ; r
        Scalar field r(g) on the 2-dimensional differentiable manifold S^2
        sage: r.display()
        r(g): S^2 → ℝ
        on U: (x, y) ↦ 2
        on V: (u, v) ↦ 2

    In dimension 2, the Riemann tensor can be expressed entirely in terms of
    the Ricci scalar `r`:

    .. MATH::

        R^i_{\ \, jlk} = \frac{r}{2} \left( \delta^i_{\ \, k} g_{jl}
            - \delta^i_{\ \, l} g_{jk} \right)

    This formula can be checked here, with the r.h.s. rewritten as
    `-r g_{j[k} \delta^i_{\ \, l]}`::

        sage: delta = M.tangent_identity_field()
        sage: riem == - r*(g*delta).antisymmetrize(2,3)
        True
    """
    _derived_objects = ('_connection', '_ricci_scalar', '_weyl',
                       '_schouten', '_cotton', '_cotton_york')

    def __init__(self, vector_field_module, name, signature=None,
                 latex_name=None):
        r"""
        Construct a metric.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:               intersection_name='W', restrictions1= x>0,
            ....:               restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: XM = M.vector_field_module()
            sage: from sage.manifolds.differentiable.metric import \
            ....:                                        PseudoRiemannianMetric
            sage: g = PseudoRiemannianMetric(XM, 'g', signature=0); g
            Lorentzian metric g on the 2-dimensional differentiable
             manifold M
            sage: g[e_xy,0,0], g[e_xy,1,1] = -(1+x^2), 1+y^2
            sage: g.add_comp_by_continuation(e_uv, W, c_uv)
            sage: TestSuite(g).run(skip=['_test_category', '_test_pickling'])

        .. TODO::

            - fix _test_pickling (in the superclass TensorField)
            - add a specific parent to the metrics, to fit with the category
              framework
        """
        TensorField.__init__(self, vector_field_module, (0,2),
                             name=name, latex_name=latex_name, sym=(0,1))
        # signature:
        ndim = self._ambient_domain.dimension()
        if signature is None:
            signature = ndim
        else:
            if not isinstance(signature, (int, Integer)):
                raise TypeError("the metric signature must be an integer")
            if (signature < - ndim) or (signature > ndim):
                raise ValueError("metric signature out of range")
            if (signature+ndim) % 2 == 1:
                if ndim % 2 == 0:
                    raise ValueError("the metric signature must be even")
                else:
                    raise ValueError("the metric signature must be odd")
        self._signature = signature
        # the pair (n_+, n_-):
        self._signature_pm = ((ndim+signature)//2, (ndim-signature)//2)
        self._indic_signat = 1 - 2*(self._signature_pm[1] % 2)  # (-1)^n_-
        # Initialization of derived quantities:
        PseudoRiemannianMetric._init_derived(self)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g')
            sage: g._repr_()
            'Riemannian metric g on the 5-dimensional differentiable manifold M'
            sage: g = M.metric('g', signature=3)
            sage: g._repr_()
            'Lorentzian metric g on the 5-dimensional differentiable manifold M'
            sage: g = M.metric('g', signature=1)
            sage: g._repr_()
            'Pseudo-Riemannian metric g on the 5-dimensional differentiable manifold M'
        """
        n = self._ambient_domain.dimension()
        s = self._signature
        if s == n:
            description = "Riemannian metric "
        elif s == n-2 or s == 2-n:
            description = "Lorentzian metric "
        else:
            description = "Pseudo-Riemannian metric "
        description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` with the same
        signature.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g', signature=3)
            sage: g1 = g._new_instance(); g1
            Lorentzian metric unnamed metric on the 5-dimensional
             differentiable manifold M
            sage: type(g1) == type(g)
            True
            sage: g1.parent() is g.parent()
            True
            sage: g1.signature() == g.signature()
            True
        """
        return type(self)(self._vmodule, 'unnamed metric',
                          signature=self._signature,
                          latex_name=r'\text{unnamed metric}')

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g')
            sage: g._init_derived()
        """
        # Initialization of quantities pertaining to the mother class:
        TensorField._init_derived(self)
        # inverse metric:
        inv_name = 'inv_' + self._name
        inv_latex_name = self._latex_name + r'^{-1}'
        self._inverse = self._vmodule.tensor((2,0), name=inv_name,
                                             latex_name=inv_latex_name,
                                             sym=(0,1))
        for attr in self._derived_objects:
            self.__setattr__(attr, None)
        self._determinants = {} # determinants in various frames
        self._sqrt_abs_dets = {} # sqrt(abs(det g)) in various frames
        self._vol_forms = [] # volume form and associated tensors

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g')
            sage: g._del_derived()
        """
        # First the derived quantities from the mother class are deleted:
        TensorField._del_derived(self)
        # The inverse metric is cleared:
        self._del_inverse()
        # The connection, Ricci scalar and Weyl tensor are reset to None:
        # The Schouten, Cotton and Cotton-York tensors are reset to None:
        for attr in self._derived_objects:
            self.__setattr__(attr, None)
        # The dictionary of determinants over the various frames is cleared:
        self._determinants.clear()
        self._sqrt_abs_dets.clear()
        # The volume form and the associated tensors is deleted:
        del self._vol_forms[:]

    def _del_inverse(self):
        r"""
        Delete the inverse metric.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g')
            sage: g._del_inverse()
        """
        self._inverse._restrictions.clear()
        self._inverse._del_derived()

    def signature(self):
        r"""
        Signature of the metric.

        OUTPUT:

        - signature `S` of the metric, defined as the integer
          `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number of
          positive terms (resp. number of negative terms) in any diagonal
          writing of the metric components

        EXAMPLES:

        Signatures on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: g = M.metric('g') # if not specified, the signature is Riemannian
            sage: g.signature()
            2
            sage: h = M.metric('h', signature=0)
            sage: h.signature()
            0
        """
        return self._signature

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of the metric to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of the metric's domain (must be an
          instance of :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `V` is a subdomain of
          ``self._codomain``
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          If ``None``, the restriction of ``self._vmodule._dest_map`` to `U` is
          used.

        OUTPUT:

        - instance of :class:`PseudoRiemannianMetric` representing the
          restriction.

        EXAMPLES::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g', signature=3)
            sage: U = M.open_subset('U')
            sage: g.restrict(U)
            Lorentzian metric g on the Open subset U of the
             5-dimensional differentiable manifold M
            sage: g.restrict(U).signature()
            3

        See the top documentation of :class:`PseudoRiemannianMetric` for more
        examples.
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            # Construct the restriction at the tensor field level:
            resu = TensorField.restrict(self, subdomain, dest_map=dest_map)
            # the type is correctly handled by TensorField.restrict, i.e.
            # resu is of type self.__class__, but the signature is not handled
            # by TensorField.restrict; we have to set it here:
            resu._signature = self._signature
            resu._signature_pm = self._signature_pm
            resu._indic_signat = self._indic_signat
            # Restrictions of derived quantities:
            resu._inverse = self.inverse().restrict(subdomain)
            for attr in self._derived_objects:
                derived = self.__getattribute__(attr)
                if derived is not None:
                    resu.__setattr__(attr, derived.restrict(subdomain))
            if self._vol_forms != []:
                for eps in self._vol_forms:
                    resu._vol_forms.append(eps.restrict(subdomain))
            # NB: no initialization of resu._determinants nor
            # resu._sqrt_abs_dets
            # The restriction is ready:
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def set(self, symbiform):
        r"""
        Define the metric from a field of symmetric bilinear forms.

        INPUT:

        - ``symbiform`` -- instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          representing a field of symmetric bilinear forms

        EXAMPLES:

        Metric defined from a field of symmetric bilinear forms on a
        non-parallelizable 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: h = M.sym_bilin_form_field(name='h')
            sage: h[eU,0,0], h[eU,0,1], h[eU,1,1] = 1+x, x*y, 1-y
            sage: h.add_comp_by_continuation(eV, W, c_uv)
            sage: h.display(eU)
            h = (x + 1) dx⊗dx + x*y dx⊗dy + x*y dy⊗dx + (-y + 1) dy⊗dy
            sage: h.display(eV)
            h = (1/8*u^2 - 1/8*v^2 + 1/4*v + 1/2) du⊗du + 1/4*u du⊗dv
             + 1/4*u dv⊗du + (-1/8*u^2 + 1/8*v^2 + 1/4*v + 1/2) dv⊗dv
            sage: g = M.metric('g')
            sage: g.set(h)
            sage: g.display(eU)
            g = (x + 1) dx⊗dx + x*y dx⊗dy + x*y dy⊗dx + (-y + 1) dy⊗dy
            sage: g.display(eV)
            g = (1/8*u^2 - 1/8*v^2 + 1/4*v + 1/2) du⊗du + 1/4*u du⊗dv
             + 1/4*u dv⊗du + (-1/8*u^2 + 1/8*v^2 + 1/4*v + 1/2) dv⊗dv
        """
        if not isinstance(symbiform, TensorField):
            raise TypeError("the argument must be a tensor field")
        if symbiform._tensor_type != (0,2):
            raise TypeError("the argument must be of tensor type (0,2)")
        if symbiform._sym != ((0,1),):
            raise TypeError("the argument must be symmetric")
        if not symbiform._domain.is_subset(self._domain):
            raise TypeError("the symmetric bilinear form is not defined " +
                            "on the metric domain")
        self._del_derived()
        self._restrictions.clear()
        if isinstance(symbiform, TensorFieldParal):
            rst = self.restrict(symbiform._domain)
            rst.set(symbiform)
        else:
            for dom, symbiform_rst in symbiform._restrictions.items():
                rst = self.restrict(dom)
                rst.set(symbiform_rst)

    def inverse(self, expansion_symbol=None, order=1):
        r"""
        Return the inverse metric.

        INPUT:

        - ``expansion_symbol`` -- (default: ``None``) symbolic variable; if
          specified, the inverse will be expanded in power series with respect
          to this variable (around its zero value)
        - ``order`` -- integer (default: 1); the order of the expansion
          if ``expansion_symbol`` is not ``None``; the *order* is defined as
          the degree of the polynomial representing the truncated power series
          in ``expansion_symbol``; currently only first order inverse is
          supported

        If ``expansion_symbol`` is set, then the zeroth order metric must be
        invertible. Moreover, subsequent calls to this method will return
        a cached value, even when called with the default value (to enable
        computation of derived quantities). To reset, use :meth:`_del_derived`.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          with ``tensor_type`` = (2,0) representing the inverse metric

        EXAMPLES:

        Inverse of the standard metric on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)  # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coord.
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)  # the complement of the two poles
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: g = M.metric('g')
            sage: g[eU,1,1], g[eU,2,2] = 4/(1+x^2+y^2)^2, 4/(1+x^2+y^2)^2
            sage: g.add_comp_by_continuation(eV, W, c_uv)
            sage: ginv = g.inverse(); ginv
            Tensor field inv_g of type (2,0) on the 2-dimensional differentiable manifold S^2
            sage: ginv.display(eU)
            inv_g = (1/4*x^4 + 1/4*y^4 + 1/2*(x^2 + 1)*y^2 + 1/2*x^2 + 1/4) ∂/∂x⊗∂/∂x
             + (1/4*x^4 + 1/4*y^4 + 1/2*(x^2 + 1)*y^2 + 1/2*x^2 + 1/4) ∂/∂y⊗∂/∂y
            sage: ginv.display(eV)
            inv_g = (1/4*u^4 + 1/4*v^4 + 1/2*(u^2 + 1)*v^2 + 1/2*u^2 + 1/4) ∂/∂u⊗∂/∂u
             + (1/4*u^4 + 1/4*v^4 + 1/2*(u^2 + 1)*v^2 + 1/2*u^2 + 1/4) ∂/∂v⊗∂/∂v

        Let us check that ``ginv`` is indeed the inverse of ``g``::

            sage: s = g.contract(ginv); s  # contraction of last index of g with first index of ginv
            Tensor field of type (1,1) on the 2-dimensional differentiable manifold S^2
            sage: s == M.tangent_identity_field()
            True
        """
        # Is the inverse metric up to date?
        for dom, rst in self._restrictions.items():
            self._inverse._restrictions[dom] = rst.inverse(
                                             expansion_symbol=expansion_symbol,
                                             order=order) # forces the update
                                                          # of the restriction
        return self._inverse

    def connection(self, name=None, latex_name=None, init_coef=True):
        r"""
        Return the unique torsion-free affine connection compatible with
        ``self``.

        This is the so-called Levi-Civita connection.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Levi-Civita
          connection; if ``None``, it is formed from the metric name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Levi-Civita connection; if ``None``, it is set to ``name``, or if the
          latter is None as well, it formed from the symbol `\nabla` and the
          metric symbol
        - ``init_coef`` -- boolean (default: ``True``); determines whether the
          connection coefficients are initialized, as Christoffel symbols
          in the top charts of the domain of ``self`` (i.e. disregarding
          the subcharts)

        OUTPUT:

        - the Levi-Civita connection, as an instance of
          :class:`~sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection`

        EXAMPLES:

        Levi-Civita connection associated with the Euclidean metric on
        `\RR^3`::

            sage: M = Manifold(3, 'R^3', start_index=1)

        Let us use spherical coordinates on `\RR^3`::

            sage: U = M.open_subset('U') # the complement of the half-plane (y=0, x>=0)
            sage: c_spher.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2 , (r*sin(th))^2  # the Euclidean metric
            sage: g.connection()
            Levi-Civita connection nabla_g associated with the Riemannian
             metric g on the Open subset U of the 3-dimensional differentiable
             manifold R^3
            sage: g.connection().display()  # Nonzero connection coefficients
            Gam^r_th,th = -r
            Gam^r_ph,ph = -r*sin(th)^2
            Gam^th_r,th = 1/r
            Gam^th_th,r = 1/r
            Gam^th_ph,ph = -cos(th)*sin(th)
            Gam^ph_r,ph = 1/r
            Gam^ph_th,ph = cos(th)/sin(th)
            Gam^ph_ph,r = 1/r
            Gam^ph_ph,th = cos(th)/sin(th)

        Test of compatibility with the metric::

            sage: Dg = g.connection()(g) ; Dg
            Tensor field nabla_g(g) of type (0,3) on the Open subset U of the
             3-dimensional differentiable manifold R^3
            sage: Dg == 0
            True
            sage: Dig = g.connection()(g.inverse()) ; Dig
            Tensor field nabla_g(inv_g) of type (2,1) on the Open subset U of
             the 3-dimensional differentiable manifold R^3
            sage: Dig == 0
            True
        """
        from sage.manifolds.differentiable.levi_civita_connection import \
                                                           LeviCivitaConnection
        if self._connection is None:
            if latex_name is None:
                if name is None:
                    latex_name = r'\nabla_{' + self._latex_name + '}'
                else:
                    latex_name = name
            if name is None:
                name = 'nabla_' + self._name
            self._connection = LeviCivitaConnection(self, name,
                                                    latex_name=latex_name,
                                                    init_coef=init_coef)
        return self._connection

    def christoffel_symbols(self, chart=None):
        r"""
        Christoffel symbols of ``self`` with respect to a chart.

        INPUT:

        - ``chart`` -- (default: ``None``) chart with respect to which the
          Christoffel symbols are required; if none is provided, the
          default chart of the metric's domain is assumed.

        OUTPUT:

        - the set of Christoffel symbols in the given chart, as an instance of
          :class:`~sage.tensor.modules.comp.CompWithSym`

        EXAMPLES:

        Christoffel symbols of the flat metric on `\RR^3` with respect to
        spherical coordinates::

            sage: M = Manifold(3, 'R3', r'\RR^3', start_index=1)
            sage: U = M.open_subset('U') # the complement of the half-plane (y=0, x>=0)
            sage: X.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2, r^2*sin(th)^2
            sage: g.display()  # the standard flat metric expressed in spherical coordinates
            g = dr⊗dr + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph
            sage: Gam = g.christoffel_symbols() ; Gam
            3-indices components w.r.t. Coordinate frame (U, (∂/∂r,∂/∂th,∂/∂ph)),
             with symmetry on the index positions (1, 2)
            sage: type(Gam)
            <class 'sage.tensor.modules.comp.CompWithSym'>
            sage: Gam[:]
            [[[0, 0, 0], [0, -r, 0], [0, 0, -r*sin(th)^2]],
            [[0, 1/r, 0], [1/r, 0, 0], [0, 0, -cos(th)*sin(th)]],
            [[0, 0, 1/r], [0, 0, cos(th)/sin(th)], [1/r, cos(th)/sin(th), 0]]]
            sage: Gam[1,2,2]
            -r
            sage: Gam[2,1,2]
            1/r
            sage: Gam[3,1,3]
            1/r
            sage: Gam[3,2,3]
            cos(th)/sin(th)
            sage: Gam[2,3,3]
            -cos(th)*sin(th)

        Note that a better display of the Christoffel symbols is provided by
        the method :meth:`christoffel_symbols_display`::

            sage: g.christoffel_symbols_display()
            Gam^r_th,th = -r
            Gam^r_ph,ph = -r*sin(th)^2
            Gam^th_r,th = 1/r
            Gam^th_ph,ph = -cos(th)*sin(th)
            Gam^ph_r,ph = 1/r
            Gam^ph_th,ph = cos(th)/sin(th)
        """
        if chart is None:
            frame = self._domain._def_chart._frame
        else:
            frame = chart._frame
        return self.connection().coef(frame)

    def christoffel_symbols_display(self, chart=None, symbol=None,
                latex_symbol=None, index_labels=None, index_latex_labels=None,
                coordinate_labels=True, only_nonzero=True,
                only_nonredundant=True):
        r"""
        Display the Christoffel symbols w.r.t. to a given chart, one
        per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``chart`` -- (default: ``None``) chart with respect to which the
          Christoffel symbols are defined; if none is provided, the
          default chart of the metric's domain is assumed.
        - ``symbol`` -- (default: ``None``) string specifying the
          symbol of the connection coefficients; if ``None``, 'Gam' is used
        - ``latex_symbol`` -- (default: ``None``) string specifying the LaTeX
          symbol for the components; if ``None``, '\\Gamma' is used
        - ``index_labels`` -- (default: ``None``) list of strings representing
          the labels of each index; if ``None``, coordinate symbols are used
          except if ``coordinate_symbols`` is set to ``False``, in which case
          integer labels are used
        - ``index_latex_labels`` -- (default: ``None``) list of strings
          representing the LaTeX labels of each index; if ``None``, coordinate
          LaTeX symbols are used, except if ``coordinate_symbols`` is set to
          ``False``, in which case integer labels are used
        - ``coordinate_labels`` -- boolean (default: ``True``); if ``True``,
          coordinate symbols are used by default (instead of integers)
        - ``only_nonzero`` -- boolean (default: ``True``); if ``True``, only
          nonzero connection coefficients are displayed
        - ``only_nonredundant`` -- boolean (default: ``True``); if ``True``,
          only nonredundant (w.r.t. the symmetry of the last two indices)
          connection coefficients are displayed

        EXAMPLES:

        Christoffel symbols of the flat metric on `\RR^3` with respect to
        spherical coordinates::

            sage: M = Manifold(3, 'R3', r'\RR^3', start_index=1)
            sage: U = M.open_subset('U') # the complement of the half-plane (y=0, x>=0)
            sage: X.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2, r^2*sin(th)^2
            sage: g.display()  # the standard flat metric expressed in spherical coordinates
            g = dr⊗dr + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph
            sage: g.christoffel_symbols_display()
            Gam^r_th,th = -r
            Gam^r_ph,ph = -r*sin(th)^2
            Gam^th_r,th = 1/r
            Gam^th_ph,ph = -cos(th)*sin(th)
            Gam^ph_r,ph = 1/r
            Gam^ph_th,ph = cos(th)/sin(th)

        To list all nonzero Christoffel symbols, including those that can be
        deduced by symmetry, use ``only_nonredundant=False``::

            sage: g.christoffel_symbols_display(only_nonredundant=False)
            Gam^r_th,th = -r
            Gam^r_ph,ph = -r*sin(th)^2
            Gam^th_r,th = 1/r
            Gam^th_th,r = 1/r
            Gam^th_ph,ph = -cos(th)*sin(th)
            Gam^ph_r,ph = 1/r
            Gam^ph_th,ph = cos(th)/sin(th)
            Gam^ph_ph,r = 1/r
            Gam^ph_ph,th = cos(th)/sin(th)

        Listing all Christoffel symbols (except those that can be deduced by
        symmetry), including the vanishing one::

            sage: g.christoffel_symbols_display(only_nonzero=False)
            Gam^r_r,r = 0
            Gam^r_r,th = 0
            Gam^r_r,ph = 0
            Gam^r_th,th = -r
            Gam^r_th,ph = 0
            Gam^r_ph,ph = -r*sin(th)^2
            Gam^th_r,r = 0
            Gam^th_r,th = 1/r
            Gam^th_r,ph = 0
            Gam^th_th,th = 0
            Gam^th_th,ph = 0
            Gam^th_ph,ph = -cos(th)*sin(th)
            Gam^ph_r,r = 0
            Gam^ph_r,th = 0
            Gam^ph_r,ph = 1/r
            Gam^ph_th,th = 0
            Gam^ph_th,ph = cos(th)/sin(th)
            Gam^ph_ph,ph = 0

        Using integer labels::

            sage: g.christoffel_symbols_display(coordinate_labels=False)
            Gam^1_22 = -r
            Gam^1_33 = -r*sin(th)^2
            Gam^2_12 = 1/r
            Gam^2_33 = -cos(th)*sin(th)
            Gam^3_13 = 1/r
            Gam^3_23 = cos(th)/sin(th)
        """
        if chart is None:
            chart = self._domain.default_chart()
        return self.connection().display(frame=chart.frame(), chart=chart,
              symbol=symbol, latex_symbol=latex_symbol,
              index_labels=index_labels, index_latex_labels=index_latex_labels,
              coordinate_labels=coordinate_labels, only_nonzero=only_nonzero,
              only_nonredundant=only_nonredundant)

    def riemann(self, name=None, latex_name=None):
        r"""
        Return the Riemann curvature tensor associated with the metric.

        This method is actually a shortcut for ``self.connection().riemann()``

        The Riemann curvature tensor is the tensor field `R` of type (1,3)
        defined by

        .. MATH::

            R(\omega, u, v, w) = \left\langle \omega, \nabla_u \nabla_v w
                - \nabla_v \nabla_u w - \nabla_{[u, v]} w \right\rangle

        for any 1-form  `\omega`  and any vector fields `u`, `v` and `w`.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Riemann tensor;
          if none, it is set to "Riem(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Riemann tensor; if none, it is set to "\\mathrm{Riem}(g)", where "g"
          is the metric's name

        OUTPUT:

        - the Riemann curvature tensor `R`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`

        EXAMPLES:

        Riemann tensor of the standard metric on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of a meridian (domain of spherical coordinates)
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: a = var('a') # the sphere radius
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2] = a^2, a^2*sin(th)^2
            sage: g.display() # standard metric on the 2-sphere of radius a:
            g = a^2 dth⊗dth + a^2*sin(th)^2 dph⊗dph
            sage: g.riemann()
            Tensor field Riem(g) of type (1,3) on the Open subset U of the
             2-dimensional differentiable manifold S^2
            sage: g.riemann()[:]
            [[[[0, 0], [0, 0]], [[0, sin(th)^2], [-sin(th)^2, 0]]],
             [[[0, -1], [1, 0]], [[0, 0], [0, 0]]]]

        In dimension 2, the Riemann tensor can be expressed entirely in terms of
        the Ricci scalar `r`:

        .. MATH::

            R^i_{\ \, jlk} = \frac{r}{2} \left( \delta^i_{\ \, k} g_{jl}
                - \delta^i_{\ \, l} g_{jk} \right)

        This formula can be checked here, with the r.h.s. rewritten as
        `-r g_{j[k} \delta^i_{\ \, l]}`::

            sage: g.riemann() == \
            ....:  -g.ricci_scalar()*(g*U.tangent_identity_field()).antisymmetrize(2,3)
            True

        Using SymPy as symbolic engine::

            sage: M.set_calculus_method('sympy')
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2] = a**2, a**2*sin(th)**2
            sage: g.riemann()[:]
            [[[[0, 0], [0, 0]],
              [[0, sin(2*th)/(2*tan(th)) - cos(2*th)],
               [-sin(2*th)/(2*tan(th)) + cos(2*th), 0]]],
             [[[0, -1], [1, 0]], [[0, 0], [0, 0]]]]
        """
        return self.connection().riemann(name, latex_name)

    def ricci(self, name=None, latex_name=None):
        r"""
        Return the Ricci tensor associated with the metric.

        This method is actually a shortcut for ``self.connection().ricci()``

        The Ricci tensor is the tensor field `Ric` of type (0,2)
        defined from the Riemann curvature tensor `R` by

        .. MATH::

            Ric(u, v) = R(e^i, u, e_i, v)

        for any vector fields `u` and `v`, `(e_i)` being any vector frame and
        `(e^i)` the dual coframe.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Ricci tensor;
          if none, it is set to "Ric(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Ricci tensor; if none, it is set to "\\mathrm{Ric}(g)", where "g"
          is the metric's name

        OUTPUT:

        - the Ricci tensor `Ric`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField` of tensor
          type (0,2) and symmetric

        EXAMPLES:

        Ricci tensor of the standard metric on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of a meridian (domain of spherical coordinates)
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: a = var('a') # the sphere radius
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2] = a^2, a^2*sin(th)^2
            sage: g.display() # standard metric on the 2-sphere of radius a:
            g = a^2 dth⊗dth + a^2*sin(th)^2 dph⊗dph
            sage: g.ricci()
            Field of symmetric bilinear forms Ric(g) on the Open subset U of
             the 2-dimensional differentiable manifold S^2
            sage: g.ricci()[:]
            [        1         0]
            [        0 sin(th)^2]
            sage: g.ricci() == a^(-2) * g
            True
        """
        return self.connection().ricci(name, latex_name)

    def ricci_scalar(self, name=None, latex_name=None):
        r"""
        Return the Ricci scalar associated with the metric.

        The Ricci scalar is the scalar field `r` defined from the Ricci tensor
        `Ric` and the metric tensor `g` by

        .. MATH::

            r = g^{ij} Ric_{ij}

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Ricci scalar;
          if none, it is set to "r(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Ricci scalar; if none, it is set to "\\mathrm{r}(g)", where "g"
          is the metric's name

        OUTPUT:

        - the Ricci scalar `r`, as an instance of
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`

        EXAMPLES:

        Ricci scalar of the standard metric on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of a meridian (domain of spherical coordinates)
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: a = var('a') # the sphere radius
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2] = a^2, a^2*sin(th)^2
            sage: g.display() # standard metric on the 2-sphere of radius a:
            g = a^2 dth⊗dth + a^2*sin(th)^2 dph⊗dph
            sage: g.ricci_scalar()
            Scalar field r(g) on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: g.ricci_scalar().display() # The Ricci scalar is constant:
            r(g): U → ℝ
               (th, ph) ↦ 2/a^2

        """
        if self._ricci_scalar is None:
            resu = self.inverse().contract(0, 1, self.ricci(), 0, 1)
            if name is None:
                name = "r(" + self._name + ")"
            if latex_name is None:
                latex_name = r"\mathrm{r}\left(" + self._latex_name + \
                              r"\right)"
            resu._name = name
            resu._latex_name = latex_name
            self._ricci_scalar = resu
        return self._ricci_scalar

    def weyl(self, name=None, latex_name=None):
        r"""
        Return the Weyl conformal tensor associated with the metric.

        The Weyl conformal tensor is the tensor field `C` of type (1,3)
        defined as the trace-free part of the Riemann curvature tensor `R`

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Weyl conformal
          tensor; if ``None``, it is set to "C(g)", where "g" is the metric's
          name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Weyl conformal tensor; if ``None``, it is set to "\\mathrm{C}(g)",
          where "g" is the metric's name

        OUTPUT:

        - the Weyl conformal tensor `C`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`

        EXAMPLES:

        Checking that the Weyl tensor identically vanishes on a 3-dimensional
        manifold, for instance the hyperbolic space `H^3`::

            sage: M = Manifold(3, 'H^3', start_index=1)
            sage: U = M.open_subset('U') # the complement of the half-plane (y=0, x>=0)
            sage: X.<rh,th,ph> = U.chart(r'rh:(0,+oo):\rho th:(0,pi):\theta  ph:(0,2*pi):\phi')
            sage: g = U.metric('g')
            sage: b = var('b')
            sage: g[1,1], g[2,2], g[3,3] = b^2, (b*sinh(rh))^2, (b*sinh(rh)*sin(th))^2
            sage: g.display()  # standard metric on H^3:
            g = b^2 drh⊗drh + b^2*sinh(rh)^2 dth⊗dth
             + b^2*sin(th)^2*sinh(rh)^2 dph⊗dph
            sage: C = g.weyl() ; C
            Tensor field C(g) of type (1,3) on the Open subset U of the
             3-dimensional differentiable manifold H^3
            sage: C == 0
            True
        """
        if self._weyl is None:
            n = self._ambient_domain.dimension()
            if n < 3:
                raise ValueError("the Weyl tensor is not defined for a " +
                                 "manifold of dimension n <= 2")
            delta = self._domain.tangent_identity_field(dest_map=self._vmodule._dest_map)
            riem = self.riemann()
            ric = self.ricci()
            rscal = self.ricci_scalar()
            # First index of the Ricci tensor raised with the metric
            ricup = ric.up(self, 0)
            aux = self*ricup + ric*delta - rscal/(n-1) * self*delta
            self._weyl = riem + 2/(n-2) * aux.antisymmetrize(2,3)
            if name is None:
                name = "C(" + self._name + ")"
            if latex_name is None:
                latex_name = r"\mathrm{C}\left(" + self._latex_name + r"\right)"
            self._weyl.set_name(name=name, latex_name=latex_name)
        return self._weyl

    def schouten(self, name=None, latex_name=None):
        r"""
        Return the Schouten tensor associated with the metric.

        The Schouten tensor is the tensor field `Sc` of type (0,2) defined
        from the Ricci curvature tensor `Ric` (see :meth:`ricci`) and the
        scalar curvature `r` (see :meth:`ricci_scalar`) and the metric `g` by

        .. MATH::

            Sc(u, v) = \frac{1}{n-2}\left(Ric(u, v) + \frac{r}{2(n-1)}g(u,v)
            \right)

        for any vector fields `u` and `v`.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Schouten tensor;
          if none, it is set to "Schouten(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Schouten tensor; if none, it is set to "\\mathrm{Schouten}(g)",
          where "g" is the metric's name

        OUTPUT:

        - the Schouten tensor `Sc`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField` of tensor
          type (0,2) and symmetric

        EXAMPLES:

        Schouten tensor of the left invariant metric of Heisenberg's
        Nil group::

            sage: M = Manifold(3, 'Nil', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.riemannian_metric('g')
            sage: g[1,1], g[2,2], g[2,3], g[3,3] = 1, 1+x^2, -x, 1
            sage: g.display()
            g = dx⊗dx + (x^2 + 1) dy⊗dy - x dy⊗dz - x dz⊗dy + dz⊗dz
            sage: g.schouten()
            Field of symmetric bilinear forms Schouten(g) on the 3-dimensional
             differentiable manifold Nil
            sage: g.schouten().display()
            Schouten(g) = -3/8 dx⊗dx + (5/8*x^2 - 3/8) dy⊗dy - 5/8*x dy⊗dz
             - 5/8*x dz⊗dy + 5/8 dz⊗dz
        """
        n = self._ambient_domain.dimension()
        if n < 3:
            raise ValueError("the Schouten tensor is only defined for a " +
                             "manifold of dimension >= 3")
        if self._schouten is None:
            s = (1/(n-2))*self.ricci() - (self.ricci_scalar()/(2*(n-1)*(n-2)))*self
            name = name or 'Schouten(' + self._name + ')'
            latex_name = latex_name or r'\mathrm{Schouten}(' + self._latex_name + ')'
            s.set_name(name=name, latex_name=latex_name)
            self._schouten = s
        return self._schouten

    def cotton(self, name=None, latex_name=None):
        r"""
        Return the Cotton conformal tensor associated with the metric.
        The tensor has type (0,3) and is defined in terms of the Schouten
        tensor `S` (see :meth:`schouten`):

        .. MATH::

            C_{ijk} = (n-2) \left(\nabla_k S_{ij}
            - \nabla_j S_{ik}\right)

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Cotton conformal
          tensor; if ``None``, it is set to "Cot(g)", where "g" is the metric's
          name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Cotton conformal tensor; if ``None``, it is set to "\\mathrm{Cot}(g)",
          where "g" is the metric's name

        OUTPUT:

        - the Cotton conformal tensor `Cot`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`

        EXAMPLES:

        Checking that the Cotton tensor identically vanishes on a conformally flat
        3-dimensional manifold, for instance the hyperbolic space `H^3`::

            sage: M = Manifold(3, 'H^3', start_index=1)
            sage: U = M.open_subset('U') # the complement of the half-plane (y=0, x>=0)
            sage: X.<rh,th,ph> = U.chart(r'rh:(0,+oo):\rho th:(0,pi):\theta  ph:(0,2*pi):\phi')
            sage: g = U.metric('g')
            sage: b = var('b')
            sage: g[1,1], g[2,2], g[3,3] = b^2, (b*sinh(rh))^2, (b*sinh(rh)*sin(th))^2
            sage: g.display()  # standard metric on H^3:
            g = b^2 drh⊗drh + b^2*sinh(rh)^2 dth⊗dth
             + b^2*sin(th)^2*sinh(rh)^2 dph⊗dph
            sage: Cot = g.cotton() ; Cot # long time
            Tensor field Cot(g) of type (0,3) on the Open subset U of the
             3-dimensional differentiable manifold H^3
            sage: Cot == 0 # long time
            True
        """
        n = self._ambient_domain.dimension()
        if n < 3:
            raise ValueError("the Cotton tensor is only defined for a " +
                             "manifold of dimension >= 3")
        if self._cotton is None:
            nabla = self.connection()
            s = self.schouten()
            cot = 2*(n-2)*nabla(s).antisymmetrize(1,2)
            name = name or 'Cot(' + self._name + ')'
            latex_name = latex_name or r'\mathrm{Cot}(' + self._latex_name + ')'
            cot.set_name(name=name, latex_name=latex_name)
            self._cotton = cot
        return self._cotton

    def cotton_york(self, name=None, latex_name=None):
        r"""
        Return the Cotton-York conformal tensor associated with the metric.
        The tensor has type (0,2) and is only defined for manifolds of
        dimension 3. It is defined in terms of the Cotton tensor `C`
        (see :meth:`cotton`) or the Schouten tensor `S` (see :meth:`schouten`):

        .. MATH::

            CY_{ij} = \frac{1}{2} \epsilon^{kl}_{\ \ \, i} C_{jlk}
                    = \epsilon^{kl}_{\ \ \, i} \nabla_k S_{lj}

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Cotton-York
          tensor; if ``None``, it is set to "CY(g)", where "g" is the metric's
          name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Cotton-York tensor; if ``None``, it is set to "\\mathrm{CY}(g)",
          where "g" is the metric's name

        OUTPUT:

        - the Cotton-York conformal tensor `CY`, as an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`

        EXAMPLES:

        Compute the determinant of the Cotton-York tensor for the Heisenberg
        group with the left invariant metric::

            sage: M = Manifold(3, 'Nil', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.riemannian_metric('g')
            sage: g[1,1], g[2,2], g[2,3], g[3,3] = 1, 1+x^2, -x, 1
            sage: g.display()
            g = dx⊗dx + (x^2 + 1) dy⊗dy - x dy⊗dz - x dz⊗dy + dz⊗dz
            sage: CY = g.cotton_york() ; CY # long time
            Tensor field CY(g) of type (0,2) on the 3-dimensional
             differentiable manifold Nil
            sage: CY.display()  # long time
            CY(g) = 1/2 dx⊗dx + (-x^2 + 1/2) dy⊗dy + x dy⊗dz + x dz⊗dy - dz⊗dz
            sage: det(CY[:]) # long time
            -1/4
        """
        n = self._ambient_domain.dimension()
        if n != 3:
            raise ValueError("the Cotton-York tensor is only defined for a " +
                             "manifold of dimension 3")
        if self._cotton_york is None:
            cot = self.cotton()
            eps = self.volume_form(2)
            cy = eps.contract(0, 1, cot, 2, 1)/2
            name = name or 'CY(' + self._name + ')'
            latex_name = latex_name or r'\mathrm{CY}(' + self._latex_name + ')'
            cy.set_name(name=name, latex_name=latex_name)
            self._cotton_york = cy
        return self._cotton_york

    def determinant(self, frame=None):
        r"""
        Determinant of the metric components in the specified frame.

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with
          respect to which the components `g_{ij}` of the metric are defined;
          if ``None``, the default frame of the metric's domain is used. If a
          chart is provided instead of a frame, the associated coordinate
          frame is used

        OUTPUT:

        - the determinant `\det (g_{ij})`, as an instance of
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`

        EXAMPLES:

        Metric determinant on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[1, 2], g[2, 2] = 1+x, x*y , 1-y
            sage: g[:]
            [ x + 1    x*y]
            [   x*y -y + 1]
            sage: s = g.determinant()  # determinant in M's default frame
            sage: s.expr()
            -x^2*y^2 - (x + 1)*y + x + 1

        A shortcut is ``det()``::

            sage: g.det() == g.determinant()
            True

        The notation ``det(g)`` can be used::

            sage: det(g) == g.determinant()
            True

        Determinant in a frame different from the default's one::

            sage: Y.<u,v> = M.chart()
            sage: ch_X_Y = X.transition_map(Y, [x+y, x-y])
            sage: ch_X_Y.inverse()
            Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y))
            sage: g.comp(Y.frame())[:, Y]
            [ 1/8*u^2 - 1/8*v^2 + 1/4*v + 1/2                            1/4*u]
            [                           1/4*u -1/8*u^2 + 1/8*v^2 + 1/4*v + 1/2]
            sage: g.determinant(Y.frame()).expr()
            -1/4*x^2*y^2 - 1/4*(x + 1)*y + 1/4*x + 1/4
            sage: g.determinant(Y.frame()).expr(Y)
            -1/64*u^4 - 1/64*v^4 + 1/32*(u^2 + 2)*v^2 - 1/16*u^2 + 1/4*v + 1/4

        A chart can be passed instead of a frame::

            sage: g.determinant(X) is g.determinant(X.frame())
            True
            sage: g.determinant(Y) is g.determinant(Y.frame())
            True

        The metric determinant depends on the frame::

            sage: g.determinant(X.frame()) == g.determinant(Y.frame())
            False

        Using SymPy as symbolic engine::

            sage: M.set_calculus_method('sympy')
            sage: g = M.metric('g')
            sage: g[1,1], g[1, 2], g[2, 2] = 1+x, x*y , 1-y
            sage: s = g.determinant()  # determinant in M's default frame
            sage: s.expr()
            -x**2*y**2 + x - y*(x + 1) + 1
        """
        from sage.matrix.constructor import matrix
        dom = self._domain
        if frame is None:
            frame = dom._def_frame
        if frame in dom._atlas:
            # frame is actually a chart and is changed to the associated
            # coordinate frame:
            frame = frame._frame
        if frame not in self._determinants:
            # a new computation is necessary
            resu = frame._domain.scalar_field()
            manif = self._ambient_domain
            gg = self.comp(frame)
            i1 = manif.start_index()
            for chart in gg[[i1, i1]]._express:
                # TODO: do the computation without the 'SR' enforcement
                gm = matrix( [[ gg[i, j, chart].expr(method='SR')
                            for j in manif.irange()] for i in manif.irange()] )
                detgm = chart.simplify(gm.det(), method='SR')
                resu.add_expr(detgm, chart=chart)
            self._determinants[frame] = resu
        return self._determinants[frame]

    det = determinant

    def sqrt_abs_det(self, frame=None):
        r"""
        Square root of the absolute value of the determinant of the metric
        components in the specified frame.

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with
          respect to which the components `g_{ij}` of ``self`` are defined;
          if ``None``, the domain's default frame is used. If a chart is
          provided, the associated coordinate frame is used

        OUTPUT:

        - `\sqrt{|\det (g_{ij})|}`, as an instance of
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`

        EXAMPLES:

        Standard metric in the Euclidean space `\RR^3` with spherical
        coordinates::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: U = M.open_subset('U') # the complement of the half-plane (y=0, x>=0)
            sage: c_spher.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2, (r*sin(th))^2
            sage: g.display()
            g = dr⊗dr + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph
            sage: g.sqrt_abs_det().expr()
            r^2*sin(th)

        Metric determinant on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[1, 2], g[2, 2] = 1+x, x*y , 1-y
            sage: g[:]
            [ x + 1    x*y]
            [   x*y -y + 1]
            sage: s = g.sqrt_abs_det() ; s
            Scalar field on the 2-dimensional differentiable manifold M
            sage: s.expr()
            sqrt(-x^2*y^2 - (x + 1)*y + x + 1)

        Determinant in a frame different from the default's one::

            sage: Y.<u,v> = M.chart()
            sage: ch_X_Y = X.transition_map(Y, [x+y, x-y])
            sage: ch_X_Y.inverse()
            Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y))
            sage: g[Y.frame(),:,Y]
            [ 1/8*u^2 - 1/8*v^2 + 1/4*v + 1/2                            1/4*u]
            [                           1/4*u -1/8*u^2 + 1/8*v^2 + 1/4*v + 1/2]
            sage: g.sqrt_abs_det(Y.frame()).expr()
            1/2*sqrt(-x^2*y^2 - (x + 1)*y + x + 1)
            sage: g.sqrt_abs_det(Y.frame()).expr(Y)
            1/8*sqrt(-u^4 - v^4 + 2*(u^2 + 2)*v^2 - 4*u^2 + 16*v + 16)

        A chart can be passed instead of a frame::

            sage: g.sqrt_abs_det(Y) is g.sqrt_abs_det(Y.frame())
            True

        The metric determinant depends on the frame::

            sage: g.sqrt_abs_det(X.frame()) == g.sqrt_abs_det(Y.frame())
            False

        Using SymPy as symbolic engine::

            sage: M.set_calculus_method('sympy')
            sage: g = M.metric('g')
            sage: g[1,1], g[1, 2], g[2, 2] = 1+x, x*y , 1-y
            sage: g.sqrt_abs_det().expr()
            sqrt(-x**2*y**2 - x*y + x - y + 1)
            sage: g.sqrt_abs_det(Y.frame()).expr()
            sqrt(-x**2*y**2 - x*y + x - y + 1)/2
            sage: g.sqrt_abs_det(Y.frame()).expr(Y)
            sqrt(-u**4 + 2*u**2*v**2 - 4*u**2 - v**4 + 4*v**2 + 16*v + 16)/8
        """
        dom = self._domain
        if frame is None:
            frame = dom._def_frame
        if frame in dom._atlas:
            # frame is actually a chart and is changed to the associated
            # coordinate frame:
            frame = frame._frame
        if frame not in self._sqrt_abs_dets:
            # a new computation is necessary
            detg = self.determinant(frame)
            resu = frame._domain.scalar_field()
            for chart, funct in detg._express.items():
                x = (self._indic_signat * funct).sqrt().expr()
                resu.add_expr(x, chart=chart)
            self._sqrt_abs_dets[frame] = resu
        return self._sqrt_abs_dets[frame]

    @overload
    def volume_form(self) -> DiffForm: ...
    @overload
    def volume_form(self, contra: int) -> TensorField: ...
    def volume_form(self, contra=0):
        r"""
        Volume form (Levi-Civita tensor) `\epsilon` associated with the metric.

        The volume form `\epsilon` is an `n`-form (`n` being the manifold's
        dimension) such that for any oriented vector basis `(e_i)` which is
        orthonormal with respect to the metric, the condition

        .. MATH::

            \epsilon(e_1,\ldots,e_n) = 1

        holds.

        Notice that a volume form requires an orientable manifold with
        a preferred orientation, see
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.orientation`
        for details.

        INPUT:

        - ``contra`` -- (default: 0) number of contravariant indices of the
          returned tensor

        OUTPUT:

        - if ``contra = 0`` (default value): the volume `n`-form `\epsilon`, as
          an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
        - if ``contra = k``, with `1\leq k \leq n`, the tensor field of type
          (k,n-k) formed from `\epsilon` by raising the first k indices with
          the metric (see method
          :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.up`);
          the output is then an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`, with
          the appropriate antisymmetries, or of the subclass
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
          if `k=n`

        EXAMPLES:

        Volume form on `\RR^3` with spherical coordinates, using the standard
        orientation, which is predefined::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: U = M.open_subset('U') # the complement of the half-plane (y=0, x>=0)
            sage: c_spher.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2, (r*sin(th))^2
            sage: g.display()
            g = dr⊗dr + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph
            sage: eps = g.volume_form() ; eps
            3-form eps_g on the Open subset U of the 3-dimensional
             differentiable manifold M
            sage: eps.display()
            eps_g = r^2*sin(th) dr∧dth∧dph
            sage: eps[[1,2,3]] == g.sqrt_abs_det()
            True
            sage: latex(eps)
            \epsilon_{g}

        The tensor field of components `\epsilon^i_{\ \, jk}` (``contra=1``)::

            sage: eps1 = g.volume_form(1) ; eps1
            Tensor field of type (1,2) on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: eps1.symmetries()
            no symmetry;  antisymmetry: (1, 2)
            sage: eps1[:]
            [[[0, 0, 0], [0, 0, r^2*sin(th)], [0, -r^2*sin(th), 0]],
             [[0, 0, -sin(th)], [0, 0, 0], [sin(th), 0, 0]],
             [[0, 1/sin(th), 0], [-1/sin(th), 0, 0], [0, 0, 0]]]

        The tensor field of components `\epsilon^{ij}_{\ \ k}` (``contra=2``)::

            sage: eps2 = g.volume_form(2) ; eps2
            Tensor field of type (2,1) on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: eps2.symmetries()
            no symmetry;  antisymmetry: (0, 1)
            sage: eps2[:]
            [[[0, 0, 0], [0, 0, sin(th)], [0, -1/sin(th), 0]],
             [[0, 0, -sin(th)], [0, 0, 0], [1/(r^2*sin(th)), 0, 0]],
             [[0, 1/sin(th), 0], [-1/(r^2*sin(th)), 0, 0], [0, 0, 0]]]

        The tensor field of components `\epsilon^{ijk}` (``contra=3``)::

            sage: eps3 = g.volume_form(3) ; eps3
            3-vector field on the Open subset U of the 3-dimensional
             differentiable manifold M
            sage: eps3.tensor_type()
            (3, 0)
            sage: eps3.symmetries()
            no symmetry;  antisymmetry: (0, 1, 2)
            sage: eps3[:]
            [[[0, 0, 0], [0, 0, 1/(r^2*sin(th))], [0, -1/(r^2*sin(th)), 0]],
             [[0, 0, -1/(r^2*sin(th))], [0, 0, 0], [1/(r^2*sin(th)), 0, 0]],
             [[0, 1/(r^2*sin(th)), 0], [-1/(r^2*sin(th)), 0, 0], [0, 0, 0]]]
            sage: eps3[1,2,3]
            1/(r^2*sin(th))
            sage: eps3[[1,2,3]] * g.sqrt_abs_det() == 1
            True

        If the manifold has no predefined orientation, an orientation must be
        set before invoking ``volume_form()``. For instance let consider the
        2-sphere described by the stereographic charts from the North and
        South pole::

            sage: M = Manifold(2, 'M', structure='Riemannian')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart()  # stereographic chart from the North pole
            sage: c_uv.<u,v> = V.chart()  # stereographic chart from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: eU = c_xy.frame(); eV = c_uv.frame()
            sage: g = M.metric()
            sage: g[eU,0,0], g[eU,1,1] = 4/(1+x^2+y^2)^2, 4/(1+x^2+y^2)^2
            sage: g.add_comp_by_continuation(eV, U.intersection(V), chart=c_uv)
            sage: eps = g.volume_form()
            Traceback (most recent call last):
            ...
            ValueError: 2-dimensional Riemannian manifold M must admit an
             orientation

        Let us define the orientation of ``M`` such that ``eU`` is
        right-handed; ``eV`` is then left-handed and in order to define an
        orientation on the whole of ``M``, we introduce a vector frame on ``V``
        by swapping ``eV``'s vectors::

            sage: f = V.vector_frame('f', (eV[1], eV[0]))
            sage: M.set_orientation([eU, f])

        We have then, factorizing the components for a nicer display::

            sage: eps = g.volume_form()
            sage: eps.apply_map(factor, frame=eU, keep_other_components=True)
            sage: eps.apply_map(factor, frame=eV, keep_other_components=True)
            sage: eps.display(eU)
            eps_g = 4/(x^2 + y^2 + 1)^2 dx∧dy
            sage: eps.display(eV)
            eps_g = -4/(u^2 + v^2 + 1)^2 du∧dv

        Note the minus sign in the above expression, reflecting the fact that
        ``eV`` is left-handed with respect to the chosen orientation.
        """
        dom = self._domain
        orient = dom.orientation()
        manif = self._ambient_domain
        ndim = manif.dimension()
        if not orient:
            raise ValueError('{} must admit an orientation'.format(dom))
        if contra > ndim:
            raise ValueError('The number of contravariant indices is greater '
                             'than the manifold dimension')
        if self._vol_forms == []:
            # a new computation is necessary
            # The result is constructed on the vector field module,
            # so that dest_map is taken automatically into account:
            eps = self._vmodule.alternating_form(ndim, name='eps_'+self._name,
                                latex_name=r'\epsilon_{'+self._latex_name+r'}')
            si = manif.start_index()
            ind = tuple(range(si, si+ndim))
            for frame in orient:
                if frame.destination_map() is frame.domain().identity_map():
                    eps.add_comp(frame)[[ind]] = self.sqrt_abs_det(frame)
            self._vol_forms.append(eps)  # Levi-Civita tensor constructed
        if contra >= len(self._vol_forms):
            # Tensors related to the Levi-Civita one by index rising:
            for k in range(len(self._vol_forms), contra+1):
                epskm1 = self._vol_forms[k-1]
                epsk = epskm1.up(self, k-1)
                if k > 1:
                    # restoring the antisymmetry after the up operation:
                    epsk = epsk.antisymmetrize(*range(k))
                self._vol_forms.append(epsk)
        return self._vol_forms[contra]

    def hodge_star(self, pform: DiffForm) -> DiffForm:
        r"""
        Compute the Hodge dual of a differential form with respect to the
        metric.

        If the differential form is a `p`-form `A`, its *Hodge dual* with
        respect to the metric `g` is the
        `(n-p)`-form `*A` defined by

        .. MATH::

            *A_{i_1\ldots i_{n-p}} = \frac{1}{p!} A_{k_1\ldots k_p}
                \epsilon^{k_1\ldots k_p}_{\qquad\ i_1\ldots i_{n-p}}

        where `n` is the manifold's dimension, `\epsilon` is the volume
        `n`-form associated with `g` (see :meth:`volume_form`) and the indices
        `k_1,\ldots, k_p` are raised with `g`.

        Notice that the hodge star dual requires an orientable manifold
        with a preferred orientation, see
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.orientation`
        for details.

        INPUT:

        - ``pform`` -- a `p`-form `A`; must be an instance of
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`
          for `p=0` and of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm` or
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`
          for `p\geq 1`.

        OUTPUT:

        - the `(n-p)`-form `*A`

        EXAMPLES:

        Hodge dual of a 1-form in the Euclidean space `R^3`::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: var('Ax Ay Az')
            (Ax, Ay, Az)
            sage: a = M.one_form(Ax, Ay, Az, name='A')
            sage: sa = g.hodge_star(a) ; sa
            2-form *A on the 3-dimensional differentiable manifold M
            sage: sa.display()
            *A = Az dx∧dy - Ay dx∧dz + Ax dy∧dz
            sage: ssa = g.hodge_star(sa) ; ssa
            1-form **A on the 3-dimensional differentiable manifold M
            sage: ssa.display()
            **A = Ax dx + Ay dy + Az dz
            sage: ssa == a  # must hold for a Riemannian metric in dimension 3
            True

        Hodge dual of a 0-form (scalar field) in `R^3`::

            sage: f = M.scalar_field(function('F')(x,y,z), name='f')
            sage: sf = g.hodge_star(f) ; sf
            3-form *f on the 3-dimensional differentiable manifold M
            sage: sf.display()
            *f = F(x, y, z) dx∧dy∧dz
            sage: ssf = g.hodge_star(sf) ; ssf
            Scalar field **f on the 3-dimensional differentiable manifold M
            sage: ssf.display()
            **f: M → ℝ
               (x, y, z) ↦ F(x, y, z)
            sage: ssf == f # must hold for a Riemannian metric
            True

        Hodge dual of a 0-form in Minkowski spacetime::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = M.lorentzian_metric('g')
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: g.display()  # Minkowski metric
            g = -dt⊗dt + dx⊗dx + dy⊗dy + dz⊗dz
            sage: var('f0')
            f0
            sage: f = M.scalar_field(f0, name='f')
            sage: sf = g.hodge_star(f) ; sf
            4-form *f on the 4-dimensional differentiable manifold M
            sage: sf.display()
            *f = f0 dt∧dx∧dy∧dz
            sage: ssf = g.hodge_star(sf) ; ssf
            Scalar field **f on the 4-dimensional differentiable manifold M
            sage: ssf.display()
            **f: M → ℝ
               (t, x, y, z) ↦ -f0
            sage: ssf == -f  # must hold for a Lorentzian metric
            True

        Hodge dual of a 1-form in Minkowski spacetime::

            sage: var('At Ax Ay Az')
            (At, Ax, Ay, Az)
            sage: a = M.one_form(At, Ax, Ay, Az, name='A')
            sage: a.display()
            A = At dt + Ax dx + Ay dy + Az dz
            sage: sa = g.hodge_star(a) ; sa
            3-form *A on the 4-dimensional differentiable manifold M
            sage: sa.display()
            *A = -Az dt∧dx∧dy + Ay dt∧dx∧dz - Ax dt∧dy∧dz - At dx∧dy∧dz
            sage: ssa = g.hodge_star(sa) ; ssa
            1-form **A on the 4-dimensional differentiable manifold M
            sage: ssa.display()
            **A = At dt + Ax dx + Ay dy + Az dz
            sage: ssa == a  # must hold for a Lorentzian metric in dimension 4
            True

        Hodge dual of a 2-form in Minkowski spacetime::

            sage: F = M.diff_form(2, name='F')
            sage: var('Ex Ey Ez Bx By Bz')
            (Ex, Ey, Ez, Bx, By, Bz)
            sage: F[0,1], F[0,2], F[0,3] = -Ex, -Ey, -Ez
            sage: F[1,2], F[1,3], F[2,3] = Bz, -By, Bx
            sage: F[:]
            [  0 -Ex -Ey -Ez]
            [ Ex   0  Bz -By]
            [ Ey -Bz   0  Bx]
            [ Ez  By -Bx   0]
            sage: sF = g.hodge_star(F) ; sF
            2-form *F on the 4-dimensional differentiable manifold M
            sage: sF[:]
            [  0  Bx  By  Bz]
            [-Bx   0  Ez -Ey]
            [-By -Ez   0  Ex]
            [-Bz  Ey -Ex   0]
            sage: ssF = g.hodge_star(sF) ; ssF
            2-form **F on the 4-dimensional differentiable manifold M
            sage: ssF[:]
            [  0  Ex  Ey  Ez]
            [-Ex   0 -Bz  By]
            [-Ey  Bz   0 -Bx]
            [-Ez -By  Bx   0]
            sage: ssF.display()
            **F = Ex dt∧dx + Ey dt∧dy + Ez dt∧dz - Bz dx∧dy + By dx∧dz
             - Bx dy∧dz
            sage: F.display()
            F = -Ex dt∧dx - Ey dt∧dy - Ez dt∧dz + Bz dx∧dy - By dx∧dz
             + Bx dy∧dz
            sage: ssF == -F  # must hold for a Lorentzian metric in dimension 4
            True

        Test of the standard identity

        .. MATH::

            *(A\wedge B) = \epsilon(A^\sharp, B^\sharp, ., .)

        where `A` and `B` are any 1-forms and `A^\sharp` and `B^\sharp` the
        vectors associated to them by the metric `g` (index raising)::

            sage: var('Bt Bx By Bz')
            (Bt, Bx, By, Bz)
            sage: b = M.one_form(Bt, Bx, By, Bz, name='B')
            sage: b.display()
            B = Bt dt + Bx dx + By dy + Bz dz
            sage: epsilon = g.volume_form()
            sage: g.hodge_star(a.wedge(b)) == epsilon.contract(0,a.up(g)).contract(0,b.up(g))
            True
        """
        return pform.hodge_dual(self)


#******************************************************************************

class PseudoRiemannianMetricParal(PseudoRiemannianMetric, TensorFieldParal):
    r"""
    Pseudo-Riemannian metric with values on a parallelizable manifold.

    An instance of this class is a field of nondegenerate symmetric bilinear
    forms (metric field) along a differentiable manifold `U` with values in a
    parallelizable manifold `M` over `\RR`, via a differentiable mapping
    `\Phi: U \rightarrow M`. The standard case of a metric field *on* a
    manifold corresponds to `U=M` and `\Phi = \mathrm{Id}_M`. Other common
    cases are `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is
    then an open interval of `\RR`).

    A *metric* `g` is a field on `U`, such that at each
    point `p\in U`, `g(p)` is a bilinear map of the type:

    .. MATH::

        g(p):\ T_q M\times T_q M  \longrightarrow \RR

    where `T_q M` stands for the tangent space to manifold `M` at the point
    `q=\Phi(p)`, such that `g(p)` is symmetric:
    `\forall (u,v)\in  T_q M\times T_q M, \ g(p)(v,u) = g(p)(u,v)`
    and nondegenerate:
    `(\forall v\in T_q M,\ \ g(p)(u,v) = 0) \Longrightarrow u=0`.

    .. NOTE::

        If `M` is not parallelizable, the class :class:`PseudoRiemannianMetric`
        should be used instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `\Phi(U)\subset M`
    - ``name`` -- name given to the metric
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      single integer: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number
      of positive terms (resp. number of negative terms) in any diagonal
      writing of the metric components; if ``signature`` is ``None``, `S` is
      set to the dimension of manifold `M` (Riemannian signature)
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the metric;
      if ``None``, it is formed from ``name``

    EXAMPLES:

    Metric on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M', start_index=1)
        sage: c_xy.<x,y> = M.chart()
        sage: g = M.metric('g') ; g
        Riemannian metric g on the 2-dimensional differentiable manifold M
        sage: latex(g)
        g

    A metric is a special kind of tensor field and therefore inheritates all the
    properties from class
    :class:`~sage.manifolds.differentiable.tensorfield.TensorField`::

        sage: g.parent()
        Free module T^(0,2)(M) of type-(0,2) tensors fields on the
         2-dimensional differentiable manifold M
        sage: g.tensor_type()
        (0, 2)
        sage: g.symmetries()  # g is symmetric:
        symmetry: (0, 1);  no antisymmetry

    Setting the metric components in the manifold's default frame::

        sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-x
        sage: g[:]
        [ x + 1    x*y]
        [   x*y -x + 1]
        sage: g.display()
        g = (x + 1) dx⊗dx + x*y dx⊗dy + x*y dy⊗dx + (-x + 1) dy⊗dy

    Metric components in a frame different from the manifold's default one::

        sage: c_uv.<u,v> = M.chart()  # new chart on M
        sage: xy_to_uv = c_xy.transition_map(c_uv, [x+y, x-y]) ; xy_to_uv
        Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
        sage: uv_to_xy = xy_to_uv.inverse() ; uv_to_xy
        Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y))
        sage: M.atlas()
        [Chart (M, (x, y)), Chart (M, (u, v))]
        sage: M.frames()
        [Coordinate frame (M, (∂/∂x,∂/∂y)), Coordinate frame (M, (∂/∂u,∂/∂v))]
        sage: g[c_uv.frame(),:]  # metric components in frame c_uv.frame() expressed in M's default chart (x,y)
        [ 1/2*x*y + 1/2          1/2*x]
        [         1/2*x -1/2*x*y + 1/2]
        sage: g.display(c_uv.frame())
        g = (1/2*x*y + 1/2) du⊗du + 1/2*x du⊗dv + 1/2*x dv⊗du
         + (-1/2*x*y + 1/2) dv⊗dv
        sage: g[c_uv.frame(),:,c_uv]   # metric components in frame c_uv.frame() expressed in chart (u,v)
        [ 1/8*u^2 - 1/8*v^2 + 1/2            1/4*u + 1/4*v]
        [           1/4*u + 1/4*v -1/8*u^2 + 1/8*v^2 + 1/2]
        sage: g.display(c_uv.frame(), c_uv)
        g = (1/8*u^2 - 1/8*v^2 + 1/2) du⊗du + (1/4*u + 1/4*v) du⊗dv
         + (1/4*u + 1/4*v) dv⊗du + (-1/8*u^2 + 1/8*v^2 + 1/2) dv⊗dv

    As a shortcut of the above command, on can pass just the chart ``c_uv``
    to ``display``, the vector frame being then assumed to be the coordinate
    frame associated with the chart::

        sage: g.display(c_uv)
        g = (1/8*u^2 - 1/8*v^2 + 1/2) du⊗du + (1/4*u + 1/4*v) du⊗dv
         + (1/4*u + 1/4*v) dv⊗du + (-1/8*u^2 + 1/8*v^2 + 1/2) dv⊗dv

    The inverse metric is obtained via :meth:`inverse`::

        sage: ig = g.inverse() ; ig
        Tensor field inv_g of type (2,0) on the 2-dimensional differentiable
         manifold M
        sage: ig[:]
        [ (x - 1)/(x^2*y^2 + x^2 - 1)      x*y/(x^2*y^2 + x^2 - 1)]
        [     x*y/(x^2*y^2 + x^2 - 1) -(x + 1)/(x^2*y^2 + x^2 - 1)]
        sage: ig.display()
        inv_g = (x - 1)/(x^2*y^2 + x^2 - 1) ∂/∂x⊗∂/∂x
         + x*y/(x^2*y^2 + x^2 - 1) ∂/∂x⊗∂/∂y + x*y/(x^2*y^2 + x^2 - 1) ∂/∂y⊗∂/∂x
         - (x + 1)/(x^2*y^2 + x^2 - 1) ∂/∂y⊗∂/∂y
    """
    def __init__(self, vector_field_module, name, signature=None,
                 latex_name=None):
        r"""
        Construct a metric on a parallelizable manifold.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: from sage.manifolds.differentiable.metric import \
            ....:                                   PseudoRiemannianMetricParal
            sage: g = PseudoRiemannianMetricParal(XM, 'g', signature=0); g
            Lorentzian metric g on the 2-dimensional differentiable manifold M
            sage: g[0,0], g[1,1] = -(1+x^2), 1+y^2
            sage: TestSuite(g).run(skip='_test_category')

        .. TODO::

            - add a specific parent to the metrics, to fit with the category
              framework
        """
        TensorFieldParal.__init__(self, vector_field_module, (0,2),
                                  name=name, latex_name=latex_name, sym=(0,1))
        # signature:
        ndim = self._ambient_domain.dimension()
        if signature is None:
            signature = ndim
        else:
            if not isinstance(signature, (int, Integer)):
                raise TypeError("the metric signature must be an integer")
            if (signature < - ndim) or (signature > ndim):
                raise ValueError("metric signature out of range")
            if (signature+ndim) % 2 == 1:
                if ndim % 2 == 0:
                    raise ValueError("the metric signature must be even")
                else:
                    raise ValueError("the metric signature must be odd")
        self._signature = signature
        # the pair (n_+, n_-):
        self._signature_pm = ((ndim+signature)//2, (ndim-signature)//2)
        self._indic_signat = 1 - 2*(self._signature_pm[1] % 2)  # (-1)^n_-
        # Initialization of derived quantities:
        PseudoRiemannianMetricParal._init_derived(self)

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: g = M.metric('g')
            sage: g._init_derived()
        """
        # Initialization of quantities pertaining to the mother classes:
        TensorFieldParal._init_derived(self)
        PseudoRiemannianMetric._init_derived(self)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- boolean (default: ``True``); determines whether the
          restrictions of ``self`` to subdomains are deleted

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: g = M.metric('g')
            sage: g._del_derived(del_restrictions=False)
            sage: g._del_derived()
        """
        # The derived quantities from the mother classes are deleted:
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)
        PseudoRiemannianMetric._del_derived(self)

    def _del_inverse(self):
        r"""
        Delete the inverse metric.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: g = M.metric('g')
            sage: g._del_inverse()
        """
        self._inverse._components.clear()
        self._inverse._del_derived()

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of the metric to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of ``self._domain`` (must be an
          instance of
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `V` is a subdomain of
          ``self._codomain``
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          If ``None``, the restriction of ``self._vmodule._dest_map`` to `U` is
          used.

        OUTPUT:

        - instance of :class:`PseudoRiemannianMetricParal` representing the
          restriction.

        EXAMPLES:

        Restriction of a Lorentzian metric on `\RR^2` to the upper half plane::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: g = M.lorentzian_metric('g')
            sage: g[0,0], g[1,1] = -1, 1
            sage: U = M.open_subset('U', coord_def={X: y>0})
            sage: gU = g.restrict(U); gU
            Lorentzian metric g on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: gU.signature()
            0
            sage: gU.display()
            g = -dx⊗dx + dy⊗dy
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            # Construct the restriction at the tensor field level:
            resu = TensorFieldParal.restrict(self, subdomain, dest_map=dest_map)
            # the type is correctly handled by TensorFieldParal.restrict, i.e.
            # resu is of type self.__class__, but the signature is not handled
            # by TensorFieldParal.restrict; we have to set it here:
            resu._signature = self._signature
            resu._signature_pm = self._signature_pm
            resu._indic_signat = self._indic_signat
            # Restrictions of derived quantities:
            resu._inverse = self.inverse().restrict(subdomain)
            if self._connection is not None:
                resu._connection = self._connection.restrict(subdomain)
            if self._ricci_scalar is not None:
                resu._ricci_scalar = self._ricci_scalar.restrict(subdomain)
            if self._weyl is not None:
                resu._weyl = self._weyl.restrict(subdomain)
            if self._vol_forms != []:
                for eps in self._vol_forms:
                    resu._vol_forms.append(eps.restrict(subdomain))
            # NB: no initialization of resu._determinants nor
            # resu._sqrt_abs_dets
            # The restriction is ready:
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def set(self, symbiform):
        r"""
        Define the metric from a field of symmetric bilinear forms.

        INPUT:

        - ``symbiform`` -- instance of
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
          representing a field of symmetric bilinear forms

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: s = M.sym_bilin_form_field(name='s')
            sage: s[0,0], s[0,1], s[1,1] = 1+x^2, x*y, 1+y^2
            sage: g = M.metric('g')
            sage: g.set(s)
            sage: g.display()
            g = (x^2 + 1) dx⊗dx + x*y dx⊗dy + x*y dy⊗dx + (y^2 + 1) dy⊗dy
        """
        if not isinstance(symbiform, TensorFieldParal):
            raise TypeError("the argument must be a tensor field with " +
                            "values on a parallelizable domain")
        if symbiform._tensor_type != (0,2):
            raise TypeError("the argument must be of tensor type (0,2)")
        if symbiform._sym != ((0,1),):
            raise TypeError("the argument must be symmetric")
        if symbiform._vmodule is not self._vmodule:
            raise TypeError("the symmetric bilinear form and the metric are " +
                            "not defined on the same vector field module")
        self._del_derived()
        self._components.clear()
        for frame in symbiform._components:
            self._components[frame] = symbiform._components[frame].copy()
        for dom, symbiform_rst in symbiform._restrictions.items():
            rst = self.restrict(dom)
            rst.set(symbiform_rst)

    def inverse(self, expansion_symbol=None, order=1):
        r"""
        Return the inverse metric.

        INPUT:

        - ``expansion_symbol`` -- (default: ``None``) symbolic variable; if
          specified, the inverse will be expanded in power series with respect
          to this variable (around its zero value)
        - ``order`` -- integer (default: 1); the order of the expansion
          if ``expansion_symbol`` is not ``None``; the *order* is defined as
          the degree of the polynomial representing the truncated power series
          in ``expansion_symbol``; currently only first order inverse is
          supported

        If ``expansion_symbol`` is set, then the zeroth order metric must be
        invertible. Moreover, subsequent calls to this method will return
        a cached value, even when called with the default value (to enable
        computation of derived quantities). To reset, use :meth:`_del_derived`.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
          with ``tensor_type`` = (2,0) representing the inverse metric

        EXAMPLES:

        Inverse metric on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-x
            sage: g[:]  # components in the manifold's default frame
            [ x + 1    x*y]
            [   x*y -x + 1]
            sage: ig = g.inverse() ; ig
            Tensor field inv_g of type (2,0) on the 2-dimensional
              differentiable manifold M
            sage: ig[:]
            [ (x - 1)/(x^2*y^2 + x^2 - 1)      x*y/(x^2*y^2 + x^2 - 1)]
            [     x*y/(x^2*y^2 + x^2 - 1) -(x + 1)/(x^2*y^2 + x^2 - 1)]

        If the metric is modified, the inverse metric is automatically updated::

            sage: g[1,2] = 0 ; g[:]
            [ x + 1      0]
            [     0 -x + 1]
            sage: g.inverse()[:]
            [ 1/(x + 1)          0]
            [         0 -1/(x - 1)]

        Using SymPy as symbolic engine::

            sage: M.set_calculus_method('sympy')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-x
            sage: g[:]  # components in the manifold's default frame
            [x + 1   x*y]
            [  x*y 1 - x]
            sage: g.inverse()[:]
            [ (x - 1)/(x**2*y**2 + x**2 - 1)      x*y/(x**2*y**2 + x**2 - 1)]
            [     x*y/(x**2*y**2 + x**2 - 1) -(x + 1)/(x**2*y**2 + x**2 - 1)]

        Demonstration of the series expansion capabilities::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: C.<t,x,y,z> = M.chart()
            sage: e = var('e')
            sage: g = M.metric()
            sage: h = M.tensor_field(0, 2, sym=(0,1))
            sage: g[0, 0], g[1, 1], g[2, 2], g[3, 3] = -1, 1, 1, 1
            sage: h[0, 1], h[1, 2], h[2, 3] = 1, 1, 1
            sage: g.set(g + e*h)

        If ``e`` is a small parameter, ``g`` is a tridiagonal approximation of
        the Minkowski metric::

            sage: g[:]
            [-1  e  0  0]
            [ e  1  e  0]
            [ 0  e  1  e]
            [ 0  0  e  1]

        The inverse, truncated to first order in ``e``, is::

            sage: g.inverse(expansion_symbol=e)[:]
            [-1  e  0  0]
            [ e  1 -e  0]
            [ 0 -e  1 -e]
            [ 0  0 -e  1]

        If ``inverse()`` is called subsequently, the result will be the same.
        This allows for all computations to be made to first order::

            sage: g.inverse()[:]
            [-1  e  0  0]
            [ e  1 -e  0]
            [ 0 -e  1 -e]
            [ 0  0 -e  1]
        """
        if expansion_symbol is not None:
            if (self._inverse is not None and bool(self._inverse._components)
                and self._inverse._components.values()[0][0,0]._expansion_symbol
                    == expansion_symbol
                and self._inverse._components.values()[0][0,0]._order == order):
                return self._inverse

            if order != 1:
                raise NotImplementedError("only first order inverse is implemented")
            decompo = self.series_expansion(expansion_symbol, order)
            g0 = decompo[0]
            g1 = decompo[1]

            g0m = self._new_instance()   # needed because only metrics have
            g0m.set_comp()[:] = g0[:]    # an "inverse" method.

            contraction = g1.contract(0, g0m.inverse(), 0)
            contraction = contraction.contract(1, g0m.inverse(), 1)
            self._inverse = g0m.inverse() - expansion_symbol * contraction
            self._inverse.set_calc_order(expansion_symbol, order)
            return self._inverse

        from sage.matrix.constructor import matrix
        from sage.tensor.modules.comp import CompFullySym
        # Is the inverse metric up to date ?
        for frame in self._components:
            if frame not in self._inverse._components:
                # the computation is necessary
                fmodule = self._fmodule
                si = fmodule._sindex
                nsi = fmodule._rank + si
                dom = self._domain
                cinv = CompFullySym(fmodule._ring, frame, 2, start_index=si,
                                    output_formatter=fmodule._output_formatter)
                cinv_scal = {}  # dict. of scalars representing the components
                                # of the inverse (keys: comp. indices)
                for i in range(si, nsi):
                    for j in range(i, nsi):   # symmetry taken into account
                        cinv_scal[(i,j)] = dom.scalar_field()
                for chart in dom.top_charts():
                    # TODO: do the computation without the 'SR' enforcement
                    try:
                        gmat = matrix(
                                  [[self.comp(frame)[i, j, chart].expr(method='SR')
                                  for j in range(si, nsi)] for i in range(si, nsi)])
                        gmat_inv = gmat.inverse()
                    except (KeyError, ValueError):
                        continue
                    for i in range(si, nsi):
                        for j in range(i, nsi):
                            val = chart.simplify(gmat_inv[i-si,j-si], method='SR')
                            cinv_scal[(i,j)].add_expr(val, chart=chart)
                for i in range(si, nsi):
                    for j in range(i, nsi):
                        cinv[i,j] = cinv_scal[(i,j)]
                self._inverse._components[frame] = cinv
        return self._inverse

    def ricci_scalar(self, name=None, latex_name=None):
        r"""
        Return the metric's Ricci scalar.

        The Ricci scalar is the scalar field `r` defined from the Ricci tensor
        `Ric` and the metric tensor `g` by

        .. MATH::

            r = g^{ij} Ric_{ij}

        INPUT:

        - ``name`` -- (default: ``None``) name given to the Ricci scalar;
          if none, it is set to "r(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          Ricci scalar; if none, it is set to "\\mathrm{r}(g)", where "g"
          is the metric's name

        OUTPUT:

        - the Ricci scalar `r`, as an instance of
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`

        EXAMPLES:

        Ricci scalar of the standard metric on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of a meridian (domain of spherical coordinates)
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: a = var('a') # the sphere radius
            sage: g = U.metric('g')
            sage: g[1,1], g[2,2] = a^2, a^2*sin(th)^2
            sage: g.display() # standard metric on the 2-sphere of radius a:
            g = a^2 dth⊗dth + a^2*sin(th)^2 dph⊗dph
            sage: g.ricci_scalar()
            Scalar field r(g) on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: g.ricci_scalar().display() # The Ricci scalar is constant:
            r(g): U → ℝ
               (th, ph) ↦ 2/a^2

        """
        if self._ricci_scalar is None:
            manif = self._ambient_domain
            ric = self.ricci()
            ig = self.inverse()
            frame = ig.common_basis(ric)
            cric = ric._components[frame]
            cig = ig._components[frame]
            rsum1 = 0
            for i in manif.irange():
                rsum1 += cig[[i,i]] * cric[[i,i]]
            rsum2 = 0
            for i in manif.irange():
                for j in manif.irange(start=i+1):
                    rsum2 += cig[[i,j]] * cric[[i,j]]
            self._ricci_scalar = rsum1 + 2*rsum2
            if name is None:
                self._ricci_scalar._name = "r(" + self._name + ")"
            else:
                self._ricci_scalar._name = name
            if latex_name is None:
                self._ricci_scalar._latex_name = r"\mathrm{r}\left(" + \
                                                 self._latex_name + r"\right)"
            else:
                self._ricci_scalar._latex_name = latex_name
        return self._ricci_scalar


#****************************************************************************************************


class DegenerateMetric(TensorField):
    r"""
    Degenerate (or null or lightlike) metric with values on an open subset of a
    differentiable manifold.

    An instance of this class is a field of degenerate symmetric bilinear
    forms (metric field) along a differentiable manifold `U` with
    values on a differentiable manifold `M` over `\RR`, via a differentiable
    mapping `\Phi: U \rightarrow M`.
    The standard case of a degenerate metric field *on* a manifold corresponds to `U=M`
    and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval
    of `\RR`).

    A *degenerate metric* `g` is a field on `U`, such that at each point `p\in U`, `g(p)`
    is a bilinear map of the type:

    .. MATH::

        g(p):\ T_q M\times T_q M  \longrightarrow \RR

    where `T_q M` stands for the tangent space to the
    manifold `M` at the point `q=\Phi(p)`, such that `g(p)` is symmetric:
    `\forall (u,v)\in  T_q M\times T_q M, \ g(p)(v,u) = g(p)(u,v)`
    and degenerate:
    `\exists v\in T_q M;\ \ g(p)(u,v) = 0\ \ \forall u\in T_qM`.

    .. NOTE::

        If `M` is parallelizable, the class :class:`DegenerateMetricParal`
        should be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `\Phi(U)\subset M`
    - ``name`` -- name given to the metric
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      tuple: `S = (n_+, n_-, n_0)`, where `n_+` (resp. `n_-`, resp. `n_0`) is the
      number of positive terms (resp. negative terms, resp. zero tems) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to `(ndim-1, 0, 1)`, being `ndim` the manifold's dimension
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the metric;
      if ``None``, it is formed from ``name``

    EXAMPLES:

    Lightlike cone::

        sage: M = Manifold(3, 'M'); X.<x,y,z> = M.chart()
        sage: g = M.metric('g', signature=(2,0,1)); g
        degenerate metric g on the 3-dimensional differentiable manifold M
        sage: det(g)
        Scalar field zero on the 3-dimensional differentiable manifold M
        sage: g.parent()
        Free module T^(0,2)(M) of type-(0,2) tensors fields on the
        3-dimensional differentiable manifold M
        sage: g[0,0], g[0,1], g[0,2] = (y^2 + z^2)/(x^2 + y^2 + z^2), \
        ....: - x*y/(x^2 + y^2 + z^2), - x*z/(x^2 + y^2 + z^2)
        sage: g[1,1], g[1,2], g[2,2] = (x^2 + z^2)/(x^2 + y^2 + z^2), \
        ....: - y*z/(x^2 + y^2 + z^2), (x^2 + y^2)/(x^2 + y^2 + z^2)
        sage: g.disp()
        g = (y^2 + z^2)/(x^2 + y^2 + z^2) dx⊗dx - x*y/(x^2 + y^2 + z^2) dx⊗dy
        - x*z/(x^2 + y^2 + z^2) dx⊗dz - x*y/(x^2 + y^2 + z^2) dy⊗dx
        + (x^2 + z^2)/(x^2 + y^2 + z^2) dy⊗dy - y*z/(x^2 + y^2 + z^2) dy⊗dz
        - x*z/(x^2 + y^2 + z^2) dz⊗dx - y*z/(x^2 + y^2 + z^2) dz⊗dy
        + (x^2 + y^2)/(x^2 + y^2 + z^2) dz⊗dz

    The position vector is a lightlike vector field::

        sage: v = M.vector_field()
        sage: v[0], v[1], v[2] = x , y, z
        sage: g(v, v).disp()
        M → ℝ
        (x, y, z) ↦ 0
    """

    def __init__(self, vector_field_module, name, signature=None,
                 latex_name=None):
        r"""
        Construct a metric.

        TESTS::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: var('m'); assume(m>0)
            m
            sage: Int = M.open_subset('Int')
            sage: X.<t,r,th,ph>=Int.chart(r"t r:(0,2*m) th:(0,pi):\theta ph:(0,2*pi):\phi")
            sage: XM = M.vector_field_module(); e = X.frame()
            sage: from sage.manifolds.differentiable.metric import \
            ....:                                        DegenerateMetric
            sage: g = DegenerateMetric(XM, 'g', signature=(2,1,1)); g
            degenerate metric g on the 4-dimensional Lorentzian manifold M
            sage: g[e, 0,0], g[e, 0,1], g[e, 1,1], g[e, 2,2], \
            ....: g[e, 3,3] = -1+2*m/r, 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2
            sage: g.disp(e)
            g = (2*m/r - 1) dt⊗dt + 2*m/r dt⊗dr + 2*m/r dr⊗dt + (2*m/r + 1) dr⊗dr
            + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph
        """
        TensorField.__init__(self, vector_field_module, (0,2),
                             name=name, latex_name=latex_name, sym=(0,1))
        # signature:
        ndim = self._ambient_domain.dimension()
        if signature is None:
            signature = (ndim-1,0,1)
        else:
            try:
                for elt in signature:
                    if (elt < 0) or (not isinstance(elt, (int, Integer))):
                        raise ValueError("{} must be a positive integer".format(elt))
                    if elt > ndim:
                        raise ValueError("{} must be less than {}".format(elt,ndim))
                    sign = signature[0]+signature[1]+signature[2]
                    if sign != ndim:
                        raise ValueError("{} is different from {}".format(sign, ndim))
            except TypeError:
                raise TypeError("signature must be an iterable")
        self._signature = (signature[0],signature[1],signature[2])
        # the tuple (n_+, n_-, n_0):
        self._signature_pm = self._signature

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: g = M.metric('g', signature=(1,1,1))
            sage: g._repr_()
            'degenerate metric g on the 3-dimensional differentiable manifold M'
        """
        return self._final_repr("degenerate metric "+self._name + " ")

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` with the same
        signature.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: g = M.metric('g', signature=(1,1,1))
            sage: g1 = g._new_instance(); g1
            degenerate metric unnamed metric on the 3-dimensional differentiable manifold M
            sage: type(g1) == type(g)
            True
            sage: g1.parent() is g.parent()
            True
            sage: g1.signature() == g.signature()
            True
        """
        return type(self)(self._vmodule, 'unnamed metric',
                          signature=self._signature,
                          latex_name=r'\text{unnamed metric}')

    def signature(self):
        r"""
        Signature of the metric.

        OUTPUT:

        - signature of a degenerate metric is defined as the tuple
          `(n_+,n_-,n_0)`, where `n_+` (resp. `n_-`, resp. `n_0`) is the number of
          positive terms (resp. negative terms, resp. zero terms) eigenvalues

        EXAMPLES:

        Signatures on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M')
            sage: g = M.metric('g', signature=(1,1,1))
            sage: g.signature()
            (1, 1, 1)
            sage: M = Manifold(3, 'M', structure='degenerate_metric')
            sage: g = M.metric()
            sage: g.signature()
            (0, 2, 1)
        """
        return self._signature

    def set(self, symbiform):
        r"""
        Define the metric from a field of symmetric bilinear forms.

        INPUT:

        - ``symbiform`` -- instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          representing a field of symmetric bilinear forms

        EXAMPLES:

        Metric defined from a field of symmetric bilinear forms on a
        non-parallelizable 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: h = M.sym_bilin_form_field(name='h')
            sage: h[eU,0,0], h[eU,0,1], h[eU,1,1] = 1+x, x*y, 1-y
            sage: h.add_comp_by_continuation(eV, W, c_uv)
            sage: h.display(eU)
            h = (x + 1) dx⊗dx + x*y dx⊗dy + x*y dy⊗dx + (-y + 1) dy⊗dy
            sage: h.display(eV)
            h = (1/8*u^2 - 1/8*v^2 + 1/4*v + 1/2) du⊗du + 1/4*u du⊗dv
             + 1/4*u dv⊗du + (-1/8*u^2 + 1/8*v^2 + 1/4*v + 1/2) dv⊗dv
            sage: g = M.metric('g')
            sage: g.set(h)
            sage: g.display(eU)
            g = (x + 1) dx⊗dx + x*y dx⊗dy + x*y dy⊗dx + (-y + 1) dy⊗dy
            sage: g.display(eV)
            g = (1/8*u^2 - 1/8*v^2 + 1/4*v + 1/2) du⊗du + 1/4*u du⊗dv
             + 1/4*u dv⊗du + (-1/8*u^2 + 1/8*v^2 + 1/4*v + 1/2) dv⊗dv
        """
        if not isinstance(symbiform, TensorField):
            raise TypeError("the argument must be a tensor field")
        if symbiform._tensor_type != (0,2):
            raise TypeError("the argument must be of tensor type (0,2)")
        if symbiform._sym != ((0,1),):
            raise TypeError("the argument must be symmetric")
        if not symbiform._domain.is_subset(self._domain):
            raise TypeError("the symmetric bilinear form is not defined " +
                            "on the metric domain")
        self._restrictions.clear()
        if isinstance(symbiform, TensorFieldParal):
            rst = self.restrict(symbiform._domain)
            rst.set(symbiform)
        else:
            for dom, symbiform_rst in symbiform._restrictions.items():
                rst = self.restrict(dom)
                rst.set(symbiform_rst)

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of the metric to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of the metric's domain (must be an
          instance of :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `V` is a subdomain of
          ``self._codomain``
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          If ``None``, the restriction of ``self._vmodule._dest_map`` to `U` is
          used.

        OUTPUT:

        - instance of :class:`DegenerateMetric` representing the
          restriction.

        EXAMPLES::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g', signature=(3,1,1))
            sage: U = M.open_subset('U')
            sage: g.restrict(U)
            degenerate metric g on the Open subset U of the 5-dimensional
            differentiable manifold M
            sage: g.restrict(U).signature()
            (3, 1, 1)

        See the top documentation of :class:`DegenerateMetric` for more
        examples.
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            # Construct the restriction at the tensor field level:
            resu = TensorField.restrict(self, subdomain, dest_map=dest_map)
            # the type is correctly handled by TensorField.restrict, i.e.
            # resu is of type self.__class__, but the signature is not handled
            # by TensorField.restrict; we have to set it here:
            resu._signature = self._signature
            resu._signature_pm = self._signature_pm
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def determinant(self):
        r"""
        Determinant of a degenerate metric is always '0'.

        EXAMPLES::

            sage: S = Manifold(2, 'S')
            sage: g = S.metric('g', signature=([0,1,1]))
            sage: g.determinant()
            Scalar field zero on the 2-dimensional differentiable manifold S
        """
        return self._domain.zero_scalar_field()

    det = determinant

#****************************************************************************************


class DegenerateMetricParal(DegenerateMetric, TensorFieldParal):
    r"""
    Degenerate (or null or lightlike) metric with values on an open subset of a
    differentiable manifold.

    An instance of this class is a field of degenerate symmetric bilinear
    forms (metric field) along a differentiable manifold `U` with
    values on a differentiable manifold `M` over `\RR`, via a differentiable
    mapping `\Phi: U \rightarrow M`.
    The standard case of a degenerate metric field *on* a manifold corresponds
    to `U=M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval
    of `\RR`).

    A *degenerate metric* `g` is a field on `U`, such that at each point
    `p\in U`, `g(p)` is a bilinear map of the type:

    .. MATH::

        g(p):\ T_q M\times T_q M  \longrightarrow \RR

    where `T_q M` stands for the tangent space to the
    manifold `M` at the point `q=\Phi(p)`, such that `g(p)` is symmetric:
    `\forall (u,v)\in  T_q M\times T_q M, \ g(p)(v,u) = g(p)(u,v)`
    and degenerate:
    `\exists v\in T_q M;\ \ g(p)(u,v) = 0\ \ \forall u\in T_qM`.

    .. NOTE::

        If `M` is not parallelizable, the class :class:`DegenerateMetric`
        should be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `\Phi(U)\subset M`
    - ``name`` -- name given to the metric
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      tuple: `S = (n_+, n_-, n_0)`, where `n_+` (resp. `n_-`, resp. `n_0`) is the
      number of positive terms (resp. negative terms, resp. zero tems) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to `(ndim-1, 0, 1)`, being `ndim` the manifold's dimension
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the metric;
      if ``None``, it is formed from ``name``

    EXAMPLES:

    Lightlike cone::

        sage: M = Manifold(3, 'M'); X.<x,y,z> = M.chart()
        sage: g = M.metric('g', signature=(2,0,1)); g
        degenerate metric g on the 3-dimensional differentiable manifold M
        sage: det(g)
        Scalar field zero on the 3-dimensional differentiable manifold M
        sage: g.parent()
        Free module T^(0,2)(M) of type-(0,2) tensors fields on the
        3-dimensional differentiable manifold M
        sage: g[0,0], g[0,1], g[0,2] = (y^2 + z^2)/(x^2 + y^2 + z^2), \
        ....: - x*y/(x^2 + y^2 + z^2), - x*z/(x^2 + y^2 + z^2)
        sage: g[1,1], g[1,2], g[2,2] = (x^2 + z^2)/(x^2 + y^2 + z^2), \
        ....: - y*z/(x^2 + y^2 + z^2), (x^2 + y^2)/(x^2 + y^2 + z^2)
        sage: g.disp()
        g = (y^2 + z^2)/(x^2 + y^2 + z^2) dx⊗dx - x*y/(x^2 + y^2 + z^2) dx⊗dy
        - x*z/(x^2 + y^2 + z^2) dx⊗dz - x*y/(x^2 + y^2 + z^2) dy⊗dx
        + (x^2 + z^2)/(x^2 + y^2 + z^2) dy⊗dy - y*z/(x^2 + y^2 + z^2) dy⊗dz
        - x*z/(x^2 + y^2 + z^2) dz⊗dx - y*z/(x^2 + y^2 + z^2) dz⊗dy
        + (x^2 + y^2)/(x^2 + y^2 + z^2) dz⊗dz

    The position vector is a lightlike vector field::

        sage: v = M.vector_field()
        sage: v[0], v[1], v[2] = x , y, z
        sage: g(v, v).disp()
        M → ℝ
        (x, y, z) ↦ 0
    """

    def __init__(self, vector_field_module, name, signature=None,
                 latex_name=None):
        r"""
        Construct a metric.

        TESTS::

            sage: M = Manifold(3, 'M'); X.<x,y,z> = M.chart()
            sage: XM = M.vector_field_module()
            sage: from sage.manifolds.differentiable.metric import \
            ....:                                        DegenerateMetricParal
            sage: g = DegenerateMetricParal(XM, 'g', signature=(2,0,1)); g
            degenerate metric g on the 3-dimensional differentiable manifold M
            sage: g[0,0], g[0,1], g[0,2] = (y^2 + z^2)/(x^2 + y^2 + z^2), \
            ....: - x*y/(x^2 + y^2 + z^2), - x*z/(x^2 + y^2 + z^2)
            sage: g[1,1], g[1,2], g[2,2] = (x^2 + z^2)/(x^2 + y^2 + z^2), \
            ....: - y*z/(x^2 + y^2 + z^2), (x^2 + y^2)/(x^2 + y^2 + z^2)
            sage: g.disp()
            g = (y^2 + z^2)/(x^2 + y^2 + z^2) dx⊗dx - x*y/(x^2 + y^2 + z^2) dx⊗dy
            - x*z/(x^2 + y^2 + z^2) dx⊗dz - x*y/(x^2 + y^2 + z^2) dy⊗dx
            + (x^2 + z^2)/(x^2 + y^2 + z^2) dy⊗dy - y*z/(x^2 + y^2 + z^2) dy⊗dz
            - x*z/(x^2 + y^2 + z^2) dz⊗dx - y*z/(x^2 + y^2 + z^2) dz⊗dy
            + (x^2 + y^2)/(x^2 + y^2 + z^2) dz⊗dz
        """
        TensorFieldParal.__init__(self, vector_field_module, (0,2),
                             name=name, latex_name=latex_name, sym=(0,1))
        # signature:
        ndim = self._ambient_domain.dimension()
        if signature is None:
            signature = (ndim-1,0,1)
        else:
            try:
                for elt in signature:
                    if (elt < 0) or (not isinstance(elt, (int, Integer))):
                        raise ValueError("{} must be a positive integer".format(elt))
                    sign = signature[0]+signature[1]+signature[2]
                    if sign != ndim:
                        raise ValueError("{} is different from {}".format(sign, ndim))
            except TypeError:
                raise TypeError("signature must be an iterable")
        self._signature = (signature[0],signature[1],signature[2])
        # the tuple (n_+, n_-, n_0):
        self._signature_pm = self._signature

    def set(self, symbiform):
        r"""
        Define the metric from a field of symmetric bilinear forms.

        INPUT:

        - ``symbiform`` -- instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          representing a field of symmetric bilinear forms

        EXAMPLES:

        Metric defined from a field of symmetric bilinear forms on a
        parallelizable 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1);
            sage: X.<x,y,z> = M.chart()
            sage: dx, dy = X.coframe()[1], X.coframe()[2]
            sage: b = dx*dx + dy*dy
            sage: g = M.metric('g', signature=(1,1,1)); g
            degenerate metric g on the 3-dimensional differentiable manifold M
            sage: g.set(b)
            sage: g.display()
            g = dx⊗dx + dy⊗dy
        """
        if not isinstance(symbiform, TensorFieldParal):
            raise TypeError("the argument must be a tensor field with " +
                            "values on a parallelizable domain")
        if symbiform._tensor_type != (0,2):
            raise TypeError("the argument must be of tensor type (0,2)")
        if symbiform._sym != ((0,1),):
            raise TypeError("the argument must be symmetric")
        if symbiform._vmodule is not self._vmodule:
            raise TypeError("the symmetric bilinear form and the metric are " +
                            "not defined on the same vector field module")
        self._components.clear()
        for frame in symbiform._components:
            self._components[frame] = symbiform._components[frame].copy()
        for dom, symbiform_rst in symbiform._restrictions.items():
            rst = self.restrict(dom)
            rst.set(symbiform_rst)

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of the metric to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of the metric's domain (must be an
          instance of :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `V` is a subdomain of
          ``self._codomain``
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          If ``None``, the restriction of ``self._vmodule._dest_map`` to `U` is
          used.

        OUTPUT:

        - instance of :class:`DegenerateMetric` representing the
          restriction.

        EXAMPLES::

            sage: M = Manifold(5, 'M')
            sage: g = M.metric('g', signature=(3,1,1))
            sage: U = M.open_subset('U')
            sage: g.restrict(U)
            degenerate metric g on the Open subset U of the 5-dimensional differentiable manifold M
            sage: g.restrict(U).signature()
            (3, 1, 1)

        See the top documentation of :class:`DegenerateMetric` for more
        examples.
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            # Construct the restriction at the tensor field level:
            resu = TensorFieldParal.restrict(self, subdomain, dest_map=dest_map)
            # the type is correctly handled by TensorField.restrict, i.e.
            # resu is of type self.__class__, but the signature is not handled
            # by TensorField.restrict; we have to set it here:
            resu._signature = self._signature
            resu._signature_pm = self._signature_pm
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

#****************************************************************************************
