# sage_setup: distribution = sagemath-symbolics
r"""
Tensor Fields

The class :class:`TensorField` implements tensor fields on differentiable
manifolds. The derived class
:class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
is devoted to tensor fields with values on parallelizable manifolds.

Various derived classes of :class:`TensorField` are devoted to specific tensor
fields:

* :class:`~sage.manifolds.differentiable.vectorfield.VectorField` for vector
  fields (rank-1 contravariant tensor fields)

* :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
  for fields of tangent-space automorphisms

* :class:`~sage.manifolds.differentiable.diff_form.DiffForm` for differential
  forms (fully antisymmetric covariant tensor fields)

* :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
  for multivector fields (fully antisymmetric contravariant tensor fields)

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Travis Scrimshaw (2016): review tweaks
- Eric Gourgoulhon (2018): operators divergence, Laplacian and d'Alembertian;
  method :meth:`TensorField.along`
- Florentin Jaffredo (2018) : series expansion with respect to a given
  parameter
- Michael Jung (2019): improve treatment of the zero element; add method
  :meth:`TensorField.copy_from`
- Eric Gourgoulhon (2020): add method :meth:`TensorField.apply_map`

REFERENCES:

- [KN1963]_
- [Lee2013]_
- [ONe1983]_
"""

# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import annotations

from typing import TYPE_CHECKING, Optional, TypeVar, Union

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.element import ModuleElementWithMutability
from sage.tensor.modules.comp import CompWithSym
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.tensor_with_indices import TensorWithIndices

if TYPE_CHECKING:
    from sage.manifolds.differentiable.diff_map import DiffMap
    from sage.manifolds.differentiable.manifold import DifferentiableManifold
    from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
    from sage.manifolds.differentiable.poisson_tensor import PoissonTensorField
    from sage.manifolds.differentiable.symplectic_form import SymplecticForm
    from sage.manifolds.differentiable.vectorfield_module import VectorFieldModule
    from sage.manifolds.point import ManifoldPoint
    from sage.tensor.modules.comp import Components


TensorType = tuple[int, int]
T = TypeVar("T", bound='TensorField')


class TensorField(ModuleElementWithMutability):
    r"""
    Tensor field along a differentiable manifold.

    An instance of this class is a tensor field along a differentiable
    manifold `U` with values on a differentiable manifold `M`, via a
    differentiable map `\Phi: U \rightarrow M`. More precisely, given two
    nonnegative integers `k` and `l` and a differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *tensor field of type* `(k,l)` *along* `U` *with values on* `M` is
    a differentiable map

    .. MATH::

        t:\ U  \longrightarrow T^{(k,l)}M

    (where `T^{(k,l)}M` is the tensor bundle of type `(k,l)` over `M`) such
    that

    .. MATH::

        \forall p \in U,\ t(p) \in T^{(k,l)}(T_q M)

    i.e. `t(p)` is a tensor of type `(k,l)` on the tangent space `T_q M` at
    the point `q = \Phi(p)`, that is to say a multilinear map

    .. MATH::

        t(p):\ \underbrace{T_q^*M\times\cdots\times T_q^*M}_{k\ \; \text{times}}
        \times \underbrace{T_q M\times\cdots\times T_q M}_{l\ \; \text{times}}
        \longrightarrow K,

    where `T_q^* M` is the dual vector space to `T_q M` and `K` is the
    topological field over which the manifold `M` is defined. The integer `k+l`
    is called the *tensor rank*.

    The standard case of a tensor
    field *on* a differentiable manifold corresponds to `U=M` and
    `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval
    of `\RR`).

    If `M` is parallelizable, the class
    :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
    should be used instead.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` associated with the map `\Phi: U \rightarrow M` (cf.
      :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`)
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank
      and `l` the covariant rank
    - ``name`` -- (default: ``None``) name given to the tensor field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the tensor
      field; if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the tensor arguments: each symmetry is described by a tuple containing
      the positions of the involved arguments, with the convention
      ``position = 0`` for the first argument; for instance:

      * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
      * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
        arguments and a symmetry between the 2nd, 4th and 5th arguments.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the arguments, with the same convention as for ``sym``
    - ``parent`` -- (default: ``None``) some specific parent (e.g. exterior
      power for differential forms); if ``None``,
      ``vector_field_module.tensor_module(k,l)`` is used

    EXAMPLES:

    Tensor field of type (0,2) on the sphere `S^2`::

        sage: M = Manifold(2, 'S^2') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
        ....:                 restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: W = U.intersection(V)
        sage: t = M.tensor_field(0,2, name='t') ; t
        Tensor field t of type (0,2) on the 2-dimensional differentiable
         manifold S^2
        sage: t.parent()
        Module T^(0,2)(S^2) of type-(0,2) tensors fields on the 2-dimensional
         differentiable manifold S^2
        sage: t.parent().category()
        Category of tensor products of modules over Algebra of differentiable scalar fields
         on the 2-dimensional differentiable manifold S^2

    The parent of `t` is not a free module, for the sphere `S^2` is not
    parallelizable::

        sage: isinstance(t.parent(), FiniteRankFreeModule)
        False

    To fully define `t`, we have to specify its components in some vector
    frames defined on subsets of `S^2`; let us start by the open subset `U`::

        sage: eU = c_xy.frame()
        sage: t[eU,:] = [[1,0], [-2,3]]
        sage: t.display(eU)
        t = dx⊗dx - 2 dy⊗dx + 3 dy⊗dy

    To set the components of `t` on `V` consistently, we copy the expressions
    of the components in the common subset `W`::

        sage: eV = c_uv.frame()
        sage: eVW = eV.restrict(W)
        sage: c_uvW = c_uv.restrict(W)
        sage: t[eV,0,0] = t[eVW,0,0,c_uvW].expr()  # long time
        sage: t[eV,0,1] = t[eVW,0,1,c_uvW].expr()  # long time
        sage: t[eV,1,0] = t[eVW,1,0,c_uvW].expr()  # long time
        sage: t[eV,1,1] = t[eVW,1,1,c_uvW].expr()  # long time

    Actually, the above operation can be performed in a single line by means
    of the method
    :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.add_comp_by_continuation`::

        sage: t.add_comp_by_continuation(eV, W, chart=c_uv)  # long time

    At this stage, `t` is fully defined, having components in frames eU and eV
    and the union of the domains of eU and eV being the whole manifold::

        sage: t.display(eV)  # long time
        t = (u^4 - 4*u^3*v + 10*u^2*v^2 + 4*u*v^3 + v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du⊗du
         - 4*(u^3*v + 2*u^2*v^2 - u*v^3)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du⊗dv
         + 2*(u^4 - 2*u^3*v - 2*u^2*v^2 + 2*u*v^3 + v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv⊗du
         + (3*u^4 + 4*u^3*v - 2*u^2*v^2 - 4*u*v^3 + 3*v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv⊗dv

    Let us consider two vector fields, `a` and `b`, on `S^2`::

        sage: a = M.vector_field({eU: [-y, x]}, name='a')
        sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: a.display(eV)
        a = -v ∂/∂u + u ∂/∂v
        sage: b = M.vector_field({eU: [y, -1]}, name='b')
        sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: b.display(eV)
        b = ((2*u + 1)*v^3 + (2*u^3 - u^2)*v)/(u^2 + v^2) ∂/∂u
         - (u^4 - v^4 + 2*u*v^2)/(u^2 + v^2) ∂/∂v

    As a tensor field of type `(0,2)`, `t` acts on the pair `(a,b)`,
    resulting in a scalar field::

        sage: f = t(a,b); f
        Scalar field t(a,b) on the 2-dimensional differentiable manifold S^2
        sage: f.display()  # long time
        t(a,b): S^2 → ℝ
        on U: (x, y) ↦ -2*x*y - y^2 - 3*x
        on V: (u, v) ↦ -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)

    The vectors can be defined only on subsets of `S^2`, the domain of the
    result is then the common subset::

        sage: # long time
        sage: s = t(a.restrict(U), b) ; s
        Scalar field t(a,b) on the Open subset U of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()
        t(a,b): U → ℝ
           (x, y) ↦ -2*x*y - y^2 - 3*x
        on W: (u, v) ↦ -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)
        sage: s = t(a.restrict(U), b.restrict(W)) ; s
        Scalar field t(a,b) on the Open subset W of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()
        t(a,b): W → ℝ
           (x, y) ↦ -2*x*y - y^2 - 3*x
           (u, v) ↦ -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)

    The tensor itself can be defined only on some open subset of `S^2`,
    yielding a result whose domain is this subset::

        sage: s = t.restrict(V)(a,b); s  # long time
        Scalar field t(a,b) on the Open subset V of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()  # long time
        t(a,b): V → ℝ
           (u, v) ↦ -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)
        on W: (x, y) ↦ -2*x*y - y^2 - 3*x

    Tests regarding the multiplication by a scalar field::

        sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2),
        ....:                     c_uv: (u^2 + v^2)/(u^2 + v^2 + 1)}, name='f')
        sage: t.parent().base_ring() is f.parent()
        True
        sage: s = f*t; s  # long time
        Tensor field f*t of type (0,2) on the 2-dimensional differentiable
         manifold S^2
        sage: s[[0,0]] == f*t[[0,0]]  # long time
        True
        sage: s.restrict(U) == f.restrict(U) * t.restrict(U)  # long time
        True
        sage: s = f*t.restrict(U); s
        Tensor field f*t of type (0,2) on the Open subset U of the 2-dimensional
         differentiable manifold S^2
        sage: s.restrict(U) == f.restrict(U) * t.restrict(U)
        True

    .. RUBRIC:: Same examples with SymPy as the symbolic engine

    From now on, we ask that all symbolic calculus on manifold `M` are
    performed by SymPy::

        sage: M.set_calculus_method('sympy')

    We define the tensor `t` as above::

        sage: t = M.tensor_field(0, 2, {eU:  [[1,0], [-2,3]]}, name='t')
        sage: t.display(eU)
        t = dx⊗dx - 2 dy⊗dx + 3 dy⊗dy
        sage: t.add_comp_by_continuation(eV, W, chart=c_uv)  # long time
        sage: t.display(eV)  # long time
        t = (u**4 - 4*u**3*v + 10*u**2*v**2 + 4*u*v**3 + v**4)/(u**8 +
         4*u**6*v**2 + 6*u**4*v**4 + 4*u**2*v**6 + v**8) du⊗du +
         4*u*v*(-u**2 - 2*u*v + v**2)/(u**8 + 4*u**6*v**2 + 6*u**4*v**4
         + 4*u**2*v**6 + v**8) du⊗dv + 2*(u**4 - 2*u**3*v - 2*u**2*v**2
         + 2*u*v**3 + v**4)/(u**8 + 4*u**6*v**2 + 6*u**4*v**4 +
         4*u**2*v**6 + v**8) dv⊗du + (3*u**4 + 4*u**3*v - 2*u**2*v**2 -
         4*u*v**3 + 3*v**4)/(u**8 + 4*u**6*v**2 + 6*u**4*v**4 +
         4*u**2*v**6 + v**8) dv⊗dv

    The default coordinate representations of tensor components are now
    SymPy objects::

        sage: t[eV,1,1,c_uv].expr() # long time
        (3*u**4 + 4*u**3*v - 2*u**2*v**2 - 4*u*v**3 + 3*v**4)/(u**8 +
         4*u**6*v**2 + 6*u**4*v**4 + 4*u**2*v**6 + v**8)
        sage: type(t[eV,1,1,c_uv].expr()) # long time
        <class 'sympy.core.mul.Mul'>

    Let us consider two vector fields, `a` and `b`, on `S^2`::

        sage: a = M.vector_field({eU: [-y, x]}, name='a')
        sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: a.display(eV)
        a = -v ∂/∂u + u ∂/∂v
        sage: b = M.vector_field({eU: [y, -1]}, name='b')
        sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: b.display(eV)
        b = v*(2*u**3 - u**2 + 2*u*v**2 + v**2)/(u**2 + v**2) ∂/∂u
            + (-u**4 - 2*u*v**2 + v**4)/(u**2 + v**2) ∂/∂v

    As a tensor field of type `(0,2)`, `t` acts on the pair `(a,b)`,
    resulting in a scalar field::

        sage: f = t(a,b)
        sage: f.display()  # long time
        t(a,b): S^2 → ℝ
        on U: (x, y) ↦ -2*x*y - 3*x - y**2
        on V: (u, v) ↦ (-3*u**3 - 3*u*v**2 - 2*u*v - v**2)/(u**4 + 2*u**2*v**2 + v**4)

    The vectors can be defined only on subsets of `S^2`, the domain of the
    result is then the common subset::

        sage: s = t(a.restrict(U), b)
        sage: s.display()  # long time
        t(a,b): U → ℝ
           (x, y) ↦ -2*x*y - 3*x - y**2
        on W: (u, v) ↦ (-3*u**3 - 3*u*v**2 - 2*u*v - v**2)/(u**4 + 2*u**2*v**2 + v**4)
        sage: s = t(a.restrict(U), b.restrict(W))  # long time
        sage: s.display()  # long time
        t(a,b): W → ℝ
           (x, y) ↦ -2*x*y - 3*x - y**2
           (u, v) ↦ (-3*u**3 - 3*u*v**2 - 2*u*v - v**2)/(u**4 + 2*u**2*v**2 + v**4)

    The tensor itself can be defined only on some open subset of `S^2`,
    yielding a result whose domain is this subset::

        sage: s = t.restrict(V)(a,b)  # long time
        sage: s.display()  # long time
        t(a,b): V → ℝ
           (u, v) ↦ (-3*u**3 - 3*u*v**2 - 2*u*v - v**2)/(u**4 + 2*u**2*v**2 + v**4)
        on W: (x, y) ↦ -2*x*y - 3*x - y**2

    Tests regarding the multiplication by a scalar field::

        sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2),
        ....:                     c_uv: (u^2 + v^2)/(u^2 + v^2 + 1)}, name='f')
        sage: s = f*t # long time
        sage: s[[0,0]] == f*t[[0,0]]  # long time
        True
        sage: s.restrict(U) == f.restrict(U) * t.restrict(U)  # long time
        True
        sage: s = f*t.restrict(U)
        sage: s.restrict(U) == f.restrict(U) * t.restrict(U)
        True

    Notice that the zero tensor field is immutable, and therefore its
    components cannot be changed::

        sage: zer = M.tensor_field_module((1, 1)).zero()
        sage: zer.is_immutable()
        True
        sage: zer.set_comp()
        Traceback (most recent call last):
        ...
        ValueError: the components of an immutable element cannot be
         changed

    Other tensor fields can be declared immutable, too::

        sage: t.is_immutable()
        False
        sage: t.set_immutable()
        sage: t.is_immutable()
        True
        sage: t.set_comp()
        Traceback (most recent call last):
        ...
        ValueError: the components of an immutable element cannot be
         changed
        sage: t.set_name('b')
        Traceback (most recent call last):
        ...
        ValueError: the name of an immutable element cannot be changed
    """

    _name: Optional[str]
    _latex_name: Optional[str]
    _vmodule: VectorFieldModule
    _domain: DifferentiableManifold
    _ambient_domain: DifferentiableManifold

    def __init__(
        self,
        vector_field_module: VectorFieldModule,
        tensor_type: TensorType,
        name: Optional[str] = None,
        latex_name: Optional[str] = None,
        sym=None,
        antisym=None,
        parent=None,
    ):
        r"""
        Construct a tensor field.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``TensorField``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:          intersection_name='W', restrictions1= x>0,
            ....:          restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: XM = M.vector_field_module()
            sage: T02 = M.tensor_field_module((0,2))
            sage: t = T02.element_class(XM, (0,2), name='t'); t
            Tensor field t of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: t[e_xy,:] = [[1+x^2, x*y], [0, 1+y^2]]
            sage: t.add_comp_by_continuation(e_uv, W, c_uv)
            sage: t.display(e_xy)
            t = (x^2 + 1) dx⊗dx + x*y dx⊗dy + (y^2 + 1) dy⊗dy
            sage: t.display(e_uv)
            t = (3/16*u^2 + 1/16*v^2 + 1/2) du⊗du
             + (-1/16*u^2 + 1/4*u*v + 1/16*v^2) du⊗dv
             + (1/16*u^2 + 1/4*u*v - 1/16*v^2) dv⊗du
             + (1/16*u^2 + 3/16*v^2 + 1/2) dv⊗dv
            sage: TestSuite(t).run(skip='_test_pickling')

        Construction with ``DifferentiableManifold.tensor_field``::

            sage: t1 = M.tensor_field(0, 2, name='t'); t1
            Tensor field t of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: type(t1) == type(t)
            True
        """
        if parent is None:
            parent = vector_field_module.tensor_module(*tensor_type)
        ModuleElementWithMutability.__init__(self, parent)
        self._vmodule = vector_field_module
        self._tensor_type = tuple(tensor_type)
        self._tensor_rank = self._tensor_type[0] + self._tensor_type[1]
        self._is_zero = False
        # a priori, may be changed below or via  method __bool__()

        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain

        self._extensions_graph = {self._domain: self}
                    # dict. of known extensions of self on bigger domains,
                    # including self, with domains as keys. Its elements can be
                    # seen as incoming edges on a graph.
        self._restrictions_graph = {self._domain: self}
                    # dict. of known restrictions of self on smaller domains,
                    # including self, with domains as keys. Its elements can be
                    # seen as outgoing edges on a graph.

        self._restrictions = {} # dict. of restrictions of self on subdomains
                                # of self._domain, with the subdomains as keys
        # Treatment of symmetry declarations:
        self._sym, self._antisym = CompWithSym._canonicalize_sym_antisym(
            self._tensor_rank, sym, antisym)
        # Initialization of derived quantities:
        self._init_derived()

    # ###### Required methods for ModuleElement (beside arithmetic) #######

    def __bool__(self):
        r"""
        Return ``True`` if ``self`` is nonzero and ``False`` otherwise.

        This method is called by :meth:`is_zero`.

        EXAMPLES:

        Tensor field defined by parts on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1, 2, name='t')
            sage: tu = U.tensor_field(1, 2, name='t')
            sage: tv = V.tensor_field(1, 2, name='t')
            sage: tu[0,0,0] = 0
            sage: tv[0,0,0] = 0
            sage: t.set_restriction(tv)
            sage: t.set_restriction(tu)
            sage: bool(t)
            False
            sage: t.is_zero()  # indirect doctest
            True
            sage: tv[0,0,0] = 1
            sage: t.set_restriction(tv)
            sage: bool(t)
            True
            sage: t.is_zero()  # indirect doctest
            False
        """
        if self._is_zero:
            return False
        if any(bool(rst) for rst in self._restrictions.values()):
            self._is_zero = False
            return True
        self._is_zero = True
        return False

    # #### End of required methods for ModuleElement (beside arithmetic) #####

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t
            Tensor field t of type (1,3) on the 2-dimensional differentiable manifold M
        """
        # Special cases
        if self._tensor_type == (0,2) and self._sym == ((0,1),):
            description = "Field of symmetric bilinear forms "
            if self._name is not None:
                description += self._name + " "
        else:
        # Generic case
            description = "Tensor field "
            if self._name is not None:
                description += self._name + " "
            description += "of type ({},{}) ".format(
                                    self._tensor_type[0], self._tensor_type[1])
        return self._final_repr(description)

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._latex_()
            't'
            sage: t = M.tensor_field(1, 3, name='t', latex_name=r'\tau')
            sage: latex(t)
            \tau
        """
        if self._latex_name is None:
            return r'\text{' + str(self) + r'}'
        else:
            return self._latex_name

    def set_name(self, name: Optional[str] = None, latex_name: Optional[str] = None):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``name`` -- string (default: ``None``); name given to the tensor
          field
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to denote
          the tensor field; if ``None`` while ``name`` is provided, the LaTeX
          symbol is set to ``name``

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3); t
            Tensor field of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: t.set_name(name='t')
            sage: t
            Tensor field t of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: latex(t)
            t
            sage: t.set_name(latex_name=r'\tau')
            sage: latex(t)
            \tau
            sage: t.set_name(name='a')
            sage: t
            Tensor field a of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: latex(t)
            a
        """
        if self.is_immutable():
            raise ValueError("the name of an immutable element "
                             "cannot be changed")
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name
        for rst in self._restrictions.values():
            rst.set_name(name=name, latex_name=latex_name)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same
        vector field module, with the same tensor type and same symmetries

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t1 = t._new_instance(); t1
            Tensor field of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: type(t1) == type(t)
            True
            sage: t1.parent() is t.parent()
            True
        """
        return type(self)(self._vmodule, self._tensor_type, sym=self._sym,
                          antisym=self._antisym, parent=self.parent())

    def _final_repr(self, description: str) -> str:
        r"""
        Part of string representation common to all derived classes of
        :class:`TensorField`.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._final_repr('Tensor field t ')
            'Tensor field t on the 2-dimensional differentiable manifold M'
        """
        if self._domain == self._ambient_domain:
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} ".format(self._domain) + \
                           "with values on the {}".format(self._ambient_domain)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._init_derived()
        """
        self._lie_derivatives = {} # dict. of Lie derivatives of self (keys: id(vector))

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._del_derived()
        """
        # First deletes any reference to self in the vectors' dictionaries:
        for val in self._lie_derivatives.values():
            del val[0]._lie_der_along_self[id(self)]
        # Then clears the dictionary of Lie derivatives
        self._lie_derivatives.clear()

    def _del_restrictions(self):
        r"""
        Delete the restrictions defined on ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: t = M.tensor_field(1,2)
            sage: U = M.open_subset('U', coord_def={c_xy: x<0})
            sage: h = t.restrict(U)
            sage: t._restrictions
            {Open subset U of the 2-dimensional differentiable manifold M:
             Tensor field of type (1,2) on the Open subset U of the
             2-dimensional differentiable manifold M}
            sage: t._del_restrictions()
            sage: t._restrictions
            {}
        """
        self._restrictions.clear()
        self._extensions_graph = {self._domain: self}
        self._restrictions_graph = {self._domain: self}

    def _init_components(self, *comp, **kwargs):
        r"""
        Initialize the tensor field components in some given vector frames.

        INPUT:

        - ``comp`` -- either the components of the tensor field with respect
          to the vector frame specified by the argument ``frame`` or a
          dictionary of components, the keys of which are vector frames or
          pairs ``(f,c)`` where ``f`` is a vector frame and ``c`` a chart
        - ``frame`` -- (default: ``None``; unused if ``comp`` is a dictionary)
          vector frame in which the components are given; if ``None``, the
          default vector frame on the domain of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is a dictionary)
          coordinate chart in which the components are expressed; if ``None``,
          the default chart on the domain of ``frame`` is assumed

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t._init_components([[1+x, x*y], [-2, y^2]])
            sage: t.display()
            t = (x + 1) ∂/∂x⊗dx + x*y ∂/∂x⊗dy - 2 ∂/∂y⊗dx + y^2 ∂/∂y⊗dy
            sage: Y.<u,v> = M.chart()
            sage: t._init_components([[2*u, 3*v], [u+v, -u]], frame=Y.frame(),
            ....:                    chart=Y)
            sage: t.display(Y)
            t = 2*u ∂/∂u⊗du + 3*v ∂/∂u⊗dv + (u + v) ∂/∂v⊗du - u ∂/∂v⊗dv
            sage: t._init_components({X.frame(): [[2*x, 1-y],[0, x]]})
            sage: t.display()
            t = 2*x ∂/∂x⊗dx + (-y + 1) ∂/∂x⊗dy + x ∂/∂y⊗dy
            sage: t._init_components({(Y.frame(), Y): [[2*u, 0],[v^3, u+v]]})
            sage: t.display(Y)
            t = 2*u ∂/∂u⊗du + v^3 ∂/∂v⊗du + (u + v) ∂/∂v⊗dv

        TESTS:

        Check that :issue:`29639` is fixed::

            sage: v = M.vector_field()
            sage: v._init_components(1/2, -1)
            sage: v.display()
            1/2 ∂/∂x - ∂/∂y
        """
        comp0 = comp[0]
        self._is_zero = False  # a priori
        if isinstance(comp0, dict):
            for frame, components in comp0.items():
                chart = None
                if isinstance(frame, tuple):
                    # frame is actually a pair (frame, chart):
                    frame, chart = frame
                self.add_comp(frame)[:, chart] = components
        elif isinstance(comp0, str):
            # For compatibility with previous use of tensor_field():
            self.set_name(comp0)
        else:
            if hasattr(comp0, '__len__') and hasattr(comp0, '__getitem__'):
                # comp0 is a list/vector of components
                # otherwise comp is the tuple of components in a specific frame
                comp = comp0
            frame = kwargs.get('frame')
            chart = kwargs.get('chart')
            self.add_comp(frame)[:, chart] = comp

    #### Simple accessors ####

    def domain(self) -> DifferentiableManifold:
        r"""
        Return the manifold on which ``self`` is defined.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: t = M.tensor_field(1,2)
            sage: t.domain()
            2-dimensional differentiable manifold M
            sage: U = M.open_subset('U', coord_def={c_xy: x<0})
            sage: h = t.restrict(U)
            sage: h.domain()
            Open subset U of the 2-dimensional differentiable manifold M
        """
        return self._domain

    def base_module(self) -> VectorFieldModule:
        r"""
        Return the vector field module on which ``self`` acts as a tensor.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`

        EXAMPLES:

        The module of vector fields on the 2-sphere as a "base module"::

            sage: M = Manifold(2, 'S^2')
            sage: t = M.tensor_field(0,2)
            sage: t.base_module()
            Module X(S^2) of vector fields on the 2-dimensional differentiable
             manifold S^2
            sage: t.base_module() is M.vector_field_module()
            True
            sage: XM = M.vector_field_module()
            sage: XM.an_element().base_module() is XM
            True
        """
        return self._vmodule

    def tensor_type(self) -> TensorType:
        r"""
        Return the tensor type of ``self``.

        OUTPUT:

        - pair `(k,l)`, where `k` is the contravariant rank and `l` is
          the covariant rank

        EXAMPLES::

            sage: M = Manifold(2, 'S^2')
            sage: t = M.tensor_field(1,2)
            sage: t.tensor_type()
            (1, 2)
            sage: v = M.vector_field()
            sage: v.tensor_type()
            (1, 0)
        """
        return self._tensor_type

    def tensor_rank(self):
        r"""
        Return the tensor rank of ``self``.

        OUTPUT:

        - integer `k+l`, where `k` is the contravariant rank and `l` is
          the covariant rank

        EXAMPLES::

            sage: M = Manifold(2, 'S^2')
            sage: t = M.tensor_field(1,2)
            sage: t.tensor_rank()
            3
            sage: v = M.vector_field()
            sage: v.tensor_rank()
            1
        """
        return self._tensor_rank

    def symmetries(self):
        r"""
        Print the list of symmetries and antisymmetries.

        EXAMPLES::

            sage: M = Manifold(2, 'S^2')
            sage: t = M.tensor_field(1,2)
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: t = M.tensor_field(1,2, sym=(1,2))
            sage: t.symmetries()
            symmetry: (1, 2);  no antisymmetry
            sage: t = M.tensor_field(2,2, sym=(0,1), antisym=(2,3))
            sage: t.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
            sage: t = M.tensor_field(2,2, antisym=[(0,1),(2,3)])
            sage: t.symmetries()
            no symmetry;  antisymmetries: [(0, 1), (2, 3)]
        """
        if not self._sym:
            s = "no symmetry; "
        elif len(self._sym) == 1:
            s = "symmetry: {}; ".format(self._sym[0])
        else:
            s = "symmetries: {}; ".format(list(self._sym))
        if not self._antisym:
            a = "no antisymmetry"
        elif len(self._antisym) == 1:
            a = "antisymmetry: {}".format(self._antisym[0])
        else:
            a = "antisymmetries: {}".format(list(self._antisym))
        print(s + a)

    #### End of simple accessors #####

    def set_immutable(self):
        r"""
        Set ``self`` and all restrictions of ``self`` immutable.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1})
            sage: a = M.tensor_field(1, 1, [[1+y,x], [0,x+y]], name='a')
            sage: aU = a.restrict(U)
            sage: a.set_immutable()
            sage: aU.is_immutable()
            True
        """
        for rst in self._restrictions.values():
            rst.set_immutable()
        super().set_immutable()

    def set_restriction(self, rst: TensorField):
        r"""
        Define a restriction of ``self`` to some subdomain.

        INPUT:

        - ``rst`` -- :class:`TensorField` of the same type and symmetries
          as the current tensor field ``self``, defined on a subdomain of
          the domain of ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: t = M.tensor_field(1, 2, name='t')
            sage: s = U.tensor_field(1, 2)
            sage: s[0,0,1] = x+y
            sage: t.set_restriction(s)
            sage: t.display(c_xy.frame())
            t = (x + y) ∂/∂x⊗dx⊗dy
            sage: t.restrict(U) == s
            True

        If the restriction is defined on the very same domain, the tensor field
        becomes a copy of it (see :meth:`copy_from`)::

            sage: v = M.tensor_field(1, 2, name='v')
            sage: v.set_restriction(t)
            sage: v.restrict(U) == t.restrict(U)
            True
        """
        if self.is_immutable():
            raise ValueError("the restrictions of an immutable element "
                             "cannot be changed")
        if not isinstance(rst, TensorField):
            raise TypeError("the argument must be a tensor field")
        if not rst._domain.is_subset(self._domain):
            raise ValueError("the domain of the declared restriction is not " +
                             "a subset of the field's domain")
        if not rst._ambient_domain.is_subset(self._ambient_domain):
            raise ValueError("the ambient domain of the declared " +
                             "restriction is not a subset of the " +
                             "field's ambient domain")
        if rst._tensor_type != self._tensor_type:
            raise ValueError("the declared restriction has not the same " +
                             "tensor type as the current tensor field")
        if rst._tensor_type != self._tensor_type:
            raise ValueError("the declared restriction has not the same " +
                             "tensor type as the current tensor field")
        if rst._sym != self._sym:
            raise ValueError("the declared restriction has not the same " +
                             "symmetries as the current tensor field")
        if rst._antisym != self._antisym:
            raise ValueError("the declared restriction has not the same " +
                             "antisymmetries as the current tensor field")
        if self._domain is rst._domain:
            self.copy_from(rst)
        else:
            self._restrictions[rst._domain] = rst.copy(name=self._name,
                                                       latex_name=self._latex_name)
        self._is_zero = False  # a priori

    def restrict(
        self: T, subdomain: DifferentiableManifold, dest_map: Optional[DiffMap] = None
    ) -> T:
        r"""
        Return the restriction of ``self`` to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` --
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
          open subset `U` of the tensor field domain `S`
        - ``dest_map`` --
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          (default: ``None``); destination map `\Psi:\ U \rightarrow V`,
          where `V` is an open subset of the manifold `M` where the tensor
          field takes it values; if ``None``, the restriction of `\Phi`
          to `U` is used, `\Phi` being the differentiable map
          `S \rightarrow M` associated with the tensor field

        OUTPUT: :class:`TensorField` representing the restriction

        EXAMPLES:

        Restrictions of a vector field on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of the North pole
            sage: stereoN.<x,y> = U.chart()  # stereographic coordinates from the North pole
            sage: eN = stereoN.frame() # the associated vector frame
            sage: V =  M.open_subset('V') # the complement of the South pole
            sage: stereoS.<u,v> = V.chart()  # stereographic coordinates from the South pole
            sage: eS = stereoS.frame() # the associated vector frame
            sage: transf = stereoN.transition_map(stereoS, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:               intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:               restrictions2= u^2+v^2!=0)
            sage: inv = transf.inverse() # transformation from stereoS to stereoN
            sage: W = U.intersection(V) # the complement of the North and South poles
            sage: stereoN_W = W.atlas()[0]  # restriction of stereographic coord. from North pole to W
            sage: stereoS_W = W.atlas()[1]  # restriction of stereographic coord. from South pole to W
            sage: eN_W = stereoN_W.frame() ; eS_W = stereoS_W.frame()
            sage: v = M.vector_field({eN: [1, 0]}, name='v')
            sage: v.display()
            v = ∂/∂x
            sage: vU = v.restrict(U) ; vU
            Vector field v on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: vU.display()
            v = ∂/∂x
            sage: vU == eN[1]
            True
            sage: vW = v.restrict(W) ; vW
            Vector field v on the Open subset W of the 2-dimensional
             differentiable manifold S^2
            sage: vW.display()
            v = ∂/∂x
            sage: vW.display(eS_W, stereoS_W)
            v = (-u^2 + v^2) ∂/∂u - 2*u*v ∂/∂v
            sage: vW == eN_W[1]
            True

        At this stage, defining the restriction of ``v`` to the open
        subset ``V`` fully specifies ``v``::

            sage: v.restrict(V)[1] = vW[eS_W, 1, stereoS_W].expr()  # note that eS is the default frame on V
            sage: v.restrict(V)[2] = vW[eS_W, 2, stereoS_W].expr()
            sage: v.display(eS, stereoS)
            v = (-u^2 + v^2) ∂/∂u - 2*u*v ∂/∂v
            sage: v.restrict(U).display()
            v = ∂/∂x
            sage: v.restrict(V).display()
            v = (-u^2 + v^2) ∂/∂u - 2*u*v ∂/∂v

        The restriction of the vector field to its own domain is of course
        itself::

            sage: v.restrict(M) is v
            True
            sage: vU.restrict(U) is vU
            True
        """
        if (subdomain == self._domain
                and (dest_map is None or dest_map == self._vmodule._dest_map)):
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subset of " +
                                 "the field's domain")
            if dest_map is None:
                dest_map = self._vmodule._dest_map.restrict(subdomain)
            elif not dest_map._codomain.is_subset(self._ambient_domain):
                raise ValueError("the argument 'dest_map' is not compatible " +
                                 "with the ambient domain of " +
                                 "the {}".format(self))
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and subdomain in rst._restrictions:
                    res = rst._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and dom is not self._domain:
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    self._restrictions_graph[subdomain] = rst.restrict(subdomain)
                    return self._restrictions[subdomain]

            # Secondly one tries to get the restriction from one previously
            # defined on a larger domain:
            for dom, ext in self._extensions_graph.items():
                if subdomain in ext._restrictions:
                    res = ext._restrictions_graph[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

            # If this fails, the restriction is created from scratch:
            smodule = subdomain.vector_field_module(dest_map=dest_map)
            res = smodule.tensor(self._tensor_type, name=self._name,
                                 latex_name=self._latex_name, sym=self._sym,
                                 antisym=self._antisym,
                                 specific_type=type(self))
            res._extensions_graph.update(self._extensions_graph)
            for dom, ext in self._extensions_graph.items():
                ext._restrictions[subdomain] = res
                ext._restrictions_graph[subdomain] = res

            for dom, rst in self._restrictions.items():
                if dom.is_subset(subdomain):
                    if rst is not res:
                        res._restrictions.update(rst._restrictions)
                    res._restrictions_graph.update(rst._restrictions_graph)
                    rst._extensions_graph.update(res._extensions_graph)
            if self.is_immutable():
                res.set_immutable()  # restrictions must be immutable, too
            self._restrictions[subdomain] = res
            self._restrictions_graph[subdomain] = res
            res._extensions_graph.update(self._extensions_graph)

        return self._restrictions[subdomain]

    def _set_comp_unsafe(self, basis=None):
        r"""
        Return the components of ``self`` in a given vector frame
        for assignment. This private method invokes no security check. Use
        this method at your own risk.

        The components with respect to other frames having the same domain
        as the provided vector frame are deleted, in order to avoid any
        inconsistency. To keep them, use the method :meth:`_add_comp_unsafe`
        instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if none is provided, the components are
          assumed to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as a
          :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 2, name='t')
            sage: t._set_comp_unsafe(e_uv)
            3-indices components w.r.t. Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: t._set_comp_unsafe(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) ∂/∂v⊗du⊗dv

        Setting the components in a new frame (``e``)::

            sage: e = V.vector_frame('e')
            sage: t._set_comp_unsafe(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t._set_comp_unsafe(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0⊗e^1⊗e^1

        Since the frames ``e`` and ``e_uv`` are defined on the same domain, the
        components w.r.t. ``e_uv`` have been erased::

            sage: t.display(c_uv.frame())
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (V, (∂/∂u,∂/∂v))
        """
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst._set_comp_unsafe(basis)

    def set_comp(self, basis=None):
        r"""
        Return the components of ``self`` in a given vector frame
        for assignment.

        The components with respect to other frames having the same domain
        as the provided vector frame are deleted, in order to avoid any
        inconsistency. To keep them, use the method :meth:`add_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if none is provided, the components are
          assumed to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as a
          :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        EXAMPLES::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 2, name='t')
            sage: t.set_comp(e_uv)
            3-indices components w.r.t. Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: t.set_comp(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) ∂/∂v⊗du⊗dv

        Setting the components in a new frame (``e``)::

            sage: e = V.vector_frame('e')
            sage: t.set_comp(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t.set_comp(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0⊗e^1⊗e^1

        Since the frames ``e`` and ``e_uv`` are defined on the same domain, the
        components w.r.t. ``e_uv`` have been erased::

            sage: t.display(c_uv.frame())
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (V, (∂/∂u,∂/∂v))

        Since zero is an immutable, its components cannot be changed::

            sage: z = M.tensor_field_module((1, 1)).zero()
            sage: z.set_comp(e)[0,1] = u*v
            Traceback (most recent call last):
            ...
            ValueError: the components of an immutable element cannot be
             changed
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        self._is_zero = False  # a priori
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.set_comp(basis=basis)

    def _add_comp_unsafe(self, basis=None):
        r"""
        Return the components of ``self`` in a given vector frame
        for assignment. This private method invokes no security check. Use
        this method at your own risk.

        The components with respect to other frames having the same domain
        as the provided vector frame are kept. To delete them, use the
        method :meth:`_set_comp_unsafe` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if ``None``, the components are assumed
          to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as a
          :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 2, name='t')
            sage: t._add_comp_unsafe(e_uv)
            3-indices components w.r.t. Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: t._add_comp_unsafe(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) ∂/∂v⊗du⊗dv

        Setting the components in a new frame::

            sage: e = V.vector_frame('e')
            sage: t._add_comp_unsafe(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t._add_comp_unsafe(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0⊗e^1⊗e^1

        The components with respect to ``e_uv`` are kept::

            sage: t.display(e_uv)
            t = (u + v) ∂/∂v⊗du⊗dv
        """
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst._add_comp_unsafe(basis)

    def add_comp(self, basis=None) -> Components:
        r"""
        Return the components of ``self`` in a given vector frame
        for assignment.

        The components with respect to other frames having the same domain
        as the provided vector frame are kept. To delete them, use the
        method :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if ``None``, the components are assumed
          to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as a
          :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        EXAMPLES::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 2, name='t')
            sage: t.add_comp(e_uv)
            3-indices components w.r.t. Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: t.add_comp(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) ∂/∂v⊗du⊗dv

        Setting the components in a new frame::

            sage: e = V.vector_frame('e')
            sage: t.add_comp(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t.add_comp(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0⊗e^1⊗e^1

        The components with respect to ``e_uv`` are kept::

            sage: t.display(e_uv)
            t = (u + v) ∂/∂v⊗du⊗dv

        Since zero is a special element, its components cannot be changed::

            sage: z = M.tensor_field_module((1, 1)).zero()
            sage: z.add_comp(e_uv)[1, 1] = u^2
            Traceback (most recent call last):
            ...
            ValueError: the components of an immutable element cannot be
             changed
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        self._is_zero = False  # a priori
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.add_comp(basis=basis)

    def add_comp_by_continuation(self, frame, subdomain, chart=None):
        r"""
        Set components with respect to a vector frame by continuation of the
        coordinate expression of the components in a subframe.

        The continuation is performed by demanding that the components have
        the same coordinate expression as those on the restriction of the
        frame to a given subdomain.

        INPUT:

        - ``frame`` -- vector frame `e` in which the components are to be set
        - ``subdomain`` -- open subset of `e`'s domain in which the
          components are known or can be evaluated from other components
        - ``chart`` -- (default: ``None``) coordinate chart on `e`'s domain in
          which the extension of the expression of the components is to be
          performed; if ``None``, the default's chart of `e`'s domain is
          assumed

        EXAMPLES:

        Components of a vector field on the sphere `S^2`::

            sage: M = Manifold(2, 'S^2', start_index=1)

        The two open subsets covered by stereographic coordinates (North and South)::

            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coordinates
            sage: transf = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:             intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:             restrictions2= u^2+v^2!=0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.vector_field({eU: [x, 2+y]}, name='a')

        At this stage, the vector field has been defined only on the open
        subset ``U`` (through its components in the frame ``eU``)::

            sage: a.display(eU)
            a = x ∂/∂x + (y + 2) ∂/∂y

        The components with respect to the restriction of ``eV`` to the common
        subdomain ``W``, in terms of the ``(u,v)`` coordinates, are obtained
        by a change-of-frame formula on ``W``::

            sage: a.display(eV.restrict(W), c_uv.restrict(W))
            a = (-4*u*v - u) ∂/∂u + (2*u^2 - 2*v^2 - v) ∂/∂v

        The continuation consists in extending the definition of the vector
        field to the whole open subset ``V`` by demanding that the components
        in the frame eV have the same coordinate expression as the above one::

            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)

        We have then::

            sage: a.display(eV)
            a = (-4*u*v - u) ∂/∂u + (2*u^2 - 2*v^2 - v) ∂/∂v

        and `a` is defined on the entire manifold `S^2`.
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        dom = frame._domain
        if not dom.is_subset(self._domain):
            raise ValueError("the vector frame is not defined on a subset " +
                             "of the tensor field domain")
        if chart is None:
            chart = dom._def_chart
        sframe = frame.restrict(subdomain)
        schart = chart.restrict(subdomain)
        scomp = self.comp(sframe)
        resu = self._add_comp_unsafe(frame) # _del_derived is performed here
        for ind in resu.non_redundant_index_generator():
            resu[[ind]] = dom.scalar_field({chart: scomp[[ind]].expr(schart)})

    def add_expr_from_subdomain(self, frame, subdomain):
        r"""
        Add an expression to an existing component from a subdomain.

        INPUT:

        - ``frame`` -- vector frame `e` in which the components are to be set
        - ``subdomain`` -- open subset of `e`'s domain in which the
          components have additional expressions

        EXAMPLES:

        We are going to consider a vector field in `\RR^3` along the 2-sphere::

            sage: M = Manifold(3, 'M', structure='Riemannian')
            sage: S = Manifold(2, 'S', structure='Riemannian')
            sage: E.<X,Y,Z> = M.chart()

        Let us define ``S`` in terms of stereographic charts::

            sage: U = S.open_subset('U')
            sage: V = S.open_subset('V')
            sage: S.declare_union(U,V)
            sage: stereoN.<x,y> = U.chart()
            sage: stereoS.<xp,yp> = V.chart("xp:x' yp:y'")
            sage: stereoN_to_S = stereoN.transition_map(stereoS,
            ....:                                 (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                 intersection_name='W',
            ....:                                 restrictions1= x^2+y^2!=0,
            ....:                                 restrictions2= xp^2+yp^2!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: W = U.intersection(V)
            sage: stereoN_W = stereoN.restrict(W)
            sage: stereoS_W = stereoS.restrict(W)

        The embedding of `S^2` in `\RR^3`::

            sage: phi = S.diff_map(M, {(stereoN, E): [2*x/(1+x^2+y^2),
            ....:                                     2*y/(1+x^2+y^2),
            ....:                                     (x^2+y^2-1)/(1+x^2+y^2)],
            ....:                        (stereoS, E): [2*xp/(1+xp^2+yp^2),
            ....:                                       2*yp/(1+xp^2+yp^2),
            ....:                               (1-xp^2-yp^2)/(1+xp^2+yp^2)]},
            ....:                   name='Phi', latex_name=r'\Phi')

        To define a vector field ``v`` along ``S`` taking its values in ``M``,
        we first set the components on ``U``::

            sage: v = M.vector_field(name='v').along(phi)
            sage: vU = v.restrict(U)
            sage: vU[:] = [x,y,x**2+y**2]

        But because ``M`` is parallelizable, these components can be extended
        to ``S`` itself::

            sage: v.add_comp_by_continuation(E.frame().along(phi), U)

        One can see that ``v`` is not yet fully defined: the components
        (scalar fields) do not have values on the whole manifold::

            sage: sorted(v._components.values())[0]._comp[(0,)].display()
            S → ℝ
            on U: (x, y) ↦ x
            on W: (xp, yp) ↦ xp/(xp^2 + yp^2)

        To fix that, we first extend the components from ``W`` to ``V`` using
        :meth:`add_comp_by_continuation`::

            sage: v.add_comp_by_continuation(E.frame().along(phi).restrict(V),
            ....:                            W, stereoS)

        Then, the expression on the subdomain ``V`` is added to the
        already known components on ``S`` by::

            sage: v.add_expr_from_subdomain(E.frame().along(phi), V)

        The definition of ``v`` is now complete::

            sage: sorted(v._components.values())[0]._comp[(2,)].display()
            S → ℝ
            on U: (x, y) ↦ x^2 + y^2
            on V: (xp, yp) ↦ 1/(xp^2 + yp^2)
        """
        if self.is_immutable():
            raise ValueError("the expressions of an immutable element "
                             "cannot be changed")
        dom = frame._domain
        if not dom.is_subset(self._domain):
            raise ValueError("the vector frame is not defined on a subset " +
                             "of the tensor field domain")
        if frame not in self.restrict(frame.domain())._components:
            raise ValueError("the tensor doesn't have an expression in "
                             "the frame"+frame._repr_())
        comp = self._add_comp_unsafe(frame)  # the components stay the same
        scomp = self.restrict(subdomain).comp(frame.restrict(subdomain))
        for ind in comp.non_redundant_index_generator():
            comp[[ind]]._express.update(scomp[[ind]]._express)

        rst = self._restrictions.copy()
        self._del_derived()  # may delete restrictions
        self._restrictions = rst

    def comp(self, basis=None, from_basis=None):
        r"""
        Return the components in a given vector frame.

        If the components are not known already, they are computed by the
        tensor change-of-basis formula from components in another vector frame.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are required; if none is provided, the components are assumed to
          refer to the tensor field domain's default frame
        - ``from_basis`` -- (default: ``None``) vector frame from which the
          required components are computed, via the tensor change-of-basis
          formula, if they are not known already in the basis ``basis``

        OUTPUT:

        - components in the vector frame ``basis``, as a
          :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES:

        Components of a type-`(1,1)` tensor field defined on two
        open subsets::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: e = U.default_frame() ; e
            Coordinate frame (U, (∂/∂x,∂/∂y))
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: f = V.default_frame() ; f
            Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e,0,0] = - x + y^3
            sage: t[e,0,1] = 2+x
            sage: t[f,1,1] = - u*v
            sage: t.comp(e)
            2-indices components w.r.t. Coordinate frame (U, (∂/∂x,∂/∂y))
            sage: t.comp(e)[:]
            [y^3 - x   x + 2]
            [      0       0]
            sage: t.comp(f)
            2-indices components w.r.t. Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: t.comp(f)[:]
            [   0    0]
            [   0 -u*v]

        Since ``e`` is ``M``'s default frame, the argument ``e`` can
        be omitted::

            sage: e is M.default_frame()
            True
            sage: t.comp() is t.comp(e)
            True

        Example of computation of the components via a change of frame::

            sage: a = V.automorphism_field()
            sage: a[:] = [[1+v, -u^2], [0, 1-u]]
            sage: h = f.new_frame(a, 'h')
            sage: t.comp(h)
            2-indices components w.r.t. Vector frame (V, (h_0,h_1))
            sage: t.comp(h)[:]
            [             0 -u^3*v/(v + 1)]
            [             0           -u*v]
        """
        if basis is None:
            basis = self._domain._def_frame
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.comp(basis=basis, from_basis=from_basis)

    def display(self, frame=None, chart=None):
        r"""
        Display the tensor field in terms of its expansion with respect
        to a given vector frame.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to
          which the tensor is expanded; if ``frame`` is ``None`` and ``chart``
          is not ``None``, the coordinate frame associated with ``chart`` is
          assumed; if both ``frame`` and ``chart`` are ``None``, the default
          frame of the domain of definition of the tensor field is assumed
        - ``chart`` -- (default: ``None``) chart with respect to which the
          components of the tensor field in the selected frame are expressed;
          if ``None``, the default chart of the vector frame domain is assumed

        EXAMPLES:

        Display of a type-`(1,1)` tensor field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e_xy,:] = [[x, 1], [y, 0]]
            sage: t.add_comp_by_continuation(e_uv, W, c_uv)
            sage: t.display(e_xy)
            t = x ∂/∂x⊗dx + ∂/∂x⊗dy + y ∂/∂y⊗dx
            sage: t.display(e_uv)
            t = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv
              + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv

        Since ``e_xy`` is ``M``'s default frame, the argument ``e_xy`` can
        be omitted::

            sage: e_xy is M.default_frame()
            True
            sage: t.display()
            t = x ∂/∂x⊗dx + ∂/∂x⊗dy + y ∂/∂y⊗dx

        Similarly, since ``e_uv`` is ``V``'s default frame, the argument ``e_uv``
        can be omitted when considering the restriction of ``t`` to ``V``::

            sage: t.restrict(V).display()
            t = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv
              + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv

        If the coordinate expression of the components are to be displayed in
        a chart distinct from the default one on the considered domain, then
        the chart has to be passed as the second argument of ``display``.
        For instance, on `W = U \cap V`, two charts are available:
        ``c_xy.restrict(W)`` (the default one) and ``c_uv.restrict(W)``.
        Accordingly, one can have two views of the expansion of ``t`` in the
        *same* vector frame ``e_uv.restrict(W)``::

            sage: t.display(e_uv.restrict(W))  # W's default chart assumed
            t = (1/2*x + 1/2*y + 1/2) ∂/∂u⊗du + (1/2*x + 1/2*y - 1/2) ∂/∂u⊗dv
              + (1/2*x - 1/2*y + 1/2) ∂/∂v⊗du + (1/2*x - 1/2*y - 1/2) ∂/∂v⊗dv
            sage: t.display(e_uv.restrict(W), c_uv.restrict(W))
            t = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv
              + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv

        As a shortcut, one can pass just a chart to ``display``. It is then
        understood that the expansion is to be performed with respect to the
        coordinate frame associated with this chart. Therefore the above
        command can be abridged to::

            sage: t.display(c_uv.restrict(W))
            t = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv
              + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv

        and one has::

            sage: t.display(c_xy)
            t = x ∂/∂x⊗dx + ∂/∂x⊗dy + y ∂/∂y⊗dx
            sage: t.display(c_uv)
            t = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv
              + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv
            sage: t.display(c_xy.restrict(W))
            t = x ∂/∂x⊗dx + ∂/∂x⊗dy + y ∂/∂y⊗dx
            sage: t.restrict(W).display(c_uv.restrict(W))
            t = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv
              + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv

        One can ask for the display with respect to a frame in which ``t`` has
        not been initialized yet (this will automatically trigger the use of
        the change-of-frame formula for tensors)::

            sage: a = V.automorphism_field()
            sage: a[:] = [[1+v, -u^2], [0, 1-u]]
            sage: f = e_uv.new_frame(a, 'f')
            sage: [f[i].display() for i in M.irange()]
            [f_0 = (v + 1) ∂/∂u, f_1 = -u^2 ∂/∂u + (-u + 1) ∂/∂v]
            sage: t.display(f)
            t = -1/2*(u^2*v + 1)/(u - 1) f_0⊗f^0
              - 1/2*(2*u^3 - 5*u^2 - (u^4 + u^3 - u^2)*v + 3*u - 1)/((u - 1)*v + u - 1) f_0⊗f^1
              - 1/2*(v^2 + 2*v + 1)/(u - 1) f_1⊗f^0
              + 1/2*(u^2 + (u^2 + u - 1)*v - u + 1)/(u - 1) f_1⊗f^1

        A shortcut of ``display()`` is ``disp()``::

            sage: t.disp(e_uv)
            t = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv
              + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv
        """
        if frame is None:
            if chart is not None:
                frame = chart.frame()
            else:
                if self._vmodule._dest_map.is_identity():
                    frame = self._domain._def_frame
                else:
                    for rst in self._restrictions.values():
                        try:
                            return rst.display()
                        except ValueError:
                            pass
                if frame is None:  # should be "is still None" ;-)
                    raise ValueError("a frame must be provided for the display")
        else:
            try:
                frame0 = frame.frame()
                # if this succeeds, frame is actually not a vector frame, but
                # a coordinate chart
                if chart is None:
                    chart = frame
                frame = frame0
            except AttributeError:
                # case of a genuine vector frame
                pass
        rst = self.restrict(frame._domain, dest_map=frame._dest_map)
        return rst.display(frame, chart)

    disp = display

    def display_comp(self, frame=None, chart=None, coordinate_labels=True,
                     only_nonzero=True, only_nonredundant=False):
        r"""
        Display the tensor components with respect to a given frame,
        one per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to which
          the tensor field components are defined; if ``None``, then

          * if ``chart`` is not ``None``, the coordinate frame associated to
            ``chart`` is used
          * otherwise, the default basis of the vector field module on which
            the tensor field is defined is used

        - ``chart`` -- (default: ``None``) chart specifying the coordinate
          expression of the components; if ``None``, the default chart of the
          tensor field domain is used
        - ``coordinate_labels`` -- boolean (default: ``True``); if ``True``,
          coordinate symbols are used by default (instead of integers) as
          index labels whenever ``frame`` is a coordinate frame
        - ``only_nonzero`` -- boolean (default: ``True``); if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- boolean (default: ``False``); if ``True``,
          only nonredundant components are displayed in case of symmetries

        EXAMPLES:

        Display of the components of a type-`(1,1)` tensor field defined
        on two open subsets::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: e = U.default_frame()
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: f = V.default_frame()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e,0,0] = - x + y^3
            sage: t[e,0,1] = 2+x
            sage: t[f,1,1] = - u*v
            sage: t.display_comp(e)
            t^x_x = y^3 - x
            t^x_y = x + 2
            sage: t.display_comp(f)
            t^v_v = -u*v

        Components in a chart frame::

            sage: t.display_comp(chart=c_xy)
            t^x_x = y^3 - x
            t^x_y = x + 2
            sage: t.display_comp(chart=c_uv)
            t^v_v = -u*v

        See documentation of
        :meth:`sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.display_comp`
        for more options.
        """
        if frame is None:
            if chart is not None:
                frame = chart.frame()
            else:
                if self._vmodule._dest_map.is_identity():
                    frame = self._domain.default_frame()
                else:
                    for rst in self._restrictions.values():
                        try:
                            return rst.display_comp(chart=chart,
                                       coordinate_labels=coordinate_labels,
                                       only_nonzero=only_nonzero,
                                       only_nonredundant=only_nonredundant)
                        except ValueError:
                            pass
                if frame is None:  # should be "is still None" ;-)
                    raise ValueError("a frame must be provided for the display")
        rst = self.restrict(frame.domain(), dest_map=frame._dest_map)
        return rst.display_comp(frame=frame, chart=chart,
                                coordinate_labels=coordinate_labels,
                                only_nonzero=only_nonzero,
                                only_nonredundant=only_nonredundant)

    def __getitem__(self, args):
        r"""
        Return a component with respect to some frame.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are returned

        The frame can be passed as the first item of ``args``. If not, the
        default frame of the tensor field's domain is assumed. If ``args``
        is a string, this method acts as a shortcut for tensor contractions
        and symmetrizations, the string containing abstract indices.

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_xy = c_xy.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy, :] = [[x+y, -2], [3*y^2, x*y]]
            sage: t.__getitem__((1,0))
            3*y^2
            sage: t.__getitem__((1,1))
            x*y
            sage: t.__getitem__((e_xy,1,0))
            3*y^2
            sage: t.__getitem__(slice(None))
            [x + y    -2]
            [3*y^2   x*y]
            sage: t.__getitem__((e_xy,slice(None)))
            [x + y    -2]
            [3*y^2   x*y]
            sage: t.__getitem__('^a_a')  # trace
            Scalar field on the 2-dimensional differentiable manifold M
            sage: t.__getitem__('^a_a').display()
            M → ℝ
            on U: (x, y) ↦ (x + 1)*y + x
        """
        if isinstance(args, str): # tensor with specified indices
            return TensorWithIndices(self, args).update()
        if isinstance(args, list):  # case of [[...]] syntax
            if not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._domain._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        return self.comp(frame)[args]

    def __setitem__(self, args, value):
        r"""
        Set a component with respect to some vector frame.

        INPUT:

        - ``args`` -- list of indices; if ``[:]`` is provided, all the
          components are set; the frame can be passed as the first item
          of ``args``; if not, the default frame of the tensor field's
          domain is assumed
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_xy = c_xy.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t.__setitem__((e_xy, 0, 1), x+y^2)
            sage: t.display(e_xy)
            t = (y^2 + x) ∂/∂x⊗dy
            sage: t.__setitem__((0, 1), x+y^2)  # same as above since e_xy is the default frame on M
            sage: t.display()
            t = (y^2 + x) ∂/∂x⊗dy
            sage: t.__setitem__(slice(None), [[x+y, -2], [3*y^2, x*y]])
            sage: t.display()
            t = (x + y) ∂/∂x⊗dx - 2 ∂/∂x⊗dy + 3*y^2 ∂/∂y⊗dx + x*y ∂/∂y⊗dy
        """
        if isinstance(args, list):  # case of [[...]] syntax
            if not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._domain._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        self.set_comp(frame)[args] = value

    def copy_from(self, other):
        r"""
        Make ``self`` a copy of ``other``.

        INPUT:

        - ``other`` -- other tensor field, in the same module as ``self``

        .. NOTE::

            While the derived quantities are not copied, the name is kept.

        .. WARNING::

            All previous defined components and restrictions will be deleted!

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = M.tensor_field(1, 1, name='s')
            sage: s.copy_from(t)
            sage: s.display(e_xy)
            s = (x + y) ∂/∂x⊗dx + 2 ∂/∂y⊗dx + (-y + 1) ∂/∂y⊗dy
            sage: s == t
            True

        While the original tensor field is modified, the copy is not::

            sage: t[e_xy,0,0] = -1
            sage: t.display(e_xy)
            t = -∂/∂x⊗dx + 2 ∂/∂y⊗dx + (-y + 1) ∂/∂y⊗dy
            sage: s.display(e_xy)
            s = (x + y) ∂/∂x⊗dx + 2 ∂/∂y⊗dx + (-y + 1) ∂/∂y⊗dy
            sage: s == t
            False
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        if other not in self.parent():
            raise TypeError("the original must be an element of "
                            f"{self.parent()}")
        self._del_derived()
        self._del_restrictions() # delete restrictions
        for dom, rst in other._restrictions.items():
            self._restrictions[dom] = rst.copy(name=self._name,
                                               latex_name=self._latex_name)
        self._is_zero = other._is_zero

    def copy(self, name=None, latex_name=None):
        r"""
        Return an exact copy of ``self``.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the copy
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          copy; if none is provided, the LaTeX symbol is set to ``name``

        .. NOTE::

            The name and the derived quantities are not copied.

        EXAMPLES:

        Copy of a type-`(1,1)` tensor field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = t.copy(); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            (x + y) ∂/∂x⊗dx + 2 ∂/∂y⊗dx + (-y + 1) ∂/∂y⊗dy
            sage: s == t
            True

        If the original tensor field is modified, the copy is not::

            sage: t[e_xy,0,0] = -1
            sage: t.display(e_xy)
            t = -∂/∂x⊗dx + 2 ∂/∂y⊗dx + (-y + 1) ∂/∂y⊗dy
            sage: s.display(e_xy)
            (x + y) ∂/∂x⊗dx + 2 ∂/∂y⊗dx + (-y + 1) ∂/∂y⊗dy
            sage: s == t
            False
        """
        resu = self._new_instance()
        # set resu name
        if name is not None:
            resu._name = name
            if latex_name is None:
                resu._latex_name = name
        if latex_name is not None:
            resu._latex_name = latex_name
        # set restrictions
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = rst.copy(name=name,
                                               latex_name=latex_name)
        resu._is_zero = self._is_zero
        return resu

    def _common_subdomains(self, other):
        r"""
        Return the list of subdomains of ``self._domain`` on which
        both ``self`` and ``other`` have known restrictions.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 0]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: sorted(t._common_subdomains(t), key=str)
            [Open subset U of the 2-dimensional differentiable manifold M,
             Open subset V of the 2-dimensional differentiable manifold M,
             Open subset W of the 2-dimensional differentiable manifold M]
            sage: a = M.tensor_field(1, 1, name='a')
            sage: t._common_subdomains(a)
            []
            sage: a[e_xy, 0, 1] = 0
            sage: t._common_subdomains(a)
            [Open subset U of the 2-dimensional differentiable manifold M]
            sage: a[e_uv, 0, 0] = 0
            sage: sorted(t._common_subdomains(a), key=str)
            [Open subset U of the 2-dimensional differentiable manifold M,
             Open subset V of the 2-dimensional differentiable manifold M]
        """
        resu = []
        for dom in self._restrictions:
            if dom in other._restrictions:
                resu.append(dom)
        return resu

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a tensor field or 0

        OUTPUT: ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: t == t
            True
            sage: t == t.copy()
            True
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a.set_restriction(t.restrict(U))
            sage: t == a  # False since a has not been defined on V
            False
            sage: a.set_restriction(t.restrict(V))
            sage: t == a  # True now
            True
            sage: a[e_xy, 0, 0] = -1
            sage: t == a  # False since a has been reset on U (domain of e_xy)
            False
            sage: t.parent().zero() == 0
            True
        """
        from .mixed_form import MixedForm

        if other is self:
            return True
        if other in ZZ: # to compare with 0
            if other == 0:
                return self.is_zero()
            return False
        elif isinstance(other, MixedForm):
            # use comparison of MixedForm:
            return other == self
        elif not isinstance(other, TensorField):
            return False
        else: # other is another tensor field
            if other._vmodule != self._vmodule:
                return False
            if other._tensor_type != self._tensor_type:
                return False
            # Non-trivial open covers of the domain:
            for oc in self._domain.open_covers(trivial=False):
                resu = True
                for dom in oc:
                    try:
                        resu = resu and \
                                bool(self.restrict(dom) == other.restrict(dom))
                    except ValueError:
                        break
                else:
                    # If this point is reached, no exception has occurred; hence
                    # the result is valid and can be returned:
                    return resu
            # If this point is reached, the comparison has not been possible
            # on any open cover; we then compare the restrictions to
            # subdomains:
            if not self._restrictions:
                return False  # self is not initialized
            if len(self._restrictions) != len(other._restrictions):
                return False  # the restrictions are not on the same subdomains
            resu = True
            for dom, rst in self._restrictions.items():
                if dom in other._restrictions:
                    resu = resu and bool(rst == other._restrictions[dom])
                else:
                    return False  # the restrictions are not on the same
                                  # subdomains
            return resu

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- a tensor field or 0

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: t != t
            False
            sage: t != t.copy()
            False
            sage: t != 0
            True
        """
        return not (self == other)

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT: an exact copy of ``self``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = t.__pos__(); s
            Tensor field +t of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            +t = (x + y) ∂/∂x⊗dx + 2 ∂/∂y⊗dx + (-y + 1) ∂/∂y⊗dy
        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = + rst
        # Compose names:
        from sage.tensor.modules.format_utilities import (format_unop_txt,
                                                          format_unop_latex)
        resu._name = format_unop_txt('+', self._name)
        resu._latex_name = format_unop_latex(r'+', self._latex_name)
        return resu

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT: the tensor field `-T`, where `T` is ``self``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy, :] = [[x, -x], [y, -y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: t.display(e_xy)
            t = x ∂/∂x⊗dx - x ∂/∂x⊗dy + y ∂/∂y⊗dx - y ∂/∂y⊗dy
            sage: t.display(e_uv)
            t = u ∂/∂u⊗dv + v ∂/∂v⊗dv
            sage: s = t.__neg__(); s
            Tensor field -t of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            -t = -x ∂/∂x⊗dx + x ∂/∂x⊗dy - y ∂/∂y⊗dx + y ∂/∂y⊗dy
            sage: s.display(e_uv)
            -t = -u ∂/∂u⊗dv - v ∂/∂v⊗dv
            sage: s == -t  # indirect doctest
            True
        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = - rst
        # Compose names:
        from sage.tensor.modules.format_utilities import (format_unop_txt,
                                                          format_unop_latex)
        resu._name = format_unop_txt('-', self._name)
        resu._latex = format_unop_latex(r'-', self._latex_name)
        return resu

    ######### ModuleElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Tensor field addition.

        INPUT:

        - ``other`` -- a tensor field, in the same tensor module as ``self``

        OUTPUT:

        - the tensor field resulting from the addition of ``self``
          and ``other``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.tensor_field(1, 1, name='b')
            sage: b[e_xy,:] = [[2, y], [x, -x]]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a._add_(b); s
            Tensor field a+b of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a.display(e_xy)
            a = x ∂/∂x⊗dx + ∂/∂x⊗dy + y ∂/∂y⊗dx
            sage: b.display(e_xy)
            b = 2 ∂/∂x⊗dx + y ∂/∂x⊗dy + x ∂/∂y⊗dx - x ∂/∂y⊗dy
            sage: s.display(e_xy)
            a+b = (x + 2) ∂/∂x⊗dx + (y + 1) ∂/∂x⊗dy + (x + y) ∂/∂y⊗dx - x ∂/∂y⊗dy
            sage: s == a + b  # indirect doctest
            True
            sage: z = a.parent().zero(); z
            Tensor field zero of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a._add_(z) == a
            True
            sage: z._add_(a) == a
            True
        """
        # Case zero:
        if self._is_zero:
            return other
        if other._is_zero:
            return self
        # Generic case:
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] + other._restrictions[dom]
        some_rst = next(iter(resu_rst.values()))
        resu_sym = some_rst._sym
        resu_antisym = some_rst._antisym
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_sym,
                                    antisym=resu_antisym)
        resu._restrictions = resu_rst
        if self._name is not None and other._name is not None:
            resu._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            resu._latex_name = self._latex_name + '+' + other._latex_name
        return resu

    def _sub_(self, other):
        r"""
        Tensor field subtraction.

        INPUT:

        - ``other`` -- a tensor field, in the same tensor module as ``self``

        OUTPUT:

        - the tensor field resulting from the subtraction of ``other``
          from ``self``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.tensor_field(1, 1, name='b')
            sage: b[e_xy,:] = [[2, y], [x, -x]]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a._sub_(b); s
            Tensor field a-b of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a.display(e_xy)
            a = x ∂/∂x⊗dx + ∂/∂x⊗dy + y ∂/∂y⊗dx
            sage: b.display(e_xy)
            b = 2 ∂/∂x⊗dx + y ∂/∂x⊗dy + x ∂/∂y⊗dx - x ∂/∂y⊗dy
            sage: s.display(e_xy)
            a-b = (x - 2) ∂/∂x⊗dx + (-y + 1) ∂/∂x⊗dy + (-x + y) ∂/∂y⊗dx + x ∂/∂y⊗dy
            sage: s == a - b
            True
            sage: z = a.parent().zero()
            sage: a._sub_(z) == a
            True
            sage: z._sub_(a) == -a
            True
        """
        # Case zero:
        if self._is_zero:
            return -other
        if other._is_zero:
            return self
        # Generic case:
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] - other._restrictions[dom]
        some_rst = next(iter(resu_rst.values()))
        resu_sym = some_rst._sym
        resu_antisym = some_rst._antisym
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_sym,
                                   antisym=resu_antisym)
        resu._restrictions = resu_rst
        if self._name is not None and other._name is not None:
            resu._name = self._name + '-' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            resu._latex_name = self._latex_name + '-' + other._latex_name
        return resu

    def _rmul_(self, scalar):
        r"""
        Reflected multiplication operator: performs ``scalar * self``.

        This is actually the multiplication by an element of the ring over
        which the tensor field module is constructed.

        INPUT:

        - ``scalar`` -- scalar field in the scalar field algebra over which
          the module containing ``self`` is defined

        OUTPUT: the tensor field ``scalar * self``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2)}, name='f')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: f.display()
            f: M → ℝ
            on U: (x, y) ↦ 1/(x^2 + y^2 + 1)
            on V: (u, v) ↦ 2/(u^2 + v^2 + 2)
            sage: s = a._rmul_(f); s
            Tensor field f*a of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a.display(e_xy)
            a = x ∂/∂x⊗dx + ∂/∂x⊗dy + y ∂/∂y⊗dx
            sage: s.display(e_xy)
            f*a = x/(x^2 + y^2 + 1) ∂/∂x⊗dx + 1/(x^2 + y^2 + 1) ∂/∂x⊗dy + y/(x^2 + y^2 + 1) ∂/∂y⊗dx
            sage: a.display(e_uv)
            a = (1/2*u + 1/2) ∂/∂u⊗du + (1/2*u - 1/2) ∂/∂u⊗dv + (1/2*v + 1/2) ∂/∂v⊗du + (1/2*v - 1/2) ∂/∂v⊗dv
            sage: s.display(e_uv)
            f*a = (u + 1)/(u^2 + v^2 + 2) ∂/∂u⊗du + (u - 1)/(u^2 + v^2 + 2) ∂/∂u⊗dv + (v + 1)/(u^2 + v^2 + 2) ∂/∂v⊗du + (v - 1)/(u^2 + v^2 + 2) ∂/∂v⊗dv
            sage: s == f*a  # indirect doctest
            True
            sage: z = a.parent().zero(); z
            Tensor field zero of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a._rmul_(M.zero_scalar_field()) == z
            True
            sage: z._rmul_(f) == z
            True
        """
        # Case zero:
        if scalar._is_zero:
            return self.parent().zero()
        # Case one:
        if scalar is self._domain._one_scalar_field:
            return self
        # Generic case:
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = scalar.restrict(dom) * rst
        # Compose names:
        from sage.tensor.modules.format_utilities import (format_mul_txt,
                                                          format_mul_latex)
        resu_name = format_mul_txt(scalar._name, '*', self._name)
        resu_latex = format_mul_latex(scalar._latex_name, r' \cdot ',
                                      self._latex_name)
        resu.set_name(name=resu_name, latex_name=resu_latex)
        return resu

    ######### End of ModuleElement arithmetic operators ########

    # TODO: Move to acted_upon or _rmul_
    def __mul__(self, other: TensorField) -> TensorField:
        r"""
        Tensor product (or multiplication of the right by a scalar).

        INPUT:

        - ``other`` -- tensor field on the same manifold as ``self`` (or an
          object that can be coerced to a scalar field on the same manifold
          as ``self``)

        OUTPUT:

        - the tensor field resulting from the tensor product of ``self``
          with ``other`` (or from the product ``other * self`` if ``other``
          is a scalar)

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)

        Tensor product with another tensor field::

            sage: b = M.vector_field(name='b')
            sage: b[e_xy, :] = [x, y]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a.__mul__(b); s
            Tensor field of type (2,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            a⊗b = x^2 ∂/∂x⊗∂/∂x⊗dx + x ∂/∂x⊗∂/∂x⊗dy + x*y ∂/∂x⊗∂/∂y⊗dx
             + y ∂/∂x⊗∂/∂y⊗dy + x*y ∂/∂y⊗∂/∂x⊗dx + y^2 ∂/∂y⊗∂/∂y⊗dx
            sage: s.display(e_uv)
            a⊗b = (1/2*u^2 + 1/2*u) ∂/∂u⊗∂/∂u⊗du + (1/2*u^2 - 1/2*u) ∂/∂u⊗∂/∂u⊗dv
             + 1/2*(u + 1)*v ∂/∂u⊗∂/∂v⊗du + 1/2*(u - 1)*v ∂/∂u⊗∂/∂v⊗dv
             + (1/2*u*v + 1/2*u) ∂/∂v⊗∂/∂u⊗du + (1/2*u*v - 1/2*u) ∂/∂v⊗∂/∂u⊗dv
             + (1/2*v^2 + 1/2*v) ∂/∂v⊗∂/∂v⊗du + (1/2*v^2 - 1/2*v) ∂/∂v⊗∂/∂v⊗dv

        Multiplication on the right by a scalar field::

            sage: f = M.scalar_field({c_xy: x*y}, name='f')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: s = a.__mul__(f); s
            Tensor field f*a of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            f*a = x^2*y ∂/∂x⊗dx + x*y ∂/∂x⊗dy + x*y^2 ∂/∂y⊗dx
            sage: s.display(e_uv)
            f*a = (1/8*u^3 - 1/8*(u + 1)*v^2 + 1/8*u^2) ∂/∂u⊗du
             + (1/8*u^3 - 1/8*(u - 1)*v^2 - 1/8*u^2) ∂/∂u⊗dv
             + (1/8*u^2*v - 1/8*v^3 + 1/8*u^2 - 1/8*v^2) ∂/∂v⊗du
             + (1/8*u^2*v - 1/8*v^3 - 1/8*u^2 + 1/8*v^2) ∂/∂v⊗dv
            sage: s == f*a
            True

        Multiplication on the right by a number::

            sage: s = a.__mul__(2); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            2*x ∂/∂x⊗dx + 2 ∂/∂x⊗dy + 2*y ∂/∂y⊗dx
            sage: s.display(e_uv)
            (u + 1) ∂/∂u⊗du + (u - 1) ∂/∂u⊗dv + (v + 1) ∂/∂v⊗du
             + (v - 1) ∂/∂v⊗dv
            sage: s.restrict(U) == 2*a.restrict(U)
            True
            sage: s.restrict(V) == 2*a.restrict(V)
            True
            sage: s == 2*a
            True

        Test with SymPy as calculus engine::

            sage: M.set_calculus_method('sympy')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: s = a.__mul__(f); s
            Tensor field f*a of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            f*a = x**2*y ∂/∂x⊗dx + x*y ∂/∂x⊗dy + x*y**2 ∂/∂y⊗dx
            sage: s.display(e_uv)
            f*a = (u**3/8 + u**2/8 - u*v**2/8 - v**2/8) ∂/∂u⊗du + (u**3/8 -
            u**2/8 - u*v**2/8 + v**2/8) ∂/∂u⊗dv + (u**2*v/8 + u**2/8 -
            v**3/8 - v**2/8) ∂/∂v⊗du + (u**2*v/8 - u**2/8 - v**3/8 +
            v**2/8) ∂/∂v⊗dv
            sage: s == f*a
            True
        """
        from sage.manifolds.differentiable.mixed_form import MixedForm
        if isinstance(other, MixedForm):
            return other.parent()(self)._mul_(other)
        if not isinstance(other, TensorField):
            # Multiplication by a scalar field or a number
            return other * self
        # Tensor product:
        dom_resu = self._domain.intersection(other._domain)
        ambient_dom_resu = self._ambient_domain.intersection(other._ambient_domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the FreeModuleTensor version:
            return FreeModuleTensor.__mul__(self_r, other_r)
        dest_map = self._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu, subcodomain=ambient_dom_resu)
        vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)
        com_dom = []
        for dom in self_r._restrictions:
            if dom in other_r._restrictions:
                com_dom.append(dom)
        resu_rst = []
        for dom in com_dom:
            self_rr = self_r._restrictions[dom]
            other_rr = other_r._restrictions[dom]
            resu_rst.append(self_rr * other_rr)
        k1, l1 = self._tensor_type
        k2, l2 = other._tensor_type
        resu = vmodule.tensor((k1+k2, l1+l2),
                              sym=resu_rst[0]._sym,
                              antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst

        return resu

    def __truediv__(self, scalar) -> TensorField:
        r"""
        Division by a scalar field.

        INPUT:

        - ``scalar`` -- scalar field in the scalar field algebra over which
          the module containing ``self`` is defined

        OUTPUT: the tensor field ``scalar * self``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2)}, name='f')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: f.display()
            f: M → ℝ
            on U: (x, y) ↦ 1/(x^2 + y^2 + 1)
            on V: (u, v) ↦ 2/(u^2 + v^2 + 2)
            sage: s = a.__truediv__(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            (x^3 + x*y^2 + x) ∂/∂x⊗dx + (x^2 + y^2 + 1) ∂/∂x⊗dy
             + (y^3 + (x^2 + 1)*y) ∂/∂y⊗dx
            sage: f*s == a
            True

        Division by a number::

            sage: s = a.__truediv__(2); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            1/2*x ∂/∂x⊗dx + 1/2 ∂/∂x⊗dy + 1/2*y ∂/∂y⊗dx
            sage: s.display(e_uv)
            (1/4*u + 1/4) ∂/∂u⊗du + (1/4*u - 1/4) ∂/∂u⊗dv
             + (1/4*v + 1/4) ∂/∂v⊗du + (1/4*v - 1/4) ∂/∂v⊗dv
            sage: s == a/2
            True
            sage: 2*s == a
            True
        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = rst / scalar
        return resu

    def __call__(self, *args):
        r"""
        The tensor field acting on 1-forms and vector fields as a
        multilinear map.

        In the particular case of tensor field of type `(1,1)`, the action can
        be on a single vector field, the tensor field being identified to a
        field of tangent-space endomorphisms. The output is then a vector
        field.

        INPUT:

        - ``*args`` -- list of `k` 1-forms and `l` vector fields, ``self``
          being a tensor of type `(k,l)`

        OUTPUT:

        - either the scalar field resulting from the action of ``self`` on
          the 1-forms and vector fields passed as arguments or the vector
          field resulting from the action of ``self`` as a field of
          tangent-space endomorphisms (case of a type-(1,1) tensor field)

        TESTS:

        Action of a tensor field of type `(1,1)` on the 2-sphere::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e_xy,:] = [[x, 1], [y, 0]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: w = M.vector_field(name='w')
            sage: w[e_xy,:] = [y^2, x^2]
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: a = M.one_form(name='a')
            sage: a[e_xy,:] = [-1+y, x*y]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)

        The tensor field acting on a pair (1-form, vector field)::

            sage: s = t.__call__(a,w); s
            Scalar field t(a,w) on the 2-dimensional differentiable manifold M
            sage: s.display()
            t(a,w): M → ℝ
            on U: (x, y) ↦ x*y^4 + x*y^3 + x^2*y - x*y^2 - x^2
            on V: (u, v) ↦ 1/32*u^5 - 1/32*(3*u + 2)*v^4 + 1/32*v^5
             + 1/16*u^4 + 1/16*(u^2 + 2*u - 4)*v^3 + 1/16*(u^3 - 4)*v^2
             - 1/4*u^2 - 1/32*(3*u^4 + 4*u^3 - 8*u^2 + 16*u)*v
            sage: s.restrict(U) == t.restrict(U)(a.restrict(U), w.restrict(U))
            True
            sage: s.restrict(V) == t.restrict(V)(a.restrict(V), w.restrict(V))
            True

        The tensor field acting on vector field, as a field of tangent-space
        endomorphisms::

            sage: s = t.__call__(w); s
            Vector field t(w) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            t(w) = (x*y^2 + x^2) ∂/∂x + y^3 ∂/∂y
            sage: s.display(e_uv)
            t(w) = (1/4*u^3 + 1/4*(u + 1)*v^2 + 1/4*u^2 - 1/2*(u^2 - u)*v) ∂/∂u
             + (-1/4*(2*u - 1)*v^2 + 1/4*v^3 + 1/4*u^2 + 1/4*(u^2 + 2*u)*v) ∂/∂v
            sage: s.restrict(U) == t.restrict(U)(w.restrict(U))
            True
            sage: s.restrict(V) == t.restrict(V)(w.restrict(V))
            True
        """
        p = len(args)
        if p == 1 and self._tensor_type == (1,1):
            # type-(1,1) tensor acting as a field of tangent-space
            # endomorphisms:
            vector = args[0]
            if vector._tensor_type != (1,0):
                raise TypeError("the argument must be a vector field")
            dom_resu = self._domain.intersection(vector._domain)
            if dom_resu.is_manifestly_parallelizable():
                # call of the TensorFieldParal version:
                return self.restrict(dom_resu)(vector.restrict(dom_resu))
            if self._name is not None and vector._name is not None:
                name_resu = "{}({})".format(self._name, vector._name)
            else:
                name_resu = None
            if self._latex_name is not None and vector._latex_name is not None:
                latex_name_resu = r"{}\left({}\right)".format(self._latex_name,
                                                              vector._latex_name)
            else:
                latex_name_resu = None
            dest_map = vector._vmodule._dest_map
            dest_map_resu = dest_map.restrict(dom_resu)
            resu = dom_resu.vector_field(name=name_resu,
                                         latex_name=latex_name_resu,
                                         dest_map=dest_map_resu)
            for dom in self._common_subdomains(vector):
                if dom.is_subset(dom_resu):
                    resu._restrictions[dom] = \
                        self._restrictions[dom](vector._restrictions[dom])
            return resu
        # Generic case
        if p != self._tensor_rank:
            raise TypeError("{} arguments must be ".format(self._tensor_rank) +
                            "provided")
        # Domain of the result
        dom_resu = self._domain
        ambient_dom = self._ambient_domain
        for arg in args:
            dom_resu = dom_resu.intersection(arg._domain)
            ambient_dom = ambient_dom.intersection(arg._ambient_domain)
        self_r = self.restrict(dom_resu)
        args_r = [args[i].restrict(dom_resu) for i in range(p)]
        if ambient_dom.is_manifestly_parallelizable():
            # TensorFieldParal version
            return self_r(*args_r)
        else:
            resu = dom_resu.scalar_field()
            com_dom = []
            for dom in self_r._restrictions:
                for arg in args_r:
                    if dom not in arg._restrictions:
                        break
                else:
                    com_dom.append(dom)
            for dom in com_dom:
                self_rr = self_r._restrictions[dom]
                args_rr = [args_r[i]._restrictions[dom] for i in range(p)]
                resu_rr = self_rr(*args_rr)
                if resu_rr.is_trivial_zero():
                    for chart in resu_rr._domain._atlas:
                        resu._express[chart] = chart.zero_function()
                else:
                    for chart, expr in resu_rr._express.items():
                        resu._express[chart] = expr
            if resu.is_trivial_zero():
                return dom_resu._zero_scalar_field
            # Name of the output:
            res_name = None
            if self._name is not None:
                res_name = self._name + "("
                for i in range(p-1):
                    if args[i]._name is not None:
                        res_name += args[i]._name + ","
                    else:
                        res_name = None
                        break
                if res_name is not None:
                    if args[p-1]._name is not None:
                        res_name += args[p-1]._name + ")"
                    else:
                        res_name = None
            resu._name = res_name
            # LaTeX symbol of the output:
            res_latex = None
            if self._latex_name is not None:
                res_latex = self._latex_name + r"\left("
                for i in range(p-1):
                    if args[i]._latex_name is not None:
                        res_latex += args[i]._latex_name + ","
                    else:
                        res_latex = None
                        break
                if res_latex is not None:
                    if args[p-1]._latex_name is not None:
                        res_latex += args[p-1]._latex_name + r"\right)"
                    else:
                        res_latex = None
            resu._latex_name = res_latex
            return resu

    def trace(
        self,
        pos1=0,
        pos2=1,
        using: Optional[
            Union[PseudoRiemannianMetric, SymplecticForm, PoissonTensorField]
        ] = None,
    ):
        r"""
        Trace (contraction) on two slots of the tensor field.

        If a non-degenerate form is provided, the trace of a `(0,2)` tensor field
        is computed by first raising the last index.

        INPUT:

        - ``pos1`` -- (default: 0) position of the first index for the
          contraction, with the convention ``pos1=0`` for the first slot
        - ``pos2`` -- (default: 1) position of the second index for the
          contraction, with the same convention as for ``pos1``. The variance
          type of ``pos2`` must be opposite to that of ``pos1``
        - ``using`` -- (default: ``None``) a non-degenerate form

        OUTPUT:

        - tensor field resulting from the ``(pos1, pos2)`` contraction

        EXAMPLES:

        Trace of a type-`(1,1)` tensor field on a 2-dimensional
        non-parallelizable manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: W = U.intersection(V)
            sage: a = M.tensor_field(1,1, name='a')
            sage: a[e_xy,:] = [[1,x], [2,y]]
            sage: a.add_comp_by_continuation(e_uv, W, chart=c_uv)
            sage: s = a.trace() ; s
            Scalar field on the 2-dimensional differentiable manifold M
            sage: s.display()
            M → ℝ
            on U: (x, y) ↦ y + 1
            on V: (u, v) ↦ 1/2*u - 1/2*v + 1
            sage: s == a.trace(0,1) # explicit mention of the positions
            True

        The trace of a type-`(0,2)` tensor field using a metric::

            sage: g = M.metric('g')
            sage: g[0,0], g[0,1], g[1,1] = 1, 0, 1
            sage: g.trace(using=g).display()
            M → ℝ
            on U: (x, y) ↦ 2
            on W: (u, v) ↦ 2

        Instead of the explicit call to the method :meth:`trace`, one
        may use the index notation with Einstein convention (summation over
        repeated indices); it suffices to pass the indices as a string inside
        square brackets::

            sage: a['^i_i']
            Scalar field on the 2-dimensional differentiable manifold M
            sage: a['^i_i'] == s
            True

        Any letter can be used to denote the repeated index::

            sage: a['^b_b'] == s
            True

        Trace of a type-`(1,2)` tensor field::

            sage: b = M.tensor_field(1,2, name='b') ; b
            Tensor field b of type (1,2) on the 2-dimensional differentiable
             manifold M
            sage: b[e_xy,:] = [[[0,x+y], [y,0]], [[0,2], [3*x,-2]]]
            sage: b.add_comp_by_continuation(e_uv, W, chart=c_uv)  # long time
            sage: s = b.trace(0,1) ; s # contraction on first and second slots
            1-form on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            3*x dx + (x + y - 2) dy
            sage: s.display(e_uv)  # long time
            (5/4*u + 3/4*v - 1) du + (1/4*u + 3/4*v + 1) dv

        Use of the index notation::

            sage: b['^k_ki']
            1-form on the 2-dimensional differentiable manifold M
            sage: b['^k_ki'] == s  # long time
            True

        Indices not involved in the contraction may be replaced by dots::

            sage: b['^k_k.'] == s  # long time
            True

        The symbol ``^`` may be omitted::

            sage: b['k_k.'] == s  # long time
            True

        LaTeX notations are allowed::

            sage: b['^{k}_{ki}'] == s  # long time
            True

        Contraction on first and third slots::

            sage: s = b.trace(0,2) ; s
            1-form on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            2 dx + (y - 2) dy
            sage: s.display(e_uv)  # long time
            (1/4*u - 1/4*v) du + (-1/4*u + 1/4*v + 2) dv

        Use of index notation::

            sage: b['^k_.k'] == s  # long time
            True
        """
        if using is not None:
            if self.tensor_type() != (0, 2):
                raise ValueError(
                    "trace with respect to a non-degenerate form is only defined for type-(0,2) tensor fields"
                )
            return self.up(using, 1).trace()

        # The indices at pos1 and pos2 must be of different types:
        k_con = self._tensor_type[0]
        l_cov = self._tensor_type[1]
        if pos1 < k_con and pos2 < k_con:
            raise IndexError("contraction on two contravariant indices is " +
                             "not allowed")
        if pos1 >= k_con and pos2 >= k_con:
            raise IndexError("contraction on two covariant indices is " +
                             "not allowed")
        resu_rst = []
        for rst in self._restrictions.values():
            resu_rst.append(rst.trace(pos1, pos2))
        if (k_con, l_cov) == (1,1):
            # scalar field result
            resu = self._domain.scalar_field()
            all_zero = True
            for rst in resu_rst:
                if rst == 0:
                    for chart in rst._domain._atlas:
                        resu._express[chart] = 0
                else:
                    all_zero = False
                    for chart, funct in rst._express.items():
                        resu._express[chart] = funct
            if all_zero:
                resu = self._domain._zero_scalar_field
        else:
            # tensor field result
            resu = self._vmodule.tensor((k_con-1, l_cov-1),
                            sym=resu_rst[0]._sym, antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def contract(self, *args: Union[int, TensorField]) -> TensorField:
        r"""
        Contraction of ``self`` with another tensor field on one or
        more indices.

        INPUT:

        - ``pos1`` -- positions of the indices in the current tensor field
          involved in the contraction; ``pos1`` must be a sequence of integers,
          with 0 standing for the first index position, 1 for the second one,
          etc.; if ``pos1`` is not provided, a single contraction on the last
          index position of the tensor field is assumed
        - ``other`` -- the tensor field to contract with
        - ``pos2`` -- positions of the indices in ``other`` involved in the
          contraction, with the same conventions as for ``pos1``; if ``pos2``
          is not provided, a single contraction on the first index position of
          ``other`` is assumed

        OUTPUT:

        - tensor field resulting from the contraction at the positions
          ``pos1`` and ``pos2`` of the tensor field with ``other``

        EXAMPLES:

        Contractions of a type-`(1,1)` tensor field with a type-`(2,0)`
        one on a 2-dimensional non-parallelizable manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(1, 1, {eU: [[1, x], [0, 2]]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: b = M.tensor_field(2, 0, {eU: [[y, -1], [x+y, 2]]}, name='b')
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = a.contract(b) ; s   # contraction on last index of a and first one of b
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M

        Check 1: components with respect to the manifold's default
        frame (``eU``)::

            sage: all(bool(s[i,j] == sum(a[i,k]*b[k,j] for k in M.irange()))
            ....:     for i in M.irange() for j in M.irange())
            True

        Check 2: components with respect to the frame ``eV``::

            sage: all(bool(s[eV,i,j] == sum(a[eV,i,k]*b[eV,k,j]
            ....:                           for k in M.irange()))
            ....:      for i in M.irange() for j in M.irange())
            True

        Instead of the explicit call to the method :meth:`contract`, one
        may use the index notation with Einstein convention (summation over
        repeated indices); it suffices to pass the indices as a string inside
        square brackets::

            sage: a['^i_k']*b['^kj'] == s
            True

        Indices not involved in the contraction may be replaced by dots::

            sage: a['^._k']*b['^k.'] == s
            True

        LaTeX notation may be used::

            sage: a['^{i}_{k}']*b['^{kj}'] == s
            True

        Contraction on the last index of ``a`` and last index of ``b``::

            sage: s = a.contract(b, 1) ; s
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: a['^i_k']*b['^jk'] == s
            True

        Contraction on the first index of ``b`` and the last index of ``a``::

            sage: s = b.contract(0,a,1) ; s
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: b['^ki']*a['^j_k'] == s
            True

        The domain of the result is the intersection of the domains of
        the two tensor fields::

            sage: aU = a.restrict(U) ; bV = b.restrict(V)
            sage: s = aU.contract(b) ; s
            Tensor field of type (2,0) on the Open subset U of the
             2-dimensional differentiable manifold M
            sage: s = a.contract(bV) ; s
            Tensor field of type (2,0) on the Open subset V of the
             2-dimensional differentiable manifold M
            sage: s = aU.contract(bV) ; s
            Tensor field of type (2,0) on the Open subset W of the
             2-dimensional differentiable manifold M
            sage: s0 = a.contract(b)
            sage: s == s0.restrict(W)
            True

        The contraction can be performed on more than one index: ``c`` being a
        type-`(2,2)` tensor, contracting the indices in positions 2 and 3
        of ``c`` with respectively those in positions 0 and 1 of ``b`` is::

            sage: c = a*a ; c
            Tensor field of type (2,2) on the 2-dimensional differentiable
             manifold M
            sage: s = c.contract(2,3, b, 0,1) ; s  # long time
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M

        The same double contraction using index notation::

            sage: s == c['^.._kl']*b['^kl']  # long time
            True

        The symmetries are either conserved or destroyed by the contraction::

            sage: c = c.symmetrize(0,1).antisymmetrize(2,3)
            sage: c.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
            sage: s = b.contract(0, c, 2) ; s
            Tensor field of type (3,1) on the 2-dimensional differentiable
             manifold M
            sage: s.symmetries()
            symmetry: (1, 2);  no antisymmetry

        Case of a scalar field result::

            sage: a = M.one_form({eU: [y, 1+x]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: b = M.vector_field({eU: [x, y^2]}, name='b')
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a.display(eU)
            a = y dx + (x + 1) dy
            sage: b.display(eU)
            b = x ∂/∂x + y^2 ∂/∂y
            sage: s = a.contract(b) ; s
            Scalar field on the 2-dimensional differentiable manifold M
            sage: s.display()
            M → ℝ
            on U: (x, y) ↦ (x + 1)*y^2 + x*y
            on V: (u, v) ↦ 1/8*u^3 - 1/8*u*v^2 + 1/8*v^3 + 1/2*u^2 - 1/8*(u^2 + 4*u)*v
            sage: s == a['_i']*b['^i'] # use of index notation
            True
            sage: s == b.contract(a)
            True

        Case of a vanishing scalar field result::

            sage: b = M.vector_field({eU: [1+x, -y]}, name='b')
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = a.contract(b) ; s
            Scalar field zero on the 2-dimensional differentiable manifold M
            sage: s.display()
            zero: M → ℝ
            on U: (x, y) ↦ 0
            on V: (u, v) ↦ 0
        """
        nargs = len(args)
        for i, arg in enumerate(args):
            if isinstance(arg, TensorField):
                other = arg
                it = i
                break
        else:
            raise TypeError("a tensor field must be provided in the " +
                            "argument list")
        if it == 0:
            pos1 = (self._tensor_rank - 1,)
        else:
            pos1 = args[:it]
        if it == nargs-1:
            pos2 = (0,)
        else:
            pos2 = args[it+1:]
        ncontr = len(pos1) # number of contractions
        if len(pos2) != ncontr:
            raise IndexError("different number of indices for the contraction")
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise ValueError("incompatible ambient domains for contraction")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for contraction")
        dom_resu = self._domain.intersection(other._domain)
        ambient_dom_resu = self._ambient_domain.intersection(other._ambient_domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        k1, l1 = self._tensor_type
        k2, l2 = other._tensor_type
        tensor_type_resu = (k1 + k2 - ncontr, l1 + l2 - ncontr)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the FreeModuleTensor version:
            args = pos1 + (other_r,) + pos2
            return FreeModuleTensor.contract(self_r, *args)
        com_dom = []
        for dom in self_r._restrictions:
            if dom in other_r._restrictions:
                com_dom.append(dom)
        resu_rst = []
        for dom in com_dom:
            self_rr = self_r._restrictions[dom]
            other_rr = other_r._restrictions[dom]
            args = pos1 + (other_rr,) + pos2
            resu_rst.append(self_rr.contract(*args))
        if tensor_type_resu == (0,0):
            # scalar field result
            resu = dom_resu.scalar_field()
            all_zero = True
            for rst in resu_rst:
                if rst == 0:
                    for chart in rst._domain._atlas:
                        resu._express[chart] = 0
                else:
                    all_zero = False
                    for chart, funct in rst._express.items():
                        resu._express[chart] = funct
            if all_zero:
                resu = dom_resu._zero_scalar_field
        else:
            # tensor field result
            dest_map = self._vmodule._dest_map
            dest_map_resu = dest_map.restrict(dom_resu,
                                              subcodomain=ambient_dom_resu)
            vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)

            resu = vmodule.tensor(tensor_type_resu, sym=resu_rst[0]._sym,
                                  antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def symmetrize(self, *pos):
        r"""
        Symmetrization over some arguments.

        INPUT:

        - ``pos`` -- (default: ``None``) list of argument positions involved
          in the symmetrization (with the convention ``position=0`` for the
          first argument); if ``None``, the symmetrization is performed
          over all the arguments

        OUTPUT:

        - the symmetrized tensor field (instance of :class:`TensorField`)

        EXAMPLES:

        Symmetrization of a type-`(0,2)` tensor field on a 2-dimensional
        non-parallelizable manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(0,2, {eU: [[1,x], [2,y]]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a[eV,:]
            [ 1/4*u + 3/4 -1/4*u + 3/4]
            [ 1/4*v - 1/4 -1/4*v - 1/4]
            sage: s = a.symmetrize() ; s
            Field of symmetric bilinear forms on the 2-dimensional
             differentiable manifold M
            sage: s[eU,:]
            [        1 1/2*x + 1]
            [1/2*x + 1         y]
            sage: s[eV,:]
            [         1/4*u + 3/4 -1/8*u + 1/8*v + 1/4]
            [-1/8*u + 1/8*v + 1/4         -1/4*v - 1/4]
            sage: s == a.symmetrize(0,1)  # explicit positions
            True

        .. SEEALSO::

            For more details and examples, see
            :meth:`sage.tensor.modules.free_module_tensor.FreeModuleTensor.symmetrize`.
        """
        resu_rst = []
        for rst in self._restrictions.values():
            resu_rst.append(rst.symmetrize(*pos))
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_rst[0]._sym,
                                    antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def antisymmetrize(self, *pos):
        r"""
        Antisymmetrization over some arguments.

        INPUT:

        - ``pos`` -- (default: ``None``) list of argument positions involved
          in the antisymmetrization (with the convention ``position=0`` for
          the first argument); if ``None``, the antisymmetrization is
          performed over all the arguments

        OUTPUT:

        - the antisymmetrized tensor field (instance of :class:`TensorField`)

        EXAMPLES:

        Antisymmetrization of a type-`(0,2)` tensor field on a 2-dimensional
        non-parallelizable manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(0,2, {eU: [[1,x], [2,y]]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a[eV,:]
            [ 1/4*u + 3/4 -1/4*u + 3/4]
            [ 1/4*v - 1/4 -1/4*v - 1/4]
            sage: s = a.antisymmetrize() ; s
            2-form on the 2-dimensional differentiable manifold M
            sage: s[eU,:]
            [         0  1/2*x - 1]
            [-1/2*x + 1          0]
            sage: s[eV,:]
            [                   0 -1/8*u - 1/8*v + 1/2]
            [ 1/8*u + 1/8*v - 1/2                    0]
            sage: s == a.antisymmetrize(0,1)  # explicit positions
            True
            sage: s == a.antisymmetrize(1,0)  # the order of positions does not matter
            True

        .. SEEALSO::

            For more details and examples, see
            :meth:`sage.tensor.modules.free_module_tensor.FreeModuleTensor.antisymmetrize`.
        """
        resu_rst = []
        for rst in self._restrictions.values():
            resu_rst.append(rst.antisymmetrize(*pos))
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_rst[0]._sym,
                                    antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def lie_derivative(self, vector):
        r"""
        Lie derivative of ``self`` with respect to a vector field.

        INPUT:

        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken

        OUTPUT:

        - the tensor field that is the Lie derivative of the current tensor
          field with respect to ``vector``

        EXAMPLES:

        Lie derivative of a type-`(1,1)` tensor field along a vector field on
        a non-parallelizable 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, {e_xy: [[x, 1], [y, 0]]}, name='t')
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: w = M.vector_field({e_xy: [-y, x]}, name='w')
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: lt = t.lie_derivative(w); lt
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: lt.display(e_xy)
            ∂/∂x⊗dx - x ∂/∂x⊗dy + (-y - 1) ∂/∂y⊗dy
            sage: lt.display(e_uv)
            -1/2*u ∂/∂u⊗du + (1/2*u + 1) ∂/∂u⊗dv + (-1/2*v + 1) ∂/∂v⊗du + 1/2*v ∂/∂v⊗dv

        The result is cached::

            sage: t.lie_derivative(w) is lt
            True

        An alias is ``lie_der``::

            sage: t.lie_der(w) is t.lie_derivative(w)
            True

        Lie derivative of a vector field::

            sage: a = M.vector_field({e_xy: [1-x, x-y]}, name='a')
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: a.lie_der(w)
            Vector field on the 2-dimensional differentiable manifold M
            sage: a.lie_der(w).display(e_xy)
            x ∂/∂x + (-y - 1) ∂/∂y
            sage: a.lie_der(w).display(e_uv)
            (v - 1) ∂/∂u + (u + 1) ∂/∂v

        The Lie derivative is antisymmetric::

            sage: a.lie_der(w) == - w.lie_der(a)
            True

        and it coincides with the commutator of the two vector fields::

            sage: f = M.scalar_field({c_xy: 3*x-1, c_uv:  3/2*(u+v)-1})
            sage: a.lie_der(w)(f) == w(a(f)) - a(w(f))  # long time
            True
        """
        if vector._tensor_type != (1,0):
            raise TypeError("the argument must be a vector field")

        # The Lie derivative is cached in _lie_derivates while neither
        #    the tensor field nor ``vector`` have been modified
        if id(vector) not in self._lie_derivatives:
            # the computation must be performed:
            resu_rst = []
            for dom, rst in self._restrictions.items():
                resu_rst.append(rst.lie_derivative(vector.restrict(dom)))
            resu = self._vmodule.tensor(self._tensor_type,
                                        sym=resu_rst[0]._sym,
                                        antisym=resu_rst[0]._antisym)
            for rst in resu_rst:
                resu._restrictions[rst._domain] = rst
            self._lie_derivatives[id(vector)] = (vector, resu)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]

    lie_der = lie_derivative

    def at(self, point: ManifoldPoint) -> FreeModuleTensor:
        r"""
        Value of ``self`` at a point of its domain.

        If the current tensor field is

        .. MATH::

            t:\ U  \longrightarrow T^{(k,l)} M

        associated with the differentiable map

        .. MATH::

            \Phi:\ U \longrightarrow M,

        where `U` and `M` are two manifolds (possibly `U = M` and
        `\Phi = \mathrm{Id}_M`), then for any point `p \in U`, `t(p)`
        is a tensor on the tangent space to `M` at the point `\Phi(p)`.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` in the domain of the tensor field `U`

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          representing the tensor `t(p)` on the tangent vector space
          `T_{\Phi(p)} M`

        EXAMPLES:

        Tensor on a tangent space of a non-parallelizable 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(1, 1, {eU: [[1+y,x], [0,x+y]]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a.display(eU)
            a = (y + 1) ∂/∂x⊗dx + x ∂/∂x⊗dy + (x + y) ∂/∂y⊗dy
            sage: a.display(eV)
            a = (u + 1/2) ∂/∂u⊗du + (-1/2*u - 1/2*v + 1/2) ∂/∂u⊗dv
             + 1/2 ∂/∂v⊗du + (1/2*u - 1/2*v + 1/2) ∂/∂v⊗dv
            sage: p = M.point((2,3), chart=c_xy, name='p')
            sage: ap = a.at(p) ; ap
            Type-(1,1) tensor a on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: ap.parent()
            Free module of type-(1,1) tensors on the Tangent space at Point p
             on the 2-dimensional differentiable manifold M
            sage: ap.display(eU.at(p))
            a = 4 ∂/∂x⊗dx + 2 ∂/∂x⊗dy + 5 ∂/∂y⊗dy
            sage: ap.display(eV.at(p))
            a = 11/2 ∂/∂u⊗du - 3/2 ∂/∂u⊗dv + 1/2 ∂/∂v⊗du + 7/2 ∂/∂v⊗dv
            sage: p.coord(c_uv) # to check the above expression
            (5, -1)
        """
        if point not in self._domain:
            raise ValueError("the {} is not a point in the ".format(point) +
                             "domain of {}".format(self))
        for dom, rst in self._restrictions.items():
            if point in dom:
                return rst.at(point)

    def up(
        self,
        non_degenerate_form: Union[PseudoRiemannianMetric, SymplecticForm, PoissonTensorField],
        pos: Optional[int] = None,
    ) -> TensorField:
        r"""
        Compute a dual of the tensor field by raising some index with the
        given tensor field (usually, a pseudo-Riemannian metric, a symplectic form or a Poisson tensor).

        If `T` is the tensor field, `(k,l)` its type and `p` the position of a
        covariant index (i.e. `k\leq p < k+l`), this method called with
        ``pos`` `=p` yields the tensor field `T^\sharp` of type `(k+1,l-1)`
        whose components are

        .. MATH::

            (T^\sharp)^{a_1\ldots a_{k+1}}_{\phantom{a_1\ldots a_{k+1}}\,
            b_1 \ldots b_{l-1}} = g^{a_{k+1} i} \,
            T^{a_1\ldots a_k}_{\phantom{a_1\ldots a_k}\, b_1 \ldots b_{p-k}
            \, i \, b_{p-k+1}\ldots b_{l-1}},

        `g^{ab}` being the components of the inverse metric or the Poisson tensor, respectively.

        The reverse operation is :meth:`TensorField.down`.

        INPUT:

        - ``non_degenerate_form`` -- non-degenerate form `g`, or a Poisson tensor
        - ``pos`` -- (default: ``None``) position of the index (with the
          convention ``pos=0`` for the first index); if ``None``, the raising
          is performed over all the covariant indices, starting from the first
          one

        OUTPUT:

        - the tensor field `T^\sharp` resulting from the index raising
          operation

        EXAMPLES:

        Raising the index of a 1-form results in a vector field::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-y
            sage: w = M.one_form(-1, 2)
            sage: v = w.up(g) ; v
            Vector field on the 2-dimensional differentiable manifold M
            sage: v.display()
            ((2*x - 1)*y + 1)/(x^2*y^2 + (x + 1)*y - x - 1) ∂/∂x
             - (x*y + 2*x + 2)/(x^2*y^2 + (x + 1)*y - x - 1) ∂/∂y
            sage: ig = g.inverse(); ig[:]
            [ (y - 1)/(x^2*y^2 + (x + 1)*y - x - 1)      x*y/(x^2*y^2 + (x + 1)*y - x - 1)]
            [     x*y/(x^2*y^2 + (x + 1)*y - x - 1) -(x + 1)/(x^2*y^2 + (x + 1)*y - x - 1)]

        Using the index notation instead of :meth:`up`::

            sage: v == ig['^ab']*w['_b']
            True

        The reverse operation::

            sage: w1 = v.down(g) ; w1
            1-form on the 2-dimensional differentiable manifold M
            sage: w1.display()
            -dx + 2 dy
            sage: w1 == w
            True

        The reverse operation in index notation::

            sage: g['_ab']*v['^b'] == w
            True

        Raising the indices of a tensor field of type (0,2)::

            sage: t = M.tensor_field(0, 2, [[1,2], [3,4]])
            sage: tu0 = t.up(g, 0) ; tu0  # raising the first index
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: tu0[:]
            [  ((3*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1) 2*((2*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1)]
            [    (x*y - 3*x - 3)/(x^2*y^2 + (x + 1)*y - x - 1)   2*(x*y - 2*x - 2)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: tu0 == ig['^ac']*t['_cb'] # the same operation in index notation
            True
            sage: tuu0 = tu0.up(g) ; tuu0 # the two indices have been raised, starting from the first one
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: tuu0 == tu0['^a_c']*ig['^cb'] # the same operation in index notation
            True
            sage: tu1 = t.up(g, 1) ; tu1 # raising the second index
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: tu1 == ig['^ac']*t['_bc'] # the same operation in index notation
            True
            sage: tu1[:]
            [((2*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1) ((4*x + 3)*y - 3)/(x^2*y^2 + (x + 1)*y - x - 1)]
            [  (x*y - 2*x - 2)/(x^2*y^2 + (x + 1)*y - x - 1) (3*x*y - 4*x - 4)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: tuu1 = tu1.up(g) ; tuu1 # the two indices have been raised, starting from the second one
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: tuu1 == tu1['^a_c']*ig['^cb'] # the same operation in index notation
            True
            sage: tuu0 == tuu1 # the order of index raising is important
            False
            sage: tuu = t.up(g) ; tuu # both indices are raised, starting from the first one
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: tuu0 == tuu # the same order for index raising has been applied
            True
            sage: tuu1 == tuu # to get tuu1, indices have been raised from the last one, contrary to tuu
            False
            sage: d0tuu = tuu.down(g, 0) ; d0tuu # the first index is lowered again
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: dd0tuu = d0tuu.down(g) ; dd0tuu  # the second index is then lowered
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: d1tuu = tuu.down(g, 1) ; d1tuu # lowering operation, starting from the last index
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: dd1tuu = d1tuu.down(g) ; dd1tuu
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: ddtuu = tuu.down(g) ; ddtuu # both indices are lowered, starting from the last one
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: ddtuu == t # should be true
            True
            sage: dd0tuu == t # not true, because of the order of index lowering to get dd0tuu
            False
            sage: dd1tuu == t # should be true
            True
        """
        n_con = self._tensor_type[0] # number of contravariant indices = k
        if pos is None:
            result = self
            for p in range(n_con, self._tensor_rank):
                k = result._tensor_type[0]
                result = result.up(non_degenerate_form, k)
            return result

        if pos < n_con or pos > self._tensor_rank - 1:
            print("pos = {}".format(pos))
            raise ValueError("position out of range")

        from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
        from sage.manifolds.differentiable.symplectic_form import SymplecticForm
        from sage.manifolds.differentiable.poisson_tensor import PoissonTensorField

        if isinstance(non_degenerate_form, PseudoRiemannianMetric):
            return self.contract(pos, non_degenerate_form.inverse(), 1)
        elif isinstance(non_degenerate_form, SymplecticForm):
            return self.contract(pos, non_degenerate_form.poisson(), 1)
        elif isinstance(non_degenerate_form, PoissonTensorField):
            return self.contract(pos, non_degenerate_form, 1)
        else:
            raise ValueError("The non-degenerate form has to be a metric, a symplectic form or a Poisson tensor field")

    def down(
        self,
        non_degenerate_form: Union[PseudoRiemannianMetric, SymplecticForm],
        pos: Optional[int] = None,
    ) -> TensorField:
        r"""
        Compute a dual of the tensor field by lowering some index with a
        given non-degenerate form (pseudo-Riemannian metric or symplectic form).

        If `T` is the tensor field, `(k,l)` its type and `p` the position of a
        contravariant index (i.e. `0\leq p < k`), this method called with
        ``pos`` `=p`  yields the tensor field `T^\flat` of type `(k-1,l+1)`
        whose components are

        .. MATH::

            (T^\flat)^{a_1\ldots a_{k-1}}_{\phantom{a_1\ldots a_{k-1}}
            \, b_1 \ldots b_{l+1}} = g_{i b_1} \,
            T^{a_1\ldots a_{p} \, i \, a_{p+1}\ldots a_{k-1}}_{\phantom{a_1
            \ldots a_{p} \, i \, a_{p+1}\ldots a_{k-1}}\, b_2 \ldots b_{l+1}},

        `g_{ab}` being the components of the metric tensor or the symplectic form, respectively.

        The reverse operation is :meth:`TensorField.up`.

        INPUT:

        - ``non_degenerate_form`` -- non-degenerate form `g`
        - ``pos`` -- (default: ``None``) position of the index (with the
          convention ``pos=0`` for the first index); if ``None``, the lowering
          is performed over all the contravariant indices, starting from the
          last one

        OUTPUT:

        - the tensor field `T^\flat` resulting from the index lowering
          operation

        EXAMPLES:

        Lowering the index of a vector field results in a 1-form::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-y
            sage: v = M.vector_field(-1, 2)
            sage: w = v.down(g) ; w
            1-form on the 2-dimensional differentiable manifold M
            sage: w.display()
            (2*x*y - x - 1) dx + (-(x + 2)*y + 2) dy

        Using the index notation instead of :meth:`down`::

            sage: w == g['_ab']*v['^b']
            True

        The reverse operation::

            sage: v1 = w.up(g) ; v1
            Vector field on the 2-dimensional differentiable manifold M
            sage: v1 == v
            True

        Lowering the indices of a tensor field of type (2,0)::

            sage: t = M.tensor_field(2, 0, [[1,2], [3,4]])
            sage: td0 = t.down(g, 0) ; td0  # lowering the first index
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: td0 == g['_ac']*t['^cb'] # the same operation in index notation
            True
            sage: td0[:]
            [  3*x*y + x + 1   (x - 3)*y + 3]
            [4*x*y + 2*x + 2 2*(x - 2)*y + 4]
            sage: tdd0 = td0.down(g) ; tdd0 # the two indices have been lowered, starting from the first one
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: tdd0 == g['_ac']*td0['^c_b'] # the same operation in index notation
            True
            sage: tdd0[:]
            [      4*x^2*y^2 + x^2 + 5*(x^2 + x)*y + 2*x + 1 2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3]
            [(3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: td1 = t.down(g, 1) ; td1  # lowering the second index
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: td1 == g['_ac']*t['^bc'] # the same operation in index notation
            True
            sage: td1[:]
            [  2*x*y + x + 1   (x - 2)*y + 2]
            [4*x*y + 3*x + 3 (3*x - 4)*y + 4]
            sage: tdd1 = td1.down(g) ; tdd1 # the two indices have been lowered, starting from the second one
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: tdd1 == g['_ac']*td1['^c_b'] # the same operation in index notation
            True
            sage: tdd1[:]
            [      4*x^2*y^2 + x^2 + 5*(x^2 + x)*y + 2*x + 1 (3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2]
            [2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: tdd1 == tdd0   # the order of index lowering is important
            False
            sage: tdd = t.down(g) ; tdd  # both indices are lowered, starting from the last one
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: tdd[:]
            [      4*x^2*y^2 + x^2 + 5*(x^2 + x)*y + 2*x + 1 (3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2]
            [2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: tdd0 == tdd  # to get tdd0, indices have been lowered from the first one, contrary to tdd
            False
            sage: tdd1 == tdd  # the same order for index lowering has been applied
            True
            sage: u0tdd = tdd.up(g, 0) ; u0tdd # the first index is raised again
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: uu0tdd = u0tdd.up(g) ; uu0tdd # the second index is then raised
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: u1tdd = tdd.up(g, 1) ; u1tdd  # raising operation, starting from the last index
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: uu1tdd = u1tdd.up(g) ; uu1tdd
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: uutdd = tdd.up(g) ; uutdd  # both indices are raised, starting from the first one
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: uutdd == t  # should be true
            True
            sage: uu0tdd == t # should be true
            True
            sage: uu1tdd == t # not true, because of the order of index raising to get uu1tdd
            False
        """
        n_con = self._tensor_type[0]  # number of contravariant indices = k
        if pos is None:
            result = self
            for p in range(n_con):
                k = result._tensor_type[0]
                result = result.down(non_degenerate_form, k - 1)
            return result
        if not isinstance(pos, (int, Integer)):
            raise TypeError("the argument 'pos' must be an integer")
        if pos < 0 or pos >= n_con:
            print("pos = {}".format(pos))
            raise ValueError("position out of range")
        return non_degenerate_form.contract(0, self, pos)

    def divergence(self, metric=None):
        r"""
        Return the divergence of ``self`` (with respect to a given
        metric).

        The divergence is taken on the *last* index: if
        ``self`` is a tensor field `t` of type `(k,0)` with `k\geq 1`, the
        divergence of `t` with respect to the metric `g` is the tensor field
        of type `(k-1,0)` defined by

        .. MATH::

            (\mathrm{div}\, t)^{a_1\ldots a_{k-1}} =
            \nabla_i t^{a_1\ldots a_{k-1} i} =
            (\nabla t)^{a_1\ldots a_{k-1} i}_{\phantom{a_1\ldots a_{k-1} i}\, i},

        where `\nabla` is the Levi-Civita connection of `g` (cf.
        :class:`~sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection`).

        This definition is extended to tensor fields of type `(k,l)` with
        `k\geq 0` and `l\geq 1`, by raising the last index with the metric `g`:
        `\mathrm{div}\, t` is then the tensor field of type `(k,l-1)` defined by

        .. MATH::

            (\mathrm{div}\, t)^{a_1\ldots a_k}_{\phantom{a_1\ldots a_k}\, b_1
            \ldots b_{l-1}} = \nabla_i (g^{ij} t^{a_1\ldots a_k}_{\phantom{a_1
            \ldots a_k}\, b_1\ldots b_{l-1} j})
            = (\nabla t^\sharp)^{a_1\ldots a_k i}_{\phantom{a_1\ldots a_k i}\,
            b_1\ldots b_{l-1} i},

        where `t^\sharp` is the tensor field deduced from `t` by raising the
        last index with the metric `g` (see :meth:`up`).

        INPUT:

        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the divergence; if none is provided,
          the domain of ``self`` is supposed to be endowed with a default
          metric (i.e. is supposed to be pseudo-Riemannian manifold, see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the divergence.

        OUTPUT:

        - instance of either
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`
          if `(k,l)=(1,0)` (``self`` is a vector field) or `(k,l)=(0,1)`
          (``self`` is a 1-form) or of :class:`TensorField` if `k+l\geq 2`
          representing the divergence of ``self`` with respect to ``metric``

        EXAMPLES:

        Divergence of a vector field in the Euclidean plane::

            sage: M.<x,y> = EuclideanSpace()
            sage: v = M.vector_field(x, y, name='v')
            sage: s = v.divergence(); s
            Scalar field div(v) on the Euclidean plane E^2
            sage: s.display()
            div(v): E^2 → ℝ
               (x, y) ↦ 2

        A shortcut alias of ``divergence`` is ``div``::

            sage: v.div() == s
            True

        The function :func:`~sage.manifolds.operators.div` from the
        :mod:`~sage.manifolds.operators` module can be used instead of the
        method :meth:`divergence`::

            sage: from sage.manifolds.operators import div
            sage: div(v) == s
            True

        The divergence can be taken with respect to a metric tensor that is
        not the default one::

            sage: h = M.lorentzian_metric('h')
            sage: h[1,1], h[2,2] = -1, 1/(1+x^2+y^2)
            sage: s = v.div(h); s
            Scalar field div_h(v) on the Euclidean plane E^2
            sage: s.display()
            div_h(v): E^2 → ℝ
               (x, y) ↦ (x^2 + y^2 + 2)/(x^2 + y^2 + 1)

        The standard formula

        .. MATH::

            \mathrm{div}_h \, v = \frac{1}{\sqrt{|\det h|}}
            \frac{\partial}{\partial x^i} \left( \sqrt{|\det h|} \, v^i \right)

        is checked as follows::

            sage: sqrth = h.sqrt_abs_det().expr(); sqrth
            1/sqrt(x^2 + y^2 + 1)
            sage: s == 1/sqrth * sum( (sqrth*v[i]).diff(i) for i in M.irange())
            True

        A divergence-free vector::

            sage: w = M.vector_field(-y, x, name='w')
            sage: w.div().display()
            div(w): E^2 → ℝ
               (x, y) ↦ 0
            sage: w.div(h).display()
            div_h(w): E^2 → ℝ
               (x, y) ↦ 0

        Divergence of a type-``(2,0)`` tensor field::

            sage: t = v*w; t
            Tensor field v⊗w of type (2,0) on the Euclidean plane E^2
            sage: s = t.div(); s
            Vector field div(v⊗w) on the Euclidean plane E^2
            sage: s.display()
            div(v⊗w) = -y e_x + x e_y
        """
        n_con = self._tensor_type[0] # number of contravariant indices = k
        n_cov = self._tensor_type[1] # number of covariant indices = l
        default_metric = metric is None
        if default_metric:
            metric = self._domain.metric()
        nabla = metric.connection()
        if n_cov == 0:
            resu = nabla(self).trace(n_con-1, n_con)
        else:
            tup = self.up(metric, self._tensor_rank-1)
            resu = nabla(tup).trace(n_con, self._tensor_rank)
        if self._name is not None:
            if default_metric:
                resu._name = "div({})".format(self._name)
                resu._latex_name = r"\mathrm{div}\left(" + self._latex_name + \
                                   r"\right)"
            else:
                resu._name = "div_{}({})".format(metric._name, self._name)
                resu._latex_name = r"\mathrm{div}_{" + metric._latex_name + \
                                   r"}\left(" + self._latex_name + r"\right)"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    div = divergence

    def laplacian(self, metric=None):
        r"""
        Return the Laplacian of ``self`` with respect to a given
        metric (Laplace-Beltrami operator).

        If ``self`` is a tensor field `t` of type `(k,l)`, the Laplacian of `t`
        with respect to the metric `g` is the tensor field of type `(k,l)`
        defined by

        .. MATH::

            (\Delta t)^{a_1\ldots a_k}_{\phantom{a_1\ldots a_k}\,{b_1\ldots b_k}}
            = \nabla_i \nabla^i
            t^{a_1\ldots a_k}_{\phantom{a_1\ldots a_k}\,{b_1\ldots b_k}},

        where `\nabla` is the Levi-Civita connection of `g` (cf.
        :class:`~sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection`)
        and `\nabla^i := g^{ij} \nabla_j`. The operator
        `\Delta = \nabla_i \nabla^i` is called the *Laplace-Beltrami operator*
        of metric `g`.

        INPUT:

        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the Laplacian; if none is provided, the
          domain of ``self`` is supposed to be endowed with a default metric
          (i.e. is supposed to be pseudo-Riemannian manifold, see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the Laplacian

        OUTPUT:

        - instance of :class:`TensorField` representing the Laplacian of
          ``self``

        EXAMPLES:

        Laplacian of a vector field in the Euclidean plane::

            sage: M.<x,y> = EuclideanSpace()
            sage: v = M.vector_field(x^3 + y^2, x*y, name='v')
            sage: Dv = v.laplacian(); Dv
            Vector field Delta(v) on the Euclidean plane E^2
            sage: Dv.display()
            Delta(v) = (6*x + 2) e_x

        The function :func:`~sage.manifolds.operators.laplacian` from the
        :mod:`~sage.manifolds.operators` module can be used instead of the
        method :meth:`laplacian`::

            sage: from sage.manifolds.operators import laplacian
            sage: laplacian(v) == Dv
            True

        In the present case (Euclidean metric and Cartesian coordinates), the
        components of the Laplacian are the Laplacians of the components::

            sage: all(Dv[[i]] == laplacian(v[[i]]) for i in M.irange())
            True

        The Laplacian can be taken with respect to a metric tensor that is
        not the default one::

            sage: h = M.lorentzian_metric('h')
            sage: h[1,1], h[2,2] = -1, 1+x^2
            sage: Dv = v.laplacian(h); Dv
            Vector field Delta_h(v) on the Euclidean plane E^2
            sage: Dv.display()
            Delta_h(v) = -(8*x^5 - 2*x^4 - x^2*y^2 + 15*x^3 - 4*x^2 + 6*x
             - 2)/(x^4 + 2*x^2 + 1) e_x - 3*x^3*y/(x^4 + 2*x^2 + 1) e_y
        """
        n_con = self._tensor_type[0] # number of contravariant indices = k
        trank = self._tensor_rank    # k + l
        default_metric = metric is None
        if default_metric:
            metric = self._domain.metric()
        nabla = metric.connection()
        tmp = nabla(nabla(self).up(metric, pos=trank))
        resu = tmp.trace(n_con, trank+1)
        if self._name is not None:
            if default_metric:
                resu._name = "Delta({})".format(self._name)
                resu._latex_name = r"\Delta\left(" + self._latex_name + \
                                   r"\right)"
            else:
                resu._name = "Delta_{}({})".format(metric._name, self._name)
                resu._latex_name = r"\Delta_{" + metric._latex_name + \
                                   r"}\left(" + self._latex_name + r"\right)"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    def dalembertian(self, metric=None):
        r"""
        Return the d'Alembertian of ``self`` with respect to a given
        Lorentzian metric.

        The *d'Alembertian* of a tensor field `t` with respect to a Lorentzian
        metric `g` is nothing but the Laplace-Beltrami operator of `g` applied
        to `t` (see :meth:`laplacian`); if ``self`` a tensor field `t` of type
        `(k,l)`, the d'Alembertian of `t` with respect to `g` is then the
        tensor field of type `(k,l)` defined by

        .. MATH::

            (\Box t)^{a_1\ldots a_k}_{\phantom{a_1\ldots a_k}\,{b_1\ldots b_k}}
            = \nabla_i \nabla^i
            t^{a_1\ldots a_k}_{\phantom{a_1\ldots a_k}\,{b_1\ldots b_k}},

        where `\nabla` is the Levi-Civita connection of `g` (cf.
        :class:`~sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection`)
        and `\nabla^i := g^{ij} \nabla_j`.

        .. NOTE::

            If the metric `g` is not Lorentzian, the name *d'Alembertian* is
            not appropriate and one should use :meth:`laplacian` instead.

        INPUT:

        - ``metric`` -- (default: ``None``) the Lorentzian metric `g`
          involved in the definition of the d'Alembertian; if none is provided,
          the domain of ``self`` is supposed to be endowed with a default
          Lorentzian metric (i.e. is supposed to be Lorentzian manifold, see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the d'Alembertian

        OUTPUT:

        - instance of :class:`TensorField` representing the d'Alembertian of
          ``self``

        EXAMPLES:

        d'Alembertian of a vector field in Minkowski spacetime, representing
        the electric field of a simple plane electromagnetic wave::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: e = M.vector_field(name='e')
            sage: e[1] = cos(t-z)
            sage: e.display()  # plane wave propagating in the z direction
            e = cos(t - z) ∂/∂x
            sage: De = e.dalembertian(); De # long time
            Vector field Box(e) on the 4-dimensional Lorentzian manifold M

        The function :func:`~sage.manifolds.operators.dalembertian` from the
        :mod:`~sage.manifolds.operators` module can be used instead of the
        method :meth:`dalembertian`::

            sage: from sage.manifolds.operators import dalembertian
            sage: dalembertian(e) == De # long time
            True

        We check that the electric field obeys the wave equation::

            sage: De.display() # long time
            Box(e) = 0
        """
        default_metric = metric is None
        if default_metric:
            metric = self._domain.metric()
        nm2 = self._domain.dim() - 2
        if metric.signature() not in [nm2, -nm2]:
            raise TypeError("the {} is not a Lorentzian ".format(metric) +
                            "metric; use laplacian() instead")
        resu = self.laplacian(metric=metric)
        if self._name is not None:
            if default_metric:
                resu._name = "Box({})".format(self._name)
                resu._latex_name = r"\Box\left(" + self._latex_name + \
                                   r"\right)"
            else:
                resu._name = "Box_{}({})".format(metric._name, self._name)
                resu._latex_name = r"\Box_{" + metric._latex_name + \
                                   r"}\left(" + self._latex_name + r"\right)"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    def along(self, mapping):
        r"""
        Return the tensor field deduced from ``self`` via a differentiable map,
        the codomain of which is included in the domain of ``self``.

        More precisely, if ``self`` is a tensor field `t` on `M` and if
        `\Phi: U \rightarrow M` is a differentiable map from some
        differentiable manifold `U` to `M`, the returned object is
        a tensor field `\tilde t` along `U` with values on `M` such that

        .. MATH::

           \forall p \in U,\  \tilde t(p) = t(\Phi(p)).

        INPUT:

        - ``mapping`` -- differentiable map `\Phi: U \rightarrow M`

        OUTPUT: tensor field `\tilde t` along `U` defined above

        EXAMPLES:

        Let us consider the 2-dimensional sphere `S^2`::

            sage: M = Manifold(2, 'S^2') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)

        and the following map from the open interval `(0,5\pi/2)` to `S^2`,
        the image of it being the great circle `x=0`, `u=0`, which goes through
        the North and South poles::

            sage: I.<t> = manifolds.OpenInterval(0, 5*pi/2)
            sage: J = I.open_interval(0, 3*pi/2)
            sage: K = I.open_interval(pi, 5*pi/2)
            sage: c_J = J.canonical_chart(); c_K = K.canonical_chart()
            sage: Phi = I.diff_map(M, {(c_J, c_xy):
            ....:                      (0, sgn(pi-t)*sqrt((1+cos(t))/(1-cos(t)))),
            ....:                      (c_K, c_uv):
            ....:                      (0,  sgn(t-2*pi)*sqrt((1-cos(t))/(1+cos(t))))},
            ....:                  name='Phi')

        Let us consider a vector field on `S^2`::

            sage: eU = c_xy.frame(); eV = c_uv.frame()
            sage: w = M.vector_field(name='w')
            sage: w[eU,0] = 1
            sage: w.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: w.display(eU)
            w = ∂/∂x
            sage: w.display(eV)
            w = (-u^2 + v^2) ∂/∂u - 2*u*v ∂/∂v

        We have then::

            sage: wa = w.along(Phi); wa
            Vector field w along the Real interval (0, 5/2*pi) with values on
             the 2-dimensional differentiable manifold S^2
            sage: wa.display(eU.along(Phi))
            w = ∂/∂x
            sage: wa.display(eV.along(Phi))
            w = -(cos(t) - 1)*sgn(-2*pi + t)^2/(cos(t) + 1) ∂/∂u

        Some tests::

            sage: p = K.an_element()
            sage: wa.at(p) == w.at(Phi(p))
            True
            sage: wa.at(J(4*pi/3)) == wa.at(K(4*pi/3))
            True
            sage: wa.at(I(4*pi/3)) == wa.at(K(4*pi/3))
            True
            sage: wa.at(K(7*pi/4)) == eU[0].at(Phi(I(7*pi/4))) # since eU[0]=∂/∂x
            True
        """
        dom = self._domain
        if self._ambient_domain is not dom:
            raise ValueError("{} is not a tensor field ".format(self) +
                             "with values in the {}".format(dom))
        if mapping.codomain().is_subset(dom):
            rmapping = mapping
        else:
            rmapping = None
            for doms, rest in mapping._restrictions.items():
                if doms[1].is_subset(dom):
                    rmapping = rest
                    break
            else:
                raise ValueError("the codomain of {} is not ".format(mapping) +
                                 "included in the domain of {}".format(self))
        resu_ambient_domain = rmapping.codomain()
        if resu_ambient_domain.is_manifestly_parallelizable():
            return self.restrict(resu_ambient_domain).along(rmapping)
        dom_resu = rmapping.domain()
        vmodule = dom_resu.vector_field_module(dest_map=rmapping)
        resu = vmodule.tensor(self._tensor_type, name=self._name,
                              latex_name=self._latex_name, sym=self._sym,
                              antisym=self._antisym)
        for rdom in resu_ambient_domain._parallelizable_parts:
            if rdom in resu_ambient_domain._top_subsets:
                for chart1, chart2 in rmapping._coord_expression:
                    if chart2.domain() is rdom:
                        dom1 = chart1.domain()
                        rrmap = rmapping.restrict(dom1, subcodomain=rdom)
                        resu._restrictions[dom1] = self.restrict(rdom).along(rrmap)
        return resu

    def set_calc_order(self, symbol, order, truncate=False):
        r"""
        Trigger a series expansion with respect to a small parameter in
        computations involving the tensor field.

        This property is propagated by usual operations. The internal
        representation must be ``SR`` for this to take effect.

        If the small parameter is `\epsilon` and `T` is ``self``, the
        power series expansion to order `n` is

        .. MATH::

            T = T_0 + \epsilon T_1 + \epsilon^2 T_2 + \cdots + \epsilon^n T_n
                + O(\epsilon^{n+1}),

        where `T_0, T_1, \ldots, T_n` are `n+1` tensor fields of the same
        tensor type as ``self`` and do not depend upon `\epsilon`.

        INPUT:

        - ``symbol`` -- symbolic variable (the "small parameter" `\epsilon`)
          with respect to which the components of ``self`` are expanded in
          power series
        - ``order`` -- integer; the order `n` of the expansion, defined as the
          degree of the polynomial representing the truncated power series in
          ``symbol``
        - ``truncate`` -- boolean (default: ``False``); determines whether the
          components of ``self`` are replaced by their expansions to the
          given order

        EXAMPLES:

        Let us consider two vector fields depending on a small parameter `h`
        on a non-parallelizable manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.vector_field()
            sage: h = var('h', domain='real')
            sage: a[eU,:] = (cos(h*x), -y)
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: b = M.vector_field()
            sage: b[eU,:] = (exp(h*x), exp(h*y))
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)

        If we set the calculus order on one of the vector fields, any operation
        involving both of them is performed to that order::

            sage: a.set_calc_order(h, 2)
            sage: s = a + b
            sage: s[eU,:]
            [h*x + 2, 1/2*h^2*y^2 + h*y - y + 1]
            sage: s[eV,:]
            [1/8*(u^2 - 2*u*v + v^2)*h^2 + h*u - 1/2*u + 1/2*v + 3,
             -1/8*(u^2 - 2*u*v + v^2)*h^2 + h*v + 1/2*u - 1/2*v + 1]

        Note that the components of ``a`` have not been affected by the above
        call to ``set_calc_order``::

            sage: a[eU,:]
            [cos(h*x), -y]
            sage: a[eV,:]
            [cos(1/2*h*u)*cos(1/2*h*v) - sin(1/2*h*u)*sin(1/2*h*v) - 1/2*u + 1/2*v,
             cos(1/2*h*u)*cos(1/2*h*v) - sin(1/2*h*u)*sin(1/2*h*v) + 1/2*u - 1/2*v]

        To have ``set_calc_order`` act on them, set the optional argument
        ``truncate`` to ``True``::

            sage: a.set_calc_order(h, 2, truncate=True)
            sage: a[eU,:]
            [-1/2*h^2*x^2 + 1, -y]
            sage: a[eV,:]
            [-1/8*(u^2 + 2*u*v + v^2)*h^2 - 1/2*u + 1/2*v + 1,
             -1/8*(u^2 + 2*u*v + v^2)*h^2 + 1/2*u - 1/2*v + 1]
        """
        for rst in self._restrictions.values():
            rst.set_calc_order(symbol, order, truncate)
        self._del_derived()

    def apply_map(self, fun, frame=None, chart=None,
                  keep_other_components=False):
        r"""
        Apply a function to the coordinate expressions of all components of
        ``self`` in a given vector frame.

        This method allows operations like factorization, expansion,
        simplification or substitution to be performed on all components of
        ``self`` in a given vector frame (see examples below).

        INPUT:

        - ``fun`` -- function to be applied to the coordinate expressions of
          the components
        - ``frame`` -- (default: ``None``) vector frame defining the
          components on which the operation ``fun`` is to be performed; if
          ``None``, the default frame of the domain of ``self`` is assumed
        - ``chart`` -- (default: ``None``) coordinate chart; if specified, the
          operation ``fun`` is performed only on the coordinate expressions
          with respect to ``chart`` of the components w.r.t. ``frame``; if
          ``None``, the operation ``fun`` is performed on all available
          coordinate expressions
        - ``keep_other_components`` -- boolean (default: ``False``); determine whether
          the components with respect to vector frames distinct from ``frame``
          and having the same domain as ``frame`` are kept. If ``fun`` is
          non-destructive, ``keep_other_components`` can be set to ``True``;
          otherwise, it is advised to set it to ``False`` (the default) in
          order to avoid any inconsistency between the various sets of
          components

        EXAMPLES:

        Factorizing all components in the default frame of a vector field::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a, b = var('a b')
            sage: v = M.vector_field(x^2 - y^2, a*(b^2 - b)*x)
            sage: v.display()
            (x^2 - y^2) ∂/∂x + (b^2 - b)*a*x ∂/∂y
            sage: v.apply_map(factor)
            sage: v.display()
            (x + y)*(x - y) ∂/∂x + a*(b - 1)*b*x ∂/∂y

        Performing a substitution in all components in the default frame::

            sage: v.apply_map(lambda f: f.subs({a: 2}))
            sage: v.display()
            (x + y)*(x - y) ∂/∂x + 2*(b - 1)*b*x ∂/∂y

        Specifying the vector frame via the argument ``frame``::

            sage: P.<p, q> = M.chart()
            sage: X_to_P = X.transition_map(P, [x + 1, y - 1])
            sage: P_to_X = X_to_P.inverse()
            sage: v.display(P)
            (p^2 - q^2 - 2*p - 2*q) ∂/∂p + (-2*b^2 + 2*(b^2 - b)*p + 2*b) ∂/∂q
            sage: v.apply_map(lambda f: f.subs({b: pi}), frame=P.frame())
            sage: v.display(P)
            (p^2 - q^2 - 2*p - 2*q) ∂/∂p + (2*pi - 2*pi^2 - 2*(pi - pi^2)*p) ∂/∂q

        Note that the required operation has been performed in all charts::

            sage: v.display(P.frame(), P)
            (p^2 - q^2 - 2*p - 2*q) ∂/∂p + (2*pi - 2*pi^2 - 2*(pi - pi^2)*p) ∂/∂q
            sage: v.display(P.frame(), X)
            (x + y)*(x - y) ∂/∂p + 2*pi*(pi - 1)*x ∂/∂q

        By default, the components of ``v`` in frames distinct from the
        specified one have been deleted::

            sage: X.frame() in v._components
            False

        When requested, they are recomputed by change-of-frame formulas,
        thereby enforcing the consistency between the representations in
        various vector frames. In particular, we can check that the
        substitution of ``b`` by ``pi``, which was asked in ``P.frame()``,
        is effective in ``X.frame()`` as well::

            sage: v.display(X.frame(), X)
            (x + y)*(x - y) ∂/∂x + 2*pi*(pi - 1)*x ∂/∂y

        When the requested operation does not change the value of the tensor
        field, one can use the keyword argument ``keep_other_components=True``,
        in order to avoid the recomputation of the components in other frames::

            sage: v.apply_map(factor, keep_other_components=True)
            sage: v.display()
            (x + y)*(x - y) ∂/∂x + 2*pi*(pi - 1)*x ∂/∂y

        The components with respect to ``P.frame()`` have been kept::

            sage: P.frame() in v._components
            True

        One can restrict the operation to expressions in a given chart, via
        the argument ``chart``::

            sage: v.display(X.frame(), P)
            (p + q)*(p - q - 2) ∂/∂x + 2*pi*(pi - 1)*(p - 1) ∂/∂y
            sage: v.apply_map(expand, chart=P)
            sage: v.display(X.frame(), P)
            (p^2 - q^2 - 2*p - 2*q) ∂/∂x + (2*pi + 2*pi^2*p - 2*pi^2 - 2*pi*p) ∂/∂y
            sage: v.display(X.frame(), X)
            (x + y)*(x - y) ∂/∂x + 2*pi*(pi - 1)*x ∂/∂y

        TESTS:

        Check that the cached quantities derived from the components are
        erased::

            sage: w = M.vector_field(a*x, 0)
            sage: diff(w[[0]]).display()
            a dx
            sage: w.apply_map(lambda t: t.subs(a=-2))
            sage: w.display()
            -2*x ∂/∂x
            sage: diff(w[[0]]).display()
            -2 dx
        """
        # The dictionary of components w.r.t. frame:
        if keep_other_components:
            comps = self.comp(frame)._comp
        else:
            comps = self.set_comp(frame)._comp # set_comp() deletes the
                                               # components in other frames
        if chart:
            for scalar in comps.values():
                scalar.add_expr(fun(scalar.expr(chart=chart)), chart=chart)
        else:
            for scalar in comps.values():
                cfunc_dict = {}  # new dict of chart functions in order not to
                                 # modify scalar._express while looping on it
                for ch, fct in scalar._express.items():
                    cfunc_dict[ch] = ch.function(fun(fct.expr()))
                scalar._express = cfunc_dict
                scalar._del_derived()
