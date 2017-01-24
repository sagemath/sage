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

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [KN1963]_
- [Lee2013]_
- [ONe1983]_

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from __future__ import print_function

from sage.rings.integer import Integer
from sage.structure.element import ModuleElement
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.tensor_with_indices import TensorWithIndices
from sage.rings.integer_ring import ZZ

class TensorField(ModuleElement):
    r"""
    Tensor field along a differentiable manifold.

    An instance of this class is a tensor field along a differentiable
    manifold `U` with values on a differentiable manifold `M`, via a
    differentiable map `\Phi: U \rightarrow M`. More precisely, given two
    non-negative integers `k` and `l` and a differentiable map

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

        t(p):\ \underbrace{T_q^*M\times\cdots\times T_q^*M}_{k\ \; \mbox{times}}
        \times \underbrace{T_q M\times\cdots\times T_q M}_{l\ \; \mbox{times}}
        \longrightarrow K

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

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
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
        Category of modules over Algebra of differentiable scalar fields on the
         2-dimensional differentiable manifold S^2

    The parent of `t` is not a free module, for the sphere `S^2` is not
    parallelizable::

        sage: isinstance(t.parent(), FiniteRankFreeModule)
        False

    To fully define `t`, we have to specify its components in some vector
    frames defined on subsets of `S^2`; let us start by the open subset `U`::

        sage: eU = c_xy.frame()
        sage: t[eU,:] = [[1,0], [-2,3]]
        sage: t.display(eU)
        t = dx*dx - 2 dy*dx + 3 dy*dy

    To set the components of `t` on `V` consistently, we copy the expressions
    of the components in the common subset `W`::

        sage: eV = c_uv.frame()
        sage: eVW = eV.restrict(W)
        sage: c_uvW = c_uv.restrict(W)
        sage: t[eV,0,0] = t[eVW,0,0,c_uvW].expr()  # long time
        sage: t[eV,0,1] = t[eVW,0,1,c_uvW].expr()  # long time
        sage: t[eV,1,0] = t[eVW,1,0,c_uvW].expr()  # long time
        sage: t[eV,1,1] = t[eVW,1,1,c_uvW].expr()  # long time

    Actually, the above operation can by performed in a single line by means
    of the method
    :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.add_comp_by_continuation`::

        sage: t.add_comp_by_continuation(eV, W, chart=c_uv)  # long time

    At this stage, `t` is fully defined, having components in frames eU and eV
    and the union of the domains of eU and eV being whole manifold::

        sage: t.display(eV)  # long time
        t = (u^4 - 4*u^3*v + 10*u^2*v^2 + 4*u*v^3 + v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du*du
         - 4*(u^3*v + 2*u^2*v^2 - u*v^3)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du*dv
         + 2*(u^4 - 2*u^3*v - 2*u^2*v^2 + 2*u*v^3 + v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv*du
         + (3*u^4 + 4*u^3*v - 2*u^2*v^2 - 4*u*v^3 + 3*v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv*dv

    Let us consider two vector fields, `a` and `b`, on `S^2`::

        sage: a = M.vector_field(name='a')
        sage: a[eU,:] = [-y,x]
        sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: a.display(eV)
        a = -v d/du + u d/dv
        sage: b = M.vector_field(name='b')
        sage: b[eU,:] = [y,-1]
        sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: b.display(eV)
        b = ((2*u + 1)*v^3 + (2*u^3 - u^2)*v)/(u^2 + v^2) d/du
         - (u^4 - v^4 + 2*u*v^2)/(u^2 + v^2) d/dv

    As a tensor field of type `(0,2)`, `t` acts on the pair `(a,b)`,
    resulting in a scalar field::

        sage: f = t(a,b); f
        Scalar field t(a,b) on the 2-dimensional differentiable manifold S^2
        sage: f.display()  # long time
        t(a,b): S^2 --> R
        on U: (x, y) |--> -2*x*y - y^2 - 3*x
        on V: (u, v) |--> -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)

    The vectors can be defined only on subsets of `S^2`, the domain of the
    result is then the common subset::

        sage: s = t(a.restrict(U), b) ; s  # long time
        Scalar field t(a,b) on the Open subset U of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()  # long time
        t(a,b): U --> R
           (x, y) |--> -2*x*y - y^2 - 3*x
        on W: (u, v) |--> -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)
        sage: s = t(a.restrict(U), b.restrict(W)) ; s  # long time
        Scalar field t(a,b) on the Open subset W of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()  # long time
        t(a,b): W --> R
           (x, y) |--> -2*x*y - y^2 - 3*x
           (u, v) |--> -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)

    The tensor itself can be defined only on some open subset of `S^2`,
    yielding a result whose domain is this subset::

        sage: s = t.restrict(V)(a,b); s  # long time
        Scalar field t(a,b) on the Open subset V of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()  # long time
        t(a,b): V --> R
           (u, v) |--> -(3*u^3 + (3*u + 1)*v^2 + 2*u*v)/(u^4 + 2*u^2*v^2 + v^4)
        on W: (x, y) |--> -2*x*y - y^2 - 3*x

    Tests regarding the multiplication by a scalar field::

        sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2),
        ....:                     c_uv: (u^2 + v^2)/(u^2 + v^2 + 1)}, name='f')
        sage: t.parent().base_ring() is f.parent()
        True
        sage: s = f*t; s  # long time
        Tensor field of type (0,2) on the 2-dimensional differentiable
         manifold S^2
        sage: s[[0,0]] == f*t[[0,0]]  # long time
        True
        sage: s.restrict(U) == f.restrict(U) * t.restrict(U)  # long time
        True
        sage: s = f*t.restrict(U); s
        Tensor field of type (0,2) on the Open subset U of the 2-dimensional
         differentiable manifold S^2
        sage: s.restrict(U) == f.restrict(U) * t.restrict(U)
        True

    """
    def __init__(self, vector_field_module, tensor_type, name=None,
                 latex_name=None, sym=None, antisym=None, parent=None):
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
            t = (x^2 + 1) dx*dx + x*y dx*dy + (y^2 + 1) dy*dy
            sage: t.display(e_uv)
            t = (3/16*u^2 + 1/16*v^2 + 1/2) du*du
             + (-1/16*u^2 + 1/4*u*v + 1/16*v^2) du*dv
             + (1/16*u^2 + 1/4*u*v - 1/16*v^2) dv*du
             + (1/16*u^2 + 3/16*v^2 + 1/2) dv*dv
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
        ModuleElement.__init__(self, parent)
        self._vmodule = vector_field_module
        self._tensor_type = tuple(tensor_type)
        self._tensor_rank = self._tensor_type[0] + self._tensor_type[1]
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        self._restrictions = {} # dict. of restrictions of self on subdomains
                                # of self._domain, with the subdomains as keys
        # Treatment of symmetry declarations:
        self._sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):
                # a single symmetry is provided as a tuple -> 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) > 1:
                    for i in isym:
                        if i < 0 or i > self._tensor_rank - 1:
                            raise IndexError("invalid position: {}".format(i) +
                                 " not in [0,{}]".format(self._tensor_rank-1))
                    self._sym.append(tuple(isym))
        self._antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):
                # a single antisymmetry is provided as a tuple -> 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) > 1:
                    for i in isym:
                        if i < 0 or i > self._tensor_rank - 1:
                            raise IndexError("invalid position: {}".format(i) +
                                " not in [0,{}]".format(self._tensor_rank-1))
                    self._antisym.append(tuple(isym))
        # Final consistency check:
        index_list = []
        for isym in self._sym:
            index_list += isym
        for isym in self._antisym:
            index_list += isym
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("incompatible lists of symmetries: the same " +
                             "position appears more than once")
        # Initialization of derived quantities:
        self._init_derived()

    ####### Required methods for ModuleElement (beside arithmetic) #######

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
        return any(bool(rst) for rst in self._restrictions.values())

    __nonzero__ = __bool__  # For Python2 compatibility

    ##### End of required methods for ModuleElement (beside arithmetic) #####

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
        if self._tensor_type == (0,2) and self._sym == [(0,1)]:
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
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def set_name(self, name=None, latex_name=None):
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

    def _final_repr(self, description):
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
        for vid, val in self._lie_derivatives.items():
            del val[0]._lie_der_along_self[id(self)]
        # Then clears the dictionary of Lie derivatives
        self._lie_derivatives.clear()

    #### Simple accessors ####

    def domain(self):
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

    def base_module(self):
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

    def tensor_type(self):
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
            s = "symmetries: {}; ".format(self._sym)
        if not self._antisym:
            a = "no antisymmetry"
        elif len(self._antisym) == 1:
            a = "antisymmetry: {}".format(self._antisym[0])
        else:
            a = "antisymmetries: {}".format(self._antisym)
        print(s + a)

    #### End of simple accessors #####

    def set_restriction(self, rst):
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
            t = (x + y) d/dx*dx*dy
            sage: t.restrict(U) == s
            True

        """
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
        self._restrictions[rst._domain] = rst.copy()
        self._restrictions[rst._domain].set_name(name=self._name,
                                                 latex_name=self._latex_name)

    def restrict(self, subdomain, dest_map=None):
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

        OUTPUT:

        - :class:`TensorField` representing the restriction

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
            sage: v = M.vector_field('v')
            sage: v.set_comp(eN)[1] = 1  # given the default settings, this can be abriged to v[1] = 1
            sage: v.display()
            v = d/dx
            sage: vU = v.restrict(U) ; vU
            Vector field v on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: vU.display()
            v = d/dx
            sage: vU == eN[1]
            True
            sage: vW = v.restrict(W) ; vW
            Vector field v on the Open subset W of the 2-dimensional
             differentiable manifold S^2
            sage: vW.display()
            v = d/dx
            sage: vW.display(eS_W, stereoS_W)
            v = (-u^2 + v^2) d/du - 2*u*v d/dv
            sage: vW == eN_W[1]
            True

        At this stage, defining the restriction of ``v`` to the open
        subset ``V`` fully specifies ``v``::

            sage: v.restrict(V)[1] = vW[eS_W, 1, stereoS_W].expr()  # note that eS is the default frame on V
            sage: v.restrict(V)[2] = vW[eS_W, 2, stereoS_W].expr()
            sage: v.display(eS, stereoS)
            v = (-u^2 + v^2) d/du - 2*u*v d/dv
            sage: v.restrict(U).display()
            v = d/dx
            sage: v.restrict(V).display()
            v = (-u^2 + v^2) d/du - 2*u*v d/dv

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
                if subdomain.is_subset(dom):
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    break
            # If this fails, the restriction is created from scratch:
            else:
                smodule = subdomain.vector_field_module(dest_map=dest_map)
                self._restrictions[subdomain] = smodule.tensor(
                                                  self._tensor_type,
                                                  name=self._name,
                                                  latex_name=self._latex_name,
                                                  sym=self._sym,
                                                  antisym=self._antisym,
                                                  specific_type=type(self))
        return self._restrictions[subdomain]

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
            3-indices components w.r.t. Coordinate frame (V, (d/du,d/dv))
            sage: t.set_comp(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) d/dv*du*dv

        Setting the components in a new frame (``e``)::

            sage: e = V.vector_frame('e')
            sage: t.set_comp(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t.set_comp(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0*e^1*e^1

        Since the frames ``e`` and ``e_uv`` are defined on the same domain, the
        components w.r.t. ``e_uv`` have been erased::

            sage: t.display(c_uv.frame())
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (V, (d/du,d/dv))

        """
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.set_comp(basis)

    def add_comp(self, basis=None):
        r"""
        Return the components of ``self`` in a given vector frame
        for assignment.

        The components with respect to other frames having the same domain
        as the provided vector frame are kept. To delete them them, use the
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
            3-indices components w.r.t. Coordinate frame (V, (d/du,d/dv))
            sage: t.add_comp(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) d/dv*du*dv

        Setting the components in a new frame::

            sage: e = V.vector_frame('e')
            sage: t.add_comp(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t.add_comp(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0*e^1*e^1

        The components with respect to ``e_uv`` are kept::

            sage: t.display(e_uv)
            t = (u + v) d/dv*du*dv

        """
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.add_comp(basis)

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
            sage: # The two open subsets covered by stereographic coordinates (North and South):
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coordinates
            sage: transf = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:             intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:             restrictions2= u^2+v^2!=0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.vector_field('a')
            sage: a[eU,:] = [x, 2+y]

        At this stage, the vector field has been defined only on the open
        subset ``U`` (through its components in the frame ``eU``)::

            sage: a.display(eU)
            a = x d/dx + (y + 2) d/dy

        The components with respect to the restriction of ``eV`` to the common
        subdomain ``W``, in terms of the ``(u,v)`` coordinates, are obtained
        by a change-of-frame formula on ``W``::

            sage: a.display(eV.restrict(W), c_uv.restrict(W))
            a = (-4*u*v - u) d/du + (2*u^2 - 2*v^2 - v) d/dv

        The continuation consists in extending the definition of the vector
        field to the whole open subset ``V`` by demanding that the components
        in the frame eV have the same coordinate expression as the above one::

            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)

        We have then::

            sage: a.display(eV)
            a = (-4*u*v - u) d/du + (2*u^2 - 2*v^2 - v) d/dv

        and `a` is defined on the entire manifold `S^2`.

        """
        dom = frame._domain
        if not dom.is_subset(self._domain):
            raise ValueError("the vector frame is not defined on a subset " +
                             "of the tensor field domain")
        if chart is None:
            chart = dom._def_chart
        sframe = frame.restrict(subdomain)
        schart = chart.restrict(subdomain)
        scomp = self.comp(sframe)
        resu = self.add_comp(frame) # _del_derived is performed here
        for ind in resu.non_redundant_index_generator():
            resu[[ind]] = dom.scalar_field({chart: scomp[[ind]].expr(schart)})

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
            Coordinate frame (U, (d/dx,d/dy))
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: f = V.default_frame() ; f
            Coordinate frame (V, (d/du,d/dv))
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e,0,0] = - x + y^3
            sage: t[e,0,1] = 2+x
            sage: t[f,1,1] = - u*v
            sage: t.comp(e)
            2-indices components w.r.t. Coordinate frame (U, (d/dx,d/dy))
            sage: t.comp(e)[:]
            [y^3 - x   x + 2]
            [      0       0]
            sage: t.comp(f)
            2-indices components w.r.t. Coordinate frame (V, (d/du,d/dv))
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

    def display(self, basis=None, chart=None):
        r"""
        Display the tensor field in terms of its expansion with respect
        to a given vector frame.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame with respect to
          which the tensor is expanded; if ``None``, the default frame
          of the domain of definition of the tensor field is assumed
        - ``chart`` -- (default: ``None``) chart with respect to which the
          components of the tensor field in the selected frame are expressed;
          if ``None``, the default chart of the vector frame domain is assumed

        EXAMPLES:

        Display of a type-`(1,1)` tensor field defined on two open subsets::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: e = U.default_frame() ; e
            Coordinate frame (U, (d/dx,d/dy))
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: f = V.default_frame() ; f
            Coordinate frame (V, (d/du,d/dv))
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e,0,0] = - x + y^3
            sage: t[e,0,1] = 2+x
            sage: t[f,1,1] = - u*v
            sage: t.display(e)
            t = (y^3 - x) d/dx*dx + (x + 2) d/dx*dy
            sage: t.display(f)
            t = -u*v d/dv*dv

        Since ``e`` is ``M``'s default frame, the argument ``e`` can
        be omitted::

            sage: e is M.default_frame()
            True
            sage: t.display()
            t = (y^3 - x) d/dx*dx + (x + 2) d/dx*dy

        Similarly, since ``f`` is ``V``'s default frame, the argument ``f``
        can be omitted when considering the restriction of ``t`` to ``V``::

            sage: t.restrict(V).display()
            t = -u*v d/dv*dv

        Display with respect to a frame in which ``t`` has not been
        initialized (automatic use of a change-of-frame formula)::

            sage: a = V.automorphism_field()
            sage: a[:] = [[1+v, -u^2], [0, 1-u]]
            sage: h = f.new_frame(a, 'h')
            sage: t.display(h)
            t = -u^3*v/(v + 1) h_0*h^1 - u*v h_1*h^1

        A shortcut of ``display()`` is ``disp()``::

            sage: t.disp(h)
            t = -u^3*v/(v + 1) h_0*h^1 - u*v h_1*h^1

        """
        if basis is None:
            if self._vmodule._dest_map.is_identity():
                basis = self._domain._def_frame
            else:
                for rst in self._restrictions.values():
                    try:
                        return rst.display()
                    except ValueError:
                        pass
            if basis is None:  # should be "is still None" ;-)
                raise ValueError("a frame must be provided for the display")
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.display(basis, chart)

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
        - ``coordinate_labels`` -- (default: ``True``) boolean; if ``True``,
          coordinate symbols are used by default (instead of integers) as
          index labels whenever ``frame`` is a coordinate frame
        - ``only_nonzero`` -- (default: ``True``) boolean; if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- (default: ``False``) boolean; if ``True``,
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
            M --> R
            on U: (x, y) |--> (x + 1)*y + x

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
        Sets a component with respect to some vector frame.

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
            t = (y^2 + x) d/dx*dy
            sage: t.__setitem__((0, 1), x+y^2)  # same as above since e_xy is the default frame on M
            sage: t.display()
            t = (y^2 + x) d/dx*dy
            sage: t.__setitem__(slice(None), [[x+y, -2], [3*y^2, x*y]])
            sage: t.display()
            t = (x + y) d/dx*dx - 2 d/dx*dy + 3*y^2 d/dy*dx + x*y d/dy*dy

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

    def copy(self):
        r"""
        Return an exact copy of ``self``.

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
            Tensor field of type (1,1) on
             the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s == t
            True

        If the original tensor field is modified, the copy is not::

            sage: t[e_xy,0,0] = -1
            sage: t.display(e_xy)
            t = -d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s.display(e_xy)
            (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s == t
            False

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = rst.copy()
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

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

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
        if other is self:
            return True
        if other in ZZ: # to compare with 0
            if other == 0:
                return self.is_zero()
            return False
        elif not isinstance(other, TensorField):
            return False
        else: # other is another tensor field
            if other._vmodule != self._vmodule:
                return False
            if other._tensor_type != self._tensor_type:
                return False
            # Non-trivial open covers of the domain:
            open_covers = self._domain.open_covers()[1:]  # the open cover 0
                                                          # is trivial
            for oc in open_covers:
                resu = True
                for dom in oc:
                    try:
                        resu = resu and \
                                bool(self.restrict(dom) == other.restrict(dom))
                    except ValueError:
                        break
                else:
                    # If this point is reached, no exception has occured; hence
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

        OUTPUT:

        - an exact copy of ``self``

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
            +t = (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = + rst
        if self._name is not None:
            resu._name = '+' + self._name
        if self._latex_name is not None:
            resu._latex_name = '+' + self._latex_name
        return resu

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the tensor field `-T`, where `T` is ``self``

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
            t = x d/dx*dx - x d/dx*dy + y d/dy*dx - y d/dy*dy
            sage: t.display(e_uv)
            t = u d/du*dv + v d/dv*dv
            sage: s = t.__neg__(); s
            Tensor field -t of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            -t = -x d/dx*dx + x d/dx*dy - y d/dy*dx + y d/dy*dy
            sage: s.display(e_uv)
            -t = -u d/du*dv - v d/dv*dv
            sage: s == -t  # indirect doctest
            True

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = - rst
        if self._name is not None:
            resu._name = '-' + self._name
        if self._latex_name is not None:
            resu._latex_name = '-' + self._latex_name
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
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: b.display(e_xy)
            b = 2 d/dx*dx + y d/dx*dy + x d/dy*dx - x d/dy*dy
            sage: s.display(e_xy)
            a+b = (x + 2) d/dx*dx + (y + 1) d/dx*dy + (x + y) d/dy*dx - x d/dy*dy
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
        if other == 0:
            return +self
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] + other._restrictions[dom]
        some_rst = resu_rst.values()[0]
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
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: b.display(e_xy)
            b = 2 d/dx*dx + y d/dx*dy + x d/dy*dx - x d/dy*dy
            sage: s.display(e_xy)
            a-b = (x - 2) d/dx*dx + (-y + 1) d/dx*dy + (-x + y) d/dy*dx + x d/dy*dy
            sage: s == a - b
            True
            sage: z = a.parent().zero()
            sage: a._sub_(z) == a
            True
            sage: z._sub_(a) == -a
            True

        """
        if other == 0:
            return +self
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] - other._restrictions[dom]
        some_rst = resu_rst.values()[0]
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
        Reflected multiplication operator: performs ``scalar * self``

        This is actually the multiplication by an element of the ring over
        which the tensor field module is constructed.

        INPUT:

        - ``scalar`` -- scalar field in the scalar field algebra over which
          the module containing ``self`` is defined

        OUTPUT:

        - the tensor field ``scalar * self``

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
            f: M --> R
            on U: (x, y) |--> 1/(x^2 + y^2 + 1)
            on V: (u, v) |--> 2/(u^2 + v^2 + 2)
            sage: s = a._rmul_(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a.display(e_xy)
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: s.display(e_xy)
            x/(x^2 + y^2 + 1) d/dx*dx + 1/(x^2 + y^2 + 1) d/dx*dy + y/(x^2 + y^2 + 1) d/dy*dx
            sage: a.display(e_uv)
            a = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv
            sage: s.display(e_uv)
            (u + 1)/(u^2 + v^2 + 2) d/du*du + (u - 1)/(u^2 + v^2 + 2) d/du*dv + (v + 1)/(u^2 + v^2 + 2) d/dv*du + (v - 1)/(u^2 + v^2 + 2) d/dv*dv
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
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = scalar.restrict(dom) * rst
        return resu

    ######### End of ModuleElement arithmetic operators ########

    # TODO: Move to acted_upon or _rmul_
    def __mul__(self, other):
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
            a*b = x^2 d/dx*d/dx*dx + x d/dx*d/dx*dy + x*y d/dx*d/dy*dx
             + y d/dx*d/dy*dy + x*y d/dy*d/dx*dx + y^2 d/dy*d/dy*dx
            sage: s.display(e_uv)
            a*b = (1/2*u^2 + 1/2*u) d/du*d/du*du + (1/2*u^2 - 1/2*u) d/du*d/du*dv
             + 1/2*(u + 1)*v d/du*d/dv*du + 1/2*(u - 1)*v d/du*d/dv*dv
             + (1/2*u*v + 1/2*u) d/dv*d/du*du + (1/2*u*v - 1/2*u) d/dv*d/du*dv
             + (1/2*v^2 + 1/2*v) d/dv*d/dv*du + (1/2*v^2 - 1/2*v) d/dv*d/dv*dv

        Multiplication on the right by a scalar field::

            sage: f = M.scalar_field({c_xy: x*y}, name='f')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: s = a.__mul__(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            x^2*y d/dx*dx + x*y d/dx*dy + x*y^2 d/dy*dx
            sage: s.display(e_uv)
            (1/8*u^3 - 1/8*(u + 1)*v^2 + 1/8*u^2) d/du*du
             + (1/8*u^3 - 1/8*(u - 1)*v^2 - 1/8*u^2) d/du*dv
             + (1/8*u^2*v - 1/8*v^3 + 1/8*u^2 - 1/8*v^2) d/dv*du
             + (1/8*u^2*v - 1/8*v^3 - 1/8*u^2 + 1/8*v^2) d/dv*dv
            sage: s == f*a
            True

        Multiplication on the right by a number::

            sage: s = a.__mul__(2); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            2*x d/dx*dx + 2 d/dx*dy + 2*y d/dy*dx
            sage: s.display(e_uv)
            (u + 1) d/du*du + (u - 1) d/du*dv + (v + 1) d/dv*du
             + (v - 1) d/dv*dv
            sage: s.restrict(U) == 2*a.restrict(U)
            True
            sage: s.restrict(V) == 2*a.restrict(V)
            True
            sage: s == 2*a
            True

        """
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

    def __truediv__(self, scalar):
        r"""
        Division by a scalar field.

        INPUT:

        - ``scalar`` -- scalar field in the scalar field algebra over which
          the module containing ``self`` is defined

        OUTPUT:

        - the tensor field ``scalar * self``

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
            f: M --> R
            on U: (x, y) |--> 1/(x^2 + y^2 + 1)
            on V: (u, v) |--> 2/(u^2 + v^2 + 2)
            sage: s = a.__div__(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            (x^3 + x*y^2 + x) d/dx*dx + (x^2 + y^2 + 1) d/dx*dy
             + (y^3 + (x^2 + 1)*y) d/dy*dx
            sage: f*s == a
            True

        Division by a number::

            sage: s = a.__div__(2); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            1/2*x d/dx*dx + 1/2 d/dx*dy + 1/2*y d/dy*dx
            sage: s.display(e_uv)
            (1/4*u + 1/4) d/du*du + (1/4*u - 1/4) d/du*dv
             + (1/4*v + 1/4) d/dv*du + (1/4*v - 1/4) d/dv*dv
            sage: s == a/2
            True
            sage: 2*s == a
            True

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = rst / scalar
        return resu

    __div__ = __truediv__

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
            t(a,w): M --> R
            on U: (x, y) |--> x*y^4 + x*y^3 + x^2*y - x*y^2 - x^2
            on V: (u, v) |--> 1/32*u^5 - 1/32*(3*u + 2)*v^4 + 1/32*v^5
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
            t(w) = (x*y^2 + x^2) d/dx + y^3 d/dy
            sage: s.display(e_uv)
            t(w) = (1/4*u^3 + 1/4*(u + 1)*v^2 + 1/4*u^2 - 1/2*(u^2 - u)*v) d/du
             + (-1/4*(2*u - 1)*v^2 + 1/4*v^3 + 1/4*u^2 + 1/4*(u^2 + 2*u)*v) d/dv
            sage: s.restrict(U) == t.restrict(U)(w.restrict(U))
            True
            sage: s.restrict(V) == t.restrict(V)(w.restrict(V))
            True

        """
        p = len(args)
        if p == 1 and self._tensor_type == (1,1):
            # type-(1,1) tensor acting as a a field of tangent-space
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
                if resu_rr == 0:
                    for chart in resu_rr._domain._atlas:
                        resu._express[chart] = chart._zero_function
                else:
                    for chart, expr in resu_rr._express.items():
                        resu._express[chart] = expr
            if resu == 0:
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

    def trace(self, pos1=0, pos2=1):
        r"""
        Trace (contraction) on two slots of the tensor field.

        INPUT:

        - ``pos1`` -- (default: 0) position of the first index for the
          contraction, with the convention ``pos1=0`` for the first slot
        - ``pos2`` -- (default: 1) position of the second index for the
          contraction, with the same convention as for ``pos1``. The variance
          type of ``pos2`` must be opposite to that of ``pos1``

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
            M --> R
            on U: (x, y) |--> y + 1
            on V: (u, v) |--> 1/2*u - 1/2*v + 1
            sage: s == a.trace(0,1) # explicit mention of the positions
            True

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

    def contract(self, *args):
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
            sage: a = M.tensor_field(1,1, name='a')
            sage: a[eU,:] = [[1,x], [0,2]]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: b = M.tensor_field(2,0, name='b')
            sage: b[eU,:] = [[y,-1], [x+y,2]]
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = a.contract(b) ; s   # contraction on last index of a and first one of b
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M

        Check 1: components with respect to the manifold's default
        frame (``eU``)::

            sage: [[bool(s[i,j] == sum(a[i,k]*b[k,j] for k in M.irange()))
            ....:   for j in M.irange()] for i in M.irange()]
            [[True, True], [True, True]]

        Check 2: components with respect to the frame ``eV``::

            sage: [[bool(s[eV,i,j] == sum(a[eV,i,k]*b[eV,k,j]
            ....:                         for k in M.irange()))
            ....:   for j in M.irange()] for i in M.irange()]
            [[True, True], [True, True]]

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
            sage: s == c['^.._kl']*b['^kl']  # the same double contraction in index notation; long time
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

            sage: a = M.one_form('a')
            sage: a[eU,:] = [y, 1+x]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: b = M.vector_field('b')
            sage: b[eU,:] = [x, y^2]
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a.display(eU)
            a = y dx + (x + 1) dy
            sage: b.display(eU)
            b = x d/dx + y^2 d/dy
            sage: s = a.contract(b) ; s
            Scalar field on the 2-dimensional differentiable manifold M
            sage: s.display()
            M --> R
            on U: (x, y) |--> (x + 1)*y^2 + x*y
            on V: (u, v) |--> 1/8*u^3 - 1/8*u*v^2 + 1/8*v^3 + 1/2*u^2 - 1/8*(u^2 + 4*u)*v
            sage: s == a['_i']*b['^i'] # use of index notation
            True
            sage: s == b.contract(a)
            True

        Case of a vanishing scalar field result::

            sage: b = M.vector_field('b')
            sage: b[eU,:] = [1+x, -y]
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = a.contract(b) ; s
            Scalar field zero on the 2-dimensional differentiable manifold M
            sage: s.display()
            zero: M --> R
            on U: (x, y) |--> 0
            on V: (u, v) |--> 0

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
            sage: a = M.tensor_field(0,2, name='a')
            sage: a[eU,:] = [[1,x], [2,y]]
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
            sage: a = M.tensor_field(0,2, name='a')
            sage: a[eU,:] = [[1,x], [2,y]]
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
        the 2-sphere::

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
            sage: w[e_xy,:] = [-y, x]
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: lt = t.lie_derivative(w); lt
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: lt.display(e_xy)
            d/dx*dx - x d/dx*dy + (-y - 1) d/dy*dy
            sage: lt.display(e_uv)
            -1/2*u d/du*du + (1/2*u + 1) d/du*dv + (-1/2*v + 1) d/dv*du + 1/2*v d/dv*dv

        The result is cached::

            sage: t.lie_derivative(w) is lt
            True

        An alias is ``lie_der``::

            sage: t.lie_der(w) is t.lie_derivative(w)
            True

        Lie derivative of a vector field::

            sage: a = M.vector_field(name='a')
            sage: a[e_xy,:] = [1-x, x-y]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: a.lie_der(w)
            Vector field on the 2-dimensional differentiable manifold M
            sage: a.lie_der(w).display(e_xy)
            x d/dx + (-y - 1) d/dy
            sage: a.lie_der(w).display(e_uv)
            (v - 1) d/du + (u + 1) d/dv

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
                resu_rst.append(rst.lie_der(vector.restrict(dom)))
            resu = self._vmodule.tensor(self._tensor_type,
                                        sym=resu_rst[0]._sym,
                                        antisym=resu_rst[0]._antisym)
            for rst in resu_rst:
                resu._restrictions[rst._domain] = rst
            self._lie_derivatives[id(vector)] = (vector, resu)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]

    lie_der = lie_derivative

    def at(self, point):
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
            sage: a = M.tensor_field(1,1, name='a')
            sage: a[eU,:] = [[1+y,x], [0,x+y]]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a.display(eU)
            a = (y + 1) d/dx*dx + x d/dx*dy + (x + y) d/dy*dy
            sage: a.display(eV)
            a = (u + 1/2) d/du*du + (-1/2*u - 1/2*v + 1/2) d/du*dv
             + 1/2 d/dv*du + (1/2*u - 1/2*v + 1/2) d/dv*dv
            sage: p = M.point((2,3), chart=c_xy, name='p')
            sage: ap = a.at(p) ; ap
            Type-(1,1) tensor a on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: ap.parent()
            Free module of type-(1,1) tensors on the Tangent space at Point p
             on the 2-dimensional differentiable manifold M
            sage: ap.display(eU.at(p))
            a = 4 d/dx*dx + 2 d/dx*dy + 5 d/dy*dy
            sage: ap.display(eV.at(p))
            a = 11/2 d/du*du - 3/2 d/du*dv + 1/2 d/dv*du + 7/2 d/dv*dv
            sage: p.coord(c_uv) # to check the above expression
            (5, -1)

        """
        if point not in self._domain:
            raise ValueError("the {} is not a point in the ".format(point) +
                             "domain of {}".format(self))
        for dom, rst in self._restrictions.items():
            if point in dom:
                return rst.at(point)

    def up(self, metric, pos=None):
        r"""
        Compute a metric dual of the tensor field by raising some index with a
        given metric.

        If `T` is the tensor field, `(k,l)` its type and `p` the position of a
        covariant index (i.e. `k\leq p < k+l`), this method called with
        ``pos`` `=p` yields the tensor field `T^\sharp` of type `(k+1,l-1)`
        whose components are

        .. MATH::

            (T^\sharp)^{a_1\ldots a_{k+1}}_{\qquad\ \ b_1 \ldots b_{l-1}}
            = g^{a_{k+1} i} \,
            T^{a_1\ldots a_k}_{\qquad\   b_1 \ldots b_{p-k} \, i \, b_{p-k+1}\ldots b_{l-1}},

        `g^{ab}` being the components of the inverse metric.

        The reverse operation is :meth:`TensorField.down`

        INPUT:

        - ``metric`` -- metric `g`, as an instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
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
            sage: w = M.one_form()
            sage: w[:] = [-1, 2]
            sage: v = w.up(g) ; v
            Vector field on the 2-dimensional differentiable manifold M
            sage: v.display()
            ((2*x - 1)*y + 1)/(x^2*y^2 + (x + 1)*y - x - 1) d/dx
             - (x*y + 2*x + 2)/(x^2*y^2 + (x + 1)*y - x - 1) d/dy
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

            sage: t = M.tensor_field(0, 2)
            sage: t[:] = [[1,2], [3,4]]
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
                result = result.up(metric, k)
            return result
        if not isinstance(pos, (int, Integer)):
            raise TypeError("the argument 'pos' must be an integer")
        if pos<n_con or pos>self._tensor_rank-1:
            print("pos = {}".format(pos))
            raise ValueError("position out of range")
        return self.contract(pos, metric.inverse(), 0)

    def down(self, metric, pos=None):
        r"""
        Compute a metric dual of the tensor field by lowering some index with a
        given metric.

        If `T` is the tensor field, `(k,l)` its type and `p` the position of a
        contravariant index (i.e. `0\leq p < k`), this method called with
        ``pos`` `=p`  yields the tensor field `T^\flat` of type `(k-1,l+1)`
        whose components are

        .. MATH::

            (T^\flat)^{a_1\ldots a_{k-1}}_{\qquad\ \  b_1 \ldots b_{l+1}}
            = g_{b_1 i} \,
            T^{a_1\ldots a_{p} \, i \, a_{p+1}\ldots a_{k-1}}_{\qquad\qquad\quad\; b_2 \ldots b_{l+1}},

        `g_{ab}` being the components of the metric tensor.

        The reverse operation is :meth:`TensorField.up`

        INPUT:

        - ``metric`` -- metric `g`, as an instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
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
            sage: v = M.vector_field()
            sage: v[:] = [-1,2]
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

            sage: t = M.tensor_field(2, 0)
            sage: t[:] = [[1,2], [3,4]]
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
        n_con = self._tensor_type[0] # number of contravariant indices = k
        if pos is None:
            result = self
            for p in range(0, n_con):
                k = result._tensor_type[0]
                result = result.down(metric, k-1)
            return result
        if not isinstance(pos, (int, Integer)):
            raise TypeError("the argument 'pos' must be an integer")
        if pos<0 or pos>=n_con:
            print("pos = {}".format(pos))
            raise ValueError("position out of range")
        return metric.contract(1, self, pos)
