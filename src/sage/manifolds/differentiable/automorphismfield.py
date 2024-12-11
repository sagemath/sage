r"""
Tangent-Space Automorphism Fields

The class :class:`AutomorphismField` implements fields of automorphisms of
tangent spaces to a generic (a priori not parallelizable) differentiable
manifold, while the class :class:`AutomorphismFieldParal`
is devoted to fields of automorphisms of tangent spaces to a parallelizable
manifold. The latter play the important role of transitions between vector
frames sharing the same domain on a differentiable manifold.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks
"""

# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism


class AutomorphismField(TensorField):
    r"""
    Field of automorphisms of tangent spaces to a generic (a priori
    not parallelizable) differentiable manifold.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi: U \rightarrow M` to a differentiable manifold `M`,
    a *field of tangent-space automorphisms along* `U` *with values on*
    `M \supset\Phi(U)` is a differentiable map

    .. MATH::

        a:\ U  \longrightarrow T^{(1,1)} M,

    with `T^{(1,1)} M` being the tensor bundle of type `(1,1)` over `M`,
    such that

    .. MATH::

        \forall p \in U,\ a(p) \in \mathrm{Aut}(T_{\Phi(p)} M),

    i.e. `a(p)` is an automorphism of the tangent space to `M` at the
    point `\Phi(p)`.

    The standard case of a field of tangent-space automorphisms *on* a
    manifold corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other
    common cases are `\Phi` being an immersion and `\Phi` being a curve
    in `M` (`U` is then an open interval of `\RR`).

    .. NOTE::

        If `M` is parallelizable, then :class:`AutomorphismFieldParal`
        *must* be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `M` via the map `\Phi`
    - ``name`` -- (default: ``None``) name given to the field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the field;
      if none is provided, the LaTeX symbol is set to ``name``
    - ``is_identity`` -- boolean (default: ``False``); determines whether the
      constructed object is a field of identity automorphisms

    EXAMPLES:

    Field of tangent-space automorphisms on a non-parallelizable
    2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
        ....:                              restrictions1= x>0, restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: a = M.automorphism_field(name='a') ; a
        Field of tangent-space automorphisms a on the 2-dimensional
         differentiable manifold M
        sage: a.parent()
        General linear group of the Module X(M) of vector fields on the
         2-dimensional differentiable manifold M

    We first define the components of `a` with respect to the
    coordinate frame on `U`::

        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: a[eU,:] = [[1,x], [0,2]]

    It is equivalent to pass the components while defining `a`::

        sage: a = M.automorphism_field({eU: [[1,x], [0,2]]}, name='a')

    We then set the components with respect to the coordinate frame
    on `V` by extending the expressions of the components in the
    corresponding subframe on `W = U \cap V`::

        sage: W = U.intersection(V)
        sage: a.add_comp_by_continuation(eV, W, c_uv)

    At this stage, the automorphism field `a` is fully defined::

        sage: a.display(eU)
        a = ∂/∂x⊗dx + x ∂/∂x⊗dy + 2 ∂/∂y⊗dy
        sage: a.display(eV)
        a = (1/4*u + 1/4*v + 3/2) ∂/∂u⊗du + (-1/4*u - 1/4*v - 1/2) ∂/∂u⊗dv
         + (1/4*u + 1/4*v - 1/2) ∂/∂v⊗du + (-1/4*u - 1/4*v + 3/2) ∂/∂v⊗dv

    In particular, we may ask for its inverse on the whole manifold `M`::

        sage: ia = a.inverse() ; ia
        Field of tangent-space automorphisms a^(-1) on the 2-dimensional
         differentiable manifold M
        sage: ia.display(eU)
        a^(-1) = ∂/∂x⊗dx - 1/2*x ∂/∂x⊗dy + 1/2 ∂/∂y⊗dy
        sage: ia.display(eV)
        a^(-1) = (-1/8*u - 1/8*v + 3/4) ∂/∂u⊗du + (1/8*u + 1/8*v + 1/4) ∂/∂u⊗dv
         + (-1/8*u - 1/8*v + 1/4) ∂/∂v⊗du + (1/8*u + 1/8*v + 3/4) ∂/∂v⊗dv

    Equivalently, one can use the power minus one to get the inverse::

        sage: ia is a^(-1)
        True

    or the operator ``~``::

        sage: ia is ~a
        True
    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a field of tangent-space automorphisms on a
        non-parallelizable manifold.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``AutomorphismField``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: XM = M.vector_field_module()
            sage: GL = XM.general_linear_group()
            sage: a = GL.element_class(XM, name='a'); a
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold M
            sage: a[c_xy.frame(), :] = [[1+x^2, 0], [0, 1+y^2]]
            sage: a.add_comp_by_continuation(c_uv.frame(), U.intersection(V), c_uv)
            sage: TestSuite(a).run(skip='_test_pickling')

        Construction of the identity field::

            sage: b = GL.one(); b
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold M
            sage: TestSuite(b).run(skip='_test_pickling')

        Construction with ``DifferentiableManifold.automorphism_field``::

            sage: a1 = M.automorphism_field(name='a'); a1
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold M
            sage: type(a1) == type(a)
            True

        .. TODO::

            Fix ``_test_pickling`` (in the superclass :class:`TensorField`).
        """
        TensorField.__init__(self, vector_field_module, (1,1), name=name,
                             latex_name=latex_name,
                             parent=vector_field_module.general_linear_group())
        self._is_identity = False # a priori
        self._init_derived() # initialization of derived quantities

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: a = M.automorphism_field(name='a')
            sage: a._repr_()
            'Field of tangent-space automorphisms a on the 2-dimensional differentiable manifold M'
            sage: repr(a)  # indirect doctest
            'Field of tangent-space automorphisms a on the 2-dimensional differentiable manifold M'
            sage: a  # indirect doctest
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold M
        """
        description = "Field of tangent-space "
        if self is self.parent().one():
            description += "identity maps "
        else:
            description += "automorphisms "
            if self._name is not None:
                description += self._name + " "
        return self._final_repr(description)

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: a = M.automorphism_field(name='a')
            sage: a._init_derived()
        """
        TensorField._init_derived(self)
        self._inverse = None  # inverse not set yet

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: a = M.automorphism_field(name='a')
            sage: a._del_derived()
        """
        # First delete the derived quantities pertaining to the mother class:
        TensorField._del_derived(self)
        # then deletes the inverse automorphism:
        self._inverse = None

    def set_comp(self, basis=None):
        r"""
        Return the components of ``self`` w.r.t. a given module basis for
        assignment.

        The components with respect to other bases are deleted, in order to
        avoid any inconsistency. To keep them, use the method :meth:`add_comp`
        instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        EXAMPLES::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: a= M.automorphism_field(name='a')
            sage: a.set_comp(e_uv)
            2-indices components w.r.t. Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: a.set_comp(e_uv)[0,0] = u+v
            sage: a.set_comp(e_uv)[1,1] = u+v
            sage: a.display(e_uv)
            a = (u + v) ∂/∂u⊗du + (u + v) ∂/∂v⊗dv

        Setting the components in a new frame::

            sage: e = V.vector_frame('e')
            sage: a.set_comp(e)
            2-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: a.set_comp(e)[0,1] = u*v
            sage: a.set_comp(e)[1,0] = u*v
            sage: a.display(e)
            a = u*v e_0⊗e^1 + u*v e_1⊗e^0

        Since the frames ``e`` and ``e_uv`` are defined on the same domain, the
        components w.r.t. ``e_uv`` have been erased::

            sage: a.display(c_uv.frame())
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (V, (∂/∂u,∂/∂v))

        Since the identity map is an immutable element, its components
        cannot be changed::

            sage: id = M.tangent_identity_field()
            sage: id.add_comp(e)[0,1] = u*v
            Traceback (most recent call last):
            ...
            ValueError: the components of an immutable element cannot be
             changed
        """
        comp = super().set_comp(basis=basis)
        self._is_identity = False  # a priori
        return comp

    def add_comp(self, basis=None):
        r"""
        Return the components of ``self`` w.r.t. a given module basis for
        assignment, keeping the components w.r.t. other bases.

        To delete the components w.r.t. other bases, use the method
        :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        .. WARNING::

            If the automorphism field has already components in other bases, it
            is the user's responsibility to make sure that the components
            to be added are consistent with them.

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`;
          if such components did not exist previously, they are created

        EXAMPLES::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: a= M.automorphism_field(name='a')
            sage: a.add_comp(e_uv)
            2-indices components w.r.t. Coordinate frame (V, (∂/∂u,∂/∂v))
            sage: a.add_comp(e_uv)[0,0] = u+v
            sage: a.add_comp(e_uv)[1,1] = u+v
            sage: a.display(e_uv)
            a = (u + v) ∂/∂u⊗du + (u + v) ∂/∂v⊗dv

        Setting the components in a new frame::

            sage: e = V.vector_frame('e')
            sage: a.add_comp(e)
            2-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: a.add_comp(e)[0,1] = u*v
            sage: a.add_comp(e)[1,0] = u*v
            sage: a.display(e)
            a = u*v e_0⊗e^1 + u*v e_1⊗e^0

        The components with respect to ``e_uv`` are kept::

            sage: a.display(e_uv)
            a = (u + v) ∂/∂u⊗du + (u + v) ∂/∂v⊗dv

        Since the identity map is a special element, its components cannot be
        changed::

            sage: id = M.tangent_identity_field()
            sage: id.add_comp(e)[0,1] = u*v
            Traceback (most recent call last):
            ...
            ValueError: the components of an immutable element cannot be
             changed
        """
        comp = super().add_comp(basis=basis)
        self._is_identity = False  # a priori
        return comp

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same
        vector field module.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: a = M.automorphism_field(name='a')
            sage: a._new_instance()
            Field of tangent-space automorphisms on the 5-dimensional
             differentiable manifold M
            sage: a._new_instance().parent() is a.parent()
            True
        """
        return type(self)(self._vmodule)

    def __call__(self, *arg):
        r"""
        Redefinition of
        :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.__call__`
        to allow for a proper treatment of the identity map and of the call
        with a single argument

        TESTS:

        Field of identity maps on the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: w = M.vector_field({e_xy: [3, 1]}, name='w')
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: z = M.one_form({e_xy: [-y, x]}, name='z')
            sage: z.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: Id = M.tangent_identity_field()
            sage: s = Id(w); s
            Vector field w on the 2-dimensional differentiable manifold M
            sage: s == w
            True
            sage: s = Id(z, w); s
            Scalar field z(w) on the 2-dimensional differentiable manifold M
            sage: s == z(w)
            True

        Field of automorphisms on the 2-sphere::

            sage: a = M.automorphism_field({e_xy: [[-1, 0], [0, 1]]}, name='a')
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)

        Call with a single argument::

            sage: s = a(w); s
            Vector field a(w) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            a(w) = -3 ∂/∂x + ∂/∂y
            sage: s.display(e_uv)
            a(w) = (3*u^2 - 2*u*v - 3*v^2) ∂/∂u + (u^2 + 6*u*v - v^2) ∂/∂v
            sage: s.restrict(U) == a.restrict(U)(w.restrict(U))
            True
            sage: s.restrict(V) == a.restrict(V)(w.restrict(V))
            True
            sage: s.restrict(U) == a(w.restrict(U))
            True
            sage: s.restrict(U) == a.restrict(U)(w)
            True

        Call with two arguments::

            sage: s = a(z, w); s
            Scalar field a(z,w) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(z,w): M → ℝ
            on U: (x, y) ↦ x + 3*y
            on V: (u, v) ↦ (u + 3*v)/(u^2 + v^2)
            sage: s.restrict(U) == a.restrict(U)(z.restrict(U), w.restrict(U))
            True
            sage: s.restrict(V) == a.restrict(V)(z.restrict(V), w.restrict(V))
            True
            sage: s.restrict(U) == a(z.restrict(U), w.restrict(U))
            True
            sage: s.restrict(U) == a(z, w.restrict(U))
            True
        """
        if self._is_identity:
            if len(arg) == 1:
                # The identity map acting as such, on a vector field:
                vector = arg[0]
                if vector._tensor_type != (1,0):
                    raise TypeError("the argument must be a vector field")
                dom = self._domain.intersection(vector._domain)
                return vector.restrict(dom)
            elif len(arg) == 2:
                # self acting as a type-(1,1) tensor on a pair
                # (1-form, vector field), returning a scalar field:
                oneform = arg[0]
                vector = arg[1]
                dom = self._domain.intersection(
                                  oneform._domain).intersection(vector._domain)
                return oneform.restrict(dom)(vector.restrict(dom))
            else:
                raise TypeError("wrong number of arguments")
        # Generic case
        if len(arg) == 1:
            # The field of automorphisms acting on a vector field:
            vector = arg[0]
            if vector._tensor_type != (1,0):
                raise TypeError("the argument must be a vector field")
            dom = self._domain.intersection(vector._domain)
            vector_dom = vector.restrict(dom)
            if dom != self._domain:
                return self.restrict(dom)(vector_dom)
            resu = dom.vector_field()
            if self._name is not None and vector._name is not None:
                resu._name = self._name + "(" + vector._name + ")"
            if self._latex_name is not None and vector._latex_name is not None:
                resu._latex_name = self._latex_name + r"\left(" + \
                                   vector._latex_name + r"\right)"
            for sdom, automorph in self._restrictions.items():
                resu._restrictions[sdom] = automorph(vector_dom.restrict(sdom))
            return resu
        # Case of 2 arguments:
        return TensorField.__call__(self, *arg)

    def copy(self, name=None, latex_name=None):
        r"""
        Return an exact copy of the automorphism field ``self``.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the copy
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          copy; if none is provided, the LaTeX symbol is set to ``name``

        .. NOTE::

            The name and the derived quantities are not copied.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: Id = M.tangent_identity_field(); Id
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold M
            sage: one = Id.copy('1'); one
            Field of tangent-space automorphisms 1 on the 2-dimensional
             differentiable manifold M
        """
        copy = super().copy(name=name, latex_name=latex_name)
        copy._is_identity = self._is_identity
        return copy

    #### MultiplicativeGroupElement methods ####

    def __invert__(self):
        r"""
        Return the inverse automorphism of ``self``.

        EXAMPLES:

        Inverse of a field of tangent-space automorphisms on a
        non-parallelizable 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: W = U.intersection(V)
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:    intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.automorphism_field({eU: [[1,x], [0,2]]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, c_uv)
            sage: ia = a.inverse() ; ia
            Field of tangent-space automorphisms a^(-1) on the 2-dimensional
             differentiable manifold M
            sage: a[eU,:], ia[eU,:]
            (
            [1 x]  [     1 -1/2*x]
            [0 2], [     0    1/2]
            )
            sage: a[eV,:], ia[eV,:]
            (
            [ 1/4*u + 1/4*v + 3/2 -1/4*u - 1/4*v - 1/2]
            [ 1/4*u + 1/4*v - 1/2 -1/4*u - 1/4*v + 3/2],
            [-1/8*u - 1/8*v + 3/4  1/8*u + 1/8*v + 1/4]
            [-1/8*u - 1/8*v + 1/4  1/8*u + 1/8*v + 3/4]
            )

        Let us check that ia is indeed the inverse of a::

            sage: s = a.contract(ia)
            sage: s[eU,:], s[eV,:]
            (
            [1 0]  [1 0]
            [0 1], [0 1]
            )
            sage: s = ia.contract(a)
            sage: s[eU,:], s[eV,:]
            (
            [1 0]  [1 0]
            [0 1], [0 1]
            )

        The result is cached::

            sage: a.inverse() is ia
            True

        Instead of ``inverse()``, one can use the power minus one to get the
        inverse::

            sage: ia is a^(-1)
            True

        or the operator ``~``::

            sage: ia is ~a
            True
        """
        if self._is_identity:
            return self
        if self._inverse is None:
            from sage.tensor.modules.format_utilities import is_atomic
            if self._name is None:
                inv_name = None
            else:
                if is_atomic(self._name, ['*']):
                    inv_name = self._name + '^(-1)'
                else:
                    inv_name = '(' + self._name + ')^(-1)'
            if self._latex_name is None:
                inv_latex_name = None
            else:
                if is_atomic(self._latex_name, ['\\circ', '\\otimes']):
                    inv_latex_name = self._latex_name + r'^{-1}'
                else:
                    inv_latex_name = r'\left(' + self._latex_name + \
                                     r'\right)^{-1}'
            self._inverse = self._vmodule.automorphism(name=inv_name,
                                                       latex_name=inv_latex_name)
            for dom, rst in self._restrictions.items():
                self._inverse._restrictions[dom] = rst.inverse()
        return self._inverse

    inverse = __invert__

    def _mul_(self, other):
        r"""
        Automorphism composition.

        This implements the group law of `GL(X(U,\Phi))`, with `X(U,\Phi)`
        being the module of ``self``.

        INPUT:

        - ``other`` -- an automorphism of the same module as ``self``

        OUTPUT:

        - the automorphism resulting from the composition of ``other`` and
          ``self``

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.automorphism_field({e_xy: [[-1, 0], [0, 1]]}, name='a')
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.automorphism_field({e_uv: [[1, 0], [0, -2]]}, name='b')
            sage: b.add_comp_by_continuation(e_xy, U.intersection(V), c_xy)
            sage: s = a._mul_(b); s
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: s.display(e_xy)
            -(x^4 - 10*x^2*y^2 + y^4)/(x^4 + 2*x^2*y^2 + y^4) ∂/∂x⊗dx
             - 6*(x^3*y - x*y^3)/(x^4 + 2*x^2*y^2 + y^4) ∂/∂x⊗dy
             + 6*(x^3*y - x*y^3)/(x^4 + 2*x^2*y^2 + y^4) ∂/∂y⊗dx
             - 2*(x^4 - 4*x^2*y^2 + y^4)/(x^4 + 2*x^2*y^2 + y^4) ∂/∂y⊗dy
            sage: s.display(e_uv)
            -(u^4 - 6*u^2*v^2 + v^4)/(u^4 + 2*u^2*v^2 + v^4) ∂/∂u⊗du
             + 8*(u^3*v - u*v^3)/(u^4 + 2*u^2*v^2 + v^4) ∂/∂u⊗dv
             - 4*(u^3*v - u*v^3)/(u^4 + 2*u^2*v^2 + v^4) ∂/∂v⊗du
            - 2*(u^4 - 6*u^2*v^2 + v^4)/(u^4 + 2*u^2*v^2 + v^4) ∂/∂v⊗dv
            sage: w = M.vector_field(name='w')
            sage: w[e_xy, :] = [3, 1]
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s(w) == a(b(w))  # long time
            True
        """
        # No need for consistency check since self and other are guaranteed
        # to have the same parent. In particular, they are defined on the same
        # module.
        #
        # Special cases:
        if self._is_identity:
            return other
        if other._is_identity:
            return self
        if other is self._inverse or self is other._inverse:
            return self.parent().one()
        # General case:
        resu = type(self)(self._vmodule)
        for dom in self._common_subdomains(other):
            resu._restrictions[dom] = (self._restrictions[dom]
                                       * other._restrictions[dom])
        return resu

    #### End of MultiplicativeGroupElement methods ####

    def __mul__(self, other):
        r"""
        Redefinition of
        :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.__mul__`
        so that ``*`` dispatches either to automorphism composition or
        to the tensor product.

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.automorphism_field(name='a')
            sage: a[e_xy, :] = [[-1, 0], [0, 1]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.automorphism_field(name='b')
            sage: b[e_uv, :] = [[1, 0], [0, -2]]
            sage: b.add_comp_by_continuation(e_xy, U.intersection(V), c_xy)
            sage: w = M.vector_field(name='w')
            sage: w[e_xy, :] = [3, 1]
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a.__mul__(b); s  # automorphism composition
            Field of tangent-space automorphisms on the 2-dimensional differentiable manifold M
            sage: s(w) == a(b(w))  # long time
            True
            sage: s = a.__mul__(w); s  # tensor product
            Tensor field of type (2,1) on the 2-dimensional differentiable manifold M
        """
        if isinstance(other, AutomorphismField):
            return self._mul_(other)  # general linear group law
        else:
            return TensorField.__mul__(self, other)  # tensor product

    def __imul__(self, other):
        r"""
        Redefinition of
        :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.__imul__`

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.automorphism_field(name='a')
            sage: a[e_xy, :] = [[-1, 0], [0, 1]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.automorphism_field(name='b')
            sage: b[e_uv, :] = [[1, 0], [0, -2]]
            sage: b.add_comp_by_continuation(e_xy, U.intersection(V), c_xy)
            sage: a.__imul__(b)
            Field of tangent-space automorphisms on the 2-dimensional differentiable manifold M
            sage: s = a*b
            sage: a *= b
            sage: a == s
            True
        """
        return self.__mul__(other)

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of ``self`` to some subdomain.

        This is a redefinition of
        :meth:`sage.manifolds.differentiable.tensorfield.TensorField.restrict`
        to take into account the identity map.

        INPUT:

        - ``subdomain`` --
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          open subset `V` of ``self._domain``
        - ``dest_map`` -- (default: ``None``)
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`;
          destination map `\Phi:\ V \rightarrow N`, where `N` is a
          subdomain of ``self._codomain``; if ``None``, the restriction
          of ``self.base_module().destination_map()`` to `V` is used

        OUTPUT: a :class:`AutomorphismField` representing the restriction

        EXAMPLES:

        Restrictions of an automorphism field on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of the North pole
            sage: stereoN.<x,y> = U.chart()  # stereographic coordinates from the North pole
            sage: eN = stereoN.frame() # the associated vector frame
            sage: V =  M.open_subset('V') # the complement of the South pole
            sage: stereoS.<u,v> = V.chart()  # stereographic coordinates from the South pole
            sage: eS = stereoS.frame() # the associated vector frame
            sage: transf = stereoN.transition_map(stereoS, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                 intersection_name='W',
            ....:                                 restrictions1= x^2+y^2!=0,
            ....:                                 restrictions2= u^2+v^2!=0)
            sage: inv = transf.inverse() # transformation from stereoS to stereoN
            sage: W = U.intersection(V) # the complement of the North and South poles
            sage: stereoN_W = W.atlas()[0]  # restriction of stereo. coord. from North pole to W
            sage: stereoS_W = W.atlas()[1]  # restriction of stereo. coord. from South pole to W
            sage: eN_W = stereoN_W.frame() ; eS_W = stereoS_W.frame()
            sage: a = M.automorphism_field({eN: [[1, atan(x^2+y^2)], [0,3]]},
            ....:                          name='a')
            sage: a.add_comp_by_continuation(eS, W, chart=stereoS); a
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold S^2
            sage: a.restrict(U)
            Field of tangent-space automorphisms a on the Open subset U of the
             2-dimensional differentiable manifold S^2
            sage: a.restrict(U)[eN,:]
            [                1 arctan(x^2 + y^2)]
            [                0                 3]
            sage: a.restrict(V)
            Field of tangent-space automorphisms a on the Open subset V of the
             2-dimensional differentiable manifold S^2
            sage: a.restrict(V)[eS,:]
            [   (u^4 + 10*u^2*v^2 + v^4 + 2*(u^3*v - u*v^3)*arctan(1/(u^2 + v^2)))/(u^4 + 2*u^2*v^2 + v^4)  -(4*u^3*v - 4*u*v^3 + (u^4 - 2*u^2*v^2 + v^4)*arctan(1/(u^2 + v^2)))/(u^4 + 2*u^2*v^2 + v^4)]
            [                    4*(u^2*v^2*arctan(1/(u^2 + v^2)) - u^3*v + u*v^3)/(u^4 + 2*u^2*v^2 + v^4) (3*u^4 - 2*u^2*v^2 + 3*v^4 - 2*(u^3*v - u*v^3)*arctan(1/(u^2 + v^2)))/(u^4 + 2*u^2*v^2 + v^4)]
            sage: a.restrict(W)
            Field of tangent-space automorphisms a on the Open subset W of the
             2-dimensional differentiable manifold S^2
            sage: a.restrict(W)[eN_W,:]
            [                1 arctan(x^2 + y^2)]
            [                0                 3]

        Restrictions of the field of tangent-space identity maps::

            sage: id = M.tangent_identity_field() ; id
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold S^2
            sage: id.restrict(U)
            Field of tangent-space identity maps on the Open subset U of the
             2-dimensional differentiable manifold S^2
            sage: id.restrict(U)[eN,:]
            [1 0]
            [0 1]
            sage: id.restrict(V)
            Field of tangent-space identity maps on the Open subset V of the
             2-dimensional differentiable manifold S^2
            sage: id.restrict(V)[eS,:]
            [1 0]
            [0 1]
            sage: id.restrict(W)[eN_W,:]
            [1 0]
            [0 1]
            sage: id.restrict(W)[eS_W,:]
            [1 0]
            [0 1]
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if self is not self.parent().one():
                return TensorField.restrict(self, subdomain, dest_map=dest_map)
            # Special case of the immutable identity map:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subset of " +
                                 "the field's domain")
            if dest_map is None:
                dest_map = self._vmodule._dest_map.restrict(subdomain)
            elif not dest_map._codomain.is_subset(self._ambient_domain):
                raise ValueError("the argument 'dest_map' is not compatible " +
                                 "with the ambient domain of " +
                                 "the {}".format(self))
            smodule = subdomain.vector_field_module(dest_map=dest_map)
            self._restrictions[subdomain] = smodule.identity_map()
        return self._restrictions[subdomain]


#******************************************************************************

class AutomorphismFieldParal(FreeModuleAutomorphism, TensorFieldParal):
    r"""
    Field of tangent-space automorphisms with values on a parallelizable
    manifold.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi: U \rightarrow M` to a parallelizable manifold `M`,
    a *field of tangent-space automorphisms along* `U` *with values on*
    `M\supset\Phi(U)` is a differentiable map

    .. MATH::

        a:\ U  \longrightarrow T^{(1,1)}M

    (`T^{(1,1)}M` being the tensor bundle of type `(1,1)` over `M`) such
    that

    .. MATH::

        \forall p \in U,\ a(p) \in \mathrm{Aut}(T_{\Phi(p)} M)

    i.e. `a(p)` is an automorphism of the tangent space to `M` at the point
    `\Phi(p)`.

    The standard case of a field of tangent-space automorphisms *on* a
    manifold corresponds to `U=M` and `\Phi = \mathrm{Id}_M`. Other
    common cases are `\Phi` being an immersion and `\Phi` being a curve in `M`
    (`U` is then an open interval of `\RR`).

    .. NOTE::

        If `M` is not parallelizable, the class :class:`AutomorphismField`
        *must* be used instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `M` via the map `\Phi`
    - ``name`` -- (default: ``None``) name given to the field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the field;
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A `\pi/3`-rotation in the Euclidean 2-plane::

        sage: M = Manifold(2, 'R^2')
        sage: c_xy.<x,y> = M.chart()
        sage: rot = M.automorphism_field([[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]],
        ....:                            name='R'); rot
        Field of tangent-space automorphisms R on the 2-dimensional
         differentiable manifold R^2
        sage: rot.parent()
        General linear group of the Free module X(R^2) of vector fields on the
         2-dimensional differentiable manifold R^2

    The inverse automorphism is obtained via the method :meth:`inverse`::

        sage: inv = rot.inverse() ; inv
        Field of tangent-space automorphisms R^(-1) on the 2-dimensional
         differentiable manifold R^2
        sage: latex(inv)
        R^{-1}
        sage: inv[:]
        [1/2*sqrt(3)         1/2]
        [       -1/2 1/2*sqrt(3)]
        sage: rot[:]
        [1/2*sqrt(3)        -1/2]
        [        1/2 1/2*sqrt(3)]
        sage: inv[:] * rot[:]  # check
        [1 0]
        [0 1]

    Equivalently, one can use the power minus one to get the inverse::

        sage: inv is rot^(-1)
        True

    or the operator ``~``::

        sage: inv is ~rot
        True
    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a field of tangent-space automorphisms.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``AutomorphismFieldParal``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: GL = XM.general_linear_group()
            sage: a = GL.element_class(XM, name='a'); a
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold M
            sage: a[:] = [[1+x^2, x*y], [0, 1+y^2]]
            sage: a.parent()
            General linear group of the Free module X(M) of vector fields on
             the 2-dimensional differentiable manifold M
            sage: a.parent() is M.automorphism_field_group()
            True
            sage: TestSuite(a).run()

        Construction of the field of identity maps::

            sage: b = GL.one(); b
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold M
            sage: b[:]
            [1 0]
            [0 1]
            sage: TestSuite(b).run()
        """
        FreeModuleAutomorphism.__init__(self, vector_field_module,
                                        name=name, latex_name=latex_name)
        # TensorFieldParal attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        self._is_identity = False # a priori
        # Initialization of derived quantities:
        TensorFieldParal._init_derived(self)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.automorphism_field(name='a')
            sage: a._repr_()
            'Field of tangent-space automorphisms a on the 2-dimensional differentiable manifold M'
            sage: repr(a)  # indirect doctest
            'Field of tangent-space automorphisms a on the 2-dimensional differentiable manifold M'
            sage: a  # indirect doctest
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold M
        """
        description = "Field of tangent-space "
        if self is self.parent().one():
            description += "identity maps "
        else:
            description += "automorphisms "
            if self._name is not None:
                description += self._name + " "
        return self._final_repr(description)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- boolean (default: ``True``); determines whether the
          restrictions of ``self`` to subdomains are deleted

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.automorphism_field(name='a')
            sage: a._del_derived()
        """
        # Delete the derived quantities pertaining to the mother classes:
        FreeModuleAutomorphism._del_derived(self)
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)

     # Method _new_instance() is defined in mother class FreeModuleAutomorphism

    def __call__(self, *arg):
        r"""
        Redefinition of
        :meth:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism.__call__`
        to allow for domain treatment.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.automorphism_field([[0, 1], [-1, 0]], name='a')
            sage: v = M.vector_field(-y, x, name='v')
            sage: z = M.one_form(1+y^2, x*y, name='z')
            sage: s = a.__call__(v); s
            Vector field a(v) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(v) = x ∂/∂x + y ∂/∂y
            sage: s = a.__call__(z, v); s
            Scalar field a(z,v) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(z,v): M → ℝ
               (x, y) ↦ 2*x*y^2 + x
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: s = a.__call__(v.restrict(U)); s
            Vector field a(v) on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: s = a.__call__(z.restrict(U), v); s
            Scalar field a(z,v) on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: s.display()
            a(z,v): U → ℝ
               (x, y) ↦ 2*x*y^2 + x
        """
        if len(arg) == 1:
            # the automorphism acting as such (map of a vector field to a
            # vector field)
            vector = arg[0]
            dom = self._domain.intersection(vector._domain)
            return FreeModuleAutomorphism.__call__(self.restrict(dom),
                                                   vector.restrict(dom))
        elif len(arg) == 2:
            # the automorphism acting as a type (1,1) tensor on a pair
            # (1-form, vector field), returning a scalar field:
            oneform = arg[0]
            vector = arg[1]
            dom = self._domain.intersection(oneform._domain).intersection(
                                                                vector._domain)
            return FreeModuleAutomorphism.__call__(self.restrict(dom),
                                                   oneform.restrict(dom),
                                                   vector.restrict(dom))
        else:
            raise TypeError("wrong number of arguments")

    def __invert__(self):
        r"""
        Return the inverse automorphism of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.automorphism_field([[0, 2], [-1, 0]], name='a')
            sage: b = a.inverse(); b
            Field of tangent-space automorphisms a^(-1) on the 2-dimensional
             differentiable manifold M
            sage: b[:]
            [  0  -1]
            [1/2   0]
            sage: a[:]
            [ 0  2]
            [-1  0]

        The result is cached::

            sage: a.inverse() is b
            True

        Instead of ``inverse()``, one can use the power minus one to get the
        inverse::

            sage: b is a^(-1)
            True

        or the operator ``~``::

            sage: b is ~a
            True
        """
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        from sage.matrix.constructor import matrix
        from sage.tensor.modules.comp import Components
        if self._is_identity:
            return self
        if self._inverse is None:
            from sage.tensor.modules.format_utilities import is_atomic
            if self._name is None:
                inv_name = None
            else:
                if is_atomic(self._name, ['*']):
                    inv_name = self._name + '^(-1)'
                else:
                    inv_name = '(' + self._name + ')^(-1)'
            if self._latex_name is None:
                inv_latex_name = None
            else:
                if is_atomic(self._latex_name, ['\\circ', '\\otimes']):
                    inv_latex_name = self._latex_name + r'^{-1}'
                else:
                    inv_latex_name = r'\left(' + self._latex_name + \
                                     r'\right)^{-1}'
            fmodule = self._fmodule
            si = fmodule._sindex
            nsi = fmodule._rank + si
            self._inverse = fmodule.automorphism(name=inv_name,
                                                 latex_name=inv_latex_name)
            for frame in self._components:
                if isinstance(frame, CoordFrame):
                    chart = frame._chart
                else:
                    chart = self._domain._def_chart #!# to be improved
                try:
                    # TODO: do the computation without the 'SR' enforcement
                    mat_self = matrix(
                              [[self.comp(frame)[i, j, chart].expr(method='SR')
                              for j in range(si, nsi)] for i in range(si, nsi)])
                except (KeyError, ValueError):
                    continue
                mat_inv = mat_self.inverse()
                cinv = Components(fmodule._ring, frame, 2, start_index=si,
                                  output_formatter=fmodule._output_formatter)
                for i in range(si, nsi):
                    for j in range(si, nsi):
                        val = chart.simplify(mat_inv[i-si,j-si], method='SR')
                        cinv[i, j] = {chart: val}
                self._inverse._components[frame] = cinv
        return self._inverse

    inverse = __invert__

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of ``self`` to some subset of its domain.

        If such restriction has not been defined yet, it is constructed here.

        This is a redefinition of
        :meth:`sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.restrict`
        to take into account the identity map.

        INPUT:

        - ``subdomain`` --
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
          open subset `V` of ``self._domain``
        - ``dest_map`` -- (default: ``None``)
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          destination map `\Phi:\ V \rightarrow N`, where `N` is a subset of
          ``self._codomain``; if ``None``, the restriction of
          ``self.base_module().destination_map()`` to `V` is used

        OUTPUT: a :class:`AutomorphismFieldParal` representing the restriction

        EXAMPLES:

        Restriction of an automorphism field defined on `\RR^2` to a disk::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: D = M.open_subset('D') # the unit open disc
            sage: c_cart_D = c_cart.restrict(D, x^2+y^2<1)
            sage: a = M.automorphism_field([[1, x*y], [0, 3]], name='a'); a
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold R^2
            sage: a.restrict(D)
            Field of tangent-space automorphisms a on the Open subset D of the
             2-dimensional differentiable manifold R^2
            sage: a.restrict(D)[:]
            [  1 x*y]
            [  0   3]

        Restriction to the disk of the field of tangent-space identity maps::

            sage: id = M.tangent_identity_field() ; id
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold R^2
            sage: id.restrict(D)
            Field of tangent-space identity maps on the Open subset D of the
             2-dimensional differentiable manifold R^2
            sage: id.restrict(D)[:]
            [1 0]
            [0 1]
            sage: id.restrict(D) == D.tangent_identity_field()
            True
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not self._is_identity:
                return TensorFieldParal.restrict(self, subdomain,
                                                 dest_map=dest_map)
            # Special case of the identity map:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subset of " +
                                 "the field's domain.")
            if dest_map is None:
                dest_map = self._fmodule._dest_map.restrict(subdomain)
            elif not dest_map._codomain.is_subset(self._ambient_domain):
                raise ValueError("the argument 'dest_map' is not compatible " +
                                 "with the ambient domain of " +
                                 "the {}".format(self))
            smodule = subdomain.vector_field_module(dest_map=dest_map)
            self._restrictions[subdomain] = smodule.identity_map()
        return self._restrictions[subdomain]

    def at(self, point):
        r"""
        Value of ``self`` at a given point.

        If the current field of tangent-space automorphisms is

        .. MATH::

            a:\ U \longrightarrow T^{(1,1)} M

        associated with the differentiable map

        .. MATH::

            \Phi:\ U \longrightarrow M,

        where `U` and `M` are two manifolds (possibly `U = M` and
        `\Phi = \mathrm{Id}_M`), then for any point `p \in U`,
        `a(p)` is an automorphism of the tangent space `T_{\Phi(p)}M`.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` in the domain of the field of automorphisms `a`

        OUTPUT:

        - the automorphism `a(p)` of the tangent vector space `T_{\Phi(p)}M`

        EXAMPLES:

        Automorphism at some point of a tangent space of a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: a = M.automorphism_field([[1+exp(y), x*y], [0, 1+x^2]],
            ....:                          name='a')
            sage: a.display()
            a = (e^y + 1) ∂/∂x⊗dx + x*y ∂/∂x⊗dy + (x^2 + 1) ∂/∂y⊗dy
            sage: p = M.point((-2,3), name='p') ; p
            Point p on the 2-dimensional differentiable manifold M
            sage: ap = a.at(p) ; ap
            Automorphism a of the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: ap.display()
            a = (e^3 + 1) ∂/∂x⊗dx - 6 ∂/∂x⊗dy + 5 ∂/∂y⊗dy
            sage: ap.parent()
            General linear group of the Tangent space at Point p on the
             2-dimensional differentiable manifold M

        The identity map of the tangent space at point ``p``::

            sage: id = M.tangent_identity_field() ; id
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold M
            sage: idp = id.at(p) ; idp
            Identity map of the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: idp is M.tangent_space(p).identity_map()
            True
            sage: idp.display()
            Id = ∂/∂x⊗dx + ∂/∂y⊗dy
            sage: idp.parent()
            General linear group of the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: idp * ap == ap
            True
        """
        if point not in self._domain:
            raise TypeError("the {} is not in the domain of the {}".format(
                                                                  point, self))
        dest_map = self._fmodule._dest_map
        if dest_map.is_identity():
            amb_point = point
        else:
            amb_point = dest_map(point)  #  "ambient" point
        ts = amb_point._manifold.tangent_space(amb_point)
        if self._is_identity:
            return ts.identity_map()
        resu = ts.automorphism(name=self._name, latex_name=self._latex_name)
        for frame, comp in self._components.items():
            comp_resu = resu.add_comp(frame.at(point))
            for ind, val in comp._comp.items():
                comp_resu._comp[ind] = val(point)
        return resu
