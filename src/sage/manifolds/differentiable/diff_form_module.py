r"""
Differential Form Modules

The set `\Omega^p(U, \Phi)` of `p`-forms along a differentiable manifold `U`
with values on a differentiable manifold `M` via a differentiable map
`\Phi:\ U \rightarrow M` (possibly `U = M` and `\Phi = \mathrm{Id}_M`)
is a module over the algebra `C^k(U)` of differentiable scalar fields on `U`.
It is a free module if and only if `M` is parallelizable. Accordingly,
two classes implement `\Omega^p(U, \Phi)`:

- :class:`DiffFormModule` for differential forms with values on a generic
  (in practice, not parallelizable) differentiable manifold `M`
- :class:`DiffFormFreeModule` for differential forms with values on a
  parallelizable manifold `M`

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [KN1963]_
- [Lee2013]_

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

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.modules import Modules
from sage.tensor.modules.ext_pow_free_module import ExtPowerDualFreeModule
from sage.manifolds.differentiable.diff_form import DiffForm, DiffFormParal
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal


class DiffFormModule(UniqueRepresentation, Parent):
    r"""
    Module of differential forms of a given degree `p` (`p`-forms) along a
    differentiable manifold `U` with values on a differentiable manifold `M`.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi: U \rightarrow M` to a differentiable manifold `M`, the set
    `\Omega^p(U, \Phi)` of `p`-forms along `U` with values on `M` is
    a module over `C^k(U)`, the commutative algebra of differentiable
    scalar fields on `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).
    The standard case of `p`-forms *on* a differentiable manifold `M`
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases
    are `\Phi` being an immersion and `\Phi` being a curve in `M`
    (`U` is then an open interval of `\RR`).

    .. NOTE::

        This class implements `\Omega^p(U,\Phi)` in the case where `M` is
        not assumed to be parallelizable; the module `\Omega^p(U, \Phi)`
        is then not necessarily free. If `M` is parallelizable, the class
        :class:`DiffFormFreeModule` must be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U, \Phi)` of vector
      fields along `U` with values on `M` via the map `\Phi: U \rightarrow M`
    - ``degree`` -- positive integer; the degree `p` of the differential forms

    EXAMPLES:

    Module of 2-forms on a non-parallelizable 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
        ....:  intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: XM = M.vector_field_module() ; XM
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M
        sage: A = M.diff_form_module(2) ; A
        Module Omega^2(M) of 2-forms on the 2-dimensional differentiable
         manifold M
        sage: latex(A)
        \Omega^{2}\left(M\right)

    ``A`` is nothing but the second exterior power of the dual of ``XM``, i.e.
    we have `\Omega^{2}(M) = \Lambda^2(\mathfrak{X}(M)^*)`::

        sage: A is XM.dual_exterior_power(2)
        True

    Modules of differential forms are unique::

        sage: A is M.diff_form_module(2)
        True

    `\Omega^2(M)` is a module over the algebra `C^k(M)` of (differentiable)
    scalar fields on `M`::

        sage: A.category()
        Category of modules over Algebra of differentiable scalar fields on
         the 2-dimensional differentiable manifold M
        sage: CM = M.scalar_field_algebra() ; CM
        Algebra of differentiable scalar fields on the 2-dimensional
         differentiable manifold M
        sage: A in Modules(CM)
        True
        sage: A.base_ring() is CM
        True
        sage: A.base_module()
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M
        sage: A.base_module() is XM
        True

    Elements can be constructed from ``A()``. In particular, ``0`` yields
    the zero element of ``A``::

        sage: z = A(0) ; z
        2-form zero on the 2-dimensional differentiable manifold M
        sage: z.display(eU)
        zero = 0
        sage: z.display(eV)
        zero = 0
        sage: z is A.zero()
        True

    while non-zero elements are constructed by providing their components in a
    given vector frame::

        sage: a = A([[0,3*x],[-3*x,0]], frame=eU, name='a') ; a
        2-form a on the 2-dimensional differentiable manifold M
        sage: a.add_comp_by_continuation(eV, W, c_uv) # finishes initializ. of a
        sage: a.display(eU)
        a = 3*x dx∧dy
        sage: a.display(eV)
        a = (-3/4*u - 3/4*v) du∧dv

    An alternative is to construct the 2-form from an empty list of
    components and to set the nonzero nonredundant components afterwards::

        sage: a = A([], name='a')
        sage: a[eU,0,1] = 3*x
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = 3*x dx∧dy
        sage: a.display(eV)
        a = (-3/4*u - 3/4*v) du∧dv

    The module `\Omega^1(M)` is nothing but the dual of `\mathfrak{X}(M)`
    (the module of vector fields on `M`)::

        sage: L1 = M.diff_form_module(1) ; L1
        Module Omega^1(M) of 1-forms on the 2-dimensional differentiable
         manifold M
        sage: L1 is XM.dual()
        True

    Since any tensor field of type `(0,1)` is a 1-form, there is a coercion
    map from the set `T^{(0,1)}(M)` of such tensors to `\Omega^1(M)`::

        sage: T01 = M.tensor_field_module((0,1)) ; T01
        Module T^(0,1)(M) of type-(0,1) tensors fields on the 2-dimensional
         differentiable manifold M
        sage: L1.has_coerce_map_from(T01)
        True

    There is also a coercion map in the reverse direction::

        sage: T01.has_coerce_map_from(L1)
        True

    For a degree `p \geq 2`, the coercion holds only in the direction
    `\Omega^p(M)\rightarrow T^{(0,p)}(M)`::

        sage: T02 = M.tensor_field_module((0,2)) ; T02
        Module T^(0,2)(M) of type-(0,2) tensors fields on the 2-dimensional
         differentiable manifold M
        sage: T02.has_coerce_map_from(A)
        True
        sage: A.has_coerce_map_from(T02)
        False

    The coercion map `T^{(0,1)}(M) \rightarrow \Omega^1(M)` in action::

        sage: b = T01([y,x], frame=eU, name='b') ; b
        Tensor field b of type (0,1) on the 2-dimensional differentiable
         manifold M
        sage: b.add_comp_by_continuation(eV, W, c_uv)
        sage: b.display(eU)
        b = y dx + x dy
        sage: b.display(eV)
        b = 1/2*u du - 1/2*v dv
        sage: lb = L1(b) ; lb
        1-form b on the 2-dimensional differentiable manifold M
        sage: lb.display(eU)
        b = y dx + x dy
        sage: lb.display(eV)
        b = 1/2*u du - 1/2*v dv

    The coercion map `\Omega^1(M) \rightarrow T^{(0,1)}(M)` in action::

        sage: tlb = T01(lb) ; tlb
        Tensor field b of type (0,1) on the 2-dimensional differentiable
         manifold M
        sage: tlb.display(eU)
        b = y dx + x dy
        sage: tlb.display(eV)
        b = 1/2*u du - 1/2*v dv
        sage: tlb == b
        True

    The coercion map `\Omega^2(M) \rightarrow T^{(0,2)}(M)` in action::

        sage: ta = T02(a) ; ta
        Tensor field a of type (0,2) on the 2-dimensional differentiable
         manifold M
        sage: ta.display(eU)
        a = 3*x dx⊗dy - 3*x dy⊗dx
        sage: a.display(eU)
        a = 3*x dx∧dy
        sage: ta.display(eV)
        a = (-3/4*u - 3/4*v) du⊗dv + (3/4*u + 3/4*v) dv⊗du
        sage: a.display(eV)
        a = (-3/4*u - 3/4*v) du∧dv

    There is also coercion to subdomains, which is nothing but the restriction
    of the differential form to some subset of its domain::

        sage: L2U = U.diff_form_module(2) ; L2U
        Free module Omega^2(U) of 2-forms on the Open subset U of the
         2-dimensional differentiable manifold M
        sage: L2U.has_coerce_map_from(A)
        True
        sage: a_U = L2U(a) ; a_U
        2-form a on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: a_U.display(eU)
        a = 3*x dx∧dy

    """
    Element = DiffForm

    def __init__(self, vector_field_module, degree):
        r"""
        Construction a module of differential forms.

        TESTS:

        Module of 2-forms on a non-parallelizable 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: from sage.manifolds.differentiable.diff_form_module import \
            ....:                                                DiffFormModule
            sage: A = DiffFormModule(M.vector_field_module(), 2) ; A
            Module Omega^2(M) of 2-forms on the 2-dimensional differentiable
             manifold M
            sage: TestSuite(A).run(skip='_test_elements')

        In the above test suite, ``_test_elements`` is skipped because of the
        ``_test_pickling`` error of the elements (to be fixed in
        :class:`sage.manifolds.differentiable.tensorfield.TensorField`)

        """
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "Omega^{}(".format(degree) + domain._name
        latex_name = r"\Omega^{{{}}}\left({}".format(degree, domain._latex_name)
        if dest_map is not domain.identity_map():
            dm_name = dest_map._name
            dm_latex_name = dest_map._latex_name
            if dm_name is None:
                dm_name = "unnamed map"
            if dm_latex_name is None:
                dm_latex_name = r"\mathrm{unnamed\; map}"
            name += "," + dm_name
            latex_name += "," + dm_latex_name
        self._name = name + ")"
        self._latex_name = latex_name + r"\right)"
        self._vmodule = vector_field_module
        self._degree = degree
        # the member self._ring is created for efficiency (to avoid calls to
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring, category=Modules(self._ring))
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: self._zero_element is not constructed here, since no element
        # can be constructed here, to avoid some infinite recursion.

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a differential form.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: A = M.diff_form_module(2)
            sage: a = A([[0, x*y], [-x*y, 0]], name='a'); a
            2-form a on the 2-dimensional differentiable manifold M
            sage: a.display(c_xy.frame())
            a = x*y dx∧dy
            sage: A(0) is A.zero()
            True

        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, (DiffForm, DiffFormParal)):
            # coercion by domain restriction
            if (self._degree == comp._tensor_type[1]
                   and self._domain.is_subset(comp._domain)
                   and self._ambient_domain.is_subset(comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise TypeError("cannot convert the {} ".format(comp) +
                                "to an element of {}".format(self))
        if isinstance(comp, TensorField):
            # coercion of a tensor of type (0,1) to a linear form
            tensor = comp # for readability
            if (tensor.tensor_type() == (0,1) and self._degree == 1
                     and tensor._vmodule is self._vmodule):
                resu = self.element_class(self._vmodule, 1, name=tensor._name,
                                          latex_name=tensor._latex_name)
                for dom, rst in tensor._restrictions.items():
                    resu._restrictions[dom] = dom.diff_form_module(1)(rst)
                return resu
            else:
                raise TypeError("cannot convert the {} ".format(tensor) +
                                "to an element of {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._vmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) differential form.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: A = M.diff_form_module(2)
            sage: A._an_element_()
            2-form on the 2-dimensional differentiable manifold M

        """
        resu = self.element_class(self._vmodule, self._degree)
        for oc in self._domain.open_covers(trivial=False):
            # the first non-trivial open cover is selected
            for dom in oc:
                vmodule_dom = dom.vector_field_module(
                                         dest_map=self._dest_map.restrict(dom))
                dmodule_dom = vmodule_dom.dual_exterior_power(self._degree)
                resu.set_restriction(dmodule_dom._an_element_())
            return resu
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: A1 = M.diff_form_module(1)
            sage: A1._coerce_map_from_(M.tensor_field_module((0,1)))
            True
            sage: A2 = M.diff_form_module(2)
            sage: A2._coerce_map_from_(M.tensor_field_module((0,2)))
            False
            sage: U = M.open_subset('U')
            sage: A2U = U.diff_form_module(2)
            sage: A2U._coerce_map_from_(A2)
            True
            sage: A2._coerce_map_from_(A2U)
            False

        """
        if isinstance(other, (DiffFormModule, DiffFormFreeModule)):
            # coercion by domain restriction
            return (self._degree == other._degree
                    and self._domain.is_subset(other._domain)
                    and self._ambient_domain.is_subset(other._ambient_domain))

        from sage.manifolds.differentiable.tensorfield_module import TensorFieldModule
        if isinstance(other, TensorFieldModule):
            # coercion of a type-(0,1) tensor to a linear form
            return (self._vmodule is other._vmodule and self._degree == 1
                    and other.tensor_type() == (0,1))

        return False

    @cached_method
    def zero(self):
        """
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: A2 = M.diff_form_module(2)
            sage: A2.zero()
            2-form zero on the 3-dimensional differentiable manifold M

        """
        zero = self._element_constructor_(name='zero', latex_name='0')
        for frame in self._domain._frames:
            if self._dest_map.restrict(frame._domain) == frame._dest_map:
                zero.add_comp(frame)
                # (since new components are initialized to zero)
        zero._is_zero = True  # This element is certainly zero
        zero.set_immutable()
        return zero

    #### End of Parent methods

    def _repr_(self):
        r"""
        Return a string representation of the object.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: A2 = M.diff_form_module(2)
            sage: A2
            Module Omega^2(M) of 2-forms on
             the 3-dimensional differentiable manifold M

        """
        description = "Module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-forms ".format(self._degree)
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {}".format(
                                            self._domain, self._ambient_domain)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(3, 'M', latex_name=r'\mathcal{M}')
            sage: A2 = M.diff_form_module(2)
            sage: A2._latex_()
            '\\Omega^{2}\\left(\\mathcal{M}\\right)'
            sage: latex(A2)  # indirect doctest
            \Omega^{2}\left(\mathcal{M}\right)

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def base_module(self):
        r"""
        Return the vector field module on which the differential form module
        ``self`` is constructed.

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`
          representing the module on which ``self`` is defined

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: A2 = M.diff_form_module(2) ; A2
            Module Omega^2(M) of 2-forms on the 3-dimensional differentiable
             manifold M
            sage: A2.base_module()
            Module X(M) of vector fields on the 3-dimensional differentiable
             manifold M
            sage: A2.base_module() is M.vector_field_module()
            True
            sage: U = M.open_subset('U')
            sage: A2U = U.diff_form_module(2) ; A2U
            Module Omega^2(U) of 2-forms on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: A2U.base_module()
            Module X(U) of vector fields on the Open subset U of the
             3-dimensional differentiable manifold M

        """
        return self._vmodule

    def degree(self):
        r"""
        Return the degree of the differential forms in ``self``.

        OUTPUT:

        - integer `p` such that ``self`` is a set of `p`-forms

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: M.diff_form_module(1).degree()
            1
            sage: M.diff_form_module(2).degree()
            2
            sage: M.diff_form_module(3).degree()
            3

        """
        return self._degree

# *****************************************************************************

class DiffFormFreeModule(ExtPowerDualFreeModule):
    r"""
    Free module of differential forms of a given degree `p` (`p`-forms) along
    a differentiable manifold `U` with values on a parallelizable manifold `M`.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi:\; U \rightarrow M` to a parallelizable manifold `M` of dimension
    `n`, the set `\Omega^p(U, \Phi)` of `p`-forms along `U` with values on `M`
    is a free module of rank `\binom{n}{p}` over `C^k(U)`, the commutative
    algebra of differentiable scalar fields on `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).
    The standard case of `p`-forms *on* a differentiable manifold `M`
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are
    `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then an
    open interval of `\RR`).

    .. NOTE::

        This class implements `\Omega^p(U, \Phi)` in the case where `M` is
        parallelizable; `\Omega^p(U, \Phi)` is then a *free* module. If `M`
        is not parallelizable, the class :class:`DiffFormModule` must be used
        instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` associated with the map `\Phi: U \rightarrow V`
    - ``degree`` -- positive integer; the degree `p` of the differential forms

    EXAMPLES:

    Free module of 2-forms on a parallelizable 3-dimensional manifold::

        sage: M = Manifold(3, 'M')
        sage: X.<x,y,z> = M.chart()
        sage: XM = M.vector_field_module() ; XM
        Free module X(M) of vector fields on the 3-dimensional differentiable
         manifold M
        sage: A = M.diff_form_module(2) ; A
        Free module Omega^2(M) of 2-forms on the 3-dimensional differentiable
         manifold M
        sage: latex(A)
        \Omega^{2}\left(M\right)

    ``A`` is nothing but the second exterior power of the dual of ``XM``, i.e.
    we have `\Omega^{2}(M) = \Lambda^2(\mathfrak{X}(M)^*)` (see
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerDualFreeModule`)::

        sage: A is XM.dual_exterior_power(2)
        True

    `\Omega^{2}(M)` is a module over the algebra `C^k(M)` of (differentiable)
    scalar fields on `M`::

        sage: A.category()
        Category of finite dimensional modules over Algebra of differentiable
         scalar fields on the 3-dimensional differentiable manifold M
        sage: CM = M.scalar_field_algebra() ; CM
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M
        sage: A in Modules(CM)
        True
        sage: A.base_ring()
        Algebra of differentiable scalar fields on
         the 3-dimensional differentiable manifold M
        sage: A.base_module()
        Free module X(M) of vector fields on
         the 3-dimensional differentiable manifold M
        sage: A.base_module() is XM
        True
        sage: A.rank()
        3

    Elements can be constructed from `A`. In particular, ``0`` yields
    the zero element of `A`::

        sage: A(0)
        2-form zero on the 3-dimensional differentiable manifold M
        sage: A(0) is A.zero()
        True

    while non-zero elements are constructed by providing their components
    in a given vector frame::

        sage: comp = [[0,3*x,-z],[-3*x,0,4],[z,-4,0]]
        sage: a = A(comp, frame=X.frame(), name='a') ; a
        2-form a on the 3-dimensional differentiable manifold M
        sage: a.display()
        a = 3*x dx∧dy - z dx∧dz + 4 dy∧dz

    An alternative is to construct the 2-form from an empty list of
    components and to set the nonzero nonredundant components afterwards::

        sage: a = A([], name='a')
        sage: a[0,1] = 3*x  # component in the manifold's default frame
        sage: a[0,2] = -z
        sage: a[1,2] = 4
        sage: a.display()
        a = 3*x dx∧dy - z dx∧dz + 4 dy∧dz

    The module `\Omega^1(M)` is nothing but the dual of `\mathfrak{X}(M)`
    (the free module of vector fields on `M`)::

        sage: L1 = M.diff_form_module(1) ; L1
        Free module Omega^1(M) of 1-forms on the 3-dimensional differentiable
         manifold M
        sage: L1 is XM.dual()
        True

    Since any tensor field of type `(0,1)` is a 1-form, there is a coercion
    map from the set `T^{(0,1)}(M)` of such tensors to `\Omega^1(M)`::

        sage: T01 = M.tensor_field_module((0,1)) ; T01
        Free module T^(0,1)(M) of type-(0,1) tensors fields on the
         3-dimensional differentiable manifold M
        sage: L1.has_coerce_map_from(T01)
        True

    There is also a coercion map in the reverse direction::

        sage: T01.has_coerce_map_from(L1)
        True

    For a degree `p \geq 2`, the coercion holds only in the direction
    `\Omega^p(M) \rightarrow T^{(0,p)}(M)`::

        sage: T02 = M.tensor_field_module((0,2)); T02
        Free module T^(0,2)(M) of type-(0,2) tensors fields on the
         3-dimensional differentiable manifold M
        sage: T02.has_coerce_map_from(A)
        True
        sage: A.has_coerce_map_from(T02)
        False

    The coercion map `T^{(0,1)}(M) \rightarrow \Omega^1(M)` in action::

        sage: b = T01([-x,2,3*y], name='b'); b
        Tensor field b of type (0,1) on the 3-dimensional differentiable
         manifold M
        sage: b.display()
        b = -x dx + 2 dy + 3*y dz
        sage: lb = L1(b) ; lb
        1-form b on the 3-dimensional differentiable manifold M
        sage: lb.display()
        b = -x dx + 2 dy + 3*y dz

    The coercion map `\Omega^1(M) \rightarrow T^{(0,1)}(M)` in action::

        sage: tlb = T01(lb); tlb
        Tensor field b of type (0,1) on
         the 3-dimensional differentiable manifold M
        sage: tlb == b
        True

    The coercion map `\Omega^2(M) \rightarrow T^{(0,2)}(M)` in action::

        sage: T02 = M.tensor_field_module((0,2)) ; T02
        Free module T^(0,2)(M) of type-(0,2) tensors fields on the
         3-dimensional differentiable manifold M
        sage: ta = T02(a) ; ta
        Tensor field a of type (0,2) on the 3-dimensional differentiable
         manifold M
        sage: ta.display()
        a = 3*x dx⊗dy - z dx⊗dz - 3*x dy⊗dx + 4 dy⊗dz + z dz⊗dx - 4 dz⊗dy
        sage: a.display()
        a = 3*x dx∧dy - z dx∧dz + 4 dy∧dz
        sage: ta.symmetries()  # the antisymmetry is preserved
        no symmetry;  antisymmetry: (0, 1)

    There is also coercion to subdomains, which is nothing but the
    restriction of the differential form to some subset of its domain::

        sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1})
        sage: B = U.diff_form_module(2) ; B
        Free module Omega^2(U) of 2-forms on the Open subset U of the
         3-dimensional differentiable manifold M
        sage: B.has_coerce_map_from(A)
        True
        sage: a_U = B(a) ; a_U
        2-form a on the Open subset U of the 3-dimensional differentiable
         manifold M
        sage: a_U.display()
        a = 3*x dx∧dy - z dx∧dz + 4 dy∧dz

    """

    Element = DiffFormParal

    def __init__(self, vector_field_module, degree):
        r"""
        Construct a free module of differential forms.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: from sage.manifolds.differentiable.diff_form_module import DiffFormFreeModule
            sage: A = DiffFormFreeModule(M.vector_field_module(), 2) ; A
            Free module Omega^2(M) of 2-forms on
             the 3-dimensional differentiable manifold M
            sage: TestSuite(A).run()

        """
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "Omega^{}(".format(degree) + domain._name
        latex_name = r"\Omega^{{{}}}\left({}".format(degree, domain._latex_name)
        if dest_map is not domain.identity_map():
            dm_name = dest_map._name
            dm_latex_name = dest_map._latex_name
            if dm_name is None:
                dm_name = "unnamed map"
            if dm_latex_name is None:
                dm_latex_name = r"\mathrm{unnamed\; map}"
            name += "," + dm_name
            latex_name += "," + dm_latex_name
        name += ")"
        latex_name += r"\right)"
        ExtPowerDualFreeModule.__init__(self, vector_field_module, degree,
                                        name=name, latex_name=latex_name)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a differential form.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: A = M.diff_form_module(2)
            sage: a = A([[0, x], [-x, 0]], name='a'); a
            2-form a on the 2-dimensional differentiable manifold M
            sage: a.display()
            a = x dx∧dy
            sage: A(0) is A.zero()
            True

        Check that :trac:`27658` is fixed::

            sage: f = M.scalar_field(x)
            sage: f in A
            False

        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, (DiffForm, DiffFormParal)):
            # coercion by domain restriction
            if (self._degree == comp._tensor_type[1]
                    and self._domain.is_subset(comp._domain)
                    and self._ambient_domain.is_subset(comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise TypeError("cannot convert the {} ".format(comp) +
                                "to a differential form in {}".format(self))
        if isinstance(comp, TensorFieldParal):
            # coercion of a tensor of type (0,1) to a linear form
            tensor = comp # for readability
            if (tensor.tensor_type() == (0,1) and self._degree == 1
                     and tensor._fmodule is self._fmodule):
                resu = self.element_class(self._fmodule, 1, name=tensor._name,
                                          latex_name=tensor._latex_name)
                for frame, comp in tensor._components.items():
                    resu._components[frame] = comp.copy()
                return resu
            else:
                raise TypeError("cannot convert the {} ".format(tensor) +
                                "to an element of {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._fmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    # Rem: _an_element_ is declared in the superclass ExtPowerDualFreeModule

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: A2 = M.diff_form_module(2)
            sage: U = M.open_subset('U', coord_def = {X: z<0})
            sage: A2U = U.diff_form_module(2)
            sage: A2U._coerce_map_from_(A2)
            True
            sage: A2._coerce_map_from_(A2U)
            False
            sage: A1 = M.diff_form_module(1)
            sage: A2U._coerce_map_from_(A1)
            False
            sage: A1._coerce_map_from_(M.tensor_field_module((0,1)))
            True
            sage: A1._coerce_map_from_(M.tensor_field_module((1,0)))
            False

        """
        if isinstance(other, (DiffFormModule, DiffFormFreeModule)):
            # coercion by domain restriction
            return (self._degree == other._degree
                    and self._domain.is_subset(other._domain)
                    and self._ambient_domain.is_subset(other._ambient_domain))

        from sage.manifolds.differentiable.tensorfield_module import TensorFieldFreeModule
        if isinstance(other, TensorFieldFreeModule):
            # coercion of a type-(0,1) tensor to a linear form
            return (self._fmodule is other._fmodule and self._degree == 1
                    and other.tensor_type() == (0,1))
        return False

    #### End of Parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: A = M.diff_form_module(2)
            sage: A
            Free module Omega^2(M) of 2-forms on
             the 3-dimensional differentiable manifold M

        """
        description = "Free module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-forms ".format(self._degree)
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {}".format(
                                            self._domain, self._ambient_domain)
        return description

