r"""
Multivector Field Modules

The set `A^p(U, \Phi)` of `p`-vector fields along a differentiable
manifold `U` with values on a differentiable manifold `M` via a
differentiable map `\Phi:\ U \rightarrow M` (possibly `U = M` and
`\Phi = \mathrm{Id}_M`) is a module over the algebra `C^k(U)` of
differentiable scalar fields on `U`. It is a free module if and only if
`M` is parallelizable. Accordingly, two classes implement
`A^p(U,\Phi)`:

- :class:`MultivectorModule` for `p`-vector fields with values on a
  generic (in practice, not parallelizable) differentiable manifold `M`
- :class:`MultivectorFreeModule` for `p`-vector fields with values on a
  parallelizable manifold `M`

AUTHORS:

- Eric Gourgoulhon (2017): initial version

REFERENCES:

- \R. L. Bishop and S. L. Goldberg (1980) [BG1980]_
- \C.-M. Marle (1997) [Mar1997]_
"""
#******************************************************************************
#       Copyright (C) 2017 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.modules import Modules
from sage.manifolds.differentiable.multivectorfield import (
    MultivectorField,
    MultivectorFieldParal,
)
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule


class MultivectorModule(UniqueRepresentation, Parent):
    r"""
    Module of multivector fields of a given degree `p` (`p`-vector
    fields) along a differentiable manifold `U` with values on a
    differentiable manifold `M`.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi: U \rightarrow M` to a differentiable manifold `M`, the set
    `A^p(U, \Phi)` of `p`-vector fields (i.e. alternating tensor fields
    of type `(p,0)`) along `U` with values on `M` is a module over
    `C^k(U)`, the commutative algebra of differentiable scalar fields on
    `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).
    The standard case of `p`-vector fields *on* a differentiable
    manifold `M` corresponds to `U = M` and `\Phi = \mathrm{Id}_M`.
    Other common cases are `\Phi` being an immersion and `\Phi` being a
    curve in `M` (`U` is then an open interval of `\RR`).

    .. NOTE::

        This class implements `A^p(U,\Phi)` in the case where `M` is
        not assumed to be parallelizable; the module `A^p(U, \Phi)`
        is then not necessarily free. If `M` is parallelizable, the
        class :class:`MultivectorFreeModule` must be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U, \Phi)` of vector
      fields along `U` with values on `M` via the map
      `\Phi: U \rightarrow M`
    - ``degree`` -- positive integer; the degree `p` of the multivector
      fields

    EXAMPLES:

    Module of 2-vector fields on a non-parallelizable 2-dimensional
    manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
        ....:  intersection_name='W', restrictions1= x>0,
        ....:  restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: XM = M.vector_field_module() ; XM
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M
        sage: A = M.multivector_module(2) ; A
        Module A^2(M) of 2-vector fields on the 2-dimensional
         differentiable manifold M
        sage: latex(A)
        A^{2}\left(M\right)

    ``A`` is nothing but the second exterior power of ``XM``, i.e.
    we have `A^{2}(M) = \Lambda^2(\mathfrak{X}(M))`::

        sage: A is XM.exterior_power(2)
        True

    Modules of multivector fields are unique::

        sage: A is M.multivector_module(2)
        True

    `A^2(M)` is a module over the algebra `C^k(M)` of (differentiable)
    scalar fields on `M`::

        sage: A.category()
        Category of modules over Algebra of differentiable scalar fields
         on the 2-dimensional differentiable manifold M
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

    Elements can be constructed from ``A()``. In particular, ``0``
    yields the zero element of ``A``::

        sage: z = A(0) ; z
        2-vector field zero on the 2-dimensional differentiable
         manifold M
        sage: z.display(eU)
        zero = 0
        sage: z.display(eV)
        zero = 0
        sage: z is A.zero()
        True

    while nonzero elements are constructed by providing their
    components in a given vector frame::

        sage: a = A([[0,3*x],[-3*x,0]], frame=eU, name='a') ; a
        2-vector field a on the 2-dimensional differentiable manifold M
        sage: a.add_comp_by_continuation(eV, W, c_uv) # finishes initializ. of a
        sage: a.display(eU)
        a = 3*x ∂/∂x∧∂/∂y
        sage: a.display(eV)
        a = (-3*u - 3*v) ∂/∂u∧∂/∂v

    An alternative is to construct the 2-vector field from an empty list
    of components and to set the nonzero nonredundant components
    afterwards::

        sage: a = A([], name='a')
        sage: a[eU,0,1] = 3*x
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = 3*x ∂/∂x∧∂/∂y
        sage: a.display(eV)
        a = (-3*u - 3*v) ∂/∂u∧∂/∂v

    The module `A^1(M)` is nothing but the dual of `\mathfrak{X}(M)`
    (the module of vector fields on `M`)::

        sage: A1 = M.multivector_module(1) ; A1
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M
        sage: A1 is XM
        True

    There is a coercion map `A^p(M)\rightarrow T^{(p,0)}(M)`::

        sage: T20 = M.tensor_field_module((2,0)) ; T20
        Module T^(2,0)(M) of type-(2,0) tensors fields on the
         2-dimensional differentiable manifold M
        sage: T20.has_coerce_map_from(A)
        True

    but of course not in the reverse direction, since not all contravariant
    tensor field is alternating::

        sage: A.has_coerce_map_from(T20)
        False

    The coercion map `A^2(M) \rightarrow T^{(2,0)}(M)` in action::

        sage: ta = T20(a) ; ta
        Tensor field a of type (2,0) on the 2-dimensional differentiable
         manifold M
        sage: ta.display(eU)
        a = 3*x ∂/∂x⊗∂/∂y - 3*x ∂/∂y⊗∂/∂x
        sage: a.display(eU)
        a = 3*x ∂/∂x∧∂/∂y
        sage: ta.display(eV)
        a = (-3*u - 3*v) ∂/∂u⊗∂/∂v + (3*u + 3*v) ∂/∂v⊗∂/∂u
        sage: a.display(eV)
        a = (-3*u - 3*v) ∂/∂u∧∂/∂v

    There is also coercion to subdomains, which is nothing but the
    restriction of the multivector field to some subset of its domain::

        sage: A2U = U.multivector_module(2) ; A2U
        Free module A^2(U) of 2-vector fields on the Open subset U of
         the 2-dimensional differentiable manifold M
        sage: A2U.has_coerce_map_from(A)
        True
        sage: a_U = A2U(a) ; a_U
        2-vector field a on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: a_U.display(eU)
        a = 3*x ∂/∂x∧∂/∂y
    """
    Element = MultivectorField

    def __init__(self, vector_field_module, degree):
        r"""
        Construction a module of multivector fields.

        TESTS:

        Module of 2-vector fields on a non-parallelizable 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:             intersection_name='W', restrictions1= x>0,
            ....:             restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: from sage.manifolds.differentiable.multivector_module import \
            ....:                                      MultivectorModule
            sage: A = MultivectorModule(M.vector_field_module(), 2) ; A
            Module A^2(M) of 2-vector fields on the 2-dimensional
             differentiable manifold M
            sage: TestSuite(A).run(skip='_test_elements')

        In the above test suite, ``_test_elements`` is skipped because
        of the ``_test_pickling`` error of the elements (to be fixed in
        :class:`sage.manifolds.differentiable.tensorfield.TensorField`)
        """
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "A^{}(".format(degree) + domain._name
        latex_name = r"A^{{{}}}\left({}".format(degree,
                                                domain._latex_name)
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
        # the member self._ring is created for efficiency (to avoid
        # calls to self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring,
                        category=Modules(self._ring))
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: self._zero_element is not constructed here, since no
        # element can be constructed here, to avoid some infinite
        # recursion.

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a multivector field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: A = M.multivector_module(2)
            sage: a = A([[0, x*y], [-x*y, 0]], name='a'); a
            2-vector field a on the 2-dimensional differentiable
             manifold M
            sage: a.display(c_xy.frame())
            a = x*y ∂/∂x∧∂/∂y
            sage: A(0) is A.zero()
            True
        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, (MultivectorField, MultivectorFieldParal)):
            # coercion by domain restriction
            if (self._degree == comp._tensor_type[0]
                   and self._domain.is_subset(comp._domain)
                   and self._ambient_domain.is_subset(
                                                 comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise TypeError("cannot convert the {} ".format(comp) +
                                "to an element of {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._vmodule, self._degree,
                                  name=name, latex_name=latex_name)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) multivector field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: A = M.multivector_module(2)
            sage: A._an_element_()
            2-vector field on the 2-dimensional differentiable
             manifold M
        """
        resu = self.element_class(self._vmodule, self._degree)
        for oc in self._domain.open_covers(trivial=False):
            # the first non-trivial open cover is selected
            for dom in oc:
                vmodule_dom = dom.vector_field_module(
                                  dest_map=self._dest_map.restrict(dom))
                dmodule_dom = vmodule_dom.exterior_power(self._degree)
                resu.set_restriction(dmodule_dom._an_element_())
            return resu
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: A2 = M.multivector_module(2)
            sage: A2._coerce_map_from_(M.tensor_field_module((2,0)))
            False
            sage: U = M.open_subset('U')
            sage: A2U = U.multivector_module(2)
            sage: A2U._coerce_map_from_(A2)
            True
            sage: A2._coerce_map_from_(A2U)
            False
        """
        if isinstance(other, (MultivectorModule, MultivectorFreeModule)):
            # coercion by domain restriction
            return (self._degree == other._degree
                    and self._domain.is_subset(other._domain)
                    and self._ambient_domain.is_subset(
                                                 other._ambient_domain))
        return False

    @cached_method
    def zero(self):
        """
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: A2 = M.multivector_module(2)
            sage: A2.zero()
            2-vector field zero on the 3-dimensional differentiable
             manifold M
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
            sage: A2 = M.multivector_module(2)
            sage: A2
            Module A^2(M) of 2-vector fields on the 3-dimensional
             differentiable manifold M
        """
        description = "Module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-vector fields ".format(self._degree)
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
            sage: A2 = M.multivector_module(2)
            sage: A2._latex_()
            'A^{2}\\left(\\mathcal{M}\\right)'
            sage: latex(A2)  # indirect doctest
            A^{2}\left(\mathcal{M}\right)
        """
        if self._latex_name is None:
            return r'\text{' + str(self) + r'}'
        else:
            return self._latex_name

    def base_module(self):
        r"""
        Return the vector field module on which the multivector field
        module ``self`` is constructed.

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`
          representing the module on which ``self`` is defined

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: A2 = M.multivector_module(2) ; A2
            Module A^2(M) of 2-vector fields on the 3-dimensional
             differentiable manifold M
            sage: A2.base_module()
            Module X(M) of vector fields on the 3-dimensional
             differentiable manifold M
            sage: A2.base_module() is M.vector_field_module()
            True
            sage: U = M.open_subset('U')
            sage: A2U = U.multivector_module(2) ; A2U
            Module A^2(U) of 2-vector fields on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: A2U.base_module()
            Module X(U) of vector fields on the Open subset U of the
             3-dimensional differentiable manifold M
        """
        return self._vmodule

    def degree(self):
        r"""
        Return the degree of the multivector fields in ``self``.

        OUTPUT: integer `p` such that ``self`` is a set of `p`-vector fields

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: M.multivector_module(2).degree()
            2
            sage: M.multivector_module(3).degree()
            3
        """
        return self._degree

#***********************************************************************


class MultivectorFreeModule(ExtPowerFreeModule):
    r"""
    Free module of multivector fields of a given degree `p` (`p`-vector
    fields) along a differentiable manifold `U` with values on a
    parallelizable manifold `M`.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi:\; U \rightarrow M` to a parallelizable manifold `M` of dimension
    `n`, the set `A^p(U, \Phi)` of `p`-vector fields (i.e. alternating tensor
    fields of type `(p,0)`) along `U` with values on `M` is a free module
    of rank `\binom{n}{p}` over `C^k(U)`, the commutative algebra of
    differentiable scalar fields on `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).
    The standard case of `p`-vector fields *on* a differentiable
    manifold `M` corresponds to `U = M` and `\Phi = \mathrm{Id}_M`.
    Other common cases are `\Phi` being an immersion and `\Phi` being a
    curve in `M` (`U` is then an open interval of `\RR`).

    .. NOTE::

        This class implements `A^p(U, \Phi)` in the case where `M` is
        parallelizable; `A^p(U, \Phi)` is then a *free* module. If `M`
        is not parallelizable, the class :class:`MultivectorModule` must
        be used instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U,\Phi)` of
      vector fields along `U` associated with the map
      `\Phi: U \rightarrow V`
    - ``degree`` -- positive integer; the degree `p` of the multivector
      fields

    EXAMPLES:

    Free module of 2-vector fields on a parallelizable 3-dimensional
    manifold::

        sage: M = Manifold(3, 'M')
        sage: X.<x,y,z> = M.chart()
        sage: XM = M.vector_field_module() ; XM
        Free module X(M) of vector fields on the 3-dimensional
         differentiable manifold M
        sage: A = M.multivector_module(2) ; A
        Free module A^2(M) of 2-vector fields on the 3-dimensional
         differentiable manifold M
        sage: latex(A)
        A^{2}\left(M\right)

    ``A`` is nothing but the second exterior power of ``XM``, i.e. we
    have `A^{2}(M) = \Lambda^2(\mathfrak{X}(M))` (see
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule`)::

        sage: A is XM.exterior_power(2)
        True

    `A^{2}(M)` is a module over the algebra `C^k(M)` of (differentiable)
    scalar fields on `M`::

        sage: A.category()
        Category of finite dimensional modules over Algebra of
         differentiable scalar fields on the 3-dimensional
         differentiable manifold M
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
        2-vector field zero on the 3-dimensional differentiable
         manifold M
        sage: A(0) is A.zero()
        True

    while nonzero elements are constructed by providing their
    components in a given vector frame::

        sage: comp = [[0,3*x,-z],[-3*x,0,4],[z,-4,0]]
        sage: a = A(comp, frame=X.frame(), name='a') ; a
        2-vector field a on the 3-dimensional differentiable manifold M
        sage: a.display()
        a = 3*x ∂/∂x∧∂/∂y - z ∂/∂x∧∂/∂z + 4 ∂/∂y∧∂/∂z

    An alternative is to construct the 2-vector field from an empty list
    of components and to set the nonzero nonredundant components
    afterwards::

        sage: a = A([], name='a')
        sage: a[0,1] = 3*x  # component in the manifold's default frame
        sage: a[0,2] = -z
        sage: a[1,2] = 4
        sage: a.display()
        a = 3*x ∂/∂x∧∂/∂y - z ∂/∂x∧∂/∂z + 4 ∂/∂y∧∂/∂z

    The module `A^1(M)` is nothing but `\mathfrak{X}(M)` (the free module
    of vector fields on `M`)::

        sage: A1 = M.multivector_module(1) ; A1
        Free module X(M) of vector fields on the 3-dimensional
         differentiable manifold M
        sage: A1 is XM
        True

    There is a coercion map `A^p(M) \rightarrow T^{(p,0)}(M)`::

        sage: T20 = M.tensor_field_module((2,0)); T20
        Free module T^(2,0)(M) of type-(2,0) tensors fields on the
         3-dimensional differentiable manifold M
        sage: T20.has_coerce_map_from(A)
        True

    but of course not in the reverse direction, since not all contravariant
    tensor field is alternating::

        sage: A.has_coerce_map_from(T20)
        False

    The coercion map `A^2(M) \rightarrow T^{(2,0)}(M)` in action::

        sage: T20 = M.tensor_field_module((2,0)) ; T20
        Free module T^(2,0)(M) of type-(2,0) tensors fields on the
         3-dimensional differentiable manifold M
        sage: ta = T20(a) ; ta
        Tensor field a of type (2,0) on the 3-dimensional differentiable
         manifold M
        sage: ta.display()
        a = 3*x ∂/∂x⊗∂/∂y - z ∂/∂x⊗∂/∂z - 3*x ∂/∂y⊗∂/∂x + 4 ∂/∂y⊗∂/∂z
         + z ∂/∂z⊗∂/∂x - 4 ∂/∂z⊗∂/∂y
        sage: a.display()
        a = 3*x ∂/∂x∧∂/∂y - z ∂/∂x∧∂/∂z + 4 ∂/∂y∧∂/∂z
        sage: ta.symmetries()  # the antisymmetry is preserved
        no symmetry;  antisymmetry: (0, 1)

    There is also coercion to subdomains, which is nothing but the
    restriction of the multivector field to some subset of its domain::

        sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1})
        sage: B = U.multivector_module(2) ; B
        Free module A^2(U) of 2-vector fields on the Open subset U of the
         3-dimensional differentiable manifold M
        sage: B.has_coerce_map_from(A)
        True
        sage: a_U = B(a) ; a_U
        2-vector field a on the Open subset U of the 3-dimensional
         differentiable manifold M
        sage: a_U.display()
        a = 3*x ∂/∂x∧∂/∂y - z ∂/∂x∧∂/∂z + 4 ∂/∂y∧∂/∂z
    """

    Element = MultivectorFieldParal

    def __init__(self, vector_field_module, degree):
        r"""
        Construct a free module of multivector fields.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: from sage.manifolds.differentiable.multivector_module \
            ....:      import MultivectorFreeModule
            sage: A = MultivectorFreeModule(M.vector_field_module(), 2)
            sage: A
            Free module A^2(M) of 2-vector fields on the 3-dimensional
             differentiable manifold M
            sage: TestSuite(A).run()
        """
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "A^{}(".format(degree) + domain._name
        latex_name = r"A^{{{}}}\left({}".format(degree,
                                                domain._latex_name)
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
        ExtPowerFreeModule.__init__(self, vector_field_module, degree,
                                    name=name, latex_name=latex_name)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a multivector field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: A = M.multivector_module(2)
            sage: a = A([[0, x], [-x, 0]], name='a'); a
            2-vector field a on the 2-dimensional differentiable
             manifold M
            sage: a.display()
            a = x ∂/∂x∧∂/∂y
            sage: A(0) is A.zero()
            True
        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, (MultivectorField, MultivectorFieldParal)):
            # coercion by domain restriction
            if (self._degree == comp._tensor_type[0]
                    and self._domain.is_subset(comp._domain)
                    and self._ambient_domain.is_subset(
                                                 comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise TypeError("cannot convert the {} ".format(comp) +
                                "to a multivector field in {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._fmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    # Rem: _an_element_ is declared in the superclass ExtPowerFreeModule

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: A2 = M.multivector_module(2)
            sage: U = M.open_subset('U', coord_def = {X: z<0})
            sage: A2U = U.multivector_module(2)
            sage: A2U._coerce_map_from_(A2)
            True
            sage: A2._coerce_map_from_(A2U)
            False
            sage: A1 = M.multivector_module(1)
            sage: A2U._coerce_map_from_(A1)
            False
            sage: A2._coerce_map_from_(M.tensor_field_module((2,0)))
            False
        """
        if isinstance(other, (MultivectorModule, MultivectorFreeModule)):
            # coercion by domain restriction
            return (self._degree == other._degree
                    and self._domain.is_subset(other._domain)
                    and self._ambient_domain.is_subset(
                                                 other._ambient_domain))
        return False

    #### End of Parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: A = M.multivector_module(2)
            sage: A
            Free module A^2(M) of 2-vector fields on
             the 3-dimensional differentiable manifold M
        """
        description = "Free module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-vector fields ".format(self._degree)
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {}".format(
                                     self._domain, self._ambient_domain)
        return description
