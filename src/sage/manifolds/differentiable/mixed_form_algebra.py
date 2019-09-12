r"""
Graded Algebra of Mixed Differential Forms

Let `M` and `N` be differentiable manifolds and `\varphi: M \to N` a
differentiable map. The space of *mixed differential forms along* `\varphi`,
denoted by `\Omega^*(M,\varphi)`, is given by the direct sum
`\bigoplus^n_{j=0} \Omega^j(M,\varphi)` of differential form modules, where
`n=\dim(N)`. With the wedge product, `\Omega^*(M,\varphi)` inherits the
structure of a graded algebra.

AUTHORS:

- Michael Jung (2019) : initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.categories.graded_algebras import GradedAlgebras
from sage.structure.unique_representation import UniqueRepresentation
from sage.symbolic.ring import ZZ, SR
from sage.manifolds.differentiable.mixed_form import MixedForm

class MixedFormAlgebra(Parent, UniqueRepresentation):
    r"""
    An instance of this class represents the graded algebra of mixed form. That
    is, if `\varphi: M \to N` is a differentiable map between two differentiable
    manifolds `M` and `N`, the *graded algebra of mixed forms*
    `\Omega^*(M,\varphi)` *along* `\varphi` is defined via the direct sum
    `\bigoplus^{n}_{j=0} \Omega^j(M,\varphi)` consisting of differential form
    modules
    (cf. :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormModule`),
    where `n` is the dimension of `N`. Hence, `\Omega^*(M,\varphi)` is a module
    over `C^k(M)` and a vector space over `\RR` or `\CC`. Furthermore notice,
    that

    .. MATH::

        \Omega^*(M,\varphi) \cong C^k \left( \bigoplus^n_{j=0} \Lambda^j(\varphi^*T^*N) \right),

    where `C^k` denotes the global section functor for differentiable sections
    of order `k` here.

    The wedge product induces a multiplication on `\Omega^*(M,\varphi)` and
    gives it the structure of a graded algebra since

    .. MATH::

        \Omega^k(M,\varphi) \wedge \Omega^l(M,\varphi) \subset \Omega^{k+l}(M,\varphi).

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(M,\varphi)` of vector
      fields along `M` associated with the map `\varphi: M \rightarrow N`

    EXAMPLES:

    Graded algebra of mixed forms on a 3-dimensional manifold::

        sage: M = Manifold(3, 'M')
        sage: X.<x,y,z> = M.chart()
        sage: Omega = M.mixed_form_algebra(); Omega
        Graded algebra Omega^*(M) of mixed differential forms on the
         3-dimensional differentiable manifold M
        sage: Omega.category()
        Category of graded algebras over Symbolic Ring
        sage: Omega.base_ring()
        Symbolic Ring
        sage: Omega.vector_field_module()
        Free module X(M) of vector fields on the 3-dimensional differentiable
         manifold M

    Elements can be created from scratch::

        sage: A = Omega(0); A
        Mixed differential form zero on the 3-dimensional differentiable
         manifold M
        sage: A is Omega.zero()
        True
        sage: B = Omega(1); B
        Mixed differential form one on the 3-dimensional differentiable
         manifold M
        sage: B is Omega.one()
        True
        sage: C = Omega([2,0,0,0]); C
        Mixed differential form on the 3-dimensional differentiable manifold M

    There are some important coercions implemented::

        sage: Omega0 = M.scalar_field_algebra(); Omega0
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M
        sage: Omega.has_coerce_map_from(Omega0)
        True
        sage: Omega2 = M.diff_form_module(2); Omega2
        Free module Omega^2(M) of 2-forms on the 3-dimensional differentiable
         manifold M
        sage: Omega.has_coerce_map_from(Omega2)
        True

    Restrictions induce coercions as well::

        sage: U = M.open_subset('U'); U
        Open subset U of the 3-dimensional differentiable manifold M
        sage: OmegaU = U.mixed_form_algebra(); OmegaU
        Graded algebra Omega^*(U) of mixed differential forms on the Open subset
         U of the 3-dimensional differentiable manifold M
        sage: OmegaU.has_coerce_map_from(Omega)
        True

    """
    Element = MixedForm

    def __init__(self, vector_field_module):
        r"""
        Construct a graded algebra of mixed forms.

        TESTS:

        Graded algebra of mixed forms on a non-parallelizable 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                   intersection_name='W', restrictions1= x>0,
            ....:                   restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: from sage.manifolds.differentiable.mixed_form_algebra import (
            ....:                                              MixedFormAlgebra)
            sage: A = MixedFormAlgebra(M.vector_field_module())
            sage: TestSuite(A).run()

        """
        if vector_field_module is None:
            raise ValueError("underlying vector field module must be provided")
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        # Set name and latex_name:
        name = "Omega^*(" + domain._name
        latex_name = r"\Omega^*\left(" + domain._latex_name
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
        # Add this algebra to the category of graded algebras:
        base_field = domain.base_field()
        if domain.base_field_type() in ['real', 'complex']:
            base_field = SR
        Parent.__init__(self, base=base_field,
                        category=GradedAlgebras(base_field))
        # Define attributes:
        self._domain = domain
        self._ambient_domain = vector_field_module._ambient_domain
        self._dest_map = dest_map
        self._vmodule = vector_field_module
        self._max_deg = vector_field_module._ambient_domain.dim()

    def _element_constructor_(self, comp=None, name=None, latex_name=None):
        r"""
        Construct a mixed form.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: A = M.mixed_form_algebra()
            sage: a = A([x,0,0], name='a'); a
            Mixed differential form a on the 2-dimensional differentiable
             manifold M

        """
        if comp is None:
            return self.element_class(self, name=name, latex_name=latex_name)
        # Treat one and zero separately for cache:
        elif comp in ZZ and comp == 0:
            return self.zero()
        elif comp in ZZ and comp == 1:
            return self.one()
        # Prepare list:
        if isinstance(comp, list):
            if len(comp) != self._max_deg + 1:
                raise IndexError(
                    "input list must have length {}".format(self._max_deg + 1))
            comp_list = comp
        elif isinstance(comp, self.Element):
            comp_list = comp[:]
        elif comp in self._domain.scalar_field_algebra():
            comp_list = [0] * (self._max_deg + 1)
            comp_list[0] = comp
        else:
            comp_list = [0] * (self._max_deg + 1)
            # Try...except?
            deg = comp.degree()
            if deg <= self._max_deg:
                comp_list[deg] = comp
        # Use already existing coercions:
        comp_list = [self._domain.diff_form_module(j,
                            self._dest_map)(comp_list[j])
                     for j in range(self._max_deg + 1)]
        # Now, define names:
        try:
            if name is None and comp._name is not None:
                name = comp._name
            if latex_name is None and comp._latex_name is not None:
                latex_name = comp._latex_name
        except AttributeError:
            # AttributeError? Then comp might be in SR:
            if comp in SR:
                if name is None:
                    name = repr(comp)
                if latex_name is None:
                    from sage.misc.latex import latex
                    latex_name = latex(comp)
        return self.element_class(self, comp_list, name, latex_name)

    def _an_element_(self):
        r"""
        Construct some (unnamed) mixed form.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: A = M.mixed_form_algebra()
            sage: A._an_element_()
            Mixed differential form on the 2-dimensional differentiable
             manifold M

        """
        resu_comp = [self._domain.diff_form_module(j,
                            self._dest_map)._an_element_()
                     for j in range(self._max_deg + 1)]
        return self.element_class(self, comp=resu_comp)

    def _coerce_map_from_(self, S):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: A = M.mixed_form_algebra()
            sage: A._coerce_map_from_(M.diff_form_module(0))
            True
            sage: A._coerce_map_from_(M.diff_form_module(1))
            True
            sage: A._coerce_map_from_(M.diff_form_module(2))
            True
            sage: A._coerce_map_from_(M.diff_form_module(3))
            True
            sage: U = M.open_subset('U')
            sage: AU = U.mixed_form_algebra()
            sage: AU._coerce_map_from_(A)
            True
            sage: A._coerce_map_from_(AU)
            False

        """
        if isinstance(S, self.__class__):
            # coercion by domain restriction
            return (self._domain.is_subset(S._domain) and
                   self._ambient_domain.is_subset(S._ambient_domain))
        # Test scalar_field_algebra separately to ensure coercion from SR:
        if self._domain.scalar_field_algebra().has_coerce_map_from(S):
            return True
        # This is tricky, we need to check the degree first:
        try:
            deg = S.degree()
            if self._domain.diff_form_module(deg,
                                        self._dest_map).has_coerce_map_from(S):
                return True
        except (NotImplementedError, AttributeError, TypeError):
            pass
        return False

    @cached_method
    def zero(self):
        r"""
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: A = M.mixed_form_algebra()
            sage: A.zero()
            Mixed differential form zero on the 3-dimensional differentiable
             manifold M

        """
        dmap = self._dest_map
        # Prepare and fill comp list:
        resu_comp = list()
        resu_comp.append(self._domain.scalar_field_algebra().zero())
        resu_comp.extend([self._domain.diff_form_module(j, dest_map=dmap).zero()
                          for j in range(1, self._max_deg + 1)])
        return self.element_class(self, comp=resu_comp, name='zero',
                                  latex_name='0')

    @cached_method
    def one(self):
        r"""
        Return the one of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: A = M.mixed_form_algebra()
            sage: A.one()
            Mixed differential form one on the 3-dimensional differentiable
             manifold M

        """
        dmap = self._dest_map
        # Prepare and fill comp list:
        resu_comp = list()
        resu_comp.append(self._domain.scalar_field_algebra().one())
        resu_comp.extend([self._domain.diff_form_module(j, dest_map=dmap).zero()
                          for j in range(1, self._max_deg + 1)])
        return self.element_class(self, comp=resu_comp, name='one',
                                          latex_name='1')

    def vector_field_module(self):
        r"""
        Return the underlying vector field module.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: N = Manifold(3, 'N')
            sage: Phi = M.diff_map(N, name='Phi'); Phi
            Differentiable map Phi from the 2-dimensional differentiable manifold M
             to the 3-dimensional differentiable manifold N
            sage: A = M.mixed_form_algebra(Phi); A
            Graded algebra Omega^*(M,Phi) of mixed differential forms along the
             2-dimensional differentiable manifold M mapped into the 3-dimensional
             differentiable manifold N via Phi
            sage: A.vector_field_module()
            Module X(M,Phi) of vector fields along the 2-dimensional differentiable
             manifold M mapped into the 3-dimensional differentiable manifold N

        """
        return self._vmodule

    def _repr_(self):
        r"""
        Return a string representation of the object.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: A = M.mixed_form_algebra(); A
            Graded algebra Omega^*(M) of mixed differential forms on the
             3-dimensional differentiable manifold M

        """
        description = ("Graded algebra " + self._name +
                      " of mixed differential forms ")
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {} ".format(
                self._domain, self._ambient_domain)
            if self._dest_map._name is None:
                dm_name = "unnamed map"
            else:
                dm_name = self._dest_map._name
            description += "via " + dm_name
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(3, 'M', latex_name=r'\mathcal{M}')
            sage: A = M.mixed_form_algebra()
            sage: A._latex_()
            '\\Omega^*\\left(\\mathcal{M}\\right)'
            sage: latex(A)  # indirect doctest
            \Omega^*\left(\mathcal{M}\right)

        """
        return self._latex_name

