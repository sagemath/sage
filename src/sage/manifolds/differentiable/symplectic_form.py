# sage_setup: distribution = sagemath-symbolics
r"""
Symplectic structures

The class :class:`SymplecticForm` implements symplectic structures
on differentiable manifolds over `\RR`. The derived class
:class:`SymplecticFormParal` is devoted to symplectic forms on a
parallelizable manifold.

AUTHORS:

- Tobias Diez (2021) : initial version

REFERENCES:

- [AM1990]_
- [RS2012]_
"""
# *****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import annotations
from typing import Union, Optional

from sage.symbolic.expression import Expression
from sage.manifolds.differentiable.diff_form import DiffForm, DiffFormParal
from sage.manifolds.differentiable.diff_map import DiffMap
from sage.manifolds.differentiable.vectorfield_module import VectorFieldModule
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.differentiable.scalarfield import DiffScalarField
from sage.manifolds.differentiable.vectorfield import VectorField
from sage.manifolds.differentiable.poisson_tensor import PoissonTensorField


class SymplecticForm(DiffForm):
    r"""
    A symplectic form on a differentiable manifold.

    An instance of this class is a closed nondegenerate differential `2`-form `\omega`
    on a differentiable manifold `M` over `\RR`.

    In particular, at each point `m \in M`, `\omega_m` is a bilinear map of the type:

    .. MATH::

        \omega_m:\ T_m M \times T_m M  \to \RR,

    where `T_m M` stands for the tangent space to the
    manifold `M` at the point `m`, such that `\omega_m` is skew-symmetric:
    `\forall u,v \in T_m M, \ \omega_m(v,u) = - \omega_m(u,v)`
    and nondegenerate:
    `(\forall v \in T_m M,\ \ \omega_m(u,v) = 0) \Longrightarrow u=0`.

    .. NOTE::

        If `M` is parallelizable, the class :class:`SymplecticFormParal`
        should be used instead.

    INPUT:

    - ``manifold`` -- module `\mathfrak{X}(M)` of vector fields on the
      manifold `M`, or the manifold `M` itself
    - ``name`` -- (default: ``omega``) name given to the symplectic form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      symplectic form; if ``None``, it is formed from ``name``

    EXAMPLES:

    A symplectic form on the 2-sphere::

        sage: M.<x,y> = manifolds.Sphere(2, coordinates='stereographic')
        sage: stereoN = M.stereographic_coordinates(pole='north')
        sage: stereoS = M.stereographic_coordinates(pole='south')
        sage: omega = M.symplectic_form(name='omega', latex_name=r'\omega')
        sage: omega
        Symplectic form omega on the 2-sphere S^2 of radius 1 smoothly embedded
         in the Euclidean space E^3

    ``omega`` is initialized by providing its single nonvanishing component
    w.r.t. the vector frame associated to ``stereoN``, which is the default
    frame on ``M``::

        sage: omega[1, 2] = 1/(1 + x^2 + y^2)^2

    The components w.r.t. the vector frame associated to ``stereoS`` are
    obtained thanks to the method
    :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.add_comp_by_continuation`::

        sage: omega.add_comp_by_continuation(stereoS.frame(),
        ....:                  stereoS.domain().intersection(stereoN.domain()))
        sage: omega.display()
        omega = (x^2 + y^2 + 1)^(-2) dx∧dy
        sage: omega.display(stereoS)
        omega = -1/(xp^4 + yp^4 + 2*(xp^2 + 1)*yp^2 + 2*xp^2 + 1) dxp∧dyp

    ``omega`` is an exact 2-form (this is trivial here, since ``M`` is
    2-dimensional)::

        sage: diff(omega).display()
        domega = 0
    """

    _name: str
    _latex_name: str
    _dim_half: int
    _poisson: Optional[PoissonTensorField]
    _vol_form: Optional[DiffForm]
    _restrictions: dict[DifferentiableManifold, SymplecticForm]

    def __init__(
        self,
        manifold: Union[DifferentiableManifold, VectorFieldModule],
        name: Optional[str] = None,
        latex_name: Optional[str] = None,
    ):
        r"""
        Construct a symplectic form.

        TESTS::

            sage: from sage.manifolds.differentiable.symplectic_form import SymplecticForm
            sage: M = manifolds.Sphere(2, coordinates='stereographic')
            sage: omega = SymplecticForm(M, name='omega', latex_name=r'\omega')
            sage: omega
            Symplectic form omega on the 2-sphere S^2 of radius 1 smoothly
             embedded in the Euclidean space E^3
        """
        try:
            vector_field_module = manifold.vector_field_module()
        except AttributeError:
            vector_field_module = manifold

        if name is None:
            name = "omega"
            if latex_name is None:
                latex_name = r"\omega"

        if latex_name is None:
            latex_name = name

        DiffForm.__init__(
            self, vector_field_module, 2, name=name, latex_name=latex_name
        )

        # Check that manifold is even dimensional
        dim = self._ambient_domain.dimension()
        if dim % 2 == 1:
            raise ValueError(
                f"the dimension of the manifold must be even but it is {dim}"
            )
        self._dim_half = dim // 2

        # Initialization of derived quantities
        SymplecticForm._init_derived(self)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M.<q, p> = EuclideanSpace(2)
            sage: omega = M.symplectic_form('omega', r'\omega'); omega
            Symplectic form omega on the Euclidean plane E^2
        """
        return self._final_repr(f"Symplectic form {self._name} ")

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``.

        TESTS::

            sage: M.<q, p> = EuclideanSpace(2)
            sage: omega = M.symplectic_form('omega', r'\omega')._new_instance(); omega
            Symplectic form unnamed symplectic form on the Euclidean plane E^2
        """
        return type(self)(
            self._vmodule,
            "unnamed symplectic form",
            latex_name=r"\text{unnamed symplectic form}",
        )

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M.<q, p> = EuclideanSpace(2)
            sage: omega = M.symplectic_form('omega', r'\omega')._init_derived()
        """
        # Initialization of quantities pertaining to the mother class
        DiffForm._init_derived(self)

        self._poisson = None
        self._vol_form = None

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M.<q, p> = EuclideanSpace(2)
            sage: omega = M.symplectic_form('omega', r'\omega')._del_derived()
        """
        # Delete the derived quantities from the mother class
        DiffForm._del_derived(self)

        # Clear the Poisson tensor
        if self._poisson is not None:
            self._poisson._restrictions.clear()
            self._poisson._del_derived()
            self._poisson = None

        # Delete the volume form
        if self._vol_form is not None:
            self._vol_form._restrictions.clear()
            self._vol_form._del_derived()
            self._vol_form = None

    def restrict(
        self, subdomain: DifferentiableManifold, dest_map: Optional[DiffMap] = None
    ) -> DiffForm:
        r"""
        Return the restriction of the symplectic form to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of the symplectic form's domain
        - ``dest_map`` -- (default: ``None``) smooth destination map
          `\Phi:\ U \to V`, where `V` is a subdomain of the symplectic form's domain
          If ``None``, the restriction of the initial vector field module is used.

        OUTPUT: the restricted symplectic form

        EXAMPLES::

            sage: M = Manifold(6, 'M')
            sage: omega = M.symplectic_form()
            sage: U = M.open_subset('U')
            sage: omega.restrict(U)
            2-form omega on the Open subset U of the 6-dimensional differentiable manifold M
        """
        if subdomain == self._domain:
            return self

        if subdomain not in self._restrictions:
            # Construct the restriction at the tensor field level
            restriction = DiffForm.restrict(self, subdomain, dest_map=dest_map)

            # Restrictions of derived quantities
            if self._poisson is not None:
                restriction._poisson = self._poisson.restrict(subdomain)
            if self._vol_form is not None:
                restriction._vol_form = self._vol_form.restrict(subdomain)

            # The restriction is ready
            self._restrictions[subdomain] = restriction
            return restriction
        else:
            return self._restrictions[subdomain]

    @staticmethod
    def wrap(
        form: DiffForm, name: Optional[str] = None, latex_name: Optional[str] = None
    ) -> SymplecticForm:
        r"""
        Define the symplectic form from a differential form.

        INPUT:

        - ``form`` -- differential `2`-form

        EXAMPLES:

        Volume form on the sphere as a symplectic form::

            sage: from sage.manifolds.differentiable.symplectic_form import SymplecticForm
            sage: M = manifolds.Sphere(2, coordinates='stereographic')
            sage: vol_form = M.induced_metric().volume_form()                   # long time
            sage: omega = SymplecticForm.wrap(vol_form, 'omega', r'\omega')     # long time
            sage: omega.display()                                               # long time
            omega = -4/(y1^4 + y2^4 + 2*(y1^2 + 1)*y2^2 + 2*y1^2 + 1) dy1∧dy2
        """
        if form.degree() != 2:
            raise TypeError("the argument must be a form of degree 2")

        if name is None:
            name = form._name
        if latex_name is None:
            latex_name = form._latex_name

        symplectic_form = form.base_module().symplectic_form(name, latex_name)

        for dom, rst in form._restrictions.items():
            symplectic_form._restrictions[dom] = SymplecticForm.wrap(
                rst, name, latex_name
            )

        if isinstance(form, DiffFormParal):
            for frame in form._components:
                symplectic_form._components[frame] = form._components[frame].copy()

        return symplectic_form

    def poisson(
        self, expansion_symbol: Optional[Expression] = None, order: int = 1
    ) -> PoissonTensorField:
        r"""
        Return the Poisson tensor associated with the symplectic form.

        INPUT:

        - ``expansion_symbol`` -- (default: ``None``) symbolic variable; if
          specified, the inverse will be expanded in power series with respect
          to this variable (around its zero value)
        - ``order`` -- integer (default: 1); the order of the expansion
          if ``expansion_symbol`` is not ``None``; the *order* is defined as
          the degree of the polynomial representing the truncated power series
          in ``expansion_symbol``; currently only first order inverse is
          supported

        If ``expansion_symbol`` is set, then the zeroth order symplectic form must be
        invertible. Moreover, subsequent calls to this method will return
        a cached value, even when called with the default value (to enable
        computation of derived quantities). To reset, use :meth:`_del_derived`.

        OUTPUT:

        - the Poisson tensor, as an instance of
          :meth:`~sage.manifolds.differentiable.poisson_tensor.PoissonTensorField`

        EXAMPLES:

        Poisson tensor of `2`-dimensional symplectic vector space::

            sage: M = manifolds.StandardSymplecticSpace(2)
            sage: omega = M.symplectic_form()
            sage: poisson = omega.poisson(); poisson
            2-vector field poisson_omega on the Standard symplectic space R2
            sage: poisson.display()
            poisson_omega = -e_q∧e_p
        """
        if self._poisson is None:
            # Initialize the Poisson tensor
            poisson_name = f"poisson_{self._name}"
            poisson_latex_name = f"{self._latex_name}^{{-1}}"
            self._poisson = self._vmodule.poisson_tensor(
                poisson_name, poisson_latex_name
            )

        # Update the Poisson tensor
        # TODO: Should this be done instead when a new restriction is added?
        for domain, restriction in self._restrictions.items():
            # Forces the update of the restriction
            self._poisson._restrictions[domain] = restriction.poisson(
                expansion_symbol=expansion_symbol, order=order
            )
        return self._poisson

    def hamiltonian_vector_field(self, function: DiffScalarField) -> VectorField:
        r"""
        The Hamiltonian vector field `X_f` generated by a function `f: M \to \RR`.

        The Hamiltonian vector field is defined by

        .. MATH::

            X_f \lrcorner \omega + df = 0.

        INPUT:

        - ``function`` -- the function generating the Hamiltonian vector field

        EXAMPLES::

            sage: M = manifolds.StandardSymplecticSpace(2)
            sage: omega = M.symplectic_form()
            sage: f = M.scalar_field({ chart: function('f')(*chart[:]) for chart in M.atlas() }, name='f')
            sage: f.display()
            f: R2 → ℝ
               (q, p) ↦ f(q, p)
            sage: Xf = omega.hamiltonian_vector_field(f)
            sage: Xf.display()
            Xf = d(f)/dp e_q - d(f)/dq e_p
        """
        return self.poisson().hamiltonian_vector_field(function)

    def flat(self, vector_field: VectorField) -> DiffForm:
        r"""
        Return the image of the given differential form under the
        map `\omega^\flat: T M \to T^*M` defined by

        .. MATH::

            <\omega^\flat(X), Y> = \omega_m (X, Y)

        for all `X, Y \in T_m M`.

        In indices, `X_i = \omega_{ji} X^j`.

        INPUT:

        - ``vector_field`` -- the vector field to calculate its flat of

        EXAMPLES::

            sage: M = manifolds.StandardSymplecticSpace(2)
            sage: omega = M.symplectic_form()
            sage: X = M.vector_field_module().an_element()
            sage: X.set_name('X')
            sage: X.display()
            X = 2 e_q + 2 e_p
            sage: omega.flat(X).display()
            X_flat = 2 dq - 2 dp
        """
        form = vector_field.down(self)
        form.set_name(
            vector_field._name + "_flat", vector_field._latex_name + "^\\flat"
        )
        return form

    def sharp(self, form: DiffForm) -> VectorField:
        r"""
        Return the image of the given differential form under the map
        `\omega^\sharp: T^* M \to TM` defined by

        .. MATH::

            \omega (\omega^\sharp(\alpha), X) = \alpha(X)

        for all `X \in T_m M` and `\alpha \in T^*_m M`.
        The sharp map is inverse to the flat map.

        In indices, `\alpha^i = \varpi^{ij} \alpha_j`, where `\varpi` is
        the Poisson tensor associated with the symplectic form.

        INPUT:

        - ``form`` -- the differential form to calculate its sharp of

        EXAMPLES::

            sage: M = manifolds.StandardSymplecticSpace(2)
            sage: omega = M.symplectic_form()
            sage: X = M.vector_field_module().an_element()
            sage: alpha = omega.flat(X)
            sage: alpha.set_name('alpha')
            sage: alpha.display()
            alpha = 2 dq - 2 dp
            sage: omega.sharp(alpha).display()
            alpha_sharp = 2 e_q + 2 e_p
        """
        return self.poisson().sharp(form)

    def poisson_bracket(
        self, f: DiffScalarField, g: DiffScalarField
    ) -> DiffScalarField:
        r"""
        Return the Poisson bracket

        .. MATH::

            \{f, g\} = \omega(X_f, X_g)

        of the given functions.

        INPUT:

        - ``f`` -- function inserted in the first slot
        - ``g`` -- function inserted in the second slot

        EXAMPLES::

            sage: M.<q, p> = EuclideanSpace(2)
            sage: poisson = M.poisson_tensor('varpi')
            sage: poisson.set_comp()[1,2] = -1
            sage: f = M.scalar_field({ chart: function('f')(*chart[:]) for chart in M.atlas() }, name='f')
            sage: g = M.scalar_field({ chart: function('g')(*chart[:]) for chart in M.atlas() }, name='g')
            sage: poisson.poisson_bracket(f, g).display()
            poisson(f, g): E^2 → ℝ
               (q, p) ↦ d(f)/dp*d(g)/dq - d(f)/dq*d(g)/dp
        """
        return self.poisson().poisson_bracket(f, g)

    def volume_form(self, contra: int = 0) -> TensorField:
        r"""
        Liouville volume form `\frac{1}{n!}\omega^n` associated with the
        symplectic form `\omega`, where `2n` is the dimension of the manifold.

        INPUT:

        - ``contra`` -- (default: 0) number of contravariant indices of the
          returned tensor

        OUTPUT:

        - if ``contra = 0``: volume form associated with the symplectic form
        - if ``contra = k``, with `1\leq k \leq n`, the tensor field of type
          (k,n-k) formed from `\epsilon` by raising the first k indices with
          the symplectic form (see method
          :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.up`)

        EXAMPLES:

        Volume form on `\RR^4`::

            sage: M = manifolds.StandardSymplecticSpace(4)
            sage: omega = M.symplectic_form()
            sage: vol = omega.volume_form() ; vol
            4-form mu_omega on the Standard symplectic space R4
            sage: vol.display()
            mu_omega = dq1∧dp1∧dq2∧dp2
        """
        if self._vol_form is None:
            from sage.functions.other import factorial

            vol_form = self
            for _ in range(1, self._dim_half):
                vol_form = vol_form.wedge(self)
            vol_form = vol_form / factorial(self._dim_half)

            volume_name = f"mu_{self._name}"
            volume_latex_name = f"\\mu_{self._latex_name}"
            vol_form.set_name(volume_name, volume_latex_name)
            self._vol_form = vol_form

        result = self._vol_form
        for k in range(0, contra):
            result = result.up(self, k)
        if contra > 1:
            # restoring the antisymmetry after the up operation:
            result = result.antisymmetrize(*range(contra))
        return result

    def hodge_star(self, pform: DiffForm) -> DiffForm:
        r"""
        Compute the Hodge dual of a differential form with respect to the
        symplectic form.

        See :meth:`~sage.manifolds.differentiable.diff_form.DiffForm.hodge_dual`
        for the definition and more details.

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

        Hodge dual of any form on the symplectic vector space `R^2`::

            sage: M = manifolds.StandardSymplecticSpace(2)
            sage: omega = M.symplectic_form()
            sage: a = M.one_form(1, 0, name='a')
            sage: omega.hodge_star(a).display()
            *a = dq
            sage: b = M.one_form(0, 1, name='b')
            sage: omega.hodge_star(b).display()
            *b = dp
            sage: f = M.scalar_field(1, name='f')
            sage: omega.hodge_star(f).display()
            *f = -dq∧dp
            sage: omega.hodge_star(omega).display()
            *omega: R2 → ℝ
               (q, p) ↦ 1
        """
        return pform.hodge_dual(self)

    def on_forms(self, first: DiffForm, second: DiffForm) -> DiffScalarField:
        r"""
        Return the contraction of the two forms with respect to the symplectic form.

        The symplectic form `\omega` gives rise to a bilinear form,
        also denoted by `\omega` on the space of `1`-forms by

        .. MATH::
            \omega(\alpha, \beta) = \omega(\alpha^\sharp, \beta^\sharp),

        where `\alpha^\sharp` is the dual of `\alpha` with respect to `\omega`, see
        :meth:`~sage.manifolds.differentiable.tensor_field.TensorField.up`.
        This bilinear form induces a bilinear form on the space of all forms determined
        by its value on decomposable elements as:

        .. MATH::
            \omega(\alpha_1 \wedge \ldots \wedge\alpha_p, \beta_1 \wedge \ldots \wedge\beta_p)
            = det(\omega(\alpha_i, \beta_j)).

        INPUT:

        - ``first`` -- a `p`-form `\alpha`
        - ``second`` -- a `p`-form `\beta`

        OUTPUT:

        - the scalar field `\omega(\alpha, \beta)`

        EXAMPLES:

            sage: M = manifolds.StandardSymplecticSpace(2)
            sage: omega = M.symplectic_form()
            sage: a = M.one_form(1, 0, name='a')
            sage: b = M.one_form(0, 1, name='b')
            sage: omega.on_forms(a, b).display()
            R2 → ℝ
            (q, p) ↦ -1
        """
        from sage.arith.misc import factorial

        if first.degree() != second.degree():
            raise ValueError("the two forms must have the same degree")

        all_positions = range(first.degree())
        return first.contract(
            *all_positions, second.up(self), *all_positions
        ) / factorial(first.degree())


class SymplecticFormParal(SymplecticForm, DiffFormParal):
    r"""
    A symplectic form on a parallelizable manifold.

    .. NOTE::

        If `M` is not parallelizable, the class :class:`SymplecticForm`
        should be used instead.

    INPUT:

    - ``manifold`` -- module `\mathfrak{X}(M)` of vector fields on the
      manifold `M`, or the manifold `M` itself
    - ``name`` -- (default: ``omega``) name given to the symplectic form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      symplectic form; if ``None``, it is formed from ``name``

    EXAMPLES:

    Standard symplectic form on `\RR^2`::

        sage: M.<q, p> = EuclideanSpace(name='R2', latex_name=r"\mathbb{R}^2")
        sage: omega = M.symplectic_form(name='omega', latex_name=r'\omega')
        sage: omega
        Symplectic form omega on the Euclidean plane R2
        sage: omega.set_comp()[1,2] = -1
        sage: omega.display()
        omega = -dq∧dp
    """
    _poisson: TensorFieldParal

    def __init__(
        self,
        manifold: Union[VectorFieldModule, DifferentiableManifold],
        name: Optional[str],
        latex_name: Optional[str] = None,
    ):
        r"""
        Construct a symplectic form.

        TESTS::

            sage: from sage.manifolds.differentiable.symplectic_form import SymplecticFormParal
            sage: M = EuclideanSpace(2)
            sage: omega = SymplecticFormParal(M, name='omega', latex_name=r'\omega')
            sage: omega
            Symplectic form omega on the Euclidean plane E^2
        """
        try:
            vector_field_module = manifold.vector_field_module()
        except AttributeError:
            vector_field_module = manifold

        if name is None:
            name = "omega"

        DiffFormParal.__init__(
            self, vector_field_module, 2, name=name, latex_name=latex_name
        )

        # Check that manifold is even dimensional
        dim = self._ambient_domain.dimension()
        if dim % 2 == 1:
            raise ValueError(
                f"the dimension of the manifold must be even but it is {dim}"
            )
        self._dim_half = dim // 2

        # Initialization of derived quantities
        SymplecticFormParal._init_derived(self)

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: from sage.manifolds.differentiable.symplectic_form import SymplecticFormParal
            sage: M = EuclideanSpace(2)
            sage: omega = SymplecticFormParal(M, 'omega', r'\omega')
            sage: omega._init_derived()
        """
        # Initialization of quantities pertaining to mother classes
        DiffFormParal._init_derived(self)
        SymplecticForm._init_derived(self)

    def _del_derived(self, del_restrictions: bool = True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- boolean (default: ``True``); determines whether the
          restrictions of ``self`` to subdomains are deleted

        TESTS::

            sage: from sage.manifolds.differentiable.symplectic_form import SymplecticFormParal
            sage: M.<q, p> = EuclideanSpace(2, "R2", r"\mathbb{R}^2", symbols=r"q:q p:p")
            sage: omega = SymplecticFormParal(M, 'omega', r'\omega')
            sage: omega._del_derived(del_restrictions=False)
            sage: omega._del_derived()
        """
        # Delete derived quantities from mother classes
        DiffFormParal._del_derived(self, del_restrictions=del_restrictions)

        # Clear the Poisson tensor
        if self._poisson is not None:
            self._poisson._components.clear()
            self._poisson._del_derived()

        SymplecticForm._del_derived(self)

    def restrict(
        self, subdomain: DifferentiableManifold, dest_map: Optional[DiffMap] = None
    ) -> SymplecticFormParal:
        r"""
        Return the restriction of the symplectic form to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of the symplectic form's domain
        - ``dest_map`` -- (default: ``None``) smooth destination map
          `\Phi:\ U \rightarrow V`, where `V` is a subdomain of the symplectic form's domain.
          If ``None``, the restriction of the initial vector field module is used.

        OUTPUT: the restricted symplectic form

        EXAMPLES:

        Restriction of the standard symplectic form on `\RR^2` to the upper half plane::

            sage: from sage.manifolds.differentiable.symplectic_form import SymplecticFormParal
            sage: M = EuclideanSpace(2, "R2", r"\mathbb{R}^2", symbols=r"q:q p:p")
            sage: X.<q, p> = M.chart()
            sage: omega = SymplecticFormParal(M, 'omega', r'\omega')
            sage: omega[1,2] = -1
            sage: U = M.open_subset('U', coord_def={X: q>0})
            sage: omegaU = omega.restrict(U); omegaU
            Symplectic form omega on the Open subset U of the Euclidean plane R2
            sage: omegaU.display()
            omega = -dq∧dp
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            # Construct the restriction at the tensor field level:
            resu = DiffFormParal.restrict(self, subdomain, dest_map=dest_map)

            self._restrictions[subdomain] = SymplecticFormParal.wrap(resu)
        return self._restrictions[subdomain]

    def poisson(
        self, expansion_symbol: Optional[Expression] = None, order: int = 1
    ) -> TensorFieldParal:
        r"""
        Return the Poisson tensor associated with the symplectic form.

        INPUT:

        - ``expansion_symbol`` -- (default: ``None``) symbolic variable; if
          specified, the inverse will be expanded in power series with respect
          to this variable (around its zero value)
        - ``order`` -- integer (default: 1); the order of the expansion
          if ``expansion_symbol`` is not ``None``; the *order* is defined as
          the degree of the polynomial representing the truncated power series
          in ``expansion_symbol``; currently only first order inverse is
          supported

        If ``expansion_symbol`` is set, then the zeroth order symplectic form must be
        invertible. Moreover, subsequent calls to this method will return
        a cached value, even when called with the default value (to enable
        computation of derived quantities). To reset, use :meth:`_del_derived`.

        OUTPUT:

        - the Poisson tensor, , as an instance of
          :meth:`~sage.manifolds.differentiable.poisson_tensor.PoissonTensorFieldParal`

        EXAMPLES:

        Poisson tensor of `2`-dimensional symplectic vector space::

            sage: from sage.manifolds.differentiable.symplectic_form import SymplecticFormParal
            sage: M.<q, p> = EuclideanSpace(2, "R2", r"\mathbb{R}^2", symbols=r"q:q p:p")
            sage: omega = SymplecticFormParal(M, 'omega', r'\omega')
            sage: omega[1,2] = -1
            sage: poisson = omega.poisson(); poisson
            2-vector field poisson_omega on the Euclidean plane R2
            sage: poisson.display()
            poisson_omega = -e_q∧e_p
        """
        super().poisson()

        if expansion_symbol is not None:
            if (
                self._poisson is not None
                and bool(self._poisson._components)
                and list(self._poisson._components.values())[0][0, 0]._expansion_symbol
                == expansion_symbol
                and list(self._poisson._components.values())[0][0, 0]._order == order
            ):
                return self._poisson

            if order != 1:
                raise NotImplementedError("only first order inverse is implemented")
            decompo = self.series_expansion(expansion_symbol, order)
            g0 = decompo[0]
            g1 = decompo[1]

            g0m = self._new_instance()  # needed because only metrics have
            g0m.set_comp()[:] = g0[:]  # an "inverse" method.

            contraction = g1.contract(0, g0m.inverse(), 0)
            contraction = contraction.contract(1, g0m.inverse(), 1)
            self._poisson = -(g0m.inverse() - expansion_symbol * contraction)
            self._poisson.set_calc_order(expansion_symbol, order)
            return self._poisson

        from sage.matrix.constructor import matrix
        from sage.tensor.modules.comp import CompFullyAntiSym

        # Is the inverse metric up to date ?
        for frame in self._components:
            if frame not in self._poisson._components:
                # the computation is necessary
                fmodule = self._fmodule
                si = fmodule._sindex
                nsi = fmodule.rank() + si
                dom = self._domain
                comp_poisson = CompFullyAntiSym(
                    fmodule._ring,
                    frame,
                    2,
                    start_index=si,
                    output_formatter=fmodule._output_formatter,
                )
                comp_poisson_scal = (
                    {}
                )  # dict. of scalars representing the components of the poisson tensor (keys: comp. indices)
                for i in fmodule.irange():
                    for j in range(i, nsi):  # symmetry taken into account
                        comp_poisson_scal[(i, j)] = dom.scalar_field()
                for chart in dom.top_charts():
                    # TODO: do the computation without the 'SR' enforcement
                    try:
                        self_matrix = matrix(
                            [
                                [
                                    self.comp(frame)[i, j, chart].expr(method='SR')
                                    for j in fmodule.irange()
                                ]
                                for i in fmodule.irange()
                            ]
                        )
                        self_matrix_inv = self_matrix.inverse()
                    except (KeyError, ValueError):
                        continue
                    for i in fmodule.irange():
                        for j in range(i, nsi):
                            val = chart.simplify(
                                -self_matrix_inv[i - si, j - si], method="SR"
                            )
                            comp_poisson_scal[(i, j)].add_expr(val, chart=chart)
                for i in range(si, nsi):
                    for j in range(i, nsi):
                        comp_poisson[i, j] = comp_poisson_scal[(i, j)]
                self._poisson._components[frame] = comp_poisson
        return self._poisson


# ****************************************************************************************
