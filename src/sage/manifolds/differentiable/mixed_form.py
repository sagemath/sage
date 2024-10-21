# sage_setup: distribution = sagemath-symbolics
r"""
Mixed Differential Forms

Let `M` and `N` be differentiable manifolds and `\varphi : M \longrightarrow N`
a differentiable map. A *mixed differential form along* `\varphi` is an element
of the graded algebra represented by
:class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`.
Its homogeneous components consist of differential forms along `\varphi`. Mixed
forms are useful to represent characteristic classes and perform computations
of such.

AUTHORS:

- Michael Jung (2019) : initial version
"""
# *****************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.element import AlgebraElement, ModuleElementWithMutability
from sage.rings.integer import Integer


class MixedForm(AlgebraElement, ModuleElementWithMutability):
    r"""
    An instance of this class is a mixed form along some differentiable map
    `\varphi: M \to N` between two differentiable manifolds `M` and `N`. More
    precisely, a mixed form `a` along `\varphi: M \to N` can be considered as a
    differentiable map

    .. MATH::

        a: M \longrightarrow \bigoplus^n_{k=0} T^{(0,k)}N,

    where `T^{(0,k)}` denotes the tensor bundle of type `(0,k)`, `\bigoplus`
    the Whitney sum and `n` the dimension of `N`, such that

    .. MATH::

        \forall x\in M, \quad a(x) \in \bigoplus^n_{k=0} \Lambda^k\left( T_{\varphi(x)}^* N \right),

    where `\Lambda^k(T^*_{\varphi(x)} N)` is the `k`-th exterior power of the
    dual of the tangent space `T_{\varphi(x)} N`. Thus, a mixed differential
    form `a` consists of homogeneous components `a_i`, `i=0,1, \dots, n`, where
    the `i`-th homogeneous component represents a differential form of
    degree `i`.

    The standard case of a mixed form *on* `M` corresponds to `M=N` with
    `\varphi = \mathrm{Id}_M`.

    INPUT:

    - ``parent`` -- graded algebra of mixed forms represented by
      :class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`
      where the mixed form ``self`` shall belong to
    - ``comp`` -- (default: ``None``) homogeneous components of the mixed form
      as a list; if none is provided, the components are set to innocent unnamed
      differential forms
    - ``name`` -- (default: ``None``) name given to the mixed form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      mixed form; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Initialize a mixed form on a 2-dimensional parallelizable differentiable
    manifold::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: e_xy = c_xy.frame()
        sage: A = M.mixed_form(name='A'); A
        Mixed differential form A on the 2-dimensional differentiable manifold M
        sage: A.parent()
        Graded algebra Omega^*(M) of mixed differential forms on the
         2-dimensional differentiable manifold M

    The default way to specify the `i`-th homogeneous component
    of a mixed form is by accessing it via ``A[i]`` or using :meth:`set_comp`::

        sage: A = M.mixed_form(name='A')
        sage: A[0].set_expr(x) # scalar field
        sage: A.set_comp(1)[0] = y*x
        sage: A.set_comp(2)[0,1] = y^2*x
        sage: A.display() # display names
        A = A_0 + A_1 + A_2
        sage: A.display_expansion() # display expansion in basis
        A = x + x*y dx + x*y^2 dx∧dy

    Another way to define the homogeneous components is using predefined
    differential forms::

        sage: f = M.scalar_field(x, name='f'); f
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: omega = M.diff_form(1, name='omega'); omega
        1-form omega on the 2-dimensional differentiable manifold M
        sage: omega[e_xy,0] = y*x; omega.display()
        omega = x*y dx
        sage: eta = M.diff_form(2, name='eta'); eta
        2-form eta on the 2-dimensional differentiable manifold M
        sage: eta[e_xy,0,1] = y^2*x; eta.display()
        eta = x*y^2 dx∧dy

    The components of a mixed form ``B`` can then be set as follows::

        sage: B = M.mixed_form(name='B')
        sage: B[:] = [f, omega, eta]; B.display() # display names
        B = f + omega + eta
        sage: B.display_expansion() # display in coordinates
        B = x + x*y dx + x*y^2 dx∧dy
        sage: B[0]
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: B[1]
        1-form omega on the 2-dimensional differentiable manifold M
        sage: B[2]
        2-form eta on the 2-dimensional differentiable manifold M

    As we can see, the names are applied. However note that the differential
    forms are different instances::

        sage: f is B[0]
        False

    Alternatively, the components can be determined from scratch::

        sage: B = M.mixed_form([f, omega, eta], name='B')
        sage: B.display()
        B = f + omega + eta

    Mixed forms are elements of an algebra so they can be added, and multiplied
    via the wedge product::

        sage: C = x*A; C
        Mixed differential form x∧A on the 2-dimensional differentiable
         manifold M
        sage: C.display_expansion()
        x∧A = x^2 + x^2*y dx + x^2*y^2 dx∧dy
        sage: D = A+C; D
        Mixed differential form A+x∧A on the 2-dimensional differentiable
         manifold M
        sage: D.display_expansion()
        A+x∧A = x^2 + x + (x^2 + x)*y dx + (x^2 + x)*y^2 dx∧dy
        sage: E = A*C; E
        Mixed differential form A∧(x∧A) on the 2-dimensional differentiable
         manifold M
        sage: E.display_expansion()
        A∧(x∧A) = x^3 + 2*x^3*y dx + 2*x^3*y^2 dx∧dy

    Coercions are fully implemented::

        sage: F = omega*A
        sage: F.display_expansion()
        omega∧A = x^2*y dx
        sage: G = omega+A
        sage: G.display_expansion()
        omega+A = x + 2*x*y dx + x*y^2 dx∧dy

    Moreover, it is possible to compute the exterior derivative of a
    mixed form::

        sage: dA = A.exterior_derivative(); dA.display()
        dA = zero + dA_0 + dA_1
        sage: dA.display_expansion()
        dA = dx - x dx∧dy

    Initialize a mixed form on a 2-dimensional non-parallelizable differentiable
    manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
        ....:                   intersection_name='W', restrictions1= x>0,
        ....:                   restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: e_xy = c_xy.frame(); e_uv = c_uv.frame() # define frames
        sage: A = M.mixed_form(name='A')
        sage: A[0].set_expr(x, c_xy)
        sage: A[0].display()
        A_0: M → ℝ
        on U: (x, y) ↦ x
        on W: (u, v) ↦ 1/2*u + 1/2*v
        sage: A[1][0] = y*x; A[1].display(e_xy)
        A_1 = x*y dx
        sage: A[2][e_uv,0,1] = u*v^2; A[2].display(e_uv)
        A_2 = u*v^2 du∧dv
        sage: A.add_comp_by_continuation(e_uv, W, c_uv)
        sage: A.display_expansion(e_uv)
        A = 1/2*u + 1/2*v + (1/8*u^2 - 1/8*v^2) du + (1/8*u^2 - 1/8*v^2) dv + u*v^2 du∧dv
        sage: A.add_comp_by_continuation(e_xy, W, c_xy)
        sage: A.display_expansion(e_xy)
        A = x + x*y dx + (-2*x^3 + 2*x^2*y + 2*x*y^2 - 2*y^3) dx∧dy

    Since zero and one are special elements, their components cannot be
    changed::

        sage: z = M.mixed_form_algebra().zero()
        sage: z[0] = 1
        Traceback (most recent call last):
        ...
        ValueError: the components of an immutable element cannot be changed
        sage: one = M.mixed_form_algebra().one()
        sage: one[0] = 0
        Traceback (most recent call last):
        ...
        ValueError: the components of an immutable element cannot be changed
    """
    def __init__(self, parent, name=None, latex_name=None):
        r"""
        Construct a mixed form.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                   intersection_name='W', restrictions1= x>0,
            ....:                   restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[e_xy,0] = y*x
            sage: omega.add_comp_by_continuation(e_uv, W, c_uv)
            sage: eta = M.diff_form(2, name='eta')
            sage: eta[e_uv,0,1] = u*v^2
            sage: eta.add_comp_by_continuation(e_xy, W, c_xy)
            sage: A = M.mixed_form_algebra()
            sage: F = A([x, omega, eta], name='F')
            sage: TestSuite(F).run(skip='_test_pickling')
        """
        if parent is None:
            raise ValueError("a parent must be provided")
        # Add this element (instance) to the set of parent algebra:
        AlgebraElement.__init__(self, parent)
        # Get the underlying vector field module:
        vmodule = parent._vmodule
        # Define attributes:
        self._vmodule = vmodule
        self._dest_map = vmodule._dest_map
        self._domain = vmodule._domain
        self._ambient_domain = vmodule._ambient_domain
        self._max_deg = vmodule._ambient_domain.dim()
        self._is_zero = False  # a priori, may be changed below or via
        # method __bool__()
        self._comp = None  # initialized on demand; see _init_comp
        # Set names:
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

    def _init_comp(self):
        r"""
        Initialize homogeneous components of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: A = M.mixed_form(name='A')
            sage: A._comp is None
            True
            sage: A._init_comp()
            sage: A._comp
            [Scalar field A_0 on the 2-dimensional differentiable manifold M,
             1-form A_1 on the 2-dimensional differentiable manifold M,
             2-form A_2 on the 2-dimensional differentiable manifold M]
        """
        self._comp = []
        for i in self.irange():
            comp_name, comp_latex_name = None, None
            if self._name is not None:
                comp_name = f"{self._name}_{i}"
            if self._latex_name is not None:
                comp_latex_name = '{' + self._latex_name + '}_{' + str(i) + '}'
            diff_form = self._domain.diff_form
            self._comp.append(diff_form(i, name=comp_name,
                                        latex_name=comp_latex_name))

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: F = M.mixed_form(name='F')
            sage: F._repr_()
            'Mixed differential form F on the 3-dimensional differentiable
             manifold M'
            sage: repr(F)  # indirect doctest
            'Mixed differential form F on the 3-dimensional differentiable
             manifold M'

        Check whether :issue:`31784` is fixed::

            sage: E3 = EuclideanSpace(3)
            sage: S2 = E3.sphere()
            sage: iota = S2.embedding()
            sage: Omega = S2.mixed_form_algebra(dest_map=iota)
            sage: Omega(1)
            Mixed differential form one along the 2-sphere S^2 of radius 1
             smoothly embedded in the Euclidean space E^3 with values on the
             Euclidean space E^3 via the map iota
        """
        desc = "Mixed differential form "
        if self._name is not None:
            desc += self._name + " "
        if self._dest_map is self._domain.identity_map():
            desc += f"on the {self._domain}"
        else:
            desc += f"along the {self._domain} with values on the {self._ambient_domain} "
            if self._dest_map._name is None:
                dm_name = "unnamed map"
            else:
                dm_name = self._dest_map._name
            desc += "via the map " + dm_name
        return desc

    def _latex_(self):
        r"""
        Return a LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(3, 'M', latex_name=r'\mathcal{M}')
            sage: omega = M.mixed_form(name='omega', latex_name=r'\omega')
            sage: omega._latex_()
            '\\omega'
            sage: latex(omega)  # indirect doctest
            \omega
        """
        if self._name is None:
            return r'\text{' + repr(self) + r'}'
        else:
            return self._latex_name

    def _new_instance(self, name=None, latex_name=None):
        r"""
        Return a new instance of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: F = M.mixed_form(name='F')
            sage: F1 = F._new_instance(); F1
            Mixed differential form on the 2-dimensional differentiable
             manifold M
            sage: type(F1) == type(F)
            True
            sage: F1.parent() is F.parent()
            True
        """
        return type(self)(self.parent(), name=name, latex_name=latex_name)

    def display_expansion(self, frame=None, chart=None, from_chart=None):
        r"""
        Display the expansion in a particular basis and chart of mixed forms.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to
          which the mixed form is expanded; if ``None``, only the names
          of the components are displayed
        - ``chart`` -- (default: ``None``) chart with respect to which the
          components of the mixed form in the selected frame are expressed;
          if ``None``, the default chart of the vector frame domain is assumed

        EXAMPLES:

        Display the expansion of a mixed form on a 2-dimensional
        non-parallelizable differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x-y, x+y),
            ....:                   intersection_name='W', restrictions1= x>0,
            ....:                   restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame() # define frames
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[e_xy,0] = x; omega.display(e_xy)
            omega = x dx
            sage: omega.add_comp_by_continuation(e_uv, W, c_uv) # continuation onto M
            sage: eta = M.diff_form(2, name='eta')
            sage: eta[e_uv,0,1] = u*v; eta.display(e_uv)
            eta = u*v du∧dv
            sage: eta.add_comp_by_continuation(e_xy, W, c_xy) # continuation onto M
            sage: F = M.mixed_form([0, omega, eta], name='F'); F
            Mixed differential form F on the 2-dimensional differentiable
             manifold M
            sage: F.display() # display names of homogeneous components
            F = zero + omega + eta
            sage: F.display_expansion(e_uv)
            F = (1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv + u*v du∧dv
            sage: F.display_expansion(e_xy)
            F = x dx + (2*x^2 - 2*y^2) dx∧dy
        """
        from sage.misc.latex import latex
        from sage.typeset.unicode_characters import unicode_wedge
        from sage.tensor.modules.format_utilities import (is_atomic,
                                                          FormattedExpansion)
        # In case, no frame is given:
        if frame is None:
            frame = self._domain._def_frame
        # In case, no chart is given:
        if chart is None:
            chart = frame._chart
        # Check names:
        if self._name is not None:
            resu_txt = self._name + " = "
        else:
            resu_txt = ""
        if self._latex_name is not None:
            resu_latex = self._latex_name + r" = "
        else:
            resu_latex = ""
        # Get terms
        terms_txt = []
        terms_latex = []
        # Scalar field term:
        if not self[0].is_trivial_zero():
            terms_txt.append(repr(self[0].expr(chart, from_chart)))
            terms_latex.append(latex(self[0].expr(chart, from_chart)))
        # Differential form terms:
        for j in self.irange(1):
            rst = self[j].restrict(frame._domain, dest_map=frame._dest_map)
            basis, format_spec = rst._preparse_display(basis=frame,
                                                       format_spec=chart)
            cobasis = basis.dual_basis()
            comp = rst.comp(basis)
            for ind in comp.non_redundant_index_generator():
                ind_arg = ind + (format_spec,)
                coef = comp[ind_arg]
                # Check whether the coefficient is zero, preferably via
                # the fast method is_trivial_zero():
                if hasattr(coef, 'is_trivial_zero'):
                    zero_coef = coef.is_trivial_zero()
                else:
                    zero_coef = coef == 0
                if not zero_coef:
                    bases_txt = []
                    bases_latex = []
                    for k in range(rst._tensor_rank):
                        bases_txt.append(cobasis[ind[k]]._name)
                        bases_latex.append(latex(cobasis[ind[k]]))
                    basis_term_txt = unicode_wedge.join(bases_txt)
                    basis_term_latex = r"\wedge ".join(bases_latex)
                    coef_txt = repr(coef)
                    if coef_txt == "1":
                        terms_txt.append(basis_term_txt)
                        terms_latex.append(basis_term_latex)
                    elif coef_txt == "-1":
                        terms_txt.append("-" + basis_term_txt)
                        terms_latex.append("-" + basis_term_latex)
                    else:
                        coef_latex = latex(coef)
                        if is_atomic(coef_txt):
                            terms_txt.append(coef_txt + " " + basis_term_txt)
                        else:
                            terms_txt.append("(" + coef_txt + ") " +
                                             basis_term_txt)
                        if is_atomic(coef_latex):
                            terms_latex.append(coef_latex + basis_term_latex)
                        else:
                            terms_latex.append(r"\left(" + coef_latex +
                                               r"\right)" + basis_term_latex)
        if not terms_txt:
            resu_txt += "0"
        else:
            resu_txt += terms_txt[0]
            for term in terms_txt[1:]:
                if term[0] == "-":
                    resu_txt += " - " + term[1:]
                else:
                    resu_txt += " + " + term
        if not terms_latex:
            resu_latex += "0"
        else:
            resu_latex += terms_latex[0]
            for term in terms_latex[1:]:
                if term[0] == "-":
                    resu_latex += term
                else:
                    resu_latex += "+" + term
        return FormattedExpansion(resu_txt, resu_latex)

    disp_exp = display_expansion
    display_exp = display_expansion

    def display(self):
        r"""
        Display the homogeneous components of the mixed form.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: f = M.scalar_field(name='f')
            sage: omega = M.diff_form(1, name='omega')
            sage: eta = M.diff_form(2, name='eta')
            sage: F = M.mixed_form([f, omega, eta], name='F'); F
            Mixed differential form F on the 2-dimensional differentiable
             manifold M
            sage: F.display() # display names of homogeneous components
            F = f + omega + eta
        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        # Mixed form name:
        if self._name is not None:
            resu_txt = self._name + " = "
        else:
            resu_txt = ""
        if self._latex_name is not None:
            resu_latex = self._latex_name + r" = "
        else:
            resu_latex = ""
        # Scalar field:
        if self[0]._name is None:
            resu_txt += "(unnamed scalar field) "
        else:
            resu_txt += self[0]._name
        if self[0]._latex_name is None:
            resu_latex += r"\text{(unnamed scalar field)}"
        else:
            resu_latex += latex(self[0])
        # Differential forms:
        for j in self.irange(1):
            if self[j]._name is None:
                resu_txt += " + (unnamed " + str(j) + "-form)"
            else:
                resu_txt += " + " + self[j]._name
            if self[j]._latex_name is None:
                resu_latex += r"+\text{(unnamed " + str(j) + r"-form)}"
            else:
                resu_latex += r"+" + latex(self[j])
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display

    def set_name(self, name=None, latex_name=None, apply_to_comp=True):
        r"""
        Redefine the string and LaTeX representation of the object.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the mixed form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          mixed form; if none is provided, the LaTeX symbol is set to ``name``
        - ``apply_to_comp`` -- boolean (default: ``True``); if ``True`` all homogeneous
          components will be renamed accordingly. If ``False`` only the mixed
          form will be renamed

        EXAMPLES:

        Rename a mixed form::

            sage: M = Manifold(4, 'M')
            sage: F = M.mixed_form(name='dummy', latex_name=r'\ugly'); F
            Mixed differential form dummy on the 4-dimensional differentiable
             manifold M
            sage: latex(F)
            \ugly
            sage: F.set_name(name='F', latex_name=r'\mathcal{F}'); F
            Mixed differential form F on the 4-dimensional differentiable
             manifold M
            sage: latex(F)
            \mathcal{F}

        If not stated otherwise, all homogeneous components are renamed
        accordingly::

            sage: F.display()
            F = F_0 + F_1 + F_2 + F_3 + F_4

        Setting the argument ``set_all`` to ``False`` prevents the renaming
        in the homogeneous components::

            sage: F.set_name(name='eta', latex_name=r'\eta', apply_to_comp=False)
            sage: F.display()
            eta = F_0 + F_1 + F_2 + F_3 + F_4

        To rename a homogeneous component individually, we simply access the
        homogeneous component and use its
        :meth:`~sage.manifolds.differentiable.tensorfield.set_name` method::

            sage: F[0].set_name(name='g'); F.display()
            eta = g + F_1 + F_2 + F_3 + F_4
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
        if apply_to_comp:
            for i in self.irange():
                comp_name, comp_latex_name = None, None
                if self._name is not None:
                    comp_name = f"{self._name}_{i}"
                if self._latex_name is not None:
                    comp_latex_name = '{' + self._latex_name + '}_{' + str(i) + '}'
                self[i].set_name(name=comp_name, latex_name=comp_latex_name)

    def __bool__(self):
        r"""
        Return ``True`` if ``self`` is nonzero and ``False`` otherwise.

        This method is called by :meth:`is_zero`.

        EXAMPLES:

        Mixed form defined by parts on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x,y> = U.chart()
            sage: V = M.open_subset('V')
            sage: c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: F = M.mixed_form(name='F')
            sage: FU = U.mixed_form(name='F')
            sage: FV = V.mixed_form(name='F')
            sage: FU[:] = [0,0,0]
            sage: FV[:] = [0,0,0]
            sage: F.set_restriction(FU)
            sage: F.set_restriction(FV)
            sage: bool(F)
            False
            sage: F.is_zero()  # indirect doctest
            True
            sage: FU[0] = 1
            sage: F.set_restriction(FU)
            sage: bool(F)
            True
            sage: F.is_zero()  # indirect doctest
            False
        """
        if self._is_zero:
            return False
        if any(bool(form) for form in self):
            self._is_zero = False
            return True
        self._is_zero = True
        return False

    def _richcmp_(self, other, op):
        r"""
        Comparison method for sage objects.

        INPUT:

        - ``other`` -- a mixed form
        - ``op`` -- comparison operator for which ``self`` and ``other`` shall
          be compared with

        TESTS::

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
            sage: f = M.scalar_field(x, name='f')
            sage: f.add_expr_by_continuation(c_uv, W)
            sage: eta = M.diff_form(1, name='eta')
            sage: eta[e_xy,0] = x+y
            sage: eta.add_comp_by_continuation(e_uv, W, c_uv)
            sage: F = M.mixed_form([f, eta, 0])
            sage: F == F
            True
            sage: F == F.copy()
            True
            sage: G = M.mixed_form()
            sage: G.set_restriction(F.restrict(U))
            sage: F == G  # False since G has not been defined on V
            False
            sage: G.set_restriction(F.restrict(V))
            sage: F == G  # ``True`` now
            True
            sage: H = M.mixed_form([f, 0, 0])
            sage: F != H  # this is fixed by issue #30108
            True
            sage: F.parent().zero() == 0
            True
        """
        from sage.structure.richcmp import op_NE, op_EQ
        if op == op_NE:
            return not self == other
        elif op == op_EQ:
            # Compare all elements separately:
            return all(self[j] == other[j] for j in self.irange())
        # Fall back on default implementation:
        return super()._richcmp_(self, other, op)

    def _add_(self, other):
        r"""
        Addition of two mixed forms.

        INPUT:

        - ``other`` -- a mixed form, in the same algebra as ``self``

        OUTPUT:

        - the mixed form resulting from the addition of ``self``
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
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: f = M.scalar_field(x, name='f')
            sage: f.add_expr_by_continuation(c_uv, W) # continuation onto M
            sage: a = M.diff_form(1, name='a')
            sage: a[e_xy,0] = x
            sage: a.add_comp_by_continuation(e_uv, W, c_uv) # continuation onto M
            sage: A = M.mixed_form([f, a, 0], name='A')
            sage: g = M.scalar_field(u, name='g', chart=c_uv)
            sage: g.add_expr_by_continuation(c_xy, W) # continuation onto M
            sage: b = M.diff_form(1, name='b')
            sage: b[e_uv,1] = v
            sage: b.add_comp_by_continuation(e_xy, W, c_xy) # continuation onto M
            sage: B = M.mixed_form([g, b, 0], name='B')
            sage: C = A._add_(B); C
            Mixed differential form A+B on the 2-dimensional differentiable
             manifold M
            sage: A.display_expansion(e_uv)
            A = 1/2*u + 1/2*v + (1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv
            sage: B.display_expansion(e_xy)
            B = x + y + (x - y) dx + (-x + y) dy
            sage: C.display_expansion(e_xy)
            A+B = 2*x + y + (2*x - y) dx + (-x + y) dy
            sage: C.display_expansion(e_uv)
            A+B = 3/2*u + 1/2*v + (1/4*u + 1/4*v) du + (1/4*u + 5/4*v) dv
            sage: C == A + B  # indirect doctest
            True
            sage: Z = A.parent().zero(); Z
            Mixed differential form zero on the 2-dimensional differentiable
             manifold M
            sage: A._add_(Z) == A
            True
            sage: Z._add_(A) == A
            True
        """
        # Case zero:
        if self._is_zero:
            return other
        if other._is_zero:
            return self
        # Generic case:
        resu = self._new_instance()
        resu._comp = [self[j] + other[j] for j in self.irange()]
        # Compose name:
        if self._name is not None and other._name is not None:
            resu._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            resu._latex_name = self._latex_name + '+' + other._latex_name
        return resu

    def _sub_(self, other):
        r"""
        Subtraction of two mixed forms.

        INPUT:

        - ``other`` -- a mixed form, in the same algebra as ``self``

        OUTPUT:

        - the mixed form resulting from the subtraction of ``self``
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
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: f = M.scalar_field(x, name='f')
            sage: f.add_expr_by_continuation(c_uv, W) # continuation onto M
            sage: a = M.diff_form(1, name='a')
            sage: a[e_xy,0] = x
            sage: a.add_comp_by_continuation(e_uv, W, c_uv) # continuation onto M
            sage: A = M.mixed_form([f, a, 0], name='A')
            sage: g = M.scalar_field(u, name='g', chart=c_uv)
            sage: g.add_expr_by_continuation(c_xy, W) # continuation onto M
            sage: b = M.diff_form(1, name='b')
            sage: b[e_uv,1] = v
            sage: b.add_comp_by_continuation(e_xy, W, c_xy) # continuation onto M
            sage: B = M.mixed_form([g, b, 0], name='B')
            sage: C = A._sub_(B); C
            Mixed differential form A-B on the 2-dimensional differentiable
             manifold M
            sage: A.display_expansion(e_uv)
            A = 1/2*u + 1/2*v + (1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv
            sage: B.display_expansion(e_xy)
            B = x + y + (x - y) dx + (-x + y) dy
            sage: C.display_expansion(e_xy)
            A-B = -y + y dx + (x - y) dy
            sage: C.display_expansion(e_uv)
            A-B = -1/2*u + 1/2*v + (1/4*u + 1/4*v) du + (1/4*u - 3/4*v) dv
            sage: C == A - B  # indirect doctest
            True
            sage: Z = A.parent().zero(); Z
            Mixed differential form zero on the 2-dimensional differentiable
             manifold M
            sage: A._sub_(Z) == A
            True
            sage: Z._sub_(A) == -A
            True
        """
        # Case zero:
        if self._is_zero:
            return -other
        if other._is_zero:
            return self
        # Generic case:
        resu = self._new_instance()
        resu._comp = [self[j] - other[j] for j in self.irange()]
        # Compose name:
        from sage.tensor.modules.format_utilities import is_atomic

        if self._name is not None and other._name is not None:
            sname = self._name
            oname = other._name
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            resu._name = sname + '-' + oname
        if self._latex_name is not None and other._latex_name is not None:
            slname = self._latex_name
            olname = other._latex_name
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            resu._latex_name = slname + '-' + olname
        return resu

    def wedge(self, other):
        r"""
        Wedge product on the graded algebra of mixed forms.

        More precisely, the wedge product is a bilinear map

        .. MATH::

            \wedge: \Omega^k(M,\varphi) \times \Omega^l(M,\varphi) \to \Omega^{k+l}(M,\varphi),

        where `\Omega^k(M,\varphi)` denotes the space of differential forms of
        degree `k` along `\varphi`. By bilinear extension, this induces a map

        .. MATH::

            \wedge: \Omega^*(M,\varphi) \times \Omega^*(M,\varphi) \to \Omega^*(M,\varphi) ``

        and equips `\Omega^*(M,\varphi)` with a multiplication such that it
        becomes a graded algebra.

        INPUT:

        - ``other`` -- mixed form in the same algebra as ``self``

        OUTPUT:

        - the mixed form resulting from the wedge product of ``self``
          with ``other``

        EXAMPLES:

        Initialize a mixed form on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: f = M.scalar_field(x, name='f')
            sage: f.display()
            f: M → ℝ
               (x, y, z) ↦ x
            sage: g = M.scalar_field(y, name='g')
            sage: g.display()
            g: M → ℝ
               (x, y, z) ↦ y
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[0] = x
            sage: omega.display()
            omega = x dx
            sage: eta = M.diff_form(1, name='eta')
            sage: eta[1] = y
            sage: eta.display()
            eta = y dy
            sage: mu = M.diff_form(2, name='mu')
            sage: mu[0,2] = z
            sage: mu.display()
            mu = z dx∧dz
            sage: A = M.mixed_form([f, omega, mu, 0], name='A')
            sage: A.display_expansion()
            A = x + x dx + z dx∧dz
            sage: B = M.mixed_form([g, eta, mu, 0], name='B')
            sage: B.display_expansion()
            B = y + y dy + z dx∧dz

        The wedge product of ``A`` and ``B`` yields::

            sage: C = A.wedge(B); C
            Mixed differential form A∧B on the 3-dimensional differentiable
             manifold M
            sage: C.display_expansion()
            A∧B = x*y + x*y dx + x*y dy + x*y dx∧dy + (x + y)*z dx∧dz - y*z dx∧dy∧dz
            sage: D = B.wedge(A); D # Don't even try, it's not commutative!
            Mixed differential form B∧A on the 3-dimensional differentiable
             manifold M
            sage: D.display_expansion() # I told you so!
            B∧A = x*y + x*y dx + x*y dy - x*y dx∧dy + (x + y)*z dx∧dz - y*z dx∧dy∧dz

        Alternatively, the multiplication symbol can be used::

            sage: A*B
            Mixed differential form A∧B on the 3-dimensional differentiable
             manifold M
            sage: A*B == C
            True

        Yet, the multiplication includes coercions::

            sage: E = x*A; E.display_expansion()
            x∧A = x^2 + x^2 dx + x*z dx∧dz
            sage: F = A*eta; F.display_expansion()
            A∧eta = x*y dy + x*y dx∧dy - y*z dx∧dy∧dz
        """
        # Case zero:
        if self._is_zero or other._is_zero:
            return self.parent().zero()
        # Case one:
        if self is self.parent().one():
            return other
        if other is self.parent().one():
            return self
        # Generic case:
        resu = self._new_instance()
        resu._comp = [sum(self[k].wedge(other[j - k]) for k in range(j + 1))
                      for j in self.irange()]
        # Compose name:
        from sage.typeset.unicode_characters import unicode_wedge
        from sage.tensor.modules.format_utilities import (format_mul_txt,
                                                          format_mul_latex)
        resu._name = format_mul_txt(self._name, unicode_wedge, other._name)
        resu._latex_name = format_mul_latex(self._latex_name, r'\wedge ',
                                            other._latex_name)
        return resu

    _mul_ = wedge

    def _lmul_(self, other):
        r"""
        Scalar multiplication operator: return ``number * self`` or
        ``self * number``.

        INPUT:

        - ``other`` -- an element of the symbolic ring

        OUTPUT:

        - the mixed form resulting from the wedge product of ``self``
          with ``other``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[0] = y*x; omega.display()
            omega = x*y dx
            sage: F = M.mixed_form([0, omega, 0], name='F')
            sage: A = x*F*y; A
            Mixed differential form y∧(x∧F) on the 2-dimensional
             differentiable manifold M
            sage: A.display_expansion()
            y∧(x∧F) = x^2*y^2 dx
        """
        try:
            if other.is_trivial_zero():
                return self.parent().zero()
            if (other - 1).is_trivial_zero():
                return self
        except AttributeError:
            # in case base ring is not SR:
            if other == 0:
                return self.parent().zero()
            if other == 1:
                return self
        resu = self._new_instance()
        resu._comp = [other * form for form in self]
        # Compose name:
        from sage.misc.latex import latex
        from sage.typeset.unicode_characters import unicode_wedge
        from sage.tensor.modules.format_utilities import (format_mul_txt,
                                                          format_mul_latex)
        resu._name = format_mul_txt(repr(other), unicode_wedge, self._name)
        resu._latex_name = format_mul_latex(latex(other), r'\wedge ',
                                            self._latex_name)
        return resu

    @cached_method
    def exterior_derivative(self):
        r"""
        Compute the exterior derivative of ``self``.

        More precisely, the *exterior derivative* on `\Omega^k(M,\varphi)` is a
        linear map

        .. MATH::

            \mathrm{d}_{k} : \Omega^k(M,\varphi) \to \Omega^{k+1}(M,\varphi),

        where `\Omega^k(M,\varphi)` denotes the space of differential forms of
        degree `k` along `\varphi`
        (see :meth:`~sage.manifolds.differentiable.diff_form.DiffForm.exterior_derivative`
        for further information). By linear extension, this induces a map on
        `\Omega^*(M,\varphi)`:

        .. MATH::

            \mathrm{d}: \Omega^*(M,\varphi) \to \Omega^*(M,\varphi).

        OUTPUT:

        - a :class:`MixedForm` representing the exterior
          derivative of the mixed form

        EXAMPLES:

        Exterior derivative of a mixed form on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: f = M.scalar_field(z^2, name='f')
            sage: f.disp()
            f: M → ℝ
                (x, y, z) ↦ z^2
            sage: a = M.diff_form(2, 'a')
            sage: a[1,2], a[1,3], a[2,3] = z+y^2, z+x, x^2
            sage: a.disp()
            a = (y^2 + z) dx∧dy + (x + z) dx∧dz + x^2 dy∧dz
            sage: F = M.mixed_form([f, 0, a, 0], name='F'); F.display()
            F = f + zero + a + zero
            sage: dF = F.exterior_derivative()
            sage: dF.display()
            dF = zero + df + dzero + da
            sage: dF = F.exterior_derivative()
            sage: dF.display_expansion()
            dF = 2*z dz + (2*x + 1) dx∧dy∧dz

        Due to long calculation times, the result is cached::

            sage: F.exterior_derivative() is dF
            True
        """
        resu = self._new_instance()
        resu[0] = self._domain.zero_scalar_field()
        resu[1:] = [self[j].exterior_derivative()
                    for j in range(self._max_deg)]
        # Compose name:
        from sage.tensor.modules.format_utilities import (format_unop_txt,
                                                          format_unop_latex)
        resu._name = format_unop_txt('d', self._name)
        resu._latex_name = format_unop_latex(r'\mathrm{d}', self._latex_name)
        return resu

    derivative = exterior_derivative

    def copy(self, name=None, latex_name=None):
        r"""
        Return an exact copy of ``self``.

        .. NOTE::

            The name and names of the components are not copied.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the copy
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          copy; if none is provided, the LaTeX symbol is set to ``name``

        EXAMPLES:

        Initialize a 2-dimensional manifold and differential forms::

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
            sage: f = M.scalar_field(x, name='f', chart=c_xy)
            sage: f.add_expr_by_continuation(c_uv, W)
            sage: f.display()
            f: M → ℝ
            on U: (x, y) ↦ x
            on V: (u, v) ↦ 1/2*u + 1/2*v
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[e_xy,0] = x
            sage: omega.add_comp_by_continuation(e_uv, W, c_uv)
            sage: omega.display()
            omega = x dx
            sage: A = M.mixed_form([f, omega, 0], name='A'); A.display()
            A = f + omega + zero
            sage: A.display_expansion(e_uv)
            A = 1/2*u + 1/2*v + (1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv

        An exact copy is made. The copy is an entirely new instance and has a
        different name, but has the very same values::

            sage: B = A.copy(); B.display()
            (unnamed scalar field) + (unnamed 1-form) + (unnamed 2-form)
            sage: B.display_expansion(e_uv)
            1/2*u + 1/2*v + (1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv
            sage: A == B
            True
            sage: A is B
            False
        """
        resu = self._new_instance()
        resu._comp = [form.copy() for form in self]
        resu.set_name(name=name, latex_name=latex_name)
        resu._is_zero = self._is_zero  # a priori

        return resu

    def __setitem__(self, index, values):
        r"""
        Set a component with respect to some vector frame.

        - ``index`` -- list of indices; if ``[:]`` is provided, all the
          components are set
        - ``values`` -- the values to be set

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x, name='f')
            sage: a = M.diff_form(1, name='a')
            sage: a[0] = y
            sage: b = M.diff_form(2, name='b')
            sage: b[0,1] = x*y
            sage: A = M.mixed_form([f, 0, 0], name='A'); A.display()
            A = f + zero + zero
            sage: A[1:3] = [a, b]; A.display()
            A = f + a + b
            sage: A.display_expansion()
            A = x + y dx + x*y dx∧dy
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        if isinstance(index, (int, Integer)):
            start, stop, step = index, index + 1, 1
        elif isinstance(index, slice):
            start, stop, step = index.indices(self._max_deg + 1)
        else:
            raise TypeError("index must be int, Integer or slice")
        if isinstance(values, list):
            form_list = values
        else:
            form_list = [values]
        if len(form_list) != len(range(start, stop, step)):
            raise IndexError("either input or index out of range")
        for deg, j in zip(range(start, stop, step), range(len(form_list))):
            dmodule = self._domain.diff_form_module(deg, self._dest_map)
            form = dmodule(form_list[j])
            self[deg].copy_from(form)  # keep the names
            self[deg].set_name(name=form._name, latex_name=form._latex_name)
        self._is_zero = False  # a priori

    def __getitem__(self, deg):
        r"""
        Return a component with respect to some frame.

        INPUT:

        - ``deg`` -- slice of which degrees shall be returned

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(name='f')
            sage: a = M.diff_form(1, name='a')
            sage: b = M.diff_form(2, name='b')
            sage: A = M.mixed_form([f, a, b], name='A'); A.display()
            A = f + a + b
            sage: A.__getitem__(0)
            Scalar field f on the 2-dimensional differentiable manifold M
            sage: A.__getitem__(1)
            1-form a on the 2-dimensional differentiable manifold M
            sage: A.__getitem__(2)
            2-form b on the 2-dimensional differentiable manifold M
            sage: A.__getitem__(slice(0,3,1))
            [Scalar field f on the 2-dimensional differentiable manifold M,
             1-form a on the 2-dimensional differentiable manifold M,
             2-form b on the 2-dimensional differentiable manifold M]
        """
        if self._comp is None:
            self._init_comp()
        return self._comp[deg]

    def set_restriction(self, rst):
        r"""
        Set a (component-wise) restriction of ``self`` to some subdomain.

        INPUT:

        - ``rst`` -- :class:`MixedForm` of the same type as ``self``, defined on
          a subdomain of the domain of ``self``

        EXAMPLES:

        Initialize the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()

        And define some forms on the subset ``U``::

            sage: f = U.scalar_field(x, name='f', chart=c_xy)
            sage: omega = U.diff_form(1, name='omega')
            sage: omega[e_xy,0] = y
            sage: AU = U.mixed_form([f, omega, 0], name='A'); AU
            Mixed differential form A on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: AU.display_expansion(e_xy)
            A = x + y dx

        A mixed form on ``M`` can be specified by some mixed form on a subset::

            sage: A = M.mixed_form(name='A'); A
            Mixed differential form A on the 2-dimensional differentiable
             manifold M
            sage: A.set_restriction(AU)
            sage: A.display_expansion(e_xy)
            A = x + y dx
            sage: A.add_comp_by_continuation(e_uv, W, c_uv)
            sage: A.display_expansion(e_uv)
            A = u/(u^2 + v^2) - (u^2*v - v^3)/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) du - 2*u*v^2/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) dv
            sage: A.restrict(U) == AU
            True
        """
        if not isinstance(rst, MixedForm):
            raise TypeError("the argument must be a mixed form")
        if not rst._domain.is_subset(self._domain):
            raise ValueError("the specified domain is not a subset of "
                             "the domain of definition of the mixed form")
        for j in self.irange():
            self[j].set_restriction(rst[j])
        self._is_zero = False  # a priori

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of ``self`` to some subdomain.

        INPUT:

        - ``subdomain`` --
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
          open subset `U` of the domain of ``self``
        - ``dest_map`` --
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          (default: ``None``); destination map `\Psi:\ U \rightarrow V`,
          where `V` is an open subset of the manifold `N` where the mixed form
          takes it values; if ``None``, the restriction of `\Phi` to `U` is
          used, `\Phi` being the differentiable map `S \rightarrow M` associated
          with the mixed form

        OUTPUT: :class:`MixedForm` representing the restriction

        EXAMPLES:

        Initialize the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()

        And predefine some forms::

            sage: f = M.scalar_field(x^2, name='f', chart=c_xy)
            sage: f.add_expr_by_continuation(c_uv, W)
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[e_xy,0] = y^2
            sage: omega.add_comp_by_continuation(e_uv, W, c_uv)
            sage: eta = M.diff_form(2, name='eta')
            sage: eta[e_xy,0,1] = x^2*y^2
            sage: eta.add_comp_by_continuation(e_uv, W, c_uv)

        Now, a mixed form can be restricted to some subdomain::

            sage: F = M.mixed_form([f, omega, eta], name='F')
            sage: FV = F.restrict(V); FV
            Mixed differential form F on the Open subset V of the 2-dimensional
             differentiable manifold M
            sage: FV[:]
            [Scalar field f on the Open subset V of the 2-dimensional
             differentiable manifold M,
             1-form omega on the Open subset V of the 2-dimensional
             differentiable manifold M,
             2-form eta on the Open subset V of the 2-dimensional
             differentiable manifold M]
            sage: FV.display_expansion(e_uv)
            F = u^2/(u^4 + 2*u^2*v^2 + v^4) - (u^2*v^2 - v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du - 2*u*v^3/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv - u^2*v^2/(u^12 + 6*u^10*v^2 + 15*u^8*v^4 + 20*u^6*v^6 + 15*u^4*v^8 + 6*u^2*v^10 + v^12) du∧dv
        """
        resu = type(self)(subdomain.mixed_form_algebra(dest_map=dest_map),
                          name=self._name, latex_name=self._latex_name)
        resu[0] = self[0].restrict(subdomain)
        resu[1:] = [self[j].restrict(subdomain, dest_map)
                    for j in self.irange(1)]
        return resu

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

        Mixed form defined by differential forms with components on different
        parts of the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: F = M.mixed_form(name='F') # No predefined components, here
            sage: F[0] = M.scalar_field(x, name='f')
            sage: F[1] = M.diff_form(1, {e_xy: [x,0]}, name='omega')
            sage: F[2].set_name(name='eta')
            sage: F[2][e_uv,0,1] = u*v
            sage: F.add_comp_by_continuation(e_uv, W, c_uv)
            sage: F.add_comp_by_continuation(e_xy, W, c_xy) # Now, F is fully defined
            sage: F.display_expansion(e_xy)
            F = x + x dx - x*y/(x^8 + 4*x^6*y^2 + 6*x^4*y^4 + 4*x^2*y^6 + y^8) dx∧dy
            sage: F.display_expansion(e_uv)
            F = u/(u^2 + v^2) - (u^3 - u*v^2)/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) du - 2*u^2*v/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) dv + u*v du∧dv
        """
        if chart is None:
            chart = frame._chart
        self[0].add_expr_by_continuation(chart, subdomain)
        for j in self.irange(1):
            self[j].add_comp_by_continuation(frame, subdomain, chart)

    def irange(self, start=None):
        r"""
        Single index generator.

        INPUT:

        - ``start`` -- (default: ``None``) initial value `i_0` of the index
          between 0 and `n`, where `n` is the manifold's dimension; if none is
          provided, the value 0 is assumed

        OUTPUT:

        - an iterable index, starting from `i_0` and ending at
          `n`, where `n` is the manifold's dimension

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: a = M.mixed_form(name='a')
            sage: list(a.irange())
            [0, 1, 2, 3]
            sage: list(a.irange(2))
            [2, 3]
        """
        return self.parent().irange(start=start)

    def set_comp(self, i):
        r"""
        Return the `i`-th homogeneous component for assignment.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: A = M.mixed_form(name='A')
            sage: A.set_comp(0).set_expr(x^2) # scalar field
            sage: A.set_comp(1)[:] = [-y, x]
            sage: A.set_comp(2)[0,1] = x-y
            sage: A.display()
            A = A_0 + A_1 + A_2
            sage: A.display_expansion()
            A = x^2 - y dx + x dy + (x - y) dx∧dy
        """
        return self[i]

    def set_immutable(self):
        r"""
        Set ``self`` and homogeneous components of ``self`` immutable.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2, name='f')
            sage: A = M.mixed_form([f, 0, 0], name='A')
            sage: A.set_immutable()
            sage: A.is_immutable()
            True
            sage: A[0].is_immutable()
            True
            sage: f.is_immutable()
            False
        """
        for form in self:
            form.set_immutable()
        super().set_immutable()
