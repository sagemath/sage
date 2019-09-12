r"""
Mixed Differential Forms

Let `M` and `N` be differentiable manifolds and `\varphi : M \longrightarrow N` a
differentiable map. A *mixed differential form along* `\varphi` is an element of
the graded algebra represented by
:class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`.
Its homogeneous components consist of differential forms along `\varphi`. Mixed
forms are useful to represent characteristic classes and perform computations
of such.

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
from sage.structure.element import AlgebraElement
from sage.rings.integer import Integer
from sage.structure.richcmp import richcmp

class MixedForm(AlgebraElement):
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
    dual of the tangent space `T_{\varphi(x)} N`.

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
        sage: F = M.mixed_form(name='F'); F
        Mixed differential form F on the 2-dimensional differentiable manifold M
        sage: F.parent()
        Graded algebra Omega^*(M) of mixed differential forms on the
         2-dimensional differentiable manifold M

    To define the homogenous components of a mixed form, it is convenient to
    define some differential forms first::

        sage: f = M.scalar_field(x, name='f'); f
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: omega = M.diff_form(1, name='omega', latex_name=r'\omega'); omega
        1-form omega on the 2-dimensional differentiable manifold M
        sage: omega[c_xy.frame(),0] = y*x; omega.disp()
        omega = x*y dx
        sage: eta = M.diff_form(2, name='eta', latex_name=r'\eta'); eta
        2-form eta on the 2-dimensional differentiable manifold M
        sage: eta[c_xy.frame(),0,1] = y^2*x; eta.disp()
        eta = x*y^2 dx/\dy

    The components of the mixed form ``F`` can be manipulated very easily::

        sage: F[:] = [f, omega, eta]; F.disp() # display names
        F = f + omega + eta
        sage: F.disp(c_xy.frame()) # display in coordinates
        F = [x] + [x*y dx] + [x*y^2 dx/\dy]
        sage: F[0]
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: F[0] is f
        True
        sage: F[1]
        1-form omega on the 2-dimensional differentiable manifold M
        sage: F[1] is omega
        True
        sage: F[2]
        2-form eta on the 2-dimensional differentiable manifold M
        sage: F[2] is eta
        True

    Alternatively, the components can be determined from scratch::

        sage: G = M.mixed_form(name='G', comp=[f, omega, eta])
        sage: G == F
        True

    Mixed forms are elements of an algebra, so they can be added, and multiplied
    via the wedge product::

        sage: xF = x*F; xF
        Mixed differential form x/\F on the 2-dimensional differentiable
         manifold M
        sage: xF.disp(c_xy.frame())
        x/\F = [x^2] + [x^2*y dx] + [x^2*y^2 dx/\dy]
        sage: FpxF = F+xF; FpxF
        Mixed differential form F+x/\F on the 2-dimensional differentiable
         manifold M
        sage: FpxF.disp(c_xy.frame())
        F+x/\F = [x^2 + x] + [(x^2 + x)*y dx] + [(x^2 + x)*y^2 dx/\dy]
        sage: FxF = F*xF; FxF
        Mixed differential form F/\(x/\F) on the 2-dimensional differentiable
         manifold M
        sage: FxF.disp(c_xy.frame())
        F/\(x/\F) = [x^3] + [2*x^3*y dx] + [2*x^3*y^2 dx/\dy]

    Coercions are fully implemented::

        sage: omegaF = omega*F
        sage: omegaF.disp(c_xy.frame())
        omega/\F = [0] + [x^2*y dx] + [0]
        sage: omegapF = omega+F
        sage: omegapF.disp(c_xy.frame())
        omega+F = [x] + [2*x*y dx] + [x*y^2 dx/\dy]

    Moreover, it is possible to compute the exterior derivative of a
    mixed form::

        sage: dF = F.exterior_derivative(); dF.disp()
        dF = zero + df + domega
        sage: dF.disp(c_xy.frame())
        dF = [0] + [dx] + [-x dx/\dy]

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
        sage: e_xy = c_xy.frame(); e_uv = c_uv.frame() # define frames
        sage: omega = M.diff_form(1, name='omega', latex_name=r'\omega')
        sage: omega[c_xy.frame(),0] = y*x; omega.disp(e_xy)
        omega = x*y dx
        sage: eta = M.diff_form(2, name='eta', latex_name=r'\eta')
        sage: eta[c_uv.frame(),0,1] = u*v^2; eta.disp(e_uv)
        eta = u*v^2 du/\dv
        sage: F = M.mixed_form(name='F', comp=[x, omega, eta]); F
        Mixed differential form F on the 2-dimensional differentiable manifold M
        sage: F.add_comp_by_continuation(e_uv, V.intersection(U), c_uv)
        sage: F.disp(e_uv)
        F = [1/2*u + 1/2*v] + [(1/8*u^2 - 1/8*v^2) du + (1/8*u^2 - 1/8*v^2) dv]
         + [u*v^2 du/\dv]
        sage: F.add_comp_by_continuation(e_xy, V.intersection(U), c_xy)
        sage: F.disp(e_xy)
        F = [x] + [x*y dx] + [(-2*x^3 + 2*x^2*y + 2*x*y^2 - 2*y^3) dx/\dy]

    """
    def __init__(self, parent, comp=None, name=None, latex_name=None):
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
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[c_xy.frame(),0] = y*x
            sage: eta = M.diff_form(2, name='eta')
            sage: eta[c_uv.frame(),0,1] = u*v^2
            sage: A = M.mixed_form_algebra()
            sage: f = A([x, omega, eta], name='F')
            sage: f.add_comp_by_continuation(e_uv, V.intersection(U), c_uv)
            sage: f.add_comp_by_continuation(e_xy, V.intersection(U), c_xy)
            sage: TestSuite(f).run(skip='_test_pickling')

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
        # Set components:
        if comp is None:
            self._comp = [self._domain.diff_form(j)
                          for j in range(self._max_deg + 1)]
        else:
            self._comp = comp
        # Set names:
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

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

        """
        description = "Mixed differential form "
        if self._name is not None:
            description += self._name + " "
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} with values on the {} ".format(
                self._domain, self._ambient_domain)
            if self._dest_map._name is None:
                dm_name = "unnamed map"
            else:
                dm_name = self._dest_map._name
            description += "via the map " + dm_name
        return description

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
            return r'\mbox{' + repr(self) + r'}'
        else:
            return self._latex_name

    def display(self, basis=None, chart=None, from_chart=None):
        r"""
        Display the components of mixed forms.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame with respect to
          which the mixed form is expanded; if ``None``, only the names
          of the components are displayed
        - ``chart`` -- (default: ``None``) chart with respect to which the
          components of the mixed form in the selected frame are expressed;
          if ``None``, the default chart of the vector frame domain is assumed

        EXAMPLES:

        Display a mixed form on a 2-dimensional non-parallelizable
        differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x-y, x+y),
            ....:                   intersection_name='W', restrictions1= x>0,
            ....:                   restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame() # define frames
            sage: omega = M.diff_form(1, name='omega', latex_name=r'\omega')
            sage: omega[c_xy.frame(),0] = x; omega.disp(e_xy)
            omega = x dx
            sage: eta = M.diff_form(2, name='eta', latex_name=r'\eta')
            sage: eta[c_uv.frame(),0,1] = u*v; eta.disp(e_uv)
            eta = u*v du/\dv
            sage: F = M.mixed_form(name='F', comp=[0, omega, eta]); F
            Mixed differential form F on the 2-dimensional differentiable
             manifold M
            sage: F.disp() # display names of homogenous components
            F = zero + omega + eta
            sage: F.add_comp_by_continuation(e_uv, V.intersection(U), c_uv)
            sage: F.disp(e_uv)
            F = [0] + [(1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv] + [u*v du/\dv]
            sage: F.add_comp_by_continuation(e_xy, V.intersection(U), c_xy)
            sage: F.disp(e_xy)
            F = [0] + [x dx] + [(2*x^2 - 2*y^2) dx/\dy]

        """
        from sage.misc.latex import latex, LatexExpr
        from sage.tensor.modules.format_utilities import FormattedExpansion

        def _display_form(rst, basis, chart):
            r"""
            Display coordinate expression of a single form ``rst`` (without
            equality sign).

            This helper method is invoked by :meth:`display`.

            TESTS::

                sage: M = Manifold(2, 'M')
                sage: c_xy.<x,y> = M.chart()
                sage: omega = M.diff_form(1, name='omega')
                sage: omega[c_xy.frame(),0] = x^2
                sage: F = M.mixed_form(comp=[0,omega,0])
                sage: F.disp(c_xy.frame())
                [0] + [x^2 dx] + [0]

            """
            from sage.tensor.modules.format_utilities import is_atomic

            cobasis = basis.dual_basis()
            comp = rst.comp(basis)
            terms_txt = []
            terms_latex = []
            for ind in comp.non_redundant_index_generator():
                ind_arg = ind + (chart,)
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
                    basis_term_txt = "/\\".join(bases_txt)
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
                            terms_txt.append(
                                "(" + coef_txt + ") " + basis_term_txt)
                        if is_atomic(coef_latex):
                            terms_latex.append(coef_latex + basis_term_latex)
                        else:
                            terms_latex.append(r"\left(" + coef_latex +
                                               r"\right)" + basis_term_latex)
            if not terms_txt:
                resu_txt = "0"
            else:
                resu_txt = terms_txt[0]
                for term in terms_txt[1:]:
                    if term[0] == "-":
                        resu_txt += " - " + term[1:]
                    else:
                        resu_txt += " + " + term
            if not terms_latex:
                resu_latex = r"0"
            else:
                resu_latex = terms_latex[0]
                for term in terms_latex[1:]:
                    if term[0] == "-":
                        resu_latex += term
                    else:
                        resu_latex += "+" + term
            return FormattedExpansion(resu_txt, resu_latex)
        plain_txt = ""
        latex_txt = r""
        if self._name is not None:
            plain_txt += self._name + " = "
        if self._latex_name is not None:
            latex_txt += latex(self) + LatexExpr(r"=")
        if basis is None:
            if self[0]._name is None:
                plain_txt += "(unnamed scalar field) "
            else:
                plain_txt += self[0]._name
            if self[0]._latex_name is None:
                latex_txt += LatexExpr(r"\mbox{(unnamed scalar field)}")
            else:
                latex_txt += latex(self[0])
            for j in range(1, self._max_deg + 1):
                if self[j]._name is None:
                    plain_txt += " + (unnamed " + str(j) + "-form)"
                else:
                    plain_txt += " + " + self[j]._name
                if self[j]._latex_name is None:
                    latex_txt += LatexExpr(
                        r"+\mbox{(unnamed " + str(j) + r"-form)}")
                else:
                    latex_txt += LatexExpr(r"+") + latex(self[j])
            return FormattedExpansion(plain_txt, latex_txt)
        # In case no other chart is given:
        if chart is None:
            chart = basis._chart
        plain_txt += "[" + repr(self[0].expr(chart, from_chart)) + "]"
        latex_txt += LatexExpr(r"\left[") + latex(
            self[0].expr(chart, from_chart)) \
                     + LatexExpr(r"\right]_0")
        for j in range(1, self._max_deg + 1):
            rst = self[j].restrict(basis._domain,
                                         dest_map=basis._dest_map)
            format_txt = _display_form(rst, basis, chart)
            plain_txt += " + [" + repr(format_txt) + "]"
            latex_txt += (LatexExpr(r"+ \left[") + latex(format_txt)
                         + LatexExpr(r"\right]_{}".format(j)))
        return FormattedExpansion(plain_txt, latex_txt)

    disp = display

    def set_name(self, name=None, latex_name=None):
        r"""
        Redefine the string and LaTeX representation of the object.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the mixed form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          mixed form; if none is provided, the LaTeX symbol is set to ``name``

        EXAMPLES::

            sage: M = Manifold(4, 'M')
            sage: F = M.mixed_form(name='dummy', latex_name=r'\ugly'); F
            Mixed differential form dummy on the 4-dimensional differentiable
             manifold M
            sage: latex(F)
            \ugly
            sage: F.set_name(name='fancy', latex_name=r'\eta'); F
            Mixed differential form fancy on the 4-dimensional differentiable
             manifold M
            sage: latex(F)
            \eta

        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name

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
        return any(bool(form) for form in self._comp)

    __nonzero__ = __bool__  # For Python2 compatibility

    def _richcmp_(self, other, op):
        r"""
        Comparison method for sage objects.

        INPUT:

        - ``other`` -- a mixed form
        - ``op`` -- comparison operator for which ``self`` and ``other`` shall
          be compared with.

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
            sage: f = M.scalar_field(x, name='f')
            sage: eta = M.diff_form(1, name='eta')
            sage: eta[e_xy,0] = x+y
            sage: F = M.mixed_form(comp=[f, eta, 0])
            sage: F.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: F == F
            True
            sage: F == F.copy()
            True
            sage: G = M.mixed_form()
            sage: G.set_restriction(F.restrict(U))
            sage: F == G  # False since G has not been defined on V
            False
            sage: G.set_restriction(F.restrict(V))
            sage: F == G  # True now
            True
            sage: F.parent().zero() == 0
            True

        """
        # Compare all elements separately:
        for j in range(self._max_deg + 1):
            if not richcmp(self[j], other[j], op):
                return False
        return True

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
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: f = M.scalar_field(x, name='f')
            sage: a = M.diff_form(1, name='a')
            sage: a[e_xy,0] = x
            sage: A = M.mixed_form(name='A', comp=[f, a, 0])
            sage: A.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: g = M.scalar_field(u, name='g', chart=c_uv)
            sage: b = M.diff_form(1, name='b')
            sage: b[e_uv,1] = v
            sage: B = M.mixed_form(name='B', comp=[g, b, 0])
            sage: B.add_comp_by_continuation(e_xy, V.intersection(U), c_xy)
            sage: C = A._add_(B); C
            Mixed differential form A+B on the 2-dimensional differentiable
             manifold M
            sage: A.display(e_uv)
            A = [1/2*u + 1/2*v] + [(1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv]
             + [0]
            sage: B.display(e_xy)
            B = [x + y] + [(x - y) dx + (-x + y) dy] + [0]
            sage: C.display(e_xy)
            A+B = [2*x + y] + [(2*x - y) dx + (-x + y) dy] + [0]
            sage: C.display(e_uv)
            A+B = [3/2*u + 1/2*v] + [(1/4*u + 1/4*v) du + (1/4*u + 5/4*v) dv]
             + [0]
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
        resu_comp = [self[j] + other[j]
                     for j in range(self._max_deg + 1)]
        resu = type(self)(self.parent(), comp=resu_comp)
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
        sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
        sage: f = M.scalar_field(x, name='f')
        sage: a = M.diff_form(1, name='a')
        sage: a[e_xy,0] = x
        sage: A = M.mixed_form(name='A', comp=[f, a, 0])
        sage: A.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
        sage: g = M.scalar_field(u, name='g', chart=c_uv)
        sage: b = M.diff_form(1, name='b')
        sage: b[e_uv,1] = v
        sage: B = M.mixed_form(name='B', comp=[g, b, 0])
        sage: B.add_comp_by_continuation(e_xy, V.intersection(U), c_xy)
        sage: C = A._sub_(B); C
        Mixed differential form A-B on the 2-dimensional differentiable
         manifold M
        sage: A.display(e_uv)
        A = [1/2*u + 1/2*v] + [(1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv] + [0]
        sage: B.display(e_xy)
        B = [x + y] + [(x - y) dx + (-x + y) dy] + [0]
        sage: C.display(e_xy)
        A-B = [-y] + [y dx + (x - y) dy] + [0]
        sage: C.display(e_uv)
        A-B = [-1/2*u + 1/2*v] + [(1/4*u + 1/4*v) du + (1/4*u - 3/4*v) dv] + [0]
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
        resu_comp = [self[j] - other[j]
                     for j in range(self._max_deg + 1)]
        resu = type(self)(self.parent(), comp=resu_comp)
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
            sage: f.disp()
            f: M --> R
               (x, y, z) |--> x
            sage: g = M.scalar_field(y, name='g')
            sage: g.disp()
            g: M --> R
               (x, y, z) |--> y
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[c_xyz.frame(),0] = x
            sage: omega.disp()
            omega = x dx
            sage: eta = M.diff_form(1, name='eta')
            sage: eta[c_xyz.frame(),1] = y
            sage: eta.disp()
            eta = y dy
            sage: mu = M.diff_form(2, name='mu')
            sage: mu[c_xyz.frame(),[0,2]] = z
            sage: mu.disp()
            mu = z dx/\dz
            sage: A = M.mixed_form(name='A', comp=[f, omega, mu, 0])
            sage: A.disp(c_xyz.frame())
            A = [x] + [x dx] + [z dx/\dz] + [0]
            sage: B = M.mixed_form(name='B', comp=[g, eta, mu, 0])
            sage: B.disp(c_xyz.frame())
            B = [y] + [y dy] + [z dx/\dz] + [0]

        The wedge product of ``A`` and ``B`` yields::

            sage: AwB = A.wedge(B); AwB
            Mixed differential form A/\B on the 3-dimensional differentiable
             manifold M
            sage: AwB.disp(c_xyz.frame())
            A/\B = [x*y] + [x*y dx + x*y dy] + [x*y dx/\dy + (x + y)*z dx/\dz] +
             [-y*z dx/\dy/\dz]
            sage: BwA = B.wedge(A); BwA # Don't even try, it's not commutative!
            Mixed differential form B/\A on the 3-dimensional differentiable
             manifold M
            sage: BwA.disp(c_xyz.frame()) # I told you so!
            B/\A = [x*y] + [x*y dx + x*y dy] + [-x*y dx/\dy + (x + y)*z dx/\dz]
             + [-y*z dx/\dy/\dz]

        The multiplication symbol may be used instead as well::

            sage: A*B
            Mixed differential form A/\B on the 3-dimensional differentiable
             manifold M
            sage: A*B == AwB
            True

        Yet, the multiplication includes coercions::

            sage: xA = x*A; xA.disp(c_xyz.frame())
            x/\A = [x^2] + [x^2 dx] + [x*z dx/\dz] + [0]
            sage: Aeta = A*eta; Aeta.disp(c_xyz.frame())
            A/\eta = [0] + [x*y dy] + [x*y dx/\dy] + [-y*z dx/\dy/\dz]

        """
        resu_comp = [None] * (self._max_deg + 1)
        resu_comp[0] = self[0] * other[0]
        for j in range(1, self._max_deg + 1):
            resu_comp[j] = self[j] * other[0]
            resu_comp[j] = resu_comp[j] + self[0] * other[j]
            for k in range(1, j):
                resu_comp[j] = (resu_comp[j] +
                               self[k].wedge(other[j - k]))
        resu = type(self)(self.parent(), comp=resu_comp)
        # Compose name:
        from sage.tensor.modules.format_utilities import (format_mul_txt,
            format_mul_latex)
        resu._name = format_mul_txt(self._name, '/\\', other._name)
        resu._latex_name = format_mul_latex(self._latex_name, r'\wedge ',
                                            other._latex_name)
        return resu

    _mul_ = wedge

    def _lmul_(self, other):
        r"""
        Scalar multiplication operator: return ``number * self`` or
        ``self * number``.

        INPUT:

        - ``num`` -- an element of the symbolic ring

        OUTPUT:

        - the mixed form resulting from the wedge product of ``self``
          with ``other``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[c_xy.frame(),0] = y*x; omega.disp()
            omega = x*y dx
            sage: F = M.mixed_form(name='F', comp=[0, omega, 0])
            sage: xFy = x*F*y; xFy
            Mixed differential form y/\(x/\F) on the 2-dimensional
             differentiable manifold M
            sage: xFy.disp(c_xy.frame())
            y/\(x/\F) = [0] + [x^2*y^2 dx] + [0]

        """
        resu_comp = [other * form for form in self._comp]
        resu = type(self)(self.parent(), comp=resu_comp)
        # Compose name:
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import (format_mul_txt,
            format_mul_latex)
        resu._name = format_mul_txt(repr(other), '/\\', self._name)
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
            f: M --> R
                (x, y, z) |--> z^2
            sage: a = M.diff_form(2, 'a')
            sage: a[1,2], a[1,3], a[2,3] = z+y^2, z+x, x^2
            sage: a.disp()
            a = (y^2 + z) dx/\dy + (x + z) dx/\dz + x^2 dy/\dz
            sage: F = M.mixed_form(name='F', comp=[f, 0, a, 0]); F.disp()
            F = f + zero + a + zero
            sage: dF = F.exterior_derivative()
            sage: dF.disp()
            dF = zero + df + dzero + da
            sage: dF = F.exterior_derivative()
            sage: dF.disp(c_xyz.frame())
            dF = [0] + [2*z dz] + [0] + [(2*x + 1) dx/\dy/\dz]

        Due to long calculation times, the result is cached, i.e. is not
        recomputed unless ``F`` is changed::

            sage: F.exterior_derivative() is dF
            True

        """
        resu_comp = list()
        resu_comp.append(self._domain.scalar_field_algebra().zero())
        resu_comp.append(self[0].differential())
        resu_comp.extend([self[j].exterior_derivative()
                          for j in range(1, self._max_deg)])
        resu = type(self)(self.parent(), comp=resu_comp)
        # Compose name:
        from sage.tensor.modules.format_utilities import (format_unop_txt,
                                                          format_unop_latex)
        resu._name = format_unop_txt('d', self._name)
        resu._latex_name = format_unop_latex(r'\mathrm{d}', self._latex_name)

        return resu

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        .. NOTE::

            The name and names of the components are not copied.

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
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: f = M.scalar_field(x, name='f', chart=c_xy); f.disp()
            f: M --> R
                on U: (x, y) |--> x
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[e_xy,0] = x; omega.disp()
            omega = x dx
            sage: A = M.mixed_form(name='A', comp=[f, omega, 0]); A.disp()
            A = f + omega + zero
            sage: A.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: A.disp(e_uv)
            A = [1/2*u + 1/2*v] + [(1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv] +
             [0]

        An exact copy is made. The copy is an entirely new instance and has a
        different name, but has the very same values::

            sage: B = A.copy(); B.disp()
            f + (unnamed 1-form) + (unnamed 2-form)
            sage: B.disp(e_uv)
            [1/2*u + 1/2*v] + [(1/4*u + 1/4*v) du + (1/4*u + 1/4*v) dv] + [0]
            sage: A == B
            True
            sage: A is B
            False

        Notice, that changes in the differential forms usually cause changes in
        the original instance. But for the copy of a mixed form, the components
        are copied as well::

            sage: omega[e_xy,0] = y; omega.disp()
            omega = y dx
            sage: A.disp(e_xy)
            A = [x] + [y dx] + [0]
            sage: B.disp(e_xy)
            [x] + [x dx] + [0]

        """
        resu_comp = [form.copy() for form in self._comp]
        return type(self)(self.parent(), comp=resu_comp)

    def __setitem__(self, index, values):
        r"""
        Sets a component with respect to some vector frame.

        - ``index`` -- list of indices; if ``[:]`` is provided, all the
          components are set
        - ``values`` -- the values to be set

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x, name='f')
            sage: a = M.diff_form(1, name='a')
            sage: a[c_xy.frame(),0] = y
            sage: b = M.diff_form(2, name='b')
            sage: b[c_xy.frame(),0,1] = x*y
            sage: A = M.mixed_form(name='A', comp=[f, 0, 0]); A.disp()
            A = f + zero + zero
            sage: A[1:3] = [a, b]; A.disp()
            A = f + a + b
            sage: A.disp(c_xy.frame())
            A = [x] + [y dx] + [x*y dx/\dy]

        """
        if isinstance(index, (int, Integer)):
            start = index
            stop = index + 1
            step = 1
        elif isinstance(index, slice):
            start, stop, step = index.indices(len(self._comp))
        else:
            raise TypeError("index must be int, Integer or slice")
        if isinstance(values, list):
            form_list = values
        else:
            form_list = [values]
        if len(form_list) != len(range(start, stop, step)):
            raise IndexError("either input or index out of range")
        for deg, j in zip(range(start, stop, step), range(len(form_list))):
            self._comp[deg] = self._domain.diff_form_module(deg,
                                            self._dest_map)(form_list[j])

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
            sage: A = M.mixed_form(name='A', comp=[f, a, b]); A.disp()
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
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()

        And define some forms on the subset ``U``::

            sage: f = U.scalar_field(x, name='f', chart=c_xy)
            sage: omega = U.diff_form(1, name='omega')
            sage: omega[e_xy,0] = y
            sage: AU = U.mixed_form(name='A', comp=[f, omega, 0]); AU
            Mixed differential form A on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: AU.disp(e_xy)
            A = [x] + [y dx] + [0]

        A mixed form on ``M`` can be specified by some mixed form on a subset::

            sage: A = M.mixed_form(name='A'); A
            Mixed differential form A on the 2-dimensional differentiable
             manifold M
            sage: A.set_restriction(AU)
            sage: A.disp(e_xy)
            A = [x] + [y dx] + [0]
            sage: A.add_comp_by_continuation(e_uv, V.intersection(U), c_uv)
            sage: A.disp(e_uv)
            A = [u/(u^2 + v^2)] + [-(u^2*v - v^3)/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 +
             v^6) du - 2*u*v^2/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) dv] + [0]
            sage: A.restrict(U) == AU
            True

        """
        if not isinstance(rst, MixedForm):
            raise TypeError("the argument must be a mixed form")
        subdomain = rst._domain
        if not subdomain.is_subset(self._domain):
            raise ValueError("the specified domain is not a subset of " +
                             "the domain of definition of the mixed form")
        # set restrictions for scalar field:
        self[0]._restrictions[subdomain] = rst[0]
        self[0]._is_zero = False  # a priori
        for chart, expr in rst[0]._express.items():
            self[0]._express[chart] = expr # automatic continuation to chart dom
        for j in range(1, self._max_deg + 1):
            self[j].set_restriction(rst[j])

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

        OUTPUT:

        - :class:`MixedForm` representing the restriction

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
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()

        And predefine some forms::

            sage: f = M.scalar_field(x^2, name='f', chart=c_xy)
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[e_xy,0] = y^2
            sage: eta = M.diff_form(2, name='eta')
            sage: eta[e_xy,0,1] = x^2*y^2

        Now, a mixed form can be restricted to some subdomain::

            sage: F = M.mixed_form(name='F', comp=[f, omega, eta])
            sage: F.add_comp_by_continuation(e_uv, V.intersection(U), c_uv)
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
            sage: FV.disp(e_uv)
            F = [u^2/(u^4 + 2*u^2*v^2 + v^4)] + [-(u^2*v^2 - v^4)/(u^8 +
             4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du - 2*u*v^3/(u^8 +
             4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv] + [-u^2*v^2/(u^12 +
             6*u^10*v^2 + 15*u^8*v^4 + 20*u^6*v^6 + 15*u^4*v^8 + 6*u^2*v^10 +
             v^12) du/\dv]

        """
        resu_comp = list()
        resu_comp.append(self[0].restrict(subdomain))
        resu_comp.extend([self[j].restrict(subdomain, dest_map)
                          for j in range(1, self._max_deg + 1)])
        return type(self)(subdomain.mixed_form_algebra(dest_map=dest_map),
                              resu_comp, self._name, self._latex_name)

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
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: f = M.scalar_field(x, name='f')
            sage: omega = M.diff_form(1, name='omega')
            sage: omega[e_xy,0] = x
            sage: eta = M.diff_form(2, name='eta')
            sage: eta[e_uv,0,1] = u*v
            sage: F = M.mixed_form(name='F', comp=[f, omega, eta])
            sage: F.add_comp_by_continuation(e_xy, U.intersection(V), c_xy)
            sage: F.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: F.disp(e_xy)
            F = [x] + [x dx] + [-x*y/(x^8 + 4*x^6*y^2 + 6*x^4*y^4 + 4*x^2*y^6 +
             y^8) dx/\dy]
            sage: F.disp(e_uv)
            F = [u/(u^2 + v^2)] + [-(u^3 - u*v^2)/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 +
             v^6) du - 2*u^2*v/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) dv] +
             [u*v du/\dv]

        """
        if chart is None:
            chart = frame._chart
        self[0].add_expr_by_continuation(chart, subdomain)
        for j in range(1, self._max_deg + 1):
            self[j].add_comp_by_continuation(frame, subdomain, chart)
