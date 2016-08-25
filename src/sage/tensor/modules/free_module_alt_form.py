r"""
Alternating forms on free modules

Given a free module `M` of finite rank over a commutative ring `R`
and a positive integer `p`, an *alternating form of degree* `p` on `M` is
a map

.. MATH::

    a:\ \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
    \longrightarrow R

that (i) is multilinear and (ii) vanishes whenever any of two of its
arguments are equal.
An alternating form of degree `p` is a tensor on `M` of type `(0,p)`.

Alternating forms are implemented via the class :class:`FreeModuleAltForm`,
which is a subclass of the generic tensor class
:class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

REFERENCES:

- Chap. 23 of R. Godement: *Algebra*, Hermann (Paris) / Houghton Mifflin
  (Boston) (1968)
- Chap. 15 of S. Lang: *Algebra*, 3rd ed., Springer (New York) (2002)

"""
from __future__ import absolute_import
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.comp import Components, CompFullyAntiSym

class FreeModuleAltForm(FreeModuleTensor):
    r"""
    Alternating form on a free module of finite rank over a commutative ring.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank over a commutative ring
      `R`, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``degree`` -- positive integer; the degree `p` of the
      alternating form (i.e. its tensor rank)
    - ``name`` -- (default: ``None``) string; name given to the alternating
      form
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      alternating form; if none is provided, ``name`` is used

    EXAMPLES:

    Alternating form of degree 2 on a rank-3 free module::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: e = M.basis('e')
        sage: a = M.alternating_form(2, name='a') ; a
        Alternating form a of degree 2 on the
         Rank-3 free module M over the Integer Ring
        sage: type(a)
        <class 'sage.tensor.modules.free_module_alt_form.ExtPowerFreeModule_with_category.element_class'>
        sage: a.parent()
        2nd exterior power of the dual of the Rank-3 free module M over the Integer Ring
        sage: a[1,2], a[2,3] = 4, -3
        sage: a.display(e)
        a = 4 e^1/\e^2 - 3 e^2/\e^3

    The alternating form acting on the basis elements::

        sage: a(e[1],e[2])
        4
        sage: a(e[1],e[3])
        0
        sage: a(e[2],e[3])
        -3
        sage: a(e[2],e[1])
        -4

    An alternating form of degree 1 is a linear form::

        sage: b = M.linear_form('b') ; b
        Linear form b on the Rank-3 free module M over the Integer Ring
        sage: b[:] = [2,-1,3]  # components w.r.t. the module's default basis (e)

    A linear form is a tensor of type (0,1)::

        sage: b.tensor_type()
        (0, 1)

    It is an element of the dual module::

        sage: b.parent()
        Dual of the Rank-3 free module M over the Integer Ring
        sage: b.parent() is M.dual()
        True

    The members of a dual basis are linear forms::

        sage: e.dual_basis()[1]
        Linear form e^1 on the Rank-3 free module M over the Integer Ring
        sage: e.dual_basis()[2]
        Linear form e^2 on the Rank-3 free module M over the Integer Ring
        sage: e.dual_basis()[3]
        Linear form e^3 on the Rank-3 free module M over the Integer Ring

    Any linear form is expanded onto them::

        sage: b.display(e)
        b = 2 e^1 - e^2 + 3 e^3

    In the above example, an equivalent writing would have been
    ``b.display()``, since the basis ``e`` is the module's default basis.
    A linear form maps module elements to ring elements::

        sage: v = M([1,1,1])
        sage: b(v)
        4
        sage: b(v) in M.base_ring()
        True

    Test of linearity::

        sage: u = M([-5,-2,7])
        sage: b(3*u - 4*v) == 3*b(u) - 4*b(v)
        True

    The standard tensor operations apply to alternating forms, like the
    extraction of components with respect to a given basis::

        sage: a[e,1,2]
        4
        sage: a[1,2]  # since e is the module's default basis
        4
        sage: all( a[i,j] == - a[j,i] for i in {1,2,3} for j in {1,2,3} )
        True

    the tensor product::

        sage: c = b*b ; c
        Symmetric bilinear form  b*b on the Rank-3 free module M over the
         Integer Ring
        sage: c.parent()
        Free module of type-(0,2) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: c.display(e)
        b*b = 4 e^1*e^1 - 2 e^1*e^2 + 6 e^1*e^3 - 2 e^2*e^1 + e^2*e^2
         - 3 e^2*e^3 + 6 e^3*e^1 - 3 e^3*e^2 + 9 e^3*e^3

    the contractions::

        sage: s = a.contract(v) ; s
        Linear form on the Rank-3 free module M over the Integer Ring
        sage: s.parent()
        Dual of the Rank-3 free module M over the Integer Ring
        sage: s.display(e)
        4 e^1 - 7 e^2 + 3 e^3

    or tensor arithmetics::

        sage: s = 3*a + c ; s
        Type-(0,2) tensor on the Rank-3 free module M over the Integer Ring
        sage: s.parent()
        Free module of type-(0,2) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: s.display(e)
        4 e^1*e^1 + 10 e^1*e^2 + 6 e^1*e^3 - 14 e^2*e^1 + e^2*e^2
         - 12 e^2*e^3 + 6 e^3*e^1 + 6 e^3*e^2 + 9 e^3*e^3

    Note that tensor arithmetics preserves the alternating character if both
    operands are alternating::

        sage: s = a - 2*a ; s
        Alternating form of degree 2 on the Rank-3 free module M over the
         Integer Ring
        sage: s.parent() # note the difference with s = 3*a + c above
        2nd exterior power of the dual of the Rank-3 free module M over the
         Integer Ring
        sage: s == -a
        True

    An operation specific to alternating forms is of course the exterior
    product::

        sage: s = a.wedge(b) ; s
        Alternating form a/\b of degree 3 on the Rank-3 free module M over the
         Integer Ring
        sage: s.parent()
        3rd exterior power of the dual of the Rank-3 free module M over the
         Integer Ring
        sage: s.display(e)
        a/\b = 6 e^1/\e^2/\e^3
        sage: s[1,2,3] == a[1,2]*b[3] + a[2,3]*b[1] + a[3,1]*b[2]
        True

    The exterior product is nilpotent on linear forms::

        sage: s = b.wedge(b) ; s
        Alternating form b/\b of degree 2 on the Rank-3 free module M over the
         Integer Ring
        sage: s.display(e)
        b/\b = 0

    """
    def __init__(self, fmodule, degree, name=None, latex_name=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = FreeModuleAltForm(M, 2, name='a')
            sage: a[e,0,1] = 2
            sage: TestSuite(a).run(skip="_test_category") # see below

        In the above test suite, _test_category fails because a is not an
        instance of a.parent().category().element_class. Actually alternating
        forms must be constructed via ExtPowerFreeModule.element_class and
        not by a direct call to FreeModuleAltForm::

            sage: a1 = M.dual_exterior_power(2).element_class(M, 2, name='a')
            sage: a1[e,0,1] = 2
            sage: TestSuite(a1).run()

        """
        FreeModuleTensor.__init__(self, fmodule, (0,degree), name=name,
                                  latex_name=latex_name, antisym=range(degree),
                                  parent=fmodule.dual_exterior_power(degree))
        FreeModuleAltForm._init_derived(self) # initialization of derived
                                              # quantities

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.alternating_form(1)
            Linear form on the Rank-3 free module M over the Integer Ring
            sage: M.alternating_form(1, name='a')
            Linear form a on the Rank-3 free module M over the Integer Ring
            sage: M.alternating_form(2)
            Alternating form of degree 2 on the
             Rank-3 free module M over the Integer Ring
            sage: M.alternating_form(2, name='a')
            Alternating form a of degree 2 on the
             Rank-3 free module M over the Integer Ring
        """
        if self._tensor_rank == 1:
            description = "Linear form "
            if self._name is not None:
                description += self._name + " "
        else:
            description = "Alternating form "
            if self._name is not None:
                description += self._name + " "
            description += "of degree {} ".format(self._tensor_rank)
        description += "on the {}".format(self._fmodule)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_form(2)
            sage: a._init_derived()

        """
        FreeModuleTensor._init_derived(self)

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_form(2)
            sage: a._del_derived()

        """
        FreeModuleTensor._del_derived(self)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``, on the same module
        and of the same degree.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_form(1)
            sage: a._new_instance()
            Linear form on the Rank-3 free module M over the Integer Ring
            sage: a._new_instance().parent() is a.parent()
            True
            sage: b = M.alternating_form(2, name='b')
            sage: b._new_instance()
            Alternating form of degree 2 on the
             Rank-3 free module M over the Integer Ring
            sage: b._new_instance().parent() is b.parent()
            True

        """
        return self.__class__(self._fmodule, self._tensor_rank)

    def _new_comp(self, basis):
        r"""
        Create some (uninitialized) components of ``self`` in a given basis.

        This method, which is already implemented in
        :meth:`FreeModuleTensor._new_comp`, is redefined here for efficiency.

        INPUT:

        - ``basis`` -- basis of the free module on which ``self`` is defined

        OUTPUT:

        - an instance of :class:`~sage.tensor.modules.comp.CompFullyAntiSym`,
          or of :class:`~sage.tensor.modules.comp.Components` if
          the degree of ``self`` is 1.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.alternating_form(2, name='a')
            sage: a._new_comp(e)
            Fully antisymmetric 2-indices components w.r.t. Basis (e_0,e_1,e_2)
             on the Rank-3 free module M over the Integer Ring
            sage: a = M.alternating_form(1)
            sage: a._new_comp(e)
            1-index components w.r.t. Basis (e_0,e_1,e_2) on the Rank-3 free
             module M over the Integer Ring

        """
        fmodule = self._fmodule  # the base free module
        if self._tensor_rank == 1:
            return Components(fmodule._ring, basis, 1,
                              start_index=fmodule._sindex,
                              output_formatter=fmodule._output_formatter)

        return CompFullyAntiSym(fmodule._ring, basis, self._tensor_rank,
                                start_index=fmodule._sindex,
                                output_formatter=fmodule._output_formatter)

    def degree(self):
        r"""
        Return the degree of ``self``.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_form(2, name='a')
            sage: a.degree()
            2

        """
        return self._tensor_rank


    def display(self, basis=None, format_spec=None):
        r"""
        Display the alternating form ``self`` in terms of its expansion w.r.t.
        a given module basis.

        The expansion is actually performed onto exterior products of elements
        of the cobasis (dual basis) associated with ``basis`` (see examples
        below). The output is either text-formatted (console mode) or
        LaTeX-formatted (notebook mode).

        INPUT:

        - ``basis`` -- (default: ``None``) basis of the free module with
          respect to which the alternating form is expanded; if none is
          provided, the module's default basis is assumed
        - ``format_spec`` -- (default: ``None``) format specification passed
          to ``self._fmodule._output_formatter`` to format the output

        EXAMPLES:

        Display of an alternating form of degree 1 (linear form) on a rank-3
        free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.dual_basis()
            Dual basis (e^0,e^1,e^2) on the Rank-3 free module M over the Integer Ring
            sage: a = M.linear_form('a', latex_name=r'\alpha')
            sage: a[:] = [1,-3,4]
            sage: a.display(e)
            a = e^0 - 3 e^1 + 4 e^2
            sage: a.display()  # a shortcut since e is M's default basis
            a = e^0 - 3 e^1 + 4 e^2
            sage: latex(a.display())  # display in the notebook
            \alpha = e^0 -3 e^1 + 4 e^2

        A shortcut is ``disp()``::

            sage: a.disp()
            a = e^0 - 3 e^1 + 4 e^2

        Display of an alternating form of degree 2 on a rank-3 free module::

            sage: b = M.alternating_form(2, 'b', latex_name=r'\beta')
            sage: b[0,1], b[0,2], b[1,2] = 3, 2, -1
            sage: b.display()
            b = 3 e^0/\e^1 + 2 e^0/\e^2 - e^1/\e^2
            sage: latex(b.display())  # display in the notebook
            \beta = 3 e^0\wedge e^1 + 2 e^0\wedge e^2 -e^1\wedge e^2

        Display of an alternating form of degree 3 on a rank-3 free module::

            sage: c = M.alternating_form(3, 'c')
            sage: c[0,1,2] = 4
            sage: c.display()
            c = 4 e^0/\e^1/\e^2
            sage: latex(c.display())
            c = 4 e^0\wedge e^1\wedge e^2

        Display of a vanishing alternating form::

            sage: c[0,1,2] = 0  # the only independent component set to zero
            sage: c.is_zero()
            True
            sage: c.display()
            c = 0
            sage: latex(c.display())
            c = 0
            sage: c[0,1,2] = 4  # value restored for what follows

        Display in a basis which is not the default one::

            sage: aut = M.automorphism(matrix=[[0,1,0], [0,0,-1], [1,0,0]],
            ....:                      basis=e)
            sage: f = e.new_basis(aut, 'f')
            sage: a.display(f)
            a = 4 f^0 + f^1 + 3 f^2
            sage: a.disp(f)     # shortcut notation
            a = 4 f^0 + f^1 + 3 f^2
            sage: b.display(f)
            b = -2 f^0/\f^1 - f^0/\f^2 - 3 f^1/\f^2
            sage: c.display(f)
            c = -4 f^0/\f^1/\f^2

        The output format can be set via the argument ``output_formatter``
        passed at the module construction::

            sage: N = FiniteRankFreeModule(QQ, 3, name='N', start_index=1,
            ....:                   output_formatter=Rational.numerical_approx)
            sage: e = N.basis('e')
            sage: b = N.alternating_form(2, 'b')
            sage: b[1,2], b[1,3], b[2,3] = 1/3, 5/2, 4
            sage: b.display()  # default format (53 bits of precision)
            b = 0.333333333333333 e^1/\e^2 + 2.50000000000000 e^1/\e^3
             + 4.00000000000000 e^2/\e^3

        The output format is then controled by the argument ``format_spec`` of
        the method :meth:`display`::

            sage: b.display(format_spec=10)  # 10 bits of precision
            b = 0.33 e^1/\e^2 + 2.5 e^1/\e^3 + 4.0 e^2/\e^3

        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import is_atomic, \
                                                         FormattedExpansion
        if basis is None:
            basis = self._fmodule._def_basis
        cobasis = basis.dual_basis()
        comp = self.comp(basis)
        terms_txt = []
        terms_latex = []
        for ind in comp.non_redundant_index_generator():
            ind_arg = ind + (format_spec,)
            coef = comp[ind_arg]
            if coef != 0:
                bases_txt = []
                bases_latex = []
                for k in range(self._tensor_rank):
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
                        terms_txt.append("(" + coef_txt + ") " +
                                         basis_term_txt)
                    if is_atomic(coef_latex):
                        terms_latex.append(coef_latex + basis_term_latex)
                    else:
                        terms_latex.append(r"\left(" + coef_latex + \
                                           r"\right)" + basis_term_latex)
        if not terms_txt:
            expansion_txt = "0"
        else:
            expansion_txt = terms_txt[0]
            for term in terms_txt[1:]:
                if term[0] == "-":
                    expansion_txt += " - " + term[1:]
                else:
                    expansion_txt += " + " + term
        if not terms_latex:
            expansion_latex = "0"
        else:
            expansion_latex = terms_latex[0]
            for term in terms_latex[1:]:
                if term[0] == "-":
                    expansion_latex += term
                else:
                    expansion_latex += "+" + term
        if self._name is None:
            resu_txt = expansion_txt
        else:
            resu_txt = self._name + " = " + expansion_txt
        if self._latex_name is None:
            resu_latex = expansion_latex
        else:
            resu_latex = latex(self) + " = " + expansion_latex
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display


    def wedge(self, other):
        r"""
        Exterior product of ``self`` with the alternating form ``other``.

        INPUT:

        - ``other`` -- an alternating form

        OUTPUT:

        - instance of :class:`FreeModuleAltForm` representing the exterior
          product ``self/\other``

        EXAMPLES:

        Exterior product of two linear forms::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.linear_form('A')
            sage: a[:] = [1,-3,4]
            sage: b = M.linear_form('B')
            sage: b[:] = [2,-1,2]
            sage: c = a.wedge(b) ; c
            Alternating form A/\B of degree 2 on the Rank-3 free module M
             over the Integer Ring
            sage: c.display()
            A/\B = 5 e^0/\e^1 - 6 e^0/\e^2 - 2 e^1/\e^2
            sage: latex(c)
            A\wedge B
            sage: latex(c.display())
            A\wedge B = 5 e^0\wedge e^1 -6 e^0\wedge e^2 -2 e^1\wedge e^2

        Test of the computation::

            sage: a.wedge(b) == a*b - b*a
            True

        Exterior product of a linear form and an alternating form of degree 2::

            sage: d = M.linear_form('D')
            sage: d[:] = [-1,2,4]
            sage: s = d.wedge(c) ; s
            Alternating form D/\A/\B of degree 3 on the Rank-3 free module M
             over the Integer Ring
            sage: s.display()
            D/\A/\B = 34 e^0/\e^1/\e^2

        Test of the computation::

            sage: s[0,1,2] == d[0]*c[1,2] + d[1]*c[2,0] + d[2]*c[0,1]
            True

        Let us check that the exterior product is associative::

            sage: d.wedge(a.wedge(b)) == (d.wedge(a)).wedge(b)
            True

        and that it is graded anticommutative::

            sage: a.wedge(b) == - b.wedge(a)
            True
            sage: d.wedge(c) == c.wedge(d)
            True

        """
        from .format_utilities import is_atomic
        if not isinstance(other, FreeModuleAltForm):
            raise TypeError("the second argument for the exterior product " +
                            "must be an alternating form")
        if other._tensor_rank == 0:
            return other*self
        if self._tensor_rank == 0:
            return self*other
        fmodule = self._fmodule
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("no common basis for the exterior product")
        rank_r = self._tensor_rank + other._tensor_rank
        cmp_s = self._components[basis]
        cmp_o = other._components[basis]
        cmp_r = CompFullyAntiSym(fmodule._ring, basis, rank_r,
                                 start_index=fmodule._sindex,
                                 output_formatter=fmodule._output_formatter)
        for ind_s, val_s in cmp_s._comp.iteritems():
            for ind_o, val_o in cmp_o._comp.iteritems():
                ind_r = ind_s + ind_o
                if len(ind_r) == len(set(ind_r)): # all indices are different
                    cmp_r[[ind_r]] += val_s * val_o
        result = fmodule.alternating_form(rank_r)
        result._components[basis] = cmp_r
        if self._name is not None and other._name is not None:
            sname = self._name
            oname = other._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            result._name = sname + '/\\' + oname
        if self._latex_name is not None and other._latex_name is not None:
            slname = self._latex_name
            olname = other._latex_name
            if not is_atomic(slname):
                slname = '(' + slname + ')'
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            result._latex_name = slname + r'\wedge ' + olname
        return result
