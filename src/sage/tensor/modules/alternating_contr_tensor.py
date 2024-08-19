r"""
Alternating contravariant tensors on free modules

Given a free module `M` of finite rank over a commutative ring `R`
and a positive integer `p`, an *alternating contravariant tensor of
degree* `p` is a map

.. MATH::

    a:\ \underbrace{M^*\times\cdots\times M^*}_{p\ \; \mbox{times}}
    \longrightarrow R

that (i) is multilinear and (ii) vanishes whenever any of two of its
arguments are equal (`M^*` stands for the dual of `M`).
`a` is an element of the `p`-th  exterior power of `M`, `\Lambda^p(M)`.

Alternating contravariant tensors are implemented via the class
:class:`AlternatingContrTensor`, which is a subclass of the generic tensor
class :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`.

AUTHORS:

- Eric Gourgoulhon (2017): initial version

REFERENCES:

- Chap. 23 of R. Godement : *Algebra* [God1968]_
- Chap. 15 of S. Lang : *Algebra* [Lan2002]_
"""
#******************************************************************************
#       Copyright (C) 2017 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.comp import Components, CompFullyAntiSym

class AlternatingContrTensor(FreeModuleTensor):
    r"""
    Alternating contravariant tensor on a free module of finite rank
    over a commutative ring.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank over a commutative
      ring `R`, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``degree`` -- positive integer; the degree `p` of the alternating
      contravariant tensor (i.e. the tensor rank)
    - ``name`` -- (default: ``None``) string; name given to the
      alternating contravariant tensor
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the alternating contravariant tensor; if none is provided,
      ``name`` is used

    EXAMPLES:

    Alternating contravariant tensor of degree 2 on a rank-3 free module::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: e = M.basis('e')
        sage: a = M.alternating_contravariant_tensor(2, name='a') ; a
        Alternating contravariant tensor a of degree 2 on the Rank-3
         free module M over the Integer Ring
        sage: type(a)
        <class 'sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule_with_category.element_class'>
        sage: a.parent()
        2nd exterior power of the Rank-3 free module M over the Integer Ring
        sage: a[1,2], a[2,3] = 4, -3
        sage: a.display(e)
        a = 4 e_1∧e_2 - 3 e_2∧e_3

    The alternating contravariant tensor acting on the dual basis elements::

        sage: f = e.dual_basis(); f
        Dual basis (e^1,e^2,e^3) on the Rank-3 free module M over the
         Integer Ring
        sage: a(f[1],f[2])
        4
        sage: a(f[1],f[3])
        0
        sage: a(f[2],f[3])
        -3
        sage: a(f[2],f[1])
        -4

    An alternating contravariant tensor of degree 1 is an element
    of the module `M`::

        sage: b = M.alternating_contravariant_tensor(1, name='b') ; b
        Element b of the Rank-3 free module M over the Integer Ring
        sage: b[:] = [2,-1,3]  # components w.r.t. the module's default basis (e)
        sage: b.parent() is M
        True

    The standard tensor operations apply to alternating contravariant
    tensors, like the extraction of components with respect to a
    given basis::

        sage: a[e,1,2]
        4
        sage: a[1,2]  # since e is the module's default basis
        4
        sage: all( a[i,j] == - a[j,i] for i in {1,2,3} for j in {1,2,3} )
        True

    the tensor product::

        sage: c = b*b ; c
        Type-(2,0) tensor b⊗b on the Rank-3 free module M over the
         Integer Ring
        sage: c.symmetries()
        symmetry: (0, 1); no antisymmetry
        sage: c.parent()
        Free module of type-(2,0) tensors on the Rank-3 free module M
         over the Integer Ring
        sage: c.display(e)
        b⊗b = 4 e_1⊗e_1 - 2 e_1⊗e_2 + 6 e_1⊗e_3 - 2 e_2⊗e_1 + e_2⊗e_2
         - 3 e_2⊗e_3 + 6 e_3⊗e_1 - 3 e_3⊗e_2 + 9 e_3⊗e_3

    the contractions::

        sage: w = f[1] + f[2] + f[3]  # a linear form
        sage: s = a.contract(w) ; s
        Element of the Rank-3 free module M over the Integer Ring
        sage: s.display(e)
        4 e_1 - 7 e_2 + 3 e_3

    or tensor arithmetics::

        sage: s = 3*a + c ; s
        Type-(2,0) tensor on the Rank-3 free module M over the Integer Ring
        sage: s.parent()
        Free module of type-(2,0) tensors on the Rank-3 free module M
         over the Integer Ring
        sage: s.display(e)
        4 e_1⊗e_1 + 10 e_1⊗e_2 + 6 e_1⊗e_3 - 14 e_2⊗e_1 + e_2⊗e_2
         - 12 e_2⊗e_3 + 6 e_3⊗e_1 + 6 e_3⊗e_2 + 9 e_3⊗e_3

    Note that tensor arithmetics preserves the alternating character if
    both operands are alternating::

        sage: s = a - 2*a ; s
        Alternating contravariant tensor of degree 2 on the Rank-3 free
         module M over the Integer Ring
        sage: s.parent() # note the difference with s = 3*a + c above
        2nd exterior power of the Rank-3 free module M over the Integer
         Ring
        sage: s == -a
        True

    An operation specific to alternating contravariant tensors is of
    course the exterior product::

        sage: s = a.wedge(b) ; s
        Alternating contravariant tensor a∧b of degree 3 on the Rank-3 free
         module M over the Integer Ring
        sage: s.parent()
        3rd exterior power of the Rank-3 free module M over the Integer Ring
        sage: s.display(e)
        a∧b = 6 e_1∧e_2∧e_3
        sage: s[1,2,3] == a[1,2]*b[3] + a[2,3]*b[1] + a[3,1]*b[2]
        True

    The exterior product is nilpotent on module elements::

        sage: s = b.wedge(b) ; s
        Alternating contravariant tensor b∧b of degree 2 on the Rank-3 free
         module M over the Integer Ring
        sage: s.display(e)
        b∧b = 0
    """
    def __init__(self, fmodule, degree, name=None, latex_name=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.tensor.modules.alternating_contr_tensor import AlternatingContrTensor
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = AlternatingContrTensor(M, 2, name='a')
            sage: a[e,0,1] = 2
            sage: TestSuite(a).run(skip='_test_category') # see below

        In the above test suite, _test_category fails because a is not an
        instance of a.parent().category().element_class. Actually alternating
        tensors must be constructed via ExtPowerFreeModule.element_class and
        not by a direct call to AlternatingContrTensor::

            sage: a1 = M.exterior_power(2).element_class(M, 2, name='a')
            sage: a1[e,0,1] = 2
            sage: TestSuite(a1).run()
        """
        FreeModuleTensor.__init__(self, fmodule, (degree,0), name=name,
                                  latex_name=latex_name,
                                  antisym=range(degree),
                                  parent=fmodule.exterior_power(degree))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.alternating_contravariant_tensor(1)
            Element of the Rank-3 free module M over the Integer Ring
            sage: M.alternating_contravariant_tensor(1, name='a')
            Element a of the Rank-3 free module M over the Integer Ring
            sage: M.alternating_contravariant_tensor(2)
            Alternating contravariant tensor of degree 2 on the Rank-3 free
             module M over the Integer Ring
            sage: M.alternating_contravariant_tensor(2, name='a')
            Alternating contravariant tensor a of degree 2 on the Rank-3 free
             module M over the Integer Ring
        """
        if self._tensor_rank == 1:
            description = "Element "
            if self._name is not None:
                description += self._name + " "
            description += "of the {}".format(self._fmodule)
        else:
            description = "Alternating contravariant tensor "
            if self._name is not None:
                description += self._name + " "
            description += "of degree {} ".format(self._tensor_rank)
            description += "on the {}".format(self._fmodule)
        return description

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``, on the same module
        and of the same degree.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_contravariant_tensor(2, name='a')
            sage: a._new_instance()
            Alternating contravariant tensor of degree 2 on the Rank-3 free
             module M over the Integer Ring
            sage: a._new_instance().parent() is a.parent()
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
            sage: a = M.alternating_contravariant_tensor(2, name='a')
            sage: a._new_comp(e)
            Fully antisymmetric 2-indices components w.r.t. Basis (e_0,e_1,e_2)
             on the Rank-3 free module M over the Integer Ring
            sage: a = M.alternating_contravariant_tensor(1)
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

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_contravariant_tensor(2, name='a')
            sage: a.degree()
            2
        """
        return self._tensor_rank

    def display(self, basis=None, format_spec=None):
        r"""
        Display the alternating contravariant tensor ``self`` in terms
        of its expansion w.r.t. a given module basis.

        The expansion is actually performed onto exterior products of
        elements of ``basis`` (see examples below). The output is either
        text-formatted (console mode) or LaTeX-formatted (notebook mode).

        INPUT:

        - ``basis`` -- (default: ``None``) basis of the free module with
          respect to which ``self`` is expanded; if none is provided,
          the module's default basis is assumed
        - ``format_spec`` -- (default: ``None``) format specification
          passed to ``self._fmodule._output_formatter`` to format the
          output

        EXAMPLES:

        Display of an alternating contravariant tensor of degree 2 on a rank-3
        free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.alternating_contravariant_tensor(2, 'a', latex_name=r'\alpha')
            sage: a[0,1], a[0,2], a[1,2] = 3, 2, -1
            sage: a.display()
            a = 3 e_0∧e_1 + 2 e_0∧e_2 - e_1∧e_2
            sage: latex(a.display())  # display in the notebook
            \alpha = 3 e_{0}\wedge e_{1} + 2 e_{0}\wedge e_{2} -e_{1}\wedge e_{2}

        Display of an alternating contravariant tensor of degree 3 on a rank-3
        free module::

            sage: b = M.alternating_contravariant_tensor(3, 'b')
            sage: b[0,1,2] = 4
            sage: b.display()
            b = 4 e_0∧e_1∧e_2
            sage: latex(b.display())
            b = 4 e_{0}\wedge e_{1}\wedge e_{2}

        Display of a vanishing alternating contravariant tensor::

            sage: b[0,1,2] = 0  # the only independent component set to zero
            sage: b.is_zero()
            True
            sage: b.display()
            b = 0
            sage: latex(b.display())
            b = 0
            sage: b[0,1,2] = 4  # value restored for what follows

        Display in a basis which is not the default one::

            sage: aut = M.automorphism(matrix=[[0,1,0], [0,0,-1], [1,0,0]],
            ....:                      basis=e)
            sage: f = e.new_basis(aut, 'f')
            sage: a.display(f)
            a = -2 f_0∧f_1 - f_0∧f_2 - 3 f_1∧f_2
            sage: a.disp(f)     # shortcut notation
            a = -2 f_0∧f_1 - f_0∧f_2 - 3 f_1∧f_2
            sage: b.display(f)
            b = -4 f_0∧f_1∧f_2

        The output format can be set via the argument ``output_formatter``
        passed at the module construction::

            sage: N = FiniteRankFreeModule(QQ, 3, name='N', start_index=1,
            ....:                   output_formatter=Rational.numerical_approx)
            sage: e = N.basis('e')
            sage: a = N.alternating_contravariant_tensor(2, 'a')
            sage: a[1,2], a[1,3], a[2,3] = 1/3, 5/2, 4
            sage: a.display()  # default format (53 bits of precision)
            a = 0.333333333333333 e_1∧e_2 + 2.50000000000000 e_1∧e_3
             + 4.00000000000000 e_2∧e_3

        The output format is then controlled by the argument ``format_spec`` of
        the method :meth:`display`::

            sage: a.display(format_spec=10)  # 10 bits of precision
            a = 0.33 e_1∧e_2 + 2.5 e_1∧e_3 + 4.0 e_2∧e_3
        """
        from sage.misc.latex import latex
        from sage.typeset.unicode_characters import unicode_wedge
        from .format_utilities import is_atomic, FormattedExpansion
        basis, format_spec = self._preparse_display(basis=basis,
                                                    format_spec=format_spec)
        comp = self.comp(basis)
        terms_txt = []
        terms_latex = []
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
                for k in range(self._tensor_rank):
                    bases_txt.append(basis[ind[k]]._name)
                    bases_latex.append(latex(basis[ind[k]]))
                basis_term_txt = unicode_wedge.join(bases_txt)
                basis_term_latex = r'\wedge '.join(bases_latex)
                coef_txt = repr(coef)
                if coef_txt == '1':
                    terms_txt.append(basis_term_txt)
                    terms_latex.append(basis_term_latex)
                elif coef_txt == '-1':
                    terms_txt.append('-' + basis_term_txt)
                    terms_latex.append('-' + basis_term_latex)
                else:
                    coef_latex = latex(coef)
                    if is_atomic(coef_txt):
                        terms_txt.append(coef_txt + ' ' + basis_term_txt)
                    else:
                        terms_txt.append('(' + coef_txt + ') ' +
                                         basis_term_txt)
                    if is_atomic(coef_latex):
                        terms_latex.append(coef_latex + basis_term_latex)
                    else:
                        terms_latex.append(r'\left(' + coef_latex +
                                           r'\right)' + basis_term_latex)
        if not terms_txt:
            expansion_txt = '0'
        else:
            expansion_txt = terms_txt[0]
            for term in terms_txt[1:]:
                if term[0] == '-':
                    expansion_txt += ' - ' + term[1:]
                else:
                    expansion_txt += ' + ' + term
        if not terms_latex:
            expansion_latex = '0'
        else:
            expansion_latex = terms_latex[0]
            for term in terms_latex[1:]:
                if term[0] == '-':
                    expansion_latex += term
                else:
                    expansion_latex += '+' + term
        if self._name is None:
            resu_txt = expansion_txt
        else:
            resu_txt = self._name + ' = ' + expansion_txt
        if self._latex_name is None:
            resu_latex = expansion_latex
        else:
            resu_latex = latex(self) + ' = ' + expansion_latex
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display

    def wedge(self, other):
        r"""
        Exterior product of ``self`` with the alternating contravariant
        tensor ``other``.

        INPUT:

        - ``other`` -- an alternating contravariant tensor

        OUTPUT:

        - instance of :class:`AlternatingContrTensor` representing the
          exterior product ``self ∧ other``

        EXAMPLES:

        Exterior product of two module elements::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M([1,-3,4], basis=e, name='A')
            sage: b = M([2,-1,2], basis=e, name='B')
            sage: c = a.wedge(b) ; c
            Alternating contravariant tensor A∧B of degree 2 on the Rank-3
             free module M over the Integer Ring
            sage: c.display()
            A∧B = 5 e_0∧e_1 - 6 e_0∧e_2 - 2 e_1∧e_2
            sage: latex(c)
            A\wedge B
            sage: latex(c.display())
            A\wedge B = 5 e_{0}\wedge e_{1} -6 e_{0}\wedge e_{2}
             -2 e_{1}\wedge e_{2}

        Test of the computation::

            sage: a.wedge(b) == a*b - b*a
            True

        Exterior product of a module element and an alternating contravariant
        tensor of degree 2::

            sage: d = M([-1,2,4], basis=e, name='D')
            sage: s = d.wedge(c) ; s
            Alternating contravariant tensor D∧A∧B of degree 3 on the Rank-3
             free module M over the Integer Ring
            sage: s.display()
            D∧A∧B = 34 e_0∧e_1∧e_2

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
        from sage.typeset.unicode_characters import unicode_wedge
        from .format_utilities import is_atomic
        if not isinstance(other, AlternatingContrTensor):
            raise TypeError("the second argument for the exterior product " +
                            "must be an alternating contravariant tensor")
        if other._tensor_rank == 0:
            return other*self
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
        for ind_s, val_s in cmp_s._comp.items():
            for ind_o, val_o in cmp_o._comp.items():
                ind_r = ind_s + ind_o
                if len(ind_r) == len(set(ind_r)): # all indices are different
                    cmp_r[[ind_r]] += val_s * val_o
        result = fmodule.alternating_contravariant_tensor(rank_r)
        result._components[basis] = cmp_r
        if self._name is not None and other._name is not None:
            sname = self._name
            oname = other._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            result._name = sname + unicode_wedge + oname
        if self._latex_name is not None and other._latex_name is not None:
            slname = self._latex_name
            olname = other._latex_name
            if not is_atomic(slname):
                slname = r'\left(' + slname + r'\right)'
            if not is_atomic(olname):
                olname = r'\left(' + olname + r'\right)'
            result._latex_name = slname + r'\wedge ' + olname
        return result

    def interior_product(self, form):
        r"""
        Interior product with an alternating form.

        If ``self`` is an alternating contravariant tensor `A` of degree `p`
        and `B` is an alternating form of degree `q\geq p` on the same free
        module, the interior product of `A` by `B` is the alternating form
        `\iota_A B` of degree `q-p` defined by

        .. MATH::

            (\iota_A B)_{i_1\ldots i_{q-p}} = A^{k_1\ldots k_p}
                            B_{k_1\ldots k_p i_1\ldots i_{q-p}}

        .. NOTE::

            ``A.interior_product(B)`` yields the same result as
            ``A.contract(0,..., p-1, B, 0,..., p-1)`` (cf.
            :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.contract`),
            but ``interior_product`` is more efficient, the alternating
            character of `A` being not used to reduce the computation in
            :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.contract`

        INPUT:

        - ``form`` -- alternating form `B` (instance of
          :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`);
          the degree of `B` must be at least equal to the degree of ``self``

        OUTPUT:

        - element of the base ring (case `p=q`) or
          :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`
          (case `p<q`) representing the interior product `\iota_A B`, where `A`
          is ``self``

        .. SEEALSO::

            :meth:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm.interior_product`
            for the interior product of an alternating form by an alternating
            contravariant tensor

        EXAMPLES:

        Let us consider a rank-4 free module::

            sage: M = FiniteRankFreeModule(ZZ, 4, name='M', start_index=1)
            sage: e = M.basis('e')

        and various interior products on it, starting with a module element
        (``p=1``) and a linear form (``q=1``)::

            sage: a = M([-2,1,2,3], basis=e, name='A')
            sage: b = M.linear_form(name='B')
            sage: b[:] = [2, 0, -3, 4]
            sage: c = a.interior_product(b); c
            2
            sage: c == a.contract(b)
            True

        Case  ``p=1`` and ``q=3``::

            sage: b = M.alternating_form(3, name='B')
            sage: b[1,2,3], b[1,2,4], b[1,3,4], b[2,3,4] = 3, -1, 2, 5
            sage: c = a.interior_product(b); c
            Alternating form i_A B of degree 2 on the Rank-4 free module M over the Integer Ring
            sage: c.display()
            i_A B = 3 e^1∧e^2 + 3 e^1∧e^3 - 3 e^1∧e^4 + 9 e^2∧e^3 - 8 e^2∧e^4 + e^3∧e^4
            sage: latex(c)
            \iota_{A} B
            sage: c == a.contract(b)
            True

        Case  ``p=2`` and ``q=3``::

            sage: a = M.alternating_contravariant_tensor(2, name='A')
            sage: a[1,2], a[1,3], a[1,4] = 2, -5, 3
            sage: a[2,3], a[2,4], a[3,4] = -1, 4, 2
            sage: c = a.interior_product(b); c
            Linear form i_A B on the Rank-4 free module M over the Integer Ring
            sage: c.display()
            i_A B = -6 e^1 + 56 e^2 - 40 e^3 - 34 e^4
            sage: c == a.contract(0, 1, b, 0, 1)  # contraction on all indices of a
            True

        Case  ``p=2`` and ``q=4``::

            sage: b = M.alternating_form(4, name='B')
            sage: b[1,2,3,4] = 5
            sage: c = a.interior_product(b); c
            Alternating form i_A B of degree 2 on the Rank-4 free module M over the Integer Ring
            sage: c.display()
            i_A B = 20 e^1∧e^2 - 40 e^1∧e^3 - 10 e^1∧e^4 + 30 e^2∧e^3 + 50 e^2∧e^4 + 20 e^3∧e^4
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case  ``p=2`` and ``q=2``::

            sage: b = M.alternating_form(2)
            sage: b[1,2], b[1,3], b[1,4] = 6, 0, -2
            sage: b[2,3], b[2,4], b[3,4] = 2, 3, 4
            sage: c = a.interior_product(b); c
            48
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case  ``p=3`` and ``q=3``::

            sage: a = M.alternating_contravariant_tensor(3, name='A')
            sage: a[1,2,3], a[1,2,4], a[1,3,4], a[2,3,4] = -3, 2, 8, -5
            sage: b = M.alternating_form(3, name='B')
            sage: b[1,2,3], b[1,2,4], b[1,3,4], b[2,3,4] = 3, -1, 2, 5
            sage: c = a.interior_product(b); c
            -120
            sage: c == a.contract(0, 1, 2, b, 0, 1, 2)
            True

        Case  ``p=3`` and ``q=4``::

            sage: b = M.alternating_form(4, name='B')
            sage: b[1,2,3,4] = 5
            sage: c = a.interior_product(b); c
            Linear form i_A B on the Rank-4 free module M over the Integer Ring
            sage: c.display()
            i_A B = 150 e^1 + 240 e^2 - 60 e^3 - 90 e^4
            sage: c == a.contract(0, 1, 2, b, 0, 1, 2)
            True

        Case  ``p=4`` and ``q=4``::

            sage: a = M.alternating_contravariant_tensor(4, name='A')
            sage: a[1,2,3,4] = -2
            sage: c = a.interior_product(b); c
            -240
            sage: c == a.contract(0, 1, 2, 3, b, 0, 1, 2, 3)
            True
        """
        from .format_utilities import is_atomic
        from .free_module_alt_form import FreeModuleAltForm
        if not isinstance(form, FreeModuleAltForm):
            raise TypeError("{} is not an alternating form".format(form))
        p_res = form._tensor_rank - self._tensor_rank  # degree of the result
        if self._tensor_rank == 1:
            # Case p = 1:
            res = self.contract(form)  # contract() deals efficiently with
                                       # the antisymmetry for p = 1
        else:
            # Case p > 1:
            if form._fmodule != self._fmodule:
                raise ValueError("{} is not defined on the same ".format(form) +
                                 "module as the {}".format(self))
            if form._tensor_rank < self._tensor_rank:
                raise ValueError("the degree of the {} is lower ".format(form) +
                                 "than that of the {}".format(self))
            # Interior product at the component level:
            basis = self.common_basis(form)
            if basis is None:
                raise ValueError("no common basis for the interior product")
            comp = self._components[basis].interior_product(
                                                       form._components[basis])
            if p_res == 0:
                res = comp  # result is a scalar
            else:
                res = self._fmodule.tensor_from_comp((0, p_res), comp)
        # Name of the result
        res_name = None
        if self._name is not None and form._name is not None:
            sname = self._name
            oname = form._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            res_name = 'i_' + sname + ' ' + oname
        res_latex_name = None
        if self._latex_name is not None and form._latex_name is not None:
            slname = self._latex_name
            olname = form._latex_name
            if not is_atomic(olname):
                olname = r'\left(' + olname + r'\right)'
            res_latex_name = r'\iota_{' + slname + '} ' + olname
        if res_name:
            try:  # there is no guarantee that the result has set_name
                  # and is mutable
                res.set_name(res_name, latex_name=res_latex_name)
            except (AttributeError, TypeError, ValueError):
                pass
        return res
