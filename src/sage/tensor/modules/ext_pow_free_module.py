r"""
Exterior powers of free modules

Given a free module `M` of finite rank over a commutative ring `R`
and a positive integer `p`, the `p`-*th exterior power of* `M`
is the set `\Lambda^p(M)` of all alternating contravariant tensors of
degree `p` on `M`, i.e. of all multilinear maps

.. MATH::

    \underbrace{M^*\times\cdots\times M^*}_{p\ \; \mbox{times}}
    \longrightarrow R

that vanish whenever any of two of their arguments are equal
(`M^*` stands for the dual of `M`).
Note that `\Lambda^1(M) = M`. The exterior power
`\Lambda^p(M)` is a free module of rank `\binom{n}{p}` over `R`,
where `n` is the rank of `M`.

Similarly, the `p`-*th exterior power of the dual of* `M`
is the set `\Lambda^p(M^*)` of all alternating forms of degree `p` on
`M`, i.e. of all multilinear maps

.. MATH::

    \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
    \longrightarrow R

that vanish whenever any of two of their arguments are equal.
Note that `\Lambda^1(M^*) = M^*` (the dual of `M`). The exterior power
`\Lambda^p(M^*)` is a free module of rank `\binom{n}{p}` over `R`,
where `n` is the rank of `M`.

The class :class:`ExtPowerFreeModule` implements `\Lambda^p(M)`, while
the class :class:`ExtPowerDualFreeModule` implements `\Lambda^p(M^*)`.

AUTHORS:

- Eric Gourgoulhon: initial version, regarding `\Lambda^p(M^*)` only
  (2015); add class for `\Lambda^p(M)` (2017)


REFERENCES:

- \K. Conrad: *Exterior powers* [Con2013]_
- Chap. 19 of S. Lang: *Algebra* [Lan2002]_
"""
#******************************************************************************
#       Copyright (C) 2017 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule_abstract
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.alternating_contr_tensor import AlternatingContrTensor
from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm

class ExtPowerFreeModule(FiniteRankFreeModule_abstract):
    r"""
    Exterior power of a free module of finite rank over a commutative
    ring.

    Given a free module `M` of finite rank over a commutative ring `R`
    and a positive integer `p`, the `p`-*th exterior power of* `M` is
    the set `\Lambda^p(M)` of all alternating contravariant tensors of
    degree `p` on `M`, i.e. of all multilinear maps

    .. MATH::

        \underbrace{M^*\times\cdots\times M^*}_{p\ \; \mbox{times}}
        \longrightarrow R

    that vanish whenever any of two of their arguments are equal.
    Note that `\Lambda^1(M) = M`.

    `\Lambda^p(M)` is a free module of rank `\binom{n}{p}` over
    `R`, where `n` is the rank of `M`.
    Accordingly, the class :class:`ExtPowerFreeModule` inherits from the
    class
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule_abstract`.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``degree`` -- positive integer; the degree `p` of the alternating
      elements
    - ``name`` -- (default: ``None``) string; name given to `\Lambda^p(M)`
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote `\Lambda^p(M)`

    EXAMPLES:

    2nd exterior power of the dual of a free `\ZZ`-module of rank 3::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
        sage: A = ExtPowerFreeModule(M, 2) ; A
        2nd exterior power of the Rank-3 free module M over the
         Integer Ring

    Instead of importing ExtPowerFreeModule in the global name space, it is
    recommended to use the module's method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.exterior_power`::

        sage: A = M.exterior_power(2) ; A
        2nd exterior power of the Rank-3 free module M over the
         Integer Ring
        sage: latex(A)
        \Lambda^{2}\left(M\right)

    ``A`` is a module (actually a free module) over `\ZZ`::

        sage: A.category()
        Category of finite dimensional modules over Integer Ring
        sage: A in Modules(ZZ)
        True
        sage: A.rank()
        3
        sage: A.base_ring()
        Integer Ring
        sage: A.base_module()
        Rank-3 free module M over the Integer Ring

    ``A`` is a *parent* object, whose elements are alternating
    contravariant tensors, represented by instances of the class
    :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`::

        sage: a = A.an_element() ; a
        Alternating contravariant tensor of degree 2 on the Rank-3 free
         module M over the Integer Ring
        sage: a.display() # expansion with respect to M's default basis (e)
        e_0∧e_1
        sage: from sage.tensor.modules.alternating_contr_tensor import AlternatingContrTensor
        sage: isinstance(a, AlternatingContrTensor)
        True
        sage: a in A
        True
        sage: A.is_parent_of(a)
        True

    Elements can be constructed from ``A``. In particular, 0 yields
    the zero element of ``A``::

        sage: A(0)
        Alternating contravariant tensor zero of degree 2 on the Rank-3
         free module M over the Integer Ring
        sage: A(0) is A.zero()
        True

    while nonzero elements are constructed by providing their components in a
    given basis::

        sage: e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: comp = [[0,3,-1],[-3,0,4],[1,-4,0]]
        sage: a = A(comp, basis=e, name='a') ; a
        Alternating contravariant tensor a of degree 2 on the Rank-3
         free module M over the Integer Ring
        sage: a.display(e)
        a = 3 e_0∧e_1 - e_0∧e_2 + 4 e_1∧e_2

    An alternative is to construct the alternating contravariant tensor from an
    empty list of components and to set the nonzero components afterwards::

        sage: a = A([], name='a')
        sage: a.set_comp(e)[0,1] = 3
        sage: a.set_comp(e)[0,2] = -1
        sage: a.set_comp(e)[1,2] = 4
        sage: a.display(e)
        a = 3 e_0∧e_1 - e_0∧e_2 + 4 e_1∧e_2

    The exterior powers are unique::

        sage: A is M.exterior_power(2)
        True

    The exterior power `\Lambda^1(M)` is nothing but `M`::

        sage: M.exterior_power(1) is M
        True

    For a degree `p\geq 2`, there is a coercion
    `\Lambda^p(M)\rightarrow T^{(p,0)}(M)`::

        sage: T20 = M.tensor_module(2,0) ; T20
        Free module of type-(2,0) tensors on the Rank-3 free module M
         over the Integer Ring
        sage: T20.has_coerce_map_from(A)
        True

    Of course, there is no coercion in the reverse direction::

        sage: A.has_coerce_map_from(T20)
        False

    The coercion map `\Lambda^2(M)\rightarrow T^{(2,0)}(M)` in action::

        sage: ta = T20(a) ; ta
        Type-(2,0) tensor a on the Rank-3 free module M over the Integer Ring
        sage: ta.display(e)
        a = 3 e_0⊗e_1 - e_0⊗e_2 - 3 e_1⊗e_0 + 4 e_1⊗e_2 + e_2⊗e_0 - 4 e_2⊗e_1
        sage: a.display(e)
        a = 3 e_0∧e_1 - e_0∧e_2 + 4 e_1∧e_2
        sage: ta.symmetries()  # the antisymmetry is of course preserved
        no symmetry;  antisymmetry: (0, 1)
        sage: ta == a  # equality as type-(2,0) tensors
        True
    """

    Element = AlternatingContrTensor

    def __init__(self, fmodule, degree, name=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = ExtPowerFreeModule(M, 2) ; A
            2nd exterior power of the Rank-3 free module M over the
             Integer Ring
            sage: TestSuite(A).run()
        """
        from sage.arith.misc import binomial
        from sage.typeset.unicode_characters import unicode_bigwedge
        self._fmodule = fmodule
        self._degree = ZZ(degree)
        rank = binomial(fmodule._rank, degree)
        if name is None and fmodule._name is not None:
            name = unicode_bigwedge + r'^{}('.format(degree) \
                   + fmodule._name + ')'
        if latex_name is None and fmodule._latex_name is not None:
            latex_name = r'\Lambda^{' + str(degree) + r'}\left(' \
                         + fmodule._latex_name + r'\right)'
        super().__init__(fmodule._ring, rank,
                         name=name, latex_name=latex_name)
        fmodule._all_modules.add(self)

    def construction(self):
        r"""
        Return the functorial construction of ``self``.

        This implementation just returns ``None``, as no functorial construction is implemented.

        TESTS::

            sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = ExtPowerFreeModule(M, 2)
            sage: A.construction() is None
            True
        """
        # No construction until https://github.com/sagemath/sage/issues/30242
        # makes this a quotient of TensorFreeModule
        return None

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an alternating contravariant tensor.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.exterior_power(2)
            sage: a = A._element_constructor_(0) ; a
            Alternating contravariant tensor zero of degree 2 on the
             Rank-3 free module M over the Integer Ring
            sage: a = A._element_constructor_([], name='a') ; a
            Alternating contravariant tensor a of degree 2 on the Rank-3
             free module M over the Integer Ring
            sage: a[e,0,2], a[e,1,2] = 3, -1
            sage: a.display()
            a = 3 e_0∧e_2 - e_1∧e_2
        """
        if isinstance(comp, (int, Integer)) and comp == 0:
            return self.zero()
        resu = self.element_class(self._fmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) alternating contravariant tensor.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 4, name='M')
            sage: e = M.basis('e')
            sage: a = M.exterior_power(2)._an_element_() ; a
            Alternating contravariant tensor of degree 2 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e_0∧e_1
            sage: a = M.exterior_power(3)._an_element_() ; a
            Alternating contravariant tensor of degree 3 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e_0∧e_1∧e_2
            sage: a = M.exterior_power(4)._an_element_() ; a
            Alternating contravariant tensor of degree 4 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e_0∧e_1∧e_2∧e_3

        TESTS:

        When the base module has no default basis, a default
        basis will be set for it::

            sage: M2 = FiniteRankFreeModule(QQ, 4, name='M2')
            sage: a = M2.exterior_power(2)._an_element_(); a
            Alternating contravariant tensor of degree 2
            on the 4-dimensional vector space M2 over the Rational Field
            sage: a + a
            Alternating contravariant tensor of degree 2
            on the 4-dimensional vector space M2 over the Rational Field
            sage: M2.default_basis()
            Basis (e_0,e_1,e_2,e_3) on the 4-dimensional vector space M2 over the Rational Field
        """
        resu = self.element_class(self._fmodule, self._degree)
        # Make sure that the base module has a default basis
        self._fmodule.an_element()
        sindex = self._fmodule._sindex
        ind = [sindex + i for i in range(resu._tensor_rank)]
        resu.set_comp()[ind] = self._fmodule._ring.an_element()
        return resu

    #### End of parent methods

    @cached_method
    def zero(self):
        r"""
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.exterior_power(2)
            sage: A.zero()
            Alternating contravariant tensor zero of degree 2 on the Rank-3 free
             module M over the Integer Ring
            sage: A(0) is A.zero()
            True
        """
        resu = self._element_constructor_(name='zero', latex_name='0')
        for basis in self._fmodule._known_bases:
            resu._add_comp_unsafe(basis)
            # (since new components are initialized to zero)
        resu._is_zero = True # This element is certainly zero
        resu.set_immutable()
        return resu

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: M.exterior_power(2)._repr_()
            '2nd exterior power of the Rank-5 free module M over the Integer Ring'
            sage: M.exterior_power(3)._repr_()
            '3rd exterior power of the Rank-5 free module M over the Integer Ring'
            sage: M.exterior_power(4)._repr_()
            '4th exterior power of the Rank-5 free module M over the Integer Ring'
            sage: M.exterior_power(5)._repr_()
            '5th exterior power of the Rank-5 free module M over the Integer Ring'
            sage: M.exterior_power(21)._repr_()
            '21st exterior power of the Rank-5 free module M over the Integer Ring'
        """
        description = "{}".format(self._degree.ordinal_str())
        description += " exterior power of the {}".format(self._fmodule)
        return description

    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the
          free module on which the exterior power is defined.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: A = M.exterior_power(2)
            sage: A.base_module()
            Rank-5 free module M over the Integer Ring
            sage: A.base_module() is M
            True
        """
        return self._fmodule

    def degree(self):
        r"""
        Return the degree of ``self``.

        OUTPUT:

        - integer `p` such that ``self`` is the exterior power
          `\Lambda^p(M)`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: A = M.exterior_power(2)
            sage: A.degree()
            2
            sage: M.exterior_power(4).degree()
            4
        """
        return self._degree

#***********************************************************************


class ExtPowerDualFreeModule(FiniteRankFreeModule_abstract):
    r"""
    Exterior power of the dual of a free module of finite rank
    over a commutative ring.

    Given a free module `M` of finite rank over a commutative ring `R`
    and a positive integer `p`, the `p`-*th exterior power of the dual of*
    `M` is the set `\Lambda^p(M^*)` of all alternating forms of degree
    `p` on `M`, i.e. of all multilinear maps

    .. MATH::

        \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
        \longrightarrow R

    that vanish whenever any of two of their arguments are equal.
    Note that `\Lambda^1(M^*) = M^*` (the dual of `M`).

    `\Lambda^p(M^*)` is a free module of rank `\binom{n}{p}` over
    `R`, where `n` is the rank of `M`.
    Accordingly, the class :class:`ExtPowerDualFreeModule` inherits from
    the class
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule_abstract`.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``degree`` -- positive integer; the degree `p` of the alternating
      forms
    - ``name`` -- (default: ``None``) string; name given to `\Lambda^p(M^*)`
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote `\Lambda^p(M^*)`

    EXAMPLES:

    2nd exterior power of the dual of a free `\ZZ`-module of rank 3::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerDualFreeModule
        sage: A = ExtPowerDualFreeModule(M, 2) ; A
        2nd exterior power of the dual of the Rank-3 free module M over the
         Integer Ring

    Instead of importing ExtPowerDualFreeModule in the global name space,
    it is recommended to use the module's method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.dual_exterior_power`::

        sage: A = M.dual_exterior_power(2) ; A
        2nd exterior power of the dual of the Rank-3 free module M over the
         Integer Ring
        sage: latex(A)
        \Lambda^{2}\left(M^*\right)

    ``A`` is a module (actually a free module) over `\ZZ`::

        sage: A.category()
        Category of finite dimensional modules over Integer Ring
        sage: A in Modules(ZZ)
        True
        sage: A.rank()
        3
        sage: A.base_ring()
        Integer Ring
        sage: A.base_module()
        Rank-3 free module M over the Integer Ring

    ``A`` is a *parent* object, whose elements are alternating forms,
    represented by instances of the class
    :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`::

        sage: a = A.an_element() ; a
        Alternating form of degree 2 on the Rank-3 free module M over the
         Integer Ring
        sage: a.display() # expansion with respect to M's default basis (e)
        e^0∧e^1
        sage: from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
        sage: isinstance(a, FreeModuleAltForm)
        True
        sage: a in A
        True
        sage: A.is_parent_of(a)
        True

    Elements can be constructed from ``A``. In particular, 0 yields
    the zero element of ``A``::

        sage: A(0)
        Alternating form zero of degree 2 on the Rank-3 free module M over the
         Integer Ring
        sage: A(0) is A.zero()
        True

    while nonzero elements are constructed by providing their components in a
    given basis::

        sage: e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: comp = [[0,3,-1],[-3,0,4],[1,-4,0]]
        sage: a = A(comp, basis=e, name='a') ; a
        Alternating form a of degree 2 on the Rank-3 free module M over the
         Integer Ring
        sage: a.display(e)
        a = 3 e^0∧e^1 - e^0∧e^2 + 4 e^1∧e^2

    An alternative is to construct the alternating form from an empty list of
    components and to set the nonzero components afterwards::

        sage: a = A([], name='a')
        sage: a.set_comp(e)[0,1] = 3
        sage: a.set_comp(e)[0,2] = -1
        sage: a.set_comp(e)[1,2] = 4
        sage: a.display(e)
        a = 3 e^0∧e^1 - e^0∧e^2 + 4 e^1∧e^2

    The exterior powers are unique::

        sage: A is M.dual_exterior_power(2)
        True

    The exterior power `\Lambda^1(M^*)` is nothing but `M^*`::

        sage: M.dual_exterior_power(1) is M.dual()
        True
        sage: M.dual()
        Dual of the Rank-3 free module M over the Integer Ring
        sage: latex(M.dual())
        M^*

    It also coincides with the module of type-`(0,1)` tensors::

        sage: M.dual_exterior_power(1) is M.tensor_module(0,1)
        True

    For a degree `p\geq 2`, there is a coercion map
    `\Lambda^p(M^*)\rightarrow T^{(0,p)}(M)`::

        sage: T02 = M.tensor_module(0,2) ; T02
        Free module of type-(0,2) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: T02.has_coerce_map_from(A)
        True
        sage: A.has_coerce_map_from(T02)
        False

    The coercion map `\Lambda^2(M^*)\rightarrow T^{(0,2)}(M)` in action::

        sage: ta = T02(a) ; ta
        Type-(0,2) tensor a on the Rank-3 free module M over the Integer Ring
        sage: ta.display(e)
        a = 3 e^0⊗e^1 - e^0⊗e^2 - 3 e^1⊗e^0 + 4 e^1⊗e^2 + e^2⊗e^0 - 4 e^2⊗e^1
        sage: a.display(e)
        a = 3 e^0∧e^1 - e^0∧e^2 + 4 e^1∧e^2
        sage: ta.symmetries() # the antisymmetry is of course preserved
        no symmetry;  antisymmetry: (0, 1)
    """

    Element = FreeModuleAltForm

    def __init__(self, fmodule, degree, name=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerDualFreeModule
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = ExtPowerDualFreeModule(M, 2) ; A
            2nd exterior power of the dual of the Rank-3 free module M over
             the Integer Ring
            sage: TestSuite(A).run()
        """
        from sage.arith.misc import binomial
        from sage.typeset.unicode_characters import unicode_bigwedge
        self._fmodule = fmodule
        self._degree = ZZ(degree)
        rank = binomial(fmodule._rank, degree)
        if name is None and fmodule._name is not None:
            name = unicode_bigwedge + r'^{}('.format(degree) \
                   + fmodule._name + '*)'
        if latex_name is None and fmodule._latex_name is not None:
            latex_name = r'\Lambda^{' + str(degree) + r'}\left(' \
                         + fmodule._latex_name + r'^*\right)'
        super().__init__(fmodule._ring, rank, name=name,
                         latex_name=latex_name)
        fmodule._all_modules.add(self)

    def construction(self):
        r"""
        Return the functorial construction of ``self``.

        This implementation just returns ``None``, as no functorial construction is implemented.

        TESTS::

            sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerDualFreeModule
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = ExtPowerDualFreeModule(M, 2)
            sage: A.construction() is None
            True
        """
        # No construction until https://github.com/sagemath/sage/issues/30242
        # makes this a quotient of TensorFreeModule
        return None

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an alternating form.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.dual_exterior_power(2)
            sage: a = A._element_constructor_(0) ; a
            Alternating form zero of degree 2 on the Rank-3 free module M over
             the Integer Ring
            sage: a = A._element_constructor_([], name='a') ; a
            Alternating form a of degree 2 on the Rank-3 free module M over
             the Integer Ring
            sage: a[e,0,2], a[e,1,2] = 3, -1
            sage: a.display()
            a = 3 e^0∧e^2 - e^1∧e^2
        """
        if isinstance(comp, (int, Integer)) and comp == 0:
            return self.zero()
        if isinstance(comp, FreeModuleTensor):
            # coercion of a tensor of type (0,1) to a linear form
            tensor = comp # for readability
            if tensor.tensor_type() == (0,1) and self._degree == 1 and \
                                         tensor.base_module() is self._fmodule:
                resu = self.element_class(self._fmodule, 1, name=tensor._name,
                                          latex_name=tensor._latex_name)
                for basis, comp in tensor._components.items():
                    resu._components[basis] = comp.copy()
                return resu
            else:
                raise TypeError("cannot coerce the {} ".format(tensor) +
                                "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._fmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) alternating form.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 4, name='M')
            sage: e = M.basis('e')
            sage: a = M.dual_exterior_power(2)._an_element_() ; a
            Alternating form of degree 2 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e^0∧e^1
            sage: a = M.dual_exterior_power(3)._an_element_() ; a
            Alternating form of degree 3 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e^0∧e^1∧e^2
            sage: a = M.dual_exterior_power(4)._an_element_() ; a
            Alternating form of degree 4 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e^0∧e^1∧e^2∧e^3

        TESTS:

        When the base module has no default basis, a default
        basis will be set for it::

            sage: M2 = FiniteRankFreeModule(QQ, 4, name='M2')
            sage: a = M2.dual_exterior_power(2)._an_element_(); a
            Alternating form of degree 2 on the 4-dimensional vector space M2 over the Rational Field
            sage: a + a
            Alternating form of degree 2 on the 4-dimensional vector space M2 over the Rational Field
            sage: M2.default_basis()
            Basis (e_0,e_1,e_2,e_3) on the 4-dimensional vector space M2 over the Rational Field
        """
        resu = self.element_class(self._fmodule, self._degree)
        # Make sure that the base module has a default basis
        self._fmodule.an_element()
        sindex = self._fmodule._sindex
        ind = [sindex + i for i in range(resu._tensor_rank)]
        resu.set_comp()[ind] = self._fmodule._ring.an_element()
        return resu

    #### End of parent methods

    @cached_method
    def zero(self):
        r"""
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.dual_exterior_power(2)
            sage: A.zero()
            Alternating form zero of degree 2 on the Rank-3 free module M over
             the Integer Ring
            sage: A(0) is A.zero()
            True
        """
        resu = self._element_constructor_(name='zero', latex_name='0')
        for basis in self._fmodule._known_bases:
            resu._components[basis] = resu._new_comp(basis)
            # (since new components are initialized to zero)
        resu._is_zero = True # This element is certainly zero
        resu.set_immutable()
        return resu

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: M.dual_exterior_power(2)._repr_()
            '2nd exterior power of the dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(3)._repr_()
            '3rd exterior power of the dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(4)._repr_()
            '4th exterior power of the dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(5)._repr_()
            '5th exterior power of the dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(21)._repr_()
            '21st exterior power of the dual of the Rank-5 free module M over the Integer Ring'
        """
        description = "{}".format(self._degree.ordinal_str())
        description += " exterior power of the dual of the {}".format(
                                                                 self._fmodule)
        return description

    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module on which the exterior power is defined.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: A = M.dual_exterior_power(2)
            sage: A.base_module()
            Rank-5 free module M over the Integer Ring
            sage: A.base_module() is M
            True
        """
        return self._fmodule

    def degree(self):
        r"""
        Return the degree of ``self``.

        OUTPUT:

        - integer `p` such that ``self`` is the exterior power `\Lambda^p(M^*)`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: A = M.dual_exterior_power(2)
            sage: A.degree()
            2
            sage: M.dual_exterior_power(4).degree()
            4
        """
        return self._degree
