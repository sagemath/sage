r"""
Tensor products of free modules

The class :class:`TensorFreeModule` implements tensor products of the type

.. MATH::

    T^{(k,l)}(M) = \underbrace{M\otimes\cdots\otimes M}_{k\ \; \mbox{times}}
    \otimes \underbrace{M^*\otimes\cdots\otimes M^*}_{l\ \; \mbox{times}},

where `M` is a free module of finite rank over a commutative ring `R` and
`M^*=\mathrm{Hom}_R(M,R)` is the dual of `M`.
Note that `T^{(1,0)}(M) = M` and  `T^{(0,1)}(M) = M^*`.

Thanks to the canonical isomorphism `M^{**} \simeq M` (which holds since `M`
is a free module of finite rank), `T^{(k,l)}(M)` can be identified with the
set of tensors of type `(k,l)` defined as multilinear maps

.. MATH::

    \underbrace{M^*\times\cdots\times M^*}_{k\ \; \mbox{times}}
    \times \underbrace{M\times\cdots\times M}_{l\ \; \mbox{times}}
    \longrightarrow R

Accordingly, :class:`TensorFreeModule` is a Sage *parent* class, whose
*element* class is
:class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`.

`T^{(k,l)}(M)` is itself a free module over `R`, of rank `n^{k+l}`, `n`
being the rank of `M`. Accordingly the class :class:`TensorFreeModule`
inherits from the class
:class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`.

.. TODO::

    implement more general tensor products, i.e. tensor product of the type
    `M_1\otimes\cdots\otimes M_n`, where the `M_i`'s are `n` free modules of
    finite rank over the same ring `R`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

REFERENCES:

- \K. Conrad: *Tensor products* [Con2015]_
- Chap. 21 (Exer. 4) of R. Godement: *Algebra* [God1968]_
- Chap. 16 of S. Lang: *Algebra* [Lan2002]_

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.alternating_contr_tensor import AlternatingContrTensor
from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
from sage.tensor.modules.free_module_morphism import \
                                                   FiniteRankFreeModuleMorphism
from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism

class TensorFreeModule(FiniteRankFreeModule):
    r"""
    Class for the free modules over a commutative ring `R` that are
    tensor products of a given free module `M` over `R` with itself and its
    dual `M^*`:

    .. MATH::

        T^{(k,l)}(M) = \underbrace{M\otimes\cdots\otimes M}_{k\ \; \mbox{times}}
        \otimes \underbrace{M^*\otimes\cdots\otimes M^*}_{l\ \; \mbox{times}}

    As recalled above, `T^{(k,l)}(M)` can be canonically identified with the
    set of tensors of type `(k,l)` on `M`.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank over a commutative ring
      `R`, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``tensor_type`` -- pair ``(k, l)`` with ``k`` being the contravariant
      rank and ``l`` the covariant rank
    - ``name`` -- (default: ``None``) string; name given to the tensor module
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      tensor module; if none is provided, it is set to ``name``

    EXAMPLES:

    Set of tensors of type `(1,2)` on a free `\ZZ`-module of rank 3::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: from sage.tensor.modules.tensor_free_module import TensorFreeModule
        sage: T = TensorFreeModule(M, (1,2)) ; T
        Free module of type-(1,2) tensors on the
         Rank-3 free module M over the Integer Ring

    Instead of importing TensorFreeModule in the global name space, it is
    recommended to use the module's method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.tensor_module`::

        sage: T = M.tensor_module(1,2) ; T
        Free module of type-(1,2) tensors on the
         Rank-3 free module M over the Integer Ring
        sage: latex(T)
        T^{(1, 2)}\left(M\right)

    The module ``M`` itself is considered as the set of tensors of
    type `(1,0)`::

        sage: M is M.tensor_module(1,0)
        True

    ``T`` is a module (actually a free module) over `\ZZ`::

        sage: T.category()
        Category of finite dimensional modules over Integer Ring
        sage: T in Modules(ZZ)
        True
        sage: T.rank()
        27
        sage: T.base_ring()
        Integer Ring
        sage: T.base_module()
        Rank-3 free module M over the Integer Ring

    ``T`` is a *parent* object, whose elements are instances of
    :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`::

        sage: t = T.an_element() ; t
        Type-(1,2) tensor on the Rank-3 free module M over the Integer Ring
        sage: from sage.tensor.modules.free_module_tensor import FreeModuleTensor
        sage: isinstance(t, FreeModuleTensor)
        True
        sage: t in T
        True
        sage: T.is_parent_of(t)
        True

    Elements can be constructed from ``T``. In particular, 0 yields
    the zero element of ``T``::

        sage: T(0)
        Type-(1,2) tensor zero on the Rank-3 free module M over the Integer Ring
        sage: T(0) is T.zero()
        True

    while non-zero elements are constructed by providing their components in
    a given basis::

        sage: e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: comp = [[[i-j+k for k in range(3)] for j in range(3)] for i in range(3)]
        sage: t = T(comp, basis=e, name='t') ; t
        Type-(1,2) tensor t on the Rank-3 free module M over the Integer Ring
        sage: t.comp(e)[:]
        [[[0, 1, 2], [-1, 0, 1], [-2, -1, 0]],
         [[1, 2, 3], [0, 1, 2], [-1, 0, 1]],
         [[2, 3, 4], [1, 2, 3], [0, 1, 2]]]
        sage: t.display(e)
        t = e_0⊗e^0⊗e^1 + 2 e_0⊗e^0⊗e^2 - e_0⊗e^1⊗e^0 + e_0⊗e^1⊗e^2
         - 2 e_0⊗e^2⊗e^0 - e_0⊗e^2⊗e^1 + e_1⊗e^0⊗e^0 + 2 e_1⊗e^0⊗e^1
         + 3 e_1⊗e^0⊗e^2 + e_1⊗e^1⊗e^1 + 2 e_1⊗e^1⊗e^2 - e_1⊗e^2⊗e^0
         + e_1⊗e^2⊗e^2 + 2 e_2⊗e^0⊗e^0 + 3 e_2⊗e^0⊗e^1 + 4 e_2⊗e^0⊗e^2
         + e_2⊗e^1⊗e^0 + 2 e_2⊗e^1⊗e^1 + 3 e_2⊗e^1⊗e^2 + e_2⊗e^2⊗e^1
         + 2 e_2⊗e^2⊗e^2

    An alternative is to construct the tensor from an empty list of components
    and to set the nonzero components afterwards::

        sage: t = T([], name='t')
        sage: t.set_comp(e)[0,1,1] = -3
        sage: t.set_comp(e)[2,0,1] = 4
        sage: t.display(e)
        t = -3 e_0⊗e^1⊗e^1 + 4 e_2⊗e^0⊗e^1

    See the documentation of
    :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
    for the full list of arguments that can be provided to the __call__
    operator. For instance, to construct a tensor symmetric with respect to the
    last two indices::

        sage: t = T([], name='t', sym=(1,2))
        sage: t.set_comp(e)[0,1,1] = -3
        sage: t.set_comp(e)[2,0,1] = 4
        sage: t.display(e)  # notice that t^2_{10} has be set equal to t^2_{01} by symmetry
        t = -3 e_0⊗e^1⊗e^1 + 4 e_2⊗e^0⊗e^1 + 4 e_2⊗e^1⊗e^0

    The tensor modules over a given module `M` are unique::

        sage: T is M.tensor_module(1,2)
        True

    There is a coercion map from `\Lambda^p(M^*)`, the set of alternating
    forms of degree `p`, to `T^{(0,p)}(M)`::

        sage: L2 = M.dual_exterior_power(2) ; L2
        2nd exterior power of the dual of the Rank-3 free module M over the
         Integer Ring
        sage: T02 = M.tensor_module(0,2) ; T02
        Free module of type-(0,2) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: T02.has_coerce_map_from(L2)
        True

    Of course, for `p\geq 2`, there is no coercion in the reverse direction,
    since not every tensor of type `(0,p)` is alternating::

        sage: L2.has_coerce_map_from(T02)
        False

    The coercion map `\Lambda^2(M^*)\rightarrow T^{(0,2)}(M)` in action::

        sage: a = M.alternating_form(2, name='a') ; a
        Alternating form a of degree 2 on the Rank-3 free module M over the
         Integer Ring
        sage: a[0,1], a[1,2] = 4, -3
        sage: a.display(e)
        a = 4 e^0∧e^1 - 3 e^1∧e^2
        sage: a.parent() is L2
        True
        sage: ta = T02(a) ; ta
        Type-(0,2) tensor a on the Rank-3 free module M over the Integer Ring
        sage: ta.display(e)
        a = 4 e^0⊗e^1 - 4 e^1⊗e^0 - 3 e^1⊗e^2 + 3 e^2⊗e^1
        sage: ta.symmetries() # the antisymmetry is of course preserved
        no symmetry;  antisymmetry: (0, 1)

    For the degree `p=1`, there is a coercion in both directions::

        sage: L1 = M.dual_exterior_power(1) ; L1
        Dual of the Rank-3 free module M over the Integer Ring
        sage: T01 = M.tensor_module(0,1) ; T01
        Free module of type-(0,1) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: T01.has_coerce_map_from(L1)
        True
        sage: L1.has_coerce_map_from(T01)
        True

    The coercion map `\Lambda^1(M^*)\rightarrow T^{(0,1)}(M)` in action::

        sage: a = M.linear_form('a')
        sage: a[:] = -2, 4, 1 ; a.display(e)
        a = -2 e^0 + 4 e^1 + e^2
        sage: a.parent() is L1
        True
        sage: ta = T01(a) ; ta
        Type-(0,1) tensor a on the Rank-3 free module M over the Integer Ring
        sage: ta.display(e)
        a = -2 e^0 + 4 e^1 + e^2

    The coercion map `T^{(0,1)}(M) \rightarrow \Lambda^1(M^*)` in action::

        sage: ta.parent() is T01
        True
        sage: lta = L1(ta) ; lta
        Linear form a on the Rank-3 free module M over the Integer Ring
        sage: lta.display(e)
        a = -2 e^0 + 4 e^1 + e^2
        sage: lta == a
        True

    There is a canonical identification between tensors of type `(1,1)` and
    endomorphisms of module `M`. Accordingly, coercion maps have been
    implemented between `T^{(1,1)}(M)` and `\mathrm{End}(M)` (the module of
    all endomorphisms of `M`, see
    :class:`~sage.tensor.modules.free_module_homset.FreeModuleHomset`)::

        sage: T11 = M.tensor_module(1,1) ; T11
        Free module of type-(1,1) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: End(M)
        Set of Morphisms from Rank-3 free module M over the Integer Ring
         to Rank-3 free module M over the Integer Ring
         in Category of finite dimensional modules over Integer Ring
        sage: T11.has_coerce_map_from(End(M))
        True
        sage: End(M).has_coerce_map_from(T11)
        True

    The coercion map `\mathrm{End}(M)\rightarrow T^{(1,1)}(M)` in action::

        sage: phi = End(M).an_element() ; phi
        Generic endomorphism of Rank-3 free module M over the Integer Ring
        sage: phi.matrix(e)
        [1 1 1]
        [1 1 1]
        [1 1 1]
        sage: tphi = T11(phi) ; tphi # image of phi by the coercion map
        Type-(1,1) tensor on the Rank-3 free module M over the Integer Ring
        sage: tphi[:]
        [1 1 1]
        [1 1 1]
        [1 1 1]
        sage: t = M.tensor((1,1))
        sage: t[0,0], t[1,1], t[2,2] = -1,-2,-3
        sage: t[:]
        [-1  0  0]
        [ 0 -2  0]
        [ 0  0 -3]
        sage: s = t + phi ; s  # phi is coerced to a type-(1,1) tensor prior to the addition
        Type-(1,1) tensor on the Rank-3 free module M over the Integer Ring
        sage: s[:]
        [ 0  1  1]
        [ 1 -1  1]
        [ 1  1 -2]

    The coercion map `T^{(1,1)}(M) \rightarrow \mathrm{End}(M)` in action::

        sage: phi1 = End(M)(tphi) ; phi1
        Generic endomorphism of Rank-3 free module M over the Integer Ring
        sage: phi1 == phi
        True
        sage: s = phi + t ; s  # t is coerced to an endomorphism prior to the addition
        Generic endomorphism of Rank-3 free module M over the Integer Ring
        sage: s.matrix(e)
        [ 0  1  1]
        [ 1 -1  1]
        [ 1  1 -2]

    There is a coercion `\mathrm{GL}(M)\rightarrow T^{(1,1)}(M)`, i.e. from
    automorphisms of `M` to type-`(1,1)` tensors on `M`::

        sage: GL = M.general_linear_group() ; GL
        General linear group of the Rank-3 free module M over the Integer Ring
        sage: T11.has_coerce_map_from(GL)
        True

    The coercion map `\mathrm{GL}(M)\rightarrow T^{(1,1)}(M)` in action::

        sage: a = GL.an_element() ; a
        Automorphism of the Rank-3 free module M over the Integer Ring
        sage: a.matrix(e)
        [ 1  0  0]
        [ 0 -1  0]
        [ 0  0  1]
        sage: ta = T11(a) ; ta
        Type-(1,1) tensor on the Rank-3 free module M over the Integer Ring
        sage: ta.display(e)
        e_0⊗e^0 - e_1⊗e^1 + e_2⊗e^2
        sage: a.display(e)
        e_0⊗e^0 - e_1⊗e^1 + e_2⊗e^2

    Of course, there is no coercion in the reverse direction, since not
    every type-`(1,1)` tensor is invertible::

        sage: GL.has_coerce_map_from(T11)
        False

    """

    Element = FreeModuleTensor

    def __init__(self, fmodule, tensor_type, name=None, latex_name=None):
        r"""
        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: T = M.tensor_module(2, 3)
            sage: TestSuite(T).run()

        """
        self._fmodule = fmodule
        self._tensor_type = tuple(tensor_type)
        rank = pow(fmodule._rank, tensor_type[0] + tensor_type[1])
        if self._tensor_type == (0,1):  # case of the dual
            if name is None and fmodule._name is not None:
                name = fmodule._name + '*'
            if latex_name is None and fmodule._latex_name is not None:
                latex_name = fmodule._latex_name + r'^*'
        else:
            if name is None and fmodule._name is not None:
                name = 'T^' + str(self._tensor_type) + '(' + fmodule._name + \
                       ')'
            if latex_name is None and fmodule._latex_name is not None:
                latex_name = r'T^{' + str(self._tensor_type) + r'}\left(' + \
                             fmodule._latex_name + r'\right)'
        FiniteRankFreeModule.__init__(self, fmodule._ring, rank, name=name,
                                      latex_name=latex_name,
                                      start_index=fmodule._sindex,
                                      output_formatter=fmodule._output_formatter)
        fmodule._all_modules.add(self)

    #### Parent Methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: T = M.tensor_module(1,1)
            sage: T._element_constructor_(0) is T.zero()
            True
            sage: e = M.basis('e')
            sage: t = T._element_constructor_(comp=[[2,0],[1/2,-3]], basis=e,
            ....:                             name='t') ; t
            Type-(1,1) tensor t on the 2-dimensional vector space M over the
             Rational Field
            sage: t.display()
            t = 2 e_0⊗e^0 + 1/2 e_1⊗e^0 - 3 e_1⊗e^1
            sage: t.parent()
            Free module of type-(1,1) tensors on the 2-dimensional vector
             space M over the Rational Field
            sage: t.parent() is T
            True

        """
        from sage.rings.integer import Integer
        if isinstance(comp, (int, Integer)) and comp == 0:
            return self.zero()
        if isinstance(comp, FiniteRankFreeModuleMorphism):
            # coercion of an endomorphism to a type-(1,1) tensor:
            endo = comp  # for readability
            if self._tensor_type == (1,1) and endo.is_endomorphism() and \
                                                self._fmodule is endo.domain():
                resu = self.element_class(self._fmodule, (1,1),
                                          name=endo._name,
                                          latex_name=endo._latex_name)
                for basis, mat in endo._matrices.items():
                    resu.add_comp(basis[0])[:] = mat
            else:
                raise TypeError("cannot coerce the {}".format(endo) +
                                " to an element of {}".format(self))
        elif isinstance(comp, AlternatingContrTensor):
            # coercion of an alternating contravariant tensor of degree
            # p to a type-(p,0) tensor:
            tensor = comp # for readability
            p = tensor.degree()
            if self._tensor_type != (p,0) or \
                                    self._fmodule != tensor.base_module():
                raise TypeError("cannot coerce the {}".format(tensor) +
                                " to an element of {}".format(self))
            if p == 1:
                asym = None
            else:
                asym = range(p)
            resu = self.element_class(self._fmodule, (p,0),
                                      name=tensor._name,
                                      latex_name=tensor._latex_name,
                                      antisym=asym)
            for basis, comp in tensor._components.items():
                resu._components[basis] = comp.copy()
        elif isinstance(comp, FreeModuleAltForm):
            # coercion of an alternating form to a type-(0,p) tensor:
            form = comp # for readability
            p = form.degree()
            if self._tensor_type != (0,p) or \
                                           self._fmodule != form.base_module():
                raise TypeError("cannot coerce the {}".format(form) +
                                " to an element of {}".format(self))
            if p == 1:
                asym = None
            else:
                asym = range(p)
            resu = self.element_class(self._fmodule, (0,p), name=form._name,
                                      latex_name=form._latex_name,
                                      antisym=asym)
            for basis, comp in form._components.items():
                resu._components[basis] = comp.copy()
        elif isinstance(comp, FreeModuleAutomorphism):
            # coercion of an automorphism to a type-(1,1) tensor:
            autom = comp # for readability
            if self._tensor_type != (1,1) or \
                                          self._fmodule != autom.base_module():
                raise TypeError("cannot coerce the {}".format(autom) +
                                " to an element of {}".format(self))
            resu = self.element_class(self._fmodule, (1,1), name=autom._name,
                                      latex_name=autom._latex_name)
            for basis, comp in autom._components.items():
                resu._components[basis] = comp.copy()
        else:
            # Standard construction:
            resu = self.element_class(self._fmodule, self._tensor_type,
                                      name=name, latex_name=latex_name,
                                      sym=sym, antisym=antisym)
            if comp:
                resu.set_comp(basis)[:] = comp
        return resu

    @cached_method
    def zero(self):
        r"""
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: T11 = M.tensor_module(1,1)
            sage: T11.zero()
            Type-(1,1) tensor zero on the Rank-3 free module M over the Integer
             Ring

        The zero element is cached::

            sage: T11.zero() is T11(0)
            True

        """
        resu = self._element_constructor_(name='zero', latex_name='0')
        for basis in self._fmodule._known_bases:
            resu._add_comp_unsafe(basis)
            # (since new components are initialized to zero)
        resu._is_zero = True # This element is certainly zero
        resu.set_immutable()
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) element of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: T = M.tensor_module(1,1)
            sage: t = T._an_element_() ; t
            Type-(1,1) tensor on the 2-dimensional vector space M over the
             Rational Field
            sage: t.display()
            1/2 e_0⊗e^0
            sage: t.parent() is T
            True
            sage: M.tensor_module(2,3)._an_element_().display()
            1/2 e_0⊗e_0⊗e^0⊗e^0⊗e^0

        """
        resu = self.element_class(self._fmodule, self._tensor_type)
        # Make sure that the base module has a default basis
        self._fmodule.an_element()
        sindex = self._fmodule._sindex
        ind = [sindex for i in range(resu._tensor_rank)]
        resu.set_comp()[ind] = self._fmodule._ring.an_element()
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        EXAMPLES:

        Sets of module endomorphisms coerce to type-`(1,1)` tensor modules::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M.tensor_module(1,1)._coerce_map_from_(End(M))
            True

        but not to tensor modules of other types::

            sage: M.tensor_module(0,1)._coerce_map_from_(End(M))
            False

        and not to type-`(1,1)` tensor modules defined on another free module::

            sage: N = FiniteRankFreeModule(ZZ, 3, name='N')
            sage: f = N.basis('f')
            sage: M.tensor_module(1,1)._coerce_map_from_(End(N))
            False

        There is no coercion if the module morphisms are not endomorphisms::

            sage: M.tensor_module(1,1)._coerce_map_from_(Hom(M,N))
            False

        Coercion from alternating contravariant tensors::

            sage: M.tensor_module(2,0)._coerce_map_from_(M.exterior_power(2))
            True
            sage: M.tensor_module(2,0)._coerce_map_from_(M.exterior_power(3))
            False
            sage: M.tensor_module(2,0)._coerce_map_from_(N.exterior_power(2))
            False

        Coercion from alternating forms::

            sage: M.tensor_module(0,1)._coerce_map_from_(M.dual_exterior_power(1))
            True
            sage: M.tensor_module(0,2)._coerce_map_from_(M.dual_exterior_power(2))
            True
            sage: M.tensor_module(0,2)._coerce_map_from_(M.dual_exterior_power(3))
            False
            sage: M.tensor_module(0,2)._coerce_map_from_(N.dual_exterior_power(2))
            False

        """
        from .free_module_homset import FreeModuleHomset
        from .ext_pow_free_module import (ExtPowerFreeModule,
                                          ExtPowerDualFreeModule)
        from .free_module_linear_group import FreeModuleLinearGroup
        if isinstance(other, FreeModuleHomset):
            # Coercion of an endomorphism to a type-(1,1) tensor:
            if self._tensor_type == (1,1):
                return other.is_endomorphism_set() and \
                                         self._fmodule is other.domain()
            else:
                return False
        if isinstance(other, ExtPowerFreeModule):
            # Coercion of an alternating contravariant tensor to a
            # type-(p,0) tensor:
            return self._tensor_type == (other.degree(), 0) and \
                                    self._fmodule is other.base_module()
        if isinstance(other, ExtPowerDualFreeModule):
            # Coercion of an alternating form to a type-(0,p) tensor:
            return self._tensor_type == (0, other.degree()) and \
                                    self._fmodule is other.base_module()
        if isinstance(other, FreeModuleLinearGroup):
            # Coercion of an automorphism to a type-(1,1) tensor:
            return self._tensor_type == (1,1) and \
                                    self._fmodule is other.base_module()
        return False

    #### End of parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: M.tensor_module(1,1)
            Free module of type-(1,1) tensors on the 2-dimensional vector space
             M over the Rational Field
            sage: M.tensor_module(0,1)
            Free module of type-(0,1) tensors on the 2-dimensional vector space
             M over the Rational Field

        """
        description = "Free module of type-({},{}) tensors on the {}".format(
                     self._tensor_type[0], self._tensor_type[1], self._fmodule)
        return description

    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module on which the tensor module is defined.

        EXAMPLES:

        Base module of a type-`(1,2)` tensor module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: T = M.tensor_module(1,2)
            sage: T.base_module()
            Rank-3 free module M over the Integer Ring
            sage: T.base_module() is M
            True

        """
        return self._fmodule

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``.

        OUTPUT:

        - pair `(k,l)` such that ``self`` is the module tensor product
          `T^{(k,l)}(M)`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: T = M.tensor_module(1,2)
            sage: T.tensor_type()
            (1, 2)

        """
        return self._tensor_type
