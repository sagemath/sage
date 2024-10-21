# sage_setup: distribution = sagemath-modules
r"""
Free submodules of tensor modules defined by monoterm symmetries

AUTHORS:

- Matthias Koeppe (2020-2022): initial version
"""

# ******************************************************************************
#       Copyright (C) 2020-2022 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ******************************************************************************

import itertools

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_set import DisjointSet
from sage.typeset.unicode_characters import unicode_otimes

from .comp import CompFullySym, CompFullyAntiSym, CompWithSym
from .tensor_free_module import TensorFreeModule
from .finite_rank_free_module import FiniteRankFreeModule_abstract


class TensorFreeSubmodule_sym(TensorFreeModule):
    r"""
    Class for free submodules of tensor products of free modules
    that are defined by some monoterm symmetries.

    EXAMPLES::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: T60M = M.tensor_module(6, 0); T60M
        Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring
        sage: T60M._name
        'T^(6, 0)(M)'
        sage: latex(T60M)
        T^{(6, 0)}\left(M\right)
        sage: T40Sym45M = M.tensor_module(6, 0, sym=((4, 5))); T40Sym45M
        Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
         with symmetry on the index positions (4, 5)
        sage: T40Sym45M._name
        'T^{0,1,2,3}(M)⊗Sym^{4,5}(M)'
        sage: latex(T40Sym45M)
        T^{\{0,1,2,3\}}(M) \otimes \mathrm{Sym}^{\{4,5\}}(M)
        sage: Sym0123x45M = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5))); Sym0123x45M
        Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
         with symmetry on the index positions (0, 1, 2, 3),
         with symmetry on the index positions (4, 5)
        sage: Sym0123x45M._name
        'Sym^{0,1,2,3}(M)⊗Sym^{4,5}(M)'
        sage: latex(Sym0123x45M)
        \mathrm{Sym}^{\{0,1,2,3\}}(M) \otimes \mathrm{Sym}^{\{4,5\}}(M)
        sage: Sym012x345M = M.tensor_module(6, 0, sym=((0, 1, 2), (3, 4, 5))); Sym012x345M
        Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
         with symmetry on the index positions (0, 1, 2),
         with symmetry on the index positions (3, 4, 5)
        sage: Sym012x345M._name
        'Sym^{0,1,2}(M)⊗Sym^{3,4,5}(M)'
        sage: latex(Sym012x345M)
        \mathrm{Sym}^{\{0,1,2\}}(M) \otimes \mathrm{Sym}^{\{3,4,5\}}(M)
        sage: Sym012345M = M.tensor_module(6, 0, sym=((0, 1, 2, 3, 4, 5))); Sym012345M
        Free module of fully symmetric type-(6,0) tensors
         on the Rank-3 free module M over the Integer Ring
        sage: Sym012345M._name
        'Sym^6(M)'
        sage: latex(Sym012345M)
        \mathrm{Sym}^6(M)

    Canonical injections from submodules are coercions::

        sage: Sym0123x45M.has_coerce_map_from(Sym012345M)
        True
        sage: T60M.has_coerce_map_from(Sym0123x45M)
        True
        sage: t = e[0] * e[0] * e[0] * e[0] * e[0] * e[0]
        sage: t.parent()
        Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring
        sage: Sym012345M(t) is t
        False

    TESTS::

        sage: T = M.tensor_module(4, 4, sym=((0, 1)), antisym=((4, 5))); T
        Free module of type-(4,4) tensors on the Rank-3 free module M over the Integer Ring,
         with symmetry on the index positions (0, 1),
         with antisymmetry on the index positions (4, 5)
        sage: T._name
        'T^{2,3}(M)⊗T^{6,7}(M*)⊗Sym^{0,1}(M)⊗ASym^{4,5}(M*)'
        sage: latex(T)
        T^{\{2,3\}}(M) \otimes T^{\{6,7\}}(M^*) \otimes \mathrm{Sym}^{\{0,1\}}(M) \otimes \mathrm{ASym}^{\{4,5\}}(M^*)
    """
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None,
                 sym=None, antisym=None, *, category=None, ambient=None):
        r"""
        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: Sym0123x45M = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: TestSuite(Sym0123x45M).run()
        """
        self._fmodule = fmodule
        self._tensor_type = tuple(tensor_type)
        if ambient is None:
            ambient = fmodule.tensor_module(*tensor_type)
        self._ambient_module = ambient
        self._sym = sym
        self._antisym = antisym
        basis_sym = self._basis_sym()
        rank = len(list(basis_sym.non_redundant_index_generator()))

        if name is None and fmodule._name is not None:
            all_indices = tuple(range(tensor_type[0] + tensor_type[1]))
            if isinstance(basis_sym, CompFullySym):
                sym = [all_indices]
                antisym = []
            elif isinstance(basis_sym, CompFullyAntiSym):
                sym = []
                antisym = [all_indices]
            elif isinstance(basis_sym, CompWithSym):
                sym = basis_sym._sym
                antisym = basis_sym._antisym
            else:
                sym = antisym = []
            nosym_0 = [i for i in range(tensor_type[0])
                       if not any(i in s for s in sym) and not any(i in s for s in antisym)]
            nosym_1 = [i for i in range(tensor_type[0], tensor_type[0] + tensor_type[1])
                       if not any(i in s for s in sym) and not any(i in s for s in antisym)]
            nosym = [s for s in [nosym_0, nosym_1] if s]

            def power_name(op, s, latex=False):
                if s[0] < tensor_type[0]:
                    assert all(i < tensor_type[0] for i in s)
                    base = fmodule
                    full = tensor_type[0]
                else:
                    assert all(i >= tensor_type[0] for i in s)
                    base = fmodule.dual()
                    full = tensor_type[1]
                if len(s) == full:
                    superscript = str(full)
                else:
                    superscript = ','.join(str(i) for i in s)
                    if latex:
                        superscript = r'\{' + superscript + r'\}'
                    else:
                        superscript = '{' + superscript + '}'
                if latex:
                    if len(superscript) != 1:
                        superscript = '{' + superscript + '}'
                    if len(base._latex_name) > 3:
                        return op + '^' + superscript + r'\left(' + base._latex_name + r'\right)'
                    else:
                        return op + '^' + superscript + '(' + base._latex_name + ')'
                else:
                    return op + '^' + superscript + '(' + base._name + ')'

            name = unicode_otimes.join(itertools.chain(
                (power_name('T', s, latex=False) for s in nosym),
                (power_name('Sym', s, latex=False) for s in sym),
                (power_name('ASym', s, latex=False) for s in antisym)))
            latex_name = r' \otimes '.join(itertools.chain(
                (power_name('T', s, latex=True) for s in nosym),
                (power_name(r'\mathrm{Sym}', s, latex=True) for s in sym),
                (power_name(r'\mathrm{ASym}', s, latex=True) for s in antisym)))

        category = fmodule.category().TensorProducts().FiniteDimensional().Subobjects().or_subcategory(category)
        # Skip TensorFreeModule.__init__
        FiniteRankFreeModule_abstract.__init__(self, fmodule._ring, rank, name=name,
                                               latex_name=latex_name,
                                               category=category, ambient=ambient)

    def construction(self):
        # TODO: Define the symmetry group and its action (https://github.com/sagemath/sage/issues/34495),
        # return the construction functor for invariant subobjects.
        r"""
        Return the functorial construction of ``self``.

        This implementation just returns ``None``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym2M = M.tensor_module(2, 0, sym=range(2)); Sym2M
            Free module of fully symmetric type-(2,0) tensors on the Rank-3 free module M over the Integer Ring
            sage: Sym2M.construction() is None
            True
        """
        return None

    @cached_method
    def _basis_sym(self):
        r"""
        Return an instance of :class:`~sage.tensor.modules.comp.Components`.

        In the current implementation of :class:`~sage.tensor.modules.tensor_free_submodule.TensorFreeSubmodule_sym`,
        it encodes the prescribed symmetry of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym2M = M.tensor_module(2, 0, sym=range(2)); Sym2M
            Free module of fully symmetric type-(2,0) tensors on the Rank-3 free module M over the Integer Ring
            sage: c = Sym2M._basis_sym(); c
            Fully symmetric 2-indices components w.r.t. (0, 1, 2)
        """
        frame = tuple(self.base_module().irange())
        # Need to call _element_constructor_ explicitly, or the passed arguments are dropped
        tensor = self.ambient()._element_constructor_(sym=self._sym, antisym=self._antisym)
        return tensor._new_comp(frame)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym2M = M.tensor_module(2, 0, sym=range(2)); Sym2M
            Free module of fully symmetric type-(2,0) tensors
             on the Rank-3 free module M over the Integer Ring
        """
        prefix, suffix = self._basis_sym()._repr_symmetry()
        return "Free module of {}type-({},{}) tensors on the {}{}".format(
            prefix.lower(), self._tensor_type[0], self._tensor_type[1], self._fmodule, suffix)

    def _is_symmetry_coarsening_of(self, coarser_comp, finer_comp):
        r"""
        Return whether ``coarser_comp`` has coarser symmetry than ``finer_comp``.

        INPUT:

        - ``coarser_comp``, ``finer_comp`` -- :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: T60M = M.tensor_module(6, 0)
            sage: Sym0123x45M = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: ten0123x45M = Sym0123x45M.an_element(); ten0123x45M
            Type-(6,0) tensor on the Rank-3 free module M over the Integer Ring
            sage: ten0123x45M.parent()
            Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
             with symmetry on the index positions (0, 1, 2, 3),
             with symmetry on the index positions (4, 5)
            sage: com0123x45M = ten0123x45M._components[e]; com0123x45M
            6-indices components w.r.t. Basis (e_0,e_1,e_2)
             on the Rank-3 free module M over the Integer Ring,
             with symmetry on the index positions (0, 1, 2, 3),
             with symmetry on the index positions (4, 5)
            sage: Sym012x345M = M.tensor_module(6, 0, sym=((0, 1, 2), (3, 4, 5)))
            sage: com012x345M = Sym012x345M.an_element()._components[e]; com012x345M
            6-indices components w.r.t. Basis (e_0,e_1,e_2)
             on the Rank-3 free module M over the Integer Ring,
             with symmetry on the index positions (0, 1, 2),
             with symmetry on the index positions (3, 4, 5)
            sage: Sym012345M  = M.tensor_module(6, 0, sym=((0, 1, 2, 3, 4, 5)))
            sage: com012345M  = Sym012345M.an_element()._components[e]; com012345M
            Fully symmetric 6-indices components w.r.t. Basis (e_0,e_1,e_2)
             on the Rank-3 free module M over the Integer Ring
            sage: Sym0123x45M._is_symmetry_coarsening_of(com0123x45M, com012x345M)
            False
            sage: Sym0123x45M._is_symmetry_coarsening_of(com012345M, com012x345M)
            True
        """
        self_tensor_type = self.tensor_type()

        def sym_antisym(comp):
            if isinstance(comp, tuple):
                sym, antisym = tuple
                if sym is None:
                    sym = []
                if antisym is None:
                    antisym = []
                return sym, antisym
            # Similar code is in Component.contract, should refactor.
            try:
                return comp._sym, comp._antisym
            except AttributeError:
                return [], []

        def is_coarsening_of(self_sym_list, other_sym_list):
            # Use the union-find data structure
            S = DisjointSet(self_tensor_type[0] + self_tensor_type[1])
            for index_set in self_sym_list:
                i = index_set[0]
                for j in index_set[1:]:
                    S.union(i, j)
            for index_set in other_sym_list:
                i = S.find(index_set[0])
                for j in index_set[1:]:
                    if S.find(j) != i:
                        return False
            return True

        finer_sym, finer_antisym = sym_antisym(finer_comp)
        if not finer_sym and not finer_antisym:
            return True
        coarser_sym, coarser_antisym = sym_antisym(coarser_comp)
        if not is_coarsening_of(coarser_sym, finer_sym):
            return False
        if not is_coarsening_of(coarser_antisym, finer_antisym):
            return False
        return True

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None, sym=None, antisym=None):
        r"""
        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: T60M = M.tensor_module(6, 0)
            sage: Sym0123x45M = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: Sym0123x45M(e[0]*e[0]*e[0]*e[0]*e[1]*e[2])
            Traceback (most recent call last):
            ...
            ValueError: this tensor does not have the symmetries of
             Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
              with symmetry on the index positions (0, 1, 2, 3),
              with symmetry on the index positions (4, 5)
            sage: t = Sym0123x45M(e[0]*e[0]*e[0]*e[0]*e[1]*e[2] + e[0]*e[0]*e[0]*e[0]*e[2]*e[1]); t.disp()
            e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2 + e_0⊗e_0⊗e_0⊗e_0⊗e_2⊗e_1
            sage: t.parent()._name
            'Sym^{0,1,2,3}(M)⊗Sym^{4,5}(M)'
        """
        if sym is not None or antisym is not None:
            # Refuse to create a tensor with finer symmetries
            # than those defining the subspace
            if not self._is_symmetry_coarsening_of((sym, antisym), self._basis_sym()):
                raise ValueError(f"cannot create a tensor with symmetries {sym=}, {antisym=} "
                                 f"as an element of {self}")

        if sym is None:
            sym = self._basis_sym()._sym
        if antisym is None:
            antisym = self._basis_sym()._antisym

        resu = super()._element_constructor_(comp=comp,
                                             basis=basis, name=name,
                                             latex_name=latex_name,
                                             sym=sym, antisym=antisym)
        if not resu._components:
            # fast path for zero tensor
            return resu

        try:
            if self.reduce(resu):
                raise ValueError(f"this tensor does not have the symmetries of {self}")
        except TypeError:
            # Averaging over the orbits of a tensor that does not have the required
            # symmetries can lead to "TypeError: no conversion of this rational to integer"
            raise ValueError(f"this tensor does not have the symmetries of {self}")

        return resu

    def is_submodule(self, other):
        r"""
        Return ``True`` if ``self`` is a submodule of ``other``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: T60M = M.tensor_module(6, 0)
            sage: Sym0123x45M = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: Sym012x345M = M.tensor_module(6, 0, sym=((0, 1, 2), (3, 4, 5)))
            sage: Sym012345M  = M.tensor_module(6, 0, sym=((0, 1, 2, 3, 4, 5)))
            sage: Sym012345M.is_submodule(Sym012345M)
            True
            sage: Sym012345M.is_submodule(Sym0123x45M)
            True
            sage: Sym0123x45M.is_submodule(Sym012345M)
            False
            sage: Sym012x345M.is_submodule(Sym0123x45M)
            False
            sage: all(S.is_submodule(T60M) for S in (Sym0123x45M, Sym012x345M, Sym012345M))
            True
        """
        if super().is_submodule(other):
            return True
        self_base_module = self.base_module()
        self_tensor_type = self.tensor_type()
        try:
            other_base_module = other.base_module()
            other_tensor_type = other.tensor_type()
        except AttributeError:
            return False
        if self_base_module != other_base_module:
            return False
        if self_tensor_type != other_tensor_type:
            return False

        other_comp = other._basis_sym()
        return self._is_symmetry_coarsening_of(self._basis_sym(), other_comp)

    @lazy_attribute
    def lift(self):
        r"""
        The lift (embedding) map from ``self`` to the ambient space.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym0123x45M = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: Sym0123x45M.lift
            Generic morphism:
              From: Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
                     with symmetry on the index positions (0, 1, 2, 3),
                     with symmetry on the index positions (4, 5)
              To:   Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring
        """
        return self.module_morphism(function=lambda x: x, codomain=self.ambient())

    @lazy_attribute
    def reduce(self):
        r"""
        The reduce map.

        This map reduces elements of the ambient space modulo this
        submodule.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: X = M.tensor_module(6, 0)
            sage: Y = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: Y.reduce
            Generic endomorphism of
             Free module of type-(6,0) tensors on the 3-dimensional vector space M over the Rational Field
            sage: t = e[0]*e[0]*e[0]*e[0]*e[1]*e[2]; t.disp()
            e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2 = e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2
            sage: r = Y.reduce(t); r
            Type-(6,0) tensor on the 3-dimensional vector space M over the Rational Field
            sage: r.disp()
            1/2 e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2 - 1/2 e_0⊗e_0⊗e_0⊗e_0⊗e_2⊗e_1
            sage: r.parent()._name
            'T^(6, 0)(M)'

        If the base ring is not a field, this may fail::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: X = M.tensor_module(6, 0)
            sage: Y = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: Y.reduce
            Generic endomorphism of
             Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring
            sage: t = e[0]*e[0]*e[0]*e[0]*e[1]*e[2]; t.disp()
            e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2 = e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2
            sage: Y.reduce(t)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: X = M.tensor_module(6, 0)
            sage: Y = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: all(Y.reduce(u.lift()) == 0 for u in Y.basis('e'))
            True
        """
        sym = self._basis_sym()._sym
        antisym = self._basis_sym()._antisym

        def _reduce_element(x):
            if not x._components:
                # zero tensor - methods symmetrize, antisymmetrize are broken
                return x
            # TODO: Implement a fast symmetry check, either as part of the symmetrize/antisymmetrize methods,
            #       or as a separate method
            symmetrized = x
            for s in sym:
                symmetrized = symmetrized.symmetrize(*s)
            for s in antisym:
                symmetrized = symmetrized.antisymmetrize(*s)
            return x - symmetrized

        return self.ambient().module_morphism(function=_reduce_element, codomain=self.ambient())

    @lazy_attribute
    def retract(self):
        r"""
        The retract map from the ambient space.

        This is a partial map, which gives an error for elements not in the subspace.

        Calling this map on elements of the ambient space is the same as calling the
        element constructor of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: X = M.tensor_module(6, 0)
            sage: Y = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: e_Y = Y.basis('e')
            sage: Y.retract
            Generic morphism:
              From: Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring
              To:   Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
                     with symmetry on the index positions (0, 1, 2, 3),
                     with symmetry on the index positions (4, 5)

            sage: t = e[0]*e[0]*e[0]*e[0]*e[1]*e[2]; t.disp()
            e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2 = e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2
            sage: Y.retract(t)
            Traceback (most recent call last):
            ...
            ValueError: this tensor does not have the symmetries of
             Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring,
              with symmetry on the index positions (0, 1, 2, 3),
              with symmetry on the index positions (4, 5)
            sage: t = e[0]*e[0]*e[0]*e[0]*e[1]*e[2] + e[0]*e[0]*e[0]*e[0]*e[2]*e[1]
            sage: y = Y.retract(t); y
            Type-(6,0) tensor on the Rank-3 free module M over the Integer Ring
            sage: y.disp()
            e_0⊗e_0⊗e_0⊗e_0⊗e_1⊗e_2 + e_0⊗e_0⊗e_0⊗e_0⊗e_2⊗e_1
            sage: y.parent()._name
            'Sym^{0,1,2,3}(M)⊗Sym^{4,5}(M)'

        TESTS::

            sage: all(Y.retract(u.lift()) == u for u in e_Y)
            True
        """
        return self.ambient().module_morphism(function=lambda x: self(x), codomain=self)
