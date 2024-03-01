r"""
Base classes for reflexive modules
"""

from sage.misc.abstract_method import abstract_method
from sage.structure.parent import Parent


class ReflexiveModule_abstract(Parent):
    r"""
    Abstract base class for reflexive modules.

    An `R`-module `M` is *reflexive* if the natural map from `M` to its double
    dual `M^{**}` is an isomorphism.

    In the category of `R`-modules, the dual module `M^*` is
    the `R`-module of linear functionals `\phi:\ M \longrightarrow R`.
    However, we do not make the assumption that the dual module
    (obtained by :meth:`dual`) is in the category :class:`Homsets`.

    We identify the double dual `M^{**}` with `M`.

    Tensor products of reflexive modules are reflexive. We identify all
    tensor products of `k` copies of `M` and `l` copies of `M^*` and
    denote it by `T^{(k,l)}(M)`. The :meth:`tensor_type` of such a tensor
    product is the pair `(k, l)`, and `M` is called its :meth:`base_module`.

    There are three abstract subclasses:

    - :class:`ReflexiveModule_base` is the base class for implementations
      of base modules `M`.

    - :class:`ReflexiveModule_dual` is the base class for implementations
      of duals `M^*`.

    - :class:`ReflexiveModule_tensor` is the base class for implementations
      of tensor modules `T^{(k,l)}(M)`.

    TESTS::

        sage: from sage.tensor.modules.reflexive_module import (
        ....:     ReflexiveModule_abstract, ReflexiveModule_base,
        ....:     ReflexiveModule_dual, ReflexiveModule_tensor)
        sage: M = FiniteRankFreeModule(ZZ, 3)
        sage: isinstance(M, ReflexiveModule_abstract)
        True
        sage: isinstance(M, ReflexiveModule_base)
        True
        sage: isinstance(M.dual(), ReflexiveModule_abstract)
        True
        sage: isinstance(M.dual(), ReflexiveModule_dual)
        True
        sage: isinstance(M.tensor_module(1, 1), ReflexiveModule_abstract)
        True
        sage: isinstance(M.tensor_module(1, 1), ReflexiveModule_tensor)
        True
    """

    @abstract_method(optional=True)
    def tensor_type(self):
        r"""
        Return the tensor type of ``self``.

        OUTPUT:

        - pair `(k,l)` such that ``self`` is the module tensor product
          `T^{(k,l)}(M)`, where `M` is the :meth:`base_module` of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: T = M.tensor_module(1, 2)
            sage: T.tensor_type()
            (1, 2)
        """

    @abstract_method
    def base_module(self):
        r"""
        Return the module on which ``self`` is constructed.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: M.base_module() is M
            True
            sage: M.dual().base_module() is M
            True
            sage: M.tensor_module(1, 2).base_module() is M
            True
        """

    def dual(self):
        r"""
        Return the dual module.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: M.dual()
            Dual of the Rank-3 free module over the Integer Ring
            sage: M.dual().dual()
            Rank-3 free module over the Integer Ring
            sage: M.tensor_module(1, 2)
            Free module of type-(1,2) tensors on the Rank-3 free module over the Integer Ring
            sage: M.tensor_module(1, 2).dual()
            Free module of type-(2,1) tensors on the Rank-3 free module over the Integer Ring
        """
        k, l = self.tensor_type()
        return self.base_module().tensor_module(l, k)

    def tensor(self, *args, **kwds):
        # Until https://github.com/sagemath/sage/issues/30373 is done,
        # TensorProductFunctor._functor_name is "tensor", so here we delegate.
        r"""
        Return the tensor product of ``self`` and ``others``.

        This method is invoked when :class:`~sage.categories.tensor.TensorProductFunctor`
        is applied to parents.

        It just delegates to :meth:`tensor_product`.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2); M
            2-dimensional vector space over the Rational Field
            sage: M20 = M.tensor_module(2, 0); M20
            Free module of type-(2,0) tensors on the 2-dimensional vector space over the Rational Field
            sage: tensor([M20, M20])
            Free module of type-(4,0) tensors on the 2-dimensional vector space over the Rational Field
        """
        return self.tensor_product(*args, **kwds)

    def tensor_power(self, n):
        r"""
        Return the ``n``-fold tensor product of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2)
            sage: M.tensor_power(3)
            Free module of type-(3,0) tensors on the 2-dimensional vector space over the Rational Field
            sage: M.tensor_module(1,2).tensor_power(3)
            Free module of type-(3,6) tensors on the 2-dimensional vector space over the Rational Field
        """
        tensor_type = self.tensor_type()
        return self.base_module().tensor_module(n * tensor_type[0], n * tensor_type[1])

    def tensor_product(self, *others):
        r"""
        Return the tensor product of ``self`` and ``others``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2)
            sage: M.tensor_product(M)
            Free module of type-(2,0) tensors on the 2-dimensional vector space over the Rational Field
            sage: M.tensor_product(M.dual())
            Free module of type-(1,1) tensors on the 2-dimensional vector space over the Rational Field
            sage: M.dual().tensor_product(M, M.dual())
            Free module of type-(1,2) tensors on the 2-dimensional vector space over the Rational Field
            sage: M.tensor_product(M.tensor_module(1,2))
            Free module of type-(2,2) tensors on the 2-dimensional vector space over the Rational Field
            sage: M.tensor_module(1,2).tensor_product(M)
            Free module of type-(2,2) tensors on the 2-dimensional vector space over the Rational Field
            sage: M.tensor_module(1,1).tensor_product(M.tensor_module(1,2))
            Free module of type-(2,3) tensors on the 2-dimensional vector space over the Rational Field

            sage: Sym2M = M.tensor_module(2, 0, sym=range(2)); Sym2M
            Free module of fully symmetric type-(2,0) tensors on the 2-dimensional vector space over the Rational Field
            sage: Sym01x23M = Sym2M.tensor_product(Sym2M); Sym01x23M
            Free module of type-(4,0) tensors on the 2-dimensional vector space over the Rational Field,
             with symmetry on the index positions (0, 1), with symmetry on the index positions (2, 3)
            sage: Sym01x23M._index_maps
            ((0, 1), (2, 3))

            sage: N = M.tensor_module(3, 3, sym=[1, 2], antisym=[3, 4]); N
            Free module of type-(3,3) tensors on the 2-dimensional vector space over the Rational Field,
             with symmetry on the index positions (1, 2),
             with antisymmetry on the index positions (3, 4)
            sage: NxN = N.tensor_product(N); NxN
            Free module of type-(6,6) tensors on the 2-dimensional vector space over the Rational Field,
             with symmetry on the index positions (1, 2), with symmetry on the index positions (4, 5),
             with antisymmetry on the index positions (6, 7), with antisymmetry on the index positions (9, 10)
            sage: NxN._index_maps
            ((0, 1, 2, 6, 7, 8), (3, 4, 5, 9, 10, 11))
        """
        from sage.modules.free_module_element import vector
        from .comp import CompFullySym, CompFullyAntiSym, CompWithSym

        base_module = self.base_module()
        if not all(module.base_module() == base_module for module in others):
            raise NotImplementedError('all factors must be tensor modules over the same base module')
        factors = [self] + list(others)
        result_tensor_type = sum(vector(factor.tensor_type()) for factor in factors)
        result_sym = []
        result_antisym = []
        # Keep track of reordering of the contravariant and covariant indices
        # (compatible with FreeModuleTensor.__mul__)
        index_maps = []
        running_indices = vector([0, result_tensor_type[0]])
        for factor in factors:
            tensor_type = factor.tensor_type()
            index_map = tuple(i + running_indices[0] for i in range(tensor_type[0]))
            index_map += tuple(i + running_indices[1] for i in range(tensor_type[1]))
            index_maps.append(index_map)

            if tensor_type[0] + tensor_type[1] > 1:
                basis_sym = factor._basis_sym()
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

                def map_isym(isym):
                    return tuple(index_map[i] for i in isym)

                result_sym.extend(tuple(index_map[i] for i in isym) for isym in sym)
                result_antisym.extend(tuple(index_map[i] for i in isym) for isym in antisym)

            running_indices += vector(tensor_type)

        result = base_module.tensor_module(*result_tensor_type,
                                           sym=result_sym, antisym=result_antisym)
        result._index_maps = tuple(index_maps)
        return result


class ReflexiveModule_base(ReflexiveModule_abstract):
    r"""
    Abstract base class for reflexive modules that are base modules.

    TESTS::

        sage: from sage.tensor.modules.reflexive_module import ReflexiveModule_base
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: isinstance(M, ReflexiveModule_base)
        True
    """

    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed, namely ``self`` itself.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.base_module() is M
            True

            sage: M = Manifold(2, 'M')                                                  # needs sage.symbolic
            sage: XM = M.vector_field_module()                                          # needs sage.symbolic
            sage: XM.base_module() is XM                                                # needs sage.symbolic
            True
        """
        return self

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``, the pair `(1, 0)`.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: M.tensor_type()
            (1, 0)

            sage: M = Manifold(2, 'M')                                                  # needs sage.symbolic
            sage: XM = M.vector_field_module()                                          # needs sage.symbolic
            sage: XM.tensor_type()                                                      # needs sage.symbolic
            (1, 0)
        """
        return (1, 0)

    def dual(self):
        r"""
        Return the dual module.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.dual()
            Dual of the Rank-3 free module M over the Integer Ring
        """
        return self.tensor_module(0, 1)

    @abstract_method
    def tensor_module(self, k, l, **kwds):
        r"""
        Return the module of all tensors of type `(k, l)` defined on ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: M.tensor_module(1, 2)
            Free module of type-(1,2) tensors on the Rank-3 free module over the Integer Ring
        """


class ReflexiveModule_dual(ReflexiveModule_abstract):
    r"""
    Abstract base class for reflexive modules that are the duals of base modules.

    TESTS::

        sage: from sage.tensor.modules.reflexive_module import ReflexiveModule_dual
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: isinstance(M.dual(), ReflexiveModule_dual)
        True
    """

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.dual().tensor_type()
            (0, 1)
        """
        return (0, 1)

    def construction(self):
        r"""
        Return the functorial construction of ``self``.

        This implementation just returns ``None``, as no functorial construction is implemented.

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.dual()
            sage: A.construction() is None
            True
        """
        # Until https://github.com/sagemath/sage/issues/34605 is done
        return None


class ReflexiveModule_tensor(ReflexiveModule_abstract):
    r"""
    Abstract base class for reflexive modules that are tensor products of base modules.

    TESTS::

        sage: from sage.tensor.modules.reflexive_module import ReflexiveModule_tensor
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: isinstance(M.tensor_module(1, 1), ReflexiveModule_tensor)
        True
    """

    def tensor_factors(self):
        r"""
        Return the tensor factors of this tensor module.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: T = M.tensor_module(2, 3)
            sage: T.tensor_factors()
            [Rank-3 free module M over the Integer Ring,
             Rank-3 free module M over the Integer Ring,
             Dual of the Rank-3 free module M over the Integer Ring,
             Dual of the Rank-3 free module M over the Integer Ring,
             Dual of the Rank-3 free module M over the Integer Ring]
        """
        tensor_type = self.tensor_type()
        if tensor_type == (0,1):  # case of the dual
            raise NotImplementedError
        bmodule = self.base_module()
        factors = [bmodule] * tensor_type[0]
        dmodule = bmodule.dual()
        if tensor_type[1]:
            factors += [dmodule] * tensor_type[1]
        return factors
