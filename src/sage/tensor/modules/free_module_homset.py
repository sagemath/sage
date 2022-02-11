r"""
Sets of morphisms between free modules

The class :class:`FreeModuleHomset` implements sets of homomorphisms between
two free modules of finite rank over the same commutative ring.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

REFERENCES:

- Chaps. 13, 14 of R. Godement : *Algebra* [God1968]_
- Chap. 3 of S. Lang : *Algebra* [Lan2002]_

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

from sage.categories.homset import Homset
from sage.tensor.modules.free_module_morphism import FiniteRankFreeModuleMorphism
from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism
from sage.tensor.modules.free_module_tensor import FreeModuleTensor

class FreeModuleHomset(Homset):
    r"""
    Set of homomorphisms between free modules of finite rank over a
    commutative ring.

    Given two free modules `M` and `N` of respective ranks `m` and `n` over a
    commutative ring `R`, the class :class:`FreeModuleHomset` implements the
    set `\mathrm{Hom}(M,N)` of homomorphisms `M\rightarrow N`.
    The set `\mathrm{Hom}(M,N)` is actually a free module of rank `mn` over
    `R`, but this aspect is not taken into account here.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`.

    INPUT:

    - ``fmodule1`` -- free module `M` (domain of the homomorphisms), as an
      instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``fmodule2`` -- free module `N` (codomain of the homomorphisms), as an
      instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``name`` -- (default: ``None``) string; name given to the hom-set; if
      none is provided, Hom(M,N) will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      hom-set; if none is provided, `\mathrm{Hom}(M,N)` will be used

    EXAMPLES:

    Set of homomorphisms between two free modules over `\ZZ`::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
        sage: H = Hom(M,N) ; H
        Set of Morphisms from Rank-3 free module M over the Integer Ring
         to Rank-2 free module N over the Integer Ring
         in Category of finite dimensional modules over Integer Ring
        sage: type(H)
        <class 'sage.tensor.modules.free_module_homset.FreeModuleHomset_with_category_with_equality_by_id'>
        sage: H.category()
        Category of homsets of modules over Integer Ring

    Hom-sets are cached::

        sage: H is Hom(M,N)
        True

    The LaTeX formatting is::

        sage: latex(H)
        \mathrm{Hom}\left(M,N\right)

    As usual, the construction of an element is performed by the ``__call__``
    method; the argument can be the matrix representing the morphism in the
    default bases of the two modules::

        sage: e = M.basis('e')
        sage: f = N.basis('f')
        sage: phi = H([[-1,2,0], [5,1,2]]) ; phi
        Generic morphism:
          From: Rank-3 free module M over the Integer Ring
          To:   Rank-2 free module N over the Integer Ring
        sage: phi.parent() is H
        True

    An example of construction from a matrix w.r.t. bases that are not the
    default ones::

        sage: ep = M.basis('ep', latex_symbol=r"e'")
        sage: fp = N.basis('fp', latex_symbol=r"f'")
        sage: phi2 = H([[3,2,1], [1,2,3]], bases=(ep,fp)) ; phi2
        Generic morphism:
          From: Rank-3 free module M over the Integer Ring
          To:   Rank-2 free module N over the Integer Ring

    The zero element::

        sage: z = H.zero() ; z
        Generic morphism:
          From: Rank-3 free module M over the Integer Ring
          To:   Rank-2 free module N over the Integer Ring
        sage: z.matrix(e,f)
        [0 0 0]
        [0 0 0]

    The test suite for H is passed::

        sage: TestSuite(H).run()

    The set of homomorphisms `M\rightarrow M`, i.e. endomorphisms, is
    obtained by the function ``End``::

        sage: End(M)
        Set of Morphisms from Rank-3 free module M over the Integer Ring
         to Rank-3 free module M over the Integer Ring
         in Category of finite dimensional modules over Integer Ring

    ``End(M)`` is actually identical to ``Hom(M,M)``::

        sage: End(M) is Hom(M,M)
        True

    The unit of the endomorphism ring is the identity map::

        sage: End(M).one()
        Identity endomorphism of Rank-3 free module M over the Integer Ring

    whose matrix in any basis is of course the identity matrix::

        sage: End(M).one().matrix(e)
        [1 0 0]
        [0 1 0]
        [0 0 1]

    There is a canonical identification between endomorphisms of `M` and
    tensors of type `(1,1)` on `M`. Accordingly, coercion maps have been
    implemented between `\mathrm{End}(M)` and `T^{(1,1)}(M)` (the module of
    all type-`(1,1)` tensors on `M`, see
    :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule`)::

        sage: T11 = M.tensor_module(1,1) ; T11
        Free module of type-(1,1) tensors on the Rank-3 free module M over
         the Integer Ring
        sage: End(M).has_coerce_map_from(T11)
        True
        sage: T11.has_coerce_map_from(End(M))
        True

    See :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule` for
    examples of the above coercions.

    There is a coercion `\mathrm{GL}(M) \rightarrow \mathrm{End}(M)`, since
    every automorphism is an endomorphism::

        sage: GL = M.general_linear_group() ; GL
        General linear group of the Rank-3 free module M over the Integer Ring
        sage: End(M).has_coerce_map_from(GL)
        True

    Of course, there is no coercion in the reverse direction, since only
    bijective endomorphisms are automorphisms::

        sage: GL.has_coerce_map_from(End(M))
        False

    The coercion `\mathrm{GL}(M) \rightarrow \mathrm{End}(M)` in action::

        sage: a = GL.an_element() ; a
        Automorphism of the Rank-3 free module M over the Integer Ring
        sage: a.matrix(e)
        [ 1  0  0]
        [ 0 -1  0]
        [ 0  0  1]
        sage: ea = End(M)(a) ; ea
        Generic endomorphism of Rank-3 free module M over the Integer Ring
        sage: ea.matrix(e)
        [ 1  0  0]
        [ 0 -1  0]
        [ 0  0  1]

    """

    Element = FiniteRankFreeModuleMorphism

    def __init__(self, fmodule1, fmodule2, name=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_homset import FreeModuleHomset
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: FreeModuleHomset(M, N)
            Set of Morphisms from Rank-3 free module M over the Integer Ring
             to Rank-2 free module N over the Integer Ring
             in Category of finite dimensional modules over Integer Ring

            sage: H = FreeModuleHomset(M, N, name='L(M,N)',
            ....:                      latex_name=r'\mathcal{L}(M,N)')
            sage: latex(H)
            \mathcal{L}(M,N)

        """
        from .finite_rank_free_module import FiniteRankFreeModule
        if not isinstance(fmodule1, FiniteRankFreeModule):
            raise TypeError("fmodule1 = {} is not an ".format(fmodule1) +
                            "instance of FiniteRankFreeModule")
        if not isinstance(fmodule2, FiniteRankFreeModule):
            raise TypeError("fmodule2 = {} is not an ".format(fmodule2) +
                            "instance of FiniteRankFreeModule")
        if fmodule1.base_ring() != fmodule2.base_ring():
            raise TypeError("the domain and codomain are not defined over " +
                            "the same ring")
        Homset.__init__(self, fmodule1, fmodule2)
        if name is None:
            self._name = "Hom(" + fmodule1._name + "," + fmodule2._name + ")"
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = \
                    r"\mathrm{Hom}\left(" + fmodule1._latex_name + "," + \
                    fmodule2._latex_name + r"\right)"
        else:
            self._latex_name = latex_name
        self._one = None # to be set by self.one() if self is an endomorphism
                         # set (fmodule1 = fmodule2)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: H = Hom(M,N)
            sage: H._latex_()
            '\\mathrm{Hom}\\left(M,N\\right)'
            sage: latex(H)  # indirect doctest
            \mathrm{Hom}\left(M,N\right)

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def __call__(self, *args, **kwds):
        r"""
        To bypass Homset.__call__, enforcing Parent.__call__ instead.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 3, name='N')
            sage: H = Hom(M,N)
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: a = H.__call__(0) ; a
            Generic morphism:
              From: Rank-2 free module M over the Integer Ring
              To:   Rank-3 free module N over the Integer Ring
            sage: a.matrix(e,f)
            [0 0]
            [0 0]
            [0 0]
            sage: a == H.zero()
            True
            sage: a == H(0)
            True
            sage: a = H.__call__([[1,2],[3,4],[5,6]], bases=(e,f), name='a') ; a
            Generic morphism:
              From: Rank-2 free module M over the Integer Ring
              To:   Rank-3 free module N over the Integer Ring
            sage: a.matrix(e,f)
            [1 2]
            [3 4]
            [5 6]
            sage: a == H([[1,2],[3,4],[5,6]], bases=(e,f))
            True

        """
        from sage.structure.parent import Parent
        return Parent.__call__(self, *args, **kwds)

    #### Methods required for any Parent

    def _element_constructor_(self, matrix_rep, bases=None, name=None,
                              latex_name=None, is_identity=False):
        r"""
        Construct an element of ``self``, i.e. a homomorphism M --> N, where
        M is the domain of ``self`` and N its codomain.

        INPUT:

        - ``matrix_rep`` -- matrix representation of the homomorphism with
          respect to the bases ``basis1`` and ``basis2``; this entry can
          actually be any material from which a matrix of size rank(N)*rank(M)
          can be constructed
        - ``bases`` -- (default: ``None``) pair (basis_M, basis_N) defining the
          matrix representation, basis_M being a basis of module `M` and
          basis_N a basis of module `N` ; if none is provided the pair formed
          by the default bases of each module is assumed.
        - ``name`` -- (default: ``None``) string; name given to the
          homomorphism
        - ``latex_name`` -- (default: ``None``)string;  LaTeX symbol to denote
          the homomorphism; if none is provided, ``name`` will be used.
        - ``is_identity`` -- (default: ``False``) determines whether the
          constructed object is the identity endomorphism; if set to ``True``,
          then N must be M and the entry ``matrix_rep`` is not used.

        EXAMPLES:

        Construction of a homomorphism between two free `\ZZ`-modules::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: H = Hom(M,N)
            sage: phi = H._element_constructor_([[2,-1,3], [1,0,-4]], bases=(e,f),
            ....:                               name='phi', latex_name=r'\phi')
            sage: phi
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: phi.matrix(e,f)
            [ 2 -1  3]
            [ 1  0 -4]
            sage: phi == H([[2,-1,3], [1,0,-4]], bases=(e,f), name='phi',
            ....:          latex_name=r'\phi')
            True

        Construction of an endomorphism::

            sage: EM = End(M)
            sage: phi = EM._element_constructor_([[1,2,3],[4,5,6],[7,8,9]],
            ....:                                name='phi', latex_name=r'\phi')
            sage: phi
            Generic endomorphism of Rank-3 free module M over the Integer Ring
            sage: phi.matrix(e,e)
            [1 2 3]
            [4 5 6]
            [7 8 9]

        Coercion of a type-`(1,1)` tensor to an endomorphism::

            sage: a = M.tensor((1,1))
            sage: a[:] = [[1,2,3],[4,5,6],[7,8,9]]
            sage: EM = End(M)
            sage: phi_a = EM._element_constructor_(a) ; phi_a
            Generic endomorphism of Rank-3 free module M over the Integer Ring
            sage: phi_a.matrix(e,e)
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: phi_a == phi
            True
            sage: phi_a1 = EM(a) ; phi_a1  # indirect doctest
            Generic endomorphism of Rank-3 free module M over the Integer Ring
            sage: phi_a1 == phi
            True

        """
        if isinstance(matrix_rep, FreeModuleTensor):
            # coercion of a type-(1,1) tensor to an endomorphism
            # (this includes automorphisms, since the class
            #  FreeModuleAutomorphism inherits from FreeModuleTensor)
            tensor = matrix_rep # for readability
            if tensor.tensor_type() == (1,1) and \
                                         self.is_endomorphism_set() and \
                                         tensor.base_module() is self.domain():
                basis = tensor.pick_a_basis()
                tcomp = tensor.comp(basis)
                fmodule = tensor.base_module()
                mat = [[ tcomp[[i,j]] for j in fmodule.irange()] \
                                                     for i in fmodule.irange()]
                if isinstance(tensor, FreeModuleAutomorphism):
                    is_identity = tensor._is_identity
                else:
                    is_identity = False
                resu = self.element_class(self, mat, bases=(basis,basis),
                              name=tensor._name, latex_name=tensor._latex_name,
                              is_identity=is_identity)
            else:
                raise TypeError("cannot coerce the {}".format(tensor) +
                                " to an element of {}".format(self))
        else:
            # Standard construction:
            resu = self.element_class(self, matrix_rep, bases=bases, name=name,
                                      latex_name=latex_name,
                                      is_identity=is_identity)
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) element.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = Hom(M,N)._an_element_() ; phi
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: phi.matrix(e,f)
            [1 1 1]
            [1 1 1]
            sage: phi == Hom(M,N).an_element()
            True

        """
        ring = self.base_ring()
        m = self.domain().rank()
        n = self.codomain().rank()
        matrix_rep = [[ring.an_element() for i in range(m)] for j in range(n)]
        return self.element_class(self, matrix_rep)

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent.

        EXAMPLES:

        The module of type-`(1,1)` tensors coerces to ``self``, if the latter
        is some endomorphism set::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: End(M)._coerce_map_from_(M.tensor_module(1,1))
            True
            sage: End(M).has_coerce_map_from(M.tensor_module(1,1))
            True
            sage: End(M)._coerce_map_from_(M.tensor_module(1,2))
            False

        The general linear group coerces to the endomorphism ring::

            sage: End(M)._coerce_map_from_(M.general_linear_group())
            True

        """
        from sage.tensor.modules.tensor_free_module import TensorFreeModule
        from sage.tensor.modules.free_module_linear_group import \
                                                          FreeModuleLinearGroup
        if isinstance(other, TensorFreeModule):
            # Coercion of a type-(1,1) tensor to an endomorphism:
            if other.tensor_type() == (1,1):
                return self.is_endomorphism_set() and \
                                           other.base_module() is self.domain()
        if isinstance(other, FreeModuleLinearGroup):
            # Coercion of an automorphism to an endomorphism:
            return self.is_endomorphism_set() and \
                                           other.base_module() is self.domain()
        return False

    #### End of methods required for any Parent


    #### Monoid methods (case of an endomorphism set) ####

    def one(self):
        r"""
        Return the identity element of ``self`` considered as a monoid (case of
        an endomorphism set).

        This applies only when the codomain of ``self`` is equal to its domain,
        i.e. when ``self`` is of the type `\mathrm{Hom}(M,M)`.

        OUTPUT:

        - the identity element of `\mathrm{End}(M) = \mathrm{Hom}(M,M)`, as an
          instance of
          :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`

        EXAMPLES:

        Identity element of the set of endomorphisms of a free module
        over `\ZZ`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: H = End(M)
            sage: H.one()
            Identity endomorphism of Rank-3 free module M over the Integer Ring
            sage: H.one().matrix(e)
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: H.one().is_identity()
            True

        NB: mathematically, ``H.one()`` coincides with the identity map of the
        free module `M`. However the latter is considered here as an
        element of `\mathrm{GL}(M)`, the general linear group of `M`.
        Accordingly, one has to use the coercion map
        `\mathrm{GL}(M) \rightarrow \mathrm{End}(M)`
        to recover ``H.one()`` from ``M.identity_map()``::

            sage: M.identity_map()
            Identity map of the Rank-3 free module M over the Integer Ring
            sage: M.identity_map().parent()
            General linear group of the Rank-3 free module M over the Integer Ring
            sage: H.one().parent()
            Set of Morphisms from Rank-3 free module M over the Integer Ring
             to Rank-3 free module M over the Integer Ring
             in Category of finite dimensional modules over Integer Ring
            sage: H.one() == H(M.identity_map())
            True

        Conversely, one can recover ``M.identity_map()`` from ``H.one()`` by
        means of a conversion `\mathrm{End}(M)\rightarrow \mathrm{GL}(M)`::

            sage: GL = M.general_linear_group()
            sage: M.identity_map() == GL(H.one())
            True

        """
        if self._one is None:
            if self.codomain() != self.domain():
                raise TypeError("the {} is not a monoid".format(self))
            self._one = self.element_class(self, [], is_identity=True)
        return self._one

    #### End of monoid methods ####
