r"""
General linear group of a free module

The set `\mathrm{GL}(M)` of automorphisms (i.e. invertible endomorphisms) of a
free module of finite rank `M` is a group under composition of automorphisms,
named the *general linear group* of `M`. In other words, `\mathrm{GL}(M)` is
the group of units (i.e. invertible elements) of `\mathrm{End}(M)`, the
endomorphism ring of `M`.

The group `\mathrm{GL}(M)` is implemented via the class
:class:`FreeModuleLinearGroup`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Michael Jung (2019): improve treatment of the identity element

REFERENCES:

- Chap. 15 of R. Godement : *Algebra* [God1968]_
"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.groups import Groups
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism


class FreeModuleLinearGroup(UniqueRepresentation, Parent):
    r"""
    General linear group of a free module of finite rank over a commutative
    ring.

    Given a free module of finite rank `M` over a commutative ring `R`, the
    *general linear group* of `M` is the group `\mathrm{GL}(M)` of
    automorphisms (i.e. invertible endomorphisms) of `M`. It is the group of
    units (i.e. invertible elements) of `\mathrm{End}(M)`, the endomorphism
    ring of `M`.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank over a commutative ring
      `R`, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`

    EXAMPLES:

    General linear group of a free `\ZZ`-module of rank 3::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: from sage.tensor.modules.free_module_linear_group import FreeModuleLinearGroup
        sage: GL = FreeModuleLinearGroup(M) ; GL
        General linear group of the Rank-3 free module M over the Integer Ring

    Instead of importing FreeModuleLinearGroup in the global name space, it is
    recommended to use the module's method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.general_linear_group`::

        sage: GL = M.general_linear_group() ; GL
        General linear group of the Rank-3 free module M over the Integer Ring
        sage: latex(GL)
        \mathrm{GL}\left( M \right)

    As most parents, the general linear group has a unique instance::

        sage: GL is M.general_linear_group()
        True

    `\mathrm{GL}(M)` is in the category of groups::

        sage: GL.category()
        Category of groups
        sage: GL in Groups()
        True

    ``GL`` is a *parent* object, whose elements are automorphisms of `M`,
    represented by instances of the class
    :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`::

        sage: GL.Element
        <class 'sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism'>
        sage: a = GL.an_element() ; a
        Automorphism of the Rank-3 free module M over the Integer Ring
        sage: a.matrix(e)
        [ 1  0  0]
        [ 0 -1  0]
        [ 0  0  1]
        sage: a in GL
        True
        sage: GL.is_parent_of(a)
        True

    As an endomorphism, ``a`` maps elements of `M` to elements of `M`::

        sage: v = M.an_element() ; v
        Element of the Rank-3 free module M over the Integer Ring
        sage: v.display()
        e_0 + e_1 + e_2
        sage: a(v)
        Element of the Rank-3 free module M over the Integer Ring
        sage: a(v).display()
        e_0 - e_1 + e_2

    An automorphism can also be viewed as a tensor of type `(1,1)` on `M`::

        sage: a.tensor_type()
        (1, 1)
        sage: a.display(e)
        e_0⊗e^0 - e_1⊗e^1 + e_2⊗e^2
        sage: type(a)
        <class 'sage.tensor.modules.free_module_linear_group.FreeModuleLinearGroup_with_category.element_class'>

    As for any group, the identity element is obtained by the method
    :meth:`one`::

        sage: id = GL.one() ; id
        Identity map of the Rank-3 free module M over the Integer Ring
        sage: id*a == a
        True
        sage: a*id == a
        True
        sage: a*a^(-1) == id
        True
        sage: a^(-1)*a == id
        True

    The identity element is of course the identity map of the module `M`::

        sage: id(v) == v
        True
        sage: id.matrix(e)
        [1 0 0]
        [0 1 0]
        [0 0 1]

    The module's changes of basis are stored as elements of the general linear
    group::

        sage: f = M.basis('f', from_family=(-e[1], 4*e[0]+3*e[2], 7*e[0]+5*e[2]))
        sage: f
        Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring
        sage: M.change_of_basis(e,f)
        Automorphism of the Rank-3 free module M over the Integer Ring
        sage: M.change_of_basis(e,f) in GL
        True
        sage: M.change_of_basis(e,f).parent()
        General linear group of the Rank-3 free module M over the Integer Ring
        sage: M.change_of_basis(e,f).matrix(e)
        [ 0  4  7]
        [-1  0  0]
        [ 0  3  5]
        sage: M.change_of_basis(e,f) == M.change_of_basis(f,e).inverse()
        True

    Since every automorphism is an endomorphism, there is a coercion
    `\mathrm{GL}(M) \rightarrow \mathrm{End}(M)` (the endomorphism ring of
    module `M`)::

        sage: End(M).has_coerce_map_from(GL)
        True

    (see :class:`~sage.tensor.modules.free_module_homset.FreeModuleHomset` for
    details), but not in the reverse direction, since only bijective
    endomorphisms are automorphisms::

        sage: GL.has_coerce_map_from(End(M))
        False

    A bijective endomorphism can be converted to an element of
    `\mathrm{GL}(M)`::

        sage: h = M.endomorphism([[1,0,0], [0,-1,2], [0,1,-3]]) ; h
        Generic endomorphism of Rank-3 free module M over the Integer Ring
        sage: h.parent() is End(M)
        True
        sage: ah = GL(h) ; ah
        Automorphism of the Rank-3 free module M over the Integer Ring
        sage: ah.parent() is GL
        True

    As maps `M\rightarrow M`, ``ah`` and ``h`` are identical::

        sage: v  # recall
        Element of the Rank-3 free module M over the Integer Ring
        sage: ah(v) == h(v)
        True
        sage: ah.matrix(e) == h.matrix(e)
        True

    Of course, non-invertible endomorphisms cannot be converted to elements of
    `\mathrm{GL}(M)`::

        sage: GL(M.endomorphism([[0,0,0], [0,-1,2], [0,1,-3]]))
        Traceback (most recent call last):
        ...
        TypeError: the Generic endomorphism of Rank-3 free module M over the
         Integer Ring is not invertible

    Similarly, there is a coercion `\mathrm{GL}(M)\rightarrow T^{(1,1)}(M)`
    (module of type-`(1,1)` tensors)::

        sage: M.tensor_module(1,1).has_coerce_map_from(GL)
        True

    (see :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule` for
    details), but not in the reverse direction, since not every type-`(1,1)`
    tensor can be considered as an automorphism::

        sage: GL.has_coerce_map_from(M.tensor_module(1,1))
        False

    Invertible type-`(1,1)` tensors can be converted to automorphisms::

        sage: t = M.tensor((1,1), name='t')
        sage: t[e,:] = [[-1,0,0], [0,1,2], [0,1,3]]
        sage: at = GL(t) ; at
        Automorphism t of the Rank-3 free module M over the Integer Ring
        sage: at.matrix(e)
        [-1  0  0]
        [ 0  1  2]
        [ 0  1  3]
        sage: at.matrix(e) == t[e,:]
        True

    Non-invertible ones cannot::

        sage: t0 = M.tensor((1,1), name='t_0')
        sage: t0[e,0,0] = 1
        sage: t0[e,:]  # the matrix is clearly not invertible
        [1 0 0]
        [0 0 0]
        [0 0 0]
        sage: GL(t0)
        Traceback (most recent call last):
        ...
        TypeError: the Type-(1,1) tensor t_0 on the Rank-3 free module M over
         the Integer Ring is not invertible
        sage: t0[e,1,1], t0[e,2,2] = 2, 3
        sage: t0[e,:]  # the matrix is not invertible in Mat_3(ZZ)
        [1 0 0]
        [0 2 0]
        [0 0 3]
        sage: GL(t0)
        Traceback (most recent call last):
        ...
        TypeError: the Type-(1,1) tensor t_0 on the Rank-3 free module M over
         the Integer Ring is not invertible
    """

    Element = FreeModuleAutomorphism

    def __init__(self, fmodule):
        r"""
        See :class:`FreeModuleLinearGroup` for documentation and examples.

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: from sage.tensor.modules.free_module_linear_group import FreeModuleLinearGroup
            sage: GL = FreeModuleLinearGroup(M) ; GL
            General linear group of the Rank-3 free module M over the Integer Ring
            sage: GL.category()
            Category of groups
            sage: TestSuite(GL).run()
        """
        if not isinstance(fmodule, FiniteRankFreeModule):
            raise TypeError("{} is not a free module of finite rank".format(
                            fmodule))
        Parent.__init__(self, category=Groups())
        self._fmodule = fmodule

    #### Parent methods ####

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct a free module automorphism.

        INPUT:

        - ``comp`` -- (default: ``[]``) components representing the
          automorphism with respect to ``basis``; this entry can actually be
          any array of size rank(M)*rank(M) from which a matrix of elements
          of ``self`` base ring can be constructed; the *columns* of ``comp``
          must be the components w.r.t. ``basis`` of the images of the elements
          of ``basis``. If ``comp`` is ``[]``, the automorphism has to be
          initialized afterwards by method
          :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.set_comp`
          or via the operator [].
        - ``basis`` -- (default: ``None``) basis of ``self`` defining the
          matrix representation; if ``None`` the default basis of ``self`` is
          assumed.
        - ``name`` -- (default: ``None``) name given to the automorphism
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          automorphism; if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`

        EXAMPLES:

        Generic construction::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: GL = M.general_linear_group()
            sage: a = GL._element_constructor_(comp=[[1,2],[1,3]], basis=e,
            ....:                              name='a')
            sage: a
            Automorphism a of the Rank-2 free module M over the Integer Ring
            sage: a.matrix(e)
            [1 2]
            [1 3]

        Identity map constructed from integer 1::

            sage: GL._element_constructor_(1)
            Identity map of the Rank-2 free module M over the Integer Ring
            sage: GL._element_constructor_(1).matrix(e)
            [1 0]
            [0 1]

        Construction from an invertible endomorphism::

            sage: phi = M.endomorphism([[1,1], [2,3]])
            sage: a = GL._element_constructor_(phi) ; a
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: a.matrix(e)
            [1 1]
            [2 3]
            sage: a.matrix(e) == phi.matrix(e)
            True

        Construction from an invertible tensor of type `(1,1)`::

            sage: t = M.tensor((1,1), name='t')
            sage: t[e,:] = [[1,1], [2,3]]
            sage: a = GL._element_constructor_(t) ; a
            Automorphism t of the Rank-2 free module M over the Integer Ring
            sage: a.matrix(e) == t[e,:]
            True
        """
        from sage.tensor.modules.free_module_tensor import FreeModuleTensor
        from sage.tensor.modules.free_module_morphism import \
                                                   FiniteRankFreeModuleMorphism
        if comp == 1:
            return self.one()
        if isinstance(comp, FreeModuleTensor):
            tens = comp # for readability
            # Conversion of a type-(1,1) tensor to an automorphism
            if tens.tensor_type() == (1,1):
                resu = self.element_class(self._fmodule, name=tens._name,
                                          latex_name=tens._latex_name)
                for basis, comp in tens._components.items():
                    resu._components[basis] = comp.copy()
                # Check whether the tensor is invertible:
                try:
                    resu.inverse()
                except (ZeroDivisionError, TypeError):
                    raise TypeError("the {} is not invertible ".format(tens))
                return resu
            else:
                    raise TypeError("the {} cannot be converted ".format(tens)
                                    + "to an automorphism.")
        if isinstance(comp, FiniteRankFreeModuleMorphism):
            # Conversion of an endomorphism to an automorphism
            endo = comp  # for readability
            if endo.is_endomorphism() and self._fmodule is endo.domain():
                resu = self.element_class(self._fmodule, name=endo._name,
                                          latex_name=endo._latex_name)
                for basis, mat in endo._matrices.items():
                    resu.add_comp(basis[0])[:] = mat
                # Check whether the endomorphism is invertible:
                try:
                    resu.inverse()
                except (ZeroDivisionError, TypeError):
                    raise TypeError("the {} is not invertible ".format(endo))
                return resu
            else:
                raise TypeError("cannot coerce the {}".format(endo) +
                                " to an element of {}".format(self))

        # standard construction
        resu = self.element_class(self._fmodule, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some specific free module automorphism.

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: GL = M.general_linear_group()
            sage: a = GL._an_element_() ; a
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: a.matrix(e)
            [ 1  0]
            [ 0 -1]
        """
        resu = self.element_class(self._fmodule)
        # Make sure that the base module has a default basis
        self._fmodule.an_element()
        comp = resu.set_comp()
        for i in self._fmodule.irange():
            if i % 2 == 0:
                comp[[i,i]] = self._fmodule._ring.one()
            else:
                comp[[i,i]] = -(self._fmodule._ring.one())
        return resu

    #### End of parent methods ####

    #### Monoid methods ####

    @cached_method
    def one(self):
        r"""
        Return the group identity element of ``self``.

        The group identity element is nothing but the module identity map.

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          representing the identity element.

        EXAMPLES:

        Identity element of the general linear group of a rank-2 free module::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M', start_index=1)
            sage: GL = M.general_linear_group()
            sage: GL.one()
            Identity map of the Rank-2 free module M over the Integer Ring

        The identity element is cached::

            sage: GL.one() is GL.one()
            True

        Check that the element returned is indeed the neutral element for
        the group law::

            sage: e = M.basis('e')
            sage: a = GL([[3,4],[5,7]], basis=e) ; a
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: a.matrix(e)
            [3 4]
            [5 7]
            sage: GL.one() * a == a
            True
            sage: a * GL.one() == a
            True
            sage: a * a^(-1) == GL.one()
            True
            sage: a^(-1) * a == GL.one()
            True

        The unit element of `\mathrm{GL}(M)` is the identity map of `M`::

            sage: GL.one()(e[1])
            Element e_1 of the Rank-2 free module M over the Integer Ring
            sage: GL.one()(e[2])
            Element e_2 of the Rank-2 free module M over the Integer Ring

        Its matrix is the identity matrix in any basis::

            sage: GL.one().matrix(e)
            [1 0]
            [0 1]
            sage: f = M.basis('f', from_family=(e[1]+2*e[2], e[1]+3*e[2]))
            sage: GL.one().matrix(f)
            [1 0]
            [0 1]
        """
        resu = self._element_constructor_(name='Id', latex_name=r'\mathrm{Id}')
        # Initialization of the components (Kronecker delta) in some basis:
        from .comp import KroneckerDelta
        fmodule = self._fmodule
        for basis in fmodule.bases():
            resu._components[basis] = KroneckerDelta(fmodule._ring, basis,
                                    start_index=fmodule._sindex,
                                    output_formatter=fmodule._output_formatter)
        resu._is_identity = True
        resu.set_immutable()
        return resu

    #### End of monoid methods ####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: GL = M.general_linear_group()
            sage: GL._repr_()
            'General linear group of the Rank-2 free module M over the Integer Ring'
        """
        return "General linear group of the {}".format(self._fmodule)

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: GL = M.general_linear_group()
            sage: GL._latex_()
            \mathrm{GL}\left( M \right)
        """
        from sage.misc.latex import latex
        return r"\mathrm{GL}\left(" + latex(self._fmodule) + r"\right)"

    def base_module(self):
        r"""
        Return the free module of which ``self`` is the general linear group.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module of which ``self`` is the general linear group

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: GL = M.general_linear_group()
            sage: GL.base_module()
            Rank-2 free module M over the Integer Ring
            sage: GL.base_module() is M
            True
        """
        return self._fmodule
