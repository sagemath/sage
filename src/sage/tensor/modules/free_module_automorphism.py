r"""
Free module automorphisms

Given a free module `M` of finite rank over a commutative ring `R`, an
*automorphism* of `M` is a map

.. MATH::

    \phi:\ M \longrightarrow M

that is linear (i.e. is a module homomorphism) and bijective.

Automorphisms of a free module of finite rank are implemented via the class
:class:`FreeModuleAutomorphism`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Michael Jung (2019): improve treatment of the identity element

REFERENCES:

- Chaps. 15, 24 of R. Godement: *Algebra* [God1968]_
"""
# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.element import MultiplicativeGroupElement
from sage.tensor.modules.free_module_tensor import FreeModuleTensor


class FreeModuleAutomorphism(FreeModuleTensor, MultiplicativeGroupElement):
    r"""
    Automorphism of a free module of finite rank over a commutative ring.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.tensor.modules.free_module_linear_group.FreeModuleLinearGroup`.

    This class inherits from the classes
    :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
    and
    :class:`~sage.structure.element.MultiplicativeGroupElement`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank over a commutative ring
      `R`, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``name`` -- (default: ``None``) name given to the automorphism
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      automorphism; if none is provided, the LaTeX symbol is set to ``name``
    - ``is_identity`` -- boolean (default: ``False``); determines whether the
      constructed object is the identity automorphism, i.e. the identity map
      of `M` considered as an automorphism (the identity element of the
      general linear group)

    EXAMPLES:

    Automorphism of a rank-2 free module over `\ZZ`::

        sage: M = FiniteRankFreeModule(ZZ, 2, name='M', start_index=1)
        sage: a = M.automorphism(name='a', latex_name=r'\alpha') ; a
        Automorphism a of the Rank-2 free module M over the Integer Ring
        sage: a.parent()
        General linear group of the Rank-2 free module M over the Integer Ring
        sage: a.parent() is M.general_linear_group()
        True
        sage: latex(a)
        \alpha

    Setting the components of ``a`` w.r.t. a basis of module ``M``::

        sage: e = M.basis('e') ; e
        Basis (e_1,e_2) on the Rank-2 free module M over the Integer Ring
        sage: a[:] = [[1,2],[1,3]]
        sage: a.matrix(e)
        [1 2]
        [1 3]
        sage: a(e[1]).display()
        a(e_1) = e_1 + e_2
        sage: a(e[2]).display()
        a(e_2) = 2 e_1 + 3 e_2

    Actually, the components w.r.t. a given basis can be specified at the
    construction of the object::

        sage: a = M.automorphism(matrix=[[1,2],[1,3]], basis=e, name='a',
        ....:                    latex_name=r'\alpha') ; a
        Automorphism a of the Rank-2 free module M over the Integer Ring
        sage: a.matrix(e)
        [1 2]
        [1 3]

    Since e is the module's default basis, it can be omitted in the argument
    list::

        sage: a == M.automorphism(matrix=[[1,2],[1,3]], name='a',
        ....:                     latex_name=r'\alpha')
        True

    The matrix of the automorphism can be obtained in any basis::

        sage: f = M.basis('f', from_family=(3*e[1]+4*e[2], 5*e[1]+7*e[2])) ; f
        Basis (f_1,f_2) on the Rank-2 free module M over the Integer Ring
        sage: a.matrix(f)
        [2 3]
        [1 2]

    Automorphisms are tensors of type `(1,1)`::

        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2

    In particular, they can be displayed as such::

        sage: a.display(e)
        a = e_1⊗e^1 + 2 e_1⊗e^2 + e_2⊗e^1 + 3 e_2⊗e^2
        sage: a.display(f)
        a = 2 f_1⊗f^1 + 3 f_1⊗f^2 + f_2⊗f^1 + 2 f_2⊗f^2

    The automorphism acting on a module element::

        sage: v = M([-2,3], name='v') ; v
        Element v of the Rank-2 free module M over the Integer Ring
        sage: a(v)
        Element a(v) of the Rank-2 free module M over the Integer Ring
        sage: a(v).display()
        a(v) = 4 e_1 + 7 e_2

    A second automorphism of the module ``M``::

        sage: b = M.automorphism([[0,1],[-1,0]], name='b') ; b
        Automorphism b of the Rank-2 free module M over the Integer Ring
        sage: b.matrix(e)
        [ 0  1]
        [-1  0]
        sage: b(e[1]).display()
        b(e_1) = -e_2
        sage: b(e[2]).display()
        b(e_2) = e_1

    The composition of automorphisms is performed via the multiplication
    operator::

        sage: s = a*b ; s
        Automorphism of the Rank-2 free module M over the Integer Ring
        sage: s(v) == a(b(v))
        True
        sage: s.matrix(f)
        [ 11  19]
        [ -7 -12]
        sage: s.matrix(f) == a.matrix(f) * b.matrix(f)
        True

    It is not commutative::

        sage: a*b != b*a
        True

    In other words, the parent of ``a`` and ``b``, i.e. the group
    `\mathrm{GL}(M)`, is not abelian::

        sage: M.general_linear_group() in CommutativeAdditiveGroups()
        False

    The neutral element for the composition law is the module identity map::

        sage: id = M.identity_map() ; id
        Identity map of the Rank-2 free module M over the Integer Ring
        sage: id.parent()
        General linear group of the Rank-2 free module M over the Integer Ring
        sage: id(v) == v
        True
        sage: id.matrix(f)
        [1 0]
        [0 1]
        sage: id*a == a
        True
        sage: a*id == a
        True

    The inverse of an automorphism is obtained via the method :meth:`inverse`,
    or the operator ~, or the exponent -1::

        sage: a.inverse()
        Automorphism a^(-1) of the Rank-2 free module M over the Integer Ring
        sage: a.inverse() is ~a
        True
        sage: a.inverse() is a^(-1)
        True
        sage: (a^(-1)).matrix(e)
        [ 3 -2]
        [-1  1]
        sage: a*a^(-1) == id
        True
        sage: a^(-1)*a == id
        True
        sage: a^(-1)*s == b
        True
        sage: (a^(-1))(a(v)) == v
        True

    The module's changes of basis are stored as automorphisms::

        sage: M.change_of_basis(e,f)
        Automorphism of the Rank-2 free module M over the Integer Ring
        sage: M.change_of_basis(e,f).parent()
        General linear group of the Rank-2 free module M over the Integer Ring
        sage: M.change_of_basis(e,f).matrix(e)
        [3 5]
        [4 7]
        sage: M.change_of_basis(f,e) == M.change_of_basis(e,f).inverse()
        True

    The opposite of an automorphism is still an automorphism::

        sage: -a
        Automorphism -a of the Rank-2 free module M over the Integer Ring
        sage: (-a).parent()
        General linear group of the Rank-2 free module M over the Integer Ring
        sage: (-a).matrix(e) == - (a.matrix(e))
        True

    Adding two automorphisms results in a generic type-`(1,1)` tensor::

        sage: s = a + b ; s
        Type-(1,1) tensor a+b on the Rank-2 free module M over the Integer Ring
        sage: s.parent()
        Free module of type-(1,1) tensors on the Rank-2 free module M over the
         Integer Ring
        sage: a[:], b[:], s[:]
        (
        [1 2]  [ 0  1]  [1 3]
        [1 3], [-1  0], [0 3]
        )

    To get the result as an endomorphism, one has to explicitly convert it via
    the parent of endomorphisms, `\mathrm{End}(M)`::

        sage: s = End(M)(a+b) ; s
        Generic endomorphism of Rank-2 free module M over the Integer Ring
        sage: s(v) == a(v) + b(v)
        True
        sage: s.matrix(e) == a.matrix(e) + b.matrix(e)
        True
        sage: s.matrix(f) == a.matrix(f) + b.matrix(f)
        True
    """
    def __init__(self, fmodule, name=None, latex_name=None):
        r"""
        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism
            sage: a = FreeModuleAutomorphism(M, name='a')
            sage: a[e,:] = [[-1,0,0],[0,1,2],[0,1,3]]
            sage: TestSuite(a).run(skip='_test_category') # see below

        In the above test suite, _test_category fails because a is not an
        instance of a.parent().category().element_class. Actually automorphism
        must be constructed via FreeModuleLinearGroup.element_class and
        not by a direct call to FreeModuleAutomorphism::

            sage: a = M.general_linear_group().element_class(M, name='a')
            sage: a[e,:] = [[-1,0,0],[0,1,2],[0,1,3]]
            sage: TestSuite(a).run()

        Test suite on the identity map::

            sage: id = M.general_linear_group().one()
            sage: TestSuite(id).run()

        Test suite on the automorphism obtained as GL.an_element()::

            sage: b = M.general_linear_group().an_element()
            sage: TestSuite(b).run()
        """
        FreeModuleTensor.__init__(self, fmodule, (1,1), name=name,
                                  latex_name=latex_name,
                                  parent=fmodule.general_linear_group())
        # MultiplicativeGroupElement attributes:
        # - none
        # Local attributes:
        self._is_identity = False # a priori
        self._inverse = None    # inverse automorphism not set yet
        self._matrices = {}

    #### SageObject methods ####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: M.automorphism()
            Automorphism of the 3-dimensional vector space M over the Rational Field
            sage: M.automorphism(name='a')
            Automorphism a of the 3-dimensional vector space M over the Rational Field
            sage: M.identity_map()
            Identity map of the 3-dimensional vector space M over the Rational Field
        """
        if self._is_identity:
            description = "Identity map "
        else:
            description = "Automorphism "
            if self._name is not None:
                description += self._name + " "
        description += "of the {}".format(self._fmodule)
        return description

    #### End of SageObject methods ####

    #### FreeModuleTensor methods ####

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.automorphism(name='a')
            sage: a._new_instance()
            Automorphism of the Rank-3 free module M over the Integer Ring
            sage: Id = M.identity_map()
            sage: Id._new_instance()
            Automorphism of the Rank-3 free module M over the Integer Ring
        """
        return self.__class__(self._fmodule)

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(name='a')
            sage: a[e,:] = [[1,0,-1], [0,3,0], [0,0,2]]
            sage: b = a.inverse()
            sage: a._inverse
            Automorphism a^(-1) of the 3-dimensional vector space M over the
             Rational Field
            sage: a._del_derived()
            sage: a._inverse  # has been reset to None
        """
        # First delete the derived quantities pertaining to FreeModuleTensor:
        FreeModuleTensor._del_derived(self)
        # Then reset the inverse automorphism to None:
        if self._inverse is not None:
            self._inverse._inverse = None  # (it was set to self)
            self._inverse = None
        # and delete the matrices:
        self._matrices.clear()

    def set_comp(self, basis=None):
        r"""
        Return the components of ``self`` w.r.t. a given module basis for
        assignment.

        The components with respect to other bases are deleted, in order to
        avoid any inconsistency. To keep them, use the method :meth:`add_comp`
        instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        EXAMPLES:

        Setting the components of an automorphism of a rank-3 free
        `\ZZ`-module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(name='a')
            sage: a.set_comp(e)
            2-indices components w.r.t. Basis (e_0,e_1,e_2) on the Rank-3 free
             module M over the Integer Ring
            sage: a.set_comp(e)[:] = [[1,0,0],[0,1,2],[0,1,3]]
            sage: a.matrix(e)
            [1 0 0]
            [0 1 2]
            [0 1 3]

        Since ``e`` is the module's default basis, one has::

            sage: a.set_comp() is a.set_comp(e)
            True

        The method :meth:`set_comp` can be used to modify a single component::

            sage: a.set_comp(e)[0,0] = -1
            sage: a.matrix(e)
            [-1  0  0]
            [ 0  1  2]
            [ 0  1  3]

        A short cut to :meth:`set_comp` is the bracket operator, with the basis
        as first argument::

            sage: a[e,:] = [[1,0,0],[0,-1,2],[0,1,-3]]
            sage: a.matrix(e)
            [ 1  0  0]
            [ 0 -1  2]
            [ 0  1 -3]
            sage: a[e,0,0] = -1
            sage: a.matrix(e)
            [-1  0  0]
            [ 0 -1  2]
            [ 0  1 -3]

        The call to :meth:`set_comp` erases the components previously defined
        in other bases; to keep them, use the method :meth:`add_comp` instead::

            sage: f = M.basis('f', from_family=(-e[0], 3*e[1]+4*e[2],
            ....:                               5*e[1]+7*e[2])) ; f
            Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer
             Ring
            sage: a._components
            {Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer
             Ring: 2-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring}
            sage: a.set_comp(f)[:] = [[-1,0,0], [0,1,0], [0,0,-1]]

        The components w.r.t. basis ``e`` have been erased::

            sage: a._components
            {Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer
             Ring: 2-indices components w.r.t. Basis (f_0,f_1,f_2) on the
             Rank-3 free module M over the Integer Ring}

        Of course, they can be computed from those in basis ``f`` by means of
        a change-of-basis formula, via the method :meth:`comp` or
        :meth:`matrix`::

            sage: a.matrix(e)
            [ -1   0   0]
            [  0  41 -30]
            [  0  56 -41]

        For the identity map, it is not permitted to set components::

            sage: id = M.identity_map()
            sage: id.set_comp(e)
            Traceback (most recent call last):
            ...
            ValueError: the components of the identity map cannot be changed

        Indeed, the components are set automatically::

            sage: id.comp(e)
            Kronecker delta of size 3x3
            sage: id.comp(f)
            Kronecker delta of size 3x3
        """
        if self._is_identity:
            raise ValueError("the components of the identity map cannot be "
                             "changed")
        return FreeModuleTensor._set_comp_unsafe(self, basis=basis)

    def add_comp(self, basis=None):
        r"""
        Return the components of ``self`` w.r.t. a given module basis for
        assignment, keeping the components w.r.t. other bases.

        To delete the components w.r.t. other bases, use the method
        :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        .. WARNING::

            If the automorphism has already components in other bases, it
            is the user's responsibility to make sure that the components
            to be added are consistent with them.

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`;
          if such components did not exist previously, they are created

        EXAMPLES:

        Adding components to an automorphism of a rank-3 free
        `\ZZ`-module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(name='a')
            sage: a[e,:] = [[1,0,0],[0,-1,2],[0,1,-3]]
            sage: f = M.basis('f', from_family=(-e[0], 3*e[1]+4*e[2],
            ....:                               5*e[1]+7*e[2])) ; f
            Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer
             Ring
            sage: a.add_comp(f)[:] = [[1,0,0], [0, 80, 143], [0, -47, -84]]

        The components in basis ``e`` have been kept::

            sage: a._components # random (dictionary output)
            {Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer
             Ring: 2-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer
             Ring: 2-indices components w.r.t. Basis (f_0,f_1,f_2) on the
             Rank-3 free module M over the Integer Ring}

        For the identity map, it is not permitted to invoke :meth:`add_comp`::

            sage: id = M.identity_map()
            sage: id.add_comp(e)
            Traceback (most recent call last):
            ...
            ValueError: the components of the identity map cannot be changed

        Indeed, the components are set automatically::

            sage: id.comp(e)
            Kronecker delta of size 3x3
            sage: id.comp(f)
            Kronecker delta of size 3x3
        """
        if self._is_identity:
            raise ValueError("the components of the identity map cannot be "
                             "changed")
        return FreeModuleTensor._add_comp_unsafe(self, basis=basis)

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleTensor.__call__` to allow for a single
        argument (module element).

        EXAMPLES:

        Call with a single argument: return a module element::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism([[-1,0,0],[0,1,2],[0,1,3]], name='a')
            sage: v = M([2,1,4], name='v')
            sage: s = a.__call__(v) ; s
            Element a(v) of the Rank-3 free module M over the Integer Ring
            sage: s.display()
            a(v) = -2 e_1 + 9 e_2 + 13 e_3
            sage: s == a(v)
            True
            sage: s == a.contract(v)
            True

        Call with two arguments (:class:`FreeModuleTensor` behaviour): return a
        scalar::

            sage: b = M.linear_form(name='b')
            sage: b[:] = 7, 0, 2
            sage: a.__call__(b,v)
            12
            sage: a(b,v) == a.__call__(b,v)
            True
            sage: a(b,v) == s(b)
            True

        Identity map with a single argument: return a module element::

            sage: id = M.identity_map()
            sage: s = id.__call__(v) ; s
            Element v of the Rank-3 free module M over the Integer Ring
            sage: s == v
            True
            sage: s == id(v)
            True
            sage: s == id.contract(v)
            True

        Identity map with two arguments (:class:`FreeModuleTensor` behaviour):
        return a scalar::

            sage: id.__call__(b,v)
            22
            sage: id(b,v) == id.__call__(b,v)
            True
            sage: id(b,v) == b(v)
            True
        """
        from .free_module_element import FiniteRankFreeModuleElement
        if len(arg) > 1:
            # The automorphism acting as a type-(1,1) tensor on a pair
            # (linear form, module element), returning a scalar:
            if self._is_identity:
                if len(arg) != 2:
                    raise TypeError("wrong number of arguments")
                linform = arg[0]
                if linform._tensor_type != (0,1):
                    raise TypeError("the first argument must be a linear form")
                vector = arg[1]
                if not isinstance(vector, FiniteRankFreeModuleElement):
                    raise TypeError("the second argument must be a module" +
                                    " element")
                return linform(vector)
            else: # self is not the identity automorphism:
                return FreeModuleTensor.__call__(self, *arg)
        # The automorphism acting as such, on a module element, returning a
        # module element:
        vector = arg[0]
        if not isinstance(vector, FiniteRankFreeModuleElement):
            raise TypeError("the argument must be an element of a free module")
        if self._is_identity:
            return vector
        basis = self.common_basis(vector)
        t = self._components[basis]
        v = vector._components[basis]
        fmodule = self._fmodule
        result = vector._new_instance()
        for i in fmodule.irange():
            res = 0
            for j in fmodule.irange():
                res += t[[i,j]]*v[[j]]
            result.set_comp(basis)[i] = res
        # Name of the output:
        result._name = None
        if self._name is not None and vector._name is not None:
            result._name = self._name + "(" + vector._name + ")"
        # LaTeX symbol for the output:
        result._latex_name = None
        if self._latex_name is not None and vector._latex_name is not None:
            result._latex_name = self._latex_name + r"\left(" + \
                              vector._latex_name + r"\right)"
        return result

    #### End of FreeModuleTensor methods ####

    #### MultiplicativeGroupElement methods ####

    def __invert__(self):
        r"""
        Return the inverse automorphism.

        OUTPUT:

        - instance of :class:`FreeModuleAutomorphism` representing the
          automorphism that is the inverse of ``self``.

        EXAMPLES:

        Inverse of an automorphism of a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(name='a')
            sage: a[e,:] = [[1,0,0],[0,-1,2],[0,1,-3]]
            sage: a.inverse()
            Automorphism a^(-1) of the Rank-3 free module M over the Integer
             Ring
            sage: a.inverse().parent()
            General linear group of the Rank-3 free module M over the Integer
             Ring

        Check that ``a.inverse()`` is indeed the inverse automorphism::

            sage: a.inverse() * a
            Identity map of the Rank-3 free module M over the Integer Ring
            sage: a * a.inverse()
            Identity map of the Rank-3 free module M over the Integer Ring
            sage: a.inverse().inverse() == a
            True

        Another check is::

            sage: a.inverse().matrix(e)
            [ 1  0  0]
            [ 0 -3 -2]
            [ 0 -1 -1]
            sage: a.inverse().matrix(e) == (a.matrix(e))^(-1)
            True

        The inverse is cached (as long as ``a`` is not modified)::

            sage: a.inverse() is a.inverse()
            True

        If ``a`` is modified, the inverse is automatically recomputed::

            sage: a[0,0] = -1
            sage: a.matrix(e)
            [-1  0  0]
            [ 0 -1  2]
            [ 0  1 -3]
            sage: a.inverse().matrix(e)  # compare with above
            [-1  0  0]
            [ 0 -3 -2]
            [ 0 -1 -1]

        Shortcuts for :meth:`inverse` are the operator ``~`` and the exponent
        ``-1``::

            sage: ~a is a.inverse()
            True
            sage: (a)^(-1) is a.inverse()
            True

        The inverse of the identity map is of course itself::

            sage: id = M.identity_map()
            sage: id.inverse() is id
            True

        and we have::

            sage: a*(a)^(-1) == id
            True
            sage: (a)^(-1)*a == id
            True

        If the name could cause some confusion, a bracket is added around the
        element before taking the inverse::

            sage: c = M.automorphism(name='a^(-1)*b')
            sage: c[e,:] = [[1,0,0],[0,-1,1],[0,2,-1]]
            sage: c.inverse()
            Automorphism (a^(-1)*b)^(-1) of the Rank-3 free module M over the
             Integer Ring
        """
        from .comp import Components
        if self._is_identity:
            return self
        if self._inverse is None:
            from sage.tensor.modules.format_utilities import is_atomic
            if self._name is None:
                inv_name = None
            else:
                if is_atomic(self._name, ['*']):
                    inv_name = self._name + '^(-1)'
                else:
                    inv_name = '(' + self._name + ')^(-1)'
            if self._latex_name is None:
                inv_latex_name = None
            else:
                if is_atomic(self._latex_name, ['\\circ', '\\otimes']):
                    inv_latex_name = self._latex_name + r'^{-1}'
                else:
                    inv_latex_name = r'\left(' + self._latex_name + \
                                     r'\right)^{-1}'
            fmodule = self._fmodule
            si = fmodule._sindex
            nsi = fmodule._rank + si
            self._inverse = self.__class__(fmodule, inv_name, inv_latex_name)
            for basis in self._components:
                try:
                    mat = self.matrix(basis)
                except (KeyError, ValueError):
                    continue
                mat_inv = mat.inverse()
                cinv = Components(fmodule._ring, basis, 2, start_index=si,
                                  output_formatter=fmodule._output_formatter)
                for i in range(si, nsi):
                    for j in range(si, nsi):
                        cinv[i, j] = mat_inv[i-si,j-si]
                self._inverse._components[basis] = cinv
            self._inverse._inverse = self
        return self._inverse

    inverse = __invert__

    def _mul_(self, other):
        r"""
        Automorphism composition.

        This implements the group law of GL(M), M being the module of ``self``.

        INPUT:

        - ``other`` -- an automorphism of the same module as ``self``

        OUTPUT:

        - the automorphism resulting from the composition of ``other`` and
          ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,2],[1,3]])
            sage: b = M.automorphism([[3,4],[5,7]])
            sage: c = a._mul_(b) ; c
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: c.matrix()
            [13 18]
            [18 25]

        TESTS::

            sage: c.parent() is a.parent()
            True
            sage: c.matrix() == a.matrix() * b.matrix()
            True
            sage: c(e[0]) == a(b(e[0]))
            True
            sage: c(e[1]) == a(b(e[1]))
            True
            sage: a.inverse()._mul_(c) == b
            True
            sage: c._mul_(b.inverse()) == a
            True
            sage: id = M.identity_map()
            sage: id._mul_(a) == a
            True
            sage: a._mul_(id) == a
            True
            sage: a._mul_(a.inverse()) == id
            True
            sage: a.inverse()._mul_(a) == id
            True
        """
        # No need for consistency check since self and other are guaranteed
        # to have the same parent. In particular, they are defined on the same
        # free module.
        #
        # Special cases:
        if self._is_identity:
            return other
        if other._is_identity:
            return self
        if other is self._inverse or self is other._inverse:
            return self._fmodule.identity_map()
        # General case:
        fmodule = self._fmodule
        resu = self.__class__(fmodule)
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("no common basis for the composition")
        # The composition is performed as a tensor contraction of the last
        # index of self (position=1) and the first index of other (position=0):
        resu._components[basis] = self._components[basis].contract(1,
                                                    other._components[basis],0)
        return resu

    #### End of MultiplicativeGroupElement methods ####

    def __mul__(self, other):
        r"""
        Redefinition of
        :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.__mul__`
        so that * dispatches either to automorphism composition or to the
        tensor product.

        EXAMPLES:

        Automorphism composition::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,2],[1,3]])
            sage: b = M.automorphism([[3,4],[5,7]])
            sage: s = a*b ; s
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: s.matrix()
            [13 18]
            [18 25]
            sage: s.matrix() == a.matrix() * b.matrix()
            True
            sage: s(e[0]) == a(b(e[0]))
            True
            sage: s(e[1]) == a(b(e[1]))
            True
            sage: s.display()
            13 e_0⊗e^0 + 18 e_0⊗e^1 + 18 e_1⊗e^0 + 25 e_1⊗e^1

        Tensor product::

            sage: c =  M.tensor((1,1)) ; c
            Type-(1,1) tensor on the Rank-2 free module M over the Integer Ring
            sage: c[:] = [[3,4],[5,7]]
            sage: c[:] == b[:]  # c and b have the same components
            True
            sage: s = a*c ; s
            Type-(2,2) tensor on the Rank-2 free module M over the Integer Ring
            sage: s.display()
            3 e_0⊗e_0⊗e^0⊗e^0 + 4 e_0⊗e_0⊗e^0⊗e^1 + 6 e_0⊗e_0⊗e^1⊗e^0
             + 8 e_0⊗e_0⊗e^1⊗e^1 + 5 e_0⊗e_1⊗e^0⊗e^0 + 7 e_0⊗e_1⊗e^0⊗e^1
             + 10 e_0⊗e_1⊗e^1⊗e^0 + 14 e_0⊗e_1⊗e^1⊗e^1 + 3 e_1⊗e_0⊗e^0⊗e^0
             + 4 e_1⊗e_0⊗e^0⊗e^1 + 9 e_1⊗e_0⊗e^1⊗e^0 + 12 e_1⊗e_0⊗e^1⊗e^1
             + 5 e_1⊗e_1⊗e^0⊗e^0 + 7 e_1⊗e_1⊗e^0⊗e^1 + 15 e_1⊗e_1⊗e^1⊗e^0
             + 21 e_1⊗e_1⊗e^1⊗e^1

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,2],[1,3]], name='a')
            sage: b = M.automorphism([[0,1],[-1,0]], name='b')
            sage: mat_a0 = a.matrix(e)
            sage: a *= b
            sage: a.matrix(e) == mat_a0 * b.matrix(e)
            True
        """
        if isinstance(other, FreeModuleAutomorphism):
            return self._mul_(other)  # general linear group law
        else:
            return FreeModuleTensor.__mul__(self, other)  # tensor product

    def matrix(self, basis1=None, basis2=None):
        r"""
        Return the matrix of ``self`` w.r.t to a pair of bases.

        If the matrix is not known already, it is computed from the matrix in
        another pair of bases by means of the change-of-basis formula.

        INPUT:

        - ``basis1`` -- (default: ``None``) basis of the free module on which
          ``self`` is defined; if none is provided, the module's default basis
          is assumed
        - ``basis2`` -- (default: ``None``) basis of the free module on which
          ``self`` is defined; if none is provided, ``basis2`` is set to
          ``basis1``

        OUTPUT:

        - the matrix representing the automorphism ``self`` w.r.t
          to bases ``basis1`` and ``basis2``; more precisely, the columns of
          this matrix are formed by the components w.r.t. ``basis2`` of
          the images of the elements of ``basis1``.

        EXAMPLES:

        Matrices of an automorphism of a rank-3 free `\ZZ`-module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism([[-1,0,0],[0,1,2],[0,1,3]], name='a')
            sage: a.matrix(e)
            [-1  0  0]
            [ 0  1  2]
            [ 0  1  3]
            sage: a.matrix()
            [-1  0  0]
            [ 0  1  2]
            [ 0  1  3]
            sage: f = M.basis('f', from_family=(-e[2], 4*e[1]+3*e[3],  7*e[1]+5*e[3])) ; f
            Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring
            sage: a.matrix(f)
            [  1  -6 -10]
            [ -7  83 140]
            [  4 -48 -81]

        Check of the above matrix::

            sage: a(f[1]).display(f)
            a(f_1) = f_1 - 7 f_2 + 4 f_3
            sage: a(f[2]).display(f)
            a(f_2) = -6 f_1 + 83 f_2 - 48 f_3
            sage: a(f[3]).display(f)
            a(f_3) = -10 f_1 + 140 f_2 - 81 f_3

        Check of the change-of-basis formula::

            sage: P = M.change_of_basis(e,f).matrix(e)
            sage: a.matrix(f) == P^(-1) * a.matrix(e) * P
            True

        Check that the matrix of the product of two automorphisms is the
        product of their matrices::

            sage: b = M.change_of_basis(e,f) ; b
            Automorphism of the Rank-3 free module M over the Integer Ring
            sage: b.matrix(e)
            [ 0  4  7]
            [-1  0  0]
            [ 0  3  5]
            sage: (a*b).matrix(e) == a.matrix(e) * b.matrix(e)
            True

        Check that the matrix of the inverse automorphism is the inverse of the
        automorphism's matrix::

            sage: (~a).matrix(e)
            [-1  0  0]
            [ 0  3 -2]
            [ 0 -1  1]
            sage: (~a).matrix(e) == ~(a.matrix(e))
            True

        Matrices of the identity map::

            sage: id = M.identity_map()
            sage: id.matrix(e)
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: id.matrix(f)
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        from sage.matrix.constructor import matrix
        fmodule = self._fmodule
        if basis1 is None:
            basis1 = fmodule.default_basis()
        elif basis1 not in fmodule.bases():
            raise TypeError("{} is not a basis on the {}".format(basis1,
                                                                 fmodule))
        if basis2 is None:
            basis2 = basis1
        elif basis2 not in fmodule.bases():
            raise TypeError("{} is not a basis on the {}".format(basis2,
                                                                 fmodule))
        if (basis1, basis2) not in self._matrices:
            if basis2 == basis1:
                comp = self.components(basis1)
                mat = [[comp[[i,j]] for j in fmodule.irange()]
                                                     for i in fmodule.irange()]
                self._matrices[(basis1, basis1)] = matrix(mat)
            else:
                # 1/ determine the matrix w.r.t. basis1:
                self.matrix(basis1)
                # 2/ perform the change (basis1, basis1) --> (basis1, basis2):
                raise NotImplementedError("basis1 != basis2 not implemented yet")
        return self._matrices[(basis1, basis2)]

    def _some_matrix(self):
        r"""
        Return the matrix of ``self`` w.r.t. some basis.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,1],[0,2]], name='a')
            sage: a._some_matrix()
            [1 1]
            [0 2]
        """
        self.matrix()  # forces the update of the matrix in the module's default
                       # basis, to make sure that the dictionary self._matrices
                       # is not empty
        return next(iter(self._matrices.values()))

    @lazy_attribute
    def characteristic_polynomial(self):
        r"""
        Return the characteristic polynomial of ``self``.

        :meth:`characteristic_polynomial` and :meth:`charpoly` are the same method.

        INPUT:

        - ``var`` -- string (default: ``'x'``); a variable name

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,1],[0,2]], name='a')
            sage: a.matrix(e)
            [1 1]
            [0 2]
            sage: a.characteristic_polynomial()
            x^2 - 3*x + 2
            sage: a.charpoly()
            x^2 - 3*x + 2
            sage: a.charpoly('T')
            T^2 - 3*T + 2
        """
        return self._some_matrix().characteristic_polynomial

    charpoly = characteristic_polynomial

    @lazy_attribute
    def det(self):
        r"""
        Return the determinant of ``self``.

        OUTPUT:

        - element of the base ring of the module on which ``self`` is defined,
          equal to the determinant of ``self``.

        EXAMPLES:

        Determinant of an automorphism on a `\ZZ`-module of rank 2::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[4,7],[3,5]], name='a')
            sage: a.matrix(e)
            [4 7]
            [3 5]
            sage: a.det()
            -1
            sage: det(a)
            -1
            sage: ~a.det()  # determinant of the inverse automorphism
            -1
            sage: id = M.identity_map()
            sage: id.det()
            1
        """
        return self._some_matrix().det

    determinant = det

    @lazy_attribute
    def fcp(self):
        r"""
        Return the factorization of the characteristic polynomial of ``self``.

        INPUT:

        - ``var`` -- string (default: ``'x'``); a variable name

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,1],[0,2]], name='a')
            sage: a.matrix(e)
            [1 1]
            [0 2]
            sage: a.fcp()                                                               # needs sage.libs.pari
            (x - 2) * (x - 1)
            sage: a.fcp('T')                                                            # needs sage.libs.pari
            (T - 2) * (T - 1)
        """
        return self._some_matrix().fcp

    @lazy_attribute
    def minimal_polynomial(self):
        r"""
        Return the minimal polynomial of ``self``.

        :meth:`minimal_polynomial` and :meth:`minpoly` are the same method.

        INPUT:

        - ``var`` -- string (default: ``'x'``); a variable name

        EXAMPLES::

            sage: M = FiniteRankFreeModule(GF(7), 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[0,1,2], [-1,0,3], [2,4,1]], name='a')
            sage: a.minpoly()                                                           # needs sage.libs.pari
            x^3 + 6*x^2 + 6*x + 1
            sage: a.minimal_polynomial()                                                # needs sage.libs.pari
            x^3 + 6*x^2 + 6*x + 1
            sage: a.minimal_polynomial('T')                                             # needs sage.libs.pari
            T^3 + 6*T^2 + 6*T + 1
        """
        return self._some_matrix().minimal_polynomial

    minpoly = minimal_polynomial

    @lazy_attribute
    def trace(self):
        r"""
        Return the trace of ``self``.

        OUTPUT:

        - element of the base ring of the module on which ``self`` is defined,
          equal to the trace of ``self``.

        EXAMPLES:

        Trace of an automorphism on a `\ZZ`-module of rank 2::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[4,7],[3,5]], name='a')
            sage: a.matrix(e)
            [4 7]
            [3 5]
            sage: a.trace()
            9
            sage: id = M.identity_map()
            sage: id.trace()
            2
        """
        return self._some_matrix().trace
