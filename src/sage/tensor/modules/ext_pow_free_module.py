r"""
Exterior powers of dual free modules

Given a free module `M` of finite rank over a commutative ring `R`
and a positive integer `p`, the *p-th exterior power* of the dual of `M` is the
set `\Lambda^p(M^*)` of all alternating forms of degree `p` on `M`, i.e. of
all multilinear maps

.. MATH::

    \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
    \longrightarrow R

that vanish whenever any of two of their arguments are equal.
Note that `\Lambda^1(M^*) = M^*` (the dual of `M`).

`\Lambda^p(M^*)` is a free module of rank `\left({n\atop p}\right)` over `R`,
where `n` is the rank of `M`.
Accordingly, exterior powers of free modules are implemented by a class,
:class:`ExtPowerFreeModule`, which inherits from the class
:class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- K. Conrad: *Exterior powers*,
  `http://www.math.uconn.edu/~kconrad/blurbs/ <http://www.math.uconn.edu/~kconrad/blurbs/>`_
- Chap. 19 of S. Lang: *Algebra*, 3rd ed., Springer (New York) (2002)

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm

class ExtPowerFreeModule(FiniteRankFreeModule):
    r"""
    Class for the exterior powers of the dual of a free module of finite rank
    over a commutative ring.

    Given a free module `M` of finite rank over a commutative ring `R`
    and a positive integer `p`, the *p-th exterior power* of the dual of `M` is
    the set `\Lambda^p(M^*)` of all alternating forms of degree `p` on `M`,
    i.e. of all multilinear maps

    .. MATH::

        \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
        \longrightarrow R

    that vanish whenever any of two of their arguments are equal.
    Note that `\Lambda^1(M^*) = M^*` (the dual of `M`).

    `\Lambda^p(M^*)` is a free module of rank `\left({n\atop p}\right)` over
    `R`, where `n` is the rank of `M`.
    Accordingly, the class :class:`ExtPowerFreeModule` inherits from the class
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``degree`` -- positive integer; the degree `p` of the alternating forms
    - ``name`` -- (default: ``None``) string; name given to `\Lambda^p(M^*)`
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      `\Lambda^p(M^*)`

    EXAMPLES:

    2nd exterior power of the dual of a free `\ZZ`-module of rank 3::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
        sage: A = ExtPowerFreeModule(M, 2) ; A
        2nd exterior power of the dual of the Rank-3 free module M over the
         Integer Ring

    Instead of importing ExtPowerFreeModule in the global name space, it is
    recommended to use the module's method
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
        e^0/\e^1
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

    while non-zero elements are constructed by providing their components in a
    given basis::

        sage: e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: comp = [[0,3,-1],[-3,0,4],[1,-4,0]]
        sage: a = A(comp, basis=e, name='a') ; a
        Alternating form a of degree 2 on the Rank-3 free module M over the
         Integer Ring
        sage: a.display(e)
        a = 3 e^0/\e^1 - e^0/\e^2 + 4 e^1/\e^2

    An alternative is to construct the alternating form from an empty list of
    components and to set the nonzero components afterwards::

        sage: a = A([], name='a')
        sage: a.set_comp(e)[0,1] = 3
        sage: a.set_comp(e)[0,2] = -1
        sage: a.set_comp(e)[1,2] = 4
        sage: a.display(e)
        a = 3 e^0/\e^1 - e^0/\e^2 + 4 e^1/\e^2

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

    Since any tensor of type (0,1) is a linear form, there is a coercion map
    from the set `T^{(0,1)}(M)` of such tensors to `M^*`::

        sage: T01 = M.tensor_module(0,1) ; T01
        Free module of type-(0,1) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: M.dual().has_coerce_map_from(T01)
        True

    There is also a coercion map in the reverse direction::

        sage: T01.has_coerce_map_from(M.dual())
        True

    For a degree `p\geq 2`, the coercion holds only in the direction
    `\Lambda^p(M^*)\rightarrow T^{(0,p)}(M)`::

        sage: T02 = M.tensor_module(0,2) ; T02
        Free module of type-(0,2) tensors on the Rank-3 free module M over the
         Integer Ring
        sage: T02.has_coerce_map_from(A)
        True
        sage: A.has_coerce_map_from(T02)
        False

    The coercion map `T^{(0,1)}(M) \rightarrow M^*` in action::

        sage: b = T01([-2,1,4], basis=e, name='b') ; b
        Type-(0,1) tensor b on the Rank-3 free module M over the Integer Ring
        sage: b.display(e)
        b = -2 e^0 + e^1 + 4 e^2
        sage: lb = M.dual()(b) ; lb
        Linear form b on the Rank-3 free module M over the Integer Ring
        sage: lb.display(e)
        b = -2 e^0 + e^1 + 4 e^2

    The coercion map `M^* \rightarrow T^{(0,1)}(M)` in action::

        sage: tlb = T01(lb) ; tlb
        Type-(0,1) tensor b on the Rank-3 free module M over the Integer Ring
        sage: tlb == b
        True

    The coercion map `\Lambda^2(M^*)\rightarrow T^{(0,2)}(M)` in action::

        sage: ta = T02(a) ; ta
        Type-(0,2) tensor a on the Rank-3 free module M over the Integer Ring
        sage: ta.display(e)
        a = 3 e^0*e^1 - e^0*e^2 - 3 e^1*e^0 + 4 e^1*e^2 + e^2*e^0 - 4 e^2*e^1
        sage: a.display(e)
        a = 3 e^0/\e^1 - e^0/\e^2 + 4 e^1/\e^2
        sage: ta.symmetries() # the antisymmetry is of course preserved
        no symmetry;  antisymmetry: (0, 1)

    """

    Element = FreeModuleAltForm

    def __init__(self, fmodule, degree, name=None, latex_name=None):
        r"""
        TEST::

            sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: A = ExtPowerFreeModule(M, 2) ; A
            2nd exterior power of the dual of the Rank-3 free module M over
             the Integer Ring
            sage: TestSuite(A).run()

        """
        from sage.functions.other import binomial
        self._fmodule = fmodule
        self._degree = degree
        rank = binomial(fmodule._rank, degree)
        self._zero_element = 0 # provisory (to avoid infinite recursion in what
                               # follows)
        if degree == 1:  # case of the dual
            if name is None and fmodule._name is not None:
                name = fmodule._name + '*'
            if latex_name is None and fmodule._latex_name is not None:
                latex_name = fmodule._latex_name + r'^*'
        else:
            if name is None and fmodule._name is not None:
                name = '/\^{}('.format(degree) + fmodule._name + '*)'
            if latex_name is None and fmodule._latex_name is not None:
                latex_name = r'\Lambda^{' + str(degree) + r'}\left(' + \
                             fmodule._latex_name + r'^*\right)'
        FiniteRankFreeModule.__init__(self, fmodule._ring, rank, name=name,
                                      latex_name=latex_name,
                                      start_index=fmodule._sindex,
                                    output_formatter=fmodule._output_formatter)
        # Unique representation:
        if self._degree in self._fmodule._dual_exterior_powers:
            raise ValueError("the {}th exterior power of ".format(degree) +
                             "the dual of {}".format(self._fmodule) +
                             " has already been created")
        else:
            self._fmodule._dual_exterior_powers[self._degree] = self
        # Zero element
        self._zero_element = self._element_constructor_(name='zero',
                                                        latex_name='0')
        for basis in self._fmodule._known_bases:
            self._zero_element._components[basis] = \
                                            self._zero_element._new_comp(basis)
            # (since new components are initialized to zero)

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an alternating form.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.dual_exterior_power(1)
            sage: a = A._element_constructor_(0) ; a
            Linear form zero on the Rank-3 free module M over the Integer Ring
            sage: a = A._element_constructor_([2,0,-1], name='a') ; a
            Linear form a on the Rank-3 free module M over the Integer Ring
            sage: a.display()
            a = 2 e^0 - e^2
            sage: A = M.dual_exterior_power(2)
            sage: a = A._element_constructor_(0) ; a
            Alternating form zero of degree 2 on the Rank-3 free module M over
             the Integer Ring
            sage: a = A._element_constructor_([], name='a') ; a
            Alternating form a of degree 2 on the Rank-3 free module M over
             the Integer Ring
            sage: a[e,0,2], a[e,1,2] = 3, -1
            sage: a.display()
            a = 3 e^0/\e^2 - e^1/\e^2

        """
        if comp == 0:
            return self._zero_element
        if isinstance(comp, FreeModuleTensor):
            # coercion of a tensor of type (0,1) to a linear form
            tensor = comp # for readability
            if tensor.tensor_type() == (0,1) and self._degree == 1 and \
                                         tensor.base_module() is self._fmodule:
                resu = self.element_class(self._fmodule, 1, name=tensor._name,
                                          latex_name=tensor._latex_name)
                for basis, comp in tensor._components.iteritems():
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
        Construct some (unamed) alternating form.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 4, name='M')
            sage: e = M.basis('e')
            sage: a = M.dual_exterior_power(1)._an_element_() ; a
            Linear form on the 4-dimensional vector space M over the Rational
             Field
            sage: a.display()
            1/2 e^0
            sage: a = M.dual_exterior_power(2)._an_element_() ; a
            Alternating form of degree 2 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e^0/\e^1
            sage: a = M.dual_exterior_power(3)._an_element_() ; a
            Alternating form of degree 3 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e^0/\e^1/\e^2
            sage: a = M.dual_exterior_power(4)._an_element_() ; a
            Alternating form of degree 4 on the 4-dimensional vector space M
             over the Rational Field
            sage: a.display()
            1/2 e^0/\e^1/\e^2/\e^3

        """
        resu = self.element_class(self._fmodule, self._degree)
        if self._fmodule._def_basis is not None:
            sindex = self._fmodule._sindex
            ind = [sindex + i for i in range(resu._tensor_rank)]
            resu.set_comp()[ind] = self._fmodule._ring.an_element()
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        EXAMPLES:

        Sets of type-(0,1) tensors coerce to ``self`` if the degree is 1::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: L1 = M.dual_exterior_power(1) ; L1
            Dual of the Rank-3 free module M over the Integer Ring
            sage: T01 = M.tensor_module(0,1) ; T01
            Free module of type-(0,1) tensors on the Rank-3 free module M over
             the Integer Ring
            sage: L1._coerce_map_from_(T01)
            True

        Of course, coercions from other tensor types are meaningless::

            sage: L1._coerce_map_from_(M.tensor_module(1,0))
            False
            sage: L1._coerce_map_from_(M.tensor_module(0,2))
            False

        If the degree is larger than 1, there is no coercion::

            sage: L2 = M.dual_exterior_power(2) ; L2
            2nd exterior power of the dual of the Rank-3 free module M over
             the Integer Ring
            sage: L2._coerce_map_from_(M.tensor_module(0,2))
            False

        """
        from sage.tensor.modules.tensor_free_module import TensorFreeModule
        if isinstance(other, TensorFreeModule):
            # coercion of a type-(0,1) tensor to a linear form
            if self._fmodule is other._fmodule and self._degree == 1 and \
               other.tensor_type() == (0,1):
                return True
        return False

    #### End of parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: M.dual_exterior_power(1)._repr_()
            'Dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(2)._repr_()
            '2nd exterior power of the dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(3)._repr_()
            '3rd exterior power of the dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(4)._repr_()
            '4th exterior power of the dual of the Rank-5 free module M over the Integer Ring'
            sage: M.dual_exterior_power(5)._repr_()
            '5th exterior power of the dual of the Rank-5 free module M over the Integer Ring'

        """
        if self._degree == 1:
            return "Dual of the {}".format(self._fmodule)
        description = "{}".format(self._degree)
        if self._degree == 2:
            description += "nd"
        elif self._degree == 3:
            description += "rd"
        else:
            description += "th"
        description += " exterior power of the dual of the {}".format(
                                                                 self._fmodule)
        return description

    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module on which the exterior power is defined.

        EXAMPLE::

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
