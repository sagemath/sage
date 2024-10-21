# sage_setup: distribution = sagemath-modules
r"""
Elements of free modules of finite rank

The class :class:`FiniteRankFreeModuleElement` implements elements of
free modules of finite rank over a commutative ring.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Eric Gourgoulhon (2017): class :class:`FiniteRankFreeModuleElement` inherits
  from :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`

REFERENCES:

- Chap. 21 of R. Godement : *Algebra* [God1968]_
- Chap. 12 of J. M. Lee: *Introduction to Smooth Manifolds* [Lee2013]_ (only
  when the free module is a vector space)
- Chap. 2 of B. O'Neill: *Semi-Riemannian Geometry* [ONe1983]_
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

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from sage.tensor.modules.alternating_contr_tensor import AlternatingContrTensor
from sage.tensor.modules.comp import Components

if TYPE_CHECKING:
    from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
    from sage.tensor.modules.free_module_basis import FreeModuleBasis


class FiniteRankFreeModuleElement(AlternatingContrTensor):
    r"""
    Element of a free module of finite rank over a commutative ring.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`.

    The class :class:`FiniteRankFreeModuleElement` inherits from
    :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`
    because the elements of a free module `M` of finite rank over a commutative
    ring `R` are identified with tensors of type `(1,0)` on `M` via the
    canonical map

    .. MATH::

        \begin{array}{lllllll}
        \Phi: & M & \longrightarrow & M^{**}  &  &   & \\
              & v & \longmapsto & \bar v : & M^* & \longrightarrow & R \\
              &   &             &          & a   & \longmapsto & a(v)
        \end{array}

    Note that for free modules of finite rank, this map is actually an
    isomorphism, enabling the canonical identification: `M^{**}= M`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank over a commutative ring
      `R`, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``name`` -- (default: ``None``) name given to the element
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the element;
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Let us consider a rank-3 free module `M` over `\ZZ`::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e') ; e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring

    There are three ways to construct an element of the free module `M`:
    the first one (recommended) is using the free module::

        sage: v = M([2,0,-1], basis=e, name='v') ; v
        Element v of the Rank-3 free module M over the Integer Ring
        sage: v.display()  # expansion on the default basis (e)
        v = 2 e_0 - e_2
        sage: v.parent() is M
        True

    The second way is to construct a tensor of type `(1,0)` on `M` (cf. the
    canonical identification `M^{**} = M` recalled above)::

        sage: v2 = M.tensor((1,0), name='v')
        sage: v2[0], v2[2] = 2, -1 ; v2
        Element v of the Rank-3 free module M over the Integer Ring
        sage: v2.display()
        v = 2 e_0 - e_2
        sage: v2 == v
        True

    Finally, the third way is via some linear combination of the basis
    elements::

        sage: v3 = 2*e[0] - e[2]
        sage: v3.set_name('v') ; v3 # in this case, the name has to be set separately
        Element v of the Rank-3 free module M over the Integer Ring
        sage: v3.display()
        v = 2 e_0 - e_2
        sage: v3 == v
        True

    The canonical identification `M^{**} = M` is implemented by letting the
    module elements act on linear forms, providing the same result as the
    reverse operation (cf. the map `\Phi` defined above)::

        sage: a = M.linear_form(name='a')
        sage: a[:] = (2, 1, -3) ; a
        Linear form a on the Rank-3 free module M over the Integer Ring
        sage: v(a)
        7
        sage: a(v)
        7
        sage: a(v) == v(a)
        True

    .. RUBRIC:: ARITHMETIC EXAMPLES

    Addition::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e') ; e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: a = M([0,1,3], name='a') ; a
        Element a of the Rank-3 free module M over the Integer Ring
        sage: a.display()
        a = e_1 + 3 e_2
        sage: b = M([2,-2,1], name='b') ; b
        Element b of the Rank-3 free module M over the Integer Ring
        sage: b.display()
        b = 2 e_0 - 2 e_1 + e_2
        sage: s = a + b ; s
        Element a+b of the Rank-3 free module M over the Integer Ring
        sage: s.display()
        a+b = 2 e_0 - e_1 + 4 e_2
        sage: all(s[i] == a[i] + b[i] for i in M.irange())
        True

    Subtraction::

        sage: s = a - b ; s
        Element a-b of the Rank-3 free module M over the Integer Ring
        sage: s.display()
        a-b = -2 e_0 + 3 e_1 + 2 e_2
        sage: all(s[i] == a[i] - b[i] for i in M.irange())
        True

    Multiplication by a scalar::

        sage: s = 2*a ; s
        Element of the Rank-3 free module M over the Integer Ring
        sage: s.display()
        2 e_1 + 6 e_2
        sage: a.display()
        a = e_1 + 3 e_2

    Tensor product::

        sage: s = a*b ; s
        Type-(2,0) tensor a⊗b on the Rank-3 free module M over the Integer Ring
        sage: s.symmetries()
        no symmetry;  no antisymmetry
        sage: s[:]
        [ 0  0  0]
        [ 2 -2  1]
        [ 6 -6  3]
        sage: s = a*s ; s
        Type-(3,0) tensor a⊗a⊗b on the Rank-3 free module M over the Integer Ring
        sage: s[:]
        [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
         [[0, 0, 0], [2, -2, 1], [6, -6, 3]],
         [[0, 0, 0], [6, -6, 3], [18, -18, 9]]]

    Exterior product::

        sage: s = a.wedge(b) ; s
        Alternating contravariant tensor a∧b of degree 2 on the Rank-3 free
         module M over the Integer Ring
        sage: s.display()
        a∧b = -2 e_0∧e_1 - 6 e_0∧e_2 + 7 e_1∧e_2
    """

    def __init__(
        self,
        fmodule: FiniteRankFreeModule,
        name: Optional[str] = None,
        latex_name: Optional[str] = None,
    ):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: v = FiniteRankFreeModuleElement(M, name='v')
            sage: v[e,:] = (-2, 1, 3)
            sage: TestSuite(v).run(skip='_test_category') # see below

        In the above test suite, _test_category fails because v is not an
        instance of v.parent().category().element_class. Actually module
        elements must be constructed via FiniteRankFreeModule.element_class and
        not by a direct call to FiniteRankFreeModuleElement::

            sage: v1 = M.element_class(M, name='v')
            sage: v1[e,:] = (-2, 1, 3)
            sage: TestSuite(v1).run()
        """
        AlternatingContrTensor.__init__(self, fmodule, 1, name=name,
                                        latex_name=latex_name)

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M([1,-2,3], name='v')
            Element v of the Rank-3 free module M over the Integer Ring
        """
        description = "Element "
        if self._name is not None:
            description += self._name + " "
        description += "of the {}".format(self._fmodule)
        return description

    def _new_comp(self, basis: FreeModuleBasis) -> Components:
        r"""
        Create some (uninitialized) components of ``self`` in a given basis.

        This method, which is already implemented in
        :meth:`FreeModuleTensor._new_comp`, is redefined here for efficiency.

        INPUT:

        - ``basis`` -- basis of the free module on which ``self`` is defined

        OUTPUT: an instance of :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: v = M([1,-2,3], name='v')
            sage: v._new_comp(e)
            1-index components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring
            sage: type(v._new_comp(e))
            <class 'sage.tensor.modules.comp.Components'>
        """
        fmodule = self._fmodule  # the base free module
        return Components(fmodule._ring, basis, 1, start_index=fmodule._sindex,
                          output_formatter=fmodule._output_formatter)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: v = M([1,-2,3], name='v')
            sage: v._new_instance()
            Element of the Rank-3 free module M over the Integer Ring
            sage: v._new_instance().parent() is v.parent()
            True
        """
        return self.__class__(self._fmodule)
