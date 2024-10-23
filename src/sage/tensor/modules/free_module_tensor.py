r"""
Tensors on free modules

The class :class:`FreeModuleTensor` implements tensors on a free module `M`
of finite rank over a commutative ring. A *tensor of type* `(k,l)` on `M`
is a multilinear map:

.. MATH::

    \underbrace{M^*\times\cdots\times M^*}_{k\ \; \mbox{times}}
    \times \underbrace{M\times\cdots\times M}_{l\ \; \mbox{times}}
    \longrightarrow R

where `R` is the commutative ring over which the free module `M` is defined
and `M^* = \mathrm{Hom}_R(M,R)` is the dual of `M`. The integer `k + l` is
called the *tensor rank*. The set `T^{(k,l)}(M)` of tensors of type `(k,l)`
on `M` is a free module of finite rank over `R`, described by the
class :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule`.

Various derived classes of :class:`FreeModuleTensor` are devoted to specific
tensors:

* :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`
  for fully antisymmetric type-`(k, 0)` tensors *(alternating contravariant
  tensors)*;

  - :class:`~sage.tensor.modules.free_module_element.FiniteRankFreeModuleElement`
    for elements of `M`, considered as type-`(1,0)` tensors thanks to the
    canonical identification `M^{**}=M` (which holds since `M` is a free module
    of finite rank);

* :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm` for
  fully antisymmetric type-`(0, l)` tensors *(alternating forms)*;

* :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
  for type-`(1,1)` tensors representing invertible endomorphisms.

Each of these classes is a Sage *element* class, the corresponding *parent*
class being:

* :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule` for
  :class:`FreeModuleTensor`;
* :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
  for :class:`~sage.tensor.modules.free_module_element.FiniteRankFreeModuleElement`;
* :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule` for
  :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`;
* :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerDualFreeModule` for
  :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`;
* :class:`~sage.tensor.modules.free_module_linear_group.FreeModuleLinearGroup`
  for
  :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`.


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Michael Jung (2019): improve treatment of the zero element; add method
  ``copy_from``

REFERENCES:

- Chap. 21 of R. Godement : *Algebra* [God1968]_
- Chap. 12 of J. M. Lee: *Introduction to Smooth Manifolds* [Lee2013]_ (only
  when the free module is a vector space)
- Chap. 2 of B. O'Neill: *Semi-Riemannian Geometry* [ONe1983]_

EXAMPLES:

A tensor of type `(1, 1)` on a rank-3 free module over `\ZZ`::

    sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
    sage: t = M.tensor((1,1), name='t') ; t
    Type-(1,1) tensor t on the Rank-3 free module M over the Integer Ring
    sage: t.parent()
    Free module of type-(1,1) tensors on the Rank-3 free module M
     over the Integer Ring
    sage: t.parent() is M.tensor_module(1,1)
    True
    sage: t in M.tensor_module(1,1)
    True

Setting some component of the tensor in a given basis::

    sage: e = M.basis('e') ; e
    Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
    sage: t.set_comp(e)[0,0] = -3  # the component [0,0] w.r.t. basis e is set to -3

The unset components are assumed to be zero::

    sage: t.comp(e)[:]  # list of all components w.r.t. basis e
    [-3  0  0]
    [ 0  0  0]
    [ 0  0  0]
    sage: t.display(e)  # displays the expansion of t on the basis e_i⊗e^j of T^(1,1)(M)
    t = -3 e_0⊗e^0

The commands ``t.set_comp(e)`` and ``t.comp(e)`` can be abridged by providing
the basis as the first argument in the square brackets::

    sage: t[e,0,0] = -3
    sage: t[e,:]
    [-3  0  0]
    [ 0  0  0]
    [ 0  0  0]

Actually, since ``e`` is ``M``'s default basis, the mention of ``e``
can be omitted::

    sage: t[0,0] = -3
    sage: t[:]
    [-3  0  0]
    [ 0  0  0]
    [ 0  0  0]

For tensors of rank 2, the matrix of components w.r.t. a given basis is
obtained via the function ``matrix``::

    sage: matrix(t.comp(e))
    [-3  0  0]
    [ 0  0  0]
    [ 0  0  0]
    sage: matrix(t.comp(e)).parent()
    Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

Tensor components can be modified (reset) at any time::

    sage: t[0,0] = 0
    sage: t[:]
    [0 0 0]
    [0 0 0]
    [0 0 0]

Checking that ``t`` is zero::

    sage: t.is_zero()
    True
    sage: t == 0
    True
    sage: t == M.tensor_module(1,1).zero()  # the zero element of the module of all type-(1,1) tensors on M
    True

The components are managed by the class
:class:`~sage.tensor.modules.comp.Components`::

    sage: type(t.comp(e))
    <class 'sage.tensor.modules.comp.Components'>

Only nonzero components are actually stored, in the dictionary :attr:`_comp`
of class :class:`~sage.tensor.modules.comp.Components`, whose keys are
the indices::

    sage: t.comp(e)._comp
    {}
    sage: t.set_comp(e)[0,0] = -3 ; t.set_comp(e)[1,2] = 2
    sage: t.comp(e)._comp  # random output order (dictionary)
    {(0, 0): -3, (1, 2): 2}
    sage: t.display(e)
    t = -3 e_0⊗e^0 + 2 e_1⊗e^2

Further tests of the comparison operator::

    sage: t.is_zero()
    False
    sage: t == 0
    False
    sage: t == M.tensor_module(1,1).zero()
    False
    sage: t1 = t.copy()
    sage: t1 == t
    True
    sage: t1[2,0] = 4
    sage: t1 == t
    False

As a multilinear map `M^* \times M \rightarrow \ZZ`, the type-`(1,1)`
tensor ``t`` acts on pairs formed by a linear form and a module element::

    sage: a = M.linear_form(name='a') ; a[:] = (2, 1, -3) ; a
    Linear form a on the Rank-3 free module M over the Integer Ring
    sage: b = M([1,-6,2], name='b') ; b
    Element b of the Rank-3 free module M over the Integer Ring
    sage: t(a,b)
    -2
"""
# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

from sage.parallel.decorate import parallel
from sage.parallel.parallelism import Parallelism
from sage.rings.integer import Integer
from sage.structure.element import ModuleElementWithMutability
from sage.tensor.modules.comp import (
    get_components_class,
    Components,
    CompWithSym,
)
from sage.tensor.modules.tensor_with_indices import TensorWithIndices

if TYPE_CHECKING:
    from sage.symbolic.expression import Expression
    from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
    from sage.tensor.modules.free_module_basis import FreeModuleBasis
    from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
    from sage.manifolds.differentiable.poisson_tensor import PoissonTensorField
    from sage.manifolds.differentiable.symplectic_form import SymplecticForm


class FreeModuleTensor(ModuleElementWithMutability):
    r"""
    Tensor over a free module of finite rank over a commutative ring.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank over a commutative ring
      `R`, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``tensor_type`` -- pair ``(k, l)`` with ``k`` being the contravariant
      rank and ``l`` the covariant rank
    - ``name`` -- (default: ``None``) name given to the tensor
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the tensor;
      if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the tensor arguments: each symmetry is described by a tuple containing
      the positions of the involved arguments, with the convention
      ``position=0`` for the first argument. For instance:

      * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments;
      * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
        arguments and a symmetry between the 2nd, 4th and 5th arguments.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the arguments, with the same convention as for ``sym``
    - ``parent`` -- (default: ``None``) some specific parent (e.g. exterior
      power for alternating forms); if ``None``, ``fmodule.tensor_module(k,l)``
      is used

    EXAMPLES:

    A tensor of type `(1,1)` on a rank-3 free module over `\ZZ`::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: t = M.tensor((1,1), name='t') ; t
        Type-(1,1) tensor t on the Rank-3 free module M over the Integer Ring

    Tensors are *Element* objects whose parents are tensor free modules::

        sage: t.parent()
        Free module of type-(1,1) tensors on the
         Rank-3 free module M over the Integer Ring
        sage: t.parent() is M.tensor_module(1,1)
        True
    """
    _fmodule: FiniteRankFreeModule

    def __init__(
        self,
        fmodule: FiniteRankFreeModule,
        tensor_type,
        name: Optional[str] = None,
        latex_name: Optional[str] = None,
        sym=None,
        antisym=None,
        parent=None,
        implementation=None,
    ):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_tensor import FreeModuleTensor
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = FreeModuleTensor(M, (2,1), name='t', latex_name=r'\tau', sym=(0,1))
            sage: t[e,0,0,0] = -3
            sage: TestSuite(t).run(skip='_test_category') # see below

        In the above test suite, _test_category fails because t is not an
        instance of t.parent().category().element_class. Actually tensors
        must be constructed via TensorFreeModule.element_class and
        not by a direct call to FreeModuleTensor::

            sage: t1 = M.tensor_module(2,1).element_class(M, (2,1), name='t',
            ....:                                         latex_name=r'\tau',
            ....:                                         sym=(0,1))
            sage: t1[e,0,0,0] = -3
            sage: TestSuite(t1).run()
        """
        if parent is None:
            parent = fmodule.tensor_module(*tensor_type)
        ModuleElementWithMutability.__init__(self, parent)
        self._fmodule = fmodule
        self._tensor_type = tuple(tensor_type)
        self._tensor_rank = self._tensor_type[0] + self._tensor_type[1]
        self._is_zero = False # a priori, may be changed below or via
                              # method __bool__()
        self._name = name
        self._implementation = implementation
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._components: dict[FreeModuleBasis, Components] = {}  # dict. of the sets of components on various
                              # bases, with the bases as keys (initially empty)

        # Treatment of symmetry declarations:
        self._sym, self._antisym = CompWithSym._canonicalize_sym_antisym(
            self._tensor_rank, sym, antisym)

        # Initialization of derived quantities:
        FreeModuleTensor._init_derived(self)

    ####### Required methods for ModuleElement (beside arithmetic) #######

    def __bool__(self):
        r"""
        Return ``True`` if ``self`` is nonzero and ``False`` otherwise.

        This method is called by ``self.is_zero()``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((2,1))
            sage: t.add_comp(e)
            3-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring
            sage: bool(t)  # uninitialized components are zero
            False
            sage: t == 0
            True
            sage: t[e,1,0,2] = 4  # setting a nonzero component in basis e
            sage: t.display()
            4 e_1⊗e_0⊗e^2
            sage: bool(t)
            True
            sage: t == 0
            False
            sage: t[e,1,0,2] = 0
            sage: t.display()
            0
            sage: bool(t)
            False
            sage: t == 0
            True
        """
        basis = self.pick_a_basis()
        if not self._components[basis].is_zero():
            self._is_zero = False
            return True
        self._is_zero = True
        return False

    ##### End of required methods for ModuleElement (beside arithmetic) #####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: t
            Type-(2,1) tensor t on the Rank-3 free module M over the Integer Ring
        """
        # Special cases
        if self._tensor_type == (0,2) and self._sym == ((0,1),):
            description = "Symmetric bilinear form "
        else:
            # Generic case
            description = "Type-({},{}) tensor".format(
                            self._tensor_type[0], self._tensor_type[1])
        if self._name is not None:
            description += " " + self._name
        description += " on the {}".format(self._fmodule)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: t._latex_()
            't'
            sage: latex(t)
            t
            sage: t = M.tensor((2,1), name='t', latex_name=r'\tau')
            sage: t._latex_()
            '\\tau'
            sage: latex(t)
            \tau
            sage: t = M.tensor((2,1))  # unnamed tensor
            sage: t._latex_()
            '\\mbox{Type-(2,1) tensor on the Rank-3 free module M over the Integer Ring}'
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        return self._latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: t._init_derived()
        """
        pass # no derived quantities

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: t._del_derived()
        """
        pass # no derived quantities

    #### Simple accessors ####

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``.

        OUTPUT:

        - pair ``(k, l)``, where ``k`` is the contravariant rank and ``l``
          is the covariant rank

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: M.an_element().tensor_type()
            (1, 0)
            sage: t = M.tensor((2,1))
            sage: t.tensor_type()
            (2, 1)
        """
        return self._tensor_type

    def tensor_rank(self):
        r"""
        Return the tensor rank of ``self``.

        OUTPUT:

        - integer ``k+l``, where ``k`` is the contravariant rank and ``l``
          is the covariant rank

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: M.an_element().tensor_rank()
            1
            sage: t = M.tensor((2,1))
            sage: t.tensor_rank()
            3
        """
        return self._tensor_rank

    def base_module(self):
        r"""
        Return the module on which ``self`` is defined.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module on which the tensor is defined.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.an_element().base_module()
            Rank-3 free module M over the Integer Ring
            sage: t = M.tensor((2,1))
            sage: t.base_module()
            Rank-3 free module M over the Integer Ring
            sage: t.base_module() is M
            True
        """
        return self._fmodule

    def symmetries(self):
        r"""
        Print the list of symmetries and antisymmetries of ``self``.

        EXAMPLES:

        Various symmetries / antisymmetries for a rank-4 tensor::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((4,0), name='T') # no symmetry declared
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: t = M.tensor((4,0), name='T', sym=(0,1))
            sage: t.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: t = M.tensor((4,0), name='T', sym=[(0,1), (2,3)])
            sage: t.symmetries()
            symmetries: [(0, 1), (2, 3)];  no antisymmetry
            sage: t = M.tensor((4,0), name='T', sym=(0,1), antisym=(2,3))
            sage: t.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
        """
        if len(self._sym) == 0:
            s = "no symmetry; "
        elif len(self._sym) == 1:
            s = "symmetry: {}; ".format(self._sym[0])
        else:
            s = "symmetries: {}; ".format(list(self._sym))
        if len(self._antisym) == 0:
            a = "no antisymmetry"
        elif len(self._antisym) == 1:
            a = "antisymmetry: {}".format(self._antisym[0])
        else:
            a = "antisymmetries: {}".format(list(self._antisym))
        print(s+a)

    #### End of simple accessors #####

    def _preparse_display(self, basis=None, format_spec=None):
        r"""
        Helper function for display (to be used in the derived classes
        FreeModuleAltForm and AlternatingContrTensor as well)

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: f = M.basis('f')
            sage: v = M([1,-2])
            sage: v._preparse_display()
            (Basis (e_0,e_1) on the Rank-2 free module M over the Integer Ring, None)
            sage: v._preparse_display(f)
            (Basis (f_0,f_1) on the Rank-2 free module M over the Integer Ring, None)
            sage: v._preparse_display(e, 10)
            (Basis (e_0,e_1) on the Rank-2 free module M over the Integer Ring, 10)
            sage: v._preparse_display(format_spec=10)
            (Basis (e_0,e_1) on the Rank-2 free module M over the Integer Ring, 10)
        """
        if basis is None:
            basis = self._fmodule._def_basis
        return (basis, format_spec)

    def display(self, basis=None, format_spec=None):
        r"""
        Display ``self`` in terms of its expansion w.r.t. a given module basis.

        The expansion is actually performed onto tensor products of elements
        of the given basis and of elements of its dual basis (see examples
        below).
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``basis`` -- (default: ``None``) basis of the free module with
          respect to which the tensor is expanded; if none is provided,
          the module's default basis is assumed
        - ``format_spec`` -- (default: ``None``) format specification passed
          to ``self._fmodule._output_formatter`` to format the output

        EXAMPLES:

        Display of a module element (type-`(1,0)` tensor)::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M', start_index=1)
            sage: e = M.basis('e') ; e
            Basis (e_1,e_2) on the 2-dimensional vector space M over the
             Rational Field
            sage: v = M([1/3,-2], name='v')
            sage: v.display(e)
            v = 1/3 e_1 - 2 e_2
            sage: v.display()  # a shortcut since e is M's default basis
            v = 1/3 e_1 - 2 e_2
            sage: latex(v.display())  # display in the notebook
            v = \frac{1}{3} e_{1} -2 e_{2}

        A shortcut is ``disp()``::

            sage: v.disp()
            v = 1/3 e_1 - 2 e_2

        Display of a linear form (type-`(0,1)` tensor)::

            sage: de = e.dual_basis() ; de
            Dual basis (e^1,e^2) on the 2-dimensional vector space M over the
             Rational Field
            sage: w = - 3/4 * de[1] + de[2] ; w
            Linear form on the 2-dimensional vector space M over the Rational
             Field
            sage: w.set_name('w', latex_name='\\omega')
            sage: w.display()
            w = -3/4 e^1 + e^2
            sage: latex(w.display())  # display in the notebook
            \omega = -\frac{3}{4} e^{1} +e^{2}

        Display of a type-`(1,1)` tensor::

            sage: t = v*w ; t  # the type-(1,1) is formed as the tensor product of v by w
            Type-(1,1) tensor v⊗w on the 2-dimensional vector space M over the
             Rational Field
            sage: t.display()
            v⊗w = -1/4 e_1⊗e^1 + 1/3 e_1⊗e^2 + 3/2 e_2⊗e^1 - 2 e_2⊗e^2
            sage: latex(t.display())  # display in the notebook
            v\otimes \omega = -\frac{1}{4} e_{1}\otimes e^{1} +
             \frac{1}{3} e_{1}\otimes e^{2} + \frac{3}{2} e_{2}\otimes e^{1}
             -2 e_{2}\otimes e^{2}

        Display in a basis which is not the default one::

            sage: a = M.automorphism(matrix=[[1,2],[3,4]], basis=e)
            sage: f = e.new_basis(a, 'f')
            sage: v.display(f) # the components w.r.t basis f are first computed via the change-of-basis formula defined by a
            v = -8/3 f_1 + 3/2 f_2
            sage: w.display(f)
            w = 9/4 f^1 + 5/2 f^2
            sage: t.display(f)
            v⊗w = -6 f_1⊗f^1 - 20/3 f_1⊗f^2 + 27/8 f_2⊗f^1 + 15/4 f_2⊗f^2

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: t2 = v*w
            sage: t2.display(f)
            v⊗w = -6 f_1⊗f^1 - 20/3 f_1⊗f^2 + 27/8 f_2⊗f^1 + 15/4 f_2⊗f^2
            sage: t2[f,:] == t[f,:]  # check of the parallel computation
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        The output format can be set via the argument ``output_formatter``
        passed at the module construction::

            sage: N = FiniteRankFreeModule(QQ, 2, name='N', start_index=1,
            ....:                   output_formatter=Rational.numerical_approx)
            sage: e = N.basis('e')
            sage: v = N([1/3,-2], name='v')
            sage: v.display()  # default format (53 bits of precision)
            v = 0.333333333333333 e_1 - 2.00000000000000 e_2
            sage: latex(v.display())
            v = 0.333333333333333 e_{1} -2.00000000000000 e_{2}

        The output format is then controlled by the argument ``format_spec`` of
        the method :meth:`display`::

            sage: v.display(format_spec=10)  # 10 bits of precision
            v = 0.33 e_1 - 2.0 e_2

        Check that the bug reported in :issue:`22520` is fixed::

            sage: # needs sage.symbolic
            sage: M = FiniteRankFreeModule(SR, 3, name='M')
            sage: e = M.basis('e')
            sage: t = SR.var('t', domain='real')
            sage: (t*e[0]).display()
            t e_0
        """
        from sage.misc.latex import latex
        from sage.typeset.unicode_characters import unicode_otimes
        from .format_utilities import is_atomic, FormattedExpansion
        basis, format_spec = self._preparse_display(basis=basis,
                                                    format_spec=format_spec)
        cobasis = basis.dual_basis()
        comp = self.comp(basis)
        terms_txt = []
        terms_latex = []
        n_con = self._tensor_type[0]
        for ind in comp.index_generator():
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
                for k in range(n_con):
                    bases_txt.append(basis[ind[k]]._name)
                    bases_latex.append(latex(basis[ind[k]]))
                for k in range(n_con, self._tensor_rank):
                    bases_txt.append(cobasis[ind[k]]._name)
                    bases_latex.append(latex(cobasis[ind[k]]))
                basis_term_txt = unicode_otimes.join(bases_txt)
                basis_term_latex = r'\otimes '.join(bases_latex)
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
        if terms_txt == []:
            expansion_txt = '0'
        else:
            expansion_txt = terms_txt[0]
            for term in terms_txt[1:]:
                if term[0] == '-':
                    expansion_txt += ' - ' + term[1:]
                else:
                    expansion_txt += ' + ' + term
        if terms_latex == []:
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

    def display_comp(self, basis=None, format_spec=None, symbol=None,
                     latex_symbol=None, index_labels=None,
                     index_latex_labels=None, only_nonzero=True,
                     only_nonredundant=False):
        r"""
        Display the tensor components with respect to a given module
        basis, one per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``basis`` -- (default: ``None``) basis of the free module with
          respect to which the tensor components are defined; if ``None``,
          the module's default basis is assumed
        - ``format_spec`` -- (default: ``None``) format specification passed
          to ``self._fmodule._output_formatter`` to format the output
        - ``symbol`` -- (default: ``None``) string (typically a single letter)
          specifying the symbol for the components; if ``None``, the tensor
          name is used if it has been set, otherwise ``'X'`` is used
        - ``latex_symbol`` -- (default: ``None``) string specifying the LaTeX
          symbol for the components; if ``None``, the tensor LaTeX name
          is used if it has been set, otherwise ``'X'`` is used
        - ``index_labels`` -- (default: ``None``) list of strings representing
          the labels of each of the individual indices; if ``None``, integer
          labels are used
        - ``index_latex_labels`` -- (default: ``None``) list of strings
          representing the LaTeX labels of each of the individual indices; if
          ``None``, integers labels are used
        - ``only_nonzero`` -- boolean (default: ``True``); if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- boolean (default: ``False``); if ``True``,
          only nonredundant components are displayed in case of symmetries

        EXAMPLES:

        Display of the components of a type-`(2,1)` tensor on a rank 2
        vector space over `\QQ`::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: M = FiniteRankFreeModule(QQ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: t = M.tensor((2,1), name='T', sym=(0,1))
            sage: t[1,2,1], t[1,2,2], t[2,2,2] = 2/3, -1/4, 3
            sage: t.display()
            T = 2/3 e_1⊗e_2⊗e^1 - 1/4 e_1⊗e_2⊗e^2 + 2/3 e_2⊗e_1⊗e^1
                - 1/4 e_2⊗e_1⊗e^2 + 3 e_2⊗e_2⊗e^2
            sage: t.display_comp()
            T^12_1 = 2/3
            T^12_2 = -1/4
            T^21_1 = 2/3
            T^21_2 = -1/4
            T^22_2 = 3

        The LaTeX output for the notebook::

            sage: latex(t.display_comp())
            \begin{array}{lcl} {T}_{\phantom{\, 1}\phantom{\, 2}\,1}^{\,1\,2\phantom{\, 1}}
             & = & \frac{2}{3} \\ {T}_{\phantom{\, 1}\phantom{\, 2}\,2}^{\,1\,2\phantom{\, 2}}
             & = & -\frac{1}{4} \\ {T}_{\phantom{\, 2}\phantom{\, 1}\,1}^{\,2\,1\phantom{\, 1}}
             & = & \frac{2}{3} \\ {T}_{\phantom{\, 2}\phantom{\, 1}\,2}^{\,2\,1\phantom{\, 2}}
             & = & -\frac{1}{4} \\ {T}_{\phantom{\, 2}\phantom{\, 2}\,2}^{\,2\,2\phantom{\, 2}}
             & = & 3 \end{array}

        By default, only the non-vanishing components are displayed; to see
        all the components, the argument ``only_nonzero`` must be set to
        ``False``::

            sage: t.display_comp(only_nonzero=False)
            T^11_1 = 0
            T^11_2 = 0
            T^12_1 = 2/3
            T^12_2 = -1/4
            T^21_1 = 2/3
            T^21_2 = -1/4
            T^22_1 = 0
            T^22_2 = 3

        ``t`` being symmetric w.r.t. to its first two indices, one may ask to
        skip the components that can be deduced by symmetry::

            sage: t.display_comp(only_nonredundant=True)
            T^12_1 = 2/3
            T^12_2 = -1/4
            T^22_2 = 3

        The index symbols can be customized::

            sage: t.display_comp(index_labels=['x', 'y'])
            T^xy_x = 2/3
            T^xy_y = -1/4
            T^yx_x = 2/3
            T^yx_y = -1/4
            T^yy_y = 3

        Display of the components w.r.t. a basis different from the
        default one::

            sage: f = M.basis('f', from_family=(-e[1]+e[2], e[1]+e[2]))
            sage: t.display_comp(basis=f)
            T^11_1 = 29/24
            T^11_2 = 13/24
            T^12_1 = 3/4
            T^12_2 = 3/4
            T^21_1 = 3/4
            T^21_2 = 3/4
            T^22_1 = 7/24
            T^22_2 = 23/24
        """
        if basis is None:
            basis = self._fmodule._def_basis
        if symbol is None:
            if self._name is not None:
                symbol = self._name
            else:
                symbol = 'X'
        if latex_symbol is None:
            if self._latex_name is not None:
                latex_symbol = r'{' + self._latex_name + r'}'
            else:
                latex_symbol = 'X'
        index_positions = self._tensor_type[0]*'u' + self._tensor_type[1]*'d'
        return self.comp(basis).display(symbol,
                                        latex_symbol=latex_symbol,
                                        index_positions=index_positions,
                                        index_labels=index_labels,
                                        index_latex_labels=index_latex_labels,
                                        format_spec=format_spec,
                                        only_nonzero=only_nonzero,
                                        only_nonredundant=only_nonredundant)

    def set_name(self, name: Optional[str] = None, latex_name: Optional[str] = None):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``name`` -- (default: ``None``) string; name given to the tensor
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the tensor; if None while ``name`` is provided, the LaTeX symbol
          is set to ``name``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1)) ; t
            Type-(2,1) tensor on the Rank-3 free module M over the Integer Ring
            sage: t.set_name('t') ; t
            Type-(2,1) tensor t on the Rank-3 free module M over the Integer Ring
            sage: latex(t)
            t
            sage: t.set_name(latex_name=r'\tau') ; t
            Type-(2,1) tensor t on the Rank-3 free module M over the Integer Ring
            sage: latex(t)
            \tau
        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name

    def _new_instance(self):
        r"""
        Create a tensor of the same tensor type and with the same symmetries
        as ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: t._new_instance()
            Type-(2,1) tensor on the Rank-3 free module M over the Integer Ring
            sage: t._new_instance().parent() is t.parent()
            True
        """
        return self.__class__(self._fmodule, self._tensor_type, sym=self._sym,
                              antisym=self._antisym)

    def _new_comp(self, basis):
        r"""
        Create some (uninitialized) components of ``self`` w.r.t a given
        module basis.

        This method, to be called by :meth:`comp`, must be redefined by derived
        classes to adapt the output to the relevant subclass of
        :class:`~sage.tensor.modules.comp.Components`.

        INPUT:

        - ``basis`` -- basis of the free module on which ``self`` is defined

        OUTPUT:

        - an instance of :class:`~sage.tensor.modules.comp.Components`
          (or of one of its subclasses)

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: e = M.basis('e')
            sage: t._new_comp(e)
            3-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring
            sage: a = M.tensor((2,1), name='a', sym=(0,1))
            sage: a._new_comp(e)
            3-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring,
             with symmetry on the index positions (0, 1)
        """
        fmodule = self._fmodule  # the base free module
        if not self._sym and not self._antisym:
            components_class = get_components_class(implementation=self._implementation)
            return components_class(fmodule._ring, basis, self._tensor_rank,
                            start_index=fmodule._sindex, output_formatter=fmodule._output_formatter)
        for isym in self._sym:
            if len(isym) == self._tensor_rank:
                components_class = get_components_class(sym="fullysym", implementation=self._implementation)
                return components_class(fmodule._ring, basis, self._tensor_rank,
                            start_index=fmodule._sindex, output_formatter=fmodule._output_formatter)
        for isym in self._antisym:
            if len(isym) == self._tensor_rank:
                components_class = get_components_class(sym="fullyantisym", implementation=self._implementation)
                return components_class(fmodule._ring, basis, self._tensor_rank,
                            start_index=fmodule._sindex, output_formatter=fmodule._output_formatter)
        components_class = get_components_class(sym="sym", implementation=self._implementation)
        return components_class(fmodule._ring, basis, self._tensor_rank,
                           start_index=fmodule._sindex, output_formatter=fmodule._output_formatter,
                           sym=self._sym, antisym=self._antisym)

    def components(self, basis=None, from_basis=None) -> Components:
        r"""
        Return the components of ``self`` w.r.t to a given module basis.

        If the components are not known already, they are computed by the
        tensor change-of-basis formula from components in another basis.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          required; if none is provided, the components are assumed to refer
          to the module's default basis
        - ``from_basis`` -- (default: ``None``) basis from which the
          required components are computed, via the tensor change-of-basis
          formula, if they are not known already in the basis ``basis``;
          if none, a basis from which both the components and a change-of-basis
          to ``basis`` are known is selected.

        OUTPUT:

        - components in the basis ``basis``, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES:

        Components of a tensor of type-`(1,1)`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t[1,2] = -3 ; t[3,3] = 2
            sage: t.components()
            2-indices components w.r.t. Basis (e_1,e_2,e_3)
             on the Rank-3 free module M over the Integer Ring
            sage: t.components() is t.components(e)  # since e is M's default basis
            True
            sage: t.components()[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]

        A shortcut is ``t.comp()``::

            sage: t.comp() is t.components()
            True

        A direct access to the components w.r.t. the module's default basis is
        provided by the square brackets applied to the tensor itself::

            sage: t[1,2] is t.comp(e)[1,2]
            True
            sage: t[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]

        Components computed via a change-of-basis formula::

            sage: a = M.automorphism()
            sage: a[:] = [[0,0,1], [1,0,0], [0,-1,0]]
            sage: f = e.new_basis(a, 'f')
            sage: t.comp(f)
            2-indices components w.r.t. Basis (f_1,f_2,f_3)
             on the Rank-3 free module M over the Integer Ring
            sage: t.comp(f)[:]
            [ 0  0  0]
            [ 0  2  0]
            [-3  0  0]
        """
        fmodule = self._fmodule
        if basis is None:
            basis = fmodule._def_basis
        try:
            # Standard bases of tensor modules are keyed to the base module's basis,
            # not to the TensorFreeSubmoduleBasis_comp instance.
            basis = basis._base_module_basis
        except AttributeError:
            pass
        if basis not in self._components:
            # The components must be computed from
            # those in the basis from_basis
            if from_basis is None:
                for known_basis in self._components:
                    if (known_basis, basis) in self._fmodule._basis_changes \
                      and (basis, known_basis) in self._fmodule._basis_changes:
                        from_basis = known_basis
                        break
                if from_basis is None:
                    raise ValueError("no basis could be found for computing " +
                                     "the components in the {}".format(basis))
            elif from_basis not in self._components:
                raise ValueError("the tensor components are not known in " +
                                 "the {}".format(from_basis))
            (n_con, n_cov) = self._tensor_type
            pp = None
            if n_cov > 0:
                if (from_basis, basis) not in fmodule._basis_changes:
                    raise ValueError("the change-of-basis matrix from the " +
                                     "{} to the {}".format(from_basis, basis)
                                     + " has not been set")
                pp = \
                  fmodule._basis_changes[(from_basis, basis)].comp(from_basis)
                # pp not used if n_cov = 0 (pure contravariant tensor)
            ppinv = None
            if n_con > 0:
                if (basis, from_basis) not in fmodule._basis_changes:
                    raise ValueError("the change-of-basis matrix from the " +
                                     "{} to the {}".format(basis, from_basis) +
                                     " has not been set")
                ppinv = \
                  fmodule._basis_changes[(basis, from_basis)].comp(from_basis)
                # ppinv not used if n_con = 0 (pure covariant tensor)
            old_comp = self._components[from_basis]
            new_comp = self._new_comp(basis)
            rank = self._tensor_rank
            # loop on the new components:
            nproc = Parallelism().get('tensor')

            if nproc != 1:
                # Parallel computation
                lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
                ind_list = list(new_comp.non_redundant_index_generator())
                ind_step = max(1, int(len(ind_list)/nproc/2))
                local_list = lol(ind_list, ind_step)
                # list of input parameters
                listParalInput = [(old_comp, ppinv, pp, n_con, rank, ii) for ii in local_list]

                @parallel(p_iter='multiprocessing', ncpus=nproc)
                def paral_newcomp(old_comp, ppinv, pp, n_con, rank, local_list_ind):
                    partial = []
                    for ind in local_list_ind:
                        res = 0
                        # Summation on the old components multiplied by the proper
                        # change-of-basis matrix elements (tensor formula):
                        for ind_old in old_comp.index_generator():
                            t = old_comp[[ind_old]]
                            for i in range(n_con): # loop on contravariant indices
                                t *= ppinv[[ind[i], ind_old[i]]]
                            for i in range(n_con,rank):  # loop on covariant indices
                                t *= pp[[ind_old[i], ind[i]]]
                            res += t
                        partial.append([ind,res])
                    return partial

                for ii,val in paral_newcomp(listParalInput):
                    for jj in val:
                        new_comp[[jj[0]]] = jj[1]

            else:
                # Sequential computation
                for ind_new in new_comp.non_redundant_index_generator():
                    # Summation on the old components multiplied by the proper
                    # change-of-basis matrix elements (tensor formula):
                    res = 0
                    for ind_old in old_comp.index_generator():
                        t = old_comp[[ind_old]]
                        for i in range(n_con): # loop on contravariant indices
                            t *= ppinv[[ind_new[i], ind_old[i]]]
                        for i in range(n_con,rank):  # loop on covariant indices
                            t *= pp[[ind_old[i], ind_new[i]]]
                        res += t
                    new_comp[ind_new] = res

            self._components[basis] = new_comp
            # end of case where the computation was necessary
        return self._components[basis]

    comp = components

    def _set_comp_unsafe(self, basis=None):
        r"""
        Return the components of ``self`` w.r.t. a given module basis for
        assignment. This private method invokes no security check. Use
        this method at your own risk.

        The components with respect to other bases are deleted, in order to
        avoid any inconsistency. To keep them, use the method
        :meth:`_add_comp_unsafe` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        TESTS:

        Setting components of a type-`(1,1)` tensor::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t._set_comp_unsafe()[0,1] = -3
            sage: t.display()
            t = -3 e_0⊗e^1
            sage: t._set_comp_unsafe()[1,2] = 2
            sage: t.display()
            t = -3 e_0⊗e^1 + 2 e_1⊗e^2
            sage: t._set_comp_unsafe(e)
            2-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring

        Setting components in a new basis::

            sage: f =  M.basis('f')
            sage: t._set_comp_unsafe(f)[0,1] = 4
            sage: list(t._components) # the components w.r.t. basis e have been deleted
            [Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]
            sage: t.display(f)
            t = 4 f_0⊗f^1

        The components w.r.t. basis e can be deduced from those w.r.t. basis f,
        once a relation between the two bases has been set::

            sage: a = M.automorphism()
            sage: a[:] = [[0,0,1], [1,0,0], [0,-1,0]]
            sage: M.set_change_of_basis(e, f, a)
            sage: t.display(e)
            t = -4 e_1⊗e^2
            sage: sorted(t._components, key=repr)
            [Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]
        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis not in self._components:
            if basis not in self._fmodule._known_bases:
                raise ValueError("the {} has not been ".format(basis) +
                                 "defined on the {}".format(self._fmodule))
            self._components[basis] = self._new_comp(basis)
        self._del_derived() # deletes the derived quantities
        self.del_other_comp(basis)
        return self._components[basis]

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

        Setting components of a type-`(1,1)` tensor::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t.set_comp()[0,1] = -3
            sage: t.display()
            t = -3 e_0⊗e^1
            sage: t.set_comp()[1,2] = 2
            sage: t.display()
            t = -3 e_0⊗e^1 + 2 e_1⊗e^2
            sage: t.set_comp(e)
            2-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring

        Setting components in a new basis::

            sage: f =  M.basis('f')
            sage: t.set_comp(f)[0,1] = 4
            sage: list(t._components) # the components w.r.t. basis e have been deleted
            [Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]
            sage: t.display(f)
            t = 4 f_0⊗f^1

        The components w.r.t. basis e can be deduced from those w.r.t. basis f,
        once a relation between the two bases has been set::

            sage: a = M.automorphism()
            sage: a[:] = [[0,0,1], [1,0,0], [0,-1,0]]
            sage: M.set_change_of_basis(e, f, a)
            sage: t.display(e)
            t = -4 e_1⊗e^2
            sage: sorted(t._components, key=repr)
            [Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]

        Since zero is an immutable element, its components cannot be changed::

            sage: z = M.tensor_module(1, 1).zero()
            sage: z.set_comp(e)[0,1] = 1
            Traceback (most recent call last):
            ...
            ValueError: the components of an immutable element cannot be changed
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        self._is_zero = False  # a priori
        return self._set_comp_unsafe(basis)

    def _add_comp_unsafe(self, basis=None):
        r"""
        Return the components of ``self`` w.r.t. a given module basis for
        assignment, keeping the components w.r.t. other bases. This private
        method invokes no security check. Use this method at your own risk.

        To delete the components w.r.t. other bases, use the method
        :meth:`_set_comp_unsafe` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`;
          if such components did not exist previously, they are created

        TESTS:

        Setting components of a type-`(1,1)` tensor::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t._add_comp_unsafe()[0,1] = -3
            sage: t.display()
            t = -3 e_0⊗e^1
            sage: t._add_comp_unsafe()[1,2] = 2
            sage: t.display()
            t = -3 e_0⊗e^1 + 2 e_1⊗e^2
            sage: t._add_comp_unsafe(e)
            2-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring

        Adding components in a new basis::

            sage: f =  M.basis('f')
            sage: t._add_comp_unsafe(f)[0,1] = 4

        The components w.r.t. basis e have been kept::

            sage: sorted(t._components, key=repr)
            [Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]
            sage: t.display(f)
            t = 4 f_0⊗f^1
            sage: t.display(e)
            t = -3 e_0⊗e^1 + 2 e_1⊗e^2
        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis not in self._components:
            if basis not in self._fmodule._known_bases:
                raise ValueError("the {} has not been ".format(basis) +
                                 "defined on the {}".format(self._fmodule))
            self._components[basis] = self._new_comp(basis)
        self._del_derived() # deletes the derived quantities
        return self._components[basis]

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

            If the tensor has already components in other bases, it
            is the user's responsibility to make sure that the components
            to be added are consistent with them.

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`;
          if such components did not exist previously, they are created

        EXAMPLES:

        Setting components of a type-`(1,1)` tensor::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t.add_comp()[0,1] = -3
            sage: t.display()
            t = -3 e_0⊗e^1
            sage: t.add_comp()[1,2] = 2
            sage: t.display()
            t = -3 e_0⊗e^1 + 2 e_1⊗e^2
            sage: t.add_comp(e)
            2-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring

        Adding components in a new basis::

            sage: f =  M.basis('f')
            sage: t.add_comp(f)[0,1] = 4

        The components w.r.t. basis e have been kept::

            sage: sorted(t._components, key=repr)
            [Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]
            sage: t.display(f)
            t = 4 f_0⊗f^1
            sage: t.display(e)
            t = -3 e_0⊗e^1 + 2 e_1⊗e^2

        Since zero is an immutable element, its components cannot be changed::

            sage: z = M.tensor_module(1, 1).zero()
            sage: z.add_comp(e)[0,1] = 1
            Traceback (most recent call last):
            ...
            ValueError: the components of an immutable element cannot be changed
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        self._is_zero = False  # a priori
        return self._add_comp_unsafe(basis)

    def del_other_comp(self, basis=None):
        r"""
        Delete all the components but those corresponding to ``basis``.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          kept; if none the module's default basis is assumed

        EXAMPLES:

        Deleting components of a module element::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: u = M([2,1,-5])
            sage: f = M.basis('f')
            sage: u.add_comp(f)[:] = [0,4,2]
            sage: sorted(u._components, key=repr)
            [Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring,
             Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring]
            sage: u.del_other_comp(f)
            sage: list(u._components)
            [Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring]

        Let us restore the components w.r.t. e and delete those w.r.t. f::

            sage: u.add_comp(e)[:] = [2,1,-5]
            sage: sorted(u._components, key=repr)
            [Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring,
             Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring]
            sage: u.del_other_comp()  # default argument: basis = e
            sage: list(u._components)
            [Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring]
        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis not in self._components:
            raise ValueError(f"the components w.r.t. the {basis}"
                             " have not been defined")
        to_be_deleted = [other_basis for other_basis in self._components
                         if other_basis != basis]
        for other_basis in to_be_deleted:
            del self._components[other_basis]

    def __getitem__(self, args) -> Components:
        r"""
        Return a component w.r.t. some basis.

        NB: if ``args`` is a string, this method acts as a shortcut for
        tensor contractions and symmetrizations, the string containing
        abstract indices.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are returned. The basis can be passed
          as the first item of ``args``; if not, the free module's default
          basis is assumed.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: e = M.basis('e')
            sage: t.add_comp(e)
            3-indices components w.r.t. Basis (e_0,e_1,e_2) on the
             Rank-3 free module M over the Integer Ring
            sage: t.__getitem__((1,2,0)) # uninitialized components are zero
            0
            sage: t.__getitem__((e,1,2,0)) # same as above since e in the default basis
            0
            sage: t[1,2,0] = -4
            sage: t.__getitem__((e,1,2,0))
            -4
            sage: v = M([3,-5,2])
            sage: v.__getitem__(slice(None))
            [3, -5, 2]
            sage: v.__getitem__(slice(None)) == v[:]
            True
            sage: v.__getitem__((e, slice(None)))
            [3, -5, 2]
        """
        if isinstance(args, str): # tensor with specified indices
            return TensorWithIndices(self, args).update()
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], (int, Integer, slice)):
                basis = self._fmodule._def_basis
            else:
                basis = args[0]
                args = args[1:]
        else:
            if isinstance(args, (int, Integer, slice)):
                basis = self._fmodule._def_basis
            elif not isinstance(args[0], (int, Integer, slice)):
                basis = args[0]
                args = args[1:]
                if len(args) == 1:
                    args = args[0]  # to accommodate for [e,:] syntax
            else:
                basis = self._fmodule._def_basis
        return self.comp(basis)[args]

    def __setitem__(self, args, value):
        r"""
        Set a component w.r.t. some basis.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are set. The basis can be passed
          as the first item of ``args``; if not, the free module's default
          basis is assumed
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((0,2), name='t')
            sage: e = M.basis('e')
            sage: t.__setitem__((e,0,1), 5)
            sage: t.display()
            t = 5 e^0⊗e^1
            sage: t.__setitem__((0,1), 5)  # equivalent to above since e is the default basis
            sage: t.display()
            t = 5 e^0⊗e^1
            sage: t[0,1] = 5  # end-user usage
            sage: t.display()
            t = 5 e^0⊗e^1
            sage: t.__setitem__(slice(None), [[1,-2,3], [-4,5,-6], [7,-8,9]])
            sage: t[:]
            [ 1 -2  3]
            [-4  5 -6]
            [ 7 -8  9]
        """
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], (int, Integer, slice, tuple)):
                basis = self._fmodule._def_basis
            else:
                basis = args[0]
                args = args[1:]
        else:
            if isinstance(args, (int, Integer, slice)):
                basis = self._fmodule._def_basis
            elif not isinstance(args[0], (int, Integer, slice)):
                basis = args[0]
                args = args[1:]
                if len(args) == 1:
                    args = args[0]  # to accommodate for [e,:] syntax
            else:
                basis = self._fmodule._def_basis
        self.set_comp(basis)[args] = value

    def copy_from(self, other):
        r"""
        Make ``self`` to a copy from ``other``.

        INPUT:

        - ``other`` -- other tensor in the very same module from which
          ``self`` should be a copy of

        .. WARNING::

            All previous defined components will be deleted!

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t[1,2] = -3 ; t[3,3] = 2
            sage: s = M.tensor((1,1), name='s')
            sage: s.copy_from(t)
            sage: s[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]
            sage: s == t
            True

        If the original tensor is modified, the copy is not::

            sage: t[2,2] = 4
            sage: s[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]
            sage: s == t
            False
        """
        if self.is_immutable():
            raise ValueError("the components of an immutable element "
                             "cannot be changed")
        if other not in self.parent():
            raise TypeError("the original must be an element "
                            + "of {}".format(self.parent()))
        self._del_derived()
        self._components.clear()
        for basis, comp in other._components.items():
            self._components[basis] = comp.copy()
        self._is_zero = other._is_zero

    def copy(self, name=None, latex_name=None):
        r"""
        Return an exact copy of ``self``.

        The name and the derived quantities are not copied.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the copy
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          copy; if none is provided, the LaTeX symbol is set to ``name``

        EXAMPLES:

        Copy of a tensor of type `(1,1)`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t[1,2] = -3 ; t[3,3] = 2
            sage: t1 = t.copy()
            sage: t1[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]
            sage: t1 == t
            True

        If the original tensor is modified, the copy is not::

            sage: t[2,2] = 4
            sage: t1[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]
            sage: t1 == t
            False
        """
        resu = self._new_instance()
        resu.set_name(name=name, latex_name=latex_name)
        for basis, comp in self._components.items():
            resu._components[basis] = comp.copy()
        resu._is_zero = self._is_zero
        return resu

    def common_basis(self, other):
        r"""
        Find a common basis for the components of ``self`` and ``other``.

        In case of multiple common bases, the free module's default basis is
        privileged. If the current components of ``self`` and ``other``
        are all relative to different bases, a common basis is searched
        by performing a component transformation, via the transformations
        listed in ``self._fmodule._basis_changes``, still privileging
        transformations to the free module's default basis.

        INPUT:

        - ``other`` -- a tensor (instance of :class:`FreeModuleTensor`)

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing the common basis; if no common basis is found, ``None``
          is returned

        EXAMPLES:

        Common basis for the components of two module elements::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: u = M([2,1,-5])
            sage: f = M.basis('f')
            sage: M._basis_changes.clear() # to ensure that bases e and f are unrelated at this stage
            sage: v = M([0,4,2], basis=f)
            sage: u.common_basis(v)

        The above result is ``None`` since ``u`` and ``v`` have been defined
        on different bases and no connection between these bases have
        been set::

            sage: list(u._components)
            [Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring]
            sage: list(v._components)
            [Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring]

        Linking bases ``e`` and ``f`` changes the result::

            sage: a = M.automorphism()
            sage: a[:] = [[0,0,1], [1,0,0], [0,-1,0]]
            sage: M.set_change_of_basis(e, f, a)
            sage: u.common_basis(v)
            Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring

        Indeed, v is now known in basis e::

            sage: sorted(v._components, key=repr)
            [Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring,
             Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring]
        """
        # Compatibility checks:
        if not isinstance(other, FreeModuleTensor):
            raise TypeError("the argument must be a tensor on a free module")
        fmodule = self._fmodule
        if other._fmodule != fmodule:
            raise TypeError("the two tensors are not defined on the same " +
                            "free module")
        def_basis = fmodule._def_basis

        # 1/ Search for a common basis among the existing components, i.e.
        #    without performing any component transformation.
        #    -------------------------------------------------------------
        if def_basis in self._components and def_basis in other._components:
            return def_basis # the module's default basis is privileged
        for basis1 in self._components:
            if basis1 in other._components:
                return basis1

        # 2/ Search for a common basis via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least
        # one component transformation to get a common basis
        if def_basis in self._components:
            for obasis in other._components:
                if (obasis, def_basis) in fmodule._basis_changes:
                    other.comp(def_basis, from_basis=obasis)
                    return def_basis
        if def_basis in other._components:
            for sbasis in self._components:
                if (sbasis, def_basis) in fmodule._basis_changes:
                    self.comp(def_basis, from_basis=sbasis)
                    return def_basis
        # If this point is reached, then def_basis cannot be a common basis
        # via a single component transformation
        for sbasis in self._components:
            for obasis in other._components:
                if (obasis, sbasis) in fmodule._basis_changes:
                    other.comp(sbasis, from_basis=obasis)
                    return sbasis
                if (sbasis, obasis) in fmodule._basis_changes:
                    self.comp(obasis, from_basis=sbasis)
                    return obasis

        # 3/ Search for a common basis via two component transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at two
        # component transformation to get a common basis
        for sbasis in self._components:
            for obasis in other._components:
                if (sbasis, def_basis) in fmodule._basis_changes and \
                   (obasis, def_basis) in fmodule._basis_changes:
                    self.comp(def_basis, from_basis=sbasis)
                    other.comp(def_basis, from_basis=obasis)
                    return def_basis
                for basis in fmodule._known_bases:
                    if (sbasis, basis) in fmodule._basis_changes and \
                       (obasis, basis) in fmodule._basis_changes:
                        self.comp(basis, from_basis=sbasis)
                        other.comp(basis, from_basis=obasis)
                        return basis

        # If this point is reached, no common basis could be found, even at
        # the price of component transformations:
        return None

    def pick_a_basis(self):
        r"""
        Return a basis in which the tensor components are defined.

        The free module's default basis is privileged.

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing the basis

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,0), name='t')
            sage: e = M.basis('e')
            sage: t[0,1] = 4  # component set in the default basis (e)
            sage: t.pick_a_basis()
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: f = M.basis('f')
            sage: t.add_comp(f)[2,1] = -4  # the components in basis e are not erased
            sage: t.pick_a_basis()
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: t.set_comp(f)[2,1] = -4  # the components in basis e not erased
            sage: t.pick_a_basis()
            Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring
        """
        if self._fmodule._def_basis in self._components:
            return self._fmodule._def_basis  # the default basis is privileged
        else:
            # a basis is picked arbitrarily:
            return next(iter(self._components.items()))[0]

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a tensor or 0

        OUTPUT: ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,0), name='t')
            sage: e = M.basis('e')
            sage: t[0,1] = 7
            sage: t.__eq__(0)
            False
            sage: t[0,1] = 0
            sage: t.__eq__(0)
            True
            sage: a = M.tensor((0,2), name='a')
            sage: a[0,1] = 7
            sage: t[0,1] = 7
            sage: a[:], t[:]
            (
            [0 7 0]  [0 7 0]
            [0 0 0]  [0 0 0]
            [0 0 0], [0 0 0]
            )
            sage: t.__eq__(a)  # False since t and a do not have the same tensor type
            False
            sage: a = M.tensor((2,0), name='a') # same tensor type as t
            sage: a[0,1] = 7
            sage: t.__eq__(a)
            True
        """
        if self is other:
            return True

        if self._tensor_rank == 0:
            raise NotImplementedError("scalar comparison not implemented")
        if isinstance(other, (int, Integer)): # other should be 0
            if other == 0:
                return self.is_zero()
            else:
                return False
        elif not isinstance(other, FreeModuleTensor):
            return False
        else: # other is another tensor
            if other._fmodule != self._fmodule:
                return False
            if other._tensor_type != self._tensor_type:
                return False
            basis = self.common_basis(other)
            if basis is None:
                raise ValueError("no common basis for the comparison")
            return bool(self._components[basis] == other._components[basis])

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- a tensor or 0

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,0), name='t')
            sage: e = M.basis('e')
            sage: t[0,1] = 7
            sage: t.__ne__(0)
            True
            sage: t[0,1] = 0
            sage: t.__ne__(0)
            False
            sage: a = M.tensor((2,0), name='a') # same tensor type as t
            sage: a[0,1] = 7
            sage: t.__ne__(a)
            True
            sage: t[0,1] = 7
            sage: t.__ne__(a)
            False
        """
        return not self == other

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT: an exact copy of ``self``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,0), name='t')
            sage: e = M.basis('e')
            sage: t[0,1] = 7
            sage: p = t.__pos__() ; p
            Type-(2,0) tensor +t on the Rank-3 free module M over the Integer Ring
            sage: p.display()
            +t = 7 e_0⊗e_1
            sage: p == t
            True
            sage: p is t
            False
        """
        result = self._new_instance()
        for basis in self._components:
            result._components[basis] = + self._components[basis]
        if self._name is not None:
            result._name = '+' + self._name
        if self._latex_name is not None:
            result._latex_name = '+' + self._latex_name
        return result

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT: the tensor `-T`, where `T` is ``self``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((2,0), name='t')
            sage: e = M.basis('e')
            sage: t[0,1], t[1,2] = 7, -4
            sage: t.display()
            t = 7 e_0⊗e_1 - 4 e_1⊗e_2
            sage: a = t.__neg__() ; a
            Type-(2,0) tensor -t on the Rank-3 free module M over the Integer Ring
            sage: a.display()
            -t = -7 e_0⊗e_1 + 4 e_1⊗e_2
            sage: a == -t
            True
        """
        result = self._new_instance()
        for basis in self._components:
            result._components[basis] = - self._components[basis]
        if self._name is not None:
            result._name = '-' + self._name
        if self._latex_name is not None:
            result._latex_name = '-' + self._latex_name
        return result

    ######### ModuleElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Tensor addition.

        INPUT:

        - ``other`` -- a tensor, of the same type as ``self``

        OUTPUT: the tensor resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[4,0], [-2,5]]
            sage: b = M.tensor((2,0), name='b')
            sage: b[:] = [[0,1], [2,3]]
            sage: s = a._add_(b) ; s
            Type-(2,0) tensor a+b on the Rank-2 free module M over the Integer Ring
            sage: s[:]
            [4 1]
            [0 8]
            sage: a._add_(-a) == 0
            True
            sage: a._add_(a) == 2*a
            True
        """
        # No need for consistency check since self and other are guaranteed
        # to belong to the same tensor module
        #
        # Case zero:
        if self._is_zero:
            return other
        if other._is_zero:
            return self
        # Generic case:
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("no common basis for the addition")
        comp_result = self._components[basis] + other._components[basis]
        result = self._fmodule.tensor_from_comp(self._tensor_type, comp_result)
        if self._name is not None and other._name is not None:
            result._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '+' + other._latex_name
        return result

    def _sub_(self, other):
        r"""
        Tensor subtraction.

        INPUT:

        - ``other`` -- a tensor, of the same type as ``self``

        OUTPUT: the tensor resulting from the subtraction of ``other`` from ``self``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[4,0], [-2,5]]
            sage: b = M.tensor((2,0), name='b')
            sage: b[:] = [[0,1], [2,3]]
            sage: s = a._sub_(b) ; s
            Type-(2,0) tensor a-b on the Rank-2 free module M over the Integer Ring
            sage: s[:]
            [ 4 -1]
            [-4  2]
            sage: b._sub_(a) == -s
            True
            sage: a._sub_(a) == 0
            True
            sage: a._sub_(-a) == 2*a
            True

        TESTS:

        Check for when there is not a basis, but the same object::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: t == t
            True
        """
        # No need for consistency check since self and other are guaranteed
        # to belong to the same tensor module
        #
        # Case zero:
        if self._is_zero:
            return -other
        if other._is_zero:
            return self
        # Generic case:
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("no common basis for the subtraction")
        comp_result = self._components[basis] - other._components[basis]
        result = self._fmodule.tensor_from_comp(self._tensor_type, comp_result)
        if self._name is not None and other._name is not None:
            result._name = self._name + '-' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '-' + other._latex_name
        return result

    def _rmul_(self, other):
        r"""
        Multiplication on the left by ``other``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[4,0], [-2,5]]
            sage: s = a._rmul_(2) ; s
            Type-(2,0) tensor on the Rank-2 free module M over the Integer Ring
            sage: s[:]
            [ 8  0]
            [-4 10]
            sage: s == a + a
            True
            sage: a._rmul_(0)
            Type-(2,0) tensor on the Rank-2 free module M over the Integer Ring
            sage: a._rmul_(0) == 0
            True
            sage: a._rmul_(1) == a
            True
            sage: a._rmul_(-1) == -a
            True
        """
        #!# The following test is probably not necessary:
        if isinstance(other, FreeModuleTensor):
            raise NotImplementedError("left tensor product not implemented")
        # Left multiplication by a scalar:
        result = self._new_instance()
        for basis in self._components:
            result._components[basis] = other * self._components[basis]
        # If other has a name, set the name of the result:
        try:
            from .format_utilities import format_mul_txt, format_mul_latex
            result_name = format_mul_txt(other._name, '*', self._name)
            result_latex = format_mul_latex(other._latex_name, r' \cdot ',
                                            self._latex_name)
            result.set_name(name=result_name, latex_name=result_latex)
        except AttributeError:
            pass
        return result

    ######### End of ModuleElement arithmetic operators ########

    def __mul__(self, other):
        r"""
        Tensor product.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[4,0], [-2,5]]
            sage: b = M.tensor((0,2), name='b', antisym=(0,1))
            sage: b[0,1] = 3
            sage: s = a.__mul__(b) ; s
            Type-(2,2) tensor a⊗b on the Rank-2 free module M over the Integer Ring
            sage: s.symmetries()
            no symmetry;  antisymmetry: (2, 3)
            sage: s[:]
            [[[[0, 12], [-12, 0]], [[0, 0], [0, 0]]],
             [[[0, -6], [6, 0]], [[0, 15], [-15, 0]]]]
        """
        from sage.typeset.unicode_characters import unicode_otimes
        from .format_utilities import format_mul_txt, format_mul_latex
        if isinstance(other, FreeModuleTensor):
            basis = self.common_basis(other)
            if basis is None:
                raise ValueError("no common basis for the tensor product")
            comp_prov = self._components[basis] * other._components[basis]
            # Reordering of the contravariant and covariant indices:
            k1, l1 = self._tensor_type
            k2, l2 = other._tensor_type
            if l1 != 0:
                comp_result = comp_prov.swap_adjacent_indices(k1,
                                                          self._tensor_rank,
                                                          self._tensor_rank+k2)
            else:
                comp_result = comp_prov  # no reordering is necessary
            result = self._fmodule.tensor_from_comp((k1+k2, l1+l2),
                                                    comp_result)
            result._name = format_mul_txt(self._name, unicode_otimes, other._name)
            result._latex_name = format_mul_latex(self._latex_name,
                                                r'\otimes ', other._latex_name)
            return result

        # multiplication by a scalar:
        return FreeModuleTensor._rmul_(self, other)

    def __truediv__(self, other):
        r"""
        Division (by a scalar).

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[4,0], [-2,5]]
            sage: s = a.__truediv__(4) ; s
            Type-(2,0) tensor on the 2-dimensional vector space M over the
             Rational Field
            sage: s[:]
            [   1    0]
            [-1/2  5/4]
            sage: 4*s == a
            True
            sage: s == a/4
            True
        """
        result = self._new_instance()
        for basis in self._components:
            result._components[basis] = self._components[basis] / other
        return result

    def __call__(self, *args) -> Expression:
        r"""
        The tensor acting on linear forms and module elements as a multilinear
        map.

        INPUT:

        - ``*args`` -- list of `k` linear forms and `l` module elements
          with ``self`` being a tensor of type `(k, l)`

        EXAMPLES:

        Action of a type-`(2,1)` tensor::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((2,1), name='t', antisym=(0,1))
            sage: t[0,1,0], t[0,1,1] = 3, 2
            sage: t.display()
            t = 3 e_0⊗e_1⊗e^0 + 2 e_0⊗e_1⊗e^1 - 3 e_1⊗e_0⊗e^0 - 2 e_1⊗e_0⊗e^1
            sage: a = M.linear_form()
            sage: a[:] = 1, 2
            sage: b = M.linear_form()
            sage: b[:] = 3, -1
            sage: v = M([-2,1])
            sage: t.__call__(a,b,v)
            28
            sage: t(a,b,v) == t.__call__(a,b,v)
            True
            sage: t(a,b,v) == t.contract(v).contract(b).contract(a)
            True

        Action of a linear form on a vector::

            sage: a.__call__(v)
            0
            sage: a.__call__(v) == a(v)
            True
            sage: a(v) == a.contract(v)
            True
            sage: b.__call__(v)
            -7
            sage: b.__call__(v) == b(v)
            True
            sage: b(v) == b.contract(v)
            True

        Action of a vector on a linear form::

            sage: v.__call__(a)
            0
            sage: v.__call__(b)
            -7
        """
        # Consistency checks:
        p = len(args)
        if p != self._tensor_rank:
            raise TypeError(str(self._tensor_rank) +
                            " arguments must be provided")
        for i in range(self._tensor_type[0]):
            if not isinstance(args[i], FreeModuleTensor):
                raise TypeError("the argument no. " + str(i+1) +
                                " must be a linear form")
            if args[i]._tensor_type != (0,1):
                raise TypeError("the argument no. " + str(i+1) +
                                " must be a linear form")
        for i in range(self._tensor_type[0], p):
            if not isinstance(args[i], FreeModuleTensor):
                raise TypeError("the argument no. " + str(i+1) +
                                " must be a module element")
            if args[i]._tensor_type != (1,0):
                raise TypeError("the argument no. " + str(i+1) +
                                " must be a module element")
        fmodule = self._fmodule
        #
        # Specific case of a linear form acting on a vector (for efficiency):
        #
        if self._tensor_type == (0,1):
            vector = args[0]
            basis = self.common_basis(vector)
            if basis is None:
                raise ValueError("no common basis for the components")
            omega = self._components[basis]
            vv = vector._components[basis]
            resu = 0
            for i in fmodule.irange():
                resu += omega[[i]]*vv[[i]]
            # Name and LaTeX symbol of the output:
            if hasattr(resu, '_name'):
                if self._name is not None and vector._name is not None:
                    resu._name = self._name + "(" + vector._name + ")"
            if hasattr(resu, '_latex_name'):
                if self._latex_name is not None and \
                                                vector._latex_name is not None:
                    resu._latex_name = self._latex_name + r"\left(" + \
                                       vector._latex_name + r"\right)"
            return resu
        #
        # Generic case
        #
        # Search for a common basis
        basis = None
        # First try with the module's default basis
        def_basis = fmodule._def_basis
        if def_basis in self._components:
            basis = def_basis
            for arg in args:
                if def_basis not in arg._components:
                    basis = None
                    break
        if basis is None:
            # Search for another basis:
            for bas in self._components:
                basis = bas
                for arg in args:
                    if bas not in arg._components:
                        basis = None
                        break
                if basis is not None: # common basis found !
                    break
        if basis is None:
            # A last attempt to find a common basis, possibly via a
            # change-of-components transformation
            for arg in args:
                self.common_basis(arg) # to trigger some change of components
            for bas in self._components:
                basis = bas
                for arg in args:
                    if bas not in arg._components:
                        basis = None
                        break
                if basis is not None: # common basis found !
                    break
        if basis is None:
            raise ValueError("no common basis for the components")
        t = self._components[basis]
        v = [args[i]._components[basis] for i in range(p)]
        res = 0
        for ind in t.index_generator():
            prod = t[[ind]]
            for i in range(p):
                prod *= v[i][[ind[i]]]
            res += prod
        # Name of the output:
        if hasattr(res, '_name'):
            res_name = None
            if self._name is not None:
                res_name = self._name + "("
                for i in range(p-1):
                    if args[i]._name is not None:
                        res_name += args[i]._name + ","
                    else:
                        res_name = None
                        break
                if res_name is not None:
                    if args[p-1]._name is not None:
                        res_name += args[p-1]._name + ")"
                    else:
                        res_name = None
            res._name = res_name
        # LaTeX symbol of the output:
        if hasattr(res, '_latex_name'):
            res_latex = None
            if self._latex_name is not None:
                res_latex = self._latex_name + r"\left("
                for i in range(p-1):
                    if args[i]._latex_name is not None:
                        res_latex += args[i]._latex_name + ","
                    else:
                        res_latex = None
                        break
                if res_latex is not None:
                    if args[p-1]._latex_name is not None:
                        res_latex += args[p-1]._latex_name + r"\right)"
                    else:
                        res_latex = None
            res._latex_name = res_latex
        return res

    def trace(
        self,
        pos1: int = 0,
        pos2: int = 1,
        using: Optional[
            Union[PseudoRiemannianMetric, SymplecticForm, PoissonTensorField]
        ] = None,
    ):
        r"""
        Trace (contraction) on two slots of the tensor.

        If a non-degenerate form is provided, the trace of a type-`(0,2)` tensor
        is computed by first raising the last index.

        INPUT:

        - ``pos1`` -- (default: 0) position of the first index for the
          contraction, with the convention ``pos1=0`` for the first slot

        - ``pos2`` -- (default: 1) position of the second index for the
          contraction, with the same convention as for ``pos1``; the variance
          type of ``pos2`` must be opposite to that of ``pos1``

        - ``using`` -- (default: ``None``) a non-degenerate form

        OUTPUT:

        - tensor or scalar resulting from the ``(pos1, pos2)`` contraction

        EXAMPLES:

        Trace of a type-`(1,1)` tensor::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e') ; e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: a = M.tensor((1,1), name='a') ; a
            Type-(1,1) tensor a on the Rank-3 free module M over the Integer Ring
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: a.trace()
            15
            sage: a.trace(0,1)  # equivalent to above (contraction of slot 0 with slot 1)
            15
            sage: a.trace(1,0)  # the order of the slots does not matter
            15

        Instead of the explicit call to the method :meth:`trace`, one
        may use the index notation with Einstein convention (summation over
        repeated indices); it suffices to pass the indices as a string inside
        square brackets::

            sage: a['^i_i']
            15

        The letter 'i' to denote the repeated index can be replaced by any
        other letter::

            sage: a['^s_s']
            15

        Moreover, the symbol ``^`` can be omitted::

            sage: a['i_i']
            15

        The contraction on two slots having the same tensor type cannot occur::

            sage: b = M.tensor((2,0), name='b') ; b
            Type-(2,0) tensor b on the Rank-3 free module M over the Integer Ring
            sage: b[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b.trace(0,1)
            Traceback (most recent call last):
            ...
            IndexError: contraction on two contravariant indices is not allowed

        The contraction either preserves or destroys the symmetries::

            sage: b = M.alternating_form(2, 'b') ; b
            Alternating form b of degree 2 on the Rank-3 free module M
             over the Integer Ring
            sage: b[0,1], b[0,2], b[1,2] = 3, 2, 1
            sage: t = a*b ; t
            Type-(1,3) tensor a⊗b on the Rank-3 free module M
             over the Integer Ring

        By construction, ``t`` is a tensor field antisymmetric w.r.t. its
        last two slots::

            sage: t.symmetries()
            no symmetry;  antisymmetry: (2, 3)
            sage: s = t.trace(0,1) ; s   # contraction on the first two slots
            Alternating form of degree 2 on the
             Rank-3 free module M over the Integer Ring
            sage: s.symmetries()    # the antisymmetry is preserved
            no symmetry;  antisymmetry: (0, 1)
            sage: s[:]
            [  0  45  30]
            [-45   0  15]
            [-30 -15   0]
            sage: s == 15*b  # check
            True
            sage: s = t.trace(0,2) ; s   # contraction on the first and third slots
            Type-(0,2) tensor on the Rank-3 free module M over the Integer Ring
            sage: s.symmetries()  # the antisymmetry has been destroyed by the above contraction:
            no symmetry;  no antisymmetry
            sage: s[:]  # indeed:
            [-26  -4   6]
            [-31  -2   9]
            [-36   0  12]
            sage: s[:] == matrix( [[sum(t[k,i,k,j] for k in M.irange())
            ....:          for j in M.irange()] for i in M.irange()] )  # check
            True

        Use of index notation instead of :meth:`trace`::

            sage: t['^k_kij'] == t.trace(0,1)
            True
            sage: t['^k_{kij}'] == t.trace(0,1) # LaTeX notation
            True
            sage: t['^k_ikj'] == t.trace(0,2)
            True
            sage: t['^k_ijk'] == t.trace(0,3)
            True

        Index symbols not involved in the contraction may be replaced by
        dots::

            sage: t['^k_k..'] == t.trace(0,1)
            True
            sage: t['^k_.k.'] == t.trace(0,2)
            True
            sage: t['^k_..k'] == t.trace(0,3)
            True
        """
        if using is not None:
            if self.tensor_type() != (0, 2):
                raise ValueError(
                    "trace with respect to a non-degenerate form is only defined for type-(0,2) tensor"
                )
            return self.up(using, 1).trace()

        # The indices at pos1 and pos2 must be of different types:
        k_con = self._tensor_type[0]
        l_cov = self._tensor_type[1]
        if pos1 < k_con and pos2 < k_con:
            raise IndexError("contraction on two contravariant indices is " +
                             "not allowed")
        if pos1 >= k_con and pos2 >= k_con:
            raise IndexError("contraction on two covariant indices is " +
                             "not allowed")
        # Frame selection for the computation:
        if self._fmodule._def_basis in self._components:
            basis = self._fmodule._def_basis
        else: # a basis is picked arbitrarily:
            basis = self.pick_a_basis()
        resu_comp = self._components[basis].trace(pos1, pos2)
        if self._tensor_rank == 2:  # result is a scalar
            return resu_comp
        else:
            return self._fmodule.tensor_from_comp((k_con-1, l_cov-1),
                                                  resu_comp)

    def contract(self, *args):
        r"""
        Contraction on one or more indices with another tensor.

        INPUT:

        - ``pos1`` -- positions of the indices in ``self`` involved in the
          contraction; ``pos1`` must be a sequence of integers, with 0 standing
          for the first index position, 1 for the second one, etc; if ``pos1``
          is not provided, a single contraction on the last index position of
          ``self`` is assumed
        - ``other`` -- the tensor to contract with
        - ``pos2`` -- positions of the indices in ``other`` involved in the
          contraction, with the same conventions as for ``pos1``; if ``pos2``
          is not provided, a single contraction on the first index position of
          ``other`` is assumed

        OUTPUT:

        - tensor resulting from the contraction at the positions ``pos1`` and
          ``pos2`` of ``self`` with ``other``

        EXAMPLES:

        Contraction of a tensor of type `(0,1)` with a tensor of type `(1,0)`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.linear_form()  # tensor of type (0,1) is a linear form
            sage: a[:] = [-3,2,1]
            sage: b = M([2,5,-2])  # tensor of type (1,0) is a module element
            sage: s = a.contract(b) ; s
            2
            sage: s in M.base_ring()
            True
            sage: s == a[0]*b[0] + a[1]*b[1] + a[2]*b[2]  # check of the computation
            True

        The positions of the contraction indices can be set explicitly::

            sage: s == a.contract(0, b, 0)
            True
            sage: s == a.contract(0, b)
            True
            sage: s == a.contract(b, 0)
            True

        Instead of the explicit call to the method :meth:`contract`, the index
        notation can be used to specify the contraction, via Einstein
        convention (summation on repeated indices); it suffices to pass the
        indices as a string inside square brackets::

            sage: s1 = a['_i']*b['^i'] ; s1
            2
            sage: s1 == s
            True

        In the present case, performing the contraction is identical to
        applying the linear form to the module element::

            sage: a.contract(b) == a(b)
            True

        or to applying the module element, considered as a tensor of type
        `(1,0)`, to the linear form::

            sage: a.contract(b) == b(a)
            True

        We have also::

            sage: a.contract(b) == b.contract(a)
            True

        Contraction of a tensor of type `(1,1)` with a tensor of type `(1,0)`::

            sage: a = M.tensor((1,1))
            sage: a[:] = [[-1,2,3],[4,-5,6],[7,8,9]]
            sage: s = a.contract(b) ; s
            Element of the Rank-3 free module M over the Integer Ring
            sage: s.display()
            2 e_0 - 29 e_1 + 36 e_2

        Since the index positions have not been specified, the contraction
        takes place on the last position of a (i.e. no. 1) and the first
        position of ``b`` (i.e. no. 0)::

            sage: a.contract(b) == a.contract(1, b, 0)
            True
            sage: a.contract(b) == b.contract(0, a, 1)
            True
            sage: a.contract(b) == b.contract(a, 1)
            True

        Using the index notation with Einstein convention::

            sage: a['^i_j']*b['^j'] == a.contract(b)
            True

        The index ``i`` can be replaced by a dot::

            sage: a['^._j']*b['^j'] == a.contract(b)
            True

        and the symbol ``^`` may be omitted, the distinction between
        contravariant and covariant indices being the position with respect to
        the symbol ``_``::

            sage: a['._j']*b['j'] == a.contract(b)
            True

        Contraction is possible only between a contravariant index and a
        covariant one::

            sage: a.contract(0, b)
            Traceback (most recent call last):
            ...
            TypeError: contraction on two contravariant indices not permitted

        Contraction of a tensor of type `(2,1)` with a tensor of type `(0,2)`::

            sage: a = a*b ; a
            Type-(2,1) tensor on the Rank-3 free module M over the Integer Ring
            sage: b = M.tensor((0,2))
            sage: b[:] = [[-2,3,1], [0,-2,3], [4,-7,6]]
            sage: s = a.contract(1, b, 1) ; s
            Type-(1,2) tensor on the Rank-3 free module M over the Integer Ring
            sage: s[:]
            [[[-9, 16, 39], [18, -32, -78], [27, -48, -117]],
             [[36, -64, -156], [-45, 80, 195], [54, -96, -234]],
             [[63, -112, -273], [72, -128, -312], [81, -144, -351]]]

        Check of the computation::

            sage: all(s[i,j,k] == a[i,0,j]*b[k,0]+a[i,1,j]*b[k,1]+a[i,2,j]*b[k,2]
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True

        Using index notation::

            sage: a['il_j']*b['_kl'] == a.contract(1, b, 1)
            True

        LaTeX notation are allowed::

            sage: a['^{il}_j']*b['_{kl}'] == a.contract(1, b, 1)
            True

        Indices not involved in the contraction may be replaced by dots::

            sage: a['.l_.']*b['_.l'] == a.contract(1, b, 1)
            True

        The two tensors do not have to be defined on the same basis for the
        contraction to take place, reflecting the fact that the contraction is
        basis-independent::

            sage: A = M.automorphism()
            sage: A[:] =  [[0,0,1], [1,0,0], [0,-1,0]]
            sage: h = e.new_basis(A, 'h')
            sage: b.comp(h)[:]  # forces the computation of b's components w.r.t. basis h
            [-2 -3  0]
            [ 7  6 -4]
            [ 3 -1 -2]
            sage: b.del_other_comp(h)  # deletes components w.r.t. basis e
            sage: list(b._components)  # indeed:
            [Basis (h_0,h_1,h_2) on the Rank-3 free module M over the Integer Ring]
            sage: list(a._components)  # while a is known only in basis e:
            [Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring]
            sage: s1 = a.contract(1, b, 1) ; s1  # yet the computation is possible
            Type-(1,2) tensor on the Rank-3 free module M over the Integer Ring
            sage: s1 == s  # ... and yields the same result as previously:
            True

        The contraction can be performed on more than a single index; for
        instance a `2`-indices contraction of a type-`(2,1)` tensor with a
        type-`(1,2)` one is::

            sage: a  # a is a tensor of type-(2,1)
            Type-(2,1) tensor on the Rank-3 free module M over the Integer Ring
            sage: b = M([1,-1,2])*b ; b # a tensor of type (1,2)
            Type-(1,2) tensor on the Rank-3 free module M over the Integer Ring
            sage: s = a.contract(1,2,b,1,0) ; s # the double contraction
            Type-(1,1) tensor on the Rank-3 free module M over the Integer Ring
            sage: s[:]
            [ -36   30   15]
            [-252  210  105]
            [-204  170   85]
            sage: s == a['^.k_l']*b['^l_k.']  # the same thing in index notation
            True
        """
        #
        # Treatment of the input
        #
        nargs = len(args)
        for i, arg in enumerate(args):
            if isinstance(arg, FreeModuleTensor):
                other = arg
                it = i
                break
        else:
            raise TypeError("a tensor must be provided in the argument list")
        if it == 0:
            pos1 = (self._tensor_rank - 1,)
        else:
            pos1 = args[:it]
        if it == nargs-1:
            pos2 = (0,)
        else:
            pos2 = args[it+1:]
        ncontr = len(pos1) # number of contractions
        if len(pos2) != ncontr:
            raise TypeError("different number of indices for the contraction")
        k1, l1 = self._tensor_type
        k2, l2 = other._tensor_type
        for i in range(ncontr):
            p1 = pos1[i]
            p2 = pos2[i]
            if p1 < k1 and p2 < k2:
                raise TypeError("contraction on two contravariant indices " +
                                "not permitted")
            if p1 >= k1 and p2 >= k2:
                raise TypeError("contraction on two covariant indices " +
                                "not permitted")
        #
        # Contraction at the component level
        #
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("no common basis for the contraction")
        args = pos1 + (other._components[basis],) + pos2
        cmp_res = self._components[basis].contract(*args)
        if self._tensor_rank + other._tensor_rank - 2*ncontr == 0:
            # Case of scalar output:
            return cmp_res
        #
        # Reordering of the indices to have all contravariant indices first:
        #
        nb_cov_s = 0  # Number of covariant indices of self not involved in the
                      # contraction
        for pos in range(k1,k1+l1):
            if pos not in pos1:
                nb_cov_s += 1
        nb_con_o = 0  # Number of contravariant indices of other not involved
                      # in the contraction
        for pos in range(0,k2):
            if pos not in pos2:
                nb_con_o += 1
        if nb_cov_s != 0 and nb_con_o != 0:
            # some reordering is necessary:
            p2 = k1 + l1 - ncontr
            p1 = p2 - nb_cov_s
            p3 = p2 + nb_con_o
            cmp_res = cmp_res.swap_adjacent_indices(p1, p2, p3)
        type_res = (k1+k2-ncontr, l1+l2-ncontr)
        return self._fmodule.tensor_from_comp(type_res, cmp_res)

    def symmetrize(self, *pos, **kwargs):
        r"""
        Symmetrization over some arguments.

        INPUT:

        - ``pos`` -- list of argument positions involved in the
          symmetrization (with the convention ``position=0`` for the first
          argument); if none, the symmetrization is performed over all the
          arguments
        - ``basis`` -- (default: ``None``) module basis with respect to which
          the component computation is to be performed; if none, the module's
          default basis is used if the tensor field has already components
          in it; otherwise another basis w.r.t. which the tensor has
          components will be picked

        OUTPUT:

        - the symmetrized tensor (instance of :class:`FreeModuleTensor`)

        EXAMPLES:

        Symmetrization of a tensor of type `(2,0)`::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((2,0))
            sage: t[:] = [[2,1,-3],[0,-4,5],[-1,4,2]]
            sage: s = t.symmetrize() ; s
            Type-(2,0) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: t[:], s[:]
            (
            [ 2  1 -3]  [  2 1/2  -2]
            [ 0 -4  5]  [1/2  -4 9/2]
            [-1  4  2], [ -2 9/2   2]
            )
            sage: s.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: all(s[i,j] == 1/2*(t[i,j]+t[j,i])   # check:
            ....:     for i in range(3) for j in range(3))
            True

        Instead of invoking the method :meth:`symmetrize`, one may use the
        index notation with parentheses to denote the symmetrization; it
        suffices to pass the indices as a string inside square brackets::

            sage: t['(ij)']
            Type-(2,0) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: t['(ij)'].symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: t['(ij)'] == t.symmetrize()
            True

        The indices names are not significant; they can even be replaced by
        dots::

            sage: t['(..)'] == t.symmetrize()
            True

        The LaTeX notation can be used as well::

            sage: t['^{(ij)}'] == t.symmetrize()
            True

        Symmetrization of a tensor of type `(0,3)` on the first two arguments::

            sage: t = M.tensor((0,3))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]],
            ....:         [[10,-11,12], [13,14,-15], [16,17,18]],
            ....:         [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.symmetrize(0,1) ; s  # (0,1) = the first two arguments
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: s[:]
            [[[1, 2, 3], [3, -3, 9], [13, -6, -15]],
             [[3, -3, 9], [13, 14, -15], [-3, 20, 21]],
             [[13, -6, -15], [-3, 20, 21], [25, 26, -27]]]
            sage: all(s[i,j,k] == 1/2*(t[i,j,k]+t[j,i,k])   # Check:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.symmetrize(0,1) == s  # another test
            True

        Again the index notation can be used::

            sage: t['_(ij)k'] == t.symmetrize(0,1)
            True
            sage: t['_(..).'] == t.symmetrize(0,1)  # no index name
            True
            sage: t['_{(ij)k}'] == t.symmetrize(0,1)  # LaTeX notation
            True
            sage: t['_{(..).}'] == t.symmetrize(0,1)  # this also allowed
            True

        Symmetrization of a tensor of type `(0,3)` on the first and
        last arguments::

            sage: s = t.symmetrize(0,2) ; s  # (0,2) = first and last arguments
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s.symmetries()
            symmetry: (0, 2);  no antisymmetry
            sage: s[:]
            [[[1, 6, 11], [-4, 9, -8], [7, 12, 8]],
             [[6, -11, -4], [9, 14, 4], [12, 17, 22]],
             [[11, -4, -21], [-8, 4, 24], [8, 22, -27]]]
            sage: all(s[i,j,k] == 1/2*(t[i,j,k]+t[k,j,i])
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.symmetrize(0,2) == s  # another test
            True

        Symmetrization of a tensor of type `(0,3)` on the last two arguments::

            sage: s = t.symmetrize(1,2) ; s  # (1,2) = the last two arguments
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s.symmetries()
            symmetry: (1, 2);  no antisymmetry
            sage: s[:]
            [[[1, -1, 5], [-1, 5, 7], [5, 7, -9]],
             [[10, 1, 14], [1, 14, 1], [14, 1, 18]],
             [[19, -21, 2], [-21, 23, 25], [2, 25, -27]]]
            sage: all(s[i,j,k] == 1/2*(t[i,j,k]+t[i,k,j])   # Check:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.symmetrize(1,2) == s  # another test
            True

        Use of the index notation::

            sage: t['_i(jk)'] == t.symmetrize(1,2)
            True
            sage: t['_.(..)'] == t.symmetrize(1,2)
            True
            sage: t['_{i(jk)}'] == t.symmetrize(1,2)  # LaTeX notation
            True

        Full symmetrization of a tensor of type `(0,3)`::

            sage: s = t.symmetrize() ; s
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s.symmetries()
            symmetry: (0, 1, 2);  no antisymmetry
            sage: s[:]
            [[[1, 8/3, 29/3], [8/3, 7/3, 0], [29/3, 0, -5/3]],
             [[8/3, 7/3, 0], [7/3, 14, 25/3], [0, 25/3, 68/3]],
             [[29/3, 0, -5/3], [0, 25/3, 68/3], [-5/3, 68/3, -27]]]
            sage: all(s[i,j,k] == 1/6*(t[i,j,k]+t[i,k,j]+t[j,k,i]+t[j,i,k]+t[k,i,j]+t[k,j,i])  # Check:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.symmetrize() == s  # another test
            True

        Index notation for the full symmetrization::

            sage: t['_(ijk)'] == t.symmetrize()
            True
            sage: t['_{(ijk)}'] == t.symmetrize()  # LaTeX notation
            True

        Symmetrization can be performed only on arguments on the same type::

            sage: t = M.tensor((1,2))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]],
            ....:         [[10,-11,12], [13,14,-15], [16,17,18]],
            ....:         [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.symmetrize(0,1)
            Traceback (most recent call last):
            ...
            TypeError: 0 is a contravariant position, while 1 is a covariant position;
            symmetrization is meaningful only on tensor arguments of the same type
            sage: s = t.symmetrize(1,2) # OK: both 1 and 2 are covariant positions

        The order of positions does not matter::

            sage: t.symmetrize(2,1) == t.symmetrize(1,2)
            True

        Use of the index notation::

            sage: t['^i_(jk)'] == t.symmetrize(1,2)
            True
            sage: t['^._(..)'] ==  t.symmetrize(1,2)
            True

        The character ``^`` can be skipped, the character ``_`` being
        sufficient to separate contravariant indices from covariant ones::

            sage: t['i_(jk)'] == t.symmetrize(1,2)
            True

        The LaTeX notation can be employed::

            sage: t['^{i}_{(jk)}'] == t.symmetrize(1,2)
            True
        """
        if not pos:
            pos = range(self._tensor_rank)
        # check whether the symmetrization is possible:
        pos_cov = self._tensor_type[0]   # first covariant position
        pos0 = pos[0]
        if pos0 < pos_cov:  # pos0 is a contravariant position
            for k in range(1,len(pos)):
                if pos[k] >= pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a contravariant position, while " +
                        str(pos[k]) + " is a covariant position; \n"
                        "symmetrization is meaningful only on tensor " +
                        "arguments of the same type")
        else:  # pos0 is a covariant position
            for k in range(1,len(pos)):
                if pos[k] < pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a covariant position, while " +
                        str(pos[k]) + " is a contravariant position; \n"
                        "symmetrization is meaningful only on tensor " +
                        "arguments of the same type")
        if 'basis' in kwargs:
            basis = kwargs['basis']
        else:
            basis = self.pick_a_basis()
        res_comp = self._components[basis].symmetrize(*pos)
        return self._fmodule.tensor_from_comp(self._tensor_type, res_comp)

    def antisymmetrize(self, *pos, **kwargs):
        r"""
        Antisymmetrization over some arguments.

        INPUT:

        - ``pos`` -- list of argument positions involved in the
          antisymmetrization (with the convention ``position=0`` for the first
          argument); if none, the antisymmetrization is performed over all the
          arguments
        - ``basis`` -- (default: ``None``) module basis with respect to which
          the component computation is to be performed; if none, the module's
          default basis is used if the tensor field has already components
          in it; otherwise another basis w.r.t. which the tensor has
          components will be picked

        OUTPUT:

        - the antisymmetrized tensor (instance of :class:`FreeModuleTensor`)

        EXAMPLES:

        Antisymmetrization of a tensor of type `(2,0)`::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((2,0))
            sage: t[:] = [[1,-2,3], [4,5,6], [7,8,-9]]
            sage: s = t.antisymmetrize() ; s
            Alternating contravariant tensor of degree 2 on the 3-dimensional
             vector space M over the Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 1)
            sage: t[:], s[:]
            (
            [ 1 -2  3]  [ 0 -3 -2]
            [ 4  5  6]  [ 3  0 -1]
            [ 7  8 -9], [ 2  1  0]
            )
            sage: all(s[i,j] == 1/2*(t[i,j]-t[j,i])   # Check:
            ....:     for i in range(3) for j in range(3))
            True
            sage: s.antisymmetrize() == s  # another test
            True
            sage: t.antisymmetrize() == t.antisymmetrize(0,1)
            True

        Antisymmetrization of a tensor of type `(0, 3)` over the first two
        arguments::

            sage: t = M.tensor((0,3))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]],
            ....:         [[10,-11,12], [13,14,-15], [16,17,18]],
            ....:         [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.antisymmetrize(0,1) ; s  # (0,1) = the first two arguments
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 1)
            sage: s[:]
            [[[0, 0, 0], [-7, 8, -3], [-6, 14, 6]],
             [[7, -8, 3], [0, 0, 0], [19, -3, -3]],
             [[6, -14, -6], [-19, 3, 3], [0, 0, 0]]]
            sage: all(s[i,j,k] == 1/2*(t[i,j,k]-t[j,i,k])   # Check:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.antisymmetrize(0,1) == s  # another test
            True
            sage: s.symmetrize(0,1) == 0  # of course
            True

        Instead of invoking the method :meth:`antisymmetrize`, one can use
        the index notation with square brackets denoting the
        antisymmetrization; it suffices to pass the indices as a string
        inside square brackets::

            sage: s1 = t['_[ij]k'] ; s1
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s1.symmetries()
            no symmetry;  antisymmetry: (0, 1)
            sage: s1 == s
            True

        The LaTeX notation is recognized::

            sage: t['_{[ij]k}'] == s
            True

        Note that in the index notation, the name of the indices is irrelevant;
        they can even be replaced by dots::

            sage: t['_[..].'] == s
            True

        Antisymmetrization of a tensor of type `(0,3)` over the first and last
        arguments::

            sage: s = t.antisymmetrize(0,2) ; s  # (0,2) = first and last arguments
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 2)
            sage: s[:]
            [[[0, -4, -8], [0, -4, 14], [0, -4, -17]],
             [[4, 0, 16], [4, 0, -19], [4, 0, -4]],
             [[8, -16, 0], [-14, 19, 0], [17, 4, 0]]]
            sage: all(s[i,j,k] == 1/2*(t[i,j,k]-t[k,j,i])   # Check:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.antisymmetrize(0,2) == s  # another test
            True
            sage: s.symmetrize(0,2) == 0  # of course
            True
            sage: s.symmetrize(0,1) == 0  # no reason for this to hold
            False

        Antisymmetrization of a tensor of type `(0,3)` over the last two
        arguments::

            sage: s = t.antisymmetrize(1,2) ; s  # (1,2) = the last two arguments
            Type-(0,3) tensor on the 3-dimensional vector space M over the
             Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (1, 2)
            sage: s[:]
            [[[0, 3, -2], [-3, 0, -1], [2, 1, 0]],
             [[0, -12, -2], [12, 0, -16], [2, 16, 0]],
             [[0, 1, -23], [-1, 0, -1], [23, 1, 0]]]
            sage: all(s[i,j,k] == 1/2*(t[i,j,k]-t[i,k,j])   # Check:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.antisymmetrize(1,2) == s  # another test
            True
            sage: s.symmetrize(1,2) == 0  # of course
            True

        The index notation can be used instead of the explicit call to
        :meth:`antisymmetrize`::

            sage: t['_i[jk]'] == t.antisymmetrize(1,2)
            True

        Full antisymmetrization of a tensor of type `(0,3)`::

            sage: s = t.antisymmetrize() ; s
            Alternating form of degree 3 on the 3-dimensional vector space M
             over the Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 1, 2)
            sage: s[:]
            [[[0, 0, 0], [0, 0, 2/3], [0, -2/3, 0]],
             [[0, 0, -2/3], [0, 0, 0], [2/3, 0, 0]],
             [[0, 2/3, 0], [-2/3, 0, 0], [0, 0, 0]]]
            sage: all(s[i,j,k] == 1/6*(t[i,j,k]-t[i,k,j]+t[j,k,i]-t[j,i,k]
            ....:                      +t[k,i,j]-t[k,j,i])
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s.antisymmetrize() == s  # another test
            True
            sage: s.symmetrize(0,1) == 0  # of course
            True
            sage: s.symmetrize(0,2) == 0  # of course
            True
            sage: s.symmetrize(1,2) == 0  # of course
            True
            sage: t.antisymmetrize() == t.antisymmetrize(0,1,2)
            True

        The index notation can be used instead of the explicit call to
        :meth:`antisymmetrize`::

            sage: t['_[ijk]'] == t.antisymmetrize()
            True
            sage: t['_[abc]'] == t.antisymmetrize()
            True
            sage: t['_[...]'] == t.antisymmetrize()
            True
            sage: t['_{[ijk]}'] == t.antisymmetrize() # LaTeX notation
            True

        Antisymmetrization can be performed only on arguments on the same type::

            sage: t = M.tensor((1,2))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]],
            ....:         [[10,-11,12], [13,14,-15], [16,17,18]],
            ....:         [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.antisymmetrize(0,1)
            Traceback (most recent call last):
            ...
            TypeError: 0 is a contravariant position, while 1 is a covariant position;
            antisymmetrization is meaningful only on tensor arguments of the same type
            sage: s = t.antisymmetrize(1,2) # OK: both 1 and 2 are covariant positions

        The order of positions does not matter::

            sage: t.antisymmetrize(2,1) == t.antisymmetrize(1,2)
            True

        Again, the index notation can be used::

            sage: t['^i_[jk]'] == t.antisymmetrize(1,2)
            True
            sage: t['^i_{[jk]}'] == t.antisymmetrize(1,2)  # LaTeX notation
            True

        The character '^' can be skipped::

            sage: t['i_[jk]'] == t.antisymmetrize(1,2)
            True
        """
        if not pos:
            pos = range(self._tensor_rank)
        # check whether the antisymmetrization is possible:
        pos_cov = self._tensor_type[0]   # first covariant position
        pos0 = pos[0]
        if pos0 < pos_cov:  # pos0 is a contravariant position
            for k in range(1,len(pos)):
                if pos[k] >= pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a contravariant position, while " +
                        str(pos[k]) + " is a covariant position; \n"
                        "antisymmetrization is meaningful only on tensor " +
                        "arguments of the same type")
        else:  # pos0 is a covariant position
            for k in range(1,len(pos)):
                if pos[k] < pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a covariant position, while " +
                        str(pos[k]) + " is a contravariant position; \n"
                        "antisymmetrization is meaningful only on tensor " +
                        "arguments of the same type")
        if 'basis' in kwargs:
            basis = kwargs['basis']
        else:
            basis = self.pick_a_basis()
        res_comp = self._components[basis].antisymmetrize(*pos)
        return self._fmodule.tensor_from_comp(self._tensor_type, res_comp)
