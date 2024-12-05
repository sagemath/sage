r"""
Free module bases

The class :class:`FreeModuleBasis` implements bases on a free module `M` of
finite rank over a commutative ring,
while the class :class:`FreeModuleCoBasis` implements the dual bases (i.e.
bases of the dual module `M^*`).

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Travis Scrimshaw (2016): ABC Basis_abstract and list functionality for bases
  (:issue:`20770`)
- Eric Gourgoulhon (2018): some refactoring and more functionalities in the
  choice of symbols for basis elements (:issue:`24792`)

REFERENCES:

- Chap. 10 of R. Godement : *Algebra* [God1968]_
- Chap. 3 of S. Lang : *Algebra* [Lan2002]_
"""
#******************************************************************************
#       Copyright (C) 2015, 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.rings.integer_ring import ZZ
from sage.sets.family import AbstractFamily
from sage.structure.unique_representation import UniqueRepresentation


class Basis_abstract(UniqueRepresentation, AbstractFamily):
    """
    Abstract base class for (dual) bases of free modules.

    A basis is an :class:`~sage.sets.family.AbstractFamily`, hence like
    :class:`collections.abc.Mapping` subclasses such as :class:`dict`, it is
    an associative :class:`Container`, providing methods :meth:`keys`,
    :meth:`values`, and :meth:`items`. Thus, ``e[i]`` returns the element
    of the basis ``e`` indexed by the key ``i``. However, in contrast to
    :class:`Mapping` subclasses, not the :meth:`keys` but the
    :meth:`values` are considered the elements.

    EXAMPLES::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: e = M.basis('e'); e
        Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
        sage: list(e)
        [Element e_1 of the Rank-3 free module M over the Integer Ring,
        Element e_2 of the Rank-3 free module M over the Integer Ring,
        Element e_3 of the Rank-3 free module M over the Integer Ring]
        sage: e.category()
        Category of facade finite enumerated sets
        sage: list(e.keys())
        [1, 2, 3]
        sage: list(e.values())
        [Element e_1 of the Rank-3 free module M over the Integer Ring,
        Element e_2 of the Rank-3 free module M over the Integer Ring,
        Element e_3 of the Rank-3 free module M over the Integer Ring]
        sage: list(e.items())
        [(1, Element e_1 of the Rank-3 free module M over the Integer Ring),
        (2, Element e_2 of the Rank-3 free module M over the Integer Ring),
        (3, Element e_3 of the Rank-3 free module M over the Integer Ring)]
    """
    def __init__(self, fmodule, symbol, latex_symbol, indices, latex_indices):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._fmodule is M
            True
        """
        self._fmodule = fmodule
        self._symbol = symbol
        self._latex_symbol = latex_symbol
        self._indices = indices
        self._latex_indices = latex_indices
        super().__init__(category=FiniteEnumeratedSets(), facade=fmodule)

    def keys(self):
        """
        Return the keys (indices) of the family.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: list(e.keys())
            [0, 1, 2]
        """
        return self._fmodule.irange()

    def values(self):
        """
        Return the basis elements of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: list(e.values())
            [Element e_0 of the Rank-3 free module M over the Integer Ring,
             Element e_1 of the Rank-3 free module M over the Integer Ring,
             Element e_2 of the Rank-3 free module M over the Integer Ring]
        """
        return self._vec

    def _element_constructor_(self, x):
        """
        Test whether ``x`` is an element of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e(e[1])
            Element e_1 of the Rank-3 free module M over the Integer Ring
            sage: f = M.basis('f')
            sage: e(f[1])
            Traceback (most recent call last):
            ...
            ValueError: no common basis for the comparison
        """
        if x in self.values():
            return x
        raise ValueError(f'{x} is not in {self}')

    def __iter__(self):
        r"""
        Return an iterator for the basis elements of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: list(e)
            [Element e_0 of the Rank-3 free module M over the Integer Ring,
             Element e_1 of the Rank-3 free module M over the Integer Ring,
             Element e_2 of the Rank-3 free module M over the Integer Ring]
            sage: ed = e.dual_basis()
            sage: list(ed)
            [Linear form e^0 on the Rank-3 free module M over the Integer Ring,
             Linear form e^1 on the Rank-3 free module M over the Integer Ring,
             Linear form e^2 on the Rank-3 free module M over the Integer Ring]

        An example with indices starting at 1 instead of 0::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M1',
            ....:                          start_index=1)
            sage: e = M.basis('e')
            sage: list(e)
            [Element e_1 of the Rank-3 free module M1 over the Integer Ring,
             Element e_2 of the Rank-3 free module M1 over the Integer Ring,
             Element e_3 of the Rank-3 free module M1 over the Integer Ring]
        """
        yield from self.values()

    def _test_iter_len(self, **options):
        r"""
        Test that __iter__ and __len__ work correctly.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._test_iter_len()
        """
        tester = self._tester(**options)
        g = iter(self)
        b = list(g)
        for x in b:
            tester.assertIn(x, self.free_module())
        tester.assertEqual(len(b), len(self))
        tester.assertEqual(len(b), self.free_module().rank())

    def __len__(self):
        r"""
        Return the basis length, i.e. the rank of the free module.

        .. NOTE::

            the method ``__len__()`` is required for the basis to act as a
            "frame" in the class :class:`~sage.tensor.modules.comp.Components`.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.__len__()
            3
            sage: len(e)
            3
        """
        return self._fmodule._rank

    def cardinality(self):
        r"""
        Return the basis length, i.e. the rank of the free module.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.cardinality()
            3
        """
        return ZZ(self._fmodule._rank)

    def __getitem__(self, index):
        r"""
        Return the basis element corresponding to a given index.

        INPUT:

        - ``index`` -- the index of the basis element

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.__getitem__(0)
            Element e_0 of the Rank-3 free module M over the Integer Ring
            sage: e.__getitem__(1)
            Element e_1 of the Rank-3 free module M over the Integer Ring
            sage: e.__getitem__(2)
            Element e_2 of the Rank-3 free module M over the Integer Ring
            sage: e[1] is e.__getitem__(1)
            True
            sage: e[1].parent() is M
            True
            sage: e[:]
            (Element e_0 of the Rank-3 free module M over the Integer Ring,
             Element e_1 of the Rank-3 free module M over the Integer Ring,
             Element e_2 of the Rank-3 free module M over the Integer Ring)
            sage: f = e.dual_basis()
            sage: f[0]
            Linear form e^0 on the Rank-3 free module M over the Integer Ring

        Examples with ``start_index=1``::

            sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e1 = M1.basis('e')
            sage: e1.__getitem__(1)
            Element e_1 of the Rank-3 free module M over the Integer Ring
            sage: e1.__getitem__(2)
            Element e_2 of the Rank-3 free module M over the Integer Ring
            sage: e1.__getitem__(3)
            Element e_3 of the Rank-3 free module M over the Integer Ring
        """
        si = self._fmodule._sindex
        if isinstance(index, slice):
            start, stop = index.start, index.stop
            if start is not None:
                start -= si
            if stop is not None:
                stop -= si
            return self._vec[start:stop:index.step]
        n = self._fmodule._rank
        i = index - si
        if i < 0 or i > n-1:
            raise IndexError("out of range: {} not in [{},{}]".format(i+si, si,
                                                                      n-1+si))
        return self._vec[i]

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._latex_()
            '\\left(e_{0},e_{1},e_{2}\\right)'
            sage: latex(e)
            \left(e_{0},e_{1},e_{2}\right)
            sage: f = M.basis('eps', latex_symbol=r'\epsilon')
            sage: f._latex_()
            '\\left(\\epsilon_{0},\\epsilon_{1},\\epsilon_{2}\\right)'
            sage: latex(f)
            \left(\epsilon_{0},\epsilon_{1},\epsilon_{2}\right)

        ::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = e.dual_basis()
            sage: f._latex_()
            '\\left(e^{0},e^{1},e^{2}\\right)'
        """
        return self._latex_name

    def free_module(self):
        """
        Return the free module of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: e.free_module() is M
            True
        """
        return self._fmodule

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='down'):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the elements of ``self``, or a list of strings,
          representing the individual symbols of the elements of ``self``
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the elements of ``self``,
          or a list of strings, representing the individual LaTeX symbols of
          the elements of ``self``; if ``None``, ``symbol`` is used in place
          of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the elements of ``self``; if ``None``, the indices will be generated
          as integers within the range declared on the free module on which
          ``self`` is defined
        - ``latex_indices`` -- (default: ``None``) list of strings representing
          the indices for the LaTeX symbols of the elements of ``self``; if
          ``None``, ``indices`` is used instead
        - ``index_position`` -- (default: ``'down'``) determines the position
          of the indices labelling the individual elements of ``self``; can be
          either ``'down'`` or ``'up'``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e'); e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: e.set_name('f'); e
            Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring
            sage: e.set_name(['a', 'b', 'c']); e
            Basis (a,b,c) on the Rank-3 free module M over the Integer Ring
            sage: e.set_name('e', indices=['x', 'y', 'z']); e
            Basis (e_x,e_y,e_z) on the Rank-3 free module M over the Integer Ring
            sage: e.set_name('e', index_position='up'); e
            Basis (e^0,e^1,e^2) on the Rank-3 free module M over the Integer Ring
            sage: latex(e)
            \left(e^{0},e^{1},e^{2}\right)
            sage: e.set_name('e', latex_symbol=r'\epsilon'); e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: latex(e)
            \left(\epsilon_{0},\epsilon_{1},\epsilon_{2}\right)
            sage: e.set_name('e', latex_symbol=[r'\alpha', r'\beta', r'\gamma'])
            sage: latex(e)
            \left(\alpha,\beta,\gamma\right)
            sage: e.set_name('e', latex_symbol='E',
            ....:            latex_indices=[r'\alpha', r'\beta', r'\gamma'])
            sage: latex(e)
            \left(E_{\alpha},E_{\beta},E_{\gamma}\right)
            sage: e.set_name('e') # back to the default
        """
        n = self._fmodule._rank
        if index_position == "down":
            pos = "_"
        else:
            pos = "^"
        if latex_symbol is None:
            latex_symbol = symbol
        self._symbol = symbol
        self._latex_symbol = latex_symbol
        self._indices = indices
        self._latex_indices = latex_indices
        # Text symbols:
        if isinstance(symbol, (list, tuple)):
            if len(symbol) != n:
                raise ValueError("symbol must contain {} strings".format(n))
            if len(set(symbol)) != n:
                raise ValueError("the individual symbols must be different")
        else:
            if indices is None:
                indices = [str(i) for i in self._fmodule.irange()]
            elif len(indices) != n:
                raise ValueError("indices must contain {} elements".format(n))
            symbol0 = symbol + pos
            symbol = [symbol0 + i for i in indices]
        # LaTeX symbols:
        if isinstance(latex_symbol, (list, tuple)):
            if len(latex_symbol) != n:
                raise ValueError(
                              "latex_symbol must contain {} strings".format(n))
            if len(set(latex_symbol)) != n:
                raise ValueError("the individual symbols must be different")
        else:
            if latex_indices is None:
                if indices is None:
                    latex_indices = [str(i) for i in self._fmodule.irange()]
                else:
                    latex_indices = indices
            elif len(latex_indices) != n:
                raise ValueError(
                            "latex_indices must contain {} elements".format(n))
            symbol0 = latex_symbol + pos
            latex_symbol = [symbol0 + "{" + i + "}" for i in latex_indices]
        # Setting the names
        self._name = "(" + ",".join(symbol) + ")"
        self._latex_name = r"\left(" + ",".join(latex_symbol) + r"\right)"
        for i in range(n):
            self._vec[i].set_name(symbol[i], latex_name=latex_symbol[i])

#******************************************************************************


class FreeModuleCoBasis(Basis_abstract):
    r"""
    Dual basis of a free module over a commutative ring.

    INPUT:

    - ``basis`` -- basis of a free module `M` of which ``self`` is the dual
      (must be an instance of :class:`FreeModuleBasis`)
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the elements of the cobasis, or a tuple of strings,
      representing the individual symbols of the elements of the cobasis
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the elements of the cobasis,
      or a tuple of strings, representing the individual LaTeX symbols of
      the elements of the cobasis; if ``None``, ``symbol`` is used in place
      of ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling
      the elements of the cobasis; if ``None``, the indices will be generated
      as integers within the range declared on the free module on which the
      cobasis is defined
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the elements of the cobasis; if
      ``None``, ``indices`` is used instead

    EXAMPLES:

    Dual basis on a rank-3 free module::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: e = M.basis('e') ; e
        Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
        sage: from sage.tensor.modules.free_module_basis import FreeModuleCoBasis
        sage: f = FreeModuleCoBasis(e, 'f') ; f
        Dual basis (f^1,f^2,f^3) on the Rank-3 free module M over the Integer Ring

    Instead of importing ``FreeModuleCoBasis`` in the global name space, it is
    recommended to use the method
    :meth:`~sage.tensor.modules.free_module_basis.FreeModuleBasis.dual_basis`
    of the basis ``e``::

        sage: f = e.dual_basis() ; f
        Dual basis (e^1,e^2,e^3) on the Rank-3 free module M over the Integer Ring

    Let us check that the elements of ``f`` are in the dual of ``M``::

        sage: f[1]
        Linear form e^1 on the Rank-3 free module M over the Integer Ring
        sage: f[1] in M.dual()
        True

    and that ``f`` is indeed the dual of ``e``::

        sage: f[1](e[1]), f[1](e[2]), f[1](e[3])
        (1, 0, 0)
        sage: f[2](e[1]), f[2](e[2]), f[2](e[3])
        (0, 1, 0)
        sage: f[3](e[1]), f[3](e[2]), f[3](e[3])
        (0, 0, 1)

    TESTS::

        sage: TestSuite(f).run()
    """
    def __init__(self, basis, symbol, latex_symbol=None, indices=None,
                 latex_indices=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_basis import FreeModuleCoBasis
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = FreeModuleCoBasis(e, 'f')
            sage: TestSuite(f).run()
        """
        self._basis = basis
        Basis_abstract.__init__(self, basis._fmodule, symbol, latex_symbol,
                                indices, latex_indices)
        # The individual linear forms:
        vl = []
        fmodule = self._fmodule
        ring_one = fmodule._ring.one()
        for i in fmodule.irange():
            v = fmodule.linear_form()
            v.set_comp(basis)[i] = ring_one
            vl.append(v)
        self._vec = tuple(vl)
        # The names:
        self.set_name(symbol, latex_symbol=latex_symbol, indices=indices,
                      latex_indices=latex_indices, index_position='up')

    def _test_iter_len(self, **options):
        r"""
        Test that __iter__ and __len__ work correctly.

        This method overrides ``Basis_abstract`` so that containment
        of elements in the dual of ``self.free_module()`` is tested instead.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = e.dual_basis()
            sage: f._test_iter_len()
        """
        tester = self._tester(**options)
        g = iter(self)
        b = list(g)
        for x in b:
            tester.assertIn(x, self.free_module().dual())
        tester.assertEqual(len(b), len(self))
        tester.assertEqual(len(b), self.free_module().rank())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = e.dual_basis()
            sage: f
            Dual basis (e^0,e^1,e^2) on the
             Rank-3 free module M over the Integer Ring
        """
        return "Dual basis {} on the {}".format(self._name, self._fmodule)

#******************************************************************************


class FreeModuleBasis(Basis_abstract):
    r"""
    Basis of a free module over a commutative ring `R`.

    INPUT:

    - ``fmodule`` -- free module `M` (as an instance of
      :class:`FiniteRankFreeModule`)
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the elements of the basis, or a tuple of strings,
      representing the individual symbols of the elements of the basis
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the elements of the basis,
      or a tuple of strings, representing the individual LaTeX symbols of
      the elements of the basis; if ``None``, ``symbol`` is used in place
      of ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling
      the elements of the basis; if ``None``, the indices will be generated
      as integers within the range declared on ``fmodule``
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the elements of the basis; if
      ``None``, ``indices`` is used instead
    - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
      dual basis; if ``None``, ``symbol`` must be a string and is used
      for the common base of the symbols of the elements of the dual basis
    - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
      but for the dual basis

    EXAMPLES:

    A basis on a rank-3 free module over `\ZZ`::

        sage: M0 = FiniteRankFreeModule(ZZ, 3, name='M_0')
        sage: from sage.tensor.modules.free_module_basis import FreeModuleBasis
        sage: e = FreeModuleBasis(M0, 'e') ; e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M_0 over the Integer Ring

    Instead of importing ``FreeModuleBasis`` in the global name space, it is
    recommended to use the module's method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.basis`::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e') ; e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring

    The individual elements constituting the basis are accessed via the
    square bracket operator::

        sage: e[0]
        Element e_0 of the Rank-3 free module M over the Integer Ring
        sage: e[0] in M
        True

    The slice operator ``:`` can be used to access to more than one element::

        sage: e[0:2]
        (Element e_0 of the Rank-3 free module M over the Integer Ring,
         Element e_1 of the Rank-3 free module M over the Integer Ring)
        sage: e[:]
        (Element e_0 of the Rank-3 free module M over the Integer Ring,
         Element e_1 of the Rank-3 free module M over the Integer Ring,
         Element e_2 of the Rank-3 free module M over the Integer Ring)

    The LaTeX symbol can be set explicitly::

        sage: latex(e)
        \left(e_{0},e_{1},e_{2}\right)
        sage: eps = M.basis('eps', latex_symbol=r'\epsilon') ; eps
        Basis (eps_0,eps_1,eps_2) on the Rank-3 free module M over the Integer
         Ring
        sage: latex(eps)
        \left(\epsilon_{0},\epsilon_{1},\epsilon_{2}\right)

    The individual elements of the basis are labelled according the
    parameter ``start_index`` provided at the free module construction::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: e = M.basis('e') ; e
        Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
        sage: e[1]
        Element e_1 of the Rank-3 free module M over the Integer Ring

    It is also possible to fully customize the labels, via the argument
    ``indices``::

        sage: f = M.basis('f', indices=('x', 'y', 'z')); f
        Basis (f_x,f_y,f_z) on the Rank-3 free module M over the Integer Ring
        sage: f[1]
        Element f_x of the Rank-3 free module M over the Integer Ring

    The symbol of each element of the basis can also be freely chosen, by
    providing a tuple of symbols as the first argument of ``basis``; it is then
    mandatory to specify some symbols for the dual basis as well::

        sage: g = M.basis(('a', 'b', 'c'), symbol_dual=('A', 'B', 'C')); g
        Basis (a,b,c) on the Rank-3 free module M over the Integer Ring
        sage: g[1]
        Element a of the Rank-3 free module M over the Integer Ring
        sage: g.dual_basis()[1]
        Linear form A on the Rank-3 free module M over the Integer Ring

    TESTS::

        sage: TestSuite(e).run()
        sage: TestSuite(f).run()
        sage: TestSuite(g).run()
    """
    # The following class attribute must be redefined by any derived class:
    _cobasis_class = FreeModuleCoBasis

    @staticmethod
    def __classcall_private__(cls, fmodule, symbol, latex_symbol=None,
                              indices=None, latex_indices=None,
                              symbol_dual=None, latex_symbol_dual=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: from sage.tensor.modules.free_module_basis import FreeModuleBasis
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = FreeModuleBasis(M, 'e', latex_symbol='e')
            sage: e is FreeModuleBasis(M, 'e')
            True
        """
        if latex_symbol is None:
            latex_symbol = symbol
        # Only tuples are valid entries for the unique representation of
        # FreeModuleBasis:
        if isinstance(symbol, list):
            symbol = tuple(symbol)
        if isinstance(latex_symbol, list):
            latex_symbol = tuple(latex_symbol)
        if isinstance(indices, list):
            indices = tuple(indices)
        if isinstance(latex_indices, list):
            latex_indices = tuple(latex_indices)
        if isinstance(symbol_dual, list):
            symbol_dual = tuple(symbol_dual)
        if isinstance(latex_symbol_dual, list):
            latex_symbol_dual = tuple(latex_symbol_dual)
        return super().__classcall__(cls, fmodule, symbol,
                                           latex_symbol=latex_symbol,
                                           indices=indices,
                                           latex_indices=latex_indices,
                                           symbol_dual=symbol_dual,
                                           latex_symbol_dual=latex_symbol_dual)

    def __init__(self, fmodule, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, symbol_dual=None, latex_symbol_dual=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: from sage.tensor.modules.free_module_basis import FreeModuleBasis
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = FreeModuleBasis(M, 'e', latex_symbol=r'\epsilon')
            sage: TestSuite(e).run()
        """
        Basis_abstract.__init__(self, fmodule, symbol, latex_symbol, indices,
                                latex_indices)
        # The basis is added to the module list of bases
        fmodule._known_bases.append(self)
        # The individual vectors:
        vl = []
        ring_one = fmodule._ring.one()
        for i in fmodule.irange():
            v = fmodule.element_class(fmodule)
            v.set_comp(self)[i] = ring_one
            vl.append(v)
        self._vec = tuple(vl)
        # The names:
        self.set_name(symbol, latex_symbol=latex_symbol, indices=indices,
                      latex_indices=latex_indices, index_position='down')
        # The first defined basis is considered as the default one:
        if fmodule._def_basis is None:
            fmodule._def_basis = self
        # Initialization of the components w.r.t the current basis of the zero
        # elements of all tensor modules constructed up to now (including the
        # base module itself, since it is considered as a type-(1,0) tensor
        # module):
        for t in fmodule._all_modules:
            t.zero()._add_comp_unsafe(self)
            # (since new components are initialized to zero)
        # Initialization of the components w.r.t. the current basis of the
        # identity map of the general linear group:
        if fmodule._general_linear_group is not None:
            from .comp import KroneckerDelta
            gl = fmodule._general_linear_group
            gl.one()._components[self] = KroneckerDelta(fmodule._ring, self,
                                    start_index=fmodule._sindex,
                                    output_formatter=fmodule._output_formatter)
        # The dual basis:
        self._symbol_dual = symbol_dual
        self._latex_symbol_dual = latex_symbol_dual
        if symbol_dual is None:
            if isinstance(symbol, (list, tuple)):
                raise ValueError("symbol_dual must be provided")
            else:
                symbol_dual = symbol
        elif latex_symbol_dual is None:
            latex_symbol_dual = symbol_dual
        if latex_symbol_dual is None:
            latex_symbol_dual = latex_symbol
        self._dual_basis = type(self)._cobasis_class(self, symbol_dual,
                                                latex_symbol=latex_symbol_dual,
                                                indices=indices,
                                                latex_indices=latex_indices)

    ###### Methods to be redefined by derived classes of FreeModuleBasis ######

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e1 = M1.basis('e')
            sage: e1
            Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
        """
        return "Basis {} on the {}".format(self._name, self._fmodule)

    def _new_instance(self, symbol, latex_symbol=None, indices=None,
                      latex_indices=None, symbol_dual=None,
                      latex_symbol_dual=None):
        r"""
        Construct a new basis on the same module as ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the elements of the basis, or a tuple of strings,
          representing the individual symbols of the elements of the basis
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the elements of the basis,
          or a tuple of strings, representing the individual LaTeX symbols of
          the elements of the basis; if ``None``, ``symbol`` is used in place
          of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the elements of the basis; if ``None``, the indices will be generated
          as integers within the range declared on the free module on which the
          ``self`` is defined
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the elements of the
          basis; if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual basis; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual basis
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual basis

        OUTPUT: instance of :class:`FreeModuleBasis`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._new_instance('f')
            Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring
            sage: e._new_instance(('a', 'b', 'c'), symbol_dual=('A', 'B', 'C'))
            Basis (a,b,c) on the Rank-3 free module M over the Integer Ring
            sage: _.dual_basis()
            Dual basis (A,B,C) on the Rank-3 free module M over the Integer Ring
            sage: e._new_instance('E', indices=('x', 'y', 'z'))
            Basis (E_x,E_y,E_z) on the Rank-3 free module M over the Integer Ring
            sage: _.dual_basis()
            Dual basis (E^x,E^y,E^z) on the Rank-3 free module M over the Integer Ring
        """
        return FreeModuleBasis(self._fmodule, symbol, latex_symbol=latex_symbol,
                               indices=indices, latex_indices=latex_indices,
                               symbol_dual=symbol_dual,
                               latex_symbol_dual=latex_symbol_dual)

    ###### End of methods to be redefined by derived classes ######

    def _init_from_family(self, family):
        r"""
        Identify ``self`` to a linearly independent spanning family of
        module elements.

        INPUT:

        - ``family`` -- a family of elements of ``self.free_module()`` that are
          linearly independent and spanning ``self.free_module()``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: f1 = e[1] + e[2]
            sage: f2 = e[1] - e[2]
            sage: f = M.basis('f')
            sage: f._init_from_family([f1, f2])
            sage: (f[1], f[2]) == (f1, f2)
            True
            sage: f[1].display()
            f_1 = e_1 + e_2
            sage: f[2].display()
            f_2 = e_1 - e_2
            sage: e[1].display(f)
            e_1 = 1/2 f_1 + 1/2 f_2
            sage: e[2].display(f)
            e_2 = 1/2 f_1 - 1/2 f_2
        """
        fmodule = self._fmodule
        n = fmodule.rank()
        if len(family) != n:
            raise ValueError("the size of the family must be {}".format(n))
        # Copy of the components of each element of the family:
        for i, ff in enumerate(family):
            if ff not in fmodule:
                raise TypeError("{} is not an element of {}".format(ff,
                                                                    fmodule))
            vs = self._vec[i]
            for basis, comp in ff._components.items():
                vs._components[basis] = comp.copy()
        # Automorphisms relating the family to previously defined bases are
        # constructed from the components of the family elements and are
        # registered as changes of basis to ``self``:
        ff0 = family[0]
        for basis in ff0._components:
            try:
                comps = [ff.components(basis) for ff in family]
            except ValueError:
                continue
            aut = fmodule.automorphism()
            mat = [[c[[i]] for c in comps] for i in fmodule.irange()]
            aut.add_comp(basis)[:] = mat
            aut.add_comp(self)[:] = mat
            fmodule.set_change_of_basis(basis, self, aut)

    def module(self):
        r"""
        Return the free module on which the basis is defined.

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
          representing the free module of which ``self`` is a basis

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.module()
            Rank-3 free module M over the Integer Ring
            sage: e.module() is M
            True
        """
        return self._fmodule

    def dual_basis(self):
        r"""
        Return the basis dual to ``self``.

        OUTPUT:

        - instance of :class:`FreeModuleCoBasis` representing the dual of
          ``self``

        EXAMPLES:

        Dual basis on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e') ; e
            Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
            sage: f = e.dual_basis() ; f
            Dual basis (e^1,e^2,e^3) on the Rank-3 free module M over the Integer Ring

        Let us check that the elements of f are elements of the dual of M::

            sage: f[1] in M.dual()
            True
            sage: f[1]
            Linear form e^1 on the Rank-3 free module M over the Integer Ring

        and that f is indeed the dual of e::

            sage: f[1](e[1]), f[1](e[2]), f[1](e[3])
            (1, 0, 0)
            sage: f[2](e[1]), f[2](e[2]), f[2](e[3])
            (0, 1, 0)
            sage: f[3](e[1]), f[3](e[2]), f[3](e[3])
            (0, 0, 1)
        """
        return self._dual_basis

    def new_basis(self, change_of_basis, symbol, latex_symbol=None,
                  indices=None, latex_indices=None, symbol_dual=None,
                  latex_symbol_dual=None):
        r"""
        Define a new module basis from ``self``.

        The new basis is defined by means of a module automorphism.

        INPUT:

        - ``change_of_basis`` -- instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          describing the automorphism `P` that relates the current basis
          `(e_i)` (described by ``self``) to the new basis `(n_i)` according
          to `n_i = P(e_i)`
        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the elements of the basis, or a tuple of strings,
          representing the individual symbols of the elements of the basis
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the elements of the basis,
          or a tuple of strings, representing the individual LaTeX symbols of
          the elements of the basis; if ``None``, ``symbol`` is used in place
          of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the elements of the basis; if ``None``, the indices will be generated
          as integers within the range declared on the free module on which
          ``self`` is defined
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the elements of the
          basis; if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual basis; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual basis
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual basis

        OUTPUT:

        - the new basis `(n_i)`, as an instance of :class:`FreeModuleBasis`

        EXAMPLES:

        Change of basis on a vector space of dimension 2::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism()
            sage: a[:] = [[1, 2], [-1, 3]]
            sage: f = e.new_basis(a, 'f') ; f
            Basis (f_1,f_2) on the 2-dimensional vector space M over the
             Rational Field
            sage: f[1].display()
            f_1 = e_1 - e_2
            sage: f[2].display()
            f_2 = 2 e_1 + 3 e_2
            sage: e[1].display(f)
            e_1 = 3/5 f_1 + 1/5 f_2
            sage: e[2].display(f)
            e_2 = -2/5 f_1 + 1/5 f_2

        Use of some keyword arguments::

            sage: b = e.new_basis(a, 'b', indices=('x', 'y'),
            ....:                 symbol_dual=('A', 'B'))
            sage: b
            Basis (b_x,b_y) on the 2-dimensional vector space M over the
             Rational Field
            sage: b.dual_basis()
            Dual basis (A,B) on the 2-dimensional vector space M over the
             Rational Field
        """
        from .free_module_automorphism import FreeModuleAutomorphism
        if not isinstance(change_of_basis, FreeModuleAutomorphism):
            raise TypeError("the argument change_of_basis must be some " +
                            "instance of FreeModuleAutomorphism")
        fmodule = self._fmodule
        # self._new_instance used instead of FreeModuleBasis for a correct
        # construction in case of derived classes:
        the_new_basis = self._new_instance(symbol, latex_symbol=latex_symbol,
                                           indices=indices,
                                           latex_indices=latex_indices,
                                           symbol_dual=symbol_dual,
                                           latex_symbol_dual=latex_symbol_dual)
        transf = change_of_basis.copy()
        inv_transf = change_of_basis.inverse().copy()
        si = fmodule._sindex
        # Components of the new basis vectors in the old basis:
        for i in fmodule.irange():
            for j in fmodule.irange():
                the_new_basis._vec[i-si].add_comp(self)[[j]] = \
                                                  transf.comp(self)[[j,i]]
        # Components of the new dual-basis elements in the old dual basis:
        for i in fmodule.irange():
            for j in fmodule.irange():
                the_new_basis._dual_basis._vec[i-si].add_comp(self)[[j]] = \
                                              inv_transf.comp(self)[[i,j]]
        # The components of the transformation and its inverse are the same in
        # the two bases:
        for i in fmodule.irange():
            for j in fmodule.irange():
                transf.add_comp(the_new_basis)[[i,j]] = transf.comp(self)[[i,j]]
                inv_transf.add_comp(the_new_basis)[[i,j]] = \
                                              inv_transf.comp(self)[[i,j]]
        # Components of the old basis vectors in the new basis:
        for i in fmodule.irange():
            for j in fmodule.irange():
                self._vec[i-si].add_comp(the_new_basis)[[j]] = \
                                                   inv_transf.comp(self)[[j,i]]
        # Components of the old dual-basis elements in the new cobasis:
        for i in fmodule.irange():
            for j in fmodule.irange():
                self._dual_basis._vec[i-si].add_comp(the_new_basis)[[j]] = \
                                                       transf.comp(self)[[i,j]]
        # The automorphism and its inverse are added to the module's dictionary
        # of changes of bases:
        fmodule._basis_changes[(self, the_new_basis)] = transf
        fmodule._basis_changes[(the_new_basis, self)] = inv_transf

        return the_new_basis
