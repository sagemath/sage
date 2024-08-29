r"""
Free module morphisms

The class :class:`FiniteRankFreeModuleMorphism` implements homomorphisms
between two free modules of finite rank over the same commutative ring.

The subclass :class:`FiniteRankFreeModuleEndomorphism` implements the
special case of endomorphisms.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Matthias Koeppe (2024): add subclass :class:`FiniteRankFreeModuleEndomorphism`

REFERENCES:

- Chap. 13, 14 of R. Godement : *Algebra* [God1968]_
- Chap. 3 of S. Lang : *Algebra* [Lan2002]_
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

from typing import TYPE_CHECKING

from sage.categories.morphism import Morphism
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer import Integer

if TYPE_CHECKING:
    from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement


class FiniteRankFreeModuleMorphism(Morphism):
    r"""
    Homomorphism between free modules of finite rank over a commutative ring.

    An instance of this class is a homomorphism

    .. MATH::

        \phi:\ M \longrightarrow N,

    where `M` and `N` are two free modules of finite rank over the same
    commutative ring `R`.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.tensor.modules.free_module_homset.FreeModuleHomset`.

    For the special case of endomorphisms (`M=N`), use the subclass
    :class:`FiniteRankFreeModuleEndomorphism`.

    INPUT:

    - ``parent`` -- Hom-set Hom(M,N) to which the homomorphism belongs
    - ``matrix_rep`` -- matrix representation of the homomorphism with
      respect to the bases ``bases``; this entry can actually
      be any material from which a matrix of size rank(N)*rank(M) of
      elements of `R` can be constructed; the *columns* of the matrix give
      the images of the basis of `M` (see the convention in the example below)
    - ``bases`` -- (default: ``None``) pair (basis_M, basis_N) defining the
      matrix representation, basis_M being a basis of module `M` and
      basis_N a basis of module `N` ; if None the pair formed by the
      default bases of each module is assumed.
    - ``name`` -- (default: ``None``) string; name given to the homomorphism
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      homomorphism. If ``None``, ``name`` will be used.

    EXAMPLES:

    A homomorphism between two free modules over `\ZZ` is constructed
    as an element of the corresponding hom-set, by means of the function
    ``__call__``::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
        sage: e = M.basis('e') ; f = N.basis('f')
        sage: H = Hom(M,N) ; H
        Set of Morphisms from Rank-3 free module M over the Integer Ring
         to Rank-2 free module N over the Integer Ring
         in Category of finite dimensional modules over Integer Ring
        sage: phi = H([[2,-1,3], [1,0,-4]], name='phi', latex_name=r'\phi') ; phi
        Generic morphism:
          From: Rank-3 free module M over the Integer Ring
          To:   Rank-2 free module N over the Integer Ring

    Since no bases have been specified in the argument list, the provided
    matrix is relative to the default bases of modules M and N, so that
    the above is equivalent to::

        sage: phi = H([[2,-1,3], [1,0,-4]], bases=(e,f), name='phi',
        ....:         latex_name=r'\phi') ; phi
        Generic morphism:
          From: Rank-3 free module M over the Integer Ring
          To:   Rank-2 free module N over the Integer Ring

    An alternative way to construct a homomorphism is to call the method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.hom`
    on the domain::

        sage: phi = M.hom(N, [[2,-1,3], [1,0,-4]], bases=(e,f), name='phi',
        ....:             latex_name=r'\phi') ; phi
        Generic morphism:
          From: Rank-3 free module M over the Integer Ring
          To:   Rank-2 free module N over the Integer Ring

    The parent of a homomorphism is of course the corresponding hom-set::

        sage: phi.parent() is H
        True
        sage: phi.parent() is Hom(M,N)
        True

    Due to Sage's category scheme, the actual class of the homomorphism ``phi``
    is a derived class of :class:`FiniteRankFreeModuleMorphism`::

        sage: type(phi)
        <class 'sage.tensor.modules.free_module_homset.FreeModuleHomset_with_category_with_equality_by_id.element_class'>
        sage: isinstance(phi, sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism)
        True

    The domain and codomain of the homomorphism are returned respectively by
    the methods ``domain()`` and ``codomain()``, which are implemented as
    Sage's constant functions::

        sage: phi.domain()
        Rank-3 free module M over the Integer Ring
        sage: phi.codomain()
        Rank-2 free module N over the Integer Ring
        sage: type(phi.domain)
        <... 'sage.misc.constant_function.ConstantFunction'>

    The matrix of the homomorphism with respect to a pair of bases is
    returned by the method :meth:`matrix`::

        sage: phi.matrix(e,f)
        [ 2 -1  3]
        [ 1  0 -4]

    The convention is that the columns of this matrix give the components of
    the images of the elements of basis ``e`` w.r.t basis ``f``::

        sage: phi(e[0]).display()
        phi(e_0) = 2 f_0 + f_1
        sage: phi(e[1]).display()
        phi(e_1) = -f_0
        sage: phi(e[2]).display()
        phi(e_2) = 3 f_0 - 4 f_1

    Test of the module homomorphism laws::

        sage: phi(M.zero()) == N.zero()
        True
        sage: u = M([1,2,3], basis=e, name='u') ; u.display()
        u = e_0 + 2 e_1 + 3 e_2
        sage: v = M([-2,1,4], basis=e, name='v') ; v.display()
        v = -2 e_0 + e_1 + 4 e_2
        sage: phi(u).display()
        phi(u) = 9 f_0 - 11 f_1
        sage: phi(v).display()
        phi(v) = 7 f_0 - 18 f_1
        sage: phi(3*u + v).display()
        34 f_0 - 51 f_1
        sage: phi(3*u + v) == 3*phi(u) + phi(v)
        True
    """
    def __init__(self, parent, matrix_rep, bases=None, name=None,
                 latex_name=None, is_identity=False):
        r"""
        TESTS:

        Generic homomorphism::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: from sage.tensor.modules.free_module_morphism import FiniteRankFreeModuleMorphism
            sage: phi = FiniteRankFreeModuleMorphism(Hom(M,N), [[1,0,-3], [2,1,4]],
            ....:                                    name='phi', latex_name=r'\phi')
            sage: phi
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: phi.matrix(e,f)
            [ 1  0 -3]
            [ 2  1  4]
            sage: latex(phi)
            \phi
        """
        if is_identity:
            raise TypeError('use the subclass FiniteRankFreeModuleEndomorphism for the identity morphis')
        from sage.matrix.constructor import matrix
        from sage.misc.constant_function import ConstantFunction
        Morphism.__init__(self, parent)
        fmodule1 = parent.domain()
        fmodule2 = parent.codomain()
        if bases is None:
            def_basis1 = fmodule1.default_basis()
            if def_basis1 is None:
                raise ValueError("the {} has no default ".format(fmodule1) +
                                 "basis")
            def_basis2 = fmodule2.default_basis()
            if def_basis2 is None:
                raise ValueError("the {} has no default ".format(fmodule2) +
                                 "basis")
            bases = (def_basis1, def_basis2)
        else:
            bases = tuple(bases)  # insures bases is a tuple
            if len(bases) != 2:
                raise TypeError("the argument bases must contain 2 bases")
            if bases[0] not in fmodule1.bases():
                raise TypeError("{} is not a basis on the {}".format(bases[0],
                                                                     fmodule1))
            if bases[1] not in fmodule2.bases():
                raise TypeError("{} is not a basis on the {}".format(bases[1],
                                                                     fmodule2))
        ring = parent.base_ring()
        n1 = fmodule1.rank()
        n2 = fmodule2.rank()
        # Construction of a generic morphism
        self._is_identity = False
        if isinstance(matrix_rep, ConstantFunction):
            # the zero morphism
            if matrix_rep().is_zero():
                matrix_rep = 0
        self._matrices = {bases: matrix(ring, n2, n1, matrix_rep)}
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

    #
    # SageObject methods
    #

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\Phi')
            sage: phi._latex_()
            '\\Phi'
            sage: latex(phi)  # indirect doctest
            \Phi

        ::

            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='F')
            sage: phi._latex_()
            'F'

        ::

            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi._latex_()
            '\\mbox{Generic morphism:\n  From: Rank-3 free module M over the Integer Ring\n  To:   Rank-2 free module N over the Integer Ring}'
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
            return self._latex_name

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a free module morphism (or 0)

        OUTPUT: ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\phi')
            sage: psi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.__eq__(psi)
            True
            sage: phi == psi
            True
            sage: phi.__eq__(phi)
            True
            sage: phi.__eq__(+phi)
            True
            sage: psi = M.hom(N, [[1,1,0], [4,1,3]])
            sage: phi.__eq__(psi)
            False
            sage: phi.__eq__(-phi)
            False

        Comparison of homomorphisms defined on different bases::

            sage: a = M.automorphism() ; a[0,2], a[1,0], a[2,1] = 1, -1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: psi = M.hom(N, [[-2,0,-1], [-1,-2, 5]], bases=(ep,f))
            sage: phi.__eq__(psi)
            True
            sage: phi.matrix(e,f) == psi.matrix(e,f)  # check
            True

        Comparison of homomorphisms having the same matrix but defined on
        different modules::

            sage: N1 = FiniteRankFreeModule(ZZ, 2, name='N1')
            sage: f1 = N1.basis('f')
            sage: phi1 = M.hom(N1, [[-1,2,0], [5,1,2]])
            sage: phi.matrix() == phi1.matrix() # same matrix in the default bases
            True
            sage: phi.__eq__(phi1)
            False

        Comparison to zero::

            sage: phi.__eq__(0)
            False
            sage: phi = M.hom(N, 0)
            sage: phi.__eq__(0)
            True
            sage: phi == 0
            True
            sage: phi.__eq__(Hom(M,N).zero())
            True
        """
        if isinstance(other, (int, Integer)): # other should be 0
            if other == 0:
                return self.is_zero()
            else:
                return False
        elif not isinstance(other, FiniteRankFreeModuleMorphism):
            return False
        elif self.parent() != other.parent():
            return False
        else:
            bases = self._common_bases(other)
            if bases is None:
                raise ValueError("no common pair of bases has been found to " +
                                 "compare {} and {}".format(self, other))
            return bool( self.matrix(*bases) == other.matrix(*bases) )

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- a free module morphism (or 0)

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\phi')
            sage: psi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.__ne__(psi)
            False
            sage: psi = M.hom(N, [[1,1,0], [4,1,3]])
            sage: phi.__ne__(psi)
            True
            sage: phi != psi
            True
            sage: phi.__ne__('junk')
            True
            sage: Hom(M,N).zero().__ne__(0)
            False
        """
        return not self == other

    #
    # Required module methods
    #

    def __bool__(self):
        r"""
        Return ``True`` if ``self`` is nonzero and ``False`` otherwise.

        This method is called by self.is_zero().

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[2,-1,3], [1,0,-4]])
            sage: bool(phi)
            True
            sage: phi.is_zero() # indirect doctest
            False
            sage: phi = M.hom(N, 0)
            sage: bool(phi)
            False
            sage: phi.is_zero() # indirect doctest
            True
            sage: bool(Hom(M,N).zero())
            False
        """
        # Some matrix representation is picked at random:
        matrix_rep = next(iter(self._matrices.values()))
        return not matrix_rep.is_zero()

    def _add_(self, other):
        r"""
        Homomorphism addition.

        INPUT:

        - ``other`` -- a free module morphism (same parent as ``self``)

        OUTPUT: the homomorphism resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\phi')
            sage: psi = M.hom(N, [[1,1,0], [4,1,3]])
            sage: s = phi._add_(psi) ; s
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: s.matrix(e,f)
            [0 3 0]
            [9 2 5]
            sage: s.matrix(e,f) == phi.matrix(e,f) + psi.matrix(e,f)  # check
            True
            sage: s == phi + psi  # indirect doctest
            True

        Addition of homomorphisms defined on different bases::

            sage: a = M.automorphism() ; a[0,2], a[1,0], a[2,1] = 1, -1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: b = N.automorphism() ; b[0,1], b[1,0] = -1, 1
            sage: fp = f.new_basis(b, 'fp', latex_symbol="f'")
            sage: psi = M.hom(N, [[-2,0,-1], [-1,-2, 5]], bases=(ep,fp))
            sage: s = phi._add_(psi) ; s
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: s.matrix(e,f)
            [-6  1 -2]
            [ 4  3  2]
            sage: s.matrix(e,f) == phi.matrix(e,f) + psi.matrix(e,f)  # check
            True
            sage: s == phi + psi  # indirect doctest
            True

        Other tests::

            sage: phi._add_(Hom(M,N).zero()) == phi
            True
        """
        # No need for consistency checks since self and other are guaranteed
        # to have the same parents
        bases = self._common_bases(other)
        if bases is None:
            raise ValueError("no common pair of bases has been found to " +
                             "add {} and {}".format(self, other))
        # Addition at the matrix level:
        resu_mat = self._matrices[bases] + other._matrices[bases]
        if self._name is not None and other._name is not None:
            resu_name = self._name + '+' + other._name
        else:
            resu_name = None
        if self._latex_name is not None and other._latex_name is not None:
            resu_latex_name = self._latex_name + '+' + other._latex_name
        else:
            resu_latex_name = None
        return self.__class__(self.parent(), resu_mat, bases=bases,
                              name=resu_name, latex_name=resu_latex_name)

    def _sub_(self, other):
        r"""
        Homomorphism subtraction.

        INPUT:

        - ``other`` -- a free module morphism (same parent as ``self``)

        OUTPUT:

        - the homomorphism resulting from the subtraction of ``other`` from
          ``self``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\phi')
            sage: psi = M.hom(N, [[1,1,0], [4,1,3]])
            sage: s = phi._sub_(psi) ; s
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: s.matrix(e,f)
            [-2  1  0]
            [ 1  0 -1]
            sage: s.matrix(e,f) == phi.matrix(e,f) - psi.matrix(e,f)  # check
            True
            sage: s == phi - psi  # indirect doctest
            True

        Subtraction of homomorphisms defined on different bases::

            sage: a = M.automorphism() ; a[0,2], a[1,0], a[2,1] = 1, -1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: b = N.automorphism() ; b[0,1], b[1,0] = -1, 1
            sage: fp = f.new_basis(b, 'fp', latex_symbol="f'")
            sage: psi = M.hom(N, [[-2,0,-1], [-1,-2, 5]], bases=(ep,fp))
            sage: s = phi._sub_(psi) ; s
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: s.matrix(e,f)
            [ 4  3  2]
            [ 6 -1  2]
            sage: s.matrix(e,f) == phi.matrix(e,f) - psi.matrix(e,f)  # check
            True
            sage: s == phi - psi  # indirect doctest
            True

        Other tests::

            sage: phi._sub_(Hom(M,N).zero()) == phi
            True
            sage: Hom(M,N).zero()._sub_(phi) == -phi
            True
            sage: phi._sub_(phi).is_zero()
            True
        """
        # No need for consistency checks since self and other are guaranteed
        # to have the same parents
        bases = self._common_bases(other)
        if bases is None:
            raise ValueError("no common pair of bases has been found to " +
                             "subtract {} from {}".format(other, self))
        # Subtraction at the matrix level:
        resu_mat = self._matrices[bases] - other._matrices[bases]
        if self._name is not None and other._name is not None:
            resu_name = self._name + '-' + other._name
        else:
            resu_name = None
        if self._latex_name is not None and other._latex_name is not None:
            resu_latex_name = self._latex_name + '-' + other._latex_name
        else:
            resu_latex_name = None
        return self.__class__(self.parent(), resu_mat, bases=bases,
                              name=resu_name, latex_name=resu_latex_name)

    def _lmul_(self, scalar):
        r"""
        Multiplication by ``scalar``.

        INPUT:

        - ``scalar`` -- element of the ring over which the parent of ``self``
          is a module

        OUTPUT: the homomorphism resulting from the multiplication of ``self``
        by ``scalar``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\phi')
            sage: s = phi._lmul_(7) ; s
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: s.matrix(e,f)
            [-7 14  0]
            [35  7 14]
            sage: s == 7 * phi
            True
        """
        resu = self.__class__(self.parent(), 0)  # 0 = provisory value
        for bases, mat in self._matrices.items():
            resu._matrices[bases] = scalar * mat
        return resu

    #
    #  Other module methods
    #

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT: an exact copy of ``self``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\phi')
            sage: s = phi.__pos__() ; s
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: s == +phi
            True
            sage: s == phi
            True
            sage: s is phi
            False
        """
        resu = self.__class__(self.parent(), 0, is_identity=self._is_identity)
                                           # 0 = provisory value
        for bases, mat in self._matrices.items():
            resu._matrices[bases] = +mat
        if self._name is not None:
            resu._name = '+' + self._name
        if self._latex_name is not None:
            resu._latex_name = '+' + self._latex_name
        return resu

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT: the homomorphism `-f`, where `f` is ``self``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi',
            ....:             latex_name=r'\phi')
            sage: s = phi.__neg__() ; s
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring
            sage: s == -phi
            True
            sage: s.matrix()
            [ 1 -2  0]
            [-5 -1 -2]
            sage: s.matrix() == -phi.matrix()
            True
        """
        resu = self.__class__(self.parent(), 0)  # 0 = provisory value
        for bases, mat in self._matrices.items():
            resu._matrices[bases] = -mat
        if self._name is not None:
            resu._name = '-' + self._name
        if self._latex_name is not None:
            resu._latex_name = '-' + self._latex_name
        return resu

    #
    # Map methods
    #

    def _call_(
        self, element: FiniteRankFreeModuleElement
    ) -> FiniteRankFreeModuleElement:
        r"""
        Action of the homomorphism ``self`` on some free module element

        INPUT:

        - ``element`` -- element of the domain of ``self``

        OUTPUT: the image of ``element`` by ``self``

        EXAMPLES:

        Images of a homomorphism between two `\ZZ`-modules::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi', latex_name=r'\phi')
            sage: v = M([1,2,3], basis=e, name='v')
            sage: w = phi(v) ; w
            Element phi(v) of the Rank-2 free module N over the Integer Ring
            sage: w.display()
            phi(v) = 3 f_0 + 13 f_1

        TESTS::

            sage: all(w[i] == sum(phi.matrix()[i,j]*v[j] for j in range(3)) for i in range(2))
            True
            sage: phi.matrix(e,f)
            [-1  2  0]
            [ 5  1  2]
            sage: phi(e[0]).display()
            phi(e_0) = -f_0 + 5 f_1
            sage: phi(e[1]).display()
            phi(e_1) = 2 f_0 + f_1
            sage: phi(e[2]).display()
            phi(e_2) = 2 f_1

        Image of an element that is not defined on the default basis::

            sage: a = M.automorphism()
            sage: a[0,2], a[1,0], a[2,1] = 1, -1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: v = M([1,2,3], basis=ep, name='v')
            sage: w = phi(v) ; w
            Element phi(v) of the Rank-2 free module N over the Integer Ring
            sage: w.display()
            phi(v) = -5 f_0 + 10 f_1
            sage: all(w[i] == sum(phi.matrix(ep,f)[i,j]*v[ep,j] for j in range(3)) for i in range(2))
            True

        Check of homomorphism properties::

            sage: phi(M.zero()) == N.zero()
            True
        """
        if self._is_identity:
            return element
        dom = self.parent().domain()
        sindex = dom._sindex
        codom = self.parent().codomain()
        basis_codom = codom.default_basis()
        # Search for a common basis to compute the result
        for basis in element._components:
            try:
                self.matrix(basis, basis_codom)
                basis_dom = basis
                break
            except ValueError:
                continue
        else:
            raise ValueError("no common basis found to evaluate the image " +
                             "of {} by {}".format(element,self))
        # Components of the result obtained by matrix multiplication
        mat = self.matrix(basis_dom, basis_codom)
        vcomp = element._components[basis_dom]
        tresu = []
        for i in range(codom.rank()):
            s = 0
            for j in range(dom.rank()):
                s += mat[i,j] * vcomp[[j+sindex]]
            tresu.append(s)
        # Name of the result
        if self._name is not None and element._name is not None:
            resu_name = self._name + '(' + element._name + ')'
        else:
            resu_name = None
        if self._latex_name is not None and element._latex_name is not None:
            resu_latex_name = self._latex_name + r'\left(' + \
                              element._latex_name + r'\right)'
        else:
            resu_latex_name = None
        # Creation of the result
        return codom(tresu, basis=basis_codom, name=resu_name,
                     latex_name=resu_latex_name)

    def is_injective(self):
        r"""
        Determine whether ``self`` is injective.

        OUTPUT:

        - ``True`` if ``self`` is an injective homomorphism and ``False``
          otherwise

        EXAMPLES:

        Homomorphisms between two `\ZZ`-modules::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.matrix(e,f)
            [-1  2  0]
            [ 5  1  2]
            sage: phi.is_injective()
            False

        Indeed, phi has a non trivial kernel::

            sage: phi(4*e[0] + 2*e[1] - 11*e[2]).display()
            0

        An injective homomorphism::

            sage: psi = N.hom(M, [[1,-1], [0,3], [4,-5]])
            sage: psi.matrix(f,e)
            [ 1 -1]
            [ 0  3]
            [ 4 -5]
            sage: psi.is_injective()
            True

        Of course, the identity endomorphism is injective::

            sage: End(M).one().is_injective()
            True
            sage: End(N).one().is_injective()
            True
        """
        # Some matrix representation is picked at random:
        matrix_rep = next(iter(self._matrices.values()))
        return not matrix_rep.right_kernel().rank()

    def is_surjective(self):
        r"""
        Determine whether ``self`` is surjective.

        OUTPUT:

        - ``True`` if ``self`` is a surjective homomorphism and ``False``
          otherwise

        EXAMPLES:

        This method has not been implemented yet::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.is_surjective()
            Traceback (most recent call last):
            ...
            NotImplementedError: FiniteRankFreeModuleMorphism.is_surjective()
             has not been implemented yet

        except for the identity endomorphism (!)::

            sage: End(M).one().is_surjective()
            True
            sage: End(N).one().is_surjective()
            True
        """
        if self._is_identity:
            return True
        raise NotImplementedError(
                              "FiniteRankFreeModuleMorphism.is_surjective() " +
                              "has not been implemented yet")
    #
    # Morphism methods
    #

    def is_identity(self):
        r"""
        Check whether ``self`` is the identity morphism.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: phi = M.endomorphism([[1,0], [0,1]])
            sage: phi.is_identity()
            True
            sage: (phi+phi).is_identity()
            False
            sage: End(M).zero().is_identity()
            False
            sage: a = M.automorphism() ; a[0,1], a[1,0] = 1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,0], [0,1]], basis=ep)
            sage: phi.is_identity()
            True

        Example illustrating that the identity can be constructed from a
        matrix that is not the identity one, provided that it is relative to
        different bases::

            sage: phi = M.hom(M, [[0,1], [-1,0]], bases=(ep,e))
            sage: phi.is_identity()
            True

        Of course, if we ask for the matrix in a single basis, it is the
        identity matrix::

            sage: phi.matrix(e)
            [1 0]
            [0 1]
            sage: phi.matrix(ep)
            [1 0]
            [0 1]
        """
        if self._is_identity:
            return True
        # The identity must be an endomorphism:
        fmodule = self.domain()
        if fmodule != self.codomain():
            return False
        # Some basis in which ``self`` has a representation is picked at
        # random and the test is performed on the images of the basis
        # elements:
        basis = list(self._matrices)[0][0]
        for i in fmodule.irange():
            if self(basis[i]) != basis[i]:
                return False
        self._is_identity = True
        return True

    #
    # End of Morphism methods
    #

    def _modules_and_bases(self, basis1=None, basis2=None):
        r"""
        Return domain, codomain, domain basis, and codomain basis.

        This method implements default argument handling for methods
        :meth:`matrix` and :meth:`display`.

        INPUT:

        - ``basis1`` -- (default: ``None``) basis of the domain of ``self``; if
          none is provided, the domain's default basis is assumed
        - ``basis2`` -- (default: ``None``) basis of the codomain of ``self``;
          if none is provided, ``basis2`` is set to ``basis1`` if ``self`` is
          an endomorphism, otherwise, ``basis2`` is set to the codomain's
          default basis.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e'); f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi._modules_and_bases()
            (Rank-3 free module M over the Integer Ring,
             Rank-2 free module N over the Integer Ring,
             Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1) on the Rank-2 free module N over the Integer Ring)
            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi._modules_and_bases(ep)
            (Rank-3 free module M over the Integer Ring,
             Rank-2 free module N over the Integer Ring,
             Basis (ep_0,ep_1,ep_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1) on the Rank-2 free module N over the Integer Ring)
        """
        fmodule1 = self.domain()
        fmodule2 = self.codomain()
        if basis1 is None:
            basis1 = fmodule1.default_basis()
        elif basis1 not in fmodule1.bases():
            raise TypeError(str(basis1) + " is not a basis on the " +
                            str(fmodule1) + ".")
        if basis2 is None:
            if self.is_endomorphism():
                basis2 = basis1
            else:
                basis2 = fmodule2.default_basis()
        elif basis2 not in fmodule2.bases():
            raise TypeError(str(basis2) + " is not a basis on the " +
                            str(fmodule2) + ".")
        return fmodule1, fmodule2, basis1, basis2

    def matrix(self, basis1=None, basis2=None):
        r"""
        Return the matrix of ``self`` w.r.t to a pair of bases.

        If the matrix is not known already, it is computed from the matrix in
        another pair of bases by means of the change-of-basis formula.

        INPUT:

        - ``basis1`` -- (default: ``None``) basis of the domain of ``self``; if
          none is provided, the domain's default basis is assumed
        - ``basis2`` -- (default: ``None``) basis of the codomain of ``self``;
          if none is provided, ``basis2`` is set to ``basis1`` if ``self`` is
          an endomorphism, otherwise, ``basis2`` is set to the codomain's
          default basis.

        OUTPUT:

        - the matrix representing the homomorphism ``self`` w.r.t
          to bases ``basis1`` and ``basis2``; more precisely, the columns of
          this matrix are formed by the components w.r.t. ``basis2`` of
          the images of the elements of ``basis1``.

        EXAMPLES:

        Matrix of a homomorphism between two `\ZZ`-modules::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.matrix()     # default bases
            [-1  2  0]
            [ 5  1  2]
            sage: phi.matrix(e, f)  # given bases
            [-1  2  0]
            [ 5  1  2]
            sage: type(phi.matrix())                                                    # needs sage.libs.flint
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

        Matrix in bases different from those in which the homomorphism has
        been defined::

            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: b = N.automorphism(matrix=[[3,5],[4,7]], basis=f)
            sage: fp = f.new_basis(b, 'fp', latex_symbol="f'")
            sage: phi.matrix(ep, fp)
            [ 32  -1 -12]
            [-19   1   8]

        Check of the change-of-basis formula::

            sage: phi.matrix(ep, fp) == (b^(-1)).matrix(f) * phi.matrix(e,f) * a.matrix(e)
            True

        Single change of basis::

            sage: phi.matrix(ep, f)
            [ 1  2  4]
            [-5  3  8]
            sage: phi.matrix(ep,f) == phi.matrix(e,f) * a.matrix(e)
            True
            sage: phi.matrix(e, fp)
            [-32   9 -10]
            [ 19  -5   6]
            sage: phi.matrix(e, fp) == (b^(-1)).matrix(f) * phi.matrix(e,f)
            True

        Matrix of an endomorphism::

            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi.matrix(ep)
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: phi.matrix(ep,ep)  # same as above
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: phi.matrix()  # matrix w.r.t to the module's default basis
            [  1  -3   1]
            [-18  39 -18]
            [-25  54 -25]
        """
        from sage.matrix.constructor import matrix
        fmodule1, fmodule2, basis1, basis2 = self._modules_and_bases(basis1, basis2)
        if (basis1, basis2) not in self._matrices:
            if self._is_identity:
                # The identity endomorphism
                # -------------------------
                if basis1 == basis2:
                    # the matrix is the identity matrix:
                    ring = fmodule1.base_ring()
                    zero = ring.zero()
                    one = ring.one()
                    size = fmodule1.rank()
                    mat = []
                    for i in range(size):
                        row = [zero]*size
                        row[i] = one
                        mat.append(row)
                else:
                    # the matrix is the change-of-basis matrix:
                    change = fmodule1.change_of_basis(basis1, basis2)
                    mat = [[change[[i,j]] for j in fmodule1.irange()]
                                                    for i in fmodule1.irange()]
                self._matrices[(basis1, basis2)] = matrix(mat)
            else:
                # Generic homomorphism
                # --------------------
                b1_list = [bases[0] for bases in self._matrices]
                b2_list = [bases[1] for bases in self._matrices]
                if basis1 in b1_list:
                    for b2 in b2_list:
                        if (basis2, b2) in fmodule2._basis_changes:
                            nb2 = b2
                            break
                    else:
                        raise ValueError("no start basis could be found for " +
                                        "applying the change-of-basis formula")
                    change2 = fmodule2._basis_changes[(basis2, nb2)]
                    mat2 = matrix( [[change2[[i,j]] for j in fmodule2.irange()]
                                                  for i in fmodule2.irange()] )
                    self._matrices[(basis1, basis2)] = \
                                            mat2 * self._matrices[(basis1,nb2)]
                elif basis2 in b2_list:
                    for b1 in b1_list:
                        if (b1, basis1) in fmodule1._basis_changes:
                            nb1 = b1
                            break
                    else:
                        raise ValueError("no start basis could be found for " +
                                        "applying the change-of-basis formula")
                    change1 = fmodule1._basis_changes[(nb1, basis1)]
                    mat1 = matrix( [[change1[[i,j]] for j in fmodule1.irange()]
                                                  for i in fmodule1.irange()] )
                    self._matrices[(basis1, basis2)] = \
                                            self._matrices[(nb1,basis2)] * mat1
                else: # most general change-of-basis formula
                    for (b1, b2) in self._matrices:
                        if (b1, basis1) in fmodule1._basis_changes and \
                           (basis2, b2) in fmodule2._basis_changes:
                            nb1, nb2 = b1, b2
                            break
                    else:
                        raise ValueError("no start basis could be found for " +
                                        "applying the change-of-basis formula")
                    change1 = fmodule1._basis_changes[(nb1, basis1)]
                    change2 = fmodule2._basis_changes[(basis2, nb2)]
                    mat1 = matrix( [[change1[[i,j]] for j in fmodule1.irange()]
                                                  for i in fmodule1.irange()] )
                    mat2 = matrix( [[change2[[i,j]] for j in fmodule2.irange()]
                                                  for i in fmodule2.irange()] )
                    self._matrices[(basis1, basis2)] = \
                                        mat2 * self._matrices[(nb1,nb2)] * mat1
        return self._matrices[(basis1, basis2)]

    def _common_bases(self, other):
        r"""
        Return a pair of bases in which ``self`` and ``other`` have a known
        matrix representation.

        INPUT:

        - ``other`` -- another homomorphism in the same hom-set

        OUTPUT:

        - a pair of bases in which ``self`` and ``other`` have a known
          matrix representation.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: psi = M.hom(N, [[1,1,0], [4,1,3]])
            sage: phi._common_bases(psi) # matrices of phi and psi both defined on (e,f)
            (Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1) on the Rank-2 free module N over the Integer Ring)
            sage: a = M.automorphism() ; a[0,2], a[1,0], a[2,1] = 1, -1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: psi = M.hom(N, [[1,1,0], [4,1,3]], bases=(ep,f))
            sage: phi._common_bases(psi) # matrix of psi w.r.t. (e,f) computed
            (Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1) on the Rank-2 free module N over the Integer Ring)
            sage: psi = M.hom(N, [[1,1,0], [4,1,3]], bases=(ep,f))
            sage: psi._common_bases(phi) # matrix of phi w.r.t. (ep,f) computed
            (Basis (ep_0,ep_1,ep_2) on the Rank-3 free module M over the Integer Ring,
             Basis (f_0,f_1) on the Rank-2 free module N over the Integer Ring)
        """
        resu = None
        for bases in self._matrices:
            try:
                other.matrix(*bases)
                resu = bases
                break
            except ValueError:
                continue
        if resu is None:
            for bases in other._matrices:
                try:
                    self.matrix(*bases)
                    resu = bases
                    break
                except ValueError:
                    continue
        return resu

    def display(self, basis1=None, basis2=None):
        r"""
        Display ``self`` as a matrix w.r.t to a pair of bases.

        If the matrix is not known already, it is computed from the matrix in
        another pair of bases by means of the change-of-basis formula.

        INPUT:

        - ``basis1`` -- (default: ``None``) basis of the domain of ``self``; if
          none is provided, the domain's default basis is assumed
        - ``basis2`` -- (default: ``None``) basis of the codomain of ``self``;
          if none is provided, ``basis2`` is set to ``basis1`` if ``self`` is
          an endomorphism, otherwise, ``basis2`` is set to the codomain's
          default basis.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e'); f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.display()     # default bases
                e_0 e_1 e_2
            f_0⎛ -1   2   0⎞
            f_1⎝  5   1   2⎠
            sage: phi.display(e, f)  # given bases
                e_0 e_1 e_2
            f_0⎛ -1   2   0⎞
            f_1⎝  5   1   2⎠

        Matrix of an endomorphism::

            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi.display(ep)
                 ep_0 ep_1 ep_2
            ep_0⎛   1    2    3⎞
            ep_1⎜   4    5    6⎟
            ep_2⎝   7    8    9⎠
            sage: phi.display(ep, ep)  # same as above
                 ep_0 ep_1 ep_2
            ep_0⎛   1    2    3⎞
            ep_1⎜   4    5    6⎟
            ep_2⎝   7    8    9⎠
            sage: phi.display()  # matrix w.r.t to the module's default basis
                e_0 e_1 e_2
            e_0⎛  1  -3   1⎞
            e_1⎜-18  39 -18⎟
            e_2⎝-25  54 -25⎠
        """
        from sage.misc.latex import latex
        from .format_utilities import is_atomic, FormattedExpansion
        fmodule1, fmodule2, basis1, basis2 = self._modules_and_bases(basis1, basis2)
        matrix = self.matrix(basis1, basis2)
        if all(element._name for element in basis1):
            basis1_names = [element._name for element in basis1]
        else:
            basis1_names = None
        if all(element._name for element in basis2):
            basis2_names = [element._name for element in basis2]
        else:
            basis2_names = None
        resu_txt = matrix.str(unicode=True,
                              top_border=basis1_names,
                              left_border=basis2_names)
        resu_latex = latex(matrix)
        return FormattedExpansion(resu_txt, resu_latex)


class FiniteRankFreeModuleEndomorphism(FiniteRankFreeModuleMorphism):
    r"""
    Endomorphism of a free module of finite rank over a commutative ring.

    An instance of this class is an endomorphism

    .. MATH::

        \phi:\ M \longrightarrow M,

    where `M` is a free module of finite rank over a commutative ring `R`.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.tensor.modules.free_module_homset.FreeModuleEndset`.

    INPUT:

    - ``parent`` -- Hom-set Hom(M,M) to which the endomorphism belongs
    - ``matrix_rep`` -- matrix representation of the endomorphism with
      respect to the basis ``bases``; this entry can actually
      be any material from which a matrix of size rank(N)*rank(M) of
      elements of `R` can be constructed; the *columns* of the matrix give
      the images of the basis of `M` (see the convention in the example below)
    - ``bases`` -- (default: ``None``) pair ``(basis_domain, basis_codomain)``
      defining the matrix representation, ``basis_domain`` and ``basis_codomain``
      being two bases (typically the same) of the same module `M`; if ``None``,
      the default basis of `M` is used for both.
    - ``name`` -- (default: ``None``) string; name given to the endomorphism
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      endomorphism. If ``None``, ``name`` will be used.
    - ``is_identity`` -- boolean (default: ``False``); determines whether the
      constructed object is the identity endomorphism. If set to ``True``,
      then the entry ``matrix_rep`` is not used.

    EXAMPLES:

    The identity endomorphism::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: Id = End(M).one(); Id
        Identity endomorphism of Rank-3 free module M over the Integer Ring
        sage: Id.parent()
        Set of Morphisms from Rank-3 free module M over the Integer Ring
         to Rank-3 free module M over the Integer Ring
         in Category of finite dimensional modules over Integer Ring
        sage: Id.parent() is End(M)
        True

    The matrix of ``Id`` with respect to the basis ``e`` is of course the identity
    matrix::

        sage: Id.matrix(e)
        [1 0 0]
        [0 1 0]
        [0 0 1]

    The identity acting on a module element::

        sage: v = M([-2,1,4], basis=e, name='v'); v.display()
        v = -2 e_0 + e_1 + 4 e_2
        sage: Id(v) is v
        True
    """
    def __init__(self, parent, matrix_rep, bases=None, name=None,
                 latex_name=None, is_identity=False):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_morphism import FiniteRankFreeModuleEndomorphism
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')

        Generic endomorphism::

            sage: phi = FiniteRankFreeModuleEndomorphism(End(M),
            ....:                                        [[1,0,-3], [2,1,4], [7,8,9]],
            ....:                                        name='phi', latex_name=r'\phi')
            sage: phi
            Generic endomorphism of Rank-3 free module M over the Integer Ring

        Identity endomorphism::

            sage: phi = FiniteRankFreeModuleEndomorphism(End(M), 'whatever',
            ....:                                        is_identity=True)
            sage: phi
            Identity endomorphism of Rank-3 free module M over the Integer Ring
            sage: phi.matrix(e)
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: latex(phi)
            \mathrm{Id}
        """
        from sage.matrix.special import identity_matrix
        from sage.misc.constant_function import ConstantFunction

        fmodule = parent.domain()
        if bases is None:
            def_basis = fmodule.default_basis()
            if def_basis is None:
                raise ValueError("the {} has no default ".format(fmodule) +
                                 "basis")
            bases = (def_basis, def_basis)
        else:
            bases = tuple(bases)  # insures bases is a tuple
            if len(bases) != 2:
                raise TypeError("the argument bases must contain 2 bases")
            if bases[0] not in fmodule.bases():
                raise TypeError("{} is not a basis on the {}".format(bases[0],
                                                                     fmodule))
            if bases[1] not in fmodule.bases():
                raise TypeError("{} is not a basis on the {}".format(bases[1],
                                                                     fmodule))
        if not is_identity:
            # Construction of a generic endomorphism
            if isinstance(matrix_rep, ConstantFunction):
                # the zero morphism
                if matrix_rep().is_zero():
                    matrix_rep = 0
            if bases[0] == bases[1] and matrix_rep == 1:
                is_identity = True
        if is_identity:
            # Construction of the identity endomorphism
            if bases[0] != bases[1]:
                raise TypeError("the two bases must coincide for " +
                                "constructing the identity endomorphism.")
            matrix_rep = 1
            if name is None:
                name = 'Id'
            if latex_name is None and name == 'Id':
                latex_name = r'\mathrm{Id}'
        FiniteRankFreeModuleMorphism.__init__(self, parent, matrix_rep, bases,
                                              name, latex_name)
        if is_identity:
            self._is_identity = True
            self._repr_type_str = 'Identity'

    def _some_matrix(self):
        r"""
        Return the matrix of the endomorphism ``self`` w.r.t. some basis.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi._some_matrix()
            [1 2 3]
            [4 5 6]
            [7 8 9]
        """
        for (basis1, basis2), matrix in self._matrices.items():
            if basis1 == basis2:
                return matrix
        for basis in self.domain().bases():
            try:
                return self.matrix(basis)
            except ValueError:
                pass
        return self.matrix()

    @lazy_attribute
    def characteristic_polynomial(self):
        r"""
        Return the characteristic polynomial of ``self``.

        :meth:`characteristic_polynomial` and :meth:`charpoly` are the same method.

        INPUT:

        - ``var`` -- string (default: ``'x'``); a variable name

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi.matrix(e)
            [  1  -3   1]
            [-18  39 -18]
            [-25  54 -25]
            sage: phi.characteristic_polynomial()
            x^3 - 15*x^2 - 18*x
            sage: phi.charpoly()
            x^3 - 15*x^2 - 18*x
            sage: phi.charpoly('T')
            T^3 - 15*T^2 - 18*T
        """
        return self._some_matrix().characteristic_polynomial

    charpoly = characteristic_polynomial

    @lazy_attribute
    def det(self):
        r"""
        Return the determinant of ``self``.

        OUTPUT:

        - element of the base ring of the modules on which ``self`` is defined,
          equal to the determinant of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi.matrix(e)
            [  1  -3   1]
            [-18  39 -18]
            [-25  54 -25]
            sage: phi.det()
            0
            sage: det(phi)
            0
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

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi.matrix(e)
            [  1  -3   1]
            [-18  39 -18]
            [-25  54 -25]
            sage: phi.fcp()                                                             # needs sage.libs.pari
            x * (x^2 - 15*x - 18)
            sage: phi.fcp('T')                                                          # needs sage.libs.pari
            T * (T^2 - 15*T - 18)
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

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi.matrix(e)
            [  1  -3   1]
            [-18  39 -18]
            [-25  54 -25]
            sage: phi.minpoly()                                                         # needs sage.libs.pari
            x^3 - 15*x^2 - 18*x
            sage: phi.minimal_polynomial()                                              # needs sage.libs.pari
            x^3 - 15*x^2 - 18*x
            sage: phi.minimal_polynomial('T')                                           # needs sage.libs.pari
            T^3 - 15*T^2 - 18*T
        """
        return self._some_matrix().minimal_polynomial

    minpoly = minimal_polynomial

    @lazy_attribute
    def trace(self):
        r"""
        Return the trace of ``self``.

        OUTPUT:

        - element of the base ring of the modules on which ``self`` is defined,
          equal to the trace of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(matrix=[[-1,0,0],[0,1,2],[0,1,3]], basis=e)
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,2,3], [4,5,6], [7,8,9]], basis=ep)
            sage: phi.matrix(e)
            [  1  -3   1]
            [-18  39 -18]
            [-25  54 -25]
            sage: phi.trace()
            15
            sage: id = M.identity_map()
            sage: id.trace()
            3
        """
        return self._some_matrix().trace
