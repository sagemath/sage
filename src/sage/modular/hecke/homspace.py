# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
r"""
Hom spaces between Hecke modules
"""
# ****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.categories.homset import HomsetWithBase
from .morphism import HeckeModuleMorphism_matrix
from .module import HeckeModule_generic


def is_HeckeModuleHomspace(x):
    r"""
    Return ``True`` if x is a space of homomorphisms in the category of Hecke modules.

    EXAMPLES::

        sage: M = ModularForms(Gamma0(7), 4)
        sage: sage.modular.hecke.homspace.is_HeckeModuleHomspace(Hom(M, M))
        doctest:warning...
        DeprecationWarning: the function is_HeckeModuleHomspace is deprecated;
        use 'isinstance(..., HeckeModuleHomspace)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        True
        sage: sage.modular.hecke.homspace.is_HeckeModuleHomspace(Hom(M, QQ))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_HeckeModuleHomspace is deprecated; use 'isinstance(..., HeckeModuleHomspace)' instead")
    return isinstance(x, HeckeModuleHomspace)


class HeckeModuleHomspace(HomsetWithBase):
    r"""
    A space of homomorphisms between two objects in the category of Hecke
    modules over a given base ring.
    """
    def __init__(self, X, Y, category=None):
        r"""
        Create the space of homomorphisms between X and Y, which must have the
        same base ring.

        EXAMPLES::

            sage: M = ModularForms(Gamma0(7), 4)
            sage: M.Hom(M)
            Set of Morphisms from ... to ... in Category of Hecke modules over Rational Field
            sage: sage.modular.hecke.homspace.HeckeModuleHomspace(M, M.base_extend(Qp(13)))
            Traceback (most recent call last):
            ...
            TypeError: X and Y must have the same base ring
            sage: M.Hom(M) == loads(dumps(M.Hom(M)))
            True

        TESTS::

            sage: M = ModularForms(Gamma0(7), 4)
            sage: H = M.Hom(M)
            sage: TestSuite(H).run(skip='_test_elements')
        """
        if not isinstance(X, HeckeModule_generic) or not isinstance(Y, HeckeModule_generic):
            raise TypeError("X and Y must be Hecke modules")
        if X.base_ring() != Y.base_ring():
            raise TypeError("X and Y must have the same base ring")
        if category is None:
            category = X.category()
        HomsetWithBase.__init__(self, X, Y, category=category)

    def __call__(self, A, name='', **kwds):
        r"""
        Create an element of the homspace ``self`` from `A`.

        INPUT:

        - ``A`` -- one of the following:

          - an element of a Hecke algebra

          - a Hecke module morphism

          - a matrix

          - a list of elements of the codomain specifying the images
            of the basis elements of the domain.

        EXAMPLES::

            sage: M = ModularForms(Gamma0(7), 4)
            sage: H = M.Hom(M)
            sage: H(M.hecke_operator(7))
            Hecke module morphism T_7 defined by the matrix
            [ -7   0   0]
            [  0   1 240]
            [  0   0 343]
            Domain:   Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            Codomain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            sage: H(H(M.hecke_operator(7)))
            Hecke module morphism T_7 defined by the matrix
            [ -7   0   0]
            [  0   1 240]
            [  0   0 343]
            Domain:   Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            Codomain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            sage: H(matrix(QQ, 3, srange(9)))
            Hecke module morphism defined by the matrix
            [0 1 2]
            [3 4 5]
            [6 7 8]
            Domain:   Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            Codomain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...

        TESTS:

        Make sure that the element is created correctly when the codomain is
        not the full module (related to :issue:`21497`)::

            sage: M = ModularSymbols(Gamma0(3), weight=22, sign=1)
            sage: S = M.cuspidal_subspace()
            sage: H = S.Hom(S)
            sage: H(S.gens())
            Hecke module morphism defined by the matrix
            [1 0 0 0 0 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
            Domain:   Modular Symbols subspace of dimension 6 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...

            sage: H.zero() in H
            True
            sage: H.one() in H
            True
        """
        try:
            if A.parent() == self:
                A._set_parent(self)
                return A
            A = A.hecke_module_morphism()
            if A.parent() == self:
                A._set_parent(self)
                return A
            else:
                raise TypeError("unable to coerce A to self")
        except AttributeError:
            pass
        side = kwds.get("side", "left")
        if A in self.base_ring():
            dim_dom = self.domain().rank()
            dim_codom = self.codomain().rank()
            if side == "left":
                MS = MatrixSpace(self.base_ring(), dim_dom, dim_codom)
            else:
                MS = MatrixSpace(self.base_ring(), dim_codom, dim_dom)
            if self.domain() == self.codomain():
                A = A * MS.identity_matrix()
            elif A == 0:
                A = MS.zero()
            else:
                raise ValueError('scalars do not coerce to this homspace')
        elif isinstance(A, (list, tuple)):
            A = matrix([self.codomain().coordinate_vector(f) for f in A])
        if side == "right":
            A = A.transpose()
        return HeckeModuleMorphism_matrix(self, A, name, side)

    def _an_element_(self):
        """
        Return an element.

        If the domain is equal to the codomain, this returns the
        action of the Hecke operator of index 2. Otherwise, this returns zero.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: M = ModularSymbols(Gamma0(2), weight=12, sign=1)
            sage: S = M.cuspidal_subspace()
            sage: S.Hom(S).an_element()
            Hecke module morphism defined by the matrix
            [      260 -2108/135]
            [     4860      -284]
            Domain:   Modular Symbols subspace of dimension 2 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 2 of Modular Symbols space ...
        """
        if self.domain() != self.codomain():
            return self.zero()
        else:
            A = self.domain().hecke_operator(2).matrix()
            return HeckeModuleMorphism_matrix(self, A)
