r"""
Hom spaces between Hecke modules
"""
from __future__ import absolute_import

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.categories.homset
from . import morphism
from . import module

def is_HeckeModuleHomspace(x):
    r"""
    Return True if x is a space of homomorphisms in the category of Hecke modules.

    EXAMPLES::

        sage: M = ModularForms(Gamma0(7), 4)
        sage: sage.modular.hecke.homspace.is_HeckeModuleHomspace(Hom(M, M))
        True
        sage: sage.modular.hecke.homspace.is_HeckeModuleHomspace(Hom(M, QQ))
        False
    """
    return isinstance(x, HeckeModuleHomspace)

class HeckeModuleHomspace(sage.categories.homset.HomsetWithBase):
    r"""
    A space of homomorphisms between two objects in the category of Hecke
    modules over a given base ring.
    """
    def __init__(self, X, Y, category = None):
        r"""
        Create the space of homomorphisms between X and Y, which must have the
        same base ring.

        EXAMPLE::

            sage: M = ModularForms(Gamma0(7), 4)
            sage: M.Hom(M)
            Set of Morphisms from ... to ... in Category of Hecke modules over Rational Field
            sage: sage.modular.hecke.homspace.HeckeModuleHomspace(M, M.base_extend(Qp(13)))
            Traceback (most recent call last):
            ...
            TypeError: X and Y must have the same base ring
            sage: M.Hom(M) == loads(dumps(M.Hom(M)))
            True
        """
        if not module.is_HeckeModule(X) or not module.is_HeckeModule(Y):
            raise TypeError("X and Y must be Hecke modules")
        if X.base_ring() != Y.base_ring():
            raise TypeError("X and Y must have the same base ring")
        if category is None:
            category = X.category()
        sage.categories.homset.HomsetWithBase.__init__(self, X, Y, category = category)

    def __call__(self, A, name=''):
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
            Domain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            Codomain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            sage: H(H(M.hecke_operator(7)))
            Hecke module morphism T_7 defined by the matrix
            [ -7   0   0]
            [  0   1 240]
            [  0   0 343]
            Domain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            Codomain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            sage: H(matrix(QQ, 3, srange(9)))
            Hecke module morphism defined by the matrix
            [0 1 2]
            [3 4 5]
            [6 7 8]
            Domain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
            Codomain: Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) ...
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
        if isinstance(A, (list, tuple)):
            from sage.matrix.constructor import matrix
            A = matrix([f.element() for f in A])
        return morphism.HeckeModuleMorphism_matrix(self, A, name)
