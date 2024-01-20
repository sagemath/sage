# sage.doctest: needs sage.graphs          (because all doctests use SimplicialComplex)
r"""
Homspaces between chain complexes

Note that some significant functionality is lacking. Namely, the homspaces
are not actually modules over the base ring. It will be necessary to
enrich some of the structure of chain complexes for this to be naturally
available. On other hand, there are various overloaded operators. ``__mul__``
acts as composition. One can ``__add__``, and one can ``__mul__`` with a ring
element on the right.

EXAMPLES::

    sage: S = simplicial_complexes.Sphere(2)
    sage: T = simplicial_complexes.Torus()
    sage: C = S.chain_complex(augmented=True, cochain=True)
    sage: D = T.chain_complex(augmented=True, cochain=True)
    sage: G = Hom(C, D); G
    Set of Morphisms
     from Chain complex with at most 4 nonzero terms over Integer Ring
       to Chain complex with at most 4 nonzero terms over Integer Ring
       in Category of chain complexes over Integer Ring

    sage: S = simplicial_complexes.ChessboardComplex(3, 3)
    sage: H = Hom(S,S)
    sage: i = H.identity()
    sage: x = i.associated_chain_complex_morphism(augmented=True); x
    Chain complex morphism:
      From: Chain complex with at most 4 nonzero terms over Integer Ring
      To:   Chain complex with at most 4 nonzero terms over Integer Ring
    sage: x._matrix_dictionary
    {-1: [1],
      0: [1 0 0 0 0 0 0 0 0]
         [0 1 0 0 0 0 0 0 0]
         [0 0 1 0 0 0 0 0 0]
         [0 0 0 1 0 0 0 0 0]
         [0 0 0 0 1 0 0 0 0]
         [0 0 0 0 0 1 0 0 0]
         [0 0 0 0 0 0 1 0 0]
         [0 0 0 0 0 0 0 1 0]
         [0 0 0 0 0 0 0 0 1],
      1: [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
         [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
         [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
         [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
         [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]
         [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
         [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
         [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
         [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
         [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
         [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
         [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
         [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]
         [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
         [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
         [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
         [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
         [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1],
      2: [1 0 0 0 0 0]
         [0 1 0 0 0 0]
         [0 0 1 0 0 0]
         [0 0 0 1 0 0]
         [0 0 0 0 1 0]
         [0 0 0 0 0 1]}

    sage: S = simplicial_complexes.Sphere(2)
    sage: A = Hom(S, S)
    sage: i = A.identity()
    sage: x = i.associated_chain_complex_morphism(); x
    Chain complex morphism:
      From: Chain complex with at most 3 nonzero terms over Integer Ring
      To: Chain complex with at most 3 nonzero terms over Integer Ring
    sage: y = x*4
    sage: z = y*y
    sage: y + z
    Chain complex morphism:
      From: Chain complex with at most 3 nonzero terms over Integer Ring
      To: Chain complex with at most 3 nonzero terms over Integer Ring
    sage: f = x._matrix_dictionary
    sage: C = S.chain_complex()
    sage: G = Hom(C, C)
    sage: w = G(f)
    sage: w == x
    True

"""

# ****************************************************************************
# Copyright (C) 2009 D. Benjamin Antieau <d.ben.antieau@gmail.com>
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
#
# ****************************************************************************

import sage.categories.homset
from sage.homology.chain_complex_morphism import ChainComplexMorphism


def is_ChainComplexHomspace(x):
    """
    Return ``True`` if and only if ``x`` is a morphism of chain complexes.

    EXAMPLES::

        sage: from sage.homology.chain_complex_homspace import is_ChainComplexHomspace
        sage: T = SimplicialComplex([[1,2,3,4],[7,8,9]])
        sage: C = T.chain_complex(augmented=True, cochain=True)
        sage: G = Hom(C, C)
        sage: is_ChainComplexHomspace(G)
        True

    """
    return isinstance(x, ChainComplexHomspace)


class ChainComplexHomspace(sage.categories.homset.Homset):
    """
    Class of homspaces of chain complex morphisms.

    EXAMPLES::

        sage: T = SimplicialComplex([[1,2,3,4],[7,8,9]])
        sage: C = T.chain_complex(augmented=True, cochain=True)
        sage: G = Hom(C, C)
        sage: G
        Set of Morphisms
         from Chain complex with at most 5 nonzero terms over Integer Ring
           to Chain complex with at most 5 nonzero terms over Integer Ring
           in Category of chain complexes over Integer Ring

    """
    def __call__(self, f):
        """
        `f` is a dictionary of matrices in the basis of the chain complex.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(5)
            sage: H = Hom(S, S)
            sage: i = H.identity()
            sage: C = S.chain_complex()
            sage: G = Hom(C, C)
            sage: x = i.associated_chain_complex_morphism()
            sage: f = x._matrix_dictionary
            sage: y = G(f)
            sage: x == y
            True

        """
        return ChainComplexMorphism(f, self.domain(), self.codomain())
