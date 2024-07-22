"""
Space of Pseudomorphisms of free modules

AUTHORS:

    - Yossef Musleh (2024-02): initial version

"""
# ****************************************************************************
#  Copyright (C) 2024 Yossef Musleh <specialholonomy@gmail.com>
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
import sage.categories.homset
from sage.structure.element import is_Matrix
from sage.matrix.constructor import matrix, identity_matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.categories.morphism import Morphism
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.derivation', 'RingDerivation')

class FreeModulePseudoHomspace(sage.categories.homset.HomsetWithBase):
    r"""
    This class implements the space of Pseudomorphisms with a fixed twist.

    For free modules, the elements of a pseudomorphism correspond to matrices
    which define the mapping on elements of a basis.

    TESTS::

        sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
        sage: PHS = M.PseudoHom(twist)
        sage: h = PHS([[1, 2], [1, 1]])
        sage: e = M((4*F.gen()^2 + F.gen() + 2, 4*F.gen()^2 + 4*F.gen() + 4))
        sage: h(e)
        (z3, 2*z3^2 + 3*z3 + 3)

    """
    def __init__(self, domain, codomain=None, twist=None):
        r"""
        Constructs the space of pseudomorphisms with a given twist.

        INPUT:
            -  ``domain``   - the domain of the pseudomorphism; a free module

            -  ``codomain`` - the codomain of the pseudomorphism; a free
                             module  (default: None)

            -  ``twist`` - a twisting morphism, this is either a morphism or
                           a derivation (default: None)

        EXAMPLES::

            sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(twist); PHS
            Set of Pseudomorphisms from Vector space of dimension 2 over Finite Field in z3 of size 5^3 to Vector space of dimension 2 over Finite Field in z3 of size 5^3
            Twisted by the morphism Frobenius endomorphism z3 |--> z3^5 on Finite Field in z3 of size 5^3
        """
        self._domain = domain
        self._codomain = domain
        if codomain is not None:
            self._codomain = codomain
        super().__init__(self._domain, self._codomain, category=None)
        self.base_homspace = self._domain.Hom(self._codomain)
        self.twist = twist
        self.twist_morphism = None
        self.derivation = None
        if twist is None:
            return
        if (twist.domain() is not self.domain().coordinate_ring()
          or twist.codomain() is not self.codomain().coordinate_ring()):
            raise TypeError("twisting morphism domain/codomain do not match\
                            coordinate rings of the modules")
        elif isinstance(twist, Morphism):
            self.twist_morphism = twist
        elif isinstance(twist, RingDerivation):
            self.derivation = twist
        else:
            raise TypeError("twist is not a ring morphism or derivation")

    def __call__(self, A, **kwds):
        r"""
        Coerce a matrix or free module homomorphism into a pseudomorphism.

        INPUTS:
            - ``A`` - either a matrix defining the morphism or a free module
                      morphism

        EXAMPLES::

            sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(twist)
            sage: h = PHS([[1, 2], [1, 1]]); h
            Free module pseudomorphism defined by the matrix
            [1 2]
            [1 1]
            twisted by the morphism Frobenius endomorphism z3 |--> z3^5 on Finite Field in z3 of size 5^3
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
        """
        return self._element_constructor_(A, **kwds)

    def __repr__(self):
        r"""
        Returns the string representation of the pseudomorphism space.

        EXAMPLES::

            sage: Fq = GF(343); M = Fq^2; frob = Fq.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(frob); PHS
            Set of Pseudomorphisms from Vector space of dimension 2 over Finite Field in z3 of size 7^3 to Vector space of dimension 2 over Finite Field in z3 of size 7^3
            Twisted by the morphism Frobenius endomorphism z3 |--> z3^7 on Finite Field in z3 of size 7^3
        """
        r = "Set of Pseudomorphisms from {} to {} {} {}"
        morph = ""
        if self.twist_morphism is not None:
            morph = "\nTwisted by the morphism {}"
            morph = morph.format(self.twist_morphism.__repr__())
        deriv = ""
        if self.derivation is not None:
            deriv = "\nTwisted by the derivation {}"
            deriv = deriv.format(self.derivation.__repr__())
        return r.format(self.domain(), self.codomain(), morph, deriv)

    def _element_constructor_(self, A, **kwds):
        r"""
        Coerce a matrix or free module homomorphism into a pseudomorphism.

        INPUTS:
            - ``A`` - either a matrix defining the morphism or a free module
                      morphism

        TESTS::

            sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(twist)
            sage: h = PHS._element_constructor_([[1, 2], [1, 1]]); h
            Free module pseudomorphism defined by the matrix
            [1 2]
            [1 1]
            twisted by the morphism Frobenius endomorphism z3 |--> z3^5 on Finite Field in z3 of size 5^3
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 5^3

        ::

            sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(twist)
            sage: morph = M.hom(matrix([[1, 2], [1, 1]]))
            sage: phi = PHS._element_constructor_(morph, side="right"); phi
            Free module pseudomorphism defined as left-multiplication by the matrix
            [1 2]
            [1 1]
            twisted by the morphism Frobenius endomorphism z3 |--> z3^5 on Finite Field in z3 of size 5^3
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
        """
        from . import free_module_pseudomorphism as pseudo
        side = kwds.get("side", "left")
        if not self.codomain().base_ring().has_coerce_map_from(
                self.domain().base_ring()) and not A.is_zero():
            raise TypeError("nontrivial morphisms require a coercion map"
                    "from the base ring of the domain to the base ring of the"
                    "codomain")
        return pseudo.FreeModulePseudoMorphism(self, A, side=side)

    def _matrix_space(self):
        r"""
        Return the full matrix space of the underlying morphisms.

        EXAMPLES::

            sage: Fq = GF(343); M = Fq^2; frob = Fq.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(frob)
            sage: PHS._matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field in z3 of size 7^3
        """
        return self.base_homspace._matrix_space()

    def basis(self, side="left"):
        r"""
        Return a basis for the underlying matrix space.

        EXAMPLES::

            sage: Fq = GF(343); M = Fq^2; frob = Fq.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(frob)
            sage: PHS.basis()
            (Vector space morphism represented by the matrix:
             [1 0]
             [0 0]
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
             Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3,
             Vector space morphism represented by the matrix:
             [0 1]
             [0 0]
             Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
             Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3,
             Vector space morphism represented by the matrix:
             [0 0]
             [1 0]
             Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
             Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3,
             Vector space morphism represented by the matrix:
             [0 0]
             [0 1]
             Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
             Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3)
        """
        return self.base_homspace.basis(side)

    def identity(self):
        r"""
        Return the pseudomorphism corresponding to the identity transformation

        EXAMPLES::

            sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(twist)
            sage: PHS.identity()
            Free module pseudomorphism defined by the matrix
            [1 0]
            [0 1]
            twisted by the morphism Frobenius endomorphism z3 |--> z3^5 on Finite Field in z3 of size 5^3
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
        """
        return self(self.base_homspace.identity())

    def zero(self):
        r"""
        Return the zero pseudomorphism. This corresponds to the zero matrix.

        EXAMPLES::

            sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
            sage: PHS = M.PseudoHom(twist)
            sage: PHS.zero()
            Free module pseudomorphism defined by the matrix
            [0 0]
            [0 0]
            twisted by the morphism Frobenius endomorphism z3 |--> z3^5 on Finite Field in z3 of size 5^3
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
        """
        return self(self.base_homspace.zero())
