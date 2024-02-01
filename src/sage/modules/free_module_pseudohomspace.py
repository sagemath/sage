"""
Space of Pseudomorphisms of free modules

AUTHORS:

    - Yossef Musleh (2024-02): initial version

"""
# ****************************************************************************
#  Copyright (C) 2024 Yossef Musleh <jbobicus@gmail.com>
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
from sage.misc.cachefunc import cached_method
from sage.categories.morphism import Morphism
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.derivation', 'RingDerivation')

class FreeModulePseudoHomspace(sage.categories.homset.HomsetWithBase):
    r"""
    This class implements the space of Pseudomorphisms with a fixed twist.

    For free modules, the elements of a pseudomorphism correspond to matrices
    which define the mapping on elements of a basis.

    EXAMPLES::

        sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
        sage: PHS = M.PseudoHom(twist)
        sage: h = PHS([[1, 2], [1, 1]])
        sage: e = M((4*F.gen()^2 + F.gen() + 2, 4*F.gen()^2 + 4*F.gen() + 4))
        sage: h(e)
        (z3, 2*z3^2 + 3*z3 + 3)

    """
    def __init__(self, X, Y, twist=None):
        self._domain = X
        self._codomain = X
        if Y is not None:
            self._codomain = Y
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
        from . import free_module_pseudomorphism
        side = kwds.get("side", "left")
        if not is_Matrix(A):
            C = self.codomain()
            try:
                if callable(A):
                    v = [C(A(g)) for g in self.domain().gens()]
                    A = matrix([C.coordinates(a) for a in v], ncols=C.rank())
                    if side == "right":
                        A = A.transpose()
                else:
                    v = [C(a) for a in A]
                    if side == "right":
                        A = matrix([C.coordinates(a) for a in v], \
                                    ncols=C.rank()).transpose()
                    else:
                        A = matrix([C.coordinates(a) for a in v], \
                                    ncols=C.rank())
            except TypeError:
                pass
        if not self.codomain().base_ring().has_coerce_map_from(\
                self.domain().base_ring()) and not A.is_zero():
            raise TypeError("nontrivial morphisms require a coercion map \
                    from the base ring of the domain to the base ring of the \
                    codomain")
        return free_module_pseudomorphism.FreeModulePseudoMorphism(\
                                self.domain(), A, twist=self.twist, \
                                codomain = self.codomain())

    def __repr__(self):
        r"""
        Returns the string representation of the pseudomorphism space.

        EXAMPLE::

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
