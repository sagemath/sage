
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



    EXAMPLES::

        sage: F = GF(25); M = F^2; twist = F.frobenius_endomorphism(5)
        sage: PHS = F.PseudoHom(twist)

    """
    def __init__(self, X, Y, twist=None):
        self._domain = X
        self._codomain = X
        if Y is not None:
            self._codomain = Y
        if (twist.domain() is not self.domain().coordinate_ring()
          or twist.codomain() is not self.codomain().coordinate_ring()):
            raise TypeError("twisting morphism domain/codomain do not match coordinate rings of the modules")
        elif isinstance(twist, Morphism) or isinstance(twist, RingDerivation):
            self.twist = twist
        else:
            raise TypeError("twist is not a ring morphism or derivation")

    def __call__(self, A, **kwds):
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
                        A = matrix([C.coordinates(a) for a in v], ncols=C.rank()).transpose()
                    else:
                        A = matrix([C.coordinates(a) for a in v], ncols=C.rank())
            except TypeError:
                pass
        if not self.codomain().base_ring().has_coerce_map_from(self.domain().base_ring()) and not A.is_zero():
            raise TypeError("nontrivial morphisms require a coercion map from the base ring of the domain to the base ring of the codomain")
        return free_module_pseudomorphism.FreeModulePseudoMorphism(self.domain(), A, twist=self.twist, codomain = self.codomain())

    def __repr__(self):
        pass
