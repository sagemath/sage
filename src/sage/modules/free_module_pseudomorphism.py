"""
Pseudomorphisms of free modules

AUTHORS:

    - Yossef Musleh (2024-02): initial version

"""
####################################################################################
#       Copyright (C) 2024 Yossef Musleh <jbobicus@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
####################################################################################

import sage.modules.free_module as free_module
from sage.categories.morphism import Morphism
from sage.modules import free_module_homspace, matrix_morphism
from sage.structure.richcmp import rich_to_bool, richcmp
from sage.structure.sequence import Sequence
from sage.structure.all import parent
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module_morphism import FreeModuleMorphism
from sage.modules.free_module_homspace import FreeModuleHomspace, is_FreeModuleHomspace
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace

lazy_import('sage.rings.derivation', 'RingDerivation')

class FreeModulePseudoMorphism(Morphism):
    r"""
    Let `M, M'` be free modules over a ring `R`, `\theta: R \to R` a ring
    homomorphism, and `\delta: R \to R` a `\theta`-derivation, which is a map
    such that:

    `\delta(xy) = \theta(x)\delta(y) + \delta(x)y`.

    Then a pseudomorphism `f : M to M` is a map such that

    `f(x + y) = f(x) + f(y)`
    `f(\lambda x) = `\theta(\lambda)f(x) + \delta(\lambda)x`

    The pair `(\theta, \delta)` may be referred to as the *twist* of
    the morphism.

    TESTS::

        sage: V = ZZ^2
        sage: f = V.pseudohom([V.1, -2*V.0]); f
        Free module pseudomorphism defined by the matrix
        [ 0  1]
        [-2  0]
        Domain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
        Codomain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
        sage: f(V((1, 2)))
        (-4, 1)

    ::

        sage: P.<x> = ZZ[]; deriv = P.derivation()
        sage: M = P^2
        sage: f = M.pseudohom([[1, 2*x], [x, 1]], deriv, side="right"); f
        Free module pseudomorphism defined as left-multiplication by the matrix
        [  1 2*x]
        [  x   1]
        twisted by the derivation d/dx
        Domain: Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring
        Codomain: Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring
        sage: e = M((2*x^2 + 3*x + 1, x^3 + 7*x + 4))
        sage: f(e)
        (2*x^4 + 16*x^2 + 15*x + 4, 3*x^3 + 6*x^2 + 8*x + 11)
        sage: f = M.pseudohom([[1, 2], [1, 1]], deriv)
        sage: f(e)
        (x^3 + 2*x^2 + 14*x + 8, x^3 + 7*x^2 + 13*x + 13)

    ::

        sage: Fq = GF(343); M = Fq^3; N = Fq^2; frob = Fq.frobenius_endomorphism(); z = Fq.gen()
        sage: phi = M.pseudohom([[2, 3, 1], [1, 4, 6]], frob, N, side="right"); phi
        Free module pseudomorphism defined as left-multiplication by the matrix
        [2 3 1]
        [1 4 6]
        twisted by the morphism Frobenius endomorphism z3 |--> z3^7 on Finite Field in z3 of size 7^3
        Domain: Vector space of dimension 3 over Finite Field in z3 of size 7^3
        Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
        sage: elem = (4*z^2 + 4*z + 3, 2, z + 5)
        sage: phi(elem)
        (2*z3 + 1, 6*z3^2 + 4*z3 + 5)
    """
    def __init__(self, pseudohomspace, base_morphism, side="left"):
        """
        Constructs a pseudomorphism of free modules.

        INPUT:
            -  ``pseudohomspace`` - the parent space of pseudomorphisms,
                                    containing

            -  ``base_morphism`` - either a morphism or a matrix defining a
                                   morphism

            - side -- side of the vectors acted on by the matrix
                      (default: ``"left"``)

        EXAMPLES::

            sage: F = GF(25); M = F^3; twist = F.frobenius_endomorphism(5)
            sage: phi = M.pseudohom(matrix(F,3,[1..9]), twist)
            sage: type(phi)
            <class 'sage.modules.free_module_pseudomorphism.FreeModulePseudoMorphism'>

        ::

            sage: F = GF(125); M = F^2; twist = F.frobenius_endomorphism()
            sage: morph = M.hom(matrix([[1,2],[0,1]]))
            sage: phi = M.pseudohom(morph, twist, side="right"); phi
            Free module pseudomorphism defined as left-multiplication by the matrix
            [1 2]
            [0 1]
            twisted by the morphism Frobenius endomorphism z3 |--> z3^5 on Finite Field in z3 of size 5^3
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 5^3
        """
        Morphism.__init__(self, pseudohomspace)
        dom = pseudohomspace.domain()
        codom = pseudohomspace.codomain()
        rows = dom.dimension()
        cols = codom.dimension()
        if side == "right":
            rows = codom.dimension()
            cols = dom.dimension()
        matrix_space = MatrixSpace(dom.coordinate_ring(), rows, cols)
        if isinstance(base_morphism, FreeModuleMorphism):
            self._base_matrix = matrix_space(base_morphism.matrix())
        else:
            self._base_matrix = matrix_space(base_morphism)
        self.derivation = pseudohomspace.derivation
        self.twist_morphism = pseudohomspace.twist_morphism
        self.side = side

    def _call_(self, x):
        r"""
        Return the result of applying a pseudomorphism to an element of the
        free module.

        TESTS::

            sage: Fq = GF(343); M = Fq^3; frob = Fq.frobenius_endomorphism()
            sage: ph = M.pseudohom([[1, Fq.gen(), 3], [0, 1, 1], [2, 1, 1]], frob, side="right")
            sage: e = M((3*Fq.gen()^2 + 5*Fq.gen() + 2, 6*Fq.gen()^2 + 2*Fq.gen() + 2, Fq.gen() + 4))
            sage: ph(e)
            (z3^2 + 6*z3 + 2, z3^2 + 2*z3 + 1, 2*z3^2 + 4*z3)
        """
        if self.domain().is_ambient():
            x = x.element()
        else:
            x = self.domain().coordinate_vector(x)
        C = self.codomain()
        if self.twist_morphism is None:
            x_twist = x
        else:
            x_twist = self.domain()(list(map(self.twist_morphism, x)))
        if self.side == "left":
            v = x_twist * self._base_matrix
        else:
            v = self._base_matrix * x_twist
        if self.derivation is not None:
            v += self.domain()(list(map(self.derivation, x)))
        if not C.is_ambient():
            v = C.linear_combination_of_basis(v)
        return C._element_constructor_(v)

    def __repr__(self):
        r"""
        Return the string representation of a pseudomorphism.

        TESTS::

            sage: Fq = GF(343); M = Fq^3; frob = Fq.frobenius_endomorphism()
            sage: ph = M.pseudohom([[1,1,1],[2,2,2],[3,3,3]], frob); ph
            Free module pseudomorphism defined by the matrix
            [1 1 1]
            [2 2 2]
            [3 3 3]
            twisted by the morphism Frobenius endomorphism z3 |--> z3^7 on Finite Field in z3 of size 7^3
            Domain: Vector space of dimension 3 over Finite Field in z3 of size 7^3
            Codomain: Vector space of dimension 3 over Finite Field in z3 of size 7^3
        """
        r = "Free module pseudomorphism defined {}by the "\
        "matrix\n{!r}{}{}\nDomain: {}\nCodomain: {}"
        act = ""
        if self.side == "right":
            act = "as left-multiplication "
        morph = ""
        if self.twist_morphism is not None:
            morph = "\ntwisted by the morphism {}"
            morph = morph.format(self.twist_morphism.__repr__())
        deriv = ""
        if self.derivation is not None:
            deriv = "\ntwisted by the derivation {}"
            deriv = deriv.format(self.derivation.__repr__())
        return r.format(act, self.matrix(), morph, deriv,
                        self.domain(), self.codomain())

    def matrix(self):
        r"""
        Return the underlying matrix of a pseudomorphism.

        If a pseudomorphism `f` on free module `M` has matrix m acting on
        the left on elements `v \in M`, with twisting morphism `\theta`.
        Then we have

        `f(v) = m*\theta(v)`

        where `\theta` acts of the coefficients of `v` in terms of the basis
        for `m`.

        EXAMPLES::

            sage: Fq = GF(343); M = Fq^3; frob = Fq.frobenius_endomorphism()
            sage: ph = M.pseudohom([[1, 2, 3], [0, 1, 1], [2, 1, 1]], frob, side="right")
            sage: e = M((3*Fq.gen()^2 + 5*Fq.gen() + 2, 6*Fq.gen()^2 + 2*Fq.gen() + 2, Fq.gen() + 4))
            sage: ph.matrix()
            [1 2 3]
            [0 1 1]
            [2 1 1]
            sage: ph(e) == ph.matrix()*vector([frob(c) for c in e])
            True
        """
        return self._base_matrix

    def twisting_derivation(self):
        r"""
        Return the twisting derivation of the pseudomorphism.

        EXAMPLES::

            sage: P.<x> = ZZ[]; deriv = P.derivation(); M = P^2
            sage: f = M.pseudohom([[1, 2*x], [x, 1]], deriv, side="right")
            sage: f.twisting_derivation()
            d/dx
        """
        return self.derivation

    def twisting_morphism(self):
        r"""
        Return the twisting homomorphism of the pseudomorphism.

        EXAMPLES::

            sage: Fq = GF(343); M = Fq^3; frob = Fq.frobenius_endomorphism()
            sage: ph = M.pseudohom([[1, 2, 3], [0, 1, 1], [2, 1, 1]], frob, side="right")
            sage: ph.twisting_morphism()
            Frobenius endomorphism z3 |--> z3^7 on Finite Field in z3 of size 7^3
        """
        return self.twist_morphism
