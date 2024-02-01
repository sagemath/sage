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

lazy_import('sage.rings.derivation', 'RingDerivation')

class FreeModulePseudoMorphism(Morphism):
    r"""
    Let `M, M'` be free modules over a ring `R`, `\theta: R \to R` a ring
    homomorphism, and `\theta: R \to R` a derivation i.e. an additive
    map such that

    `\delta(xy) = x\delta(y) + \delta(x)y`.

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
    """
    def __init__(self, domain, base_morphism, twist=None, codomain=None, side="left"):
        """
        Constructs a pseudomorphism of free modules.

        INPUT:
            -  ``domain``   - the domain of the pseudomorphism; a free module 

            -  ``base_morphism`` - either a morphism or a matrix defining a
                                   morphism

            -  ``twist`` - a twisting morphism, this is either a morphism or 
                           a derivation (default: None)

            -  ``codomain`` - the codomain of the pseudomorphism; a free
                              module  (default: None)

            - side -- side of the vectors acted on by the matrix 
                      (default: ``"left"``)

        EXAMPLES::

            sage: F = GF(25); V = F^3; twist = F.frobenius_endomorphism(5)
            sage: phi = V.pseudohom(matrix(F,3,[1..9]), twist)
            sage: type(phi)
            <class 'sage.modules.free_module_pseudomorphism.FreeModulePseudoMorphism'>
        """
        from sage.structure.element import is_Matrix
        if is_Matrix(base_morphism):
            self.base_morphism = domain.hom(base_morphism, codomain)
        elif isinstance(base_morphism, Morphism):
            self.base_morphism = base_morphism
        else:
            self.base_morphism = domain.hom(matrix(domain.coordinate_ring(), \
                                            base_morphism), codomain)
        self.derivation = None
        self.twist_morphism = None
        if isinstance(twist, Morphism):
            self.twist_morphism = twist
        elif isinstance(twist, RingDerivation):
            self.twist_morphism = twist.parent().twisting_morphism()
            if twist:
                self.derivation = twist
            else:
                self.derivation = None
        self.side = side

    def __call__(self, x):
        r"""
        Return the result of applying a pseudomorphism to an element of the
        free module.

        TESTS::

        sage: Fq = GF(343); M = Fq^3; frob = Fq.frobenius_endomorphism()
        sage: ph = M.pseudohom([[1, 2, 3], [0, 1, 1], [2, 1, 1]], frob, side="right")
        sage: e = M((3*Fq.gen()^2 + 5*Fq.gen() + 2, 6*Fq.gen()^2 + 2*Fq.gen() + 2, Fq.gen() + 4))
        sage: ph(e)
        (3*z3 + 3, 2*z3^2, 5*z3^2 + 4)
        """
        if self.twist_morphism is None and self.derivation is None:
            return self.base_morphism(x)
        else:
            try:
                if parent(x) is not self.domain():
                    x = self.domain()(x)
            except TypeError:
                raise TypeError("%s must be coercible into %s" % (x,self.domain()))
            if self.domain().is_ambient():
                x = x.element()
            else:
                x = self.domain().coordinate_vector(x)
            C = self.codomain()
            if self.side == "left":
                v = x * self.matrix()
            else:
                v = self.matrix() * x
            if self.twist_morphism is not None:
                for i in range(len(v)):
                    v[i] *= self.twist_morphism(x[i])
            if self.derivation is not None:
                v += self.domain()(list(map(self.derivation, x)))
            if not C.is_ambient():
                v = C.linear_combination_of_basis(v)
            return C._element_constructor_(v)

    def __repr__(self):
        r"""
        Return the string representation of a pseudomorphism.

        TESTS::
        """
        r = "Free module pseudomorphism defined {}by the \
                matrix\n{!r}{}{}\nDomain: {}\nCodomain: {}"
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
        return r.format(act, self.matrix(), morph, deriv, \
                        self.domain(), self.codomain())

    def domain(self):
        r"""
        Return the domain of the pseudomorphism.
        """
        return self.base_morphism.domain()

    def codomain(self):
        r"""
        Return the codomain of the pseudomorphism.
        """
        return self.base_morphism.codomain()

    def matrix(self):
        r"""
        Return the underlying matrix of a pseudomorphism.
        """
        return self.base_morphism.matrix()

    def base_morphism(self):
        r"""
        Return the underlying morphism of a pseudomorphism. This is an element
        of the Hom space of the free module.
        """
        return self.base_morphism

    def twisting_morphism(self):
        r"""
        Return the twisting homomorphism of the pseudomorphism.
        """
        return self.twist_morphism

    def derivation(self):
        r"""
        Return the twisting derivation of the pseudomorphism.
        """
        return self.derivation
