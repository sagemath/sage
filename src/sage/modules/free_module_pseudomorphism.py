"""
Pseudomorphisms of free modules

AUTHORS:


TESTS::

    sage: V = ZZ^2; f = V.hom([V.1, -2*V.0])
"""

####################################################################################
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
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

# A matrix morphism is a morphism that is defined by multiplication by a
# matrix.  Elements of domain must either have a method "vector()" that
# returns a vector that the defining matrix can hit from the left, or
# be coercible into vector space of appropriate dimension.
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
    def __init__(self, domain, base_morphism, twist=None, codomain=None, side="left"):
        """
        INPUT:
            -  ``domain``   - the domain of the pseudomorphism; a free module 

            -  ``base_morphism`` - either a morphism or a matrix defining a morphism

            -  ``twist`` - a twisting morphism, this is either a morphism or a derivation (default: None)

            -  ``codomain`` - the codomain of the pseudomorphism; a free module  (default: None)

            - side -- side of the vectors acted on by the matrix  (default: ``"left"``)

        EXAMPLES::

            sage: V = ZZ^3; W = span([[1,2,3],[-1,2,8]], ZZ)
            sage: phi = V.hom(matrix(ZZ,3,[1..9]))
            sage: type(phi)
            <class 'sage.modules.free_module_morphism.FreeModuleMorphism'>
        """
        from sage.structure.element import is_Matrix
        if is_Matrix(base_morphism):
            self.base_morphism = domain.hom(base_morphism, codomain)
        elif isinstance(base_morphism, Morphism):
            self.base_morphism = base_morphism
        else:
            self.base_morphism = domain.hom(matrix(domain.coordinate_ring(), base_morphism), codomain)
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
            if self.twist_morphism is None:
                x_twistmorphism = x
            else:
                x_twistmorphism = self.domain()(list(map(self.twist_morphism, x)))
            C = self.codomain()
            if self.side == "left":
                v = x_twistmorphism * self.matrix()
            else:
                v = self.matrix() * x_twistmorphism
            if self.derivation is not None:
                v += self.domain()(list(map(self.derivation, x)))
            if not C.is_ambient():
                v = C.linear_combination_of_basis(v)
            return C._element_constructor_(v)

    def __repr__(self):
        r = "Free module pseudomorphism defined {}by the matrix\n{!r}{}{}\nDomain: {}\nCodomain: {}"
        act = ""
        if self.side == "right":
            act = "as left-multiplication "
        morph = ""
        if self.twist_morphism is not None:
            morph = "\nTwisted by the morphism {}"
            morph = morph.format(self.twist_morphism.__repr__())
        deriv = ""
        if self.derivation is not None:
            deriv = "\nTwisted by the derivation {}"
            deriv = deriv.format(self.derivation.__repr__())
        return r.format(act, self.matrix(), morph, deriv, self.domain(), self.codomain())

    def domain(self):
        return self.base_morphism.domain()

    def codomain(self):
        return self.base_morphism.codomain()

    def matrix(self):
        return self.base_morphism.matrix()

    def base_morphism(self):
        return self.base_morphism
