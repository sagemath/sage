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

lazy_import('sage.rings.derivation', 'RingDerivation')

class FreeModulePseudoMorphism():
    def __init__(self, morphism, twist=None, side="left"):
        """
        INPUT:

            -  ``parent`` - a homspace in a (sub) category of free modules

            -  ``A`` - matrix

            - side -- side of the vectors acted on by the matrix  (default: ``"left"``)

        EXAMPLES::

            sage: V = ZZ^3; W = span([[1,2,3],[-1,2,8]], ZZ)
            sage: phi = V.hom(matrix(ZZ,3,[1..9]))
            sage: type(phi)
            <class 'sage.modules.free_module_morphism.FreeModuleMorphism'>
        """
        self.derivation = None
        if isinstance(twist, Morphism):
            self.twist_morphism = twist
        elif isinstance(twist, RingDerivation):
            self.twist_morphism = twist.parent().twisting_morphism()
            if twist:
                self.derivation = twist
            else:
                self.derivation = None
        self.side = side
        self.base_morphism = morphism
		
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

    def domain(self):
        return self.base_morphism.domain()

    def codomain(self):
        return self.base_morphism.codomain()

    def matrix(self):
        return self.base_morphism.matrix()

    def base_morphism(self):
        return self.base_morphism
