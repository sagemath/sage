"""
Elements of Hecke modules

AUTHORS:

- William Stein
"""

#*****************************************************************************
#       Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#*****************************************************************************

import sage.modules.module_element

def is_HeckeModuleElement(x):
    """
    Return True if x is a Hecke module element, i.e., of type HeckeModuleElement.

    EXAMPLES::

        sage: sage.modular.hecke.all.is_HeckeModuleElement(0)
        False
        sage: sage.modular.hecke.all.is_HeckeModuleElement(BrandtModule(37)([1,2,3]))
        True
    """
    return isinstance(x, HeckeModuleElement)

class HeckeModuleElement(sage.modules.module_element.ModuleElement):
    """
    Element of a Hecke module.
    """
    def __init__(self, parent, x=None):
        """
        INPUT:

        - ``parent`` -- a Hecke module

        - ``x`` -- element of the free module associated to parent

        EXAMPLES::

            sage: v = sage.modular.hecke.all.HeckeModuleElement(BrandtModule(37), vector(QQ,[1,2,3])); v
            (1, 2, 3)
            sage: type(v)
            <class 'sage.modular.hecke.element.HeckeModuleElement'>

        TESTS::

            sage: v = ModularSymbols(37).0
            sage: loads(dumps(v))
            (1,0)
            sage: loads(dumps(v)) == v
            True
        """
        sage.modules.module_element.ModuleElement.__init__(self, parent)
        if x is not None:
            self.__element = x

    def _repr_(self):
        """
        Return string representation of this Hecke module element.
        The default representation is just the representation of the
        underlying vector.

        EXAMPLES::

            sage: BrandtModule(37)([0,1,-1])._repr_()
            '(0, 1, -1)'
        """
        return self.element()._repr_()

    def _compute_element(self):
        """
        Use internally to compute vector underlying this element.

        EXAMPLES::

            sage: f = EllipticCurve('11a').modular_form()
            sage: hasattr(f, '_HeckeModuleElement__element')
            False
            sage: f._compute_element()
            (1, 0)
            sage: f.element()
            (1, 0)
            sage: hasattr(f, '_HeckeModuleElement__element')
            True
        """
        # You have to define this in the derived class if you ever set
        # x=None in __init__ for your element class.
        # The main reason for this is it allows for lazy constructors who
        # compute the representation of an element (e.g., a q-expansion) in
        # terms of the basis only when needed.
        raise NotImplementedError("_compute_element *must* be defined in the derived class if element is set to None in constructor")

    def element(self):
        """
        Return underlying vector space element that defines this Hecke module element.

        EXAMPLES::

            sage: z = BrandtModule(37)([0,1,-1]).element(); z
            (0, 1, -1)
            sage: type(z)
            <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        """
        try:
            return self.__element
        except AttributeError:
            self.__element = self._compute_element()
        return self.__element

    def _vector_(self, R=None):
        """
        This makes it so vector(self) and vector(self, R) both work.

        EXAMPLES::

            sage: v = BrandtModule(37)([0,1,-1]); v
            (0, 1, -1)
            sage: type(v._vector_())
            <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
            sage: type(vector(v))
            <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
            sage: type(vector(v, GF(2)))
            <type 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>
        """
        if R is None: return self.__element
        return self.__element.change_ring(R)

    def ambient_module(self):
        """
        Return the ambient Hecke module that contains this element.

        EXAMPLES::

            sage: BrandtModule(37)([0,1,-1]).ambient_module()
            Brandt module of dimension 3 of level 37 of weight 2 over Rational Field
        """
        return self.parent().ambient_module()

    def _lmul_(self, x):
        """
        EXAMPLES::

            sage: BrandtModule(37)([0,1,-1])._lmul_(3)
            (0, 3, -3)
        """
        return self.parent()(self.element()*x)

    def _rmul_(self, x):
        """
        EXAMPLES::

            sage: BrandtModule(37)([0,1,-1])._rmul_(3)
            (0, 3, -3)
        """
        return self.parent()(x * self.element())

    def _neg_(self):
        """
        EXAMPLES::

            sage: BrandtModule(37)([0,1,-1])._neg_()
            (0, -1, 1)
        """
        return self.parent()(-self.element())

    def _pos_(self):
        """
        EXAMPLES::

            sage: BrandtModule(37)([0,1,-1])._pos_()
            (0, 1, -1)
        """
        return self

    def _sub_(self, right):
        """
        EXAMPLES::

            sage: BrandtModule(37)([0,1,-1])._sub_(BrandtModule(37)([0,1,-5]))
            (0, 0, 4)
        """
        return self.parent()(self.element() - right.element())

    def is_cuspidal(self):
        r"""
        Return True if this element is cuspidal.

        EXAMPLES::

            sage: M = ModularForms(2, 22); M.0.is_cuspidal()
            True
            sage: (M.0 + M.4).is_cuspidal()
            False
            sage: EllipticCurve('37a1').newform().is_cuspidal()
            True

        It works for modular symbols too::

            sage: M = ModularSymbols(19,2)
            sage: M.0.is_cuspidal()
            False
            sage: M.1.is_cuspidal()
            True
        """
        return (self in self.parent().ambient().cuspidal_submodule())

    def is_eisenstein(self):
        r"""
        Return True if this element is Eisenstein. This makes sense for both
        modular forms and modular symbols.

        EXAMPLES::

            sage: CuspForms(2,8).0.is_eisenstein()
            False
            sage: M = ModularForms(2,8);(M.0  + M.1).is_eisenstein()
            False
            sage: M.1.is_eisenstein()
            True
            sage: ModularSymbols(19,4).0.is_eisenstein()
            False
            sage: EllipticCurve('37a1').newform().is_eisenstein()
            False
        """
        return (self in self.parent().ambient().eisenstein_submodule())

    def is_new(self, p=None):
        r"""
        Return True if this element is p-new. If p is None, return True if the
        element is new.

        EXAMPLE::

            sage: CuspForms(22, 2).0.is_new(2)
            False
            sage: CuspForms(22, 2).0.is_new(11)
            True
            sage: CuspForms(22, 2).0.is_new()
            False
        """
        return (self in self.parent().new_submodule(p))

    def is_old(self, p=None):
        r"""
        Return True if this element is p-old. If p is None, return True if the
        element is old.

        EXAMPLE::

            sage: CuspForms(22, 2).0.is_old(11)
            False
            sage: CuspForms(22, 2).0.is_old(2)
            True
            sage: CuspForms(22, 2).0.is_old()
            True
            sage: EisensteinForms(144, 2).1.is_old()  # long time (3s on sage.math, 2011)
            False
            sage: EisensteinForms(144, 2).1.is_old(2) # not implemented
            False
        """
        return (self in self.parent().old_submodule(p))

