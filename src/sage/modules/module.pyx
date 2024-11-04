"""
Abstract base class for modules

AUTHORS:

- William Stein: initial version

- Julian Rueth (2014-05-10): category parameter for Module, doc cleanup

EXAMPLES:

A minimal example of a module::

    sage: from sage.structure.richcmp import richcmp
    sage: class MyElement(sage.structure.element.ModuleElement):
    ....:     def __init__(self, parent, x):
    ....:         self.x = x
    ....:         sage.structure.element.ModuleElement.__init__(self, parent=parent)
    ....:     def _lmul_(self, c):
    ....:         return self.parent()(c*self.x)
    ....:     def _add_(self, other):
    ....:         return self.parent()(self.x + other.x)
    ....:     def _richcmp_(self, other, op):
    ....:         return richcmp(self.x, other.x, op)
    ....:     def __hash__(self):
    ....:         return hash(self.x)
    ....:     def _repr_(self):
    ....:         return repr(self.x)

    sage: from sage.modules.module import Module
    sage: class MyModule(Module):
    ....:     Element = MyElement
    ....:     def _element_constructor_(self, x):
    ....:         if isinstance(x, MyElement): x = x.x
    ....:         return self.element_class(self, self.base_ring()(x))
    ....:     def __eq__(self, other):
    ....:         if not isinstance(other, MyModule): return False
    ....:         return self.base_ring() == other.base_ring()
    ....:     def __hash__(self):
    ....:         return hash(self.base_ring())

    sage: M = MyModule(QQ)
    sage: M(1)
    1

    sage: import __main__
    sage: __main__.MyModule = MyModule
    sage: __main__.MyElement = MyElement
    sage: TestSuite(M).run()
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
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

cdef class Module(Parent):
    """
    Generic module class.

    INPUT:

    - ``base`` -- a ring; the base ring of the module

    - ``category`` -- a category (default: ``None``), the category for this
      module. If ``None``, then this is set to the category of modules/vector
      spaces over ``base``.

    - ``names`` -- names of generators

    EXAMPLES::

        sage: from sage.modules.module import Module
        sage: M = Module(ZZ)
        sage: M.base_ring()
        Integer Ring
        sage: M.category()
        Category of modules over Integer Ring

    Normally the category is set to the category of modules over ``base``. If
    ``base`` is a field, then the category is the category of vector spaces
    over ``base``::

        sage: M_QQ = Module(QQ)
        sage: M_QQ.category()
        Category of vector spaces over Rational Field

    The ``category`` parameter can be used to set a more specific category::

        sage: N = Module(ZZ, category=FiniteDimensionalModulesWithBasis(ZZ))
        sage: N.category()
        Category of finite dimensional modules with basis over Integer Ring

     TESTS:

     We check that :issue:`8119` has been resolved::

        sage: # needs sage.modules
        sage: M = ZZ^3
        sage: h = M.__hash__()
        sage: M.rename('toto')
        sage: h == M.__hash__()
        True
    """
    def __init__(self, base, category=None, names=None):
        """
        Initialization.

        TESTS::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ)
            sage: type(M)
            <class 'sage.modules.module.Module'>
        """
        from sage.categories.modules import Modules
        if category is None:
            category = Modules(base)
        Parent.__init__(self, base=base, category=category, names=names)

    cpdef _coerce_map_from_(self, M):
        """
        Return a coercion map from `M` to ``self``, or None.

        The implementation of this method in module classes should
        observe the following guidelines:

        1. We want to relate two different *ambient* modules if and
           only if they have the same rank (which is the same as
           degree) and if there is a coercion of the base rings.

        2. Two modules embedded in other modules that have a coercion
           may inherit a coercion if they are in fact sub-modules of
           one another.

        3. Since different embeddings of one abstract module are
           related by non-identical coerce maps (in 2.), we must not
           have coercion from a sub-module to the corresponding
           abstract module, for otherwise non-commuting coercion
           triangles emerge.

        4. Quotient modules must not coerce unless their modulus `W`
           is the same. There must not be forgetful maps.

        5. Coerce maps for quotient modules are already registered on
           construction.

        TESTS:

        Make sure :issue:`3638` is fixed::

            sage: vector(ZZ,[1,2,11]) == vector(Zmod(8),[1,2,3])                        # needs sage.modules
            True

        AUTHORS:

        - Simon King (2010-12)

        - Peter Bruin (June 2014)
        """
        try:
            if (isinstance(M, Module)
                    and self.base_ring().has_coerce_map_from(M.base_ring())
                    and M.change_ring(self.base_ring()).is_submodule(self)):
                return M.hom([self._element_constructor_(x) for x in M.gens()], self)
        except (TypeError, NotImplementedError, AttributeError, ArithmeticError):
            pass
        return None

    def change_ring(self, R):
        """
        Return the base change of ``self`` to `R`.

        EXAMPLES::

            sage: from sage.modular.modform.space import ModularFormsSpace              # needs sage.modular
            sage: ModularFormsSpace(Gamma0(11), 2,                                      # needs sage.modular sage.rings.finite_rings
            ....:                   DirichletGroup(1)[0], QQ).change_ring(GF(7))
            Traceback (most recent call last):
            ...
            NotImplementedError: the method change_ring() has not yet been implemented
        """
        if R is self.base_ring():
            return self
        raise NotImplementedError('the method change_ring() has not yet been implemented')

    def base_extend(self, R):
        r"""
        Return the base extension of ``self`` to `R`.

        This is the same as ``self.change_ring(R)`` except that a
        :exc:`TypeError` is raised if there is no canonical coerce map
        from the base ring of ``self`` to `R`.

        INPUT:

        - ``R`` -- ring

        EXAMPLES::

            sage: V = ZZ^7                                                              # needs sage.modules
            sage: V.base_extend(QQ)                                                     # needs sage.modules
            Vector space of dimension 7 over Rational Field

        TESTS::

            sage: N = ModularForms(6, 4)                                                # needs sage.modular
            sage: N.base_extend(CyclotomicField(7))                                     # needs sage.modular sage.rings.number_field
            Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(6)
             of weight 4 over Cyclotomic Field of order 7 and degree 6

            sage: m = ModularForms(DirichletGroup(13).0^2,2); m                         # needs sage.modular sage.rings.number_field
            Modular Forms space of dimension 3, character [zeta6] and weight 2
             over Cyclotomic Field of order 6 and degree 2
            sage: m.base_extend(CyclotomicField(12))                                    # needs sage.modular sage.rings.number_field
            Modular Forms space of dimension 3, character [zeta6] and weight 2
             over Cyclotomic Field of order 12 and degree 4

            sage: # needs sage.modular sage.rings.number_field
            sage: chi = DirichletGroup(109, CyclotomicField(3)).0
            sage: S3 = CuspForms(chi, 2)
            sage: S9 = S3.base_extend(CyclotomicField(9)); S9
            Cuspidal subspace of dimension 8 of
             Modular Forms space of dimension 10, character [zeta3 + 1] and weight 2
              over Cyclotomic Field of order 9 and degree 6
            sage: S9.has_coerce_map_from(S3)    # not implemented
            True
            sage: S9.base_extend(CyclotomicField(3))
            Traceback (most recent call last):
            ...
            TypeError: Base extension of self (over 'Cyclotomic Field of order 9 and degree 6')
            to ring 'Cyclotomic Field of order 3 and degree 2' not defined.
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        raise TypeError("Base extension of self (over '%s') to ring '%s' not defined."
                        % (self.base_ring(), R))

    def endomorphism_ring(self):
        """
        Return the endomorphism ring of this module in its category.

        EXAMPLES::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ)
            sage: M.endomorphism_ring()
            Set of Morphisms
             from <sage.modules.module.Module object at ...>
               to <sage.modules.module.Module object at ...>
               in Category of modules over Integer Ring
        """
        from sage.categories.homset import End
        return End(self)


def is_Module(x):
    """
    Return ``True`` if ``x`` is a module, ``False`` otherwise.

    INPUT:

    - ``x`` -- anything

    EXAMPLES::

        sage: from sage.modules.module import is_Module
        sage: is_Module(10)
        doctest:warning...
        DeprecationWarning: the function is_Module is deprecated;
        use 'isinstance(..., Module)' instead
        See https://github.com/sagemath/sage/issues/37924 for details.
        False
        sage: M = FreeModule(RationalField(),30)                                        # needs sage.modules
        sage: is_Module(M)                                                              # needs sage.modules
        True
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(37924, "the function is_Module is deprecated; use 'isinstance(..., Module)' instead")
    return isinstance(x, Module)


def is_VectorSpace(x):
    """
    Return ``True`` if ``x`` is a vector space, ``False`` otherwise.

    INPUT:

    - ``x`` -- anything

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.modules.module import is_Module, is_VectorSpace
        sage: M = FreeModule(RationalField(),30)
        sage: is_VectorSpace(M)
        doctest:warning...
        DeprecationWarning: the function is_VectorSpace is deprecated;
        use 'isinstance(..., Module)' and check the base ring instead
        See https://github.com/sagemath/sage/issues/37924 for details.
        True
        sage: M = FreeModule(IntegerRing(),30)
        sage: is_Module(M)
        doctest:warning...
        DeprecationWarning: the function is_Module is deprecated;
        use 'isinstance(..., Module)' instead
        See https://github.com/sagemath/sage/issues/37924 for details.
        True
        sage: is_VectorSpace(M)
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(37924, "the function is_VectorSpace is deprecated; use 'isinstance(..., Module)' and check the base ring instead")
    try:
        return isinstance(x, Module) and x.base_ring().is_field()
    except AttributeError:
        return False
