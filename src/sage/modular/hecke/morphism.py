# sage.doctest: needs sage.libs.flint
"""
Morphisms of Hecke modules

AUTHORS:

- William Stein
"""

# ****************************************************************************
#       Sage: Open Source Mathematical Software
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import sage.misc.misc as misc
from sage.modules.matrix_morphism import MatrixMorphism
from sage.categories.morphism import Morphism

# We also define other types of Hecke-module morphisms that aren't
# specified by a matrix.  E.g., Hecke operators, or maybe morphisms on
# modular abelian varieties (which are specified by matrices, but on
# integral homology). All morphisms derive from HeckeModuleMorphism.


def is_HeckeModuleMorphism(x):
    r"""
    Return ``True`` if x is of type HeckeModuleMorphism.

    EXAMPLES::

        sage: sage.modular.hecke.morphism.is_HeckeModuleMorphism(ModularSymbols(6).hecke_operator(7).hecke_module_morphism())
        doctest:warning...
        DeprecationWarning: the function is_HeckeModuleMorphism is deprecated;
        use 'isinstance(..., HeckeModuleMorphism)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_HeckeModuleMorphism is deprecated; use 'isinstance(..., HeckeModuleMorphism)' instead")
    return isinstance(x, HeckeModuleMorphism)


def is_HeckeModuleMorphism_matrix(x):
    """

    EXAMPLES::

        sage: sage.modular.hecke.morphism.is_HeckeModuleMorphism_matrix(ModularSymbols(6).hecke_operator(7).matrix_form().hecke_module_morphism())
        doctest:warning...
        DeprecationWarning: the function is_HeckeModuleMorphism_matrix is deprecated;
        use 'isinstance(..., HeckeModuleMorphism_matrix)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_HeckeModuleMorphism_matrix is deprecated; use 'isinstance(..., HeckeModuleMorphism_matrix)' instead")
    return isinstance(x, HeckeModuleMorphism_matrix)


class HeckeModuleMorphism(Morphism):
    r"""
    Abstract base class for morphisms of Hecke modules.
    """
    pass


class HeckeModuleMorphism_matrix(MatrixMorphism, HeckeModuleMorphism):
    """
    Morphisms of Hecke modules when the morphism is given by a matrix.

    Note that care is needed when composing morphisms, because morphisms in
    Sage act on the left, but their matrices act on the right (!). So if F: A
    -> B and G : B -> C are morphisms, the composition A -> C is G*F, but its
    matrix is F.matrix() * G.matrix().

    EXAMPLES::

        sage: A = ModularForms(1, 4)
        sage: B = ModularForms(1, 16)
        sage: C = ModularForms(1, 28)
        sage: F = A.Hom(B)(matrix(QQ,1,2,srange(1, 3)))
        sage: G = B.Hom(C)(matrix(QQ,2,3,srange(1, 7)))
        sage: G * F
        Hecke module morphism defined by the matrix
        [ 9 12 15]
        Domain: Modular Forms space of dimension 1 for Modular Group SL(2,Z) ...
        Codomain: Modular Forms space of dimension 3 for Modular Group SL(2,Z) ...
        sage: F * G
        Traceback (most recent call last):
        ...
        TypeError: Incompatible composition of morphisms: domain of left morphism must be codomain of right.
    """
    def __init__(self, parent, A, name='', side='left'):
        """
        INPUT:

        -  ``parent`` -- ModularSymbolsHomspace

        - ``A`` -- matrix

        - ``name`` -- string (default: ``''``); name of the morphism
           (used for printing)

        - ``side`` -- string (default: ``'left'``)

        EXAMPLES::

            sage: M = ModularSymbols(6)
            sage: t = M.Hom(M)(matrix(QQ,3,3,srange(9)), name='spam'); t
            Hecke module morphism spam defined by the matrix
            [0 1 2]
            [3 4 5]
            [6 7 8]
            Domain: Modular Symbols space of dimension 3 for Gamma_0(6) of weight ...
            Codomain: Modular Symbols space of dimension 3 for Gamma_0(6) of weight ...
            sage: t == loads(dumps(t))
            True
        """
        if not isinstance(name, str):
            raise TypeError("name must be a string")
        self.__name = name
        MatrixMorphism.__init__(self, parent, A, side)

    def name(self, new=None):
        r"""
        Return the name of this operator, or set it to a new name.

        EXAMPLES::

            sage: M = ModularSymbols(6)
            sage: t = M.Hom(M)(matrix(QQ,3,3,srange(9)), name='spam'); t
            Hecke module morphism spam defined by ...
            sage: t.name()
            'spam'
            sage: t.name("eggs"); t
            Hecke module morphism eggs defined by ...
        """
        if new is None:
            return self.__name
        self.__name = new

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(6)
            sage: t = M.Hom(M)(matrix(QQ,3,3,srange(9))); t._repr_()
            'Hecke module morphism defined by the matrix\n[0 1 2]\n[3 4 5]\n[6 7 8]\nDomain: Modular Symbols space of dimension 3 for Gamma_0(6) of weight ...\nCodomain: Modular Symbols space of dimension 3 for Gamma_0(6) of weight ...'
            sage: t.name('spam'); t._repr_()
            'Hecke module morphism spam defined by the matrix\n[0 1 2]\n[3 4 5]\n[6 7 8]\nDomain: Modular Symbols space of dimension 3 for Gamma_0(6) of weight ...\nCodomain: Modular Symbols space of dimension 3 for Gamma_0(6) of weight ...'
        """
        name = self.__name
        if name:
            name += ' '
        return "Hecke module morphism %sdefined by the matrix\n%r\nDomain: %s\nCodomain: %s" % (name, self.matrix(), misc.strunc(self.domain()), misc.strunc(self.codomain()))

# __mul__ method removed by David Loeffler 2009-04-14 as it is an exact duplicate of sage.modules.matrix_morphism.__mul__
