# sage.doctest: needs sage.libs.flint
"""
A single element of an ambient space of modular symbols
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

import sage.modules.free_module_element
from sage.misc.repr import repr_lincomb
import sage.structure.formal_sum as formal_sum
import sage.modular.hecke.all as hecke
import sage.misc.latex as latex


_print_mode = "manin"


def is_ModularSymbolsElement(x) -> bool:
    r"""
    Return ``True`` if x is an element of a modular symbols space.

    EXAMPLES::

        sage: sage.modular.modsym.element.is_ModularSymbolsElement(ModularSymbols(11, 2).0)
        doctest:warning...
        DeprecationWarning: The function is_ModularSymbolsElement is deprecated;
        use 'isinstance(..., ModularSymbolsElement)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        True
        sage: sage.modular.modsym.element.is_ModularSymbolsElement(13)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_ModularSymbolsElement is deprecated; "
                "use 'isinstance(..., ModularSymbolsElement)' instead.")
    return isinstance(x, ModularSymbolsElement)


def set_modsym_print_mode(mode='manin'):
    r"""
    Set the mode for printing of elements of modular symbols spaces.

    INPUT:

    - ``mode`` -- string; the possibilities are as follows:

      - ``'manin'`` -- (the default) formal sums of Manin
        symbols [P(X,Y),(u,v)]

      - ``'modular'`` -- formal sums of Modular symbols
        P(X,Y)\*alpha,beta, where alpha and beta are cusps

      - ``'vector'`` -- as vectors on the basis for the
        ambient space

    OUTPUT: none

    EXAMPLES::

        sage: M = ModularSymbols(13, 8)
        sage: x = M.0 + M.1 + M.14
        sage: set_modsym_print_mode('manin'); x
        [X^5*Y,(1,11)] + [X^5*Y,(1,12)] + [X^6,(1,11)]
        sage: set_modsym_print_mode('modular'); x
        1610510*X^6*{-1/11, 0} + 893101*X^5*Y*{-1/11, 0} + 206305*X^4*Y^2*{-1/11, 0} + 25410*X^3*Y^3*{-1/11, 0} + 1760*X^2*Y^4*{-1/11, 0} + 65*X*Y^5*{-1/11, 0} - 248832*X^6*{-1/12, 0} - 103680*X^5*Y*{-1/12, 0} - 17280*X^4*Y^2*{-1/12, 0} - 1440*X^3*Y^3*{-1/12, 0} - 60*X^2*Y^4*{-1/12, 0} - X*Y^5*{-1/12, 0} + Y^6*{-1/11, 0}
        sage: set_modsym_print_mode('vector'); x
        (1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        sage: set_modsym_print_mode()
    """
    mode = str(mode).lower()
    if mode not in ['manin', 'modular', 'vector']:
        raise ValueError("mode must be one of 'manin', 'modular', or 'vector'")
    global _print_mode
    _print_mode = mode


class ModularSymbolsElement(hecke.HeckeModuleElement):
    """
    An element of a space of modular symbols.

    TESTS::

        sage: x = ModularSymbols(3, 12).cuspidal_submodule().gen(0)
        sage: x == loads(dumps(x))
        True
    """
    def __init__(self, parent, x, check=True):
        """
        INPUT:

        - ``parent`` -- a space of modular symbols

        - ``x`` -- a free module element that represents the modular
          symbol in terms of a basis for the ambient space (not in
          terms of a basis for parent!)

        EXAMPLES::

            sage: S = ModularSymbols(11, sign=1).cuspidal_submodule()
            sage: S(vector([0,1])) == S.basis()[0]
            True
            sage: S(vector([1,0]))
            Traceback (most recent call last):
            ...
            TypeError: x does not coerce to an element of this Hecke module
        """
        if check:
            from .space import ModularSymbolsSpace
            if not isinstance(parent, ModularSymbolsSpace):
                raise TypeError("parent (= %s) must be a space of modular symbols" % parent)
            if not isinstance(x, sage.modules.free_module_element.FreeModuleElement):
                raise TypeError("x must be a free module element.")
            if x.degree() != parent.degree():
                raise TypeError("x (of degree %s) must be of degree the same as the degree of the parent (of degree %s)." % (x.degree(), parent.degree()))
        hecke.HeckeModuleElement.__init__(self, parent, x)

    def _repr_(self):
        r"""
        String representation of ``self``. The output will depend on the global
        modular symbols print mode setting controlled by the function
        ``set_modsym_print_mode``.

        EXAMPLES::

            sage: M = ModularSymbols(13, 4)
            sage: set_modsym_print_mode('manin'); M.0._repr_()
            '[X^2,(0,1)]'
            sage: set_modsym_print_mode('modular'); M.0._repr_()
            'X^2*{0, Infinity}'
            sage: set_modsym_print_mode('vector'); M.0._repr_()
            '(1, 0, 0, 0, 0, 0, 0, 0)'
            sage: set_modsym_print_mode()
        """
        if _print_mode == "vector":
            return str(self.element())
        elif _print_mode == "manin":
            m = self.manin_symbol_rep()
        elif _print_mode == "modular":
            m = self.modular_symbol_rep()
        return repr_lincomb([(t, c) for c, t in m])

    def _latex_(self):
        r"""
        LaTeX representation of ``self``. The output will be determined by the
        print mode setting set using ``set_modsym_print_mode``.

        EXAMPLES::

            sage: M = ModularSymbols(11, 2)
            sage: x = M.0 + M.2; x
            (1,0) + (1,9)
            sage: set_modsym_print_mode('manin'); latex(x) # indirect doctest
            (1,0) + (1,9)
            sage: set_modsym_print_mode('modular'); latex(x) # indirect doctest
            \left\{\infty, 0\right\} + \left\{\frac{-1}{9}, 0\right\}
            sage: set_modsym_print_mode('vector'); latex(x) # indirect doctest
            \left(1,\,0,\,1\right)
            sage: set_modsym_print_mode()
        """

        if _print_mode == "vector":
            return self.element()._latex_()
        elif _print_mode == "manin":
            m = self.manin_symbol_rep()
        elif _print_mode == "modular":
            m = self.modular_symbol_rep()
        c = [x[0] for x in m]
        v = [x[1] for x in m]
        # TODO: use repr_lincomb with is_latex=True
        return latex.repr_lincomb(v, c)

    def _add_(self, right):
        r"""
        Sum of ``self`` and ``other``.

        EXAMPLES::

            sage: M = ModularSymbols(3, 12)
            sage: x = M.0; y = M.1; z = x + y; z # indirect doctest
            [X^8*Y^2,(1,2)] + [X^9*Y,(1,0)]
            sage: z.parent() is M
            True
        """
        return ModularSymbolsElement(self.parent(), self.element() + right.element(), check=False)

    def _rmul_(self, other):
        r"""
        Right-multiply ``self`` by ``other``.

        EXAMPLES::

            sage: M = ModularSymbols(3, 12)
            sage: x = M.0; z = x*3; z # indirect doctest
            3*[X^8*Y^2,(1,2)]
            sage: z.parent() is M
            True
            sage: z*Mod(1, 17)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Modular Symbols space of dimension 8 for Gamma_0(3) of weight 12 with sign 0 over Rational Field' and 'Ring of integers modulo 17'
        """
        return ModularSymbolsElement(self.parent(), self.element() * other,
                                     check=False)

    def _lmul_(self, left):
        r"""
        Left-multiply ``self`` by ``other``.

        EXAMPLES::

            sage: M = ModularSymbols(3, 12)
            sage: x = M.0; z = 3*x; z # indirect doctest
            3*[X^8*Y^2,(1,2)]
            sage: z.parent() is M
            True
            sage: Mod(1, 17)*z
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Ring of integers modulo 17' and 'Modular Symbols space of dimension 8 for Gamma_0(3) of weight 12 with sign 0 over Rational Field'
        """
        return ModularSymbolsElement(self.parent(), left * self.element(),
                                     check=False)

    def _neg_(self):
        r"""
        Multiply by -1.

        EXAMPLES::

            sage: M = ModularSymbols(3, 12)
            sage: x = M.0; z = -x; z # indirect doctest
            -[X^8*Y^2,(1,2)]
            sage: z.parent() is M
            True
        """
        return ModularSymbolsElement(self.parent(), -self.element(), check=False)

    def _sub_(self, other):
        r"""
        Subtract ``other`` from ``self``.

        EXAMPLES::

            sage: M = ModularSymbols(3, 12)
            sage: x = M.0; y = M.1; z = y-x; z # indirect doctest
            -[X^8*Y^2,(1,2)] + [X^9*Y,(1,0)]
            sage: z.parent() is M
            True
        """
        return ModularSymbolsElement(self.parent(), self.element() - other.element(), check=False)

#   this clearly hasn't worked for some time -- the method embedded_vector_space doesn't exist -- DL 2009-05-18
#    def coordinate_vector(self):
#        if self.parent().is_ambient():
#            return self.element()
#        return self.parent().embedded_vector_space().coordinate_vector(self.element())

    def list(self):
        r"""
        Return a list of the coordinates of ``self`` in terms of a basis for
        the ambient space.

        EXAMPLES::

            sage: ModularSymbols(37, 2).0.list()
            [1, 0, 0, 0, 0]
        """
        return self.element().list()

    def manin_symbol_rep(self):
        """
        Return a representation of ``self`` as a formal sum of Manin symbols.

        EXAMPLES::

            sage: x = ModularSymbols(37, 4).0
            sage: x.manin_symbol_rep()
            [X^2,(0,1)]

        The result is cached::

            sage: x.manin_symbol_rep() is x.manin_symbol_rep()
            True
        """
        try:
            return self.__manin_symbols
        except AttributeError:
            A = self.parent()
            v = self.element()
            manin_symbols = A.ambient_hecke_module().manin_symbols_basis()
            F = formal_sum.FormalSums(A.base_ring())
            ms = F([(v[i], manin_symbols[i]) for i in range(v.degree())
                    if v[i] != 0], check=False, reduce=False)
            self.__manin_symbols = ms
        return self.__manin_symbols

    def modular_symbol_rep(self):
        """
        Return a representation of ``self`` as a formal sum of modular symbols.

        EXAMPLES::

            sage: x = ModularSymbols(37, 4).0
            sage: x.modular_symbol_rep()
            X^2*{0, Infinity}

        The result is cached::

            sage: x.modular_symbol_rep() is x.modular_symbol_rep()
            True
        """
        try:
            return self.__modular_symbols
        except AttributeError:
            v = self.manin_symbol_rep()
            if v == 0:
                return v
            w = [c * x.modular_symbol_rep() for c, x in v]
            self.__modular_symbols = sum(w)
            return self.__modular_symbols
