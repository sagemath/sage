"""
Degeneracy maps
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from . import morphism

class DegeneracyMap(morphism.HeckeModuleMorphism_matrix):
    """
    A degeneracy map between Hecke modules of different levels.

    EXAMPLES:

    We construct a number of degeneracy maps::

        sage: M = ModularSymbols(33)
        sage: d = M.degeneracy_map(11)
        sage: d
        Hecke module morphism degeneracy map corresponding to f(q) |--> f(q) defined by the matrix
        [ 1  0  0]
        [ 0  0  1]
        [ 0  0 -1]
        [ 0  1 -1]
        [ 0  0  1]
        [ 0 -1  1]
        [-1  0  0]
        [-1  0  0]
        [-1  0  0]
        Domain: Modular Symbols space of dimension 9 for Gamma_0(33) of weight ...
        Codomain: Modular Symbols space of dimension 3 for Gamma_0(11) of weight ...
        sage: d.t()
        1
        sage: d = M.degeneracy_map(11,3)
        sage: d.t()
        3

    The parameter d must be a divisor of the quotient of the two levels::

        sage: d = M.degeneracy_map(11,2)
        Traceback (most recent call last):
        ...
        ValueError: The level of self (=33) must be a divisor or multiple of level (=11), and t (=2) must be a divisor of the quotient.

    Degeneracy maps can also go from lower level to higher level::

        sage: M.degeneracy_map(66,2)
        Hecke module morphism degeneracy map corresponding to f(q) |--> f(q^2) defined by the matrix
        [ 2  0  0  0  0  0  1  0  0  0  1 -1  0  0  0 -1  1  0  0  0  0  0  0  0 -1]
        [ 0  0  1 -1  0 -1  1  0 -1  2  0  0  0 -1  0  0 -1  1  2 -2  0  0  0 -1  1]
        [ 0  0  1  0  0  0  0  0  1  0  0  0  1  0  0  0 -1  1  0  0 -1  1  0  0  0]
        [ 0  0  0  0  0  0  0  0  0  2 -1  0  0  1  0  0 -1  1  0  0  1  0 -1 -1  1]
        [ 0 -1  0  0  1  0  0  0  0  0  0  1  0  0  1  1 -1  0  0 -1  0  0  0  0  0]
        [ 0  0  0  0  0  0  0  1 -1  0  0  2 -1  0  0  1  0  0  0 -1  0 -1  1 -1  1]
        [ 0  0  0  0  1 -1  0  1 -1  0  0  0  0  0 -1  2  0  0  0  0  1  0  1  0  0]
        [ 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1  0  0]
        [ 0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  1  1  1  0  0  0]
        Domain: Modular Symbols space of dimension 9 for Gamma_0(33) of weight ...
        Codomain: Modular Symbols space of dimension 25 for Gamma_0(66) of weight ...
    """
    def __init__(self, matrix, domain, codomain, t):
        r"""
        Initialise a degeneracy map.

        EXAMPLES::

            sage: D = ModularSymbols(Gamma0(100)).degeneracy_map(2,5); D
            Hecke module morphism degeneracy map corresponding to f(q) |--> f(q^5) defined by the matrix
            31 x 1 dense matrix over Rational Field
            Domain: Modular Symbols space of dimension 31 for Gamma_0(100) of weight ...
            Codomain: Modular Symbols space of dimension 1 for Gamma_0(2) of weight ...
            sage: D == loads(dumps(D))
            True
        """
        self.__t = t
        H = domain.Hom(codomain)
        if t == 1:
            pow = ""
        else:
            pow = "^%s"%t
        name = "degeneracy map corresponding to f(q) |--> f(q%s)"%(pow)
        morphism.HeckeModuleMorphism_matrix.__init__(self, H, matrix, name)

    def t(self):
        """
        Return the divisor of the quotient of the two levels
        associated to the degeneracy map.

        EXAMPLES::

            sage: M = ModularSymbols(33)
            sage: d = M.degeneracy_map(11,3)
            sage: d.t()
            3
            sage: d = M.degeneracy_map(11,1)
            sage: d.t()
            1
        """
        return self.__t
