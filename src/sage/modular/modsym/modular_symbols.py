r"""
Modular symbols {alpha, beta}

The ModularSymbol class represents a single modular symbol `X^i Y^{k-2-i} \{\alpha, \beta\}`.

AUTHOR:

- William Stein (2005, 2009)

TESTS::

    sage: s = ModularSymbols(11).2.modular_symbol_rep()[0][1]; s
    {-1/9, 0}
    sage: loads(dumps(s)) == s
    True
"""

#*****************************************************************************
#       Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005, 2009 William Stein <wstein@gmail.com>
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
import sage.modular.cusps as cusps
from sage.modular.modsym.apply import apply_to_monomial
from sage.modular.modsym.manin_symbol import ManinSymbol
from sage.structure.sage_object import SageObject
import sage.structure.formal_sum as formal_sum
from sage.structure.richcmp import richcmp_method, richcmp
from sage.rings.integer_ring import ZZ
from sage.misc.latex import latex

_C = cusps.Cusps

X, Y = ZZ['X,Y'].gens()


@richcmp_method
class ModularSymbol(SageObject):
    r"""
    The modular symbol `X^i\cdot Y^{k-2-i}\cdot \{\alpha, \beta\}`.
    """
    def __init__(self, space, i, alpha, beta):
        """
        Initialise a modular symbol.

        INPUT:

        - ``space`` -- space of Manin symbols

        - ``i`` -- integer

        - ``alpha`` -- cusp

        - ``beta`` -- cusp

        EXAMPLES::

            sage: s = ModularSymbols(11).2.modular_symbol_rep()[0][1]; s
            {-1/9, 0}
            sage: type(s)
            <class 'sage.modular.modsym.modular_symbols.ModularSymbol'>
            sage: s = ModularSymbols(11,4).2.modular_symbol_rep()[0][1]; s
            X^2*{-1/7, 0}
        """
        self.__space = space
        self.__i = i
        self.__alpha = _C(alpha)
        self.__beta = _C(beta)

    def _repr_(self):
        """
        String representation of this modular symbol.

        EXAMPLES::

            sage: s = ModularSymbols(11,4).2.modular_symbol_rep()[0][1]; s
            X^2*{-1/7, 0}
            sage: s._repr_()
            'X^2*{-1/7, 0}'
            sage: s.rename('sym')
            sage: s
            sym
        """
        if self.weight() == 2:
            polypart = ''
        else:
            polypart = str(self.polynomial_part()) + '*'
        return "%s{%s, %s}"%(polypart, self.__alpha, self.__beta)

    def __getitem__(self, j):
        r"""
        Given a modular symbols `s = X^i Y^{k-2-i}\{\alpha, \beta\}`, ``s[0]`` is `\alpha`
        and ``s[1]`` is `\beta`.

        EXAMPLES::

            sage: s = ModularSymbols(11).2.modular_symbol_rep()[0][1]; s
            {-1/9, 0}
            sage: s[0]
            -1/9
            sage: s[1]
            0
            sage: s[2]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return [self.__alpha, self.__beta][j]

    def _latex_(self):
        r"""
        Return Latex representation of this modular symbol.

        EXAMPLES::

            sage: s = ModularSymbols(11,4).2.modular_symbol_rep()[0][1]; s
            X^2*{-1/7, 0}
            sage: latex(s)                         # indirect doctest
            X^{2}\left\{\frac{-1}{7}, 0\right\}
        """
        if self.weight() == 2:
            polypart = ''
        else:
            polypart = latex(self.polynomial_part())
        return "%s\\left\\{%s, %s\\right\\}"%(polypart,
                  latex(self.__alpha), latex(self.__beta))

    def __richcmp__(self, other, op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: s = M.2.modular_symbol_rep()[0][1]
            sage: t = M.0.modular_symbol_rep()[0][1]
            sage: s, t
            ({-1/9, 0}, {Infinity, 0})
            sage: s < t
            True
            sage: t > s
            True
            sage: s == s
            True
            sage: t == t
            True
        """
        if not isinstance(other, ModularSymbol):
            return NotImplemented
        return richcmp((self.__space, -self.__i, self.__alpha, self.__beta),
                       (other.__space,-other.__i,other.__alpha,other.__beta),
                       op)

    def __hash__(self):
        """
        EXAMPLES::

            sage: s = ModularSymbols(11).2.modular_symbol_rep()[0][1]
            sage: hash(s)  # random
            -7344656798833624820
        """
        return hash((self.__space, self.__i, self.__alpha, self.__beta))

    def space(self):
        """
        The list of Manin symbols to which this symbol belongs.

        EXAMPLES::

            sage: s = ModularSymbols(11).2.modular_symbol_rep()[0][1]
            sage: s.space()
            Manin Symbol List of weight 2 for Gamma0(11)
        """
        return self.__space

    def polynomial_part(self):
        r"""
        Return the polynomial part of this symbol, i.e. for a symbol of the
        form `X^i Y^{k-2-i}\{\alpha, \beta\}`, return `X^i Y^{k-2-i}`.

        EXAMPLES::

            sage: s = ModularSymbols(11).2.modular_symbol_rep()[0][1]
            sage: s.polynomial_part()
            1
            sage: s = ModularSymbols(1,28).0.modular_symbol_rep()[0][1]; s
            X^22*Y^4*{0, Infinity}
            sage: s.polynomial_part()
            X^22*Y^4
        """
        i = self.__i
        return X**i*Y**(self.weight()-2-i)

    def i(self):
        r"""
        For a symbol of the form `X^i Y^{k-2-i}\{\alpha, \beta\}`, return `i`.

        EXAMPLES::

            sage: s = ModularSymbols(11).2.modular_symbol_rep()[0][1]
            sage: s.i()
            0
            sage: s = ModularSymbols(1,28).0.modular_symbol_rep()[0][1]; s
            X^22*Y^4*{0, Infinity}
            sage: s.i()
            22
        """
        return self.__i

    def weight(self):
        r"""
        Return the weight of the modular symbols space to which this symbol
        belongs; i.e. for a symbol of the form `X^i Y^{k-2-i}\{\alpha,
        \beta\}`, return `k`.

        EXAMPLES::

            sage: s = ModularSymbols(1,28).0.modular_symbol_rep()[0][1]
            sage: s.weight()
            28
        """
        return self.__space.weight()

    def alpha(self):
        r"""
        For a symbol of the form `X^i Y^{k-2-i}\{\alpha, \beta\}`, return `\alpha`.

        EXAMPLES::

            sage: s = ModularSymbols(11,4).1.modular_symbol_rep()[0][1]; s
            X^2*{-1/6, 0}
            sage: s.alpha()
            -1/6
            sage: type(s.alpha())
            <class 'sage.modular.cusps.Cusp'>
        """
        return self.__alpha

    def beta(self):
        r"""
        For a symbol of the form `X^i Y^{k-2-i}\{\alpha, \beta\}`, return `\beta`.

        EXAMPLES::

            sage: s = ModularSymbols(11,4).1.modular_symbol_rep()[0][1]; s
            X^2*{-1/6, 0}
            sage: s.beta()
            0
            sage: type(s.beta())
            <class 'sage.modular.cusps.Cusp'>
        """
        return self.__beta

    def apply(self, g):
        r"""
        Act on this symbol by the element `g \in {\rm GL}_2(\QQ)`.

        INPUT:

        - ``g`` -- a list ``[a,b,c,d]``, corresponding to the 2x2 matrix
          `\begin{pmatrix} a & b \\ c & d \end{pmatrix} \in {\rm GL}_2(\QQ)`.

        OUTPUT:

        - ``FormalSum`` -- a formal sum `\sum_i c_i x_i`, where `c_i` are
          scalars and `x_i` are ModularSymbol objects, such that the sum
          `\sum_i c_i x_i` is the image of this symbol under the action of g.
          No reduction is performed modulo the relations that hold in
          self.space().

        The action of `g` on symbols is by

        .. MATH::

           P(X,Y)\{\alpha, \beta\} \mapsto  P(dX-bY, -cx+aY) \{g(\alpha), g(\beta)\}.

        Note that for us we have `P=X^i Y^{k-2-i}`, which simplifies computation
        of the polynomial part slightly.

        EXAMPLES::

            sage: s = ModularSymbols(11,2).1.modular_symbol_rep()[0][1]; s
            {-1/8, 0}
            sage: a=1;b=2;c=3;d=4; s.apply([a,b,c,d])
            {15/29, 1/2}
            sage: x = -1/8;  (a*x+b)/(c*x+d)
            15/29
            sage: x = 0;  (a*x+b)/(c*x+d)
            1/2
            sage: s = ModularSymbols(11,4).1.modular_symbol_rep()[0][1]; s
            X^2*{-1/6, 0}
            sage: s.apply([a,b,c,d])
            16*X^2*{11/21, 1/2} - 16*X*Y*{11/21, 1/2} + 4*Y^2*{11/21, 1/2}
            sage: P = s.polynomial_part()
            sage: X,Y = P.parent().gens()
            sage: P(d*X-b*Y, -c*X+a*Y)
            16*X^2 - 16*X*Y + 4*Y^2
            sage: x=-1/6; (a*x+b)/(c*x+d)
            11/21
            sage: x=0; (a*x+b)/(c*x+d)
            1/2
            sage: type(s.apply([a,b,c,d]))
            <class 'sage.structure.formal_sum.FormalSum'>
        """
        space = self.__space
        i = self.__i
        k = space.weight()
        a,b,c,d = tuple(g)
        coeffs = apply_to_monomial(i, k-2, d, -b, -c, a)
        g_alpha = self.__alpha.apply(g)
        g_beta = self.__beta.apply(g)
        return formal_sum.FormalSum([(coeffs[j], ModularSymbol(space, j, g_alpha, g_beta)) \
                                     for j in reversed(range(k-1)) if coeffs[j] != 0])

    def __manin_symbol_rep(self, alpha):
        """
        Return Manin symbol representation of X^i*Y^(k-2-i){0,alpha}.

        EXAMPLES::

            sage: s = ModularSymbols(11,2).1.modular_symbol_rep()[0][1]; s
            {-1/8, 0}
            sage: s.manin_symbol_rep()          # indirect doctest
            -(1,1) - (-8,1)
            sage: M = ModularSymbols(11,2)
            sage: s = M( (1,9) ); s
            (1,9)
            sage: t = s.modular_symbol_rep()[0][1].manin_symbol_rep(); t
            -(1,1) - (-9,1)
            sage: M(t)
            (1,9)
        """
        space = self.__space
        i = self.__i
        k = space.weight()
        v = [(0,1), (1,0)]
        if not alpha.is_infinity():
            cf = alpha._rational_().continued_fraction()
            v.extend((cf.p(k),cf.q(k)) for k in range(len(cf)))
        sign = 1
        z = formal_sum.FormalSum(0)
        for j in range(1,len(v)):
            c = sign*v[j][1]
            d = v[j-1][1]
            coeffs = apply_to_monomial(i, k-2, sign*v[j][0], v[j-1][0],
                                       sign*v[j][1], v[j-1][1])
            w = [(coeffs[j], ManinSymbol(space, (j, c, d)))
                 for j in range(k-1) if coeffs[j] != 0]
            z += formal_sum.FormalSum(w)
            sign *= -1
        return z

    def manin_symbol_rep(self):
        """
        Return a representation of ``self`` as a formal sum of Manin symbols.

        The result is not cached.

        EXAMPLES::

            sage: M = ModularSymbols(11,4)
            sage: s = M.1.modular_symbol_rep()[0][1]; s
            X^2*{-1/6, 0}
            sage: s.manin_symbol_rep()
            -2*[X*Y,(-1,0)] - [X^2,(-1,0)] - [Y^2,(1,1)] - [X^2,(-6,1)]
            sage: M(s.manin_symbol_rep()) == M([2,-1/6,0])
            True
        """
        alpha = self.__alpha
        beta = self.__beta
        return -1*self.__manin_symbol_rep(alpha) + self.__manin_symbol_rep(beta)

