"""
Elements of Free Monoids

AUTHORS:

- David Kohel (2005-09-29)

Elements of free monoids are represented internally as lists of
pairs of integers.
"""

# ****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer import Integer
from sage.structure.element import MonoidElement
from sage.structure.richcmp import richcmp, richcmp_not_equal
from sage.rings.semirings.non_negative_integer_semiring import NN


def is_FreeMonoidElement(x):
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_FreeMonoidElement is deprecated; "
                "use 'isinstance(..., FreeMonoidElement)' instead.")
    return isinstance(x, FreeMonoidElement)


class FreeMonoidElement(MonoidElement):
    """
    Element of a free monoid.

    EXAMPLES::

        sage: a = FreeMonoid(5, 'a').gens()
        sage: x = a[0]*a[1]*a[4]**3
        sage: x**3
        a0*a1*a4^3*a0*a1*a4^3*a0*a1*a4^3
        sage: x**0
        1
        sage: x**(-1)
        Traceback (most recent call last):
        ...
        NotImplementedError
    """
    def __init__(self, F, x, check=True):
        """
        Create the element `x` of the FreeMonoid `F`.

        This should typically be called by a FreeMonoid.
        """
        MonoidElement.__init__(self, F)
        if isinstance(x, (int, Integer)):
            if x == 1:
                self._element_list = []
            else:
                raise TypeError("argument x (= %s) is of the wrong type" % x)
        elif isinstance(x, list):
            if check:
                x2 = []
                for v in x:
                    if not (isinstance(v, tuple) and len(v) == 2):
                        raise TypeError("x (= %s) must be a list of 2-tuples or 1" % x)
                    if not (isinstance(v[0], (int, Integer)) and
                            isinstance(v[1], (int, Integer))):
                        raise TypeError("x (= %s) must be a list of 2-tuples of integers or 1" % x)
                    if len(x2) > 0 and v[0] == x2[len(x2)-1][0]:
                        x2[len(x2)-1] = (v[0], v[1]+x2[len(x2)-1][1])
                    else:
                        x2.append(v)
                self._element_list = x2
            else:
                self._element_list = list(x)  # make copy, so user can't accidentally change monoid.

        else:
            # TODO: should have some other checks here...
            raise TypeError("argument x (= %s) is of the wrong type" % x)

    def __hash__(self):
        r"""
        TESTS::

            sage: R.<x,y> = FreeMonoid(2)
            sage: hash(x) == hash(((0, 1),))
            True
            sage: hash(y) == hash(((1, 1),))
            True
            sage: hash(x*y) == hash(((0, 1), (1, 1)))
            True
        """
        return hash(tuple(self._element_list))

    def __iter__(self):
        """
        Return an iterator which yields tuples of variable and exponent.

        EXAMPLES::

            sage: a = FreeMonoid(5, 'a').gens()
            sage: list(a[0]*a[1]*a[4]**3*a[0])
            [(a0, 1), (a1, 1), (a4, 3), (a0, 1)]
        """
        gens = self.parent().gens()
        return ((gens[index], exponent)
                for (index, exponent) in self._element_list)

    def _repr_(self):
        s = ""
        v = self._element_list
        x = self.parent().variable_names()
        for i in range(len(v)):
            if len(s) > 0:
                s += "*"
            g = x[int(v[i][0])]
            e = v[i][1]
            if e == 1:
                s += "%s" % g
            else:
                s += f"{g}^{e}"
        if len(s) == 0:
            s = "1"
        return s

    def _latex_(self) -> str:
        r"""
        Return latex representation of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(3, 'a')
            sage: z = F([(0,5),(1,2),(0,10),(0,2),(1,2)])
            sage: z._latex_()
            'a_{0}^{5}a_{1}^{2}a_{0}^{12}a_{1}^{2}'
            sage: F.<alpha,beta,gamma> = FreeMonoid(3)
            sage: latex(alpha*beta*gamma)
            \alpha \beta \gamma

        Check that :issue:`14509` is fixed::

            sage: # needs sage.symbolic
            sage: K.< alpha,b > = FreeAlgebra(SR)
            sage: latex(alpha*b)
            \alpha b
            sage: latex(b*alpha)
            b \alpha
            sage: "%s" % latex(alpha*b)
            '\\alpha b'
        """
        s = ""
        v = self._element_list
        x = self.parent().latex_variable_names()
        for i in range(len(v)):
            g = x[int(v[i][0])]
            e = v[i][1]
            if e == 1:
                s += f"{g} "
            else:
                s += f"{g}^{{{e}}}"
        s = s.rstrip(" ")  # strip the trailing whitespace caused by adding a space after each element name
        if len(s) == 0:
            s = "1"
        return s

    def __call__(self, *x, **kwds):
        """
        EXAMPLES::

            sage: M.<x,y> = FreeMonoid(2)
            sage: (x*y).substitute(x=1)
            y

            sage: M.<a> = FreeMonoid(1)
            sage: a.substitute(a=5)
            5

            sage: M.<x,y,z> = FreeMonoid(3)
            sage: (x*y).subs(x=1,y=2,z=14)
            2
            sage: (x*y).subs({x:z,y:z})
            z^2

        It is still possible to substitute elements
        that have no common parent::

            sage: M1 = MatrixSpace(ZZ,1,2)                                              # needs sage.modules
            sage: M2 = MatrixSpace(ZZ,2,1)                                              # needs sage.modules
            sage: (x*y).subs({x: M1([1,2]), y: M2([3,4])})                              # needs sage.modules
            [11]

        TESTS::

            sage: M.<x,y> = FreeMonoid(2)
            sage: (x*y)(QQ(4),QQ(5)).parent()
            Rational Field

        The codomain is by default the first parent::

            sage: M.one()(QQ(4),QQ(5)).parent()
            Rational Field

        unless there is no variable and no substitution::

            sage: M = FreeMonoid(0, [])
            sage: M.one()().parent()
            Free monoid on 0 generators ()

        AUTHORS:

        - Joel B. Mohler (2007-10-27)
        """
        if kwds and x:
            raise ValueError("must not specify both a keyword and positional argument")

        P = self.parent()

        if kwds:
            x = self.gens()
            gens_dict = {name: i for i, name in enumerate(P.variable_names())}
            for key, value in kwds.items():
                if key in gens_dict:
                    x[gens_dict[key]] = value

        if x and isinstance(x[0], tuple):
            x = x[0]

        if len(x) != self.parent().ngens():
            raise ValueError("must specify as many values as generators in parent")

        # if no substitution, do nothing
        if not x:
            return self

        try:
            # This will land in the parent of the first element
            result = x[0].parent().one()
        except (AttributeError, TypeError):
            # unless the parent has no unit
            result = NN.one()
        for var_index, exponent in self._element_list:
            replacement = x[var_index]
            if exponent > 1:
                result *= replacement ** exponent
            elif exponent == 1:
                result *= replacement
        return result

    def _mul_(self, y):
        """
        Multiply two elements ``self`` and ``y`` of the
        free monoid.

        EXAMPLES::

            sage: a = FreeMonoid(5, 'a').gens()
            sage: x = a[0] * a[1] * a[4]**3
            sage: y = a[4] * a[0] * a[1]
            sage: x*y
            a0*a1*a4^4*a0*a1
        """
        M = self.parent()
        z = M(1)
        x_elt = self._element_list
        y_elt = y._element_list
        if not x_elt:
            z._element_list = y_elt
        elif not y_elt:
            z._element_list = x_elt
        else:
            k = len(x_elt)-1
            if x_elt[k][0] != y_elt[0][0]:
                z._element_list = x_elt + y_elt
            else:
                m = (y_elt[0][0], x_elt[k][1]+y_elt[0][1])
                z._element_list = x_elt[:k] + [m] + y_elt[1:]
        return z

    def __invert__(self):
        """
        EXAMPLES::

            sage: a = FreeMonoid(5, 'a').gens()
            sage: x = a[0]*a[1]*a[4]**3
            sage: x**(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __len__(self) -> int:
        """
        Return the degree of the monoid element ``self``, where each
        generator of the free monoid is given degree `1`.

        For example, the length of the identity is `0`, and the
        length of `x_0^2x_1` is `3`.

        EXAMPLES::

            sage: F = FreeMonoid(3, 'a')
            sage: z = F(1)
            sage: len(z)
            0
            sage: a = F.gens()
            sage: len(a[0]**2 * a[1])
            3
        """
        return sum(x[1] for x in self._element_list)

    def _richcmp_(self, other, op) -> bool:
        """
        Compare two free monoid elements with the same parents.

        The ordering is first by increasing length, then lexicographically
        on the underlying word.

        EXAMPLES::

            sage: S = FreeMonoid(3, 'a')
            sage: (x,y,z) = S.gens()
            sage: x * y < y * x
            True

            sage: a = FreeMonoid(5, 'a').gens()
            sage: x = a[0]*a[1]*a[4]**3
            sage: x < x
            False
            sage: x == x
            True
            sage: x >= x*x
            False
        """
        m = sum(i for x, i in self._element_list)
        n = sum(i for x, i in other._element_list)
        if m != n:
            return richcmp_not_equal(m, n, op)
        v = tuple([x for x, i in self._element_list for j in range(i)])
        w = tuple([x for x, i in other._element_list for j in range(i)])
        return richcmp(v, w, op)

    def _acted_upon_(self, x, self_on_left):
        """
        Return the action of the integer 1 on this element.

        EXAMPLES::

            sage: M.<x,y,z>=FreeMonoid(3)
            sage: 1*x
            x
        """
        if x == 1:
            return self
        return None

    def to_word(self, alph=None):
        """
        Return ``self`` as a word.

        INPUT:

        - ``alph`` -- (optional) the alphabet which the result should
          be specified in

        EXAMPLES::

            sage: M.<x,y,z> = FreeMonoid(3)
            sage: a = x * x * y * x
            sage: w = a.to_word(); w
            word: xxyx
            sage: w.to_monoid_element() == a
            True

        .. SEEALSO::

            :meth:`to_list`
        """
        from sage.combinat.words.finite_word import Words
        gens = self.parent().gens()
        if alph is None:
            alph = gens
        alph = [str(c) for c in alph]
        W = Words(alph, infinite=False)
        return W(sum([[alph[gens.index(i[0])]] * i[1] for i in self], []))

    def to_list(self, indices=False) -> list:
        r"""
        Return ``self`` as a list of generators.

        If ``self`` equals `x_{i_1} x_{i_2} \cdots x_{i_n}`, with
        `x_{i_1}, x_{i_2}, \ldots, x_{i_n}` being some of the
        generators of the free monoid, then this method returns
        the list `[x_{i_1}, x_{i_2}, \ldots, x_{i_n}]`.

        If the optional argument ``indices`` is set to ``True``,
        then the list `[i_1, i_2, \ldots, i_n]` is returned instead.

        EXAMPLES::

            sage: M.<x,y,z> = FreeMonoid(3)
            sage: a = x * x * y * x
            sage: w = a.to_list(); w
            [x, x, y, x]
            sage: M.prod(w) == a
            True
            sage: w = a.to_list(indices=True); w
            [0, 0, 1, 0]
            sage: a = M.one()
            sage: a.to_list()
            []

        .. SEEALSO::

            :meth:`to_word`
        """
        if not indices:
            return sum(([i[0]] * i[1] for i in list(self)), [])
        gens = self.parent().gens()
        return sum(([gens.index(i[0])] * i[1] for i in list(self)), [])
