r"""
The modular group `\SL_2(\ZZ)`

AUTHORS:

- Niles Johnson (2010-08): :issue:`3893`: ``random_element()`` should pass on ``*args`` and ``**kwds``.
"""

################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
################################################################################

from sage.arith.misc import gcd
from sage.modular.cusps import Cusp
from sage.modular.modsym.p1list import lift_to_sl2z
from sage.rings.integer_ring import ZZ

from .congroup_gamma0 import Gamma0_class
from .arithgroup_element import ArithmeticSubgroupElement


def is_SL2Z(x):
    r"""
    Return ``True`` if x is the modular group `\SL_2(\ZZ)`.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import is_SL2Z
        sage: is_SL2Z(SL2Z)
        doctest:warning...
        DeprecationWarning: The function is_SL2Z is deprecated; use 'isinstance(..., SL2Z_class)' instead.
        See https://github.com/sagemath/sage/issues/38035 for details.
        True
        sage: is_SL2Z(Gamma0(6))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38035, "The function is_SL2Z is deprecated; use 'isinstance(..., SL2Z_class)' instead.")
    return isinstance(x, SL2Z_class)


class SL2Z_class(Gamma0_class):
    r"""
    The full modular group `\SL_2(\ZZ)`, regarded as a congruence
    subgroup of itself.
    """

    def __init__(self):
        r"""
        The modular group `\SL_2(\Z)`.

        EXAMPLES::

            sage: G = SL2Z; G
            Modular Group SL(2,Z)
            sage: G.gens()
            (
            [ 0 -1]  [1 1]
            [ 1  0], [0 1]
            )
            sage: G.0
            [ 0 -1]
            [ 1  0]
            sage: G.1
            [1 1]
            [0 1]
            sage: latex(G)
            \mbox{\rm SL}_2(\Bold{Z})
            sage: G([1,-1,0,1])
            [ 1 -1]
            [ 0  1]
            sage: TestSuite(G).run()
            sage: SL2Z.0 * SL2Z.1
            [ 0 -1]
            [ 1  1]
            sage: SL2Z is loads(dumps(SL2Z))
            True
        """
        Gamma0_class.__init__(self, 1)

    def __reduce__(self):
        """
        Used for pickling ``self``.

        EXAMPLES::

            sage: SL2Z.__reduce__()
            (<function _SL2Z_ref at ...>, ())
        """
        return _SL2Z_ref, ()

    def _element_constructor_(self, x, check=True):
        r"""
        Create an element of ``self`` from x, which must be something that can be
        coerced into a 2x2 integer matrix. If ``check=True`` (the default), check
        that x really has determinant 1.

        EXAMPLES::

            sage: SL2Z([1,0,0,1]) # indirect doctest
            [1 0]
            [0 1]
            sage: SL2Z([2, 0, 0, 2], check=False) # don't do this!
            [2 0]
            [0 2]
            sage: SL2Z([1, QQ, False])
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce <class 'sage.rings.rational_field.RationalField_with_category'> to an integer
        """
        return ArithmeticSubgroupElement(self, x, check=check)

    def _contains_sl2(self,a,b,c,d):
        r"""
        Test whether [a,b,c,d] is an element of self, where a,b,c,d are integers with `ad-bc=1`. In other words, always return True.

        EXAMPLES::

            sage: [8,7,9,8] in SL2Z # indirect doctest
            True
        """
        return True

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: SL2Z._repr_()
            'Modular Group SL(2,Z)'
        """
        return "Modular Group SL(2,Z)"

    def _latex_(self):
        r"""
        Return the \LaTeX representation of ``self``.

        EXAMPLES::

            sage: SL2Z._latex_()
            '\\mbox{\\rm SL}_2(\\Bold{Z})'
            sage: latex(SL2Z)
            \mbox{\rm SL}_2(\Bold{Z})
        """
        return "\\mbox{\\rm SL}_2(%s)" % (ZZ._latex_())

    def is_subgroup(self, right):
        """
        Return ``True`` if ``self`` is a subgroup of ``right``.

        EXAMPLES::

            sage: SL2Z.is_subgroup(SL2Z)
            True
            sage: SL2Z.is_subgroup(Gamma1(1))
            True
            sage: SL2Z.is_subgroup(Gamma0(6))
            False
        """
        return right.level() == 1

    def reduce_cusp(self, c):
        r"""
        Return the unique reduced cusp equivalent to c under the
        action of ``self``. Always returns Infinity, since there is only
        one equivalence class of cusps for `SL_2(Z)`.

        EXAMPLES::

            sage: SL2Z.reduce_cusp(Cusps(-1/4))
            Infinity
        """
        return Cusp(1,0)

    def random_element(self, bound=100, *args, **kwds):
        r"""
        Return a random element of `\SL_2(\ZZ)` with entries whose
        absolute value is strictly less than bound (default: 100).
        Additional arguments and keywords are passed to the random_element
        method of ZZ.

        (Algorithm: Generate a random pair of integers at most bound. If they
        are not coprime, throw them away and start again. If they are, find an
        element of `\SL_2(\ZZ)` whose bottom row is that, and
        left-multiply it by `\begin{pmatrix} 1 & w \\ 0 & 1\end{pmatrix}` for
        an integer `w` randomly chosen from a small enough range that the
        answer still has entries at most bound.)

        It is, unfortunately, not true that all elements of SL2Z with entries <
        bound appear with equal probability; those with larger bottom rows are
        favoured, because there are fewer valid possibilities for w.

        EXAMPLES::

            sage: s = SL2Z.random_element()
            sage: s.parent() is SL2Z
            True
            sage: all(a in range(-99, 100) for a in list(s))
            True
            sage: S = set()
            sage: while len(S) < 180:
            ....:     s = SL2Z.random_element(5)
            ....:     assert all(a in range(-4, 5) for a in list(s))
            ....:     assert s.parent() is SL2Z
            ....:     assert s in SL2Z
            ....:     S.add(s)

        Passes extra positional or keyword arguments through::

            sage: SL2Z.random_element(5, distribution='1/n').parent() is SL2Z
            True
        """
        if bound <= 1:
            raise ValueError("bound must be greater than 1")
        c = ZZ.random_element(1-bound, bound, *args, **kwds)
        d = ZZ.random_element(1-bound, bound, *args, **kwds)
        if gcd(c,d) != 1: # try again
            return self.random_element(bound, *args, **kwds)
        else:
            a,b,c,d = lift_to_sl2z(c,d,0)
            whi = bound
            wlo = bound
            if c > 0:
                whi = min(whi, ((bound - a)/ZZ(c)).ceil())
                wlo = min(wlo, ((bound + a)/ZZ(c)).ceil())
            elif c < 0:
                whi = min(whi, ((bound + a)/ZZ(-c)).ceil())
                wlo = min(wlo, ((bound - a)/ZZ(-c)).ceil())

            if d > 0:
                whi = min(whi, ((bound - b)/ZZ(d)).ceil())
                wlo = min(wlo, ((bound + b)/ZZ(d)).ceil())
            elif d < 0:
                whi = min(whi, ((bound + b)/ZZ(-d)).ceil())
                wlo = min(wlo, ((bound - b)/ZZ(-d)).ceil())

            w = ZZ.random_element(1-wlo, whi, *args, **kwds)
            a += c*w
            b += d*w
            return self([a,b,c,d])


SL2Z = SL2Z_class()


def _SL2Z_ref():
    """
    Return SL2Z. (Used for pickling SL2Z.).

    EXAMPLES::

        sage: sage.modular.arithgroup.congroup_sl2z._SL2Z_ref()
        Modular Group SL(2,Z)
        sage: sage.modular.arithgroup.congroup_sl2z._SL2Z_ref() is SL2Z
        True
    """
    return SL2Z
