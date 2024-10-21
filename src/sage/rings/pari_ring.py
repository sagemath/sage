# sage_setup: distribution = sagemath-pari
"""
Ring of pari objects

AUTHORS:

- William Stein (2004): Initial version.
- Simon King (2011-08-24): Use UniqueRepresentation, element_class and
  proper initialisation of elements.
"""
# ****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import sage.libs.pari.all as pari
from sage.categories.rings import Rings
from sage.structure.parent import Parent
from sage.structure.element import RingElement
from sage.structure.richcmp import richcmp
from sage.misc.fast_methods import Singleton


class Pari(RingElement):
    """
    Element of Pari pseudo-ring.
    """
    def __init__(self, x, parent=None) -> None:
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: f = R('x^3 + 1/2')
            sage: f
            x^3 + 1/2
            sage: type(f)
            <class 'sage.rings.pari_ring.PariRing_with_category.element_class'>
            sage: loads(f.dumps()) == f
            True
        """
        if parent is None:
            parent = _inst
        RingElement.__init__(self, parent)
        self.__x = pari.pari(x)

    def __repr__(self) -> str:
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3); a
            3
        """
        return str(self.__x)

    def _add_(self, other):
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: b = R(11)
            sage: a = R(3)
            sage: a + b
            14
        """
        return self.__class__(self.__x + other.__x, parent=_inst)

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3)
            sage: b = R(11)
            sage: b - a
            8
        """
        return self.__class__(self.__x - other.__x, parent=_inst)

    def _mul_(self, other):
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3)
            sage: b = R(11)
            sage: b * a
            33
        """
        return self.__class__(self.__x * other.__x, parent=_inst)

    def _div_(self, other):
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3)
            sage: b = R(11)
            sage: b / a
            11/3
        """
        return self.__x * (~other.__x)

    def __neg__(self):
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3)
            sage: -a
            -3
        """
        return self.__class__(-self.__x, parent=_inst)

    def __pow__(self, other):
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3)
            sage: a^2
            9
        """
        if other not in PariRing():
            other = Pari(other)
        return self.__class__(self.__x ** other.__x, parent=_inst)

    def __invert__(self):
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3)
            sage: ~a
            1/3
        """
        return self.__class__(~self.__x, parent=_inst)

    def _richcmp_(self, other, op) -> bool:
        """
        EXAMPLES::

            sage: R = PariRing()
            sage: a = R(3)
            sage: b = R(11)
            sage: a < b
            True
            sage: a == b
            False
            sage: a > b
            False
        """
        return richcmp(self.__x, other.__x, op)

    def __int__(self) -> int:
        return int(self.__x)


class PariRing(Singleton, Parent):
    """
    EXAMPLES::

        sage: R = PariRing(); R
        Pseudoring of all PARI objects.
        sage: loads(R.dumps()) is R
        True
    """
    Element = Pari

    def __init__(self):
        Parent.__init__(self, self, category=Rings())

    def __repr__(self) -> str:
        return 'Pseudoring of all PARI objects.'

    def _element_constructor_(self, x):
        if isinstance(x, Pari):
            return x
        return self.element_class(x, parent=self)

    def is_field(self, proof=True) -> bool:
        return False

    def characteristic(self):
        raise RuntimeError("Not defined.")

    def random_element(self, x=None, y=None, distribution=None):
        """
        Return a random integer in Pari.

        .. NOTE::

            The given arguments are passed to ``ZZ.random_element(...)``.

        INPUT:

        - `x`, `y` -- optional integers, that are lower and upper bound
          for the result. If only `x` is provided, then the result is
          between 0 and `x-1`, inclusive. If both are provided, then the
          result is between `x` and `y-1`, inclusive.

        - ``distribution`` -- (optional) string, so that ``ZZ`` can make sense
          of it as a probability distribution

        EXAMPLES::

            sage: R = PariRing()
            sage: R.random_element().parent() is R
            True
            sage: R(5) <= R.random_element(5,13) < R(13)
            True
            sage: R.random_element(distribution='1/n').parent() is R
            True
        """
        from sage.rings.integer_ring import ZZ
        return self(ZZ.random_element(x, y, distribution))

    def zeta(self):
        """
        Return -1.

        EXAMPLES::

            sage: R = PariRing()
            sage: R.zeta()
            -1
        """
        return self(-1)


_inst = PariRing()
