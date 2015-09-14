"""
Abelian Monoid Elements

AUTHORS:

- David Kohel (2005-09)

EXAMPLES:

Recall the example from abelian monoids.

::

    sage: F = FreeAbelianMonoid(5,names = list("abcde"))
    sage: (a,b,c,d,e) = F.gens()
    sage: a*b^2*e*d
    a*b^2*d*e
    sage: x = b^2*e*d*a^7
    sage: x
    a^7*b^2*d*e
    sage: x.list()
    [7, 2, 0, 1, 1]

It is important to note that lists are mutable and the returned
list is not a copy. As a result, reassignment of an element of the
list changes the object.

::

    sage: x.list()[0] = 0
    sage: x
    b^2*d*e
"""

#*****************************************************************************
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.integer import Integer
from sage.structure.element import MonoidElement

def is_FreeAbelianMonoidElement(x):
    r"""
    Queries whether ``x`` is an object of type ``FreeAbelianMonoidElement``.

    INPUT:

    - ``x`` -- an object.

    OUTPUT:

    - ``True`` if ``x`` is an object of type ``FreeAbelianMonoidElement``;
      ``False`` otherwise.
    """
    return isinstance(x, FreeAbelianMonoidElement)

class FreeAbelianMonoidElement(MonoidElement):
    def __init__(self, F, x):
        """
        Create the element x of the FreeAbelianMonoid F.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F
            Free abelian monoid on 5 generators (a, b, c, d, e)
            sage: F(1)
            1
            sage: a, b, c, d, e = F.gens()
            sage: a^2 * b^3 * a^2 * b^4
            a^4*b^7
            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: a in F
            True
            sage: a*b in F
            True
        """
        MonoidElement.__init__(self, F)
        self.__repr = None
        n = F.ngens()
        if isinstance(x, (int, long, Integer)) and x == 1:
            self.__element_vector = [ 0 for i in range(n) ]
        elif isinstance(x, list):
            if len(x) != n:
                raise IndexError("Argument length (= %s) must be %s."%(len(x), n))
            self.__element_vector = x
        else:
            raise TypeError("Argument x (= %s) is of wrong type."%x)

    def __repr__(self):
        s = ""
        A = self.parent()
        n = A.ngens()
        x = A.variable_names()
        v = self.__element_vector
        for i in range(n):
            if v[i] == 0:
                continue
            elif v[i] == 1:
                if len(s) > 0: s += "*"
                s += "%s"%x[i]
            else:
                if len(s) > 0: s += "*"
                s += "%s^%s"%(x[i],v[i])
        if len(s) == 0: s = "1"
        return s

    def __mul__(self, y):
        if not isinstance(y, FreeAbelianMonoidElement):
            raise TypeError("Argument y (= %s) is of wrong type."%y)
        M = self.parent()
        z = M(1)
        xelt = self.__element_vector
        yelt = y.__element_vector
        z.__element_vector = [ xelt[i]+yelt[i] for i in range(len(xelt)) ]
        return z

    def __pow__(self, n):
        """
        Raises self to the power of n.

        AUTHORS:

        - Tom Boothby (2007-08): Replaced O(log n) time, O(n) space
          algorithm with O(1) time and space"algorithm".

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5,names = list("abcde"))
            sage: (a,b,c,d,e) = F.gens()
            sage: x = a*b^2*e*d; x
            a*b^2*d*e
            sage: x^3
            a^3*b^6*d^3*e^3
            sage: x^0
            1
        """

        if not isinstance(n, (int, long, Integer)):
            raise TypeError("Argument n (= %s) must be an integer."%n)
        if n < 0:
            raise IndexError("Argument n (= %s) must be positive."%n)
        elif n == 1:
            return self
        z = self.parent()(1)
        if n == 0:
            return z
        else:
            z.__element_vector = [i*n for i in self.__element_vector]
            return z


    def list(self):
        """
        Return (a reference to) the underlying list used to represent this
        element. If this is a monoid in an abelian monoid on `n`
        generators, then this is a list of nonnegative integers of length
        `n`.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: (a, b, c, d, e) = F.gens()
            sage: a.list()
            [1, 0, 0, 0, 0]
        """
        return self.__element_vector


