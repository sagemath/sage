# sage_setup: distribution = sagemath-combinat
"""
Abelian Monoid Elements

AUTHORS:

- David Kohel (2005-09)

EXAMPLES:

Recall the example from abelian monoids::

    sage: F = FreeAbelianMonoid(5,names = list("abcde"))
    sage: (a,b,c,d,e) = F.gens()
    sage: a*b^2*e*d
    a*b^2*d*e
    sage: x = b^2*e*d*a^7
    sage: x
    a^7*b^2*d*e
    sage: x.list()
    [7, 2, 0, 1, 1]

The list is a copy, so changing the list does not change the element::

    sage: x.list()[0] = 0
    sage: x
    a^7*b^2*d*e
"""

# ***************************************************************************
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_on, sig_off
from sage.structure.richcmp cimport rich_to_bool
from sage.rings.integer cimport Integer, _Integer_from_mpz
from sage.libs.gmp.mpz cimport *


def is_FreeAbelianMonoidElement(x):
    r"""
    Queries whether ``x`` is an object of type ``FreeAbelianMonoidElement``.

    INPUT:

    - ``x`` -- an object

    OUTPUT:

    - ``True`` if ``x`` is an object of type ``FreeAbelianMonoidElement``;
      ``False`` otherwise.
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38184,
                       "The function is_FreeAbelianMonoidElement is deprecated; "
                       "use 'isinstance(..., FreeAbelianMonoidElement)' instead.")
    return isinstance(x, FreeAbelianMonoidElement)


cdef class FreeAbelianMonoidElement(MonoidElement):
    cdef int _init(self, Py_ssize_t n, Parent parent) except -1:
        """
        Initialize the C data structures in this vector. After calling
        this, ``self`` is the identity element represented as a zero vector
        of degree ``n`` with parent ``parent``.

        Only if you call ``__new__`` without a ``parent`` argument, it
        is needed to call this function manually. The only reason to do
        that is for efficiency: calling ``__new__`` without any
        additional arguments besides the type and then calling ``_init``
        is (slightly) more efficient than calling ``__new__`` with a
        ``parent`` argument.
        """
        # Assign variables only when the array is fully initialized
        cdef mpz_t* entries = <mpz_t*>check_allocarray(n, sizeof(mpz_t))
        cdef Py_ssize_t i
        sig_on()
        for i in range(n):
            mpz_init(entries[i])
        sig_off()
        self._element_vector = entries
        self._n = n
        self._parent = parent

    def __cinit__(self, parent=None, x=None):
        """
        C level initialization of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F(1)
            1
        """
        self._element_vector = NULL
        if parent is None:
            self._n = 0
            return
        self._init(parent.ngens(), <Parent?>parent)

    def __init__(self, parent, x):
        r"""
        Create the element ``x`` of the FreeAbelianMonoid ``parent``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F
            Free abelian monoid on 5 generators (a, b, c, d, e)
            sage: F(1)
            1
            sage: F(2)
            Traceback (most recent call last):
            ...
            TypeError: argument x (= 2) is of the wrong type
            sage: F(int(1))
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
        cdef Py_ssize_t i
        cdef Integer z

        if isinstance(x, int):
            if (<int> x) == 1:
                return  # the _element_vector is already set to the 0 vector
        elif isinstance(x, Integer):
            if mpz_cmp_si((<Integer> x).value, 1) == 0:
                return  # the _element_vector is already set to the 0 vector
        elif isinstance(x, (list, tuple)):
            if len(x) != self._n:
                raise IndexError("argument length (= %s) must be %s" % (len(x), self._n))
            for i in range(self._n):
                z = Integer(x[i])
                mpz_set(self._element_vector[i], z.value)
            return
        raise TypeError("argument x (= %s) is of the wrong type" % x)

    def __dealloc__(self):
        """
        Deallocate ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: x = a^2 * b^3
            sage: del x
        """
        cdef Py_ssize_t i
        if self._element_vector:
            for i in range(self._n):
                mpz_clear(self._element_vector[i])
            sig_free(self._element_vector)

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: x = a^2 * b^3
            sage: copy(x) == x
            True
            sage: copy(x) is x
            False
        """
        cdef FreeAbelianMonoidElement y
        y = self._new_c()
        cdef Py_ssize_t i
        for i in range(self._n):
            mpz_set(y._element_vector[i], self._element_vector[i])
        return y

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: x = a^2 * b^3
            sage: loads(dumps(x)) == x
            True
        """
        return (self._parent, (self.list(),))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F(1)
            1
            sage: a, b, c, d, e = F.gens()
            sage: a^2 * b^3 * a^2 * b^4
            a^4*b^7
        """
        cdef str s = ""
        A = self._parent
        cdef tuple x = A.variable_names()
        cdef mpz_t *v = self._element_vector
        cdef Integer val
        for i in range(self._n):
            val = _Integer_from_mpz(v[i])
            if val == 0:
                continue
            elif val == 1:
                if s:
                    s += "*"
                s += "%s" % x[i]
            else:
                if s:
                    s += "*"
                s += "%s^%s" % (x[i], val)
        if not s:
            s = "1"
        return s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(2, 'alpha,beta')
            sage: latex(F(1))
            1
            sage: a, b = F.gens()
            sage: latex(a^2 * b^3 * a^2 * b^4)
            \alpha^{4} \beta^{7}
        """
        cdef str s = ""
        A = self._parent
        cdef list x = A.latex_variable_names()
        cdef mpz_t *v = self._element_vector
        cdef Integer val
        for i in range(self._n):
            val = _Integer_from_mpz(v[i])
            if val == 0:
                continue
            elif val == 1:
                if s:
                    s += " "
                s += "%s" % x[i]
            else:
                if s:
                    s += " "
                s += "%s^{%s}" % (x[i], val)
        if not s:
            s = "1"
        return s

    cpdef _richcmp_(left, right, int op):
        """
        Rich comparison.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F(1)
            1
            sage: a, b, c, d, e = F.gens()
            sage: x = a^2 * b^3
            sage: F(1) < x
            True
            sage: x > b
            True
            sage: x <= a^4
            True
            sage: x != a*b
            True
            sage: a*b == b*a
            True
            sage: x > a^3*b^2
            False
        """
        cdef Py_ssize_t i
        cdef int c
        for i in range(left._n):
            c = mpz_cmp(left._element_vector[i], (<FreeAbelianMonoidElement>right)._element_vector[i])
            if c < 0:
                return rich_to_bool(op, -1)
            elif c > 0:
                return rich_to_bool(op, 1)
        return rich_to_bool(op, 0)

    def __mul__(self, y):
        """
        Multiply ``self`` with ``y``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: b * a^2 * b^3
            a^2*b^4
        """
        if not isinstance(y, FreeAbelianMonoidElement):
            raise TypeError("argument y (= %s) is of wrong type" % y)
        cdef FreeAbelianMonoidElement s, z, r
        s = self
        r = y
        z = s._new_c()
        cdef Py_ssize_t i
        for i in range(s._n):
            mpz_add(z._element_vector[i], s._element_vector[i], r._element_vector[i])
        return z

    def __pow__(self, n, modulus):
        """
        Raise ``self`` to the power of ``n``.

        AUTHORS:

        - Tom Boothby (2007-08): Replaced O(log n) time, O(n) space
          algorithm with O(1) time and space "algorithm".

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
        if modulus is not None:
            raise NotImplementedError("modulus for exponents not implemented")
        cdef Integer val
        if isinstance(n, int):
            val = Integer(n)
        else:
            if not isinstance(n, Integer):
                raise TypeError("argument n (= %s) must be an integer" % (n,))
            val = <Integer> n
        if val < 0:
            raise IndexError("argument n (= %s) must be positive" % val)
        elif val == 1:
            return self
        cdef FreeAbelianMonoidElement s, z
        s = self
        z = s._new_c()
        cdef Py_ssize_t i
        for i in range(s._n):
            mpz_mul(z._element_vector[i], s._element_vector[i], val.value)
        return z

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5,names = list("abcde"))
            sage: (a,b,c,d,e) = F.gens()
            sage: x = a*b^2*e*d
            sage: hash(x) == hash(x)
            True
        """
        return hash(tuple(self.list()))

    def list(self):
        r"""
        Return the underlying list used to represent ``self``.

        If this is a monoid in an abelian monoid on `n` generators,
        then this is a list of nonnegative integers of length `n`.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: (a, b, c, d, e) = F.gens()
            sage: a.list()
            [1, 0, 0, 0, 0]
        """
        cdef Py_ssize_t i
        return [_Integer_from_mpz(self._element_vector[i]) for i in range(self._n)]
