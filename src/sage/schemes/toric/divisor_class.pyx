# sage_setup: distribution = sagemath-polyhedra
# sage.doctest: needs sage.geometry.polyhedron sage.graphs
r"""
Toric rational divisor classes

This module is a part of the framework for :mod:`toric varieties
<sage.schemes.toric.variety>`.

AUTHORS:

- Volker Braun and Andrey Novoseltsev (2010-09-05): initial version.

TESTS:

Toric rational divisor classes are elements of the rational class group of a
toric variety, represented as rational vectors in some basis::

    sage: dP6 = toric_varieties.dP6()
    sage: Cl = dP6.rational_class_group()
    sage: D = Cl([1, -2, 3, -4])
    sage: D
    Divisor class [1, -2, 3, -4]
    sage: E = Cl([1/2, -2/3, 3/4, -4/5])
    sage: E
    Divisor class [1/2, -2/3, 3/4, -4/5]

They behave much like ordinary vectors::

    sage: D + E
    Divisor class [3/2, -8/3, 15/4, -24/5]
    sage: 2 * D
    Divisor class [2, -4, 6, -8]
    sage: E / 10
    Divisor class [1/20, -1/15, 3/40, -2/25]
    sage: D * E
    Traceback (most recent call last):
    ...
    TypeError: cannot multiply two divisor classes

The only special method is :meth:`~ToricRationalDivisorClass.lift` to get a
divisor representing a divisor class::

    sage: D.lift()
    -3*V(x) - 9*V(u) + 7*V(z) + 3*V(w)
    sage: E.lift()
    -3/10*V(x) - 133/60*V(u) + 31/20*V(z) + 3/4*V(w)
"""


#*****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.mpq cimport *

from sage.misc.latex import latex
from sage.modules.free_module_element import vector
from sage.modules.vector_rational_dense cimport Vector_rational_dense
from sage.rings.rational_field import QQ
from sage.rings.rational cimport Rational
from sage.structure.element cimport Vector
from sage.structure.element import Vector


def is_ToricRationalDivisorClass(x):
    r"""
    Check if ``x`` is a toric rational divisor class.

    INPUT:

    - ``x`` -- anything

    OUTPUT: ``True`` if ``x`` is a toric rational divisor class, ``False`` otherwise

    EXAMPLES::

        sage: from sage.schemes.toric.divisor_class import is_ToricRationalDivisorClass
        sage: is_ToricRationalDivisorClass(1)
        False
        sage: dP6 = toric_varieties.dP6()
        sage: D = dP6.rational_class_group().gen(0); D
        Divisor class [1, 0, 0, 0]
        sage: is_ToricRationalDivisorClass(D)
        True
    """
    return isinstance(x, ToricRationalDivisorClass)


cdef class ToricRationalDivisorClass(Vector_rational_dense):
    r"""
    Create a toric rational divisor class.

    .. WARNING::

        You probably should not construct divisor classes explicitly.

    INPUT:

    - same as for
      :class:`~sage.modules.vector_rational_dense.Vector_rational_dense`.

    OUTPUT: toric rational divisor class

    TESTS::

        sage: dP6 = toric_varieties.dP6()
        sage: Cl = dP6.rational_class_group()
        sage: D = dP6.divisor(2)
        sage: Cl(D)
        Divisor class [0, 0, 1, 0]
    """

    def __reduce__(self):
        """
        Prepare ``self`` for pickling.

        TESTS::

            sage: dP6 = toric_varieties.dP6()
            sage: Cl = dP6.rational_class_group()
            sage: D = Cl([1, -2, 3, -4])
            sage: D
            Divisor class [1, -2, 3, -4]
            sage: loads(dumps(D))
            Divisor class [1, -2, 3, -4]
        """
        return (_ToricRationalDivisorClass_unpickle_v1,
                (self._parent, list(self), self._degree,
                 not self._is_immutable))

    cpdef _act_on_(self, other, bint self_on_left):
        """
        Act on ``other``.

        INPUT:

        - ``other`` -- something that
          :class:`~sage.modules.vector_rational_dense.Vector_rational_dense`
          can act on *except* for another toric rational divisor class.

        OUTPUT:

        - standard output for ``self`` acting as a rational vector on
          ``other`` if the latter one is not a toric rational divisor class.

        TESTS::

            sage: dP6 = toric_varieties.dP6()
            sage: Cl = dP6.rational_class_group()
            sage: D = Cl([1, -2, 3, -4])
            sage: D
            Divisor class [1, -2, 3, -4]
            sage: D * D
            Traceback (most recent call last):
            ...
            TypeError: cannot multiply two divisor classes

        We test standard behaviour::

            sage: v = vector([1, 2, 3, 4])
            sage: v * D # indirect doctest
            -10
            sage: D * v # indirect doctest
            -10
            sage: v = vector([1, 2/3, 4/5, 6/7])
            sage: v * D # indirect doctest
            -143/105
            sage: D * v # indirect doctest
            -143/105
            sage: A = matrix(4, range(16))
            sage: A * D # indirect doctest
            (-8, -16, -24, -32)
            sage: D * A # indirect doctest
            (-32, -34, -36, -38)
            sage: B = A / 3
            sage: B * D # indirect doctest
            (-8/3, -16/3, -8, -32/3)
            sage: D * B # indirect doctest
            (-32/3, -34/3, -12, -38/3)
        """
        # If we don't treat vectors separately, they get converted into
        # divisor classes where multiplication is prohibited on purpose.
        if isinstance(other, Vector_rational_dense):
            return Vector_rational_dense._dot_product_(self, other)
        cdef Vector v
        if isinstance(other, Vector) and not is_ToricRationalDivisorClass(other):
            try:
                v = vector(QQ, other)
                if v._degree == self._degree:
                    return Vector_rational_dense._dot_product_(self, v)
            except TypeError:
                pass
        # Now let the standard framework work...
        return Vector_rational_dense._act_on_(self, other, self_on_left)

    cpdef _dot_product_(self, Vector right):
        r"""
        Raise a :exc:`TypeError` exception.

        Dot product is not defined on toric rational divisor classes.

        INPUT:

        - ``right`` -- vector

        OUTPUT: a :exc:`TypeError` exception is raised

        TESTS::

            sage: c = toric_varieties.dP8().rational_class_group().gens()
            sage: c[0]._dot_product_(c[1])
            Traceback (most recent call last):
            ...
            TypeError: cannot multiply two divisor classes
            sage: c[0] * c[1]      # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: cannot multiply two divisor classes
        """
        raise TypeError("cannot multiply two divisor classes")

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: D = toric_varieties.dP6().divisor(0).divisor_class()
            sage: print(D._latex_())
            \left[ 1, 0, 0, 0 \right]_{\mathop{Cl}_{\QQ}\left(\mathbb{P}_{\Delta^{2}_{9}}\right)}
        """
        return r"\left[ %s \right]_{%s}" % (
                    ", ".join(latex(e) for e in self), latex(self.parent()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: toric_varieties.dP6().divisor(0).divisor_class()._repr_()
            'Divisor class [1, 0, 0, 0]'
        """
        return 'Divisor class %s' % list(self)

    def lift(self):
        r"""
        Return a divisor representing this divisor class.

        OUTPUT: an instance of :class:`ToricDivisor` representing ``self``

        EXAMPLES::

            sage: X = toric_varieties.Cube_nonpolyhedral()
            sage: D = X.divisor([0,1,2,3,4,5,6,7]); D
            V(z1) + 2*V(z2) + 3*V(z3) + 4*V(z4) + 5*V(z5) + 6*V(z6) + 7*V(z7)
            sage: D.divisor_class()
            Divisor class [29, 6, 8, 10, 0]
            sage: Dequiv = D.divisor_class().lift(); Dequiv
            15*V(z1) - 11*V(z2) - 9*V(z5) + 19*V(z6) + 10*V(z7)
            sage: Dequiv == D
            False
            sage: Dequiv.divisor_class() == D.divisor_class()
            True
        """
        Cl = self.parent()
        return Cl._variety.divisor(Cl._lift_matrix * self)


def _ToricRationalDivisorClass_unpickle_v1(parent, entries,
                                           degree, is_mutable):
    """
    Unpickle a :class:`toric rational divisor class
    <ToricRationalDivisorClass>`.

    INPUT:

    - ``parent`` -- rational divisor class group of a toric variety

    - ``entries`` -- list of rationals specifying the divisor class

    - ``degree`` -- integer; dimension of the ``parent``

    - ``is_mutable`` -- boolean, whether the divisor class is mutable

    OUTPUT: :class:`toric rational divisor class <ToricRationalDivisorClass>`

    TESTS::

        sage: dP6 = toric_varieties.dP6()
        sage: Cl = dP6.rational_class_group()
        sage: D = Cl([1, -2, 3, -4])
        sage: D
        Divisor class [1, -2, 3, -4]
        sage: loads(dumps(D))   # indirect test
        Divisor class [1, -2, 3, -4]
        sage: from sage.schemes.toric.divisor_class import (
        ....:     _ToricRationalDivisorClass_unpickle_v1)
        sage: _ToricRationalDivisorClass_unpickle_v1(
        ....:    Cl, [1, -2, 3, -4], 4, True)
        Divisor class [1, -2, 3, -4]
    """
    cdef ToricRationalDivisorClass v
    v = ToricRationalDivisorClass.__new__(ToricRationalDivisorClass)
    v._init(degree, parent)
    cdef Rational z
    cdef Py_ssize_t i
    for i in range(degree):
        z = Rational(entries[i])
        mpq_set(v._entries[i], z.value)
    v._is_immutable = not is_mutable
    return v
