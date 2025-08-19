# sage.doctest: needs sage.modules
"""
Octonion Algebras

AUTHORS:

- Travis Scrimshaw (2023-05-06): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.functional import sqrt
from sage.structure.element cimport Element, parent
from sage.structure.parent cimport Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp cimport richcmp
from sage.structure.category_object import normalize_names
from sage.modules.free_module import FreeModule
from sage.categories.magmatic_algebras import MagmaticAlgebras
from sage.categories.rings import Rings
from sage.categories.metric_spaces import MetricSpaces


cdef class Octonion_generic(AlgebraElement):
    r"""
    An octonion with generic parameters.
    """
    def __init__(self, parent, v):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: TestSuite(elt).run()
            sage: elt = sum(O.basis())
            sage: TestSuite(elt).run()
            sage: TestSuite(O.zero()).run()

            sage: O = OctonionAlgebra(Zmod(8), 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: TestSuite(elt).run()
            sage: elt = sum(O.basis())
            sage: TestSuite(elt).run()
            sage: TestSuite(O.zero()).run()
        """
        v.set_immutable()
        self.vec = v
        super().__init__(parent)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: elt
            2 + 3*i + 4*j + 5*k + 6*l + 7*li + 8*lj + 9*lk
        """
        from sage.algebras.weyl_algebra import repr_from_monomials
        data = [p for p in enumerate(self.vec) if p[1]]
        return repr_from_monomials(data, self._parent._repr_term)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: latex(elt)
            2 + 3 i + 4 j + 5 k + 6 l + 7 li + 8 lj + 9 lk
        """
        from sage.algebras.weyl_algebra import repr_from_monomials
        data = [p for p in enumerate(self.vec) if p[1]]
        return repr_from_monomials(data, self._parent._repr_term, True)

    def __bool__(self):
        r"""
        Return if ``self`` is nonzero or not.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: bool(elt)
            True
            sage: bool(O.zero())
            False
        """
        return bool(self.vec)

    def __reduce__(self):
        r"""
        For pickling.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: x = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: x.__reduce__()
            (<class 'sage.algebras.octonion_algebra.Octonion_generic'>,
             (Octonion algebra over Rational Field with parameters (1, 3, 7),
              (2, 3, 4, 5, 6, 7, 8, 9)))
        """
        return (self.__class__, (self._parent, self.vec))

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare ``self`` to ``other`` with type ``op``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: x = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: y = sum(O.basis())
            sage: x != y
            True
            sage: x == y
            False
            sage: x < y
            False
        """
        return richcmp(self.vec, (<Octonion_generic> other).vec, op)

    def __hash__(self):
        r"""
        Return a hash of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: hash(elt) == hash(elt.vector())
            True
        """
        return hash(self.vec)

    cpdef _add_(self, other):
        r"""
        Return ``self`` plus ``other``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: x = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: y = sum(O.basis())
            sage: x + y
            3 + 4*i + 5*j + 6*k + 7*l + 8*li + 9*lj + 10*lk
        """
        return self.__class__(self._parent, self.vec + (<Octonion_generic> other).vec)

    cpdef _sub_(self, other):
        r"""
        Return ``self`` minus ``other``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: x = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: y = sum(O.basis())
            sage: x - y
            1 + 2*i + 3*j + 4*k + 5*l + 6*li + 7*lj + 8*lk
        """
        return self.__class__(self._parent, self.vec - (<Octonion_generic> other).vec)

    def __neg__(self):
        r"""
        Return the negative of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: x = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: -x
            -2 - 3*i - 4*j - 5*k - 6*l - 7*li - 8*lj - 9*lk
            sage: y = sum(O.basis())
            sage: -y
            -1 - i - j - k - l - li - lj - lk
        """
        return self.__class__(self._parent, -self.vec)

    cpdef _lmul_(self, Element other):
        r"""
        Return ``self * other`` for a scalar ``other``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: i, j, k, l = O.gens()
            sage: elt = 2 * i + 3 * j + k * 4 + l * 5
            sage: elt * 5
            10*i + 15*j + 20*k + 25*l
        """
        return self.__class__(self._parent, self.vec * other)

    cpdef _rmul_(self, Element other):
        r"""
        Return ``self * other`` for a scalar ``other``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: i, j, k, l = O.gens()
            sage: elt = 2 * i + 3 * j + k * 4 + l * 5
            sage: 5 * elt
            10*i + 15*j + 20*k + 25*l
        """
        return self.__class__(self._parent, other * self.vec)

    cpdef _mul_(self, other):
        r"""
        Return ``self`` multiplied by ``other``.

        INPUT:

        - ``other`` -- element of the octonion algebra as ``self``

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: table([[r * c for c in O.basis()] for r in O.basis()])
              1    i     j     k     l     li    lj    lk
              i    -1    -k    j     -li   l     lk    -lj
              j    k     -1    -i    -lj   -lk   l     li
              k    -j    i     -1    -lk   lj    -li   l
              l    li    lj    lk    -1    -i    -j    -k
              li   -l    lk    -lj   i     -1    k     -j
              lj   -lk   -l    li    j     -k    -1    i
              lk   lj    -li   -l    k     j     -i    -1

            sage: SO = OctonionAlgebra(QQ, c=1)
            sage: table([[r * c for c in SO.basis()] for r in SO.basis()])
              1    i     j     k     l     li    lj    lk
              i    -1    -k    j     -li   l     lk    -lj
              j    k     -1    -i    -lj   -lk   l     li
              k    -j    i     -1    -lk   lj    -li   l
              l    li    lj    lk    1     i     j     k
              li   -l    lk    -lj   -i    1     -k    j
              lj   -lk   -l    li    -j    k     1     -i
              lk   lj    -li   -l    -k    -j    i     1
        """
        cdef Parent P = self._parent
        cdef list ret = [P._base.zero()] * 8
        cdef tuple table = <tuple> P._mult_table
        cdef tuple row
        cdef int i, j, k
        cdef FreeModuleElement rv = (<Octonion_generic> other).vec
        for i in range(8):
            row = <tuple> table[i]
            cl = self.vec.get_unsafe(i)
            for j in range(8):
                k = <int> row[j][0]
                coeff = row[j][1]
                cr = rv.get_unsafe(j)
                ret[k] += cl * cr * coeff
        return self.__class__(P, P._module(ret))

    cpdef _div_(self, other):
        """
        Return ``self`` divided by ``other``.

        EXAMPLES::

            sage: O = OctonionAlgebra(ZZ)
            sage: x = O([0, 0, 2, 0, 2, 0, 0, 0])
            sage: y = O([1, 0, 0, 0, 1, 0, 0, 0])
            sage: y.quadratic_form()
            2
            sage: x * y.conjugate()
            2 + 2*j + 2*l + 2*lj
            sage: x / y
            1 + j + l + lj
        """
        cdef Octonion_generic temp = <Octonion_generic> self * (<Octonion_generic> other).conjugate()
        qf = (<Octonion_generic> other).quadratic_form()
        return self.__class__(self._parent, temp.vec / qf)

    def is_unit(self):
        r"""
        Return if ``self`` is a unit or not.

        EXAMPLES::

            sage: O = OctonionAlgebra(ZZ)
            sage: x = O([1, 0, 1, 0, 0, 0, 0, 0])
            sage: x.quadratic_form()
            2
            sage: x.is_unit()
            False
            sage: O([1, 0, -1, 0, 0, 0, 0, 0]).is_unit()
            False
            sage: x = O([1, 0, 0, 0, 1, 0, 0, 0])
            sage: x.quadratic_form()
            2
            sage: x.is_unit()
            False
            sage: x = O([0, 0, 0, 0, 1, 0, 0, 0])
            sage: x.quadratic_form()
            1
            sage: x.is_unit()
            True

            sage: O = OctonionAlgebra(ZZ, -1, 1, 2)
            sage: x = O([1, 0, 1, 0, 0, 0, 0, 0])
            sage: x.quadratic_form()
            0
            sage: x.is_unit()
            False
            sage: O([1, 0, -1, 0, 0, 0, 0, 0]).is_unit()
            False
            sage: x = O([1, 0, 0, 0, 1, 0, 0, 0])
            sage: x.quadratic_form()
            -1
            sage: x.is_unit()
            True
            sage: x = O([0, 0, 0, 0, 1, 0, 0, 0])
            sage: x.quadratic_form()
            -2
            sage: x.is_unit()
            False
        """
        return self.quadratic_form().is_unit()

    def __invert__(self):
        r"""
        Return the (multiplicative) inverse of ``self``.

        EXAMPLES::

            sage: O137 = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(O137.basis())
            sage: elt.quadratic_form()
            0
            sage: ~elt
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
            sage: ~O137.zero()
            Traceback (most recent call last):
            ...
            ZeroDivisionError
        """
        if not self.vec:
            raise ZeroDivisionError
        return self.quadratic_form().inverse_of_unit() * self.conjugate()

    cpdef Octonion_generic conjugate(self):
        r"""
        Return the conjugate of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: elt = sum(O.basis()); elt
            1 + i + j + k + l + li + lj + lk
            sage: elt.conjugate()
            1 - i - j - k - l - li - lj - lk
        """
        cdef FreeModuleElement v = <FreeModuleElement> -self.vec
        v.set_unsafe(0, -v.get_unsafe(0))
        return self.__class__(self._parent, v)

    cpdef quadratic_form(self):
        r"""
        Return the quadratic form of ``self``.

        The octonion algebra has a distinguished quadratic form given
        by `N(x) = x x^*`, where `x^*` is the :meth:`conjugate` of `x`.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(O.basis())
            sage: elt.quadratic_form()
            0
            sage: elt * elt.conjugate()
            0
        """
        cdef int i
        cdef tuple table = self._parent._mult_table
        ret = self.vec.get_unsafe(0) ** 2
        for i in range(1, 8):
            ret += -(<tuple> table[i])[i][1] * self.vec.get_unsafe(i) ** 2
        return ret

    cpdef norm(self):
        r"""
        Return the norm of ``self``.

        The norm of an octonion `x` is `\lVert x \rVert = \sqrt{x x^*}`,
        where `x^*` is the :meth:`conjugate` of `x`.

        .. SEEALSO::

            This is the square root of :meth:`quadratic_form()`.

        .. WARNING::

            If any of the parameters `a, b, c \not> 0`, then this is not
            an actual norm.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: elt.norm()                                                            # needs sage.symbolic
            2*sqrt(-61)
            sage: elt = sum(O.basis())
            sage: elt.norm()
            0
        """
        return sqrt(self.quadratic_form())

    cpdef abs(self):
        r"""
        Return the absolute value of ``self``.

        This is equal to the :meth:`norm`.

        .. WARNING::

            If any of the parameters `a, b, c \not> 0`, then this does
            not define a metric.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: elt.abs()                                                             # needs sage.symbolic
            2*sqrt(-61)
            sage: elt = sum(O.basis())
            sage: elt.abs()
            0
        """
        return self.norm()

    cpdef real_part(self):
        r"""
        Return the real part of ``self``.

        OUTPUT: the real part of ``self`` as an element in the base ring

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2)); elt
            2 + 3*i + 4*j + 5*k + 6*l + 7*li + 8*lj + 9*lk
            sage: r = elt.real_part(); r
            2
            sage: r.parent() is QQ
            True
        """
        return self.vec.get_unsafe(0)

    cpdef Octonion_generic imag_part(self):
        r"""
        Return the imginary part of ``self``.

        OUTPUT: the imaginary part of ``self`` as an element in the octonion algebra

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2)); elt
            2 + 3*i + 4*j + 5*k + 6*l + 7*li + 8*lj + 9*lk
            sage: elt.imag_part()
            3*i + 4*j + 5*k + 6*l + 7*li + 8*lj + 9*lk

        TESTS::

            sage: O = OctonionAlgebra(QQ)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: elt.imag_part()
            3*i + 4*j + 5*k + 6*l + 7*li + 8*lj + 9*lk
            sage: elt
            2 + 3*i + 4*j + 5*k + 6*l + 7*li + 8*lj + 9*lk
        """
        cdef FreeModuleElement v = <FreeModuleElement> self.vec.__copy__()
        v.set_unsafe(0, self._parent._base.zero())
        return self.__class__(self._parent, v)

    def _vector_(self, new_base_ring=None):
        r"""
        Return ``self`` as a vector in `R^8`, where `R` is the base
        ring of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: elt.vector()
            (2, 3, 4, 5, 6, 7, 8, 9)
        """
        if new_base_ring is None:
            return self.vec
        return self.vec.change_ring(new_base_ring)

    vector = _vector_

    def monomial_coefficients(self, copy=False):
        """
        Return ``self`` as a ``dict`` with keys being indices
        for the basis and the values being the corresponding
        nonzero coefficients.

        INPUT:

        - ``copy`` -- ignored

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: x = O([2/7, 0, 0, 0, 2/3, 0, -5, 0])
            sage: x.monomial_coefficients()
            {0: 2/7, 4: 2/3, 6: -5}
        """
        return self.vec.dict()

    dict = monomial_coefficients


cdef class Octonion(Octonion_generic):
    r"""
    An octonion.

    This is an element of the octonion algebra with parameters
    `a = b = c = -1`, which is a classical octonion number.
    """
    cpdef quadratic_form(self):
        r"""
        Return the quadratic form of ``self``.

        The octonion algebra has a distinguished quadratic form given
        by `N(x) = x x^*`, where `x^*` is the :meth:`conjugate` of `x`.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: elt = sum(O.basis()); elt
            1 + i + j + k + l + li + lj + lk
            sage: elt.quadratic_form()
            8
            sage: elt * elt.conjugate()
            8
        """
        return self.vec * self.vec

    cpdef norm(self):
        r"""
        Return the norm of ``self``.

        The norm of an octonion `x` is `\lVert x \rVert = \sqrt{x x^*}`,
        where `x^*` is the :meth:`conjugate` of `x`.

        .. SEEALSO::

            This is the square root of :meth:`quadratic_form()`.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: elt = sum(i * b for i, b in enumerate(O.basis(), start=2))
            sage: elt.norm()                                                            # needs sage.symbolic
            2*sqrt(71)
            sage: elt = sum(O.basis())
            sage: elt.norm()                                                            # needs sage.symbolic
            2*sqrt(2)
        """
        return self.vec.norm()


class OctonionAlgebra(UniqueRepresentation, Parent):
    r"""
    The octonion algebra.

    Let `R` be a commutative ring of characteristic not equal to `2`. The
    *octonion algebra* with parameters `a, b, c` is a non-associative
    non-commutative unital 8-dimensional `R`-algebra that is a deformation
    of the usual octonions, which are when `a = b = c = -1`. The octonions
    were originally constructed by Graves and independently discovered by
    Cayley (due to being published first, these are sometimes called
    the Cayley numbers) and can also be built from the Cayley-Dickson
    construction with the :class:`quaternions <QuaternionAlgebra>`.

    We use the multiplication table from [Scha1996]_. The octonion
    algebra `\mathbf{O}_{a,b,c}(R)` is a composition (Hurwitz) algebra,
    which means it is also an alternative algebra as it satisfies
    `x^2 y = (x x) y = x (x y)` and  `y x^2 = y (x x) = (y x) x`
    for all `x, y \in \mathbf{O}_{a,b,c}`.

    EXAMPLES:

    We first create the classical octonions and perform some basic
    computations::

        sage: O = OctonionAlgebra(QQ)
        sage: O
        Octonion algebra over Rational Field
        sage: i, j, k, l = O.gens()
        sage: i * j * k
        1
        sage: k * j * i
        -1
        sage: (i * k) * l
        -lj
        sage: i * (k * l)
        lj
        sage: elt = sum(O.basis())
        sage: elt^2
        -6 + 2*i + 2*j + 2*k + 2*l + 2*li + 2*lj + 2*lk
        sage: prod(O.basis())
        1
        sage: (i + l)^2
        -2
        sage: (1 + l) * (1 + l).conjugate()
        2
        sage: S = O.some_elements()
        sage: B = O.basis()
        sage: S.extend(x * (i + j/2 - 5*k/3) for x in O.some_elements())
        sage: all((x * x) * y == (x * (x * y)) for x in S for y in S)
        True
        sage: all(y * (x * x) == (y * x) * x for x in S for y in S)
        True
        sage: all((x + x.conjugate()) / 2 == x.real_part() for x in S)
        True
        sage: all((x - x.conjugate()) / 2 == x.imag_part() for x in S)
        True
        sage: all(sum((b*x)*b for b in B) == -6 * x.conjugate() for x in S)
        True

    We construct the (rescaled) `E_8` lattice as the integral octonions,
    which we verify by constructing `240` shortest length elements in
    the lattice (see also :wikipedia:`E8_lattice#Integral_octonions`)::

        sage: m = (i + j + k + l) / 2
        sage: basis = [i, j, i*j, i^2, m, i * m, j * m, (i * j) * m]
        sage: basis
        [i,
         j,
         -k,
         -1,
         1/2*i + 1/2*j + 1/2*k + 1/2*l,
         -1/2 + 1/2*j - 1/2*k - 1/2*li,
         -1/2 - 1/2*i + 1/2*k - 1/2*lj,
         1/2 - 1/2*i + 1/2*j + 1/2*lk]
        sage: matrix([vector(b) for b in basis]).rank()
        8
        sage: [b.norm() for b in basis]
        [1, 1, 1, 1, 1, 1, 1, 1]
        sage: roots = set(basis)
        sage: roots.update(-b for b in basis)
        sage: new_roots = set(roots)  # make a copy
        sage: while new_roots:
        ....:     prev_roots = new_roots
        ....:     new_roots = set()
        ....:     for a in prev_roots:
        ....:         for b in roots:
        ....:             c = a + b
        ....:             if c.quadratic_form() != 1 or c in roots:
        ....:                 continue
        ....:             new_roots.update([c, -c])
        ....:         roots.update(new_roots)
        sage: len(roots)
        240

    A classical construction of the Lie algebra of type `G_2` is
    the Lie algebra of all derivations of `\mathbf{O}` (as the
    automorphism group is the Lie group of type `G_2`). We verify
    that the derivations have the correct dimension::

        sage: len(O.derivations_basis())
        14

    We can construct the split octonions by taking the parameter `c = 1`::

        sage: SO = OctonionAlgebra(QQ, c=1)
        sage: SO
        Octonion algebra over Rational Field with parameters (-1, -1, 1)
        sage: i, j, k, l = SO.gens()
        sage: i^2 == j^2 == k^2 == -1
        True
        sage: l^2
        1
        sage: (i + l)^2
        0
        sage: (1 + l) * (1 + l).conjugate()
        0

    REFERENCES:

    - [Scha1996]_
    - :wikipedia:`octonion`
    - :wikipedia:`Split-octonion`
    - :wikipedia:`Hurwitz's_theorem_(composition_algebras)`
    """
    @staticmethod
    def __classcall_private__(cls, R, a=-1, b=-1, c=-1, names=('i', 'j', 'k', 'l')):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: O1 = OctonionAlgebra(QQ, 1)
            sage: O2 = OctonionAlgebra(QQ, QQ.one(), -1, QQ['x'](-1))
            sage: O3 = algebras.Octonion(QQ, QQ.one(), c=int(-1), names='ijkl')
            sage: O1 is O2
            True
            sage: O1 is O3
            True

            sage: OctonionAlgebra(ExteriorAlgebra(QQ, 3))
            Traceback (most recent call last):
            ...
            ValueError: the base ring must be a commutative ring
            sage: OctonionAlgebra(GF(2)['x','y','z'])
            Traceback (most recent call last):
            ...
            ValueError: the characteristic must not be 2
        """
        if R not in Rings().Commutative():
            raise ValueError("the base ring must be a commutative ring")
        if R.one() + R.one() == R.zero():
            raise ValueError("the characteristic must not be 2")
        a = R(a)
        b = R(b)
        c = R(c)
        names = normalize_names(4, names)
        return super().__classcall__(cls, R, a, b, c, names)

    def __init__(self, R, a, b, c, names):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: TestSuite(O).run()                                                    # needs sage.symbolic

            sage: O = OctonionAlgebra(QQ, 1, 3, 7)
            sage: TestSuite(O).run()

            sage: O = OctonionAlgebra(GF(3), 1, 2, 2)
            sage: TestSuite(O).run()

            sage: O = OctonionAlgebra(Zmod(6), -1, 2, 3)
            sage: TestSuite(O).run()

            sage: R.<a, b, c> = QQ[]
            sage: O = OctonionAlgebra(R, a, b, c)
            sage: TestSuite(O).run()
        """
        self._params = (a, b, c)
        self._module = FreeModule(R, 8)
        cat = MagmaticAlgebras(R.category()).Unital().WithBasis().FiniteDimensional()
        if a == b == c == -1:
            self.Element = Octonion
            if R in MetricSpaces():
                cat &= MetricSpaces()
        Parent.__init__(self, base=R, category=cat, names=names)

        # setup the multiplication table
        d = a * b
        e = a * c
        f = b * c
        g = a * b * c
        self._mult_table = (
          ((0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 1)),
          ((1, 1), (0, a), (3, -1), (2, -a), (5, -1), (4, -a), (7, 1), (6, a)),
          ((2, 1), (3, 1), (0, b), (1, b), (6, -1), (7, -1), (4, -b), (5, -b)),
          ((3, 1), (2, a), (1, -b), (0, -d), (7, -1), (6, -a), (5, b), (4, d)),
          ((4, 1), (5, 1), (6, 1), (7, 1), (0, c), (1, c), (2, c), (3, c)),
          ((5, 1), (4, a), (7, 1), (6, a), (1, -c), (0, -e), (3, -c), (2, -e)),
          ((6, 1), (7, -1), (4, b), (5, -b), (2, -c), (3, c), (0, -f), (1, f)),
          ((7, 1), (6, -a), (5, b), (4, -d), (3, -c), (2, e), (1, -f), (0, g)),
        )

    def _test_alternative(self, **options):
        r"""
        Test that ``self`` is an alternative algebra.

        An algebra `A` is *alternative* if for all `x, y \in A`, we have
        `(x * x) * y = x * (x * y)` and `(y * x) * x = y * (x * x)`.

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: O = OctonionAlgebra(R, a, b, c)
            sage: O._test_alternative()
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(S, 2, tester._max_runs):
            tester.assertEqual((x * x) * y, x * (x * y))
            tester.assertEqual(y * (x * x), (y * x) * x)

    def _test_hurwitz(self, **options):
        r"""
        Test that ``self`` is an Hurwitz algebra.

        An algebra `A` is *Hurwitz* if there exists a nondegenerate quadratic
        form `N` such that `N(x y) = N(x) N(y)` for all `x, y \in A`.

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: O = OctonionAlgebra(R, a, b, c)
            sage: O._test_hurwitz()
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(S, 2, tester._max_runs):
            tester.assertEqual((x * y).quadratic_form(), x.quadratic_form() * y.quadratic_form())

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: O = OctonionAlgebra(QQ)
            sage: [O._repr_term(i) for i in range(8)]
            ['1', 'i', 'j', 'k', 'l', 'li', 'lj', 'lk']
        """
        data = ['1', 'i', 'j', 'k', 'l', 'li', 'lj', 'lk']
        return data[m]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: OctonionAlgebra(QQ)
            Octonion algebra over Rational Field
            sage: OctonionAlgebra(QQ, 1, 3, 7)
            Octonion algebra over Rational Field with parameters (1, 3, 7)
        """
        ret = f"Octonion algebra over {self.base_ring()}"
        a, b, c = self._params
        if a != -1 or b != -1 or c != -1:
            ret += f" with parameters ({a}, {b}, {c})"
        return ret

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(OctonionAlgebra(QQ))
            \Bold{O}\left(\Bold{Q}\right)
            sage: latex(OctonionAlgebra(QQ, 1, 3, 7))
            \Bold{O}_{1, 3, 7}\left(\Bold{Q}\right)
        """
        from sage.misc.latex import latex
        ret = r"\Bold{O}"
        a, b, c = self._params
        if a != -1 or b != -1 or c != -1:
            ret += r"_{{{}, {}, {}}}".format(latex(a), latex(b), latex(c))
        return ret + r"\left({}\right)".format(latex(self.base_ring()))

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: O137 = OctonionAlgebra(QQ, 1, 3, 7)
            sage: x = O137([1, 5, 2, 3, 6, -2, 3, -3]); x
            1 + 5*i + 2*j + 3*k + 6*l - 2*li + 3*lj - 3*lk
            sage: O = OctonionAlgebra(ZZ)
            sage: y = O([1, 5, 2, 3, 6, -2, 3, -3])
            sage: O(x)
            1 + 5*i + 2*j + 3*k + 6*l - 2*li + 3*lj - 3*lk
            sage: O137(y)
            1 + 5*i + 2*j + 3*k + 6*l - 2*li + 3*lj - 3*lk

            sage: xp = O137([1, 5, 2, 3, 6, -2/3, 3, -3]); xp
            1 + 5*i + 2*j + 3*k + 6*l - 2/3*li + 3*lj - 3*lk
            sage: O(xp)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: O([1, 5, 2, 3, 6, -2, 3, -3/4])
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
        """
        if isinstance(x, Octonion_generic):
            v = self._module((<Octonion_generic> x).vec)
        elif self.base_ring().has_coerce_map_from(parent(x)):
            R = self.base_ring()
            v = self._module([R(x)] + [R.zero()]*7)
        else:
            v = self._module(x)
        return self.element_class(self, v)

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(ZZ)
            sage: O.some_elements()
            [2, 1, i, j, k, l, li, lj, lk,
             2 + 3*i + 4*j + 5*k + 6*l + 7*li + 8*lj + 9*lk,
             -2*j + 3*k - li - lj,
             8 - 7*i + 2*j + 13*k - 18*l + 45*li - 40*lj + 5*lk]
            sage: O = OctonionAlgebra(Zmod(6))
            sage: O.some_elements()
            [2, 1, i, j, k, l, li, lj, lk,
             2 + 3*i + 4*j + 5*k + li + 2*lj + 3*lk,
             4*j + 3*k + 5*li + 5*lj,
             2 + 5*i + 2*j + k + 3*li + 2*lj + 5*lk]
        """
        elts = [self.an_element()]
        elts.extend(self.basis())
        elt = sum(i * b for i, b in enumerate(self.basis(), start=2))
        elts.append(elt)
        elt2 = self([0, 0, -2, 3, 0, -1, -1, 0])
        elts.append(elt2)
        elts.append(elt * elt2)
        return elts

    @cached_method
    def one_basis(self):
        r"""
        Return the index for the basis element of `1`.

        EXAMPLES::

            sage: O = OctonionAlgebra(ZZ)
            sage: O.one_basis()
            0
        """
        return 0

    @cached_method
    def gens(self) -> tuple:
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(ZZ)
            sage: O.gens()
            (i, j, k, l)
        """
        e = self._module.basis()
        return tuple([self.element_class(self, v) for v in e[1:5]])

    @cached_method
    def basis(self):
        r"""
        Return the basis of ``self``.

        EXAMPLES::

            sage: O = OctonionAlgebra(ZZ)
            sage: O.basis()
            Family (1, i, j, k, l, li, lj, lk)
        """
        e = self._module.basis()
        from sage.sets.family import Family
        return Family([self.element_class(self, v) for v in e])

    Element = Octonion_generic
