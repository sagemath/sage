"""
Symplectic Linear Groups

EXAMPLES::

    sage: G = Sp(4, GF(7));  G
    Symplectic Group of degree 4 over Finite Field of size 7

    sage: # needs sage.libs.gap
    sage: g = prod(G.gens());  g
    [3 0 3 0]
    [1 0 0 0]
    [0 1 0 1]
    [0 2 0 0]
    sage: m = g.matrix()
    sage: m * G.invariant_form() * m.transpose() == G.invariant_form()
    True
    sage: G.order()
    276595200

AUTHORS:

- David Joyner (2006-03): initial version, modified from
  special_linear (by W. Stein)

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.

- Sebastian Oehms (2018-8) add option for user defined invariant bilinear
  form and bug-fix in
  :meth:`~sage.groups.matrix_gps.symplectic.SymplecticMatrixGroup_generic.invariant_form`
  (see :issue:`26028`)
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#                     2013 Volker Braun <vbraun.name@gmail.com>
#                     2018 Sebastian Oehms
#                     2018 Travis Scrimshaw
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, normalize_args_invariant_form,
    NamedMatrixGroup_generic)


###############################################################################
# Symplectic Group
###############################################################################

def Sp(n, R, var='a', invariant_form=None):
    r"""
    Return the symplectic group.

    The special linear group `GL( d, R )` consists of all `d \times d`
    matrices that are invertible over the ring `R` with determinant one.

    .. NOTE::

        This group is also available via ``groups.matrix.Sp()``.

    INPUT:

    - ``n`` -- positive integer

    - ``R`` -- ring or an integer; if an integer is specified, the
      corresponding finite field is used

    - ``var`` -- (default: ``'a'``) variable used to
      represent generator of the finite field, if needed

    - ``invariant_form`` -- (optional) instances being accepted by
      the matrix-constructor which define a `n \times n` square matrix
      over ``R`` describing the alternating form to be kept invariant
      by the symplectic group

    EXAMPLES::

        sage: Sp(4, 5)
        Symplectic Group of degree 4 over Finite Field of size 5

        sage: Sp(4, IntegerModRing(15))
        Symplectic Group of degree 4 over Ring of integers modulo 15

        sage: Sp(3, GF(7))
        Traceback (most recent call last):
        ...
        ValueError: the degree must be even

    Using the ``invariant_form`` option::

        sage: m = matrix(QQ, 4,4, [[0, 0, 1, 0], [0, 0, 0, 2], [-1, 0, 0, 0], [0, -2, 0, 0]])
        sage: Sp4m = Sp(4, QQ, invariant_form=m)
        sage: Sp4 = Sp(4, QQ)
        sage: Sp4 == Sp4m
        False
        sage: Sp4.invariant_form()
        [ 0  0  0  1]
        [ 0  0  1  0]
        [ 0 -1  0  0]
        [-1  0  0  0]
        sage: Sp4m.invariant_form()
        [ 0  0  1  0]
        [ 0  0  0  2]
        [-1  0  0  0]
        [ 0 -2  0  0]
        sage: pm = Permutation([2,1,4,3]).to_matrix()
        sage: g = Sp4(pm); g in Sp4; g
        True
        [0 1 0 0]
        [1 0 0 0]
        [0 0 0 1]
        [0 0 1 0]
        sage: Sp4m(pm)
        Traceback (most recent call last):
        ...
        TypeError: matrix must be symplectic with respect to the alternating form
        [ 0  0  1  0]
        [ 0  0  0  2]
        [-1  0  0  0]
        [ 0 -2  0  0]

        sage: Sp(4,3, invariant_form=[[0,0,0,1],[0,0,1,0],[0,2,0,0], [2,0,0,0]])
        Traceback (most recent call last):
        ...
        NotImplementedError: invariant_form for finite groups is fixed by GAP

    TESTS::

        sage: TestSuite(Sp4).run()
        sage: TestSuite(Sp4m).run()
        sage: groups.matrix.Sp(2, 3)                                                    # needs sage.modules sage.rings.finite_rings
        Symplectic Group of degree 2 over Finite Field of size 3

        sage: # needs sage.libs.gap
        sage: G = Sp(4,5)
        sage: TestSuite(G).run()
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    if degree % 2:
        raise ValueError('the degree must be even')

    if invariant_form is not None:
        if isinstance(ring, FiniteField):
            raise NotImplementedError("invariant_form for finite groups is fixed by GAP")

        invariant_form = normalize_args_invariant_form(ring, degree, invariant_form)
        if not invariant_form.is_alternating():
            raise ValueError("invariant_form must be alternating")

        name = 'Symplectic Group of degree {0} over {1} with respect to alternating bilinear form\n{2}'.format(
                                                degree, ring, invariant_form)
        ltx = r'\text{{Sp}}_{{{0}}}({1})\text{{ with respect to alternating bilinear form}}{2}'.format(
                                    degree, latex(ring), latex(invariant_form))
    else:
        name = 'Symplectic Group of degree {0} over {1}'.format(degree, ring)
        ltx = r'\text{{Sp}}_{{{0}}}({1})'.format(degree, latex(ring))

    try:
        from .symplectic_gap import SymplecticMatrixGroup_gap
    except ImportError:
        pass
    else:
        try:
            cmd = 'Sp({0}, {1})'.format(degree, ring._gap_init_())
            return SymplecticMatrixGroup_gap(degree, ring, True, name, ltx, cmd)
        except ValueError:
            pass

    return SymplecticMatrixGroup_generic(degree, ring, True, name, ltx, invariant_form=invariant_form)


class SymplecticMatrixGroup_generic(NamedMatrixGroup_generic):
    r"""
    Symplectic Group over arbitrary rings.

    EXAMPLES::

        sage: Sp43 = Sp(4,3); Sp43
        Symplectic Group of degree 4 over Finite Field of size 3
        sage: latex(Sp43)
        \text{Sp}_{4}(\Bold{F}_{3})

        sage: Sp4m = Sp(4, QQ, invariant_form=(0, 0, 1, 0,  0, 0, 0, 2,
        ....:                                  -1, 0, 0, 0, 0, -2, 0, 0)); Sp4m
        Symplectic Group of degree 4 over Rational Field
         with respect to alternating bilinear form
        [ 0  0  1  0]
        [ 0  0  0  2]
        [-1  0  0  0]
        [ 0 -2  0  0]
        sage: latex(Sp4m)
        \text{Sp}_{4}(\Bold{Q})\text{ with respect to alternating bilinear form}\left(\begin{array}{rrrr}
        0 & 0 & 1 & 0 \\
        0 & 0 & 0 & 2 \\
        -1 & 0 & 0 & 0 \\
        0 & -2 & 0 & 0
        \end{array}\right)
    """

    @cached_method
    def invariant_form(self):
        """
        Return the quadratic form preserved by the symplectic group.

        OUTPUT: a matrix

        EXAMPLES::

            sage: Sp(4, QQ).invariant_form()
            [ 0  0  0  1]
            [ 0  0  1  0]
            [ 0 -1  0  0]
            [-1  0  0  0]
        """
        if self._invariant_form is not None:
            return self._invariant_form

        R = self.base_ring()
        d = self.degree()
        from sage.matrix.constructor import zero_matrix
        m = zero_matrix(R, d)
        for i in range(d):
            m[i, d-i-1] = 1 if i < d/2 else -1
        m.set_immutable()
        return m

    def _check_matrix(self, x, *args):
        """
        Check whether the matrix ``x`` is symplectic.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = Sp(4, GF(5))
            sage: G._check_matrix(G.an_element().matrix())
        """
        F = self.invariant_form()
        if x * F * x.transpose() != F:
            raise TypeError('matrix must be symplectic with respect to the alternating form\n{}'.format(F))
