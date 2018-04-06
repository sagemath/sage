r"""
Unitary Groups `GU(n,q)` and `SU(n,q)`

These are `n \times n` unitary matrices with entries in
`GF(q^2)`.

EXAMPLES::

    sage: G = SU(3,5)
    sage: G.order()
    378000
    sage: G
    Special Unitary Group of degree 3 over Finite Field in a of size 5^2
    sage: G.gens()
    (
    [      a       0       0]  [4*a   4   1]
    [      0 2*a + 2       0]  [  4   4   0]
    [      0       0     3*a], [  1   0   0]
    )
    sage: G.base_ring()
    Finite Field in a of size 5^2

AUTHORS:

- David Joyner (2006-03): initial version, modified from
  special_linear (by W. Stein)

- David Joyner (2006-05): minor additions (examples, _latex_, __str__,
  gens)

- William Stein (2006-12): rewrite

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.
"""

#*********************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*********************************************************************************

from sage.rings.all import ZZ, GF
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.misc.latex import latex
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, NamedMatrixGroup_generic, NamedMatrixGroup_gap )


def finite_field_sqrt(ring):
    """
    Helper function.

    INPUT:

    A ring.

    OUTPUT:

    Integer q such that ``ring`` is the finite field with `q^2` elements.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.unitary import finite_field_sqrt
        sage: finite_field_sqrt(GF(4, 'a'))
        2
    """
    if not is_FiniteField(ring):
        raise ValueError('not a finite field')
    q, rem = ring.cardinality().sqrtrem()
    if rem:
        raise ValueError('cardinality not a square')
    return q


###############################################################################
# General Unitary Group
###############################################################################

def GU(n, R, var='a'):
    r"""
    Return the general unitary group.

    The general unitary group `GU( d, R )` consists of all `d \times
    d` matrices that preserve a nondegenerate sesquilinear form over
    the ring `R`.

    .. note::

        For a finite field the matrices that preserve a sesquilinear
        form over `F_q` live over `F_{q^2}`. So ``GU(n,q)`` for
        integer ``q`` constructs the matrix group over the base ring
        ``GF(q^2)``.

    .. note::

        This group is also available via ``groups.matrix.GU()``.

    INPUT:

    - ``n`` -- a positive integer.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``var`` -- variable used to represent generator of the finite
      field, if needed.

    OUTPUT:

    Return the general unitary group.

    EXAMPLES::

        sage: G = GU(3, 7); G
        General Unitary Group of degree 3 over Finite Field in a of size 7^2
        sage: G.gens()
        (
        [  a   0   0]  [6*a   6   1]
        [  0   1   0]  [  6   6   0]
        [  0   0 5*a], [  1   0   0]
        )
        sage: GU(2,QQ)
        General Unitary Group of degree 2 over Rational Field

        sage: G = GU(3, 5, var='beta')
        sage: G.base_ring()
        Finite Field in beta of size 5^2
        sage: G.gens()
        (
        [  beta      0      0]  [4*beta      4      1]
        [     0      1      0]  [     4      4      0]
        [     0      0 3*beta], [     1      0      0]
        )

    TESTS::

        sage: groups.matrix.GU(2, 3)
        General Unitary Group of degree 2 over Finite Field in a of size 3^2
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    if is_FiniteField(ring):
        q = ring.cardinality()
        ring = GF(q ** 2, name=var)
    name = 'General Unitary Group of degree {0} over {1}'.format(degree, ring)
    ltx  = r'\text{{GU}}_{{{0}}}({1})'.format(degree, latex(ring))
    if is_FiniteField(ring):
        cmd  = 'GU({0}, {1})'.format(degree, q)
        return UnitaryMatrixGroup_gap(degree, ring, False, name, ltx, cmd)
    else:
        return UnitaryMatrixGroup_generic(degree, ring, False, name, ltx)



###############################################################################
# Special Unitary Group
###############################################################################

def SU(n, R, var='a'):
    """
    The special unitary group `SU( d, R )` consists of all `d \times d`
    matrices that preserve a nondegenerate sesquilinear form over the
    ring `R` and have determinant one.

    .. note::

        For a finite field the matrices that preserve a sesquilinear
        form over `F_q` live over `F_{q^2}`. So ``SU(n,q)`` for
        integer ``q`` constructs the matrix group over the base ring
        ``GF(q^2)``.

    .. note::

        This group is also available via ``groups.matrix.SU()``.

    INPUT:

    - ``n`` -- a positive integer.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``var`` -- variable used to represent generator of the finite
      field, if needed.

    OUTPUT:

    Return the special unitary group.

    EXAMPLES::

        sage: SU(3,5)
        Special Unitary Group of degree 3 over Finite Field in a of size 5^2
        sage: SU(3, GF(5))
        Special Unitary Group of degree 3 over Finite Field in a of size 5^2
        sage: SU(3,QQ)
        Special Unitary Group of degree 3 over Rational Field

    TESTS::

        sage: groups.matrix.SU(2, 3)
        Special Unitary Group of degree 2 over Finite Field in a of size 3^2
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    if is_FiniteField(ring):
        q = ring.cardinality()
        ring = GF(q ** 2, name=var)
    name = 'Special Unitary Group of degree {0} over {1}'.format(degree, ring)
    ltx  = r'\text{{SU}}_{{{0}}}({1})'.format(degree, latex(ring))
    if is_FiniteField(ring):
        cmd  = 'SU({0}, {1})'.format(degree, q)
        return UnitaryMatrixGroup_gap(degree, ring, True, name, ltx, cmd)
    else:
        return UnitaryMatrixGroup_generic(degree, ring, True, name, ltx)


########################################################################
# Unitary Group class
########################################################################

class UnitaryMatrixGroup_generic(NamedMatrixGroup_generic):
    r"""
    General Unitary Group over arbitrary rings.

    EXAMPLES::

        sage: G = GU(3, GF(7)); G
        General Unitary Group of degree 3 over Finite Field in a of size 7^2
        sage: latex(G)
        \text{GU}_{3}(\Bold{F}_{7^{2}})

        sage: G = SU(3, GF(5));  G
        Special Unitary Group of degree 3 over Finite Field in a of size 5^2
        sage: latex(G)
        \text{SU}_{3}(\Bold{F}_{5^{2}})
    """

    def _check_matrix(self, x, *args):
        """a
        Check whether the matrix ``x`` is unitary.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = GU(2, GF(5))
            sage: G._check_matrix(G.an_element().matrix())
            sage: G = SU(2, GF(5))
            sage: G._check_matrix(G.an_element().matrix())
        """
        if self._special and x.determinant() != 1:
            raise TypeError('matrix must have determinant one')
        if not x.is_unitary():
            raise TypeError('matrix must be unitary')


class UnitaryMatrixGroup_gap(UnitaryMatrixGroup_generic, NamedMatrixGroup_gap):
    pass



