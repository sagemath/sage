r"""
Orthogonal Linear Groups

The general orthogonal group `GO(n,R)` consists of all `n\times n`
matrices over the ring `R` preserving an `n`-ary positive definite
quadratic form. In cases where there are muliple non-isomorphic
quadratic forms, additional data needs to be specified to
disambiguate. The special orthogonal group is the normal subgroup of
matrices of determinant one.

In characteristics different from 2, a quadratic form is equivalent to
a bilinear symmetric form. Furthermore, over the real numbers a
positive definite quadratic form is equivalent to the diagonal
quadratic form, equivalent to the bilinear symmetric form defined by
the identity matrix. Hence, the orthogonal group `GO(n,\RR)` is the
group of orthogonal matrices in the usual sense.

In the case of a finite field and if the degree `n` is even, then
there are two inequivalent quadratic forms and a third parameter ``e``
must be specified to disambiguate these two possibilities. The index
of `SO(e,d,q)` in `GO(e,d,q)` is `2` if `q` is odd, but `SO(e,d,q) =
GO(e,d,q)` if `q` is even.)

.. warning::

   GAP and Sage use different notations:

   * GAP notation: The optional ``e`` comes first, that is, ``GO([e,]
     d, q)``, ``SO([e,] d, q)``.

   * Sage notation: The optional ``e`` comes last, the standard Python
     convention: ``GO(d, GF(q), e=0)``, ``SO( d, GF(q), e=0)``.

EXAMPLES::

    sage: GO(3,7)
    General Orthogonal Group of degree 3 over Finite Field of size 7

    sage: G = SO( 4, GF(7), 1); G
    Special Orthogonal Group of degree 4 and form parameter 1 over Finite Field of size 7
    sage: G.random_element()   # random
    [4 3 5 2]
    [6 6 4 0]
    [0 4 6 0]
    [4 4 5 1]

TESTS::

    sage: G = GO(3, GF(5))
    sage: latex(G)
    \text{GO}_{3}(\Bold{F}_{5})
    sage: G = SO(3, GF(5))
    sage: latex(G)
    \text{SO}_{3}(\Bold{F}_{5})
    sage: G = SO(4, GF(5), 1)
    sage: latex(G)
    \text{SO}_{4}(\Bold{F}_{5}, +)

AUTHORS:

- David Joyner (2006-03): initial version

- David Joyner (2006-05): added examples, _latex_, __str__, gens,
  as_matrix_group

- William Stein (2006-12-09): rewrite

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, NamedMatrixGroup_generic, NamedMatrixGroup_gap )



def normalize_args_e(degree, ring, e):
    """
    Normalize the arguments that relate the choice of quadratic form
    for special orthogonal groups over finite fields.

    INPUT:

    - ``degree`` -- integer. The degree of the affine group, that is,
      the dimension of the affine space the group is acting on.

    - ``ring`` -- a ring. The base ring of the affine space.

    - ``e`` -- integer, one of `+1`, `0`, `-1`.  Only relevant for
      finite fields and if the degree is even. A parameter that
      distinguishes inequivalent invariant forms.

    OUTPUT:

    The integer ``e`` with values required by GAP.

    TESTS::

        sage: from sage.groups.matrix_gps.orthogonal import normalize_args_e
        sage: normalize_args_e(2, GF(3), +1)
        1
        sage: normalize_args_e(3, GF(3), 0)
        0
        sage: normalize_args_e(3, GF(3), +1)
        0
        sage: normalize_args_e(2, GF(3), 0)
        Traceback (most recent call last):
        ...
        ValueError: must have e=-1 or e=1 for even degree
    """
    if is_FiniteField(ring) and degree%2 == 0:
        if e not in (-1, +1):
            raise ValueError('must have e=-1 or e=1 for even degree')
    else:
        e = 0
    return ZZ(e)





########################################################################
# General Orthogonal Group
########################################################################

def GO(n, R, e=0, var='a'):
    """
    Return the general orthogonal group.

    The general orthogonal group `GO(n,R)` consists of all `n\times n`
    matrices over the ring `R` preserving an `n`-ary positive definite
    quadratic form. In cases where there are muliple non-isomorphic
    quadratic forms, additional data needs to be specified to
    disambiguate.

    In the case of a finite field and if the degree `n` is even, then
    there are two inequivalent quadratic forms and a third parameter
    ``e`` must be specified to disambiguate these two possibilities.

    .. note::

        This group is also available via ``groups.matrix.GO()``.

   INPUT:

    - ``n`` -- integer. The degree.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``e`` -- ``+1`` or ``-1``, and ignored by default. Only relevant
      for finite fields and if the degree is even. A parameter that
      distinguishes inequivalent invariant forms.

    OUTPUT:

    The general orthogonal group of given degree, base ring, and
    choice of invariant form.

    EXAMPLES::

        sage: GO( 3, GF(7))
        General Orthogonal Group of degree 3 over Finite Field of size 7
        sage: GO( 3, GF(7)).order()
        672
        sage: GO( 3, GF(7)).gens()
        (
        [3 0 0]  [0 1 0]
        [0 5 0]  [1 6 6]
        [0 0 1], [0 2 1]
        )

    TESTS::

        sage: groups.matrix.GO(2, 3, e=-1)
        General Orthogonal Group of degree 2 and form parameter -1 over Finite Field of size 3
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    e = normalize_args_e(degree, ring, e)
    if e == 0:
        name = 'General Orthogonal Group of degree {0} over {1}'.format(degree, ring)
        ltx  = r'\text{{GO}}_{{{0}}}({1})'.format(degree, latex(ring))
    else:
        name = 'General Orthogonal Group of degree' + \
            ' {0} and form parameter {1} over {2}'.format(degree, e, ring)
        ltx  = r'\text{{GO}}_{{{0}}}({1}, {2})'.format(degree, latex(ring), '+' if e == 1 else '-')
    if is_FiniteField(ring):
        cmd  = 'GO({0}, {1}, {2})'.format(e, degree, ring.characteristic())
        return OrthogonalMatrixGroup_gap(degree, ring, False, name, ltx, cmd)
    else:
        return OrthogonalMatrixGroup_generic(degree, ring, False, name, ltx)



########################################################################
# Special Orthogonal Group
########################################################################

def SO(n, R, e=None, var='a'):
    """
    Return the special orthogonal group.

    The special orthogonal group `GO(n,R)` consists of all `n\times n`
    matrices with determint one over the ring `R` preserving an
    `n`-ary positive definite quadratic form. In cases where there are
    muliple non-isomorphic quadratic forms, additional data needs to
    be specified to disambiguate.

    .. note::

        This group is also available via ``groups.matrix.SO()``.

    INPUT:

    - ``n`` -- integer. The degree.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``e`` -- ``+1`` or ``-1``, and ignored by default. Only relevant
      for finite fields and if the degree is even. A parameter that
      distinguishes inequivalent invariant forms.

    OUTPUT:

    The special orthogonal group of given degree, base ring, and choice of
    invariant form.

    EXAMPLES::

        sage: G = SO(3,GF(5))
        sage: G
        Special Orthogonal Group of degree 3 over Finite Field of size 5

        sage: G = SO(3,GF(5))
        sage: G.gens()
        (
        [2 0 0]  [3 2 3]  [1 4 4]
        [0 3 0]  [0 2 0]  [4 0 0]
        [0 0 1], [0 3 1], [2 0 4]
        )
        sage: G = SO(3,GF(5))
        sage: G.as_matrix_group()
        Matrix group over Finite Field of size 5 with 3 generators (
        [2 0 0]  [3 2 3]  [1 4 4]
        [0 3 0]  [0 2 0]  [4 0 0]
        [0 0 1], [0 3 1], [2 0 4]
        )

    TESTS::

        sage: groups.matrix.SO(2, 3, e=1)
        Special Orthogonal Group of degree 2 and form parameter 1 over Finite Field of size 3
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    e = normalize_args_e(degree, ring, e)
    if e == 0:
        name = 'Special Orthogonal Group of degree {0} over {1}'.format(degree, ring)
        ltx  = r'\text{{SO}}_{{{0}}}({1})'.format(degree, latex(ring))
    else:
        name = 'Special Orthogonal Group of degree' + \
            ' {0} and form parameter {1} over {2}'.format(degree, e, ring)
        ltx  = r'\text{{SO}}_{{{0}}}({1}, {2})'.format(degree, latex(ring), '+' if e == 1 else '-')
    if is_FiniteField(ring):
        cmd  = 'SO({0}, {1}, {2})'.format(e, degree, ring.characteristic())
        return OrthogonalMatrixGroup_gap(degree, ring, True, name, ltx, cmd)
    else:
        return OrthogonalMatrixGroup_generic(degree, ring, True, name, ltx)



########################################################################
# Orthogonal Group class
########################################################################

class OrthogonalMatrixGroup_generic(NamedMatrixGroup_generic):

    @cached_method
    def invariant_bilinear_form(self):
        """
        Return the symmetric bilinear form preserved by the orthogonal
        group.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: GO(2,3,+1).invariant_bilinear_form()
            [0 1]
            [1 0]
            sage: GO(2,3,-1).invariant_bilinear_form()
            [2 1]
            [1 1]
        """
        from sage.matrix.constructor import identity_matrix
        m = identity_matrix(self.base_ring(), self.degree())
        m.set_immutable()
        return m

    @cached_method
    def invariant_quadratic_form(self):
        """
        Return the quadratic form preserved by the orthogonal group.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: GO(2,3,+1).invariant_quadratic_form()
            [0 1]
            [0 0]
            sage: GO(2,3,-1).invariant_quadratic_form()
            [1 1]
            [0 2]
        """
        from sage.matrix.constructor import identity_matrix
        m = identity_matrix(self.base_ring(), self.degree())
        m.set_immutable()
        return m

    def _check_matrix(self, x, *args):
        """a
        Check whether the matrix ``x`` is symplectic.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = GO(4, GF(5), +1)
            sage: G._check_matrix(G.an_element().matrix())
        """
        if self._special and x.determinant() != 1:
            raise TypeError('matrix must have determinant one')
        F = self.invariant_bilinear_form()
        if x * F * x.transpose() != F:
            raise TypeError('matrix must be orthogonal with respect to the invariant form')
        # TODO: check that quadratic form is preserved in characteristic two


class OrthogonalMatrixGroup_gap(OrthogonalMatrixGroup_generic, NamedMatrixGroup_gap):

    @cached_method
    def invariant_bilinear_form(self):
        """
        Return the symmetric bilinear form preserved by the orthogonal
        group.

        OUTPUT:

        A matrix `M` such that, for every group element g, the
        identity `g m g^T = m` holds. In characteristic different from
        two, this uniquely determines the orthogonal group.

        EXAMPLES::

            sage: G = GO(4, GF(7), -1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]

            sage: G = GO(4, GF(7), +1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 6 0]
            [0 0 0 2]

            sage: G = GO(4, QQ)
            sage: G.invariant_bilinear_form()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = SO(4, GF(7), -1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]
        """
        m = self.gap().InvariantBilinearForm()['matrix'].matrix()
        m.set_immutable()
        return m

    @cached_method
    def invariant_quadratic_form(self):
        """
        Return the quadratic form preserved by the orthogonal group.

        OUTPUT:

        The matrix `Q` defining "orthogonal" as follows. The matrix
        determines a quadratic form `q` on the natural vector space
        `V`, on which `G` acts, by `q(v) = v Q v^t`. A matrix `M' is
        an element of the orthogonal group if `q(v) = q(v M)` for all
        `v \in V`.

        EXAMPLES::

            sage: G = GO(4, GF(7), -1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = GO(4, GF(7), +1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 3 0]
            [0 0 0 1]

            sage: G = GO(4, QQ)
            sage: G.invariant_quadratic_form()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = SO(4, GF(7), -1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        m = self.gap().InvariantQuadraticForm()['matrix'].matrix()
        m.set_immutable()
        return m



