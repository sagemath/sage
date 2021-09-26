"""
Linear Groups

EXAMPLES::

    sage: GL(4,QQ)
    General Linear Group of degree 4 over Rational Field
    sage: GL(1,ZZ)
    General Linear Group of degree 1 over Integer Ring
    sage: GL(100,RR)
    General Linear Group of degree 100 over Real Field with 53 bits of precision
    sage: GL(3,GF(49,'a'))
    General Linear Group of degree 3 over Finite Field in a of size 7^2

    sage: SL(2, ZZ)
    Special Linear Group of degree 2 over Integer Ring
    sage: G = SL(2,GF(3)); G
    Special Linear Group of degree 2 over Finite Field of size 3
    sage: G.is_finite()
    True
    sage: G.conjugacy_classes_representatives()
    (
    [1 0]  [0 2]  [0 1]  [2 0]  [0 2]  [0 1]  [0 2]
    [0 1], [1 1], [2 1], [0 2], [1 2], [2 2], [1 0]
    )
    sage: G = SL(6,GF(5))
    sage: G.gens()
    (
    [2 0 0 0 0 0]  [4 0 0 0 0 1]
    [0 3 0 0 0 0]  [4 0 0 0 0 0]
    [0 0 1 0 0 0]  [0 4 0 0 0 0]
    [0 0 0 1 0 0]  [0 0 4 0 0 0]
    [0 0 0 0 1 0]  [0 0 0 4 0 0]
    [0 0 0 0 0 1], [0 0 0 0 4 0]
    )

AUTHORS:

- William Stein: initial version

- David Joyner: degree, base_ring, random, order methods; examples

- David Joyner (2006-05): added center, more examples, renamed random
  attributes, bug fixes.

- William Stein (2006-12): total rewrite

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.

REFERENCES: See [KL1990]_ and [Car1972]_.
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.latex import latex
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, NamedMatrixGroup_generic, NamedMatrixGroup_gap )
from sage.categories.fields import Fields
from sage.categories.groups import Groups
from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap


###############################################################################
# General Linear Group
###############################################################################

def GL(n, R, var='a'):
    r"""
    Return the general linear group.

    The general linear group `GL( d, R )` consists of all `d \times d`
    matrices that are invertible over the ring `R`.

    .. NOTE::

        This group is also available via ``groups.matrix.GL()``.

    INPUT:

    - ``n`` -- a positive integer.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``var`` -- variable used to represent generator of the finite
      field, if needed.

    EXAMPLES::

        sage: G = GL(6,GF(5))
        sage: G.order()
        11064475422000000000000000
        sage: G.base_ring()
        Finite Field of size 5
        sage: G.category()
        Category of finite groups
        sage: TestSuite(G).run()

        sage: G = GL(6, QQ)
        sage: G.category()
        Category of infinite groups
        sage: TestSuite(G).run()

    Here is the Cayley graph of (relatively small) finite General Linear Group::

        sage: g = GL(2,3)
        sage: d = g.cayley_graph(); d
        Digraph on 48 vertices
        sage: d.plot(color_by_label=True, vertex_size=0.03, vertex_labels=False)  # long time
        Graphics object consisting of 144 graphics primitives
        sage: d.plot3d(color_by_label=True)  # long time
        Graphics3d Object

    ::

        sage: F = GF(3); MS = MatrixSpace(F,2,2)
        sage: gens = [MS([[2,0],[0,1]]), MS([[2,1],[2,0]])]
        sage: G = MatrixGroup(gens)
        sage: G.order()
        48
        sage: G.cardinality()
        48
        sage: H = GL(2,F)
        sage: H.order()
        48
        sage: H == G
        True
        sage: H.gens() == G.gens()
        True
        sage: H.as_matrix_group() == H
        True
        sage: H.gens()
        (
        [2 0]  [2 1]
        [0 1], [2 0]
        )

    TESTS::

        sage: groups.matrix.GL(2, 3)
        General Linear Group of degree 2 over Finite Field of size 3
        sage: groups.matrix.GL(1, ZZ).category()
        Category of groups
        sage: groups.matrix.GL(1, QQ).category()
        Category of infinite groups
    """
    degree, ring = normalize_args_vectorspace(n, R, var='a')
    try:
        if ring.is_finite():
            cat = Groups().Finite()
        elif n > 1 or ring in Fields():
            cat = Groups().Infinite()
        else:
            cat = Groups()
    except AttributeError:
        cat = Groups()
    name = 'General Linear Group of degree {0} over {1}'.format(degree, ring)
    ltx = 'GL({0}, {1})'.format(degree, latex(ring))
    try:
        cmd = 'GL({0}, {1})'.format(degree, ring._gap_init_())
        return LinearMatrixGroup_gap(degree, ring, False, name, ltx, cmd,
                                     category=cat)
    except ValueError:
        return LinearMatrixGroup_generic(degree, ring, False, name, ltx,
                                         category=cat)



###############################################################################
# Special Linear Group
###############################################################################

def SL(n, R, var='a'):
    r"""
    Return the special linear group.

    The special linear group `SL( d, R )` consists of all `d \times d`
    matrices that are invertible over the ring `R` with determinant
    one.

    .. note::

        This group is also available via ``groups.matrix.SL()``.

   INPUT:

    - ``n`` -- a positive integer.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``var`` -- variable used to represent generator of the finite
      field, if needed.

    EXAMPLES::

        sage: SL(3, GF(2))
        Special Linear Group of degree 3 over Finite Field of size 2
        sage: G = SL(15, GF(7)); G
        Special Linear Group of degree 15 over Finite Field of size 7
        sage: G.category()
        Category of finite groups
        sage: G.order()
        1956712595698146962015219062429586341124018007182049478916067369638713066737882363393519966343657677430907011270206265834819092046250232049187967718149558134226774650845658791865745408000000
        sage: len(G.gens())
        2
        sage: G = SL(2, ZZ); G
        Special Linear Group of degree 2 over Integer Ring
        sage: G.category()
        Category of infinite groups
        sage: G.gens()
        (
        [ 0  1]  [1 1]
        [-1  0], [0 1]
        )

    Next we compute generators for `\mathrm{SL}_3(\ZZ)` ::

        sage: G = SL(3,ZZ); G
        Special Linear Group of degree 3 over Integer Ring
        sage: G.gens()
        (
        [0 1 0]  [ 0  1  0]  [1 1 0]
        [0 0 1]  [-1  0  0]  [0 1 0]
        [1 0 0], [ 0  0  1], [0 0 1]
        )
        sage: TestSuite(G).run()

    TESTS::

        sage: groups.matrix.SL(2, 3)
        Special Linear Group of degree 2 over Finite Field of size 3
    """
    degree, ring = normalize_args_vectorspace(n, R, var='a')
    try:
        if ring.is_finite() or n == 1:
            cat = Groups().Finite()
        else:
            cat = Groups().Infinite()
    except AttributeError:
        cat = Groups()
    name = 'Special Linear Group of degree {0} over {1}'.format(degree, ring)
    ltx  = 'SL({0}, {1})'.format(degree, latex(ring))
    try:
        cmd  = 'SL({0}, {1})'.format(degree, ring._gap_init_())
        return LinearMatrixGroup_gap(degree, ring, True, name, ltx, cmd,
                                     category=cat)
    except ValueError:
        return LinearMatrixGroup_generic(degree, ring, True, name, ltx,
                                         category=cat)



########################################################################
# Linear Matrix Group class
########################################################################

class LinearMatrixGroup_generic(NamedMatrixGroup_generic):

    def _check_matrix(self, x, *args):
        """a
        Check whether the matrix ``x`` is special linear.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = SL(2,GF(5))
            sage: G._check_matrix(G.an_element().matrix())
        """
        if self._special:
            if x.determinant() != 1:
                raise TypeError('matrix must have determinant one')
        else:
            if x.determinant() == 0:
                raise TypeError('matrix must non-zero determinant')


class LinearMatrixGroup_gap(NamedMatrixGroup_gap, LinearMatrixGroup_generic, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The general or special linear group in GAP.

    TESTS:

    Check that :trac:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap
        sage: G = GL(3,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """
    pass
