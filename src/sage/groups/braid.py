"""
Braid groups

Braid groups are implemented as a particular case of finitely presented groups,
but with a lot of specific methods for braids.

A braid group can be created by giving the number of strands, and the name
of the generators::

    sage: BraidGroup(3)
    Braid group on 3 strands
    sage: BraidGroup(3,'a')
    Braid group on 3 strands
    sage: BraidGroup(3,'a').gens()
    (a0, a1)
    sage: BraidGroup(3,'a,b').gens()
    (a, b)

The elements can be created by operating with the generators, or by passing
a list with the indices of the letters to the group::

    sage: B.<s0,s1,s2> = BraidGroup(4)
    sage: s0*s1*s0
    s0*s1*s0
    sage: B([1,2,1])
    s0*s1*s0

The mapping class action of the braid group over the free group is
also implemented, see :class:`MappingClassGroupAction` for an
explanation. This action is left multiplication of a free group
element by a braid::

    sage: B.<b0,b1,b2> = BraidGroup()
    sage: F.<f0,f1,f2,f3> = FreeGroup()
    sage: B.strands() == F.rank()   # necessary for the action to be defined
    True
    sage: f1 * b1
    f1*f2*f1^-1
    sage: f0 * b1
    f0
    sage: f1 * b1
    f1*f2*f1^-1
    sage: f1^-1 * b1
    f1*f2^-1*f1^-1

AUTHORS:

- Miguel Angel Marco Buzunariz
- Volker Braun
- Søren Fuglede Jørgensen
- Robert Lipshitz
- Thierry Monteil: add a ``__hash__`` method consistent with the word
  problem to ensure correct Cayley graph computations.
- Sebastian Oehms (July and Nov 2018): add other versions for
  burau_matrix (unitary + simple, see :issue:`25760` and :issue:`26657`)
- Moritz Firsching (Sept 2021): Colored Jones polynomial
- Sebastian Oehms (May 2022): add :meth:`links_gould_polynomial`
"""

##############################################################################
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
##############################################################################

from itertools import combinations

from sage.algebras.free_algebra import FreeAlgebra
from sage.categories.action import Action
from sage.categories.groups import Groups
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.subset import Subsets
from sage.features.sagemath import sage__libs__braiding
from sage.functions.generalized import sign
from sage.groups.artin import FiniteTypeArtinGroup, FiniteTypeArtinGroupElement
from sage.groups.finitely_presented import (
    FinitelyPresentedGroup,
    GroupMorphismWithGensImages,
)
from sage.groups.free_group import FreeGroup, is_FreeGroup
from sage.groups.perm_gps.permgroup_named import (SymmetricGroup,
                                                  SymmetricGroupElement)
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.sets.set import Set
from sage.structure.element import Expression
from sage.structure.richcmp import rich_to_bool, richcmp

lazy_import('sage.libs.braiding',
            ['leftnormalform', 'rightnormalform', 'centralizer',
             'supersummitset', 'greatestcommondivisor',
             'leastcommonmultiple', 'conjugatingbraid', 'ultrasummitset',
             'thurston_type', 'rigidity', 'sliding_circuits', 'send_to_sss',
             'send_to_uss', 'send_to_sc', 'trajectory', 'cyclic_slidings'],
            feature=sage__libs__braiding())
lazy_import('sage.knots.knot', 'Knot')


class Braid(FiniteTypeArtinGroupElement):
    """
    An element of a braid group.

    It is a particular case of element of a finitely presented group.

    EXAMPLES::

        sage: B.<s0,s1,s2> = BraidGroup(4)
        sage: B
        Braid group on 4 strands
        sage: s0*s1/s2/s1
        s0*s1*s2^-1*s1^-1
        sage: B((1, 2, -3, -2))
        s0*s1*s2^-1*s1^-1
    """
    def _richcmp_(self, other, op):
        """
        Compare ``self`` and ``other``.

        TESTS::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1])
            sage: c = B([2, 1, 2])
            sage: b == c #indirect doctest
            True
            sage: b < c^(-1)
            True
            sage: B([]) == B.one()
            True
        """
        if self.Tietze() == other.Tietze():
            return rich_to_bool(op, 0)
        nfself = [i.Tietze() for i in self.left_normal_form()]
        nfother = [i.Tietze() for i in other.left_normal_form()]
        return richcmp(nfself, nfother, op)

    def __hash__(self):
        r"""
        Return a hash value for ``self``.

        EXAMPLES::

            sage: B.<s0,s1,s2> = BraidGroup(4)
            sage: hash(s0*s2) == hash(s2*s0)
            True
            sage: hash(s0*s1) == hash(s1*s0)
            False
        """
        return hash(tuple(i.Tietze() for i in self.left_normal_form()))

    def strands(self):
        """
        Return the number of strands in the braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, -1, 3, -2])
            sage: b.strands()
            4
        """
        return self.parent().strands()

    def components_in_closure(self):
        """
        Return the number of components of the trace closure of the braid.

        OUTPUT: a positive integer

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: b = B([1, -3])  # Three disjoint unknots
            sage: b.components_in_closure()
            3
            sage: b = B([1, 2, 3, 4])  # The unknot
            sage: b.components_in_closure()
            1
            sage: B = BraidGroup(4)
            sage: K11n42 = B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2])
            sage: K11n42.components_in_closure()
            1
        """
        cycles = self.permutation().to_cycles(singletons=False)
        return self.strands() - sum(len(c)-1 for c in cycles)

    def burau_matrix(self, var='t', reduced=False):
        r"""
        Return the Burau matrix of the braid.

        INPUT:

        - ``var`` -- string (default: ``'t'``); the name of the
          variable in the entries of the matrix
        - ``reduced`` -- boolean (default: ``False``); whether to
          return the reduced or unreduced Burau representation, can
          be one of the following:

          * ``True`` or ``'increasing'`` -- returns the reduced form using
            the basis given by `e_1 - e_i` for `2 \leq i \leq n`
          * ``'unitary'`` -- the unitary form according to Squier [Squ1984]_
          * ``'simple'`` -- returns the reduced form using the basis given
            by simple roots `e_i - e_{i+1}`, which yields the matrices
            given on the Wikipedia page

        OUTPUT:

        The Burau matrix of the braid. It is a matrix whose entries
        are Laurent polynomials in the variable ``var``. If ``reduced``
        is ``True``, return the matrix for the reduced Burau representation
        instead in the format specified. If ``reduced`` is ``'unitary'``,
        a triple ``M, Madj, H`` is returned, where ``M`` is the Burau matrix
        in the unitary form, ``Madj`` the adjoined to ``M`` and ``H``
        the hermitian form.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.inject_variables()
            Defining s0, s1, s2
            sage: b = s0*s1/s2/s1
            sage: b.burau_matrix()
            [       1 - t            0      t - t^2          t^2]
            [           1            0            0            0]
            [           0            0            1            0]
            [           0         t^-2 -t^-2 + t^-1    -t^-1 + 1]
            sage: s2.burau_matrix('x')
            [    1     0     0     0]
            [    0     1     0     0]
            [    0     0 1 - x     x]
            [    0     0     1     0]
            sage: s0.burau_matrix(reduced=True)
            [-t  0  0]
            [-t  1  0]
            [-t  0  1]

        Using the different reduced forms::

            sage: b.burau_matrix(reduced='simple')
            [    1 - t -t^-1 + 1        -1]
            [        1 -t^-1 + 1        -1]
            [        1     -t^-1         0]

            sage: M, Madj, H = b.burau_matrix(reduced='unitary')
            sage: M
            [-t^-2 + 1         t       t^2]
            [ t^-1 - t   1 - t^2      -t^3]
            [    -t^-2     -t^-1         0]
            sage: Madj
            [  1 - t^2 -t^-1 + t      -t^2]
            [     t^-1 -t^-2 + 1        -t]
            [     t^-2     -t^-3         0]
            sage: H
            [t^-1 + t       -1        0]
            [      -1 t^-1 + t       -1]
            [       0       -1 t^-1 + t]
            sage: M * H * Madj == H
            True

        REFERENCES:

        - :wikipedia:`Burau_representation`
        - [Squ1984]_
        """
        R = LaurentPolynomialRing(ZZ, var)
        t = R.gen()
        n = self.strands()
        if not reduced:
            M = identity_matrix(R, n)
            for i in self.Tietze():
                A = identity_matrix(R, n)
                if i > 0:
                    A[i-1, i-1] = 1-t
                    A[i, i] = 0
                    A[i, i-1] = 1
                    A[i-1, i] = t
                if i < 0:
                    A[-1-i, -1-i] = 0
                    A[-i, -i] = 1-t**(-1)
                    A[-1-i, -i] = 1
                    A[-i, -1-i] = t**(-1)
                M = M * A

        else:
            if reduced is True or reduced == "increasing":
                M = identity_matrix(R, n - 1)
                for j in self.Tietze():
                    A = identity_matrix(R, n - 1)
                    if j > 1:
                        i = j - 1
                        A[i-1, i-1] = 1 - t
                        A[i, i] = 0
                        A[i, i-1] = 1
                        A[i-1, i] = t
                    if j < -1:
                        i = j + 1
                        A[-1-i, -1-i] = 0
                        A[-i, -i] = 1 - t**-1
                        A[-1-i, -i] = 1
                        A[-i, -1-i] = t**-1
                    if j == 1:
                        for k in range(n - 1):
                            A[k, 0] = -t
                    if j == -1:
                        A[0, 0] = -t**-1
                        for k in range(1, n - 1):
                            A[k, 0] = -1
                    M = M * A

            elif reduced in ["simple", "unitary"]:
                M = identity_matrix(R, n - 1)
                for j in self.Tietze():
                    A = identity_matrix(R, n-1)
                    if j > 0:
                        A[j-1, j-1] = -t
                        if j > 1:
                            A[j-1, j-2] = t
                        if j < n-1:
                            A[j-1, j] = 1
                    if j < 0:
                        A[-j-1, -j-1] = -t**(-1)
                        if -j > 1:
                            A[-j-1, -j-2] = 1
                        if -j < n - 1:
                            A[-j-1, -j] = t**(-1)
                    M = M * A

            else:
                raise ValueError("invalid reduced type")

            if reduced == "unitary":
                # note: the roles of Madj and M are exchanged with respect
                # to the Squier paper in order to match the convention in
                # sage for instance in :meth:`_check_matrix` of
                # :class:`UnitaryMatrixGroup_generic`

                t_sq = R.hom([t**2], codomain=R)
                Madj = matrix(R, n - 1, n - 1,
                              lambda i, j: t**(j - i) * t_sq(M[i, j]))

                t_inv = R.hom([t**(-1)], codomain=R)
                M = matrix(R, n - 1, n - 1,
                           lambda i, j: t_inv(Madj[j, i]))

                # We see if the hermitian form has been cached
                #   in the parent
                H = self.parent()._hermitian_form
                if H is None:
                    # Defining the hermitian form
                    H = (t + t**(-1)) * identity_matrix(R, n - 1)
                    for i in range(n-2):
                        H[i, i + 1] = -1
                        H[i + 1, i] = -1
                    self.parent()._hermitian_form = H

                return M, Madj, H

        return M

    def alexander_polynomial(self, var='t', normalized=True):
        r"""
        Return the Alexander polynomial of the closure of the braid.

        INPUT:

        - ``var`` -- string (default: ``'t'``); the name of the
          variable in the entries of the matrix
        - ``normalized`` -- boolean (default: ``True``); whether to
          return the normalized Alexander polynomial

        OUTPUT: the Alexander polynomial of the braid closure of the braid

        This is computed using the reduced Burau representation. The
        unnormalized Alexander polynomial is a Laurent polynomial,
        which is only well-defined up to multiplication by plus or
        minus times a power of `t`.

        We normalize the polynomial by dividing by the largest power
        of `t` and then if the resulting constant coefficient
        is negative, we multiply by `-1`.

        EXAMPLES:

        We first construct the trefoil::

            sage: B = BraidGroup(3)
            sage: b = B([1,2,1,2])
            sage: b.alexander_polynomial(normalized=False)
            1 - t + t^2
            sage: b.alexander_polynomial()
            t^-2 - t^-1 + 1

        Next we construct the figure 8 knot::

            sage: b = B([-1,2,-1,2])
            sage: b.alexander_polynomial(normalized=False)
            -t^-2 + 3*t^-1 - 1
            sage: b.alexander_polynomial()
            t^-2 - 3*t^-1 + 1

        Our last example is the Kinoshita-Terasaka knot::

            sage: B = BraidGroup(4)
            sage: b = B([1,1,1,3,3,2,-3,-1,-1,2,-1,-3,-2])
            sage: b.alexander_polynomial(normalized=False)
            -t^-1
            sage: b.alexander_polynomial()
            1

        REFERENCES:

        - :wikipedia:`Alexander_polynomial`
        """
        n = self.strands()
        p = (self.burau_matrix(reduced=True) - identity_matrix(n - 1)).det()
        K, t = LaurentPolynomialRing(ZZ, var).objgen()
        if p == 0:
            return K.zero()
        qn = sum(t**i for i in range(n))
        p //= qn
        if normalized:
            p *= t**(-p.degree())
            if p.constant_coefficient() < 0:
                p = -p
        return p

    def permutation(self, W=None):
        """
        Return the permutation induced by the braid in its strands.

        INPUT:

        - ``W`` -- (optional) the permutation group to project
          ``self`` to; the default is ``self.parent().coxeter_group()``

        OUTPUT: the image of ``self`` under the natural projection map to ``W``

        EXAMPLES::

            sage: B.<s0,s1,s2> = BraidGroup()
            sage: S = SymmetricGroup(4)
            sage: b = s0*s1/s2/s1
            sage: c0 = b.permutation(W=S); c0
            (1,4,2)
            sage: c1 = b.permutation(W=Permutations(4)); c1
            [4, 1, 3, 2]
            sage: c1 == b.permutation()
            True

        The canonical section from the symmetric group to the braid group
        (sending a permutation to its associated permutation braid)
        can be recovered::

            sage: B(c0)
            s0*s1*s2*s1
            sage: B(c0) == B(c1)
            True
        """
        return self.coxeter_group_element(W)

    def plot(self, color='rainbow', orientation='bottom-top', gap=0.05,
             aspect_ratio=1, axes=False, **kwds):
        """
        Plot the braid.

        The following options are available:

        - ``color`` -- (default: ``'rainbow'``) the color of the
          strands. Possible values are:

            * ``'rainbow'``, uses :meth:`~sage.plot.colors.rainbow`
              according to the number of strands.

            * a valid color name for :meth:`~sage.plot.bezier_path`
              and :meth:`~sage.plot.line`. Used for all strands.

            * a list or a tuple of colors for each individual strand.

        - ``orientation`` -- (default: ``'bottom-top'``) determines how
          the braid is printed. The possible values are:

            * ``'bottom-top'``, the braid is printed from bottom to top

            * ``'top-bottom'``, the braid is printed from top to bottom

            * ``'left-right'``, the braid is printed from left to right

        - ``gap`` -- floating point number (default: 0.05); determines
          the size of the gap left when a strand goes under another

        - ``aspect_ratio`` -- floating point number (default:
          ``1``); the aspect ratio

        - ``**kwds`` -- other keyword options that are passed to
          :meth:`~sage.plot.bezier_path` and :meth:`~sage.plot.line`

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1])
            sage: b = B([2, 1, 2, 1])
            sage: c = b * a / b
            sage: d = a.conjugating_braid(c)
            sage: d * c / d == a
            True
            sage: d
            s1*s0
            sage: d * a / d == c
            False
            sage: B = BraidGroup(4, 's')
            sage: b = B([1, 2, 3, 1, 2, 1])
            sage: b.plot()                                                              # needs sage.plot
            Graphics object consisting of 30 graphics primitives
            sage: b.plot(color=["red", "blue", "red", "blue"])                          # needs sage.plot
            Graphics object consisting of 30 graphics primitives

            sage: B.<s,t> = BraidGroup(3)
            sage: b = t^-1*s^2
            sage: b.plot(orientation='left-right', color='red')                         # needs sage.plot
            Graphics object consisting of 12 graphics primitives
        """
        from sage.plot.bezier_path import bezier_path
        from sage.plot.colors import rainbow
        from sage.plot.plot import Graphics, line
        if orientation == 'top-bottom':
            orx = 0
            ory = -1
            nx = 1
            ny = 0
        elif orientation == 'left-right':
            orx = 1
            ory = 0
            nx = 0
            ny = -1
        elif orientation == 'bottom-top':
            orx = 0
            ory = 1
            nx = 1
            ny = 0
        else:
            raise ValueError('unknown value for "orientation"')
        n = self.strands()
        if isinstance(color, (list, tuple)):
            if len(color) != n:
                raise TypeError(f"color (={color}) must contain exactly {n} colors")
            col = list(color)
        elif color == "rainbow":
            col = rainbow(n)
        else:
            col = [color]*n
        braid = self.Tietze()
        a = Graphics()
        op = gap
        for i, m in enumerate(braid):
            for j in range(n):
                if m == j+1:
                    a += bezier_path([[(j*nx+i*orx, i*ory+j*ny), (j*nx+orx*(i+0.25), j*ny+ory*(i+0.25)),
                                       (nx*(j+0.5)+orx*(i+0.5), ny*(j+0.5)+ory*(i+0.5))],
                                      [(nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)),
                                       (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col[j], **kwds)
                elif m == j:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j-0.5+4*op)+orx*(i+0.5-2*op), ny*(j-0.5+4*op)+ory*(i+0.5-2*op)),
                                       (nx*(j-0.5+2*op)+orx*(i+0.5-op), ny*(j-0.5+2*op)+ory*(i+0.5-op))]],
                                     color=col[j], **kwds)
                    a += bezier_path([[(nx*(j-0.5-2*op)+orx*(i+0.5+op), ny*(j-0.5-2*op)+ory*(i+0.5+op)),
                                       (nx*(j-0.5-4*op)+orx*(i+0.5+2*op), ny*(j-0.5-4*op)+ory*(i+0.5+2*op)),
                                       (nx*(j-1)+orx*(i+0.75), ny*(j-1)+ory*(i+0.75)),
                                       (nx*(j-1)+orx*(i+1), ny*(j-1)+ory*(i+1))]], color=col[j], **kwds)
                    col[j], col[j-1] = col[j-1], col[j]
                elif -m == j+1:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j+0.5-4*op)+orx*(i+0.5-2*op), ny*(j+0.5-4*op)+ory*(i+0.5-2*op)),
                                       (nx*(j+0.5-2*op)+orx*(i+0.5-op), ny*(j+0.5-2*op)+ory*(i+0.5-op))]],
                                     color=col[j], **kwds)
                    a += bezier_path([[(nx*(j+0.5+2*op)+orx*(i+0.5+op), ny*(j+0.5+2*op)+ory*(i+0.5+op)),
                                       (nx*(j+0.5+4*op)+orx*(i+0.5+2*op), ny*(j+0.5+4*op)+ory*(i+0.5+2*op)),
                                       (nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)),
                                       (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col[j], **kwds)
                elif -m == j:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j-0.5)+orx*(i+0.5), ny*(j-0.5)+ory*(i+0.5))],
                                      [(nx*(j-1)+orx*(i+0.75), ny*(j-1)+ory*(i+0.75)),
                                       (nx*(j-1)+orx*(i+1), ny*(j-1)+ory*(i+1))]], color=col[j], **kwds)
                    col[j], col[j-1] = col[j-1], col[j]
                else:
                    a += line([(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+1), ny*j+ory*(i+1))], color=col[j], **kwds)
        a.set_aspect_ratio(aspect_ratio)
        a.axes(axes)
        return a

    def plot3d(self, color='rainbow'):
        """
        Plot the braid in 3d.

        The following option is available:

        - ``color`` -- (default: ``'rainbow'``) the color of the
          strands. Possible values are:

            * ``'rainbow'``, uses :meth:`~sage.plot.colors.rainbow`
              according to the number of strands.

            * a valid color name for :meth:`~sage.plot.plot3d.bezier3d`.
              Used for all strands.

            * a list or a tuple of colors for each individual strand.

        EXAMPLES::

            sage: B = BraidGroup(4, 's')
            sage: b = B([1, 2, 3, 1, 2, 1])
            sage: b.plot3d()                                                            # needs sage.plot sage.symbolic
            Graphics3d Object
            sage: b.plot3d(color='red')                                                 # needs sage.plot sage.symbolic
            Graphics3d Object
            sage: b.plot3d(color=["red", "blue", "red", "blue"])                        # needs sage.plot sage.symbolic
            Graphics3d Object
        """
        from sage.plot.colors import rainbow
        from sage.plot.plot3d.shapes2 import bezier3d
        b = []
        n = self.strands()
        if isinstance(color, (list, tuple)):
            if len(color) != n:
                raise TypeError("color (=%s) must contain exactly %d colors" % (color, n))
            col = list(color)
        elif color == "rainbow":
            col = rainbow(n)
        else:
            col = [color]*n
        braid = self.Tietze()

        for i, m in enumerate(braid):
            for j in range(n):
                if m == j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j+0.5, i+0.5)],
                                       [(0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]], color=col[j]))
                elif -m == j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j+0.5, i+0.5)],
                                       [(-0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]], color=col[j]))
                elif m == j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j-0.5, i+0.5)],
                                       [(-0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]], color=col[j]))
                    col[j], col[j-1] = col[j-1], col[j]
                elif -m == j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j-0.5, i+0.5)],
                                       [(0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]], color=col[j]))
                    col[j], col[j-1] = col[j-1], col[j]
                else:
                    b.append(bezier3d([[(0, j, i), (0, j, i+1)]], color=col[j]))
        return sum(b)

    def LKB_matrix(self, variables='x,y'):
        r"""
        Return the Lawrence-Krammer-Bigelow representation matrix.

        The matrix is expressed in the basis `\{e_{i, j} \mid 1\leq i
        < j \leq n\}`, where the indices are ordered
        lexicographically.  It is a matrix whose entries are in the
        ring of Laurent polynomials on the given variables.  By
        default, the variables are ``'x'`` and ``'y'``.

        INPUT:

        - ``variables`` -- string (default: ``'x,y'``); a string
          containing the names of the variables, separated by a comma

        OUTPUT: the matrix corresponding to the Lawrence-Krammer-Bigelow
        representation of the braid

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, 1])
            sage: b.LKB_matrix()
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: c = B([2, 1, 2])
            sage: c.LKB_matrix()
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]

        REFERENCES:

        - [Big2003]_
        """
        return self.parent()._LKB_matrix_(self.Tietze(), variab=variables)

    def TL_matrix(self, drain_size, variab=None, sparse=True):
        r"""
        Return the matrix representation of the Temperley--Lieb--Jones
        representation of the braid in a certain basis.

        The basis is given by non-intersecting pairings of `(n+d)` points,
        where `n` is the number of strands, `d` is given by ``drain_size``,
        and the pairings satisfy certain rules. See
        :meth:`~sage.groups.braid.BraidGroup_class.TL_basis_with_drain()`
        for details.

        We use the convention that the eigenvalues of the standard generators
        are `1` and `-A^4`, where `A` is a variable of a Laurent
        polynomial ring.

        When `d = n - 2` and the variables are picked appropriately, the
        resulting representation is equivalent to the reduced Burau
        representation.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands
          (both inclusive)

        - ``variab`` -- variable (default: ``None``); the variable in the
          entries of the matrices; if ``None``, then use a default variable
          in `\ZZ[A,A^{-1}]`

        - ``sparse`` -- boolean (default: ``True``); whether or not the
          result should be given as a sparse matrix

        OUTPUT: the matrix of the TL representation of the braid

        The parameter ``sparse`` can be set to ``False`` if it is
        expected that the resulting matrix will not be sparse. We
        currently make no attempt at guessing this.

        EXAMPLES:

        Let us calculate a few examples for `B_4` with `d = 0`::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, -3])
            sage: b.TL_matrix(0)
            [1 - A^4   -A^-2]
            [   -A^6       0]
            sage: R.<x> = LaurentPolynomialRing(GF(2))
            sage: b.TL_matrix(0, variab=x)
            [1 + x^4    x^-2]
            [    x^6       0]
            sage: b = B([])
            sage: b.TL_matrix(0)
            [1 0]
            [0 1]

        Test of one of the relations in `B_8`::

            sage: B = BraidGroup(8)
            sage: d = 0
            sage: B([4,5,4]).TL_matrix(d) == B([5,4,5]).TL_matrix(d)
            True

        An element of the kernel of the Burau representation, following
        [Big1999]_::

            sage: B = BraidGroup(6)
            sage: psi1 = B([4, -5, -2, 1])
            sage: psi2 = B([-4, 5, 5, 2, -1, -1])
            sage: w1 = psi1^(-1) * B([3]) * psi1
            sage: w2 = psi2^(-1) * B([3]) * psi2
            sage: (w1 * w2 * w1^(-1) * w2^(-1)).TL_matrix(4)
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]

        REFERENCES:

        - [Big1999]_
        - [Jon2005]_
        """
        if variab is None:
            R = LaurentPolynomialRing(ZZ, 'A')
        else:
            R = variab.parent()
        rep = self.parent().TL_representation(drain_size, variab)
        M = identity_matrix(R, self.parent().dimension_of_TL_space(drain_size),
                            sparse=sparse)
        for i in self.Tietze():
            if i > 0:
                M = M*rep[i-1][0]
            if i < 0:
                M = M*rep[-i-1][1]
        return M

    def links_gould_matrix(self, symbolics=False):
        r"""
        Return the representation matrix of ``self`` according to the R-matrix
        representation being attached to the quantum superalgebra `\mathfrak{sl}_q(2|1)`.

        See [MW2012]_, section 3 and references given there.

        INPUT:

        - ``symbolics`` -- boolean (default: ``False``); if set to ``True`` the
          coefficients will be contained in the symbolic ring. Per default they
          are elements of a quotient ring of a three variate Laurent polynomial
          ring.

        OUTPUT: the representation matrix of ``self`` over the ring according
        to the choice of the keyword ``symbolics`` (see the corresponding
        explanation)

        EXAMPLES::

            sage: Hopf = BraidGroup(2)([-1, -1])
            sage: HopfLG = Hopf.links_gould_matrix()
            sage: HopfLG.dimensions()
            (16, 16)
            sage: HopfLG.base_ring()
            Univariate Quotient Polynomial Ring in Yrbar
              over Multivariate Laurent Polynomial Ring in s0r, s1r
              over Integer Ring with modulus Yr^2 + s0r^2*s1r^2 - s0r^2 - s1r^2 + 1
            sage: HopfLGs = Hopf.links_gould_matrix(symbolics=True)                     # needs sage.symbolic
            sage: HopfLGs.base_ring()                                                   # needs sage.symbolic
            Symbolic Ring
        """
        rep = self.parent()._links_gould_representation(symbolics=symbolics)
        M = rep[0][0].parent().one()
        for i in self.Tietze():
            if i > 0:
                M = M * rep[i-1][0]
            if i < 0:
                M = M * rep[-i-1][1]
        return M

    @cached_method
    def links_gould_polynomial(self, varnames=None, use_symbolics=False):
        r"""
        Return the Links-Gould polynomial of the closure of ``self``.

        See [MW2012]_, section 3 and references given there.

        INPUT:

        - ``varnames`` -- string (default: ``'t0, t1'``)

        OUTPUT: a Laurent polynomial in the given variable names

        EXAMPLES::

            sage: Hopf = BraidGroup(2)([-1, -1])
            sage: Hopf.links_gould_polynomial()
            -1 + t1^-1 + t0^-1 - t0^-1*t1^-1
            sage: _ == Hopf.links_gould_polynomial(use_symbolics=True)
            True
            sage: Hopf.links_gould_polynomial(varnames='a, b')
            -1 + b^-1 + a^-1 - a^-1*b^-1
            sage: _ == Hopf.links_gould_polynomial(varnames='a, b', use_symbolics=True)
            True

        REFERENCES:

        - [MW2012]_
        """
        if varnames is not None:
            poly = self.links_gould_polynomial(use_symbolics=use_symbolics)
            R = LaurentPolynomialRing(ZZ, varnames)
            t0, t1 = R.gens()
            return poly(t0=t0, t1=t1)
        varnames = 't0, t1'

        rep = self.parent()._links_gould_representation(symbolics=use_symbolics)
        ln = len(rep)
        mu = rep[ln - 1]  # quantum trace factor
        M = mu * self.links_gould_matrix(symbolics=use_symbolics)
        d1, d2 = M.dimensions()
        e = d1 // 4
        B = M.base_ring()
        R = LaurentPolynomialRing(ZZ, varnames)

        # partial quantum trace according to I. Marin section 2.5
        part_trace = matrix(B, 4, 4, lambda i, j: sum(M[e * i + k, e * j + k]
                                                      for k in range(e)))
        ptemp = part_trace[0, 0]  # part_trace == psymb*M.parent().one()
        if use_symbolics:
            v1, v2 = R.variable_names()
            pstr = str(ptemp._sympy_().simplify())
            pstr = pstr.replace('t0', v1).replace('t1', v2)
            F = R.fraction_field()  # to make coercion work
            return R(F(pstr))
        else:
            ltemp = ptemp.lift().constant_coefficient()
            # Since the result of the calculation is known to be a Laurent polynomial
            # in t0 and t1 all exponents of ltemp must be divisible by 2
            L = ltemp.parent()
            lred = L({(k[0]/2, k[1]/2): v for k, v in ltemp.monomial_coefficients().items()})
            t0, t1 = R.gens()
            return lred(t0, t1)

    def tropical_coordinates(self) -> list:
        r"""
        Return the tropical coordinates of ``self`` in the braid group `B_n`.

        OUTPUT: list of `2n` tropical integers

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1])
            sage: tc = b.tropical_coordinates(); tc
            [1, 0, 0, 2, 0, 1]
            sage: tc[0].parent()
            Tropical semiring over Integer Ring

            sage: b = B([-2, -2, -1, -1, 2, 2, 1, 1])
            sage: b.tropical_coordinates()
            [1, -19, -12, 9, 0, 13]

        REFERENCES:

        - [DW2007]_
        - [Deh2011]_
        """
        coord = [0, 1] * self.strands()
        for s in self.Tietze():
            k = 2*(abs(s)-1)
            x1, y1, x2, y2 = coord[k:k+4]
            if s > 0:
                sign = 1
                z = x1 - min(y1, 0) - x2 + max(y2, 0)
                coord[k+1] = y2 - max(z, 0)
                coord[k+3] = y1 + max(z, 0)
            else:
                sign = -1
                z = x1 + min(y1, 0) - x2 - max(y2, 0)
                coord[k+1] = y2 + min(z, 0)
                coord[k+3] = y1 - min(z, 0)

            coord[k] = x1 + sign*(max(y1, 0) + max(max(y2, 0) - sign*z, 0))
            coord[k+2] = x2 + sign*(min(y2, 0) + min(min(y1, 0) + sign*z, 0))

        from sage.rings.semirings.tropical_semiring import TropicalSemiring
        T = TropicalSemiring(ZZ)
        return [T(c) for c in coord]

    def markov_trace(self, variab=None, normalized=True):
        r"""
        Return the Markov trace of the braid.

        The normalization is so that in the underlying braid group
        representation, the eigenvalues of the standard generators of
        the braid group are `1` and `-A^4`.

        INPUT:

        - ``variab`` -- variable (default: ``None``); the variable in the
          resulting polynomial; if ``None``, then use the variable `A`
          in `\ZZ[A,A^{-1}]`

        - ``normalized`` -- boolean (default: ``True``); if specified to be
          ``False``, return instead a rescaled Laurent polynomial version of
          the Markov trace

        OUTPUT:

        If ``normalized`` is ``False``, return instead the Markov trace
        of the braid, normalized by a factor of `(A^2+A^{-2})^n`. The
        result is then a Laurent polynomial in ``variab``. Otherwise it
        is a quotient of Laurent polynomials in ``variab``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, -3])
            sage: mt = b.markov_trace(); mt
            A^4/(A^12 + 3*A^8 + 3*A^4 + 1)
            sage: mt.factor()
            A^4 * (A^4 + 1)^-3

        We now give the non-normalized Markov trace::

            sage: mt = b.markov_trace(normalized=False); mt
            A^-4 + 1
            sage: mt.parent()
            Univariate Laurent Polynomial Ring in A over Integer Ring
        """
        if variab is None:
            R = LaurentPolynomialRing(ZZ, 'A')
            A = R.gens()[0]
            one = ZZ.one()
            quantum_integer = lambda d: R({i: one for i in range(-2*d, 2*d+1, 4)})
        else:
            A = variab
            quantum_integer = lambda d: (A**(2*(d+1))-A**(-2*(d+1))) // (A**2-A**(-2))

        n = self.strands()
        trace_sum = sum(quantum_integer(d) * self.TL_matrix(d, variab=variab).trace()
                        for d in range(n+1) if (n+d) % 2 == 0)

        if normalized:
            delta = A**2 + A**(-2)
            trace_sum = trace_sum / delta**n
        return trace_sum

    @lazy_attribute
    def _jones_polynomial(self):
        """
        Cached version of the Jones polynomial in a generic variable
        with the Skein normalization.

        The computation of the Jones polynomial uses the representation
        of the braid group on the Temperley--Lieb algebra. We cache the
        part of the calculation which does not depend on the choices of
        variables or normalizations.

        .. SEEALSO::

            :meth:`jones_polynomial`

        TESTS::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: b.jones_polynomial()                                                  # needs sage.symbolic
            1

            sage: B = BraidGroup(2)
            sage: b = B([])
            sage: b._jones_polynomial                                                   # needs sage.symbolic
            -A^-2 - A^2
            sage: b = B([-1, -1, -1])
            sage: b._jones_polynomial                                                   # needs sage.symbolic
            -A^-16 + A^-12 + A^-4
        """
        trace = self.markov_trace(normalized=False)
        A = trace.parent().gens()[0]
        D = A**2 + A**(-2)
        exp_sum = self.exponent_sum()
        num_comp = self.components_in_closure()
        return (-1)**(num_comp-1) * A**(2*exp_sum) * trace // D

    def jones_polynomial(self, variab=None, skein_normalization=False):
        r"""
        Return the Jones polynomial of the trace closure of the braid.

        The normalization is so that the unknot has Jones polynomial `1`. If
        ``skein_normalization`` is ``True``, the variable of the result is
        replaced by a itself to the power of `4`, so that the result
        agrees with the conventions of [Lic1997]_ (which in particular differs
        slightly from the conventions used otherwise in this class), had
        one used the conventional Kauffman bracket variable notation directly.

        If ``variab`` is ``None`` return a polynomial in the variable `A`
        or `t`, depending on the value ``skein_normalization``. In
        particular, if ``skein_normalization`` is ``False``, return the
        result in terms of the variable `t`, also used in [Lic1997]_.

        INPUT:

        - ``variab`` -- variable (default: ``None``); the variable in the
          resulting polynomial; if unspecified, use either a default variable
          in `\ZZ[A,A^{-1}]` or the variable `t` in the symbolic ring

        - ``skein_normalization`` -- boolean (default: ``False``); determines
          the variable of the resulting polynomial

        OUTPUT:

        If ``skein_normalization`` if ``False``, this returns an element
        in the symbolic ring as the Jones polynomial of the closure might
        have fractional powers when the closure of the braid is not a knot.
        Otherwise the result is a Laurent polynomial in ``variab``.

        EXAMPLES:

        The unknot::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: b.jones_polynomial()                                                  # needs sage.symbolic
            1

        Two unlinked unknots::

            sage: B = BraidGroup(2)
            sage: b = B([])
            sage: b.jones_polynomial()                                                  # needs sage.symbolic
            -sqrt(t) - 1/sqrt(t)

        The Hopf link::

            sage: B = BraidGroup(2)
            sage: b = B([-1,-1])
            sage: b.jones_polynomial()                                                  # needs sage.symbolic
            -1/sqrt(t) - 1/t^(5/2)

        Different representations of the trefoil and one of its mirror::

            sage: # needs sage.symbolic
            sage: B = BraidGroup(2)
            sage: b = B([-1, -1, -1])
            sage: b.jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: b.jones_polynomial()
            1/t + 1/t^3 - 1/t^4
            sage: B = BraidGroup(3)
            sage: b = B([-1, -2, -1, -2])
            sage: b.jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: R.<x> = LaurentPolynomialRing(GF(2))
            sage: b.jones_polynomial(skein_normalization=True, variab=x)
            x^-16 + x^-12 + x^-4
            sage: B = BraidGroup(3)
            sage: b = B([1, 2, 1, 2])
            sage: b.jones_polynomial(skein_normalization=True)
            A^4 + A^12 - A^16

        K11n42 (the mirror of the "Kinoshita-Terasaka" knot) and K11n34 (the
        mirror of the "Conway" knot)::

            sage: B = BraidGroup(4)
            sage: b11n42 = B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2])
            sage: b11n34 = B([1, 1, 2, -3, 2, -3, 1, -2, -2, -3, -3])
            sage: bool(b11n42.jones_polynomial() == b11n34.jones_polynomial())          # needs sage.symbolic
            True
        """
        if skein_normalization:
            if variab is None:
                return self._jones_polynomial
            else:
                return self._jones_polynomial(variab)
        else:
            if variab is None:
                variab = 't'
            if not isinstance(variab, Expression):
                from sage.symbolic.ring import SR
                variab = SR(variab)
            # We force the result to be in the symbolic ring because of the expand
            return self._jones_polynomial(variab**(ZZ(1)/ZZ(4))).expand()

    @cached_method
    def _enhanced_states(self):
        r"""
        Return the enhanced states of the closure of the braid diagram.

        The states are collected in a dictionary, where the dictionary
        keys are tuples of quantum and annular grading.
        Each dictionary value is itself a dictionary with the
        dictionary keys being the homological grading, and the values
        a list of enhanced states with the corresponding homology,
        quantum and annular grading.

        Each enhanced state is represented as a tuple containing:

        - A tuple with the type of smoothing made at each crossing.

        - A set with the circles marked as negative.

        - A set with the circles marked as positive.

        Each circle represented by a frozenset of tuples of the form
        (index of crossing, side where the circle passes the crossing)

        EXAMPLES::

            sage: B = BraidGroup(2)
            sage: b = B([1,1])
            sage: sorted((gr, sorted((d, [(sm,
            ....:                          sorted((sorted(A[0]), A[1]) for A in X),
            ....:                          sorted((sorted(A[0]), A[1]) for A in Y))
            ....:                         for sm, X, Y in data])
            ....:                    for d, data in v.items()))
            ....:        for gr,v in b._enhanced_states().items())
            [((0, -2),
              [(0, [((0, 0), [([(0, 1), (1, 1)], 1), ([(0, 3), (1, 3)], 1)], [])])]),
             ((2, 0),
              [(0,
                [((0, 0), [([(0, 3), (1, 3)], 1)], [([(0, 1), (1, 1)], 1)]),
                 ((0, 0), [([(0, 1), (1, 1)], 1)], [([(0, 3), (1, 3)], 1)])]),
               (1,
                [((1, 0), [([(0, 0), (0, 2), (1, 1), (1, 3)], 0)], []),
                 ((0, 1), [([(0, 1), (0, 3), (1, 0), (1, 2)], 0)], [])]),
               (2, [((1, 1), [([(0, 0), (1, 2)], 0), ([(0, 2), (1, 0)], 0)], [])])]),
             ((4, 0),
              [(1,
                [((1, 0), [], [([(0, 0), (0, 2), (1, 1), (1, 3)], 0)]),
                 ((0, 1), [], [([(0, 1), (0, 3), (1, 0), (1, 2)], 0)])]),
               (2,
                [((1, 1), [([(0, 2), (1, 0)], 0)], [([(0, 0), (1, 2)], 0)]),
                 ((1, 1), [([(0, 0), (1, 2)], 0)], [([(0, 2), (1, 0)], 0)])])]),
             ((4, 2),
              [(0, [((0, 0), [], [([(0, 1), (1, 1)], 1), ([(0, 3), (1, 3)], 1)])])]),
             ((6, 0),
              [(2, [((1, 1), [], [([(0, 0), (1, 2)], 0), ([(0, 2), (1, 0)], 0)])])])]
        """
        from sage.functions.generalized import sgn
        from sage.graphs.graph import Graph
        crossinglist = self.Tietze()
        ncross = len(crossinglist)
        writhe = 0

        # first build a "quadruply linked list", each crossing indicating its
        # previous and following neighbours
        last_crossing_in_row = [None] * self.strands()
        first_crossing_in_row = [None] * self.strands()
        crossings = [None] * ncross
        for i, cr in enumerate(crossinglist):
            writhe = writhe + sgn(cr)
            prevabove = last_crossing_in_row[abs(cr) - 1]
            prevbelow = last_crossing_in_row[abs(cr)]
            if prevabove is None:
                first_crossing_in_row[abs(cr) - 1] = i
            else:
                if abs(cr) == abs(crossings[prevabove]["cr"]):
                    crossings[prevabove]["next_above"] = i
                else:
                    crossings[prevabove]["next_below"] = i
            if prevbelow is None:
                first_crossing_in_row[abs(cr)] = i
            else:
                if abs(cr) == abs(crossings[prevbelow]["cr"]):
                    crossings[prevbelow]["next_below"] = i
                else:
                    crossings[prevbelow]["next_above"] = i
            crossings[i] = {"cr": cr,
                            "prev_above": prevabove,
                            "prev_below": prevbelow,
                            "next_above": None,
                            "next_below": None}
            last_crossing_in_row[abs(cr) - 1] = i
            last_crossing_in_row[abs(cr)] = i
        # tie up the ends of the list
        for k, i in enumerate(first_crossing_in_row):
            if i is not None:
                j = last_crossing_in_row[k]
                if abs(crossings[i]["cr"]) == k:
                    crossings[i]["prev_below"] = j
                else:
                    crossings[i]["prev_above"] = j

                if abs(crossings[j]["cr"]) == k:
                    crossings[j]["next_below"] = i
                else:
                    crossings[j]["next_above"] = i

        smoothings = []
        # generate all the resolutions
        for i in range(2**ncross):
            v = Integer(i).bits()
            v = v + [0]*(ncross - len(v))
            G = Graph()
            for j, cr in enumerate(crossings):
                if (v[j]*2-1)*sgn(cr["cr"]) == -1:  # oriented resolution
                    G.add_edge((j, cr["next_above"], abs(cr["cr"]) - 1), (j, 1))
                    G.add_edge((cr["prev_above"], j, abs(cr["cr"]) - 1), (j, 1))
                    G.add_edge((j, cr["next_below"], abs(cr["cr"])), (j, 3))
                    G.add_edge((cr["prev_below"], j, abs(cr["cr"])), (j, 3))
                else:
                    G.add_edge((j, cr["next_above"], abs(cr["cr"]) - 1), (j, 0))
                    G.add_edge((j, cr["next_below"], abs(cr["cr"])), (j, 0))
                    G.add_edge((cr["prev_above"], j, abs(cr["cr"]) - 1), (j, 2))
                    G.add_edge((cr["prev_below"], j, abs(cr["cr"])), (j, 2))
            # add loops of strands without crossing
            for k, j in enumerate(first_crossing_in_row):
                if j is None:
                    G.add_edge((ncross + k, ncross + k, k), (ncross + k, 4))
            sm = []
            for component in G.connected_components(sort=False):
                circle = set()
                trivial = 1
                # trivial switch: minus one means a circle is non-trivial.
                for vertex in component:
                    if len(vertex) == 3:
                        if vertex[1] <= vertex[0]:  # flip triviality for every looping edge
                            trivial *= -1
                    else:
                        circle.add(vertex)
                trivial = (1-trivial) // 2  # convert to 0 - trivial, 1 - non-trivial
                sm.append((frozenset(circle), trivial))
            smoothings.append((tuple(v), sm))

        states = {}
        for sm in smoothings:
            iindex = (writhe - ncross) // 2 + sum(sm[0])
            for m in range(2**len(sm[1])):
                m = [2*x-1 for x in Integer(m).bits()]
                m = m + [-1]*(len(sm[1]) - len(m))
                qagrad = (writhe + iindex + sum(m),
                          sum([x for i, x in enumerate(m) if sm[1][i][1] == 1]))
                circpos = set()
                circneg = set()
                for i, x in enumerate(m):
                    if x == 1:
                        circpos.add(sm[1][i])
                    else:
                        circneg.add(sm[1][i])

                if qagrad in states:
                    if iindex in states[qagrad]:
                        states[qagrad][iindex].append((sm[0], circneg, circpos))
                    else:
                        states[qagrad][iindex] = [(sm[0], circneg, circpos)]
                else:
                    states[qagrad] = {iindex: [(sm[0], circneg, circpos)]}
        return states

    @cached_method
    def _annular_khovanov_complex_cached(self, qagrad, ring=None):
        r"""
        Return the annular Khovanov complex of the braid.

        INPUT:

        - ``qagrad`` -- tuple of the quantum and annular grading to compute

        - ``ring`` -- (default: ``ZZ``) the coefficient ring

        OUTPUT: the annular Khovanov complex of the braid in the given grading

        .. NOTE::

            This method is intended only as the cache for
            :meth:`annular_khovanov_complex`.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: B([1,2,1,2])._annular_khovanov_complex_cached((5,-1)).homology()
            {1: Z, 2: Z, 3: 0}
        """
        from sage.homology.chain_complex import ChainComplex
        if ring is None:
            ring = ZZ
        states = self._enhanced_states()
        if qagrad in states:
            bases = states[qagrad]
        else:
            # return trivial chain complexx
            return ChainComplex()
        C_differentials = {}
        for i in bases:
            if i+1 in bases:
                m = matrix(ring, len(bases[i+1]), len(bases[i]), sparse=True)
                for ii in range(m.nrows()):
                    source = bases[i+1][ii]
                    for jj in range(m.ncols()):
                        target = bases[i][jj]
                        difs = [index for index, value in enumerate(source[0]) if value != target[0][index]]
                        if len(difs) == 1 and not (target[2].intersection(source[1]) or target[1].intersection(source[2])):
                            m[ii, jj] = (-1)**sum(target[0][:difs[0]])
            else:
                m = matrix(ring, 0, len(bases[i]), sparse=True)
            C_differentials[i] = m
        return ChainComplex(C_differentials)

    def annular_khovanov_complex(self, qagrad=None, ring=None):
        r"""
        Return the annular Khovanov complex of the closure of a braid,
        as defined in [BG2013]_.

        INPUT:

        - ``qagrad`` -- tuple of quantum and annular grading for which to compute
          the chain complex; if not specified all gradings are computed

        - ``ring`` -- (default: ``ZZ``) the coefficient ring

        OUTPUT:

        The annular Khovanov complex of the braid, given as a dictionary whose
        keys are tuples of quantum and annular grading. If ``qagrad`` is
        specified only return the chain complex of that grading.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1,-2,1,-2])
            sage: C = b.annular_khovanov_complex()
            sage: C
            {(-5, -1): Chain complex with at most 1 nonzero terms over Integer Ring,
             (-3, -3): Chain complex with at most 1 nonzero terms over Integer Ring,
             (-3, -1): Chain complex with at most 2 nonzero terms over Integer Ring,
             (-3, 1): Chain complex with at most 1 nonzero terms over Integer Ring,
             (-1, -1): Chain complex with at most 5 nonzero terms over Integer Ring,
             (-1, 1): Chain complex with at most 2 nonzero terms over Integer Ring,
             (1, -1): Chain complex with at most 2 nonzero terms over Integer Ring,
             (1, 1): Chain complex with at most 5 nonzero terms over Integer Ring,
             (3, -1): Chain complex with at most 1 nonzero terms over Integer Ring,
             (3, 1): Chain complex with at most 2 nonzero terms over Integer Ring,
             (3, 3): Chain complex with at most 1 nonzero terms over Integer Ring,
             (5, 1): Chain complex with at most 1 nonzero terms over Integer Ring}
            sage: C[1,-1].homology()
            {1: Z x Z, 2: 0}

        TESTS::

            sage: C = BraidGroup(2)([]).annular_khovanov_complex()
            sage: {qa: C[qa].homology() for qa in C}
            {(-2, -2): {0: Z}, (0, 0): {0: Z x Z}, (2, 2): {0: Z}}

            sage: BraidGroup(3)([-1]).annular_khovanov_complex((0,1), ZZ).differential()
            {-2: [],
             -1: [0]
             [1]
             [1],
             0: []}
        """
        if ring is None:
            ring = ZZ
        if qagrad is None:
            return {qa: self._annular_khovanov_complex_cached(qa, ring)
                    for qa in self._enhanced_states()}
        return self._annular_khovanov_complex_cached(qagrad, ring)

    def annular_khovanov_homology(self, qagrad=None, ring=ZZ):
        r"""
        Return the annular Khovanov homology of a closure of a braid.

        INPUT:

        - ``qagrad`` -- (optional) tuple of quantum and annular grading
          for which to compute the homology

        - ``ring`` -- (default: ``ZZ``) the coefficient ring

        OUTPUT:

        If ``qagrad`` is ``None``, return a dictionary of homologies in all
        gradings indexed by grading. If ``qagrad`` is specified, return the
        homology of that grading.

        .. NOTE::

            This is a simple wrapper around :meth:`annular_khovanov_complex`
            to compute homology from it.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1,3,-2])
            sage: b.annular_khovanov_homology()
            {(-3, -4): {0: Z},
             (-3, -2): {-1: Z},
             (-1, -2): {-1: 0, 0: Z x Z x Z, 1: 0},
             (-1, 0): {-1: Z x Z},
             (1, -2): {1: Z x Z},
             (1, 0): {-1: 0, 0: Z x Z x Z x Z, 1: 0, 2: 0},
             (1, 2): {-1: Z},
             (3, 0): {1: Z x Z x Z, 2: 0},
             (3, 2): {-1: 0, 0: Z x Z x Z, 1: 0},
             (5, 0): {2: Z},
             (5, 2): {1: Z x Z},
             (5, 4): {0: Z}}

            sage: B = BraidGroup(2)
            sage: b = B([1,1,1])
            sage: b.annular_khovanov_homology((7,0))
            {2: 0, 3: C2}

        TESTS::

            sage: b = BraidGroup(4)([1,-3])
            sage: b.annular_khovanov_homology((-4,-2))
            {-1: Z}
            sage: b.annular_khovanov_homology((0,2))
            {-1: Z}
        """
        if qagrad is None:
            C = self.annular_khovanov_complex(qagrad, ring)
            return {qa: C[qa].homology() for qa in C}
        return self.annular_khovanov_complex(qagrad, ring).homology()

    @cached_method
    def left_normal_form(self, algorithm='libbraiding'):
        r"""
        Return the left normal form of the braid.

        INPUT:

        - ``algorithm`` -- string (default: ``'artin'``); must be one of the following:

          * ``'artin'`` -- the general method for Artin groups is used
          * ``'libbraiding'`` -- the algorithm from the ``libbraiding`` package

        OUTPUT:

        A tuple of simple generators in the left normal form. The first
        element is a power of `\Delta`, and the rest are elements of the
        natural section lift from the corresponding symmetric group.

        EXAMPLES::

            sage: B = BraidGroup(6)
            sage: B.one().left_normal_form()
            (1,)
            sage: b = B([-2, 2, -4, -4, 4, -5, -1, 4, -1, 1])
            sage: L1 = b.left_normal_form(); L1
            (s0^-1*s1^-1*s0^-1*s2^-1*s1^-1*s0^-1*s3^-1*s2^-1*s1^-1*s0^-1*s4^-1*s3^-1*s2^-1*s1^-1*s0^-1,
             s0*s2*s1*s0*s3*s2*s1*s0*s4*s3*s2*s1,
             s3)
            sage: L1 == b.left_normal_form()
            True
            sage: B([1]).left_normal_form(algorithm='artin')
            (1, s0)
            sage: B([-3]).left_normal_form(algorithm='artin')
            (s0^-1*s1^-1*s0^-1*s2^-1*s1^-1*s0^-1*s3^-1*s2^-1*s1^-1*s0^-1*s4^-1*s3^-1*s2^-1*s1^-1*s0^-1,
             s0*s1*s2*s3*s4*s0*s1*s2*s3*s1*s2*s0*s1*s0)
            sage: B = BraidGroup(3)
            sage: B([1,2,-1]).left_normal_form()
            (s0^-1*s1^-1*s0^-1, s1*s0, s0*s1)
            sage: B([1,2,1]).left_normal_form()
            (s0*s1*s0,)
        """
        if algorithm == 'libbraiding':
            lnf = leftnormalform(self)
            B = self.parent()
            return tuple([B.delta()**lnf[0][0]] + [B(b) for b in lnf[1:]])
        elif algorithm == 'artin':
            return FiniteTypeArtinGroupElement.left_normal_form.f(self)
        raise ValueError("invalid algorithm")

    def _left_normal_form_coxeter(self):
        r"""
        Return the left normal form of the braid, in permutation form.

        OUTPUT: tuple whose first element is the power of `\Delta`, and the
        rest are the permutations corresponding to the simple factors

        EXAMPLES::

            sage: B = BraidGroup(12)
            sage: B([2, 2, 2, 3, 1, 2, 3, 2, 1, -2])._left_normal_form_coxeter()
            (-1,
             [12, 11, 10, 9, 8, 7, 6, 5, 2, 4, 3, 1],
             [4, 1, 3, 2, 5, 6, 7, 8, 9, 10, 11, 12],
             [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12],
             [3, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12],
             [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12])
            sage: C = BraidGroup(6)
            sage: C([2, 3, -4, 2, 3, -5, 1, -2, 3, 4, 1, -2])._left_normal_form_coxeter()
            (-2, [3, 5, 4, 2, 6, 1], [1, 6, 3, 5, 2, 4], [5, 6, 2, 4, 1, 3],
             [3, 2, 4, 1, 5, 6], [1, 5, 2, 3, 4, 6])

        .. NOTE::

            For long braids this method is slower than ``algorithm='libbraiding'``.

        .. TODO::

            Remove this method and use the default one from
            :meth:`sage.groups.artin.FiniteTypeArtinGroupElement.left_normal_form`.
        """
        delta = 0
        Delta = self.parent()._coxeter_group.long_element()
        sr = self.parent()._coxeter_group.simple_reflections()
        tz = self.Tietze()
        if not tz:
            return (0,)
        form = []
        for i in tz:
            if i > 0:
                form.append(sr[i])
            else:
                delta += 1
                form = [Delta * a * Delta for a in form]
                form.append(Delta * sr[-i])
        i = j = 0
        while j < len(form):
            while i < len(form) - j - 1:
                e = form[i].idescents(from_zero=False)
                s = form[i + 1].descents(from_zero=False)
                S = set(s).difference(set(e))
                while S:
                    a = list(S)[0]
                    form[i] = form[i] * sr[a]
                    form[i + 1] = sr[a] * form[i+1]
                    e = form[i].idescents(from_zero=False)
                    s = form[i + 1].descents(from_zero=False)
                    S = set(s).difference(set(e))
                if form[i+1].length() == 0:
                    form.pop(i+1)
                    i = 0
                else:
                    i += 1
            j += 1
            i = 0
        form = [a for a in form if a.length()]
        while form and form[0] == Delta:
            form.pop(0)
            delta -= 1
        return tuple([-delta] + form)

    def right_normal_form(self):
        r"""
        Return the right normal form of the braid.

        A tuple of simple generators in the right normal form. The last
        element is a power of `\Delta`, and the rest are elements of the
        natural section lift from the corresponding symmetric group.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1, -2, 3, 1])
            sage: b.right_normal_form()
            (s1*s0, s0*s2, 1)
        """
        rnf = rightnormalform(self)
        B = self.parent()
        return tuple([B(b) for b in rnf[:-1]] + [B.delta()**rnf[-1][0]])

    def centralizer(self) -> list:
        """
        Return a list of generators of the centralizer of the braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([2, 1, 3, 2])
            sage: b.centralizer()
            [s1*s0*s2*s1, s0*s2]
        """
        c = centralizer(self)
        B = self.parent()
        return [B._element_from_libbraiding(b) for b in c]

    def super_summit_set(self) -> list:
        """
        Return a list with the super summit set of the braid.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1, -2, -2, 1])
            sage: b.super_summit_set()
            [s0^-1*s1^-1*s0^-2*s1^2*s0^2,
             (s0^-1*s1^-1*s0^-1)^2*s1^2*s0^3*s1,
             (s0^-1*s1^-1*s0^-1)^2*s1*s0^3*s1^2,
             s0^-1*s1^-1*s0^-2*s1^-1*s0*s1^3*s0]
        """
        sss = supersummitset(self)
        B = self.parent()
        return [B._element_from_libbraiding(b) for b in sss]

    def gcd(self, other):
        """
        Return the greatest common divisor of the two braids.

        INPUT:

        - ``other`` -- the other braid with respect with the gcd is computed

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1, -2, -2, 1])
            sage: c = B([1, 2, 1])
            sage: b.gcd(c)
            s0^-1*s1^-1*s0^-2*s1^2*s0
            sage: c.gcd(b)
            s0^-1*s1^-1*s0^-2*s1^2*s0
        """
        B = self.parent()
        b = greatestcommondivisor(self, other)
        return B._element_from_libbraiding(b)

    def lcm(self, other):
        """
        Return the least common multiple of the two braids.

        INPUT:

        - ``other`` -- the other braid with respect with the lcm is computed

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1, -2, -2, 1])
            sage: c = B([1, 2, 1])
            sage: b.lcm(c)
            (s0*s1)^2*s0
        """
        B = self.parent()
        b = leastcommonmultiple(self, other)
        return B._element_from_libbraiding(b)

    def conjugating_braid(self, other):
        r"""
        Return a conjugating braid, if it exists.

        INPUT:

        - ``other`` -- a braid in the same braid group as ``self``

        OUTPUT:

        A conjugating braid. More precisely, if the output is `d`, `o` equals
        ``other``, and `s` equals ``self`` then `o = d^{-1} \cdot s \cdot d`.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: B.one().conjugating_braid(B.one())
            1
            sage: B.one().conjugating_braid(B.gen(0)) is None
            True
            sage: B.gen(0).conjugating_braid(B.gen(1))
            s1*s0
            sage: B.gen(0).conjugating_braid(B.gen(1).inverse()) is None
            True
            sage: a = B([2, 2, -1, -1])
            sage: b = B([2, 1, 2, 1])
            sage: c = b * a / b
            sage: d1 = a.conjugating_braid(c)
            sage: d1
            s1*s0
            sage: d1 * c / d1 == a
            True
            sage: d1 * a / d1 == c
            False
            sage: l = sage.groups.braid.conjugatingbraid(a,c)                           # needs sage.groups
            sage: d1 == B._element_from_libbraiding(l)                                  # needs sage.groups
            True
            sage: b = B([2, 2, 2, 2, 1])
            sage: c = b * a / b
            sage: d1 = a.conjugating_braid(c)
            sage: len(d1.Tietze())
            7
            sage: d1 * c / d1 == a
            True
            sage: d1 * a / d1 == c
            False
            sage: d1
            s1^2*s0^2*s1^2*s0
            sage: l = sage.groups.braid.conjugatingbraid(a,c)                           # needs sage.groups
            sage: d2 = B._element_from_libbraiding(l)                                   # needs sage.groups
            sage: len(d2.Tietze())                                                      # needs sage.groups
            13
            sage: c.conjugating_braid(b) is None
            True
        """
        cb = conjugatingbraid(self, other)
        if not cb:
            return None
        B = self.parent()
        cb[0][0] %= 2
        return B._element_from_libbraiding(cb)

    def is_conjugated(self, other) -> bool:
        """
        Check if the two braids are conjugated.

        INPUT:

        - ``other`` -- the other braid to check for conjugacy

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1])
            sage: b = B([2, 1, 2, 1])
            sage: c = b * a / b
            sage: c.is_conjugated(a)
            True
            sage: c.is_conjugated(b)
            False
        """
        cb = conjugatingbraid(self, other)
        return bool(cb)

    def pure_conjugating_braid(self, other):
        r"""
        Return a pure conjugating braid, i.e. a conjugating braid whose
        associated permutation is the identity, if it exists.

        INPUT:

        - ``other`` -- a braid in the same braid group as ``self``

        OUTPUT:

        A pure conjugating braid. More precisely, if the output is `d`, `o`
        equals ``other``, and `s` equals ``self`` then
        `o = d^{-1} \cdot s \cdot d`.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.one().pure_conjugating_braid(B.one())
            1
            sage: B.one().pure_conjugating_braid(B.gen(0)) is None
            True
            sage: B.gen(0).pure_conjugating_braid(B.gen(1)) is None
            True
            sage: B.gen(0).conjugating_braid(B.gen(2).inverse()) is None
            True
            sage: a = B([1, 2, 3])
            sage: b = B([3, 2,])
            sage: c = b ^ 12 * a / b ^ 12
            sage: d1 = a.conjugating_braid(c)
            sage: len(d1.Tietze())
            30
            sage: S = SymmetricGroup(4)
            sage: d1.permutation(W=S)
            (1,3)(2,4)
            sage: d1 * c / d1 == a
            True
            sage: d1 * a / d1 == c
            False
            sage: d2 = a.pure_conjugating_braid(c)
            sage: len(d2.Tietze())
            24
            sage: d2.permutation(W=S)
            ()
            sage: d2 * c / d2 == a
            True
            sage: d2
            (s0*s1*s2^2*s1*s0)^4
            sage: a.conjugating_braid(b) is None
            True
            sage: a.pure_conjugating_braid(b) is None
            True
            sage: a1 = B([1])
            sage: a2 = B([2])
            sage: a1.conjugating_braid(a2)
            s1*s0
            sage: a1.permutation(W=S)
            (1,2)
            sage: a2.permutation(W=S)
            (2,3)
            sage: a1.pure_conjugating_braid(a2) is None
            True
            sage: (a1^2).conjugating_braid(a2^2)
            s1*s0
            sage: (a1^2).pure_conjugating_braid(a2^2) is None
            True
        """
        B = self.parent()
        n = B.strands()
        S = SymmetricGroup(n)
        p1 = self.permutation(W=S)
        p2 = other.permutation(W=S)
        if p1 != p2:
            return None
        b0 = self.conjugating_braid(other)
        if b0 is None:
            return None
        p3 = b0.permutation(W=S).inverse()
        if p3.is_one():
            return b0
        LP = {a.permutation(W=S): a for a in self.centralizer()}
        if p3 not in S.subgroup(LP):
            return None
        P = p3.word_problem(list(LP), display=False, as_list=True)
        b1 = prod(LP[S(a)] ** b for a, b in P)
        b0 = b1 * b0
        n0 = len(b0.Tietze())
        L = leftnormalform(b0)
        L[0][0] %= 2
        b2 = B._element_from_libbraiding(L)
        n2 = len(b2.Tietze())
        return b2 if n2 <= n0 else b0

    def ultra_summit_set(self) -> list:
        """
        Return a list with the orbits of the ultra summit set of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: b = B([2, 1, 2, 1])
            sage: b.ultra_summit_set()
            [[s0*s1*s0^2, (s0*s1)^2]]
            sage: a.ultra_summit_set()
            [[(s0^-1*s1^-1*s0^-1)^2*s1^3*s0^2*s1^3,
              (s0^-1*s1^-1*s0^-1)^2*s1^2*s0^2*s1^4,
              (s0^-1*s1^-1*s0^-1)^2*s1*s0^2*s1^5,
              s0^-1*s1^-1*s0^-2*s1^5*s0,
              (s0^-1*s1^-1*s0^-1)^2*s1^5*s0^2*s1,
              (s0^-1*s1^-1*s0^-1)^2*s1^4*s0^2*s1^2],
             [s0^-1*s1^-1*s0^-2*s1^-1*s0^2*s1^2*s0^3,
              s0^-1*s1^-1*s0^-2*s1^-1*s0*s1^2*s0^4,
              s0^-1*s1^-1*s0^-2*s1*s0^5,
              (s0^-1*s1^-1*s0^-1)^2*s1*s0^6*s1,
              s0^-1*s1^-1*s0^-2*s1^-1*s0^4*s1^2*s0,
              s0^-1*s1^-1*s0^-2*s1^-1*s0^3*s1^2*s0^2]]
        """
        uss = ultrasummitset(self)
        B = self.parent()
        return [[B._element_from_libbraiding(i) for i in s] for s in uss]

    def thurston_type(self) -> str:
        """
        Return the thurston_type of ``self``.

        OUTPUT: one of ``'reducible'``, ``'periodic'`` or ``'pseudo-anosov'``

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1])
            sage: b.thurston_type()
            'reducible'
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.thurston_type()
            'pseudo-anosov'
            sage: c = B([2, 1, 2, 1])
            sage: c.thurston_type()
            'periodic'
        """
        return thurston_type(self)

    def is_reducible(self) -> bool:
        """
        Check whether the braid is reducible.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1])
            sage: b.is_reducible()
            True
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.is_reducible()
            False
        """
        return self.thurston_type() == 'reducible'

    def is_periodic(self) -> bool:
        """
        Check whether the braid is periodic.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: b = B([2, 1, 2, 1])
            sage: a.is_periodic()
            False
            sage: b.is_periodic()
            True
        """
        return self.thurston_type() == 'periodic'

    def is_pseudoanosov(self) -> bool:
        """
        Check if the braid is pseudo-Anosov.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: b = B([2, 1, 2, 1])
            sage: a.is_pseudoanosov()
            True
            sage: b.is_pseudoanosov()
            False
        """
        return self.thurston_type() == 'pseudo-anosov'

    def rigidity(self):
        """
        Return the rigidity of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([2, 1, 2, 1])
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.rigidity()
            6
            sage: b.rigidity()
            0
        """
        return Integer(rigidity(self))

    def sliding_circuits(self) -> list:
        """
        Return the sliding circuits of the braid.

        OUTPUT: list of sliding circuits. Each sliding circuit is itself
        a list of braids.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.sliding_circuits()
            [[(s0^-1*s1^-1*s0^-1)^2*s1^3*s0^2*s1^3],
             [s0^-1*s1^-1*s0^-2*s1^-1*s0^2*s1^2*s0^3],
             [s0^-1*s1^-1*s0^-2*s1^-1*s0^3*s1^2*s0^2],
             [(s0^-1*s1^-1*s0^-1)^2*s1^4*s0^2*s1^2],
             [(s0^-1*s1^-1*s0^-1)^2*s1^2*s0^2*s1^4],
             [s0^-1*s1^-1*s0^-2*s1^-1*s0*s1^2*s0^4],
             [(s0^-1*s1^-1*s0^-1)^2*s1^5*s0^2*s1],
             [s0^-1*s1^-1*s0^-2*s1^-1*s0^4*s1^2*s0],
             [(s0^-1*s1^-1*s0^-1)^2*s1*s0^2*s1^5],
             [s0^-1*s1^-1*s0^-2*s1*s0^5],
             [(s0^-1*s1^-1*s0^-1)^2*s1*s0^6*s1],
             [s0^-1*s1^-1*s0^-2*s1^5*s0]]
            sage: b = B([2, 1, 2, 1])
            sage: b.sliding_circuits()
            [[s0*s1*s0^2, (s0*s1)^2]]
        """
        slc = sliding_circuits(self)
        B = self.parent()
        return [[B._element_from_libbraiding(i) for i in s] for s in slc]

    def mirror_image(self):
        r"""
        Return the image of ``self`` under the mirror involution (see
        :meth:`BraidGroup_class.mirror_involution`). The link closure of
        it is mirrored to the closure of ``self`` (see the example below
        of a positive amphicheiral knot).

        EXAMPLES::

            sage: B5 = BraidGroup(5)
            sage: b  = B5((-1, 2, -3, -1, -3, 4, 2, -3, 2, 4, 2, -3)) # closure K12a_427
            sage: bm = b.mirror_image(); bm
            s0*s1^-1*s2*s0*s2*s3^-1*s1^-1*s2*s1^-1*s3^-1*s1^-1*s2
            sage: bm.is_conjugated(b)
            True
            sage: bm.is_conjugated(~b)
            False
        """
        return self.parent().mirror_involution()(self)

    def reverse(self):
        r"""
        Return the reverse of ``self`` obtained by reversing the order of the
        generators in its word. This defines an anti-involution on the braid
        group. The link closure of it has the reversed orientation (see the
        example below of a non reversible knot).

        EXAMPLES::

            sage: b  = BraidGroup(3)((1, 1, -2, 1, -2, 1, -2, -2))  # closure K8_17
            sage: br = b.reverse(); br
            s1^-1*(s1^-1*s0)^3*s0
            sage: br.is_conjugated(b)
            False
        """
        t = list(self.Tietze())
        t.reverse()
        return self.parent()(tuple(t))

    def deformed_burau_matrix(self, variab='q'):
        r"""
        Return the deformed Burau matrix of the braid.

        INPUT:

        - ``variab`` -- variable (default: ``q``); the variable in the
          resulting laurent polynomial, which is the base ring for the
          free algebra constructed

        OUTPUT: a matrix with elements in the free algebra ``self._algebra``

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, -3, -2, 3, 1])
            sage: db = b.deformed_burau_matrix(); db
            [                                ap_0*ap_5 ... bp_0*ap_1*cm_3*bp_4]
            ...
            [                           bm_2*bm_3*cp_5 ...      bm_2*am_3*bp_4]

        We check how this relates to the nondeformed Burau matrix::

            sage: def subs_gen(gen, q):
            ....:     gen_str = str(gen)
            ....:     v = q if 'p' in gen_str else 1/q
            ....:     if 'b' in gen_str:
            ....:         return v
            ....:     elif 'a' in gen_str:
            ....:         return 1 - v
            ....:     else:
            ....:         return 1
            sage: db_base = db.parent().base_ring()
            sage: q = db_base.base_ring().gen()
            sage: db_simp = db.subs({gen: subs_gen(gen, q)
            ....:                    for gen in db_base.gens()})
            sage: db_simp
            [ (1-2*q+q^2)      (q-q^2)  (q-q^2+q^3)    (q^2-q^3)]
            [       (1-q)            q            0            0]
            [           0            0        (1-q)            q]
            [      (q^-2)            0 -(q^-2-q^-1)    -(q^-1-1)]
            sage: burau = b.burau_matrix(); burau
            [1 - 2*t + t^2       t - t^2 t - t^2 + t^3     t^2 - t^3]
            [        1 - t             t             0             0]
            [            0             0         1 - t             t]
            [         t^-2             0  -t^-2 + t^-1     -t^-1 + 1]
            sage: t = burau.parent().base_ring().gen()
            sage: burau.subs({t:q}).change_ring(db_base) == db_simp
            True
        """
        R = LaurentPolynomialRing(ZZ, variab)

        n = self.strands()
        tz = self.Tietze()
        m3 = len(tz) * 3
        plus = [i for i, tzi in enumerate(tz) if tzi > 0]
        minus = [i for i, tzi in enumerate(tz) if tzi < 0]
        gens_str = [f'{s}p_{i}' for i in plus for s in 'bca']
        gens_str += [f'{s}m_{i}' for i in minus for s in 'bca']
        alg_ZZ = FreeAlgebra(ZZ, m3, gens_str)
        gen_indices = {k: i for i, k in enumerate(plus + minus)}
        gens = [alg_ZZ.gens()[k:k + 3] for k in range(0, m3, 3)]

        M = identity_matrix(alg_ZZ, n)
        for k, i in enumerate(tz):
            A = identity_matrix(alg_ZZ, n)
            b, c, a = gens[gen_indices[k]]
            # faster using row operations instead ?
            if i > 0:
                A[i-1, i-1] = a
                A[i, i] = 0
                A[i, i-1] = c
                A[i-1, i] = b
            if i < 0:
                A[-1-i, -1-i] = 0
                A[-i, -i] = a
                A[-1-i, -i] = c
                A[-i, -1-i] = b
            M = M * A

        alg_R = FreeAlgebra(R, m3, gens_str)
        return M.change_ring(alg_R)

    def _colored_jones_sum(self, N, qword):
        r"""
        Helper function to get the colored Jones polynomial.

        INPUT:

        - ``N`` -- integer; the number of colors
        - ``qword`` -- a right quantum word (possibly in unreduced form)

        EXAMPLES::

            sage: b = BraidGroup(2)([1,1,1])
            sage: db = b.deformed_burau_matrix()[1:,1:]; db
            [cp_0*ap_1*bp_2]
            sage: b._colored_jones_sum(2, db[0,0])
            1 + q - q^2
            sage: b._colored_jones_sum(3, db[0,0])
            1 + q^2 - q^5 - q^6 + q^7
            sage: b._colored_jones_sum(4, db[0,0])
            1 + q^3 - q^8 - q^10 + q^13 + q^14 - q^15
        """
        rqword = RightQuantumWord(qword).reduced_word()
        alg = qword.parent()
        result = alg.base_ring().one()
        current_word = alg.one()
        # This seemingly infinite sum is always finite if the qword comes
        # from a sum of quantum determinants; because at some point
        # the break condition will become true.
        while True:
            current_word *= rqword
            new_rqw = RightQuantumWord(current_word)
            current_word = new_rqw.reduced_word()
            new_eps = new_rqw.eps(N)
            if not new_eps:
                break
            result += new_eps
        return result

    def colored_jones_polynomial(self, N, variab=None, try_inverse=True):
        r"""
        Return the colored Jones polynomial of the trace closure of the braid.

        INPUT:

        - ``N`` -- integer; the number of colors
        - ``variab`` -- (default: `q`) the variable in the resulting
          Laurent polynomial
        - ``try_inverse`` -- boolean (default: ``True``); if ``True``,
          attempt a faster calculation by using the inverse of the braid

        ALGORITHM:

        The algorithm used is described in [HL2018]_. We follow their
        notation, but work in a suitable free algebra over a Laurent
        polynomial ring in one variable to simplify bookkeeping.

        EXAMPLES::

            sage: trefoil = BraidGroup(2)([1,1,1])
            sage: trefoil.colored_jones_polynomial(2)
            q + q^3 - q^4
            sage: trefoil.colored_jones_polynomial(4)
            q^3 + q^7 - q^10 + q^11 - q^13 - q^14 + q^15 - q^17
             + q^19 + q^20 - q^21
            sage: trefoil.inverse().colored_jones_polynomial(4)
            -q^-21 + q^-20 + q^-19 - q^-17 + q^-15 - q^-14 - q^-13
             + q^-11 - q^-10 + q^-7 + q^-3

            sage: figure_eight = BraidGroup(3)([-1, 2, -1, 2])
            sage: figure_eight.colored_jones_polynomial(2)
            q^-2 - q^-1 + 1 - q + q^2
            sage: figure_eight.colored_jones_polynomial(3, 'Q')
            Q^-6 - Q^-5 - Q^-4 + 2*Q^-3 - Q^-2 - Q^-1 + 3 - Q - Q^2
             + 2*Q^3 - Q^4 - Q^5 + Q^6
        """
        if self.components_in_closure() != 1:
            raise ValueError("the number of components must be 1")
        if not hasattr(self, '_cj_with_q'):
            # Move to the __init__ if this class adds one
            self._cj_with_q = {}
        if N in self._cj_with_q:
            cj = self._cj_with_q[N]
            if variab is None:
                return cj
            if isinstance(variab, str):
                variab = LaurentPolynomialRing(ZZ, variab).gen()
            return cj.subs(q=variab)

        db = self.deformed_burau_matrix('q')[1:, 1:]
        q = db.parent().base_ring().base_ring().gen()
        n = db.ncols()
        qword = sum((-1)**(s.cardinality() - 1)
                    * (q * db[list(s), list(s)]).quantum_determinant(q)
                    for s in Subsets(range(n)) if s)
        inverse_shorter = try_inverse
        if try_inverse:
            db_inv = self.inverse().deformed_burau_matrix('q')[1:, 1:]
            q_inv = db_inv.parent().base_ring().base_ring().gen()
            qword_inv = sum((-1)**(s.cardinality() - 1)
                            * (q_inv*db_inv[list(s), list(s)]).quantum_determinant(q_inv)
                            for s in Subsets(range(n)) if s)
            # Check if the inverse has a shorter expression at this point
            inverse_shorter = len(list(qword_inv)) < len(list(qword))
        use_inverse = try_inverse and inverse_shorter
        shorter_qword = qword_inv if use_inverse else qword
        knot = Knot(self.inverse()) if use_inverse else Knot(self)
        cj = (q**((N - 1) * (knot.writhe() - self.strands() + 1) / 2)
              * self._colored_jones_sum(N, shorter_qword))
        self._cj_with_q[N] = cj.subs({q: 1/q}) if use_inverse else cj
        return self.colored_jones_polynomial(N, variab, try_inverse)

    def super_summit_set_element(self):
        r"""
        Return an element of the braid's super summit set and the conjugating
        braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1, 2, 3, -1, 2, 1, 3])
            sage: b.super_summit_set_element()
            (s0*s2*s0*s1*s2*s1*s0, s0^-1*s1^-1*s0^-1*s2^-1*s1^-1*s0^-1*s1*s0*s2*s1*s0)
        """
        to_sss = send_to_sss(self)
        B = self.parent()
        return tuple([B._element_from_libbraiding(b) for b in to_sss])

    def ultra_summit_set_element(self):
        r"""
        Return an element of the braid's ultra summit set and the conjugating
        braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1, 2, 3, -1, 2, -1, 3])
            sage: b.ultra_summit_set_element()
            (s0*s1*s0*s2*s1, s0^-1*s1^-1*s0^-1*s2^-1*s1^-1*s0^-1*s1*s2*s1^2*s0)
        """
        to_uss = send_to_uss(self)
        B = self.parent()
        return tuple([B._element_from_libbraiding(b) for b in to_uss])

    def sliding_circuits_element(self) -> tuple:
        r"""
        Return an element of the braid's sliding circuits, and the conjugating
        braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1, 2, 3, -1, 2, -1, 3])
            sage: b.sliding_circuits_element()
            (s0*s1*s0*s2*s1, s0^2*s1*s2)
        """
        to_sc = send_to_sc(self)
        B = self.parent()
        return tuple([B._element_from_libbraiding(b) for b in to_sc])

    def trajectory(self) -> list:
        r"""
        Return the braid's trajectory.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1, 2, 3, -1, 2, -1, 3])
            sage: b.trajectory()
            [s0^-1*s1^-1*s0^-1*s2^-1*s1^-1*s2*s0*s1*s2*s1*s0^2*s1*s2^2,
             s0*s1*s2^3,
             s0*s1*s2*s1^2,
             s0*s1*s0*s2*s1]
        """
        traj = trajectory(self)
        B = self.parent()
        return [B._element_from_libbraiding(b) for b in traj]

    def cyclic_slidings(self) -> list:
        r"""
        Return the braid's cyclic slidings.

        OUTPUT: The braid's cyclic slidings. Each cyclic sliding is a list of braids.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1, 2, 3, -1, 2, 1])
            sage: b.cyclic_slidings()
            [[s0*s2*s1*s0*s1*s2, s0*s1*s2*s1*s0^2, s1*s0*s2^2*s1*s0],
             [s0*s1*s2*s1^2*s0, s0*s1*s2*s1*s0*s2, s1*s0*s2*s0*s1*s2]]
        """
        cs = cyclic_slidings(self)
        B = self.parent()
        return [[B._element_from_libbraiding(b) for b in t] for t in cs]


class RightQuantumWord:
    """
    A right quantum word as in Definition 4.1 of [HL2018]_.

    INPUT:

    - ``words`` -- an element in a suitable free algebra over a Laurent
      polynomial ring in one variable; this input does not need to be in
      reduced form, but the monomials for the input can come in any order

    EXAMPLES::

        sage: from sage.groups.braid import RightQuantumWord
        sage: fig_8 = BraidGroup(3)([-1, 2, -1, 2])
        sage: (
        ....:  bp_1, cp_1, ap_1,
        ....:  bp_3, cp_3, ap_3,
        ....:  bm_0, cm_0, am_0,
        ....:  bm_2, cm_2, am_2
        ....: ) = fig_8.deformed_burau_matrix().parent().base_ring().gens()
        sage: q = bp_1.base_ring().gen()
        sage: RightQuantumWord(ap_1*cp_1 + q**3*bm_2*bp_1*am_0*cm_0)
        The right quantum word represented by
         q*cp_1*ap_1 + q^2*bp_1*cm_0*am_0*bm_2
         reduced from ap_1*cp_1 + q^3*bm_2*bp_1*am_0*cm_0
    """
    def __init__(self, words):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.groups.braid import RightQuantumWord
            sage: fig_8 = BraidGroup(3)([-1, 2, -1, 2])
            sage: (
            ....:  bp_1, cp_1, ap_1,
            ....:  bp_3, cp_3, ap_3,
            ....:  bm_0, cm_0, am_0,
            ....:  bm_2, cm_2, am_2
            ....: ) = fig_8.deformed_burau_matrix().parent().base_ring().gens()
            sage: q = bp_1.base_ring().gen()
            sage: Q = RightQuantumWord(ap_1*cp_1 + q**3*bm_2*bp_1*am_0*cm_0)
            sage: TestSuite(Q).run(skip='_test_pickling')
        """
        self._algebra = words.parent()
        self.q = self._algebra.base_ring().gen()
        self.iq = ~self.q
        self.R = self._algebra.base_ring()
        self._unreduced_words = words
        self._gens = self._algebra._indices.gens()
        self._minus_begin = min((i for i, gen in enumerate(self._gens) if 'm' in str(gen)),
                                default=len(self._gens)) // 3
        split = ((g, str(g), i) for i, g in enumerate(self._gens))
        self._recognize = {g: (s[0], s[1] == 'm', 3 * (i // 3))
                           for g, s, i in split}

    @lazy_attribute
    def tuples(self):
        r"""
        Get a representation of the right quantum word as a ``dict``, with
        keys monomials in the free algebra represented as tuples and
        values in elements the Laurent polynomial ring in one variable.

        This is in the reduced form as outlined in Definition 4.1
        of [HL2018]_.

        OUTPUT: a dict of tuples of ints corresponding to the exponents in the
        generators with values in the algebra's base ring

        EXAMPLES::

            sage: from sage.groups.braid import RightQuantumWord
            sage: fig_8 = BraidGroup(3)([-1, 2, -1, 2])
            sage: (
            ....:  bp_1, cp_1, ap_1,
            ....:  bp_3, cp_3, ap_3,
            ....:  bm_0, cm_0, am_0,
            ....:  bm_2, cm_2, am_2
            ....: ) = fig_8.deformed_burau_matrix().parent().base_ring().gens()
            sage: q = bp_1.base_ring().gen()
            sage: qw = RightQuantumWord(ap_1*cp_1 +
            ....:                       q**3*bm_2*bp_1*am_0*cm_0)
            sage: for key, value in qw.tuples.items():
            ....:     print(key, value)
            ....:
            (0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0) q
            (1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0) q^2
        """
        from collections import defaultdict
        ret = defaultdict(self.R)
        convert = self._recognize
        q = self.q
        iq = self.iq
        for unreduced_monom, q_power in list(self._unreduced_words):
            ret_tuple = [0] * len(self._gens)
            for gen, exp in unreduced_monom:
                letter, is_minus, index = convert[gen]
                # This uses the relations in equations (4.1) and (4.2)
                # of [HL2018]_.
                if letter == 'a':  # is_a
                    ret_tuple[index + 2] += exp
                elif letter == 'b':  # is_b
                    j, k = ret_tuple[index + 1:index + 3]
                    ret_tuple[index] += exp
                    q_power *= q**(2*(k*exp + j*exp)) if is_minus else iq**(2*j*exp)
                else:  # is_c
                    k = ret_tuple[index + 2]
                    ret_tuple[index + 1] += exp
                    q_power *= iq**(k*exp) if is_minus else q**(k*exp)
            ret[tuple(ret_tuple)] += q_power
        return ret

    def reduced_word(self):
        r"""
        Return the (reduced) right quantum word.

        OUTPUT: an element in the free algebra

        EXAMPLES::

            sage: from sage.groups.braid import RightQuantumWord
            sage: fig_8 = BraidGroup(3)([-1, 2, -1, 2])
            sage: (
            ....:  bp_1, cp_1, ap_1,
            ....:  bp_3, cp_3, ap_3,
            ....:  bm_0, cm_0, am_0,
            ....:  bm_2, cm_2, am_2
            ....: ) = fig_8.deformed_burau_matrix().parent().base_ring().gens()
            sage: q = bp_1.base_ring().gen()
            sage: qw = RightQuantumWord(ap_1*cp_1 +
            ....:                            q**3*bm_2*bp_1*am_0*cm_0)
            sage: qw.reduced_word()
            q*cp_1*ap_1 + q^2*bp_1*cm_0*am_0*bm_2

        TESTS:

        Testing the equations (4.1) and (4.2) in [HL2018]_::

            sage: RightQuantumWord(ap_3*bp_3).reduced_word()
            bp_3*ap_3
            sage: RightQuantumWord(ap_3*cp_3).reduced_word()
            q*cp_3*ap_3
            sage: RightQuantumWord(cp_3*bp_3).reduced_word()
            (q^-2)*bp_3*cp_3
            sage: RightQuantumWord(am_2*bm_2).reduced_word()
            q^2*bm_2*am_2
            sage: RightQuantumWord(am_2*cm_2).reduced_word()
            (q^-1)*cm_2*am_2
            sage: RightQuantumWord(cm_2*bm_2).reduced_word()
            q^2*bm_2*cm_2

        .. TODO::

            Parallelize this function, calculating all summands in the sum
            in parallel.
        """
        M = self._algebra._indices

        def tuple_to_word(q_tuple):
            return M.prod(self._gens[i]**exp
                          for i, exp in enumerate(q_tuple))

        ret = {tuple_to_word(q_tuple): q_factor
               for q_tuple, q_factor in self.tuples.items() if q_factor}
        return self._algebra._from_dict(ret, remove_zeros=False)

    def eps(self, N):
        r"""
        Evaluate the map `\mathcal{E}_N` for a braid.

        INPUT:

        - ``N`` -- integer; the number of colors

        EXAMPLES::

            sage: from sage.groups.braid import RightQuantumWord
            sage: B = BraidGroup(3)
            sage: b = B([1,-2,1,2])
            sage: db = b.deformed_burau_matrix()
            sage: q = db.parent().base_ring().base_ring().gen()
            sage: (bp_0, cp_0, ap_0,
            ....:  bp_2, cp_2, ap_2,
            ....:  bp_3, cp_3, ap_3,
            ....:  bm_1, cm_1, am_1) = db.parent().base_ring().gens()
            sage: rqw = RightQuantumWord(
            ....:    q^3*bp_2*bp_0*ap_0 + q*ap_3*bm_1*am_1*bp_0)
            sage: rqw.eps(3)
            -q^-1 + 2*q - q^5
            sage: rqw.eps(2)
            -1 + 2*q - q^2 + q^3 - q^4

        TESTS::

            sage: rqw.eps(1)
            0

        .. TODO::

            Parallelize this function, calculating all summands in the sum
            in parallel.
        """
        def eps_monom(q_tuple):
            r"""
            Evaluate the map `\mathcal{E}_N` for a single mononial.
            """
            q = self.q
            ret_q = q**sum((N - 1 - q_tuple[3*i + 2])*q_tuple[3*i + 1]
                           for i in range(self._minus_begin))
            ret_q *= q**sum((N - 1)*(-q_tuple[rj])
                            for rj in range(self._minus_begin * 3 + 1,
                                            len(q_tuple), 3))
            ret_q *= prod(prod(1 - q**(N - 1 - q_tuple[3*i + 1] - h)
                               for h in range(q_tuple[3*i + 2]))
                          for i in range(self._minus_begin))
            ret_q *= prod(prod(1 - q**(q_tuple[3*j + 1] + k + 1 - N)
                               for k in range(q_tuple[3*j + 2]))
                          for j in range(self._minus_begin,
                                         len(q_tuple)//3))
            return ret_q

        return sum(q_factor * eps_monom(q_tuple)
                   for q_tuple, q_factor in self.tuples.items())

    def __repr__(self) -> str:
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.groups.braid import RightQuantumWord
            sage: b = BraidGroup(3)([1,2,-1,2,-1])
            sage: db = b.deformed_burau_matrix(); db[2,2]
            cp_1*am_2*bp_3
            sage: RightQuantumWord(db[2,2])
            The right quantum word represented by cp_1*bp_3*am_2 reduced from
             cp_1*am_2*bp_3
        """
        return ('The right quantum word represented by '
                + f'{str(self.reduced_word())} reduced from '
                + f'{str(self._unreduced_words)}')


class BraidGroup_class(FiniteTypeArtinGroup):
    """
    The braid group on `n` strands.

    EXAMPLES::

        sage: B1 = BraidGroup(5)
        sage: B1
        Braid group on 5 strands
        sage: B2 = BraidGroup(3)
        sage: B1==B2
        False
        sage: B2 is BraidGroup(3)
        True
    """
    Element = Braid

    def __init__(self, names) -> None:
        """
        Python constructor.

        INPUT:

        - ``names`` -- tuple of strings; the names of the generators

        TESTS::

            sage: B1 = BraidGroup(5) # indirect doctest
            sage: B1
            Braid group on 5 strands
            sage: TestSuite(B1).run()
            sage: B1.category()
            Category of infinite groups

        Check that :issue:`14081` is fixed::

            sage: BraidGroup(2)
            Braid group on 2 strands
            sage: BraidGroup(('a',))
            Braid group on 2 strands

        Check that :issue:`15505` is fixed::

            sage: B = BraidGroup(4)
            sage: B.relations()
            (s0*s1*s0*s1^-1*s0^-1*s1^-1, s0*s2*s0^-1*s2^-1, s1*s2*s1*s2^-1*s1^-1*s2^-1)
            sage: B = BraidGroup('a,b,c,d,e,f')
            sage: B.relations()
            (a*b*a*b^-1*a^-1*b^-1,
             a*c*a^-1*c^-1,
             a*d*a^-1*d^-1,
             a*e*a^-1*e^-1,
             a*f*a^-1*f^-1,
             b*c*b*c^-1*b^-1*c^-1,
             b*d*b^-1*d^-1,
             b*e*b^-1*e^-1,
             b*f*b^-1*f^-1,
             c*d*c*d^-1*c^-1*d^-1,
             c*e*c^-1*e^-1,
             c*f*c^-1*f^-1,
             d*e*d*e^-1*d^-1*e^-1,
             d*f*d^-1*f^-1,
             e*f*e*f^-1*e^-1*f^-1)

            sage: BraidGroup([])
            Traceback (most recent call last):
            ...
            ValueError: the number of strands must be at least 2
        """
        n = len(names)
        # n is the number of generators, not the number of strands
        # see issue 14081
        if n < 1:
            raise ValueError("the number of strands must be at least 2")
        free_group = FreeGroup(names)
        rels = []
        for i in range(1, n):
            rels.append(free_group([i, i + 1, i, -i - 1, -i, -i - 1]))
            rels.extend(free_group([i, j, -i, -j])
                        for j in range(i + 2, n + 1))
        cat = Groups().Infinite()
        FinitelyPresentedGroup.__init__(self, free_group, tuple(rels),
                                        category=cat)
        self._nstrands = n + 1
        self._coxeter_group = Permutations(self._nstrands)

        # For caching TL_representation()
        self._TL_representation_dict = {}

        # For caching hermitian form of unitary Burau representation
        self._hermitian_form = None

    def __reduce__(self) -> tuple:
        """
        TESTS::

            sage: B = BraidGroup(3)
            sage: B.__reduce__()
            (<class 'sage.groups.braid.BraidGroup_class'>, (('s0', 's1'),))
            sage: B = BraidGroup(3, 'sigma')
            sage: B.__reduce__()
            (<class 'sage.groups.braid.BraidGroup_class'>, (('sigma0', 'sigma1'),))
        """
        return (BraidGroup_class, (self.variable_names(), ))

    def _repr_(self) -> str:
        """
        Return a string representation.

        TESTS::

            sage: B1 = BraidGroup(5)
            sage: B1  # indirect doctest
            Braid group on 5 strands
        """
        return "Braid group on %s strands" % self._nstrands

    def cardinality(self):
        """
        Return the number of group elements.

        OUTPUT: Infinity

        TESTS::

            sage: B1 = BraidGroup(5)
            sage: B1.cardinality()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    order = cardinality

    def as_permutation_group(self):
        """
        Return an isomorphic permutation group.

        OUTPUT: this raises a :exc:`ValueError` error since braid groups
        are infinite

        TESTS::

            sage: B = BraidGroup(4, 'g')
            sage: B.as_permutation_group()
            Traceback (most recent call last):
            ...
            ValueError: the group is infinite
        """
        raise ValueError("the group is infinite")

    def strands(self):
        """
        Return the number of strands.

        OUTPUT: integer

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.strands()
            4
        """
        return self._nstrands

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: B = BraidGroup(4)
            sage: B([1, 2, 3]) # indirect doctest
            s0*s1*s2
            sage: p = Permutation([3,1,2,4]); B(p)
            s0*s1
            sage: q = SymmetricGroup(4)((1,2)); B(q)
            s0
        """
        if not isinstance(x, (tuple, list)):
            if isinstance(x, (SymmetricGroupElement, Permutation)):
                x = self._standard_lift_Tietze(x)
        return self.element_class(self, x)

    def an_element(self):
        """
        Return an element of the braid group.

        This is used both for illustration and testing purposes.

        EXAMPLES::

            sage: B = BraidGroup(2)
            sage: B.an_element()
            s
        """
        return self.gen(0)

    def some_elements(self) -> list:
        """
        Return a list of some elements of the braid group.

        This is used both for illustration and testing purposes.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: B.some_elements()
            [s0, s0*s1, (s0*s1)^3]
        """
        elements_list = [self.gen(0)]
        elements_list.append(self(range(1, self.strands())))
        elements_list.append(elements_list[-1]**self.strands())
        return elements_list

    def _standard_lift_Tietze(self, p) -> tuple:
        """
        Helper for :meth:`_standard_lift_Tietze`.

        INPUT:

        - ``p`` -- a permutation

        The standard lift of a permutation is the only braid with
        the following properties:

        - The braid induces the given permutation.

        - The braid is positive (that is, it can be written without
          using the inverses of the generators).

        - Every two strands cross each other at most once.

        OUTPUT: a shortest word that represents the braid, in Tietze list form

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: P = Permutation([5, 3, 1, 2, 4])
            sage: B._standard_lift_Tietze(P)
            (1, 2, 3, 4, 1, 2)
        """
        G = SymmetricGroup(self.strands())
        return tuple(G(p).reduced_word())

    @cached_method
    def _links_gould_representation(self, symbolics=False):
        """
        Compute the representation matrices of the generators of the R-matrix
        representation being attached the quantum superalgebra `sl_q(2|1)`.

        INPUT:

        - ``symbolics`` -- boolean (default: ``False``); if set to ``True`` the
          coefficients will be contained in the symbolic ring. Per default they
          are elements of a quotient ring of a three variate Laurent polynomial
          ring.

        OUTPUT:

        A tuple of length equal to the number `n` of strands. The first `n-1`
        items are pairs of the representation matrices of the generators and
        their inverses. The last item is the quantum trace operator of the
        `n`-fold tensorproduct of the natural module.

        TESTS::

            sage: B = BraidGroup(3)
            sage: g1, g2, mu3 = B._links_gould_representation()
            sage: R1, R1I = g1
            sage: R2, R2I = g2
            sage: R1*R2*R1 == R2*R1*R2
            True
        """
        from sage.matrix.constructor import matrix
        n = self.strands()
        d = 4  # dimension of the natural module
        from sage.matrix.special import diagonal_matrix
        if symbolics:
            from sage.misc.functional import sqrt
            from sage.symbolic.ring import SR as BR
            t0, t1 = BR.var('t0, t1')
            s0 = sqrt(t0)
            s1 = sqrt(t1)
            Y = sqrt(-(t0 - 1)*(t1 - 1))
            sparse = False
        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            LR = LaurentPolynomialRing(ZZ, 's0r, s1r')
            s0r, s1r = LR.gens()
            PR = PolynomialRing(LR, 'Yr')
            Yr = PR.gens()
            pqr = Yr**2 + (s0r**2 - 1) * (s1r**2 - 1)
            BR = PR.quotient_ring(pqr)
            s0 = BR(s0r)
            s1 = BR(s1r)
            t0 = BR(s0r**2)
            t1 = BR(s1r**2)
            Y = BR(Yr)
            sparse = True

        # degree one quantum trace operator as defined in I. Marin
        mu = diagonal_matrix([t0**(-1), - t1, - t0**(-1), t1])
        if n == 2:
            # R-Matrix taken from I. Marin
            R = matrix(BR, {(0, 0): t0, (1, 4): s0, (2, 8): s0, (3, 12): 1,
                            (4, 1): s0, (4, 4): t0 - 1, (5, 5): -1, (6, 6): t0*t1 - 1,
                            (6, 9): -s0*s1, (6, 12): -Y*s0*s1, (7, 13): s1, (8, 2): s0,
                            (8, 8): t0 - 1, (9, 6): -s0*s1, (9, 12): Y, (10, 10): -1,
                            (11, 14): s1, (12, 3): 1, (12, 6): -Y*s0*s1, (12, 9): Y,
                            (12, 12): -(t0 - 1)*(t1 - 1), (13, 7): s1, (13, 13): t1 - 1,
                            (14, 11): s1, (14, 14): t1 - 1, (15, 15): t1}, sparse=sparse)
            RI = (~t0 + ~t1)*(1 + R) - ~t0*~t1*(R + R**2) - 1

            # quantum trace operator on two fold tensor space
            E = mu.parent().one()
            mu2 = E.tensor_product(mu)
            return ([R, RI], mu2)

        from sage.matrix.matrix_space import MatrixSpace
        Ed = MatrixSpace(BR, d, d, sparse=sparse).one()
        BGsub = BraidGroup(n-1)
        if n > 3:
            BG2 = BraidGroup(2)
        else:
            BG2 = BGsub
        g1 = list(BG2._links_gould_representation(symbolics=symbolics))
        mu2 = g1.pop()
        R, RI = g1[0]
        lg_sub = list(BGsub._links_gould_representation(symbolics=symbolics))
        musub = lg_sub.pop()

        # extend former generators
        lg = [(g.tensor_product(Ed), gi.tensor_product(Ed)) for g, gi in lg_sub]
        En = MatrixSpace(BR, d**(n-2), d**(n-2), sparse=sparse).one()

        # define new  generator
        gn = En.tensor_product(R)
        gni = En.tensor_product(RI)

        # quantum trace operator on n fold tensor space
        mun = musub.tensor_product(mu)
        return tuple(lg + [(gn, gni), mun])

    @cached_method
    def _LKB_matrix_(self, braid, variab):
        """
        Compute the Lawrence-Krammer-Bigelow representation matrix.

        The variables of the matrix must be given. This actual
        computation is done in this helper method for caching
        purposes.

        INPUT:

        - ``braid`` -- tuple of integers; the Tietze list of the
          braid

        - ``variab`` -- string. The names of the variables that will
          appear in the matrix. They must be given as a string,
          separated by a comma

        OUTPUT: the LKB matrix of the braid, with respect to the variables

        TESTS::

            sage: B = BraidGroup(3)
            sage: B._LKB_matrix_((2, 1, 2), 'x, y')
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: B._LKB_matrix_((1, 2, 1), 'x, y')
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: B._LKB_matrix_((-1, -2, -1, 2, 1, 2), 'x, y')
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        n = self.strands()
        if len(braid) > 1:
            A = self._LKB_matrix_(braid[:1], variab)
            for i in braid[1:]:
                A = A*self._LKB_matrix_((i,), variab)
            return A
        n2 = [set(X) for X in combinations(range(n), 2)]
        R = LaurentPolynomialRing(ZZ, variab)
        q = R.gens()[0]
        t = R.gens()[1]
        if not braid:
            return identity_matrix(R, len(n2), sparse=True)
        A = matrix(R, len(n2), sparse=True)
        if braid[0] > 0:
            i = braid[0] - 1
            for m in range(len(n2)):
                j = min(n2[m])
                k = max(n2[m])
                if i == j-1:
                    A[n2.index(Set([i, k])), m] = q
                    A[n2.index(Set([i, j])), m] = q*q-q
                    A[n2.index(Set([j, k])), m] = 1-q
                elif i == j and not j == k-1:
                    A[n2.index(Set([j, k])), m] = 0
                    A[n2.index(Set([j+1, k])), m] = 1
                elif k-1 == i and not k-1 == j:
                    A[n2.index(Set([j, i])), m] = q
                    A[n2.index(Set([j, k])), m] = 1-q
                    A[n2.index(Set([i, k])), m] = (1-q)*q*t
                elif i == k:
                    A[n2.index(Set([j, k])), m] = 0
                    A[n2.index(Set([j, k+1])), m] = 1
                elif i == j and j == k-1:
                    A[n2.index(Set([j, k])), m] = -t*q*q
                else:
                    A[n2.index(Set([j, k])), m] = 1
            return A
        else:
            i = -braid[0]-1
            for m in range(len(n2)):
                j = min(n2[m])
                k = max(n2[m])
                if i == j-1:
                    A[n2.index(Set([j-1, k])), m] = 1
                elif i == j and not j == k-1:
                    A[n2.index(Set([j+1, k])), m] = q**(-1)
                    A[n2.index(Set([j, k])), m] = 1-q**(-1)
                    A[n2.index(Set([j, j+1])), m] = t**(-1)*q**(-1)-t**(-1)*q**(-2)
                elif k-1 == i and not k-1 == j:
                    A[n2.index(Set([j, k-1])), m] = 1
                elif i == k:
                    A[n2.index(Set([j, k+1])), m] = q**(-1)
                    A[n2.index(Set([j, k])), m] = 1-q**(-1)
                    A[n2.index(Set([k, k+1])), m] = -q**(-1)+q**(-2)
                elif i == j and j == k-1:
                    A[n2.index(Set([j, k])), m] = -t**(-1)*q**(-2)
                else:
                    A[n2.index(Set([j, k])), m] = 1
            return A

    def dimension_of_TL_space(self, drain_size):
        """
        Return the dimension of a particular Temperley--Lieb representation
        summand of ``self``.

        Following the notation of :meth:`TL_basis_with_drain`, the summand
        is the one corresponding to the number of drains being fixed to be
        ``drain_size``.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands
          (both inclusive)

        EXAMPLES:

        Calculation of the dimension of the representation of `B_8`
        corresponding to having `2` drains::

            sage: B = BraidGroup(8)
            sage: B.dimension_of_TL_space(2)
            28

        The direct sum of endomorphism spaces of these vector spaces make up
        the entire Temperley--Lieb algebra::

            sage: import sage.combinat.diagram_algebras as da                           # needs sage.combinat
            sage: B = BraidGroup(6)
            sage: dimensions = [B.dimension_of_TL_space(d)**2 for d in [0, 2, 4, 6]]
            sage: total_dim = sum(dimensions)
            sage: total_dim == len(list(da.temperley_lieb_diagrams(6)))         # long time, needs sage.combinat
            True
        """
        n = self.strands()
        if drain_size > n:
            raise ValueError("number of drains must not exceed number of strands")
        if (n + drain_size) % 2 == 1:
            raise ValueError("parity of strands and drains must agree")

        m = (n - drain_size) // 2
        return Integer(n-1).binomial(m) - Integer(n-1).binomial(m - 2)

    def TL_basis_with_drain(self, drain_size):
        """
        Return a basis of a summand of the Temperley--Lieb--Jones
        representation of ``self``.

        The basis elements are given by non-intersecting pairings of `n+d`
        points in a square with `n` points marked 'on the top' and `d` points
        'on the bottom' so that every bottom point is paired with a top point.
        Here, `n` is the number of strands of the braid group, and `d` is
        specified by ``drain_size``.

        A basis element is specified as a list of integers obtained by
        considering the pairings as obtained as the 'highest term' of
        trivalent trees marked by Jones--Wenzl projectors (see e.g. [Wan2010]_).
        In practice, this is a list of nonnegative integers whose first
        element is ``drain_size``, whose last element is `0`, and satisfying
        that consecutive integers have difference `1`. Moreover, the length
        of each basis element is `n + 1`.

        Given these rules, the list of lists is constructed recursively
        in the natural way.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands
          (both inclusive)

        OUTPUT: list of basis elements, each of which is a list of integers

        EXAMPLES:

        We calculate the basis for the appropriate vector space for `B_5` when
        `d = 3`::

            sage: B = BraidGroup(5)
            sage: B.TL_basis_with_drain(3)
            [[3, 4, 3, 2, 1, 0],
             [3, 2, 3, 2, 1, 0],
             [3, 2, 1, 2, 1, 0],
             [3, 2, 1, 0, 1, 0]]

        The number of basis elements hopefully corresponds to the general
        formula for the dimension of the representation spaces::

            sage: B = BraidGroup(10)
            sage: d = 2
            sage: B.dimension_of_TL_space(d) == len(B.TL_basis_with_drain(d))
            True
        """
        def fill_out_forest(forest, treesize):
            # The basis elements are built recursively using this function,
            # which takes a collection of partial basis elements, given in
            # terms of trivalent trees (i.e. a 'forest') and extends each of
            # the trees by one branch.
            if not forest:
                raise ValueError("forest has to start with a tree")
            if forest[0][0] + treesize % 2 == 0:
                raise ValueError("parity mismatch in forest creation")
            # Loop over all trees
            newforest = list(forest)
            for tree in forest:
                if len(tree) < treesize:
                    newtreeup = list(tree)
                    newtreedown = list(tree)
                    newforest.remove(tree)  # Cut down the original tree
                    # Add two greater trees, admissibly. We need to check two
                    # things to ensure that the tree will eventually define a
                    # basis elements: that its 'colour' is not too large, and
                    # that it is positive.
                    if tree[-1] < treesize - len(tree) + 1:
                        newtreeup.append(tree[-1] + 1)
                        newforest.append(newtreeup)
                    if tree[-1] > 0:
                        newtreedown.append(tree[-1] - 1)
                        newforest.append(newtreedown)
            # Are we there yet?
            if len(newforest[0]) == treesize:
                return newforest
            return fill_out_forest(newforest, treesize)

        n = self.strands()
        if drain_size > n:
            raise ValueError("number of drains must not exceed number of strands")
        if (n + drain_size) % 2 == 1:
            raise ValueError("parity of strands and drains must agree")

        # We can now start the process: all we know is that our basis elements
        # have a drain size of d, so we use fill_out_forest to build all basis
        # elements out of this
        basis = [[drain_size]]
        forest = fill_out_forest(basis, n-1)
        for tree in forest:
            tree.extend([1, 0])
        return forest

    @cached_method
    def _TL_action(self, drain_size):
        """
        Return a matrix representing the action of cups and caps on
        Temperley--Lieb diagrams corresponding to ``self``.

        The action space is the space of non-crossing diagrams of `n+d`
        points, where `n` is the number of strands, and `d` is specified by
        ``drain_size``. As in :meth:`TL_basis_with_drain`, we put certain
        constraints on the diagrams.

        We essentially calculate the action of the TL-algebra generators
        `e_i` on the algebra itself: the action of `e_i` on one of our basis
        diagrams is itself a basis diagram, and ``auxmat`` will store the
        index of this new basis diagram.

        In some cases, the new diagram will connect two bottom points which
        we explicitly disallow (as such a diagram is not one of our basis
        elements). In this case, the corresponding ``auxmat`` entry will
        be `-1`.

        This is used in :meth:`TL_representation` and could be included
        entirely in that method. They are split for purposes of caching.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands
          (both inclusive)

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B._TL_action(2)
            [ 0  0 -1]
            [ 1  1  1]
            [-1  2  2]
            sage: B._TL_action(0)
            [1 1]
            [0 0]
            [1 1]
            sage: B = BraidGroup(6)
            sage: B._TL_action(2)
            [ 1  1  2  3  1  2  3 -1 -1]
            [ 0  0  5  6  5  5  6  5  6]
            [ 1  1  1 -1  4  4  8  8  8]
            [ 5  2  2  2  5  5  5  7  7]
            [-1 -1  3  3  8  6  6  8  8]
        """
        n = self.strands()
        basis = self.TL_basis_with_drain(drain_size)
        auxmat = matrix(n - 1, len(basis))
        for i in range(1, n):  # For each of the e_i
            for v, tree in enumerate(basis):  # For each basis element
                if tree[i - 1] < tree[i] and tree[i + 1] < tree[i]:
                    # Here, for instance, we've created an unknot.
                    auxmat[i-1, v] = v
                if tree[i-1] > tree[i] and tree[i+1] > tree[i]:
                    newtree = list(tree)
                    newtree[i] += 2
                    auxmat[i-1, v] = basis.index(newtree)
                if tree[i-1] > tree[i] and tree[i+1] < tree[i]:
                    newtree = list(tree)
                    newtree[i-1] -= 2
                    j = 2
                    while newtree[i-j] != newtree[i] and i-j >= 0:
                        newtree[i-j] -= 2
                        j += 1
                    if newtree in basis:
                        auxmat[i-1, v] = basis.index(newtree)
                    else:
                        auxmat[i-1, v] = -1
                if tree[i-1] < tree[i] and tree[i+1] > tree[i]:
                    newtree = list(tree)
                    newtree[i+1] -= 2
                    j = 2
                    while newtree[i+j] != newtree[i] and i+j <= n:
                        newtree[i+j] -= 2
                        j += 1
                    if newtree in basis:
                        auxmat[i-1, v] = basis.index(newtree)
                    else:
                        auxmat[i-1, v] = -1
        return auxmat

    def TL_representation(self, drain_size, variab=None):
        r"""
        Return representation matrices of the Temperley--Lieb--Jones
        representation of standard braid group generators and inverses
        of ``self``.

        The basis is given by non-intersecting pairings of `(n+d)` points,
        where `n` is the number of strands, and `d` is given by
        ``drain_size``, and the pairings satisfy certain rules. See
        :meth:`TL_basis_with_drain()` for details. This basis has
        the useful property that all resulting entries can be regarded as
        Laurent polynomials.

        We use the convention that the eigenvalues of the standard generators
        are `1` and `-A^4`, where `A` is the generator of the Laurent
        polynomial ring.

        When `d = n - 2` and the variables are picked appropriately, the
        resulting representation is equivalent to the reduced Burau
        representation. When `d = n`, the resulting representation is
        trivial and 1-dimensional.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands
          (both inclusive)
        - ``variab`` -- variable (default: ``None``); the variable in the
          entries of the matrices; if ``None``, then use a default variable
          in `\ZZ[A,A^{-1}]`

        OUTPUT: list of matrices corresponding to the representations of each
        of the standard generators and their inverses

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.TL_representation(0)
            [(
              [   1    0]  [    1     0]
              [ A^2 -A^4], [ A^-2 -A^-4]
            ),
             (
              [-A^4  A^2]  [-A^-4  A^-2]
              [   0    1], [    0     1]
            ),
             (
              [   1    0]  [    1     0]
              [ A^2 -A^4], [ A^-2 -A^-4]
            )]
            sage: R.<A> = LaurentPolynomialRing(GF(2))
            sage: B.TL_representation(0, variab=A)
            [(
              [  1   0]  [   1    0]
              [A^2 A^4], [A^-2 A^-4]
            ),
             (
              [A^4 A^2]  [A^-4 A^-2]
              [  0   1], [   0    1]
            ),
             (
              [  1   0]  [   1    0]
              [A^2 A^4], [A^-2 A^-4]
            )]
            sage: B = BraidGroup(8)
            sage: B.TL_representation(8)
            [([1], [1]),
             ([1], [1]),
             ([1], [1]),
             ([1], [1]),
             ([1], [1]),
             ([1], [1]),
             ([1], [1])]
        """
        if variab is None:
            if drain_size in self._TL_representation_dict:
                return self._TL_representation_dict[drain_size]
            R = LaurentPolynomialRing(ZZ, 'A')
            A = R.gens()[0]
        else:
            R = variab.parent()
            A = variab

        n = self.strands()
        auxmat = self._TL_action(drain_size)
        dimension = auxmat.ncols()
        # The action of the sigma_i is given in terms of the actions of the
        # e_i which is what auxmat describes. Our choice of normalization means
        # that \sigma_i acts by the identity + A**2 e_i.
        rep_matrices = []  # The list which will store the actions of sigma_i

        # Store the respective powers
        Ap2 = A**2
        Apm2 = A**(-2)
        Ap4 = -A**4
        Apm4 = -A**(-4)

        for i in range(n-1):  # For each \sigma_{i+1}
            rep_mat_new = identity_matrix(R, dimension, sparse=True)
            rep_mat_new_inv = identity_matrix(R, dimension, sparse=True)
            for v in range(dimension):
                new_mat_entry = auxmat[i, v]
                if new_mat_entry == v:  # Did we create an unknot?
                    rep_mat_new[v, v] = Ap4
                    rep_mat_new_inv[v, v] = Apm4
                elif new_mat_entry >= 0:
                    rep_mat_new[new_mat_entry, v] = Ap2
                    rep_mat_new_inv[new_mat_entry, v] = Apm2
            rep_matrices.append((rep_mat_new, rep_mat_new_inv))

        if variab is None:  # Cache the result in this case
            for mat_pair in rep_matrices:
                mat_pair[0].set_immutable()
                mat_pair[1].set_immutable()
            self._TL_representation_dict[drain_size] = rep_matrices

        return rep_matrices

    def mapping_class_action(self, F):
        """
        Return the action of ``self`` in the free group F as mapping class
        group.

        This action corresponds to the action of the braid over the
        punctured disk, whose fundamental group is the free group on
        as many generators as strands.

        In Sage, this action is the result of multiplying a free group
        element with a braid. So you generally do not have to
        construct this action yourself.

        OUTPUT: a :class:`MappingClassGroupAction`

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: B.inject_variables()
            Defining s0, s1
            sage: F.<a,b,c> = FreeGroup(3)
            sage: A = B.mapping_class_action(F)
            sage: A(a,s0)
            a*b*a^-1
            sage: a * s0    # simpler notation
            a*b*a^-1
        """
        return MappingClassGroupAction(self, F)

    def _get_action_(self, S, op, self_on_left):
        """
        Let the coercion system discover actions of the braid group on free groups. ::

            sage: B.<b0,b1,b2> = BraidGroup()
            sage: F.<f0,f1,f2,f3> = FreeGroup()
            sage: f1 * b1
            f1*f2*f1^-1

            sage: from sage.structure.all import get_coercion_model
            sage: cm = get_coercion_model()
            sage: cm.explain(f1, b1, operator.mul)
            Action discovered.
                Right action by Braid group on 4 strands on Free Group on generators {f0, f1, f2, f3}
            Result lives in Free Group on generators {f0, f1, f2, f3}
            Free Group on generators {f0, f1, f2, f3}
            sage: cm.explain(b1, f1, operator.mul)
            Unknown result parent.
        """
        import operator
        if is_FreeGroup(S) and op == operator.mul and not self_on_left:
            return self.mapping_class_action(S)
        return None

    def _element_from_libbraiding(self, nf):
        """
        Return the element of ``self`` corresponding to the output
        of libbraiding.

        INPUT:

        - ``nf`` -- list of lists, as returned by libbraiding

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: B._element_from_libbraiding([[-2], [2, 1], [1, 2], [2, 1]])
            (s0^-1*s1^-1*s0^-1*s2^-1*s1^-1*s0^-1*s3^-1*s2^-1*s1^-1*s0^-1)^2*s1*s0^2*s1^2*s0
            sage: B._element_from_libbraiding([[0]])
            1
        """
        if len(nf) == 1:
            return self.delta()**nf[0][0]
        return self.delta()**nf[0][0] * prod(self(i) for i in nf[1:])

    def mirror_involution(self):
        r"""
        Return the mirror involution of ``self``.

        This automorphism maps a braid to another one by replacing
        each generator in its word by the inverse. In general this is
        different from the inverse of the braid since the order of the
        generators in the word is not reversed.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: mirr = B.mirror_involution()
            sage: b = B((1,-2,-1,3,2,1))
            sage: bm = mirr(b); bm
            s0^-1*s1*s0*s2^-1*s1^-1*s0^-1
            sage: bm == ~b
            False
            sage: bm.is_conjugated(b)
            False
            sage: bm.is_conjugated(~b)
            True
        """
        gens_mirr = [~g for g in self.gens()]
        return self.hom(gens_mirr, check=False)

    def presentation_two_generators(self, isomorphisms=False):
        r"""
        Construct a finitely presented group isomorphic to ``self`` with only two generators.

        INPUT:

        - ``isomorphism`` -- boolean (default: ``False``); if ``True``, then an isomorphism
          from ``self`` and the isomorphic group and its inverse is also returned

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: B.presentation_two_generators()
            Finitely presented group < x0, x1 | x1^3*x0^-2 >
            sage: B = BraidGroup(4)
            sage: G, hom1, hom2 = B.presentation_two_generators(isomorphisms=True)
            sage: G
            Finitely presented group < x0, x1 | x1^4*x0^-3, x0*x1*x0*x1^-2*x0^-1*x1^3*x0^-1*x1^-2 >
            sage: hom1(B.gen(0))
            x0*x1^-1
            sage: hom1(B.gen(1))
            x1*x0*x1^-2
            sage: hom1(B.gen(2))
            x1^2*x0*x1^-3
            sage: all(hom2(hom1(a)) == a for a in B.gens())
            True
            sage: all(hom2(a) == B.one() for a in G.relations())
            True
        """
        n = self.strands()
        F = FreeGroup(2, "x")
        rel = [n * (2,) + (n - 1) * (-1,)]
        rel += [(1,) + (j - 1) * (2,) + (1,) + j * (-2,) + (-1,) + (j + 1) * (2,) + (-1,) + j * (-2,)
                for j in range(2, n - 1)]
        G = F / rel
        if not isomorphisms:
            return G
        a1 = (1, -2)
        L1 = [j * (2,) + a1 + j * (-2,) for j in range(n - 1)]
        h1 = self.hom(codomain=G, im_gens=[G(a) for a in L1], check=False)
        a2 = tuple(range(1, n))
        L2 = [(1,) + a2, a2]
        h2 = G.hom(codomain=self, im_gens=[self(a) for a in L2], check=False)
        return (G, h1, h2)

    def epimorphisms(self, H) -> list:
        r"""
        Return the epimorphisms from ``self`` to ``H``, up to automorphism of `H` passing
        through the :meth:`two generator presentation
        <presentation_two_generators>` of ``self``.

        INPUT:

        - ``H`` -- another group

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: B.epimorphisms(SymmetricGroup(5))
            [Generic morphism:
            From: Braid group on 5 strands
            To:   Symmetric group of order 5! as a permutation group
            Defn: s0 |--> (1,5)
                  s1 |--> (4,5)
                  s2 |--> (3,4)
                  s3 |--> (2,3)]

        ALGORITHM:

        Uses libgap's GQuotients function.
        """
        G, hom1, hom2 = self.presentation_two_generators(isomorphisms=True)
        from sage.misc.misc_c import prod
        HomSpace = self.Hom(H)
        G0g = libgap(self)
        Gg = libgap(G)
        Hg = libgap(H)
        gquotients = Gg.GQuotients(Hg)
        hom1g = libgap.GroupHomomorphismByImagesNC(G0g, Gg,
                                                   [libgap(hom1(u))
                                                    for u in self.gens()])
        g0quotients = [hom1g * h for h in gquotients]
        res = []
        # the following closure is needed to attach a specific value of quo to
        # each function in the different morphisms
        fmap = lambda tup: (lambda a: H(prod(tup[abs(i)-1]**sign(i)
                                             for i in a.Tietze())))
        for quo in g0quotients:
            tup = tuple(H(quo.ImageElm(i.gap()).sage()) for i in self.gens())
            fhom = GroupMorphismWithGensImages(HomSpace, fmap(tup))
            res.append(fhom)
        return res


def BraidGroup(n=None, names='s'):
    """
    Construct a Braid Group.

    INPUT:

    - ``n`` -- integer or ``None`` (default). The number of
      strands. If not specified the ``names`` are counted and the
      group is assumed to have one more strand than generators.

    - ``names`` -- string or list/tuple/iterable of strings (default:
      ``'x'``); the generator names or name prefix

    EXAMPLES::

        sage: B.<a,b> = BraidGroup();  B
        Braid group on 3 strands
        sage: H = BraidGroup('a, b')
        sage: B is H
        True
        sage: BraidGroup(3)
        Braid group on 3 strands

    The entry can be either a string with the names of the generators,
    or the number of generators and the prefix of the names to be
    given. The default prefix is ``'s'`` ::

        sage: B = BraidGroup(3); B.generators()
        (s0, s1)
        sage: BraidGroup(3, 'g').generators()
        (g0, g1)

    Since the word problem for the braid groups is solvable, their Cayley graph
    can be locally obtained as follows (see :issue:`16059`)::

        sage: def ball(group, radius):
        ....:     ret = set()
        ....:     ret.add(group.one())
        ....:     for length in range(1, radius):
        ....:         for w in Words(alphabet=group.gens(), length=length):
        ....:              ret.add(prod(w))
        ....:     return ret
        sage: B = BraidGroup(4)
        sage: GB = B.cayley_graph(elements=ball(B, 4), generators=B.gens()); GB         # needs sage.combinat sage.graphs
        Digraph on 31 vertices

    Since the braid group has nontrivial relations, this graph contains less
    vertices than the one associated to the free group (which is a tree)::

        sage: F = FreeGroup(3)
        sage: GF = F.cayley_graph(elements=ball(F, 4), generators=F.gens()); GF         # needs sage.combinat sage.graphs
        Digraph on 40 vertices

    TESTS::

        sage: G1 = BraidGroup(3, 'a,b')
        sage: G2 = BraidGroup('a,b')
        sage: G3.<a,b> = BraidGroup()
        sage: G1 is G2, G2 is G3
        (True, True)
    """
    # Support Freegroup('a,b') syntax
    if n is not None:
        try:
            n = Integer(n)-1
        except TypeError:
            names = n
            n = None
    # derive n from counting names
    if n is None:
        if isinstance(names, str):
            n = len(names.split(','))
        else:
            names = list(names)
            n = len(names)
    from sage.structure.category_object import normalize_names
    names = normalize_names(n, names)
    return BraidGroup_class(names)


class MappingClassGroupAction(Action):
    r"""
    The right action of the braid group the free group as the mapping
    class group of the punctured disk.

    That is, this action is the action of the braid over the punctured
    disk, whose fundamental group is the free group on as many
    generators as strands.

    This action is defined as follows:

    .. MATH::

        x_j \cdot \sigma_i=\begin{cases}
        x_{j}\cdot x_{j+1}\cdot {x_j}^{-1} & \text{if $i=j$} \\
        x_{j-1} & \text{if $i=j-1$} \\
        x_{j} & \text{otherwise}
        \end{cases},

    where `\sigma_i` are the generators of the braid group on `n`
    strands, and `x_j` the generators of the free group of rank `n`.

    You should left multiplication of the free group element by the
    braid to compute the action. Alternatively, use the
    :meth:`~sage.groups.braid.BraidGroup_class.mapping_class_action`
    method of the braid group to construct this action.

    EXAMPLES::

        sage: B.<s0,s1,s2> = BraidGroup(4)
        sage: F.<x0,x1,x2,x3> = FreeGroup(4)
        sage: x0 * s1
        x0
        sage: x1 * s1
        x1*x2*x1^-1
        sage: x1^-1 * s1
        x1*x2^-1*x1^-1

        sage: A = B.mapping_class_action(F)
        sage: A
        Right action by Braid group on 4 strands on Free Group
        on generators {x0, x1, x2, x3}
        sage: A(x0, s1)
        x0
        sage: A(x1, s1)
        x1*x2*x1^-1
        sage: A(x1^-1, s1)
        x1*x2^-1*x1^-1
    """
    def __init__(self, G, M) -> None:
        """
        TESTS::

            sage: B = BraidGroup(3)
            sage: G = FreeGroup('a, b, c')
            sage: B.mapping_class_action(G) # indirect doctest
            Right action by Braid group on 3 strands on Free Group
            on generators {a, b, c}
        """
        import operator
        Action.__init__(self, G, M, False, operator.mul)

    def _act_(self, b, x):
        """
        Return the action of ``b`` on ``x``.

        INPUT:

        - ``b`` -- a braid

        - ``x`` -- a free group element

        OUTPUT: a new braid

        TESTS::

            sage: B = BraidGroup(3)
            sage: G = FreeGroup('a, b, c')
            sage: A = B.mapping_class_action(G)
            sage: A(G.0, B.0) # indirect doctest
            a*b*a^-1
            sage: A(G.1, B.0) # indirect doctest
            a
        """
        t = x.Tietze()
        for j in b.Tietze():
            s = []
            for i in t:
                if j == i and i > 0:
                    s += [i, i+1, -i]
                elif j == -i and i < 0:
                    s += [-i, i-1, i]
                elif j == -i and i > 0:
                    s += [i+1]
                elif j == i and i < 0:
                    s += [i-1]
                elif i > 0 and j == i-1:
                    s += [i-1]
                elif i < 0 and j == -i-1:
                    s += [i+1]
                elif i > 0 and -j == i-1:
                    s += [-i, i-1, i]
                elif i < 0 and j == i+1:
                    s += [i, i+1, -i]
                else:
                    s += [i]
            t = s
        return self.codomain()(t)
