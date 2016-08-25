# -*- coding: utf-8 -*-
"""
Braid groups

Braid groups are implemented as a particular case of finitely presented groups,
but with a lot of specific methods for braids.

A braid group can be created by giving the number of strands, and the name of the generators::

    sage: BraidGroup(3)
    Braid group on 3 strands
    sage: BraidGroup(3,'a')
    Braid group on 3 strands
    sage: BraidGroup(3,'a').gens()
    (a0, a1)
    sage: BraidGroup(3,'a,b').gens()
    (a, b)

The elements can be created by operating with the generators, or by passing a list
with the indices of the letters to the group::

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
"""

##############################################################################
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

import six
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.groups.free_group import FreeGroup, is_FreeGroup
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.matrix.constructor import identity_matrix, matrix
from sage.combinat.permutation import Permutation
from sage.categories.action import Action
from sage.sets.set import Set
from sage.groups.finitely_presented import FinitelyPresentedGroup, FinitelyPresentedGroupElement
from sage.misc.package import PackageNotFoundError

class Braid(FinitelyPresentedGroupElement):
    """
    Class that models elements of the braid group.

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
    def _cmp_(self, other):
        """
        Compare ``self`` and ``other``

        TESTS::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1])
            sage: c = B([2, 1, 2])
            sage: b == c #indirect doctest
            True
            sage: b._cmp_(c^(-1))
            -1
            sage: B([])._cmp_(B.one())
            0
        """
        if self.Tietze()==other.Tietze():
            return 0
        nfself = [i.Tietze() for i in self.left_normal_form()]
        nfother = [i.Tietze() for i in other.left_normal_form()]
        return cmp(nfself, nfother)

    __cmp__ = _cmp_

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

    def _latex_(self):
        """
        Return a LaTeX representation

        OUTPUT:

        String. A valid LaTeX math command sequence.

        TESTS::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 3, -1, 2, -3])
            sage: b._latex_()
            '\\sigma_{1}\\sigma_{2}\\sigma_{3}\\sigma_{1}^{-1}\\sigma_{2}\\sigma_{3}^{-1}'
        """
        latexrepr = ''
        for i in self.Tietze():
            if i > 0:
                latexrepr = latexrepr+"\sigma_{%s}" % i
            if i < 0:
                latexrepr = latexrepr+"\sigma_{%s}^{-1}" % (-i)
        return latexrepr

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

    def exponent_sum(self):
        """
        Return the exponent sum of the braid.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: b = B([1, 4, -3, 2])
            sage: b.exponent_sum()
            2
            sage: b = B([])
            sage: b.exponent_sum()
            0
        """
        return sum(s.sign() for s in self.Tietze())

    def components_in_closure(self):
        """
        Return the number of components of the trace closure of the braid.

        OUTPUT:

        Positive integer.

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
        """
        Return the Burau matrix of the braid.

        INPUT:

        - ``var`` -- string (default: ``'t'``); the name of the
          variable in the entries of the matrix
        - ``reduced`` -- boolean (default: ``False``); whether to
          return the reduced or unreduced Burau representation

        OUTPUT:

        The Burau matrix of the braid. It is a matrix whose entries
        are Laurent polynomials in the variable ``var``. If ``reduced``
        is ``True``, return the matrix for the reduced Burau representation
        instead.

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

        REFERENCES:

        - :wikipedia:`Burau_representation`
        """
        R = LaurentPolynomialRing(IntegerRing(), var)
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
            M = identity_matrix(R, n - 1)
            for j in self.Tietze():
                A = identity_matrix(R, n - 1)
                if j > 1:
                    i = j-1
                    A[i-1, i-1] = 1-t
                    A[i, i] = 0
                    A[i, i-1] = 1
                    A[i-1, i] = t
                if j < -1:
                    i = j+1
                    A[-1-i, -1-i] = 0
                    A[-i, -i] = 1-t**(-1)
                    A[-1-i, -i] = 1
                    A[-i, -1-i] = t**(-1)
                if j == 1:
                    for k in range(n - 1):
                        A[k,0] = -t
                if j == -1:
                    A[0,0] = -t**(-1)
                    for k in range(1, n - 1):
                        A[k,0] = -1
                M = M * A
        return M

    def alexander_polynomial(self, var='t', normalized=True):
        r"""
        Return the Alexander polynomial of the closure of the braid.

        INPUT:

        - ``var`` -- string (default: ``'t'``); the name of the
          variable in the entries of the matrix
        - ``normalized`` -- boolean (default: ``True``); whether to
          return the normalized Alexander polynomial

        OUTPUT:

        The Alexander polynomial of the braid closure of the braid.

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
        K, t = LaurentPolynomialRing(IntegerRing(), var).objgen()
        if p == 0:
            return K.zero()
        qn = sum(t ** i for i in range(n))
        p //= qn
        if normalized:
            p *= t ** (-p.degree())
            if p.constant_coefficient() < 0:
                p = -p
        return p

    def permutation(self):
        """
        Return the permutation induced by the braid in its strands.

        OUTPUT:

        A permutation.

        EXAMPLES::

            sage: B.<s0,s1,s2> = BraidGroup()
            sage: b = s0*s1/s2/s1
            sage: b.permutation()
            [4, 1, 3, 2]
            sage: b.permutation().cycle_string()
            '(1,4,2)'
        """
        per = Permutation((()))
        for i in self.Tietze():
            j = abs(i)
            per = per*Permutation(((j, j+1)))
        return per

    def plot(self, color='rainbow', orientation='bottom-top', gap=0.05, aspect_ratio=1, axes=False, **kwds):
        """
        Plot the braid

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

        - ``gap`` -- floating point number (default: 0.05). determines
          the size of the gap left when a strand goes under another.

        - ``aspect_ratio`` -- floating point number (default:
          ``1``). The aspect ratio.

        - ``**kwds`` -- other keyword options that are passed to
          :meth:`~sage.plot.bezier_path` and :meth:`~sage.plot.line`.

        EXAMPLES::

            sage: B = BraidGroup(4, 's')
            sage: b = B([1, 2, 3, 1, 2, 1])
            sage: b.plot()
            Graphics object consisting of 30 graphics primitives
            sage: b.plot(color=["red", "blue", "red", "blue"])
            Graphics object consisting of 30 graphics primitives

            sage: B.<s,t> = BraidGroup(3)
            sage: b = t^-1*s^2
            sage: b.plot(orientation="left-right", color="red")
            Graphics object consisting of 12 graphics primitives
        """
        from sage.plot.bezier_path import bezier_path
        from sage.plot.plot import Graphics, line
        from sage.plot.colors import rainbow
        if orientation=='top-bottom':
            orx = 0
            ory = -1
            nx = 1
            ny = 0
        elif orientation=='left-right':
            orx = 1
            ory = 0
            nx = 0
            ny = -1
        elif orientation=='bottom-top':
            orx = 0
            ory = 1
            nx = 1
            ny = 0
        else:
            raise ValueError('unknown value for "orientation"')
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
        a = Graphics()
        op = gap
        for i, m in enumerate(braid):
            for j in range(n):
                if m==j+1:
                    a += bezier_path([[(j*nx+i*orx, i*ory+j*ny), (j*nx+orx*(i+0.25), j*ny+ory*(i+0.25)),
                                       (nx*(j+0.5)+orx*(i+0.5), ny*(j+0.5)+ory*(i+0.5))],
                                      [(nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)),
                                       (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col[j], **kwds)
                elif m==j:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j-0.5+4*op)+orx*(i+0.5-2*op), ny*(j-0.5+4*op)+ory*(i+0.5-2*op)),
                                       (nx*(j-0.5+2*op)+orx*(i+0.5-op), ny*(j-0.5+2*op)+ory*(i+0.5-op))]],
                                     color=col[j], **kwds)
                    a += bezier_path([[(nx*(j-0.5-2*op)+orx*(i+0.5+op), ny*(j-0.5-2*op)+ory*(i+0.5+op)),
                                       (nx*(j-0.5-4*op)+orx*(i+0.5+2*op), ny*(j-0.5-4*op)+ory*(i+0.5+2*op)),
                                       (nx*(j-1)+orx*(i+0.75), ny*(j-1)+ory*(i+0.75)),
                                       (nx*(j-1)+orx*(i+1), ny*(j-1)+ory*(i+1))]], color=col[j], **kwds)
                    col[j], col[j-1] = col[j-1], col[j]
                elif -m==j+1:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j+0.5-4*op)+orx*(i+0.5-2*op), ny*(j+0.5-4*op)+ory*(i+0.5-2*op)),
                                       (nx*(j+0.5-2*op)+orx*(i+0.5-op), ny*(j+0.5-2*op)+ory*(i+0.5-op))]],
                                     color=col[j], **kwds)
                    a += bezier_path([[(nx*(j+0.5+2*op)+orx*(i+0.5+op), ny*(j+0.5+2*op)+ory*(i+0.5+op)),
                                       (nx*(j+0.5+4*op)+orx*(i+0.5+2*op), ny*(j+0.5+4*op)+ory*(i+0.5+2*op)),
                                       (nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)),
                                       (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col[j], **kwds)
                elif -m==j:
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
        Plots the braid in 3d.

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
            sage: b.plot3d()
            Graphics3d Object
            sage: b.plot3d(color="red")
            Graphics3d Object
            sage: b.plot3d(color=["red", "blue", "red", "blue"])
            Graphics3d Object
        """
        from sage.plot.plot3d.shapes2 import bezier3d
        from sage.plot.colors import rainbow
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
                if m==j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j+0.5, i+0.5)],
                                       [(0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]], color=col[j]))
                elif -m==j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j+0.5, i+0.5)],
                                       [(-0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]], color=col[j]))
                elif m==j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j-0.5, i+0.5)],
                                       [(-0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]], color=col[j]))
                    col[j], col[j-1] = col[j-1], col[j]
                elif -m==j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j-0.5, i+0.5)],
                                       [(0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]], color=col[j]))
                    col[j], col[j-1] = col[j-1], col[j]
                else:
                    b.append(bezier3d([[(0, j, i), (0, j, i+1)]], color=col[j]))
        return sum(b)

    def LKB_matrix(self, variables='x,y'):
        """
        Return the Lawrence-Krammer-Bigelow representation matrix.

        The matrix is expressed in the basis $\{e_{i, j} \mid 1\\leq i
        < j \leq n\}$, where the indices are ordered
        lexicographically.  It is a matrix whose entries are in the
        ring of Laurent polynomials on the given variables.  By
        default, the variables are ``'x'`` and ``'y'``.

        INPUT:

        - ``variables`` -- string (default: ``'x,y'``). A string
          containing the names of the variables, separated by a comma.

        OUTPUT:

        The matrix corresponding to the Lawrence-Krammer-Bigelow representation of the braid.

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

        .. [Bigelow] Bigelow, Stephen J. The Lawrence-Krammer representation.
           :arxiv:`math/0204057v1`
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

        OUTPUT:

        The matrix of the TL representation of the braid.

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
        [Big99]_::

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

        .. [Big99] Stephen J. Bigelow. The Burau representation is
           not faithful for `n = 5`. Geom. Topol., 3:397--404, 1999.
        .. [JonesNotes] Vaughan Jones. The Jones Polynomial.
           https://math.berkeley.edu/~vfr/jones.pdf
        """
        if variab is None:
            R = LaurentPolynomialRing(IntegerRing(), 'A')
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

    def tropical_coordinates(self):
        r"""
        Return the tropical coordinates of ``self`` in the braid group `B_n`.

        OUTPUT:

        - a list of `2n` tropical integers

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

        .. [Dynnikov07] \I. Dynnikov and B. Wiest, On the complexity of braids,
           J. Europ. Math. Soc. 9 (2007)
        .. [Dehornoy] \P. Dehornoy, Le probleme d'isotopie des tresses, in
           lecons de mathematiques d'aujourd'hui vol. 4
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
        T = TropicalSemiring(IntegerRing())
        return [T(_) for _ in coord]

    def markov_trace(self, variab=None, normalized=True):
        """
        Return the Markov trace of the braid.

        The normalization is so that in the underlying braid group
        representation, the eigenvalues of the standard generators of
        the braid group are `1` and `-A^4`.

        INPUT:

        - ``variab`` -- variable (default: ``None``); the variable in the
          resulting polynomial; if ``None``, then use the variable `A`
          in `\ZZ[A,A^{-1}]`

        - ``normalized`` - boolean (default: ``True``); if specified to be
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
            R = LaurentPolynomialRing(IntegerRing(), 'A')
            A = R.gens()[0]
            one = IntegerRing().one()
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
            sage: b.jones_polynomial()
            1

            sage: B = BraidGroup(2)
            sage: b = B([])
            sage: b._jones_polynomial
            -A^-2 - A^2
            sage: b = B([-1, -1, -1])
            sage: b._jones_polynomial
            -A^-16 + A^-12 + A^-4
        """
        trace = self.markov_trace(normalized=False)
        A = trace.parent().gens()[0]
        D = A**2 + A**(-2)
        exp_sum = self.exponent_sum()
        num_comp = self.components_in_closure()
        return (-1)**(num_comp-1) * A**(2*exp_sum) * trace // D

    def jones_polynomial(self, variab=None, skein_normalization=False):
        """
        Return the Jones polynomial of the trace closure of the braid.

        The normalization is so that the unknot has Jones polynomial `1`. If
        ``skein_normalization`` is ``True``, the variable of the result is
        replaced by a itself to the power of `4`, so that the result
        agrees with the conventions of [Lic]_ (which in particular differs
        slightly from the conventions used otherwise in this class), had
        one used the conventional Kauffman bracket variable notation directly.

        If ``variab`` is ``None`` return a polynomial in the variable `A`
        or `t`, depending on the value ``skein_normalization``. In
        particular, if ``skein_normalization`` is ``False``, return the
        result in terms of the variable `t`, also used in [Lic]_.

        INPUT:

        - ``variab`` -- variable (default: ``None``); the variable in the
          resulting polynomial; if unspecified, use either a default variable
          in `ZZ[A,A^{-1}]` or the variable `t` in the symbolic ring

        - ``skein_normalization`` -- boolean (default: ``False``); determines
          the variable of the resulting polynomial

        OUTPUT:

        If ``skein_normalization`` if ``False``, this returns an element
        in the symbolic ring as the Jones polynomial of the closure might
        have fractional powers when the closure of the braid is not a knot.
        Otherwise the result is a Laurant polynomial in ``variab``.

        EXAMPLES:

        The unknot::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: b.jones_polynomial()
            1

        Two unlinked unknots::

            sage: B = BraidGroup(2)
            sage: b = B([])
            sage: b.jones_polynomial()
            -sqrt(t) - 1/sqrt(t)

        The Hopf link::

            sage: B = BraidGroup(2)
            sage: b = B([-1,-1])
            sage: b.jones_polynomial()
            -1/sqrt(t) - 1/t^(5/2)

        Different representations of the trefoil and one of its mirror::

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
            sage: cmp(b11n42.jones_polynomial(), b11n34.jones_polynomial())
            0

        REFERENCES:

        .. [Lic] William B. Raymond Lickorish. An Introduction to Knot Theory,
           volume 175 of Graduate Texts in Mathematics. Springer-Verlag,
           New York, 1997. ISBN 0-387-98254-X
        """
        if skein_normalization:
            if variab is None:
                return self._jones_polynomial
            else:
                return self._jones_polynomial(variab)
        else:
            from sage.symbolic.ring import SR
            from sage.rings.integer_ring import ZZ
            if variab is None:
                variab = 't'
            # We force the result to be in the symbolic ring because of the expand
            return self._jones_polynomial(SR(variab)**(ZZ(1)/ZZ(4))).expand()

    @cached_method
    def left_normal_form(self):
        """
        Return the left normal form of the braid.

        OUTPUT:

        A tuple of braid generators in the left normal form. The first
        element is a power of $\Delta$, and the rest are permutation
        braids.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 3, -1, 2, -3])
            sage: b.left_normal_form()
            (s0^-1*s1^-1*s2^-1*s0^-1*s1^-1*s0^-1, s0*s1*s2*s1*s0, s0*s2*s1)
            sage: c = B([1])
            sage: c.left_normal_form()
            (1, s0)
        """
        lnfp = self._left_normal_form_perm_()
        a = lnfp[0]
        l = lnfp[1:]
        n = self.strands()
        delta = Permutation([n-i for i in range(n)])
        P = self.parent()
        return tuple( [P._permutation_braid(delta) ** a] +
                      [P._permutation_braid(i) for i in l] )

    def _left_normal_form_perm_(self):
        """
        Return the left normal form of the braid, in permutation form.

        OUTPUT:

        A tuple whose first element is the power of $\Delta$, and the
        rest are the permutations corresponding to the simple factors.

        EXAMPLES::

            sage: B = BraidGroup(12)
            sage: B([2, 2, 2, 3, 1, 2, 3, 2, 1, -2])._left_normal_form_perm_()
            (-1,
             [12, 11, 10, 9, 8, 7, 6, 5, 2, 4, 3, 1],
             [4, 1, 3, 2, 5, 6, 7, 8, 9, 10, 11, 12],
             [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12],
             [3, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12],
             [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12])
            sage: C=BraidGroup(6)
            sage: C([2, 3, -4, 2, 3, -5, 1, -2, 3, 4, 1, -2])._left_normal_form_perm_()
            (-2, [3, 5, 4, 2, 6, 1], [1, 6, 3, 5, 2, 4], [5, 6, 2, 4, 1, 3],
             [3, 2, 4, 1, 5, 6], [1, 5, 2, 3, 4, 6])
        """
        n = self.parent().strands()
        delta = 0
        Delta = Permutation([n-i for i in range(n)])
        l = self.Tietze()
        if l==():
            return (0,)
        form = []
        for i in l:
            if i>0:
                form.append(Permutation((i, i+1)))
            else:
                delta = delta+1
                form = [Delta*a*Delta for a in form]
                form.append(Delta*Permutation((-i, -i+1)))
        i = j = 0
        while j<len(form):
            while i<len(form)-j-1:
                e = form[i].idescents(from_zero=False)
                s = form[i + 1].descents(from_zero=False)
                S = set(s).difference(set(e))
                while S:
                    a = list(S)[0]
                    form[i] = form[i] * Permutation((a, a+1))
                    form[i + 1] = Permutation((a, a+1))*form[i+1]
                    e = form[i].idescents(from_zero=False)
                    s = form[i + 1].descents(from_zero=False)
                    S = set(s).difference(set(e))
                if form[i+1].length()==0:
                    form.pop(i+1)
                    i = 0
                else:
                    i += 1
            j += 1
            i = 0
        form = [a for a in form if a.length()]
        while form!=[] and form[0]==Delta:
            form.pop(0)
            delta = delta-1
        return tuple([-delta]+form)

    def right_normal_form(self):
        """
        Return the right normal form of the braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1, -2, 3, 1])
            sage: b.right_normal_form() # optional - libbraiding
            (s1*s0, s0*s2, 1)

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import rightnormalform
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        l = rightnormalform(self)
        B = self.parent()
        return tuple([B(b) for b in l[:-1]] + [B.Delta() ** l[-1][0]])

    def centralizer(self):
        """
        Return a list of generators of the centralizer of the braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([2, 1, 3, 2])
            sage: b.centralizer() # optional - libbraiding
            [s1*s0*s2*s1, s0*s2]

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import centralizer
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        l = centralizer(self)
        B = self.parent()
        return [B._element_from_libbraiding(b) for b in l]

    def super_summit_set(self):
        """
        Return a list with the super summit set of the braid

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1, -2, -2, 1])
            sage: b.super_summit_set() # optional - libbraiding
            [s0^-1*s1^-1*s0^-2*s1^2*s0^2,
            (s0^-1*s1^-1*s0^-1)^2*s1^2*s0^3*s1,
            (s0^-1*s1^-1*s0^-1)^2*s1*s0^3*s1^2,
            s0^-1*s1^-1*s0^-2*s1^-1*s0*s1^3*s0]

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import supersummitset
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        l = supersummitset(self)
        B = self.parent()
        return [B._element_from_libbraiding(b) for b in l]

    def gcd(self, other):
        """
        Return the greatest common divisor of the two braids.

        INPUT:

        - ``other`` -- the other braid with respect with the gcd is computed

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1, -2, -2, 1])
            sage: c = B([1, 2, 1])
            sage: b.gcd(c) # optional - libbraiding
            s0^-1*s1^-1*s0^-2*s1^2*s0
            sage: c.gcd(b) # optional - libbraiding
            s0^-1*s1^-1*s0^-2*s1^2*s0

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import greatestcommondivisor
        except ImportError:
            raise PackageNotFoundError("libbraiding")
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
            sage: b.lcm(c) # optional - libbraiding
            (s0*s1)^2*s0

        .. WARNING::

            This functionality requires the libbraiding package to be installed.
        """
        try:
            from sage.libs.braiding import leastcommonmultiple
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        B = self.parent()
        b = leastcommonmultiple(self, other)
        return B._element_from_libbraiding(b)

    def conjugating_braid(self, other):
        r"""
        Return a conjugating braid, if it exists.

        INPUT:

        - ``other`` -- the other braid to look for conjugating braid

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1])
            sage: b = B([2, 1, 2, 1])
            sage: c = b * a / b
            sage: d = a.conjugating_braid(c) # optional - libbraiding
            sage: d * c / d == a # optional - libbraiding
            True
            sage: d # optional - libbraiding
            s1*s0
            sage: d * a / d == c # optional - libbraiding
            False

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import conjugatingbraid
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        l = conjugatingbraid(self, other)
        if not l:
            return None
        else:
            return self.parent()._element_from_libbraiding(l)

    def is_conjugated(self, other):
        """
        Check if the two braids are conjugated.

        INPUT:

        - ``other`` -- the other breaid to check for conjugacy

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1])
            sage: b = B([2, 1, 2, 1])
            sage: c = b * a / b
            sage: c.is_conjugated(a) # optional - libbraiding
            True
            sage: c.is_conjugated(b) # optional - libbraiding
            False

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import conjugatingbraid
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        l = conjugatingbraid(self, other)
        return bool(l)

    def ultra_summit_set(self):
        """
        Return a list with the orbits of the ultra summit set of ``self``

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: b = B([2, 1, 2, 1])
            sage: b.ultra_summit_set() # optional - libbraiding
            [[s0*s1*s0^2, (s0*s1)^2]]
            sage: a.ultra_summit_set() # optional - libbraiding
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

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import ultrasummitset
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        uss = ultrasummitset(self)
        B = self.parent()
        return [[B._element_from_libbraiding(i) for i in s] for s in uss]

    def thurston_type(self):
        """
        Return the thurston_type of ``self``.

        OUTPUT:

        One of ``'reducible'``, ``'periodic'`` or ``'pseudo-anosov'``.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1])
            sage: b.thurston_type() # optional - libbraiding
            'reducible'
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.thurston_type() # optional - libbraiding
            'pseudo-anosov'
            sage: c = B([2, 1, 2, 1])
            sage: c.thurston_type() # optional - libbraiding
            'periodic'

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import thurston_type
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        return thurston_type(self)

    def is_reducible(self):
        """
        Check weather the braid is reducible.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, -1])
            sage: b.is_reducible() # optional - libbraiding
            True
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.is_reducible() # optional - libbraiding
            False

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        return self.thurston_type() == 'reducible'

    def is_periodic(self):
        """
        Check weather the braid is periodic.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: b = B([2, 1, 2, 1])
            sage: a.is_periodic() # optional - libbraiding
            False
            sage: b.is_periodic() # optional - libbraiding
            True

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        return self.thurston_type() == 'periodic'

    def is_pseudoanosov(self):
        """
        Check if the braid is pseudo-anosov.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: b = B([2, 1, 2, 1])
            sage: a.is_pseudoanosov() # optional - libbraiding
            True
            sage: b.is_pseudoanosov() # optional - libbraiding
            False

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        return self.thurston_type() == 'pseudo-anosov'

    def rigidity(self):
        """
        Return the rigidity of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([2, 1, 2, 1])
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.rigidity() # optional - libbraiding
            6
            sage: b.rigidity() # optional - libbraiding
            0

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import rigidity
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        return Integer(rigidity(self))

    def sliding_circuits(self):
        """
        Return the sliding circuits of the braid.

        OUTPUT:

        A list of sliding circuits. Each sliding circuit is itself
        a list of braids.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: a = B([2, 2, -1, -1, 2, 2])
            sage: a.sliding_circuits() # optional - libbraiding
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
            sage: b.sliding_circuits() # optional - libbraiding
            [[s0*s1*s0^2, (s0*s1)^2]]

        .. WARNING::

            This functionality requires the libbraiding package
            to be installed.
        """
        try:
            from sage.libs.braiding import sliding_circuits
        except ImportError:
            raise PackageNotFoundError("libbraiding")
        slc = sliding_circuits(self)
        B = self.parent()
        return [[B._element_from_libbraiding(i) for i in s] for s in slc]

class BraidGroup_class(FinitelyPresentedGroup):
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

    def __init__(self, names):
        """
        Python constructor.

        INPUT:

        - ``names`` -- a tuple of strings. The names of the
          generators.

        TESTS::

            sage: B1 = BraidGroup(5) # indirect doctest
            sage: B1
            Braid group on 5 strands
            sage: TestSuite(B1).run()


        Check that :trac:`14081` is fixed::

            sage: BraidGroup(2)
            Braid group on 2 strands
            sage: BraidGroup(('a',))
            Braid group on 2 strands

        Check that :trac:`15505` is fixed::

            sage: B=BraidGroup(4)
            sage: B.relations()
            (s0*s1*s0*s1^-1*s0^-1*s1^-1, s0*s2*s0^-1*s2^-1, s1*s2*s1*s2^-1*s1^-1*s2^-1)
            sage: B=BraidGroup('a,b,c,d,e,f')
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
        """
        n = len(names)
        if n<1: #n is the number of generators, not the number of strands (see ticket 14081)
            raise ValueError("the number of strands must be an integer bigger than one")
        free_group = FreeGroup(names)
        rels = []
        for i in range(1, n):
            rels.append(free_group([i, i+1, i, -i-1, -i, -i-1]))
            for j in range(i+2, n+1):
                rels.append(free_group([i, j, -i, -j]))
        FinitelyPresentedGroup.__init__(self, free_group, tuple(rels))
        self._nstrands_ = n+1

        # For caching TL_representation()
        self._TL_representation_dict = {}

    def __reduce__(self):
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

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        TESTS::

            sage: B1 = BraidGroup(5)
            sage: B1 # indirect doctest
            Braid group on 5 strands
        """
        return "Braid group on %s strands"%self._nstrands_

    def cardinality(self):
        """
        Return the number of group elements.

        OUTPUT:

        Infinity.

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

        OUTPUT:

        Raises a ``ValueError`` error since braid groups are infinite.

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

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.strands()
            4
        """
        return self._nstrands_

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: B = BraidGroup(4)
            sage: B([1, 2, 3]) # indirect doctest
            s0*s1*s2
        """
        return self.element_class(self, x)

    def an_element(self):
        """
        Return an element of the braid group.

        This is used both for illustration and testing purposes.

        EXAMPLES::

            sage: B=BraidGroup(2)
            sage: B.an_element()
            s
        """
        return self.gen(0)

    def some_elements(self):
        """
        Return a list of some elements of the braid group.

        This is used both for illustration and testing purposes.

        EXAMPLES::

            sage: B=BraidGroup(3)
            sage: B.some_elements()
            [s0, s0*s1, (s0*s1)^3]
        """
        elements_list = [self.gen(0)]
        elements_list.append(self(range(1,self.strands())))
        elements_list.append(elements_list[-1]**self.strands())
        return elements_list

    def _permutation_braid_Tietze(self, p):
        """
        Helper for :meth:`_permutation_braid`.

        INPUT:

        - ``p`` -- a permutation.

        OUTPUT:

        The lexicographically smallest word that represents the braid,
        in Tietze list form.

        EXAMPLES::

            sage: B=BraidGroup(5)
            sage: P=Permutation([5, 3, 1, 2, 4])
            sage: B._permutation_braid_Tietze(P)
            (1, 2, 1, 3, 2, 4)
        """
        if p.length() == 0:
            return ()
        pl = p
        l = []
        while pl.length()>0:
            i = 1
            while i<max(pl):
                if pl(i)>pl(i+1):
                    l.append(i)
                    pl = Permutation([(i, i+1)])*pl
                    i = 1
                else:
                    i = i+1
        return tuple(l)

    @cached_method
    def _permutation_braid(self, p):
        """
        Return the braid that corresponds to the given permutation.

        It is the only braid with the following properties:

        - The braid induces the given permutation.

        - The braid is positive (that is, it can be writen without using the inverses of the generators).

        - Every two strands cross each other at most once.

        INPUT:

        - ``p`` -- a permutation.

        OUTPUT:

        The braid that corresponds to the permutation.

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: P = Permutation([5, 3, 1, 2, 4])
            sage: B._permutation_braid(P)
            s0*s1*s0*s2*s1*s3
        """
        return self(self._permutation_braid_Tietze(p))

    @cached_method
    def _LKB_matrix_(self, braid, variab):
        """
        Compute the Lawrence-Krammer-Bigelow representation matrix.

        The variables of the matrix must be given. This actual
        computation is done in this helper method for caching
        purposes.

        INPUT:

        - ``braid`` -- tuple of integers. The Tietze list of the
          braid.

        - ``variab`` -- string. the names of the variables that will
          appear in the matrix. They must be given as a string,
          separated by a comma

        OUTPUT:

        The LKB matrix of the braid, with respect to the variables.

        TESTS::

            sage: B=BraidGroup(3)
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
        if len(braid)>1:
            A = self._LKB_matrix_(braid[:1], variab)
            for i in braid[1:]:
                A = A*self._LKB_matrix_((i,), variab)
            return A
        l = list(Set(range(n)).subsets(2))
        R = LaurentPolynomialRing(IntegerRing(), variab)
        q = R.gens()[0]
        t = R.gens()[1]
        if len(braid)==0:
            return identity_matrix(R, len(l), sparse=True)
        A = matrix(R, len(l), sparse=True)
        if braid[0]>0:
            i = braid[0]-1
            for m in range(len(l)):
                j = min(l[m])
                k = max(l[m])
                if i==j-1:
                    A[l.index(Set([i, k])), m] = q
                    A[l.index(Set([i, j])), m] = q*q-q
                    A[l.index(Set([j, k])), m] = 1-q
                elif i==j and not j==k-1:
                    A[l.index(Set([j, k])), m] = 0
                    A[l.index(Set([j+1, k])), m] = 1
                elif k-1==i and not k-1==j:
                    A[l.index(Set([j, i])), m] = q
                    A[l.index(Set([j, k])), m] = 1-q
                    A[l.index(Set([i, k])), m] = (1-q)*q*t
                elif i==k:
                    A[l.index(Set([j, k])), m] = 0
                    A[l.index(Set([j, k+1])), m] = 1
                elif i==j and j==k-1:
                    A[l.index(Set([j, k])), m] = -t*q*q
                else:
                    A[l.index(Set([j, k])), m] = 1
            return A
        else:
            i = -braid[0]-1
            for m in range(len(l)):
                j = min(l[m])
                k = max(l[m])
                if i==j-1:
                    A[l.index(Set([j-1, k])), m] = 1
                elif i==j and not j==k-1:
                    A[l.index(Set([j+1, k])), m] = q**(-1)
                    A[l.index(Set([j, k])), m] = 1-q**(-1)
                    A[l.index(Set([j, j+1])), m] = t**(-1)*q**(-1)-t**(-1)*q**(-2)
                elif k-1==i and not k-1==j:
                    A[l.index(Set([j, k-1])), m] = 1
                elif i==k:
                    A[l.index(Set([j, k+1])), m] = q**(-1)
                    A[l.index(Set([j, k])), m] = 1-q**(-1)
                    A[l.index(Set([k, k+1])), m] = -q**(-1)+q**(-2)
                elif i==j and j==k-1:
                    A[l.index(Set([j, k])), m] = -t**(-1)*q**(-2)
                else:
                    A[l.index(Set([j, k])), m] = 1
            return A

    def dimension_of_TL_space(self, drain_size):
        """
        Return the dimension of a particular Templerley--Lieb representation
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

            sage: import sage.combinat.diagram_algebras as da
            sage: B = BraidGroup(6)
            sage: dimensions = [B.dimension_of_TL_space(d)**2 for d in [0, 2, 4, 6]]
            sage: total_dim = sum(dimensions)
            sage: total_dim == len(list(da.temperley_lieb_diagrams(6)))
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
        trivalent trees marked by Jones--Wenzl projectors (see e.g. [Wan10]_).
        In practice, this is a list of non-negative integers whose first
        element is ``drain_size``, whose last element is `0`, and satisfying
        that consecutive integers have difference `1`. Moreover, the length
        of each basis element is `n + 1`.

        Given these rules, the list of lists is constructed recursively
        in the natural way.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands
          (both inclusive)

        OUTPUT:

        A list of basis elements, each of which is a list of integers.

        EXAMPLES:

        We calculate the basis for the appropriate vector space for `B_5` when
        `d = 3`::

            sage: B = BraidGroup(5)
            sage: B.TL_basis_with_drain(3)
            [[3, 4, 3, 2, 1, 0],
             [3, 2, 3, 2, 1, 0],
             [3, 2, 1, 2, 1, 0],
             [3, 2, 1, 0, 1, 0]]

        The number of basis elements hopefully correponds to the general
        formula for the dimension of the representation spaces::

            sage: B = BraidGroup(10)
            sage: d = 2
            sage: B.dimension_of_TL_space(d) == len(B.TL_basis_with_drain(d))
            True

        REFERENCES:

        .. [Wan10] Zhenghan Wang. Tolological quantum computation. Providence,
           RI: American Mathematical Society (AMS), 2010.
           ISBN 978-0-8218-4930-9
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
            else:
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
        auxmat = matrix(n-1, len(basis))
        for i in range(1, n):  # For each of the e_i
            for v in range(len(basis)):  # For each basis element
                tree = basis[v]
                if tree[i-1] < tree[i] and tree[i+1] < tree[i]:
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

        OUTPUT:

        A list of matrices corresponding to the representations of each
        of the standard generators and their inverses.

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
            R = LaurentPolynomialRing(IntegerRing(), 'A')
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
        Return the action of self in the free group F as mapping class group.

        This action corresponds to the action of the braid over the
        punctured disk, whose fundamental group is the free group on
        as many generators as strands.

        In Sage, this action is the result of multiplying a free group
        element with a braid. So you generally do not have to
        construct this action yourself.

        OUTPUT:

        A :class:`MappingClassGroupAction`.

        EXAMPLES ::

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
            Will try _r_action and _l_action
            Unknown result parent.
        """
        import operator
        if is_FreeGroup(S) and op == operator.mul and not self_on_left:
            return self.mapping_class_action(S)
        return None

    def Delta(self):
        r"""
        Return the `\Delta` element of the braid group.

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: B.Delta()
            s0*s1*s0*s2*s1*s0*s3*s2*s1*s0
        """
        n = self.strands()
        delta = Permutation([n-i for i in range(n)])
        return self._permutation_braid(delta)

    def _element_from_libbraiding(self, nf):
        """
        Return the element of ``self`` corresponding to the output
        of libbraiding.

        INPUT:

        - ``nf`` -- a list of lists, as returned by libbraiding

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: B._element_from_libbraiding([[-2], [2, 1], [1, 2], [2, 1]])
            (s0^-1*s1^-1*s2^-1*s3^-1*s0^-1*s1^-1*s2^-1*s0^-1*s1^-1*s0^-1)^2*s1*s0^2*s1^2*s0
            sage: B._element_from_libbraiding([[0]])
            1
        """
        if len(nf) == 1:
            return self.Delta() ** nf[0][0]
        from sage.misc.misc_c import prod
        return self.Delta() ** nf[0][0] * prod(self(i) for i in nf[1:])


def BraidGroup(n=None, names='s'):
    """
    Construct a Braid Group

    INPUT:

    - ``n`` -- integer or ``None`` (default). The number of
      strands. If not specified the ``names`` are counted and the
      group is assumed to have one more strand than generators.

    - ``names`` -- string or list/tuple/iterable of strings (default:
      ``'x'``). The generator names or name prefix.

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

        sage: B=BraidGroup(3); B.generators()
        (s0, s1)
        sage: BraidGroup(3, 'g').generators()
        (g0, g1)

    Since the word problem for the braid groups is solvable, their Cayley graph
    can be localy obtained as follows (see :trac:`16059`)::

        sage: def ball(group, radius):
        ....:     ret = set()
        ....:     ret.add(group.one())
        ....:     for length in range(1, radius):
        ....:         for w in Words(alphabet=group.gens(), length=length):
        ....:              ret.add(prod(w))
        ....:     return ret
        sage: B = BraidGroup(4)
        sage: GB = B.cayley_graph(elements=ball(B, 4), generators=B.gens()); GB
        Digraph on 31 vertices

    Since the braid group has nontrivial relations, this graph contains less
    vertices than the one associated to the free group (which is a tree)::

        sage: F = FreeGroup(3)
        sage: GF = F.cayley_graph(elements=ball(F, 4), generators=F.gens()); GF
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
        if isinstance(names, six.string_types):
            n = len(names.split(','))
        else:
            names = list(names)
            n = len(names)
    from sage.structure.category_object import normalize_names
    names = normalize_names(n, names)
    return BraidGroup_class(names)



class MappingClassGroupAction(Action):
    r"""
    The action of the braid group the free group as the mapping class
    group of the punctured disk.

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

    where $\sigma_i$ are the generators of the braid group on $n$
    strands, and $x_j$ the generators of the free group of rank $n$.

    You should left multiplication of the free group element by the
    braid to compute the action. Alternatively, use the
    :meth:`~sage.groups.braid.BraidGroup_class.mapping_class_action`
    method of the braid group to constuct this action.

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
        Right action by Braid group on 4 strands on Free Group on generators {x0, x1, x2, x3}
        sage: A(x0, s1)
        x0
        sage: A(x1, s1)
        x1*x2*x1^-1
        sage: A(x1^-1, s1)
        x1*x2^-1*x1^-1
    """
    def __init__(self, G, M, is_left=0):
        """
        TESTS::

            sage: B = BraidGroup(3)
            sage: G = FreeGroup('a, b, c')
            sage: B.mapping_class_action(G) # indirect doctest
            Right action by Braid group on 3 strands on Free Group on generators {a, b, c}
        """
        import operator
        Action.__init__(self, G, M, is_left, operator.mul)

    def _call_(self, x, b):
        """
        Return the action of ``b`` on ``x``.

        INPUT:

        - ``x`` -- a free group element.

        - ``b`` -- a braid.

        OUTPUT:

        A new braid.

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
            s=[]
            for i in t:
                if j==i and i>0:
                    s += [i, i+1, -i]
                elif j==-i and i<0:
                    s += [-i, i-1, i]
                elif j==-i and i>0:
                    s += [i+1]
                elif j==i and i<0:
                    s += [i-1]
                elif i>0 and j==i-1:
                    s += [i-1]
                elif i<0 and j==-i-1:
                    s += [i+1]
                elif i>0 and -j==i-1:
                    s += [-i, i-1, i]
                elif i<0 and j==i+1:
                    s += [i, i+1, -i]
                else:
                    s += [i]
            t = s
        return self.codomain()(t)
