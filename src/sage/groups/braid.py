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

AUTHOR:

- Miguel Angel Marco Buzunariz
- Volker Braun
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


from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.groups.free_group import FreeGroup, is_FreeGroup
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.matrix.constructor import identity_matrix, matrix
from sage.combinat.permutation import Permutation
from sage.categories.action import Action
from sage.sets.set import Set
from sage.groups.finitely_presented import FinitelyPresentedGroup, FinitelyPresentedGroupElement


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

    def __cmp__(self, other):
        """
        Compare ``self`` and ``other``

        TESTS::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1])
            sage: c = B([2, 1, 2])
            sage: b == c #indirect doctest
            True
            sage: b.__cmp__(c^(-1))
            -1
            sage: B([]).__cmp__(B.one())
            0
        """
        if self.Tietze()==other.Tietze():
            return 0
        nfself = map(lambda i: i.Tietze(), self.left_normal_form())
        nfother = map(lambda i: i.Tietze(), other.left_normal_form())
        return cmp(nfself, nfother)

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
            if i>0:
                latexrepr = latexrepr+"\sigma_{%s}"%i
            if i<0:
                latexrepr = latexrepr+"\sigma_{%s}^{-1}"%(-i)
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

    def burau_matrix(self, var='t'):
        """
        Return the Burau matrix of the braid.

        INPUT:

        - ``var`` -- string (default: ``'t'``). The name of the
          variable in the entries of the matrix.

        OUTPUT:

        The Burau matrix of the braid. It is a matrix whose entries
        are Laurent polynomials in the variable ``var``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.inject_variables()
            Defining s0, s1, s2
            sage: b=s0*s1/s2/s1
            sage: b.burau_matrix()
            [     -t + 1           0    -t^2 + t         t^2]
            [          1           0           0           0]
            [          0           0           1           0]
            [          0        t^-2 t^-1 - t^-2    1 - t^-1]
            sage: s2.burau_matrix('x')
            [     1      0      0      0]
            [     0      1      0      0]
            [     0      0 -x + 1      x]
            [     0      0      1      0]

        REFERENCES:

            http://en.wikipedia.org/wiki/Burau_representation
        """
        R = LaurentPolynomialRing(IntegerRing(), var)
        t = R.gen()
        M = identity_matrix(R, self.strands())
        for i in self.Tietze():
            A = identity_matrix(R, self.strands())
            if i>0:
                A[i-1, i-1] = 1-t
                A[i, i] = 0
                A[i, i-1] = 1
                A[i-1, i] = t
            if i<0:
                A[-1-i, -1-i] = 0
                A[-i, -i] = 1-t**(-1)
                A[-1-i, -i] = 1
                A[-i, -1-i] = t**(-1)
            M=M*A
        return M

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
            sage: b.plot(color=["red", "blue", "red", "blue"])

            sage: B.<s,t> = BraidGroup(3)
            sage: b = t^-1*s^2
            sage: b.plot(orientation="left-right", color="red")
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
            sage: b.plot3d(color="red")
            sage: b.plot3d(color=["red", "blue", "red", "blue"])
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
        Return the Lawence-Krammer-Bigelow representation matrix.

        The matrix is expressed in the basis $\{e_{i, j} \mid 1\\leq i
        < j \leq n\}$, where the indices are ordered
        lexicographically.  It is a matrix whose entries are in the
        ring of Laurent polynomials on the given variables.  By
        default, the variables are ``'x'`` and ``'y'``.

        INPUT:

        - ``variables`` -- string (default: ``'x,y'``). A string
          containing the names of the variables, separated by a comma.

        OUTPUT:

        The matrix corresponding to the Lawence-Krammer-Bigelow representation of the braid.

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

        .. [Bigelow] Bigelow, Stephen J. The Lawrence-Krammer representation. arXiv:math/0204057v1
        """
        return self.parent()._LKB_matrix_(self.Tietze(), variab=variables)

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

        .. [Dynnikov07] I. Dynnikov and B. Wiest, On the complexity of braids,
           J. Europ. Math. Soc. 9 (2007)
        .. [Dehornoy] P. Dehornoy, Le probleme d'isotopie des tresses, in
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
        return map(T, coord)

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
        return tuple( [P._permutation_braid(delta).__pow__(a)] +
                      map(lambda i:P._permutation_braid(i), l) )

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
                form = map(lambda a: Delta*a*Delta, form)
                form.append(Delta*Permutation((-i, -i+1)))
        i = j = 0
        while j<len(form):
            while i<len(form)-j-1:
                e = form[i].inverse().descents()
                s = form[i+1].descents()
                S = set(s).difference(set(e))
                while S!=set([]):
                    a = list(S)[0]
                    form[i] = form[i]*Permutation((a+1, a+2))
                    form[i+1] = Permutation((a+1, a+2))*form[i+1]
                    e = form[i].inverse().descents()
                    s = form[i+1].descents()
                    S = set(s).difference(set(e))
                if form[i+1].length()==0:
                    form.pop(i+1)
                    i = 0
                else:
                    i += 1
            j += 1
            i = 0
        form = filter(lambda a: a.length()>0, form)
        while form!=[] and form[0]==Delta:
            form.pop(0)
            delta = delta-1
        return tuple([-delta]+form)


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
        return Braid(self, x)

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
        Let the coercion system discover actions of the braid group on free groups.

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
        if is_FreeGroup(S) and op==operator.mul and not self_on_left:
            return self.mapping_class_action(S)
        return None



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
        if isinstance(names, basestring):
            n = len(names.split(','))
        else:
            names = list(names)
            n = len(names)
    from sage.structure.parent import normalize_names
    names = tuple(normalize_names(n, names))
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
