r"""
Farey Symbol for arithmetic subgroups of `{\rm PSL}_2(\ZZ)`

AUTHORS:

- Hartmut Monien (08 - 2011)

based on the *KFarey* package by Chris Kurth. Implemented as C++ module
for speed.
"""
#*****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'
include 'sage/ext/cdefs.pxi'

include "farey.pxd"

import sage.rings.arith
from sage.rings.all import CC, RR
from sage.rings.integer cimport Integer
from sage.rings.infinity import infinity
from congroup_gammaH import is_GammaH
from congroup_gamma1 import is_Gamma1
from congroup_gamma0 import is_Gamma0
from congroup_gamma import is_Gamma
from congroup_sl2z import SL2Z
from sage.modular.cusps import Cusp

from sage.plot.all import Graphics
from sage.plot.colors import to_mpl_color
from sage.plot.misc import options, rename_keyword
from sage.plot.all import hyperbolic_arc, hyperbolic_triangle, text

from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

cdef class Farey:
    r"""
    A class for calculating Farey symbols of arithmetics subgroups of
    `{\rm PSL}_2(\ZZ)`.

    The arithmetic subgroup can be either any of
    the congruence subgroups implemented in Sage, i.e. Gamma, Gamma0,
    Gamma1 and GammaH or a subgroup of `{\rm PSL}_2(\ZZ)` which is
    given by a user written helper class defining membership in that
    group.

    REFERENCES:

    - Ravi S. Kulkarni, ''An arithmetic-geometric method in the study of the
      subgroups of the modular group'', `Amer. J. Math., 113(6):1053--1133,
      1991. <http://www.jstor.org/stable/2374900>`_

    INPUT:

    - `G` - an arithmetic subgroup of `{\rm PSL}_2(\ZZ)`

    EXAMPLES:

    Create a Farey symbol for the group `\Gamma_0(11)`::

        sage: f = FareySymbol(Gamma0(11)); f
        FareySymbol(Congruence Subgroup Gamma0(11))

    Calculate the generators::

         sage: f.generators()
         [
         [1 1]  [ 7 -2]  [ 8 -3]  [-1  0]
         [0 1], [11 -3], [11 -4], [ 0 -1]
         ]

    Pickling the FareySymbol and recovering it::

         sage: f == loads(dumps(f))
         True

    Calculate the index of `\Gamma_H(33, [2, 5])` in
    `{\rm PSL}_2(\ZZ)` via FareySymbol::

         sage: FareySymbol(GammaH(33, [2, 5])).index()
         48

    Calculate the generators of `\Gamma_1(4)`::

         sage: FareySymbol(Gamma1(4)).generators()
         [
         [1 1]  [-3  1]
         [0 1], [-4  1]
         ]

    Calculate the generators of the :meth:`example
    <sage.modular.arithgroup.arithgroup_perm.HsuExample10>` of an
    index 10 arithmetic subgroup given by Tim Hsu::

         sage: from sage.modular.arithgroup.arithgroup_perm import HsuExample10
         sage: FareySymbol(HsuExample10()).generators()
         [
         [1 2]  [-2  1]  [ 4 -3]
         [0 1], [-7  3], [ 3 -2]
         ]

    Calculate the generators of the group `\Gamma' =
    \Gamma_0(8)\cap\Gamma_1(4)` using a helper class to define group membership::

         sage: class GPrime:
         ....:   def __contains__(self, M):
         ....:       return M in Gamma0(8) and M in Gamma1(4)
         ....:

         sage: FareySymbol(GPrime()).generators()
         [
         [1 1]  [ 5 -1]  [ 5 -2]
         [0 1], [16 -3], [ 8 -3]
         ]

    Calculate cusps of arithmetic subgroup defined via permutation group::

        sage: L = SymmetricGroup(4)('(1, 2, 3)')

        sage: R = SymmetricGroup(4)('(1, 2, 4)')

        sage: FareySymbol(ArithmeticSubgroup_Permutation(L, R)).cusps()
        [-1, Infinity]

    Calculate the left coset representation of `\Gamma_H(8, [3])`::

         sage: FareySymbol(GammaH(8, [3])).coset_reps()
         [
         [1 0]  [ 4 -1]  [ 3 -1]  [ 2 -1]  [ 1 -1]  [ 3 -1]  [ 2 -1]  [-1  0]
         [0 1], [ 1  0], [ 1  0], [ 1  0], [ 1  0], [ 4 -1], [ 3 -1], [ 3 -1],
         [ 1 -1]  [-1  0]  [ 0 -1]  [-1  0]
         [ 2 -1], [ 2 -1], [ 1 -1], [ 1 -1]
         ]
    """
    cdef cpp_farey *this_ptr
    cdef object group

    def __cinit__(self, group, data=None):
        r"""
        Initialize FareySymbol::

            sage: FareySymbol(Gamma0(23))
            FareySymbol(Congruence Subgroup Gamma0(23))
        """
        self.group = group
        # if data is present we want to restore
        if data is not None:
            sig_on()
            self.this_ptr = new cpp_farey(data)
            sig_off()
            return
        ## to accelerate the calculation of the FareySymbol
        ## we implement the tests for the standard congruence groups
        ## in the c++ module. For a general group the test if an element
        ## of SL2Z is in the group the python __contains__ attribute
        ## of the group is called
        cdef int p
        if hasattr(group, "level"):
            p = group.level()
        if group == SL2Z:
            sig_on()
            self.this_ptr = new cpp_farey()
            sig_off()
        elif is_Gamma0(group):
            sig_on()
            self.this_ptr = new cpp_farey(group, new is_element_Gamma0(p))
            sig_off()
        elif is_Gamma1(group):
            sig_on()
            self.this_ptr = new cpp_farey(group, new is_element_Gamma1(p))
            sig_off()
        elif is_Gamma(group):
            sig_on()
            self.this_ptr = new cpp_farey(group, new is_element_Gamma(p))
            sig_off()
        elif is_GammaH(group):
            sig_on()
            l = group._GammaH_class__H
            self.this_ptr = new cpp_farey(group, new is_element_GammaH(p, l))
            sig_off()
        else:
            sig_on()
            self.this_ptr = new cpp_farey(group)
            sig_off()

    def __deallocpp__(self):
        r"""
        Remove reference to FareySymbol::

            sage: F = FareySymbol(Gamma0(23))

            sage: del F

        """
        del self.this_ptr

    @cached_method
    def pairing_matrices_to_tietze_index(self):
        r"""
        Obtain the translation table from pairing matrices
        to generators.

        The result is cached.

        OUTPUT:

        a list where the `i`-th entry is a nonzero integer `k`,
        such that if `k > 0` then the `i`-th pairing matrix is (up to sign)
        the `(k-1)`-th generator and, if `k < 0`, then the `i`-th pairing
        matrix is (up to sign) the inverse of the `(-k-1)`-th generator.

        EXAMPLES::

            sage: F = Gamma0(40).farey_symbol()
            sage: table = F.pairing_matrices_to_tietze_index()
            sage: table[12]
            -2
            sage: F.pairing_matrices()[12]
            [  3  -1]
            [ 40 -13]
            sage: F.generators()[1]**-1
            [ -3   1]
            [-40  13]
        """
        gens_dict = {g:i+1 for i,g in enumerate(self.generators())}
        ans = []
        for pm in self.pairing_matrices():
            a, b, c, d = pm.matrix().list()
            newval = gens_dict.get(SL2Z([a, b, c, d]))
            if newval is not None:
                ans.append(newval)
                continue
            newval = gens_dict.get(SL2Z([-a, -b, -c, -d]))
            if newval is not None:
                ans.append(newval)
                continue
            newval = gens_dict.get(SL2Z([d, -b, -c, a]))
            if newval is not None:
                ans.append(-newval)
                continue
            newval = gens_dict.get(SL2Z([-d, b, c, -a]))
            if newval is not None:
                ans.append(-newval)
                continue
            raise RuntimeError("This should have not happened")
        return ans

    def word_problem(self, M, output = 'standard'):
        r"""
        Solve the word problem (up to sign) using this Farey symbol.
        See also :meth:`.syllables` and :meth:`.tietze`.

        INPUT:

        - ``M`` -- An element `M` of `{\rm SL}_2(\ZZ)`.
        - ``output`` -- (default: ``'standard'``) Should be one of ``'standard'``,
          ``'syllables'``, ``'gens'``.

        OUTPUT:

        A solution to the word problem for the matrix `M`.
        The format depends on the ``output`` parameter, as follows.

        - ``standard`` returns the so called the Tietze representation,
          consists of a tuple of nonzero integers `i`, where if `i` > 0
          then it indicates the `i`th generator (that is, ``self.generators()[0]``
          would correspond to `i` = 1), and if `i` < 0 then it indicates
          the inverse of the `i`-th generator.
        - ``syllables`` returns a tuple of tuples of the form `(i,n)`, where
          `(i,n)` represents ``self.generators()[i] ^ n``,
          whose product equals `M` up to sign.
        - ``gens`` returns tuple of tuples of the form `(g,n)`,
          `(g,n)` such that the product of the matrices `g^n`
          equals `M` up to sign.

        EXAMPLES::

            sage: F = Gamma0(30).farey_symbol()
            sage: gens = F.generators()
            sage: g = gens[3] * gens[10] * gens[8]^-1 * gens[5]
            sage: g
            [-628597   73008]
            [-692130   80387]
            sage: F.word_problem(g)
            (4, 11, -9, 6)
            sage: g = gens[3] * gens[10]^2 * gens[8]^-1 * gens[5]
            sage: g
            [-5048053   586303]
            [-5558280   645563]
            sage: F.word_problem(g, output = 'gens')
            ((
            [109 -10]
            [120 -11], 1
            ),
             (
            [ 19  -7]
            [ 30 -11], 2
            ),
             (
            [ 49  -9]
            [ 60 -11], -1
            ),
             (
            [17 -2]
            [60 -7], 1
            ))
            sage: F.word_problem(g, output = 'syllables')
            ((3, 1), (10, 2), (8, -1), (5, 1))
        """
        if output not in ['standard', 'syllables', 'gens']:
            raise ValueError('Unrecognized output format')
        cdef Integer a = M.d()
        cdef Integer b = -M.b()
        cdef Integer c = -M.c()
        cdef Integer d = M.a()
        sig_on()
        result = self.this_ptr.word_problem(a.value, b.value, c.value, d.value)
        sig_off()
        result.reverse()
        V = self.pairing_matrices_to_tietze_index()
        tietze = tuple(V[o-1] if o > 0 else -V[-o-1] for o in result)
        if output == 'standard':
            return tietze
        from itertools import groupby
        if output == 'syllables':
            return tuple((a-1,len(list(g))) if a > 0 else (-a-1,-len(list(g))) for a,g in groupby(tietze))
        else: # output == 'gens'
            gens = self.generators()
            return tuple((gens[a-1],len(list(g))) if a > 0 else (gens[-a-1],-len(list(g))) for a,g in groupby(tietze))

    def __contains__(self, M):
        r"""
        Tests if element is in the arithmetic group of the Farey symbol
        via LLT algorithm.

        EXAMPLES::

            sage: SL2Z([0, -1, 1, 0]) in FareySymbol(Gamma0(6))
            False

            sage: SL2Z([1, 1, 0, 1]) in FareySymbol(Gamma0(6))
            True
        """
        cdef Integer a = M.a()
        cdef Integer b = M.b()
        cdef Integer c = M.c()
        cdef Integer d = M.d()
        sig_on()
        result = self.this_ptr.is_element(a.value, b.value, c.value, d.value)
        sig_off()
        return result

    def __cmp__(self, other):
        r"""
        Compare self to others.

        EXAMPLES::

            sage: FareySymbol(Gamma(2)) == FareySymbol(Gamma0(7))
            False

            sage: FareySymbol(Gamma0(23)) == loads(dumps(FareySymbol(Gamma0(23))))
            True
        """
        cmp_fcts = [lambda fs: fs.coset_reps(),
                    lambda fs: fs.cusps(),
                    lambda fs: fs.fractions()]

        for cf in cmp_fcts:
            c = cmp(cf(self), cf(other))
            if c != 0:
                return c

        return c

    def __reduce__(self):
        r"""
        Serialization for pickling::

            sage: FareySymbol(Gamma0(4)).__reduce__()
            (<type 'sage.modular.arithgroup.farey_symbol.Farey'>, ...))

        """
        return Farey, (self.group, self.this_ptr.dumps())

    def __repr__(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: FareySymbol(Gamma0(23)).__repr__()
            'FareySymbol(Congruence Subgroup Gamma0(23))'
        """
        if hasattr(self.group, "_repr_"):
            return "FareySymbol(%s)" % self.group._repr_()
        elif hasattr(self.group, "__repr__"):
            return "FareySymbol(%s)" % self.group.__repr__()
        else:
            return "FareySymbol(?)"

    def _latex_(self, forced_format=None):
        r"""
        Return the LaTeX representation of ``self``.

        INPUT:

        - ``forced_format`` -- A format string ('plain' or 'xymatrix')
                               or ``None``.

        EXAMPLES::

            sage: FareySymbol(Gamma0(11))._latex_(forced_format = 'plain')
            '\\left( -\\infty\\underbrace{\\quad}_{1} 0\\underbrace{\\quad}_{2} \\frac{1}{3}\\underbrace{\\quad}_{3} \\frac{1}{2}\\underbrace{\\quad}_{2} \\frac{2}{3}\\underbrace{\\quad}_{3} 1\\underbrace{\\quad}_{1} \\infty\\right)'
            sage: FareySymbol(Gamma0(11))._latex_(forced_format = 'xymatrix')
            '\\begin{xy}\\xymatrix{& -\\infty \\ar@{-}@/_1pc/[r]_{1}& 0 \\ar@{-}@/_1pc/[r]_{2}& \\frac{1}{3} \\ar@{-}@/_1pc/[r]_{3}& \\frac{1}{2} \\ar@{-}@/_1pc/[r]_{2}& \\frac{2}{3} \\ar@{-}@/_1pc/[r]_{3}& 1 \\ar@{-}@/_1pc/[r]_{1}& \\infty }\\end{xy}'

            sage: if '\\xymatrix' in sage.misc.latex.latex.mathjax_avoid_list():
            ....:      'xymatrix' not in FareySymbol(Gamma0(11))._latex_()
            ....: else:
            ....:     'xymatrix' in FareySymbol(Gamma0(11))._latex_()
            True
        """
        if forced_format == 'plain' or \
           (forced_format is None and '\\xymatrix' in latex.mathjax_avoid_list()):
            # output not using xymatrix
            s = r'\left( -\infty'
            a = [x._latex_() for x in self.fractions()] + ['\infty']
            b = self.pairings()
            for i in xrange(len(a)):
                u = b[i]
                if u == -3:
                    u = r'\bullet'
                elif u == -2:
                    u = r'\circ'
                s += r'\underbrace{\quad}_{%s} %s' % (u, a[i])
            return s + r'\right)'
        else:
            # output using xymatrix
            s = r'\begin{xy}\xymatrix{& -\infty '
            f = [x._latex_() for x in self.fractions()]+[r'\infty']
            f.reverse()
            for p in self.pairings():
                if p >= 0:
                    s += r'\ar@{-}@/_1pc/[r]_{%s}' % p
                elif p == -2:
                    s += r'\ar@{-}@/_1pc/[r]_{\circ}'
                elif p == -3:
                    s += r'\ar@{-}@/_1pc/[r]_{\bullet}'
                s += r'& %s ' % f.pop()
            s += r'}\end{xy}'
            return s

    def index(self):
        r"""
        Return the index of the arithmetic group of the FareySymbol
        in `{\rm PSL}_2(\ZZ)`.

        EXAMPLES::

            sage: [FareySymbol(Gamma0(n)).index() for n in range(1, 16)]
            [1, 3, 4, 6, 6, 12, 8, 12, 12, 18, 12, 24, 14, 24, 24]
        """
        return self.this_ptr.index()

    def genus(self):
        r"""
        Return the genus of the arithmetic group of the FareySymbol.

        EXAMPLES::

            sage: [FareySymbol(Gamma0(n)).genus() for n in range(16, 32)]
            [0, 1, 0, 1, 1, 1, 2, 2, 1, 0, 2, 1, 2, 2, 3, 2]
        """
        return self.this_ptr.genus()

    def level(self):
        r"""
        Return the level of the arithmetic group of the FareySymbol.

        EXAMPLES::

            sage: [FareySymbol(Gamma0(n)).level() for n in range(1, 16)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        return self.this_ptr.level()

    def nu2(self):
        r"""
        Return the number of elliptic points of order two.

        EXAMPLES::

            sage: [FareySymbol(Gamma0(n)).nu2() for n in range(1, 16)]
            [1, 1, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0]
        """
        return self.this_ptr.nu2()

    def nu3(self):
        r"""
        Return the number of elliptic points of order three.

        EXAMPLES::

            sage: [FareySymbol(Gamma0(n)).nu3() for n in range(1, 16)]
            [1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0]
        """
        return self.this_ptr.nu3()

    def coset_reps(self):
        r"""
        Left coset of the arithmetic group of the FareySymbol.

        EXAMPLES:

        Calculate the left coset of `\Gamma_0(6)`::

            sage: FareySymbol(Gamma0(6)).coset_reps()
            [
            [1 0]  [ 3 -1]  [ 2 -1]  [ 1 -1]  [ 2 -1]  [ 3 -2]  [ 1 -1]  [-1  0]
            [0 1], [ 1  0], [ 1  0], [ 1  0], [ 3 -1], [ 2 -1], [ 2 -1], [ 2 -1],
            [ 1 -1]  [ 0 -1]  [-1  0]  [-2  1]
            [ 3 -2], [ 1 -1], [ 1 -1], [ 1 -1]
            ]
        """
        return self.this_ptr.get_coset()

    def generators(self):
        r"""
        Minmal set of generators of the group of the FareySymbol.

        EXAMPLES:

        Calculate the generators of `\Gamma_0(6)`::

            sage: FareySymbol(Gamma0(6)).generators()
            [
            [1 1]  [ 5 -1]  [ 7 -3]  [-1  0]
            [0 1], [ 6 -1], [12 -5], [ 0 -1]
            ]

        Calculate the generators of `{\rm SL}_2(\ZZ)`::

            sage: FareySymbol(SL2Z).generators()
            [
            [ 0 -1]  [ 0 -1]
            [ 1  0], [ 1 -1]
            ]

        The unique index 2 even subgroup and index 4 odd subgroup each get handled correctly::

            sage: FareySymbol(ArithmeticSubgroup_Permutation(S2="(1,2)", S3="()")).generators()
            [
            [ 0  1]  [-1  1]
            [-1 -1], [-1  0]
            ]
            sage: FareySymbol(ArithmeticSubgroup_Permutation(S2="(1,2, 3, 4)", S3="(1,3)(2,4)")).generators()
            [
            [ 0  1]  [-1  1]
            [-1 -1], [-1  0]
            ]
        """
        return self.this_ptr.get_generators()

    def fractions(self):
        r"""
        Fractions of the FareySymbol.

        EXAMPLES::

            sage: FareySymbol(Gamma(4)).fractions()
            [0, 1/2, 1, 3/2, 2, 5/2, 3, 7/2, 4]
        """
        return self.this_ptr.get_fractions()

    def pairings(self):
        r"""
        Pairings of the sides of the fundamental domain of the Farey symbol
        of the arithmetic group.

        The sides of the hyperbolic polygon are
        numbered 0, 1, ... from left to right. Conventions: even pairings are
        denoted by -2, odd pairings by -3 while free pairings are denoted by
        an integer number greater than zero.

        EXAMPLES:

        Odd pairings::

            sage: FareySymbol(Gamma0(7)).pairings()
            [1, -3, -3, 1]

        Even and odd pairings::

            FareySymbol(Gamma0(13)).pairings()
            [1, -3, -2, -2, -3, 1]

        Only free pairings::

            sage: FareySymbol(Gamma0(23)).pairings()
            [1, 2, 3, 5, 3, 4, 2, 4, 5, 1]
        """
        return self.this_ptr.get_pairings()

    def paired_sides(self):
        r"""
        Pairs of index of the sides of the fundamental domain of the
        Farey symbol of the arithmetic group. The sides of the
        hyperbolic polygon are numbered 0, 1, ... from left to right.

        .. image:: ../../../media/pairing.png

        EXAMPLES::

            sage: FareySymbol(Gamma0(11)).paired_sides()
            [(0, 5), (1, 3), (2, 4)]

        indicating that the side 0 is paired with 5, 1 with 3 and 2 with 4.
        """
        return self.this_ptr.get_paired_sides()

    def pairing_matrices(self):
        r"""
        Pairing matrices of the sides of the fundamental domain. The sides
        of the hyperbolic polygon are numbered 0, 1, ... from left to right.

        EXAMPLES::

            sage: FareySymbol(Gamma0(6)).pairing_matrices()
            [
            [1 1]  [ 5 -1]  [ 7 -3]  [ 5 -3]  [ 1 -1]  [-1  1]
            [0 1], [ 6 -1], [12 -5], [12 -7], [ 6 -5], [ 0 -1]
            ]
        """

        return self.this_ptr.get_pairing_matrices()

    def cusps(self):
        r"""
        Cusps of the FareySymbol.

        EXAMPLES::

            sage: FareySymbol(Gamma0(6)).cusps()
            [0, 1/3, 1/2, Infinity]
        """
        return self.this_ptr.get_cusps()+[Cusp(infinity)]

    def cusp_widths(self):
        r"""
        Cusps widths of the FareySymbol.

        EXAMPLES::

            sage: FareySymbol(Gamma0(6)).cusp_widths()
            [6, 2, 3, 1]
        """
        return self.this_ptr.get_cusp_widths()

    def cusp_class(self, c):
        r"""
        Cusp class of a cusp in the FareySymbol.

        INPUT:

        ``c`` -- a cusp

        EXAMPLES::

            sage: FareySymbol(Gamma0(12)).cusp_class(Cusp(1, 12))
            5
        """
        cusp = Cusp(c)
        cdef Integer p = cusp.numerator()
        cdef Integer q = cusp.denominator()
        sig_on()
        result = self.this_ptr.get_cusp_class(p.value, q.value)
        sig_off()
        return result

    def reduce_to_cusp(self, r):
        r"""
        Transformation of a rational number to cusp representative.

        INPUT:

        ``r`` -- a rational number

        EXAMPLES::

            sage: FareySymbol(Gamma0(12)).reduce_to_cusp(5/8)
            [ 5  -3]
            [12  -7]

        Reduce 11/17 to a cusp of for HsuExample10()::

            sage: from sage.modular.arithgroup.arithgroup_perm import HsuExample10
            sage: f = FareySymbol(HsuExample10())
            sage: f.reduce_to_cusp(11/17)
            [14 -9]
            [-3  2]
            sage: _.acton(11/17)
            1
            sage: f.cusps()[f.cusp_class(11/17)]
            1
        """
        cdef Integer p = r.numerator()
        cdef Integer q = r.denominator()
        sig_on()
        result = self.this_ptr.get_transformation_to_cusp(p.value, q.value)
        sig_off()
        return result

    @rename_keyword(rgbcolor='color')
    @options(alpha=1, fill=True, thickness=1, color='lightgray',
             color_even='white',
             zorder=2, linestyle='solid', show_pairing=True,
             tesselation='Dedekind', ymax=1)
    def fundamental_domain(self, **options):
        r"""
        Plot a fundamental domain of an arithmetic subgroup of
        `{\rm PSL}_2(\ZZ)` corresponding to the Farey symbol.

        OPTIONS:

        - ``fill`` -- boolean (default ``True``) fill the fundamental domain

        - ``linestyle`` -- string (default: 'solid') The style of the line,
          which is one of 'dashed', 'dotted', 'solid', 'dashdot', or '--',
          ':', '-', '-.', respectively

        - ``color`` -- (default: 'lightgray') fill color; fill
          color for odd part of Dedekind tesselation.

        - ``show_pairing`` -- boolean (default: ``True``) flag for pairing

        - ``tesselation`` -- (default: 'Dedekind') The type of
          hyperbolic tesselation which is one of
          'coset', 'Dedekind' or ``None`` respectively

        - ``color_even`` -- fill color for even parts of Dedekind
          tesselation (default 'white'); ignored for other tesselations

        - ``thickness`` -- float (default: `1`) the thickness of the line

        - ``ymax`` -- float (default: `1`) maximal height

        EXAMPLES:

        For example, to plot the fundamental domain of `\Gamma_0(11)`
        with pairings use the following command::

            sage: FareySymbol(Gamma0(11)).fundamental_domain()
            Graphics object consisting of 54 graphics primitives

        indicating that side 1 is paired with side 3 and side 2 is
        paired with side 4, see also :meth:`.paired_sides`.

        To plot the fundamental domain of `\Gamma(3)` without pairings
        use the following command::

            sage: FareySymbol(Gamma(3)).fundamental_domain(show_pairing=False)
            Graphics object consisting of 48 graphics primitives

        Plot the fundamental domain of `\Gamma_0(23)` showing the left
        coset representatives::

            sage: FareySymbol(Gamma0(23)).fundamental_domain(tesselation='coset')
            Graphics object consisting of 58 graphics primitives

        The same as above but with a custom linestyle::

            sage: FareySymbol(Gamma0(23)).fundamental_domain(tesselation='coset', linestyle=':', thickness='2')
            Graphics object consisting of 58 graphics primitives

        """
        from sage.plot.colors import rainbow
        I = CC(0, 1)
        w = RR(3).sqrt()
        L = 1000
        g = Graphics()
        ## show coset
        for x in self.coset_reps():
            a, b, c, d = x[1, 1], -x[0, 1], -x[1, 0], x[0, 0]
            A, B = CC(0, L), CC(0, L)
            if d != 0:
                A = b / d
            if c != 0:
                B = a / c
            C = (a*c+b*d+(a*d+b*c)/2+I*w/2)/(c*c+c*d+d*d)
            D = (a*c+b*d + CC(0, 1))/(c*c+d*d)
            if options['tesselation'] == 'Dedekind':
                g += hyperbolic_triangle(A, D, C,
                                         alpha=options['alpha'],
                                         color=options['color'],
                                         fill=options['fill'],
                                         linestyle=options['linestyle'],
                                         thickness=options['thickness'])
                g += hyperbolic_triangle(D, C, B,
                                         alpha=options['alpha'],
                                         color=options['color_even'],
                                         fill=options['fill'],
                                         linestyle=options['linestyle'],
                                         thickness=options['thickness'])
                g += hyperbolic_triangle(A, D, C, color='gray')
                g += hyperbolic_triangle(D, C, B, color='gray')
            elif options['tesselation'] == 'coset':
                g += hyperbolic_triangle(A, B, C,
                                         alpha=options['alpha'],
                                         color=options['color'],
                                         fill=options['fill'],
                                         linestyle=options['linestyle'],
                                         thickness=options['thickness'])
                g += hyperbolic_triangle(A, B, C, color='gray',
                                         linestyle=options['linestyle'],
                                         thickness=options['thickness'])
            else:
                g += hyperbolic_triangle(A, B, C,
                                         alpha=options['alpha'],
                                         color=options['color'],
                                         fill=options['fill'],
                                         linestyle=options['linestyle'],
                                         thickness=options['thickness'])
        ## show pairings
        p = self.pairings()
        x = self.fractions()
        if options['show_pairing']:
            rc = rainbow(max(p)-min(p)+1)
            if p[0] > 0:
                g += hyperbolic_arc(CC(0, L), x[0], color=rc[p[0]-min(p)],
                                    linestyle=options['linestyle'],
                                    thickness=options['thickness'])
            if p[-1] > 0:
                g += hyperbolic_arc(CC(0, L), x[-1], color=rc[p[-1]-min(p)],
                                    linestyle=options['linestyle'],
                                    thickness=options['thickness'])
            for i in range(len(x)-1):
                if p[i+1] > 0:
                    g += hyperbolic_arc(x[i], x[i+1], color=rc[p[i+1]-min(p)],
                                        linestyle=options['linestyle'],
                                        thickness=options['thickness'])
        d = g.get_minmax_data()
        g.set_axes_range(d['xmin'], d['xmax'], 0, options['ymax'])
        return g


#--- conversions ------------------------------------------------------------

cdef public long convert_to_long(n):
    cdef long m = n
    return m

cdef public object convert_to_Integer(mpz_class a):
    A = Integer()
    A.set_from_mpz(a.get_mpz_t())
    return A

cdef public object convert_to_rational(mpq_class r):
    a = Integer()
    a.set_from_mpz(r.get_num_mpz_t())
    b = Integer()
    b.set_from_mpz(r.get_den_mpz_t())
    return a/b

cdef public object convert_to_cusp(mpq_class r):
    a = Integer()
    a.set_from_mpz(r.get_num_mpz_t())
    b = Integer()
    b.set_from_mpz(r.get_den_mpz_t())
    return Cusp(a/b)

cdef public object convert_to_SL2Z(cpp_SL2Z M):
    a = convert_to_Integer(M.a())
    b = convert_to_Integer(M.b())
    c = convert_to_Integer(M.c())
    d = convert_to_Integer(M.d())
    return SL2Z([a, b, c, d])
