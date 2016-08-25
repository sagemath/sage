# -*- coding: utf-8 -*-
"""
Discrete Subgroups of `\\ZZ^n`.

AUTHORS:

- Martin Albrecht (2014-03): initial version

- Jan Pöschko (2012-08): some code in this module was taken from Jan Pöschko's
  2012 GSoC project

TESTS::

    sage: from sage.modules.free_module_integer import IntegerLattice
    sage: L = IntegerLattice(random_matrix(ZZ, 10, 10))
    sage: TestSuite(L).run()

"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2012 Jan Poeschko <jan@poeschko.com>
#       Copyright (C) 2014 Martin Albrecht <martinralbecht@googlemail.com>
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

from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.method_decorator import MethodDecorator
from sage.modules.free_module import FreeModule_submodule_with_basis_pid, FreeModule_ambient_pid
from sage.modules.free_module_element import vector
from sage.modules.vector_integer_dense import Vector_integer_dense
from sage.structure.parent import Parent
from sage.rings.number_field.number_field_element import OrderElement_absolute

def IntegerLattice(basis, lll_reduce=True):
    r"""
    Construct a new integer lattice from ``basis``.

    INPUT:

    - ``basis`` -- can be one of the following:

      - a list of vectors

      - a matrix over the integers

      - an element of an absolute order

    - ``lll_reduce`` -- (default: ``True``) run LLL reduction on the basis
      on construction.

    EXAMPLES:

    We construct a lattice from a list of rows::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: IntegerLattice([[1,0,3], [0,2,1], [0,2,7]])
        Free module of degree 3 and rank 3 over Integer Ring
        User basis matrix:
        [-2  0  0]
        [ 0  2  1]
        [ 1 -2  2]

    Sage includes a generator for hard lattices from cryptography::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: A = sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)
        sage: IntegerLattice(A)
        Free module of degree 10 and rank 10 over Integer Ring
        User basis matrix:
        [ 0  1  2  0  1  2 -1  0 -1 -1]
        [ 0  1  0 -3  0  0  0  0  3 -1]
        [ 1  1 -1  0 -3  0  0  1  2 -2]
        [-1  2 -1 -1  2 -2  1 -1  0 -1]
        [ 1  0 -4  2  0  1 -2 -1  0  0]
        [ 2  3  0  1  1  0 -2  3  0  0]
        [-2 -3 -2  0  0  1 -1  1  3 -2]
        [-3  0 -1  0 -2 -1 -2  1 -1  1]
        [ 1  4 -1  1  2  2  1  0  3  1]
        [-1 -1  0 -3 -1  2  2  3 -1  0]

    You can also construct the lattice directly::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True, lattice=True)
        Free module of degree 10 and rank 10 over Integer Ring
        User basis matrix:
        [ 0  1  2  0  1  2 -1  0 -1 -1]
        [ 0  1  0 -3  0  0  0  0  3 -1]
        [ 1  1 -1  0 -3  0  0  1  2 -2]
        [-1  2 -1 -1  2 -2  1 -1  0 -1]
        [ 1  0 -4  2  0  1 -2 -1  0  0]
        [ 2  3  0  1  1  0 -2  3  0  0]
        [-2 -3 -2  0  0  1 -1  1  3 -2]
        [-3  0 -1  0 -2 -1 -2  1 -1  1]
        [ 1  4 -1  1  2  2  1  0  3  1]
        [-1 -1  0 -3 -1  2  2  3 -1  0]

    We construct an ideal lattice from an element of an absolute order::

        sage: K.<a>  = CyclotomicField(17)
        sage: O = K.ring_of_integers()
        sage: f = O.random_element(); f
        -a^15 - a^12 - a^10 - 8*a^9 - a^8 - 4*a^7 + 3*a^6 + a^5 + 2*a^4 + 8*a^3 - a^2 + a + 1

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: IntegerLattice(f)
        Free module of degree 16 and rank 16 over Integer Ring
        User basis matrix:
        [ 1  1 -1  8  2  1  3 -4 -1 -8 -1  0 -1  0  0 -1]
        [-1  0  1  1 -1  8  2  1  3 -4 -1 -8 -1  0 -1  0]
        [ 1  1  0  1  2  2  0  9  3  2  4 -3  0 -7  0  1]
        [ 1  0  1  1  0  1  2  2  0  9  3  2  4 -3  0 -7]
        [ 2 -5 -2 -9 -2 -1 -2 -1 -1 -2 -1  0  0 -2  7  1]
        [ 1  4  0 -5 -3  0  3  5 -2  2  0 -7  4  0 -6 -2]
        [-7  4  0 -6 -2  0  1  4  0 -5 -3  0  3  5 -2  2]
        [-1  0  0 -1  0  1  1 -1  8  2  1  3 -4 -1 -8 -1]
        [-1 -1 -2  4  1  9  1 -1  0  1 -8 -1 -1 -4  3  2]
        [-1 -2  6 -6  8  1 -3  5  3  1  1  0 -2  4  3  2]
        [ 4  8  2  7 -3  2 -1  2  0 -4  0  3  6  3  0 -4]
        [ 0 -1  0  1  1 -1  8  2  1  3 -4 -1 -8 -1  0 -1]
        [ 2 -7 -1  0 -2  5  2  9  2  1  2  1  1  2  1  0]
        [-1  7 -5  9  2 -2  6  4  2  2  1 -1  5  4  3  1]
        [ 3  1 -6  0  3  1 -5  2  2  8 -4  4 -3  2 -6 -7]
        [ 2  6 -4  4  0 -1  7  0 -6  3  9  1 -3 -1  4  3]

    We construct `\ZZ^n`::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: IntegerLattice(ZZ^10)
        Free module of degree 10 and rank 10 over Integer Ring
        User basis matrix:
        [1 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 1]


    Sage also interfaces with fpLLL's lattice generator::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: from sage.libs.fplll.fplll import gen_simdioph
        sage: IntegerLattice(gen_simdioph(8, 20, 10), lll_reduce=False)
        Free module of degree 8 and rank 8 over Integer Ring
        User basis matrix:
        [   1024  829556  161099   11567  521155  769480  639201  689979]
        [      0 1048576       0       0       0       0       0       0]
        [      0       0 1048576       0       0       0       0       0]
        [      0       0       0 1048576       0       0       0       0]
        [      0       0       0       0 1048576       0       0       0]
        [      0       0       0       0       0 1048576       0       0]
        [      0       0       0       0       0       0 1048576       0]
        [      0       0       0       0       0       0       0 1048576]

    """

    if isinstance(basis, OrderElement_absolute):
        basis = basis.matrix()
    elif isinstance(basis, FreeModule_ambient_pid):
        basis = basis.basis_matrix()

    try:
        basis = matrix(ZZ, basis)
    except TypeError:
        raise NotImplementedError("only integer lattices supported")

    return FreeModule_submodule_with_basis_integer(ZZ**basis.ncols(),
                                                   basis=basis,
                                                   lll_reduce=lll_reduce)

class FreeModule_submodule_with_basis_integer(FreeModule_submodule_with_basis_pid):
    r"""
    This class represents submodules of `\ZZ^n` with a distinguished basis.

    However, most functionality in excess of standard submodules over PID
    is for these submodules considered as discrete subgroups of `\ZZ^n`, i.e.
    as lattices. That is, this class provides functions for computing LLL
    and BKZ reduced bases for this free module with respect to the standard
    Euclidean norm.

    EXAMPLE::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: L = IntegerLattice(sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)); L
        Free module of degree 10 and rank 10 over Integer Ring
        User basis matrix:
        [ 0  1  2  0  1  2 -1  0 -1 -1]
        [ 0  1  0 -3  0  0  0  0  3 -1]
        [ 1  1 -1  0 -3  0  0  1  2 -2]
        [-1  2 -1 -1  2 -2  1 -1  0 -1]
        [ 1  0 -4  2  0  1 -2 -1  0  0]
        [ 2  3  0  1  1  0 -2  3  0  0]
        [-2 -3 -2  0  0  1 -1  1  3 -2]
        [-3  0 -1  0 -2 -1 -2  1 -1  1]
        [ 1  4 -1  1  2  2  1  0  3  1]
        [-1 -1  0 -3 -1  2  2  3 -1  0]
        sage: L.shortest_vector()
        (0, 1, 2, 0, 1, 2, -1, 0, -1, -1)

    """
    def __init__(self, ambient, basis, check=True, echelonize=False,
                 echelonized_basis=None, already_echelonized=False,
                 lll_reduce=True):
        r"""
        Construct a new submodule of `\ZZ^n` with a distinguished basis.

        INPUT:

        - ``ambient`` -- ambient free module over a principal ideal domain
          `\ZZ`, i.e. `\ZZ^n`

        - ``basis`` -- either a list of vectors or a matrix over the integers

        - ``check`` -- (default: ``True``) if ``False``, correctness of
          the input will not be checked and type conversion may be omitted,
          use with care

        - ``echelonize`` -- (default:``False``) if ``True``, ``basis`` will be
          echelonized and the result will be used as the default basis of the
          constructed submodule

        - `` echelonized_basis`` -- (default: ``None``) if not ``None``, must
          be the echelonized basis spanning the same submodule as ``basis``

        - ``already_echelonized`` -- (default: ``False``) if ``True``,
          ``basis`` must be already given in the echelonized form

        - ``lll_reduce`` -- (default: ``True``) run LLL reduction on the basis
          on construction

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: IntegerLattice([[1,0,-2], [0,2,5], [0,0,7]])
            Free module of degree 3 and rank 3 over Integer Ring
            User basis matrix:
            [ 1  0 -2]
            [ 1 -2  0]
            [ 2  2  1]

            sage: IntegerLattice(random_matrix(ZZ, 5, 5, x=-2^20, y=2^20))
            Free module of degree 5 and rank 5 over Integer Ring
            User basis matrix:
            [  -7945 -381123   85872 -225065   12924]
            [-158254  120252  189195 -262144 -345323]
            [ 232388  -49556  306585  -31340  401528]
            [-353460  213748  310673  158140  172810]
            [-287787  333937 -145713 -482137  186529]

            sage: K.<a> = NumberField(x^8+1)
            sage: O = K.ring_of_integers()
            sage: f = O.random_element(); f
            a^7 - a^6 + 4*a^5 - a^4 + a^3 + 1

            sage: IntegerLattice(f)
            Free module of degree 8 and rank 8 over Integer Ring
            User basis matrix:
            [ 0  1  0  1  0  3  3  0]
            [ 1  0  0  1 -1  4 -1  1]
            [ 0  0  1  0  1  0  3  3]
            [-4  1 -1  1  0  0  1 -1]
            [ 1 -3  0  0  0  3  0 -2]
            [ 0 -1  1 -4  1 -1  1  0]
            [ 2  0 -3 -1  0 -3  0  0]
            [-1  0 -1  0 -3 -3  0  0]

        """
        basis = matrix(ZZ, basis)
        self._basis_is_LLL_reduced = False

        if lll_reduce:
            basis = matrix([v for v in basis.LLL() if v])
            self._basis_is_LLL_reduced = True

        basis.set_immutable()
        FreeModule_submodule_with_basis_pid.__init__(self,
                                                     ambient=ambient,
                                                     basis=basis,
                                                     check=check,
                                                     echelonize=echelonize,
                                                     echelonized_basis=echelonized_basis,
                                                     already_echelonized=already_echelonized)

        self._reduced_basis = basis.change_ring(ZZ)

    @property
    def reduced_basis(self):
        """
        This attribute caches the currently best known reduced basis for
        ``self``, where "best" is defined by the Euclidean norm of the
        first row vector.

        EXAMPLE::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: L = IntegerLattice(random_matrix(ZZ, 10, 10), lll_reduce=False)
            sage: L.reduced_basis
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: _ = L.LLL()
            sage: L.reduced_basis
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

        """
        return self._reduced_basis

    def LLL(self, *args, **kwds):
        r"""
        Return an LLL reduced basis for ``self``.

        A lattice basis `(b_1, b_2, ..., b_d)` is `(\delta, \eta)`-LLL-reduced
        if the two following conditions hold:

        -  For any `i > j`, we have `\lvert \mu_{i, j} \rvert \leq η`.

        -  For any `i < d`, we have
           `\delta \lvert b_i^* \rvert^2 \leq \lvert b_{i+1}^* +
           \mu_{i+1, i} b_i^* \rvert^2`,

        where `\mu_{i,j} = \langle b_i, b_j^* \rangle / \langle b_j^*,b_j^*
        \rangle` and `b_i^*` is the `i`-th vector of the Gram-Schmidt
        orthogonalisation of `(b_1, b_2, \ldots, b_d)`.

        The default reduction parameters are `\delta = 3/4` and
        `\eta = 0.501`.

        The parameters `\delta` and `\eta` must satisfy:
        `0.25 < \delta \leq 1.0` and `0.5 \leq \eta < \sqrt{\delta}`.
        Polynomial time complexity is only guaranteed for `\delta < 1`.

        INPUT:

        - ``*args`` -- passed through to
          :meth:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL`

        - ``**kwds`` -- passed through to
          :meth:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL`

        OUTPUT:

        An integer matrix which is an LLL-reduced basis for this lattice.

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: A = random_matrix(ZZ, 10, 10, x=-2000, y=2000)
            sage: L = IntegerLattice(A, lll_reduce=False); L
            Free module of degree 10 and rank 10 over Integer Ring
            User basis matrix:
            [ -645 -1037 -1775 -1619  1721 -1434  1766  1701  1669  1534]
            [ 1303   960  1998 -1838  1683 -1332   149   327  -849 -1562]
            [-1113 -1366  1379   669    54  1214 -1750  -605 -1566  1626]
            [-1367  1651   926  1731  -913   627   669 -1437  -132  1712]
            [ -549  1327 -1353    68  1479 -1803  -456  1090  -606  -317]
            [ -221 -1920 -1361  1695  1139   111 -1792  1925  -656  1992]
            [-1934   -29    88   890  1859  1820 -1912 -1614 -1724  1606]
            [ -590 -1380  1768   774   656   760  -746  -849  1977 -1576]
            [  312  -242 -1732  1594  -439 -1069   458 -1195  1715    35]
            [  391  1229 -1815   607  -413  -860  1408  1656  1651  -628]
            sage: min(v.norm().n() for v in L.reduced_basis)
            3346.57...

            sage: L.LLL()
            [ -888    53  -274   243   -19   431   710   -83   928   347]
            [  448  -330   370  -511   242  -584    -8  1220   502   183]
            [ -524  -460   402  1338  -247  -279 -1038   -28  -159  -794]
            [  166  -190  -162  1033  -340   -77 -1052  1134  -843   651]
            [  -47 -1394  1076  -132   854  -151   297  -396  -580  -220]
            [-1064   373  -706   601  -587 -1394   424   796   -22  -133]
            [-1126   398   565 -1418  -446  -890  -237  -378   252   247]
            [ -339   799   295   800   425  -605  -730 -1160   808   666]
            [  755 -1206  -918  -192 -1063   -37  -525   -75   338   400]
            [  382  -199 -1839  -482   984   -15  -695   136   682   563]
            sage: L.reduced_basis[0].norm().n()
            1613.74...

        """
        basis = self.reduced_basis
        basis = [v for v in basis.LLL(*args, **kwds) if v]
        basis = matrix(ZZ, len(basis), len(basis[0]), basis)
        basis.set_immutable()

        if self.reduced_basis[0].norm() > basis[0].norm():
            self._reduced_basis = basis
        return basis

    def BKZ(self, *args, **kwds):
        """
        Return a Block Korkine-Zolotareff reduced basis for ``self``.

        INPUT:

        - ``*args`` -- passed through to
          :meth:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.BKZ`

        - ``*kwds`` -- passed through to
          :meth:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.BKZ`

        OUTPUT:

        An integer matrix which is a BKZ-reduced basis for this lattice.

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=60, q=2^60, seed=42)
            sage: L = IntegerLattice(A, lll_reduce=False)
            sage: min(v.norm().n() for v in L.reduced_basis)
            4.17330740711759e15

            sage: L.LLL()
            60 x 60 dense matrix over Integer Ring (use the '.str()' method to see the entries)

            sage: min(v.norm().n() for v in L.reduced_basis)
            5.19615242270663

            sage: L.BKZ(block_size=10)
            60 x 60 dense matrix over Integer Ring (use the '.str()' method to see the entries)

            sage: min(v.norm().n() for v in L.reduced_basis)
            4.12310562561766

        .. NOTE::

            If ``block_size == L.rank()`` where ``L`` is this lattice, then
            this function performs Hermite-Korkine-Zolotareff (HKZ) reduction.
        """
        basis = self.reduced_basis
        basis = [v for v in basis.BKZ(*args, **kwds) if v]
        basis = matrix(ZZ, len(basis), len(basis[0]), basis)
        basis.set_immutable()

        if self.reduced_basis[0].norm() > basis[0].norm():
            self._reduced_basis = basis
        return basis

    def HKZ(self, *args, **kwds):
        r"""
        Hermite-Korkine-Zolotarev (HKZ) reduce the basis.

        A basis `B` of a lattice `L`, with orthogonalized basis `B^*` such
        that `B = M \cdot B^*` is HKZ reduced, if and only if, the following
        properties are satisfied:

        #. The basis `B` is size-reduced, i.e., all off-diagonal
           coefficients of `M` satisfy `|\mu_{i,j}| \leq 1/2`

        #. The vector `b_1` realizes the first minimum `\lambda_1(L)`.

        #. The projection of the vectors `b_2, \ldots,b_r` orthogonally to
           `b_1` form an HKZ reduced basis.

        .. NOTE::

            This is realised by calling
            :func:`sage.modules.free_module_integer.FreeModule_submodule_with_basis_integer.BKZ` with
            ``block_size == self.rank()``.

        INPUT:

        - ``*args`` -- passed through to :meth:`BKZ`

        - ``*kwds`` -- passed through to :meth:`BKZ`

        OUTPUT:

        An integer matrix which is a HKZ-reduced basis for this lattice.

        EXAMPLE::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: L = sage.crypto.gen_lattice(type='random', n=1, m=40, q=2^60, seed=42, lattice=True)
            sage: L.HKZ()
            40 x 40 dense matrix over Integer Ring (use the '.str()' method to see the entries)

            sage: L.reduced_basis[0]
            (-1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0,
             -1, 1, 0, 1, 1, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        """
        return self.BKZ(block_size=self.rank())


    @cached_method
    def volume(self):
        r"""
        Return `vol(L)` which is `\sqrt{\det(B \cdot B^T)}` for any basis `B`.

        OUTPUT:

        An integer.

        EXAMPLE::

            sage: L = sage.crypto.gen_lattice(m=10, seed=42, lattice=True)
            sage: L.volume()
            14641
        """
        if self.rank() == self.degree():
            return abs(self.reduced_basis.determinant())
        else:
            return self.gram_matrix().determinant().sqrt()

    @cached_method
    def discriminant(self):
        r"""
        Return `|\det(G)|`, i.e. the absolute value of the determinant of the
        Gram matrix `B \cdot B^T` for any basis `B`.

        OUTPUT:

        An integer.

        EXAMPLE::

            sage: L = sage.crypto.gen_lattice(m=10, seed=42, lattice=True)
            sage: L.discriminant()
            214358881
        """
        return abs(self.gram_matrix().determinant())

    @cached_method
    def is_unimodular(self):
        """
        Return ``True`` if this lattice is unimodular.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: L = IntegerLattice([[1, 0], [0, 1]])
            sage: L.is_unimodular()
            True
            sage: IntegerLattice([[2, 0], [0, 3]]).is_unimodular()
            False
        """
        return self.volume() == 1

    @cached_method
    def shortest_vector(self, update_reduced_basis=True, algorithm="fplll", *args, **kwds):
        r"""
        Return a shortest vector.

        INPUT:

        - ``update_reduced_basis`` -- (default: ``True``) set this flag if
          the found vector should be used to improve the basis 

        - ``algorithm`` -- (default: ``"fplll"``) either ``"fplll"`` or
          ``"pari"``

        - ``*args`` -- passed through to underlying implementation

        - ``**kwds`` -- passed through to underlying implementation

        OUTPUT:

        A shortest non-zero vector for this lattice.

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=30, q=2^40, seed=42)
            sage: L = IntegerLattice(A, lll_reduce=False)
            sage: min(v.norm().n() for v in L.reduced_basis)
            6.03890756700000e10

            sage: L.shortest_vector().norm().n()
            3.74165738677394

            sage: L = IntegerLattice(A, lll_reduce=False)
            sage: min(v.norm().n() for v in L.reduced_basis)
            6.03890756700000e10

            sage: L.shortest_vector(algorithm="pari").norm().n()
            3.74165738677394

            sage: L = IntegerLattice(A, lll_reduce=True)
            sage: L.shortest_vector(algorithm="pari").norm().n()
            3.74165738677394
        """
        if algorithm == "pari":
            if self._basis_is_LLL_reduced:
                B = self.basis_matrix().change_ring(ZZ)
                qf = self.gram_matrix()
            else:
                B = self.reduced_basis.LLL()
                qf = B*B.transpose()

            count, length, vectors = qf._pari_().qfminim()
            v = vectors.python().columns()[0]
            w = v*B
        elif algorithm == "fplll":
            from sage.libs.fplll.fplll import FP_LLL
            L = FP_LLL(self.reduced_basis)
            w = L.shortest_vector(*args, **kwds)
        else:
            raise ValueError("algorithm '{}' unknown".format(algorithm))

        if update_reduced_basis:
            self.update_reduced_basis(w)
        return w

    def update_reduced_basis(self, w):
        """
        Inject the vector ``w`` and run LLL to update the basis.

        INPUT:

        - ``w`` -- a vector

        OUTPUT:

        Nothing is returned but the internal state is modified.

        EXAMPLE::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=30, q=2^40, seed=42)
            sage: L = IntegerLattice(A)
            sage: B = L.reduced_basis
            sage: v = L.shortest_vector(update_reduced_basis=False)
            sage: L.update_reduced_basis(v)
            sage: bool(L.reduced_basis[0].norm() < B[0].norm())
            True
        """
        w = matrix(ZZ, w)
        L = w.stack(self.reduced_basis).LLL()
        assert(L[0] == 0)
        self._reduced_basis = L.matrix_from_rows(range(1,L.nrows()))


    @cached_method
    def voronoi_cell(self, radius=None):
        """
        Compute the Voronoi cell of a lattice, returning a Polyhedron.

        INPUT:

        - ``radius`` -- (default: automatic determination) radius of ball
          containing considered vertices

        OUTPUT:

        The Voronoi cell as a Polyhedron instance.

        The result is cached so that subsequent calls to this function
        return instantly.

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: L = IntegerLattice([[1, 0], [0, 1]])
            sage: V = L.voronoi_cell()
            sage: V.Vrepresentation()
            (A vertex at (1/2, -1/2), A vertex at (1/2, 1/2), A vertex at (-1/2, 1/2), A vertex at (-1/2, -1/2))

        The volume of the Voronoi cell is the square root of the
        discriminant of the lattice::

            sage: L = IntegerLattice(Matrix(ZZ, 4, 4, [[0,0,1,-1],[1,-1,2,1],[-6,0,3,3,],[-6,-24,-6,-5]])); L
            Free module of degree 4 and rank 4 over Integer Ring
            User basis matrix:
            [  0   0   1  -1]
            [  1  -1   2   1]
            [ -6   0   3   3]
            [ -6 -24  -6  -5]
            sage: V = L.voronoi_cell() # long time
            sage: V.volume()           # long time
            678
            sage: sqrt(L.discriminant())
            678

        Lattices not having full dimension are handled as well::

            sage: L = IntegerLattice([[2, 0, 0], [0, 2, 0]])
            sage: V = L.voronoi_cell()
            sage: V.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, -1, 0) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0)

        ALGORITHM:

        Uses parts of the algorithm from [Vit1996]_.

        REFERENCES:

        .. [Vit1996] \E. Viterbo, E. Biglieri. *Computing the Voronoi Cell
           of a Lattice: The Diamond-Cutting Algorithm*.
           IEEE Transactions on Information Theory, 1996.
        """
        if not self._basis_is_LLL_reduced:
            self.LLL()

        B = self.reduced_basis

        from .diamond_cutting import calculate_voronoi_cell
        return calculate_voronoi_cell(B, radius=radius)

    def voronoi_relevant_vectors(self):
        """
        Compute the embedded vectors inducing the Voronoi cell.

        OUTPUT:

        The list of Voronoi relevant vectors.

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: L = IntegerLattice([[3, 0], [4, 0]])
            sage: L.voronoi_relevant_vectors()
            [(-1, 0), (1, 0)]
        """
        V = self.voronoi_cell()

        def defining_point(ieq):
            """
            Compute the point defining an inequality.

            INPUT:

            - ``ieq`` -- an inequality in the form [c, a1, a2, ...]
              meaning a1 * x1 + a2 * x2 + ... ≦ c

            OUTPUT:

            The point orthogonal to the hyperplane defined by ``ieq``
            in twice the distance from the origin.
            """
            c = ieq[0]
            a = ieq[1:]
            n = sum(y ** 2 for y in a)
            return vector([2 * y * c / n for y in a])

        return [defining_point(ieq) for ieq in V.inequality_generator()]

    def closest_vector(self, t):
        """
        Compute the closest vector in the embedded lattice to a given vector.

        INPUT:

        - ``t`` -- the target vector to compute the closest vector to

        OUTPUT:

        The vector in the lattice closest to ``t``.

        EXAMPLES::

            sage: from sage.modules.free_module_integer import IntegerLattice
            sage: L = IntegerLattice([[1, 0], [0, 1]])
            sage: L.closest_vector((-6, 5/3))
            (-6, 2)

        ALGORITHM:

        Uses the algorithm from [Mic2010]_.

        REFERENCES:

        .. [Mic2010] \D. Micciancio, P. Voulgaris. *A Deterministic Single
           Exponential Time Algorithm for Most Lattice Problems based on
           Voronoi Cell Computations*.
           Proceedings of the 42nd ACM Symposium Theory of Computation, 2010.
        """
        voronoi_cell = self.voronoi_cell()

        def projection(M, v):
            Mt = M.transpose()
            P = Mt * (M * Mt) ** (-1) * M
            return P * v

        t = projection(matrix(self.reduced_basis), vector(t))

        def CVPP_2V(t, V, voronoi_cell):
            t_new = t
            while not voronoi_cell.contains(t_new.list()):
                v = max(V, key=lambda v: t_new * v / v.norm() ** 2)
                t_new = t_new - v
            return t - t_new

        V = self.voronoi_relevant_vectors()
        t = vector(t)
        p = 0
        while not (ZZ(2 ** p) * voronoi_cell).contains(t):
            p += 1
        t_new = t
        i = p
        while i >= 1:
            V_scaled = [v * (2 ** (i - 1)) for v in V]
            t_new = t_new - CVPP_2V(t_new, V_scaled, ZZ(2 ** (i - 1)) * voronoi_cell)
            i -= 1
        return t - t_new

