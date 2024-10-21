# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
r"""
Quotients of the Bruhat-Tits tree

This package contains all the functionality described and developed in [FM2014]_.
It allows for computations with fundamental domains of the Bruhat-Tits tree,
under the action of arithmetic groups arising from units in definite
quaternion algebras.

EXAMPLES:

Create the quotient attached to a maximal order of the quaternion algebra of
discriminant `13`, at the prime `p = 5`::

    sage: Y = BruhatTitsQuotient(5, 13)

We can query for its genus, as well as get it back as a graph::

    sage: Y.genus()
    5
    sage: Y.get_graph()
    Multi-graph on 2 vertices

The rest of functionality can be found in the docstrings below.

AUTHORS:

- Cameron Franc and Marc Masdeu (2011): initial version
"""

# ****************************************************************************
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from collections import deque

from sage.arith.misc import gcd, xgcd, kronecker_symbol, fundamental_discriminant
from sage.matrix.constructor import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.misc.verbose import verbose
from sage.modular.arithgroup.all import Gamma0
from sage.modular.arithgroup.congroup_gammaH import GammaH_constructor
from sage.modular.dirichlet import DirichletGroup
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
lazy_import("sage.plot.colors", "rainbow")

lazy_import('sage.algebras.quatalg.quaternion_algebra', 'QuaternionAlgebra')
lazy_import('sage.graphs.graph', 'Graph')
lazy_import('sage.libs.pari.all', 'pari')
lazy_import('sage.plot.colors', 'rainbow')
lazy_import('sage.rings.number_field.number_field', 'NumberField')
lazy_import('sage.rings.padics.factory', ['Qp', 'Zp'])


class DoubleCosetReduction(SageObject):
    r"""
    Edges in the Bruhat-Tits tree are represented by cosets of
    matrices in `GL_2`. Given a matrix `x` in `GL_2`, this
    class computes and stores the data corresponding to the
    double coset representation of `x` in terms of a fundamental
    domain of edges for the action of the arithmetic group `\Gamma`.

    More precisely:

    Initialized with an element `x` of `GL_2(\ZZ)`, finds elements
    `\gamma` in `\Gamma`, `t` and an edge `e` such that `get=x`. It
    stores these values as members ``gamma``, ``label`` and functions
    ``self.sign()``,  ``self.t()`` and ``self.igamma()``, satisfying:

    - if ``self.sign() == +1``:
      ``igamma() * edge_list[label].rep * t() == x``

    - if ``self.sign() == -1``:
      ``igamma() * edge_list[label].opposite.rep * t() == x``

    It also stores a member called power so that:

        ``p**(2*power) = gamma.reduced_norm()``

    The usual decomposition `get=x` would be:

    - g = gamma / (p ** power)

    - e = edge_list[label]

    - t' = t * p ** power

    Here usual denotes that we have rescaled gamma to have unit
    determinant, and so that the result is honestly an element
    of the arithmetic quaternion group under consideration. In
    practice we store integral multiples and keep track of the
    powers of `p`.

    INPUT:

    - ``Y`` -- BruhatTitsQuotient object in which to work
    - ``x`` -- Something coercible into a matrix in `GL_2(\ZZ)`. In
      principle we should allow elements in `GL_2(\QQ_p)`, but it is
      enough to work with integral entries
    - ``extrapow`` -- gets added to the power attribute, and it is
      used for the Hecke action

    EXAMPLES::

        sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
        sage: Y = BruhatTitsQuotient(5, 13)
        sage: x = Matrix(ZZ,2,2,[123,153,1231,1231])
        sage: d = DoubleCosetReduction(Y,x)
        sage: d.sign()
        -1
        sage: d.igamma()*Y._edge_list[d.label - len(Y.get_edge_list())].opposite.rep*d.t() == x
        True
        sage: x = Matrix(ZZ,2,2,[1423,113553,11231,12313])
        sage: d = DoubleCosetReduction(Y,x)
        sage: d.sign()
        1
        sage: d.igamma()*Y._edge_list[d.label].rep*d.t() == x
        True

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """

    def __init__(self, Y, x, extrapow=0):
        r"""
        Initialize and compute the reduction as a double coset.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
            sage: Y = BruhatTitsQuotient(5, 13)
            sage: x = Matrix(ZZ,2,2,[123,153,1231,1231])
            sage: d = DoubleCosetReduction(Y,x)
            sage: TestSuite(d).run()
        """
        e1 = Y._BT.edge(x)
        try:
            g, label, parity = Y._cached_decomps[e1]
        except KeyError:
            valuation = e1.determinant().valuation(Y._p)
            parity = valuation % 2
            v1 = Y._BT.target(e1)
            v = Y.fundom_rep(v1)
            g, e = Y._find_equivalent_edge(e1, v.entering_edges,
                                           valuation=valuation)
            label = e.label
            Y._cached_decomps[e1] = (g, label, parity)

        self._parent = Y
        self.parity = parity
        self._num_edges = len(Y.get_edge_list())
        self.label = label + parity * self._num_edges
        # The label will encode whether it is an edge or its opposite !
        self.gamma = g[0]
        self.x = x
        self.power = g[1] + extrapow
        self._t_prec = -1
        self._igamma_prec = -1

    def _repr_(self):
        r"""
        Return the representation of ``self`` as a string.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
            sage: Y = BruhatTitsQuotient(5, 13)
            sage: x = Matrix(ZZ,2,2,[123,153,1231,1231])
            sage: DoubleCosetReduction(Y,x)
            Double coset data (-1, [(4), (5), (-4), (-4)], 8)
        """
        return "Double coset data (%s, %s, %s)" % (self.sign(),
                                                   list(self.gamma), self.label)

    def __eq__(self, other):
        """
        Return ``self == other``.

        TESTS::

            sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
            sage: Y = BruhatTitsQuotient(5, 13)
            sage: x = Matrix(ZZ,2,2,[123,153,1231,1231])
            sage: d1 = DoubleCosetReduction(Y,x)
            sage: d1 == d1
            True
        """
        if self._parent != other._parent:
            return False
        if self.parity != other.parity:
            return False
        if self._num_edges != other._num_edges:
            return False
        if self.label != other.label:
            return False
        if self.gamma != other.gamma:
            return False
        if self.x != other.x:
            return False
        if self.power != other.power:
            return False
        if self._t_prec != other._t_prec:
            return False
        if self._igamma_prec != other._igamma_prec:
            return False
        return True

    def __ne__(self, other):
        """
        Return ``self != other``.

        TESTS::

            sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
            sage: Y = BruhatTitsQuotient(5, 13)
            sage: x = Matrix(ZZ,2,2,[123,153,1231,1231])
            sage: d1 = DoubleCosetReduction(Y,x)
            sage: d1 != d1
            False
        """
        return not self.__eq__(other)

    def sign(self):
        r"""
        Return the direction of the edge.

        The Bruhat-Tits quotients are directed graphs but we only store
        half the edges (we treat them more like unordered graphs).
        The sign tells whether the matrix self.x is equivalent to the
        representative in the quotient (sign = +1), or to the
        opposite of one of the representatives (sign = -1).

        OUTPUT:

        an int that is +1 or -1 according to the sign of ``self``

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
            sage: Y = BruhatTitsQuotient(3, 11)
            sage: x = Matrix(ZZ,2,2,[123,153,1231,1231])
            sage: d = DoubleCosetReduction(Y,x)
            sage: d.sign()
            -1
            sage: d.igamma()*Y._edge_list[d.label - len(Y.get_edge_list())].opposite.rep*d.t() == x
            True
            sage: x = Matrix(ZZ,2,2,[1423,113553,11231,12313])
            sage: d = DoubleCosetReduction(Y,x)
            sage: d.sign()
            1
            sage: d.igamma()*Y._edge_list[d.label].rep*d.t() == x
            True
        """
        if self.parity == 0:
            return 1
        else:
            return -1

    def igamma(self, embedding=None, scale=1):
        r"""
        Image under gamma.

        Elements of the arithmetic group can be regarded as elements
        of the global quaternion order, and hence may be represented
        exactly. This function computes the image of such an element
        under the local splitting and returns the corresponding `p`-adic
        approximation.

        INPUT:

        - ``embedding`` -- integer; or a function (default:
          none). If ``embedding`` is None, then the image of
          ``self.gamma`` under the local splitting associated to
          ``self.Y`` is used. If ``embedding`` is an integer, then
          the precision of the local splitting of self.Y is raised
          (if necessary) to be larger than this integer, and this
          new local splitting is used. If a function is passed, then
          map ``self.gamma`` under ``embedding``.
        - ``scale`` -- (default: 1) scaling factor applied to the output

        OUTPUT:

        a 2x2 matrix with `p`-adic entries encoding the image of ``self``
        under the local splitting

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
            sage: Y = BruhatTitsQuotient(7, 11)
            sage: d = DoubleCosetReduction(Y,Matrix(ZZ,2,2,[123,45,88,1]))
            sage: d.igamma()
            [6 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + O(7^5)                                   O(7^5)]
            [                                  O(7^5) 6 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + O(7^5)]
            sage: d.igamma(embedding = 7)
            [6 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + 6*7^5 + 6*7^6 + O(7^7)                                                   O(7^7)]
            [                                                  O(7^7) 6 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + 6*7^5 + 6*7^6 + O(7^7)]
        """
        Y = self._parent
        if embedding is None:
            prec = Y._prec
        else:
            try:
                # The user wants higher precision
                prec = ZZ(embedding)
            except TypeError:
                # The user knows what she is doing, so let it go
                return embedding(self.gamma)
        if prec > self._igamma_prec:
            self._igamma_prec = prec
            self._cached_igamma = Y.embed_quaternion(self.gamma, exact=False,
                                                     prec=prec)
        return scale * self._cached_igamma

    def t(self, prec=None):
        r"""
        Return the 't part' of the decomposition using the rest of the data.

        INPUT:

        - ``prec`` -- a `p`-adic precision that t will be computed
          to. Defaults to the default working precision of self.

        OUTPUT:

        a 2x2 `p`-adic matrix with entries of
        precision ``prec`` that is the 't-part' of the decomposition of
        self

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import DoubleCosetReduction
            sage: Y = BruhatTitsQuotient(5, 13)
            sage: x = Matrix(ZZ,2,2,[123,153,1231,1232])
            sage: d = DoubleCosetReduction(Y,x)
            sage: t = d.t(20)
            sage: t[1,0].valuation() > 0
            True
        """
        Y = self._parent
        if prec is None:
            prec = max([5, Y._prec])
        if self._t_prec >= prec:
            return self._cached_t
        e = Y._edge_list[self.label % self._num_edges]
        tmp_prec = prec
        while self._t_prec < prec:
            if self.parity == 0:
                self._cached_t = (self.igamma(tmp_prec) * e.rep).inverse() * self.x
                # assert self._cached_t[1, 0].valuation()>self._cached_t[1,1].valuation()
            else:
                self._cached_t = (self.igamma(tmp_prec) * e.opposite.rep).inverse() * self.x
                # assert self._cached_t[1, 0].valuation()>self._cached_t[1,1].valuation()
            tmp_prec += 1
            self._t_prec = min([xx.precision_absolute()
                                for xx in self._cached_t.list()])
        return self._cached_t


class BruhatTitsTree(SageObject, UniqueRepresentation):
    r"""
    An implementation of the Bruhat-Tits tree for `GL_2(\QQ_p)`.

    INPUT:

    - ``p`` -- a prime number. The corresponding tree is then `p+1` regular

    EXAMPLES:

    We create the tree for `GL_2(\QQ_5)`::

        sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
        sage: p = 5
        sage: T = BruhatTitsTree(p)
        sage: m = Matrix(ZZ,2,2,[p**5,p**2,p**3,1+p+p*3])
        sage: e = T.edge(m); e
        [  0  25]
        [625  21]
        sage: v0 = T.origin(e); v0
        [ 25   0]
        [ 21 125]
        sage: v1 = T.target(e); v1
        [ 25   0]
        [ 21 625]
        sage: T.origin(T.opposite(e)) == v1
        True
        sage: T.target(T.opposite(e)) == v0
        True

    A value error is raised if a prime is not passed::

        sage: T = BruhatTitsTree(4)
        Traceback (most recent call last):
        ...
        ValueError: input (4) must be prime

    AUTHORS:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self, p):
        """
        Initialize a BruhatTitsTree object for a given prime `p`.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: T = BruhatTitsTree(17)
            sage: TestSuite(T).run()
        """
        if not ZZ(p).is_prime():
            raise ValueError(f'input ({p}) must be prime')
        self._p = ZZ(p)
        self._Mat_22 = MatrixSpace(ZZ, 2, 2)
        self._mat_p001 = self._Mat_22([self._p, 0, 0, 1])

    def target(self, e, normalized=False):
        r"""
        Return the target vertex of the edge represented by the
        input matrix e.

        INPUT:

        - ``e`` -- a 2x2 matrix with integer entries

        - ``normalized`` -- boolean (default: ``False``); if True
            then the input matrix is assumed to be normalized

        OUTPUT:

        - ``e`` -- 2x2 integer matrix representing the target of
          the input edge

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: T = BruhatTitsTree(7)
            sage: T.target(Matrix(ZZ,2,2,[1,5,8,9]))
            [1 0]
            [0 1]
        """
        if normalized:
            # then the normalized target vertex is also M and we save some
            # row reductions with a simple return
            return e
        else:
            # must normalize the target vertex representative
            return self.vertex(e)

    def origin(self, e, normalized=False):
        r"""
        Return the origin vertex of the edge represented by the
        input matrix e.

        INPUT:

        - ``e`` -- a 2x2 matrix with integer entries

        - ``normalized`` -- boolean (default: ``False``); if True
          then the input matrix M is assumed to be normalized

        OUTPUT: ``e`` -- 2x2 integer matrix

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: T = BruhatTitsTree(7)
            sage: T.origin(Matrix(ZZ,2,2,[1,5,8,9]))
            [1 0]
            [1 7]
        """
        if not normalized:
            # then normalize
            x = copy(self.edge(e))
        else:
            x = copy(e)
        x.swap_columns(0, 1)
        x.rescale_col(0, self._p)
        return self.vertex(x)

    def edge(self, M):
        r"""
        Normalize a matrix to the correct normalized edge
        representative.

        INPUT:

        - ``M`` -- 2x2 integer matrix

        OUTPUT: ``newM`` -- 2x2 integer matrix

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: T = BruhatTitsTree(3)
            sage: T.edge( Matrix(ZZ,2,2,[0,-1,3,0]) )
            [0 1]
            [3 0]
        """
        p = self._p
        # M_orig = M

        def lift(a):
            """
            Naively approximate a `p`-adic integer by a positive integer.

            INPUT:

            - ``a`` -- `p`-adic integer

            OUTPUT: integer

            EXAMPLES::

                sage: x = Zp(3)(-17)                                                    # needs sage.rings.padics
                sage: lift(x)                                                           # needs sage.rings.padics
                3486784384
            """
            try:
                return ZZ(a.lift())
            except AttributeError:
                return ZZ(a)

        if M.base_ring() is not ZZ:
            M = M.apply_map(lift, R=ZZ)

        v = min([M[i, j].valuation(p) for i in range(2) for j in range(2)])

        if v != 0:
            M = p ** (-v) * M

        det = M.determinant()
        if not det:
            raise NotImplementedError("matrix must be invertible")

        m00 = M[0, 0].valuation(p)
        m01 = M[0, 1].valuation(p)

        if m00 <= m01:
            tmp = det.valuation(p) - m00
            bigpower = p ** (1 + tmp)
            r = M[0, 0]
            if r != 0:
                r /= p ** m00
            g, s, _ = xgcd(r, bigpower)
            r = (M[1, 0] * s) % bigpower
            newM = self._Mat_22([p ** m00, 0, r, bigpower / p])
        else:
            tmp = det.valuation(p) - m01
            bigpower = p ** tmp
            r = M[0, 1]
            if r != 0:
                r /= p ** m01
            g, s, _ = xgcd(r, bigpower)
            r = (ZZ(M[1, 1]) * s) % bigpower
            newM = self._Mat_22([0, p ** m01, bigpower, r])
        newM.set_immutable()
        # assert self.is_in_group(M_orig.inverse()*newM, as_edge = True)
        return newM

    def vertex(self, M):
        r"""
        Normalize a matrix to the corresponding normalized
        vertex representative

        INPUT:

        - ``M`` -- 2x2 integer matrix

        OUTPUT: a 2x2 integer matrix

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 5
            sage: T = BruhatTitsTree(p)
            sage: m = Matrix(ZZ,2,2,[p**5,p**2,p**3,1+p+p*3])
            sage: e = T.edge(m)
            sage: t = m.inverse()*e
            sage: scaling = Qp(p,20)(t.determinant()).sqrt()
            sage: t = 1/scaling * t
            sage: min([t[ii,jj].valuation(p) for ii in range(2) for jj in range(2)]) >= 0
            True
            sage: t[1,0].valuation(p) > 0
            True
        """
        p = self._p

        def lift(a):
            try:
                return ZZ(a.lift())
            except AttributeError:
                return ZZ(a)

        if M.base_ring() is not ZZ:
            M = M.apply_map(lift, R=ZZ)

        v = min(M[i, j].valuation(p) for i in range(2) for j in range(2))

        if v:
            M = p ** (-v) * M
        m00 = M[0, 0].valuation(p)
        m01 = M[0, 1].valuation(p)
        if m01 < m00:
            M = copy(M)
            M.swap_columns(0, 1)
            m00 = m01
        tmp = M.determinant().valuation(p) - m00
        bigpower = p ** tmp
        r = M[0, 0]
        if r:
            r /= p ** m00
        # r = ZZ(r) % bigpower
        g, s, _ = xgcd(r, bigpower)
        m10 = M[1, 0] % bigpower
        r = (m10 * s) % bigpower
        newM = self._Mat_22([p ** m00, 0, r, bigpower])
        newM.set_immutable()
        # assert self.is_in_group(M_orig.inverse()*newM, as_edge=False)
        return newM

    def edges_leaving_origin(self):
        r"""
        Find normalized representatives for the `p+1` edges
        leaving the origin vertex corresponding to the homothety class
        of `\ZZ_p^2`. These are cached.

        OUTPUT: list of size `p+1` of 2x2 integer matrices

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: T = BruhatTitsTree(3)
            sage: T.edges_leaving_origin()
            [
            [0 1]  [3 0]  [0 1]  [0 1]
            [3 0], [0 1], [3 1], [3 2]
            ]
        """
        try:
            return self._edges_leaving_origin
        except AttributeError:
            p = self._p
            self._edges_leaving_origin = [self.edge(self._Mat_22([0, -1, p, 0]))]
            self._edges_leaving_origin.extend([self.edge(self._Mat_22([p, i, 0, 1])) for i in range(p)])
            return self._edges_leaving_origin

    def edge_between_vertices(self, v1, v2, normalized=False):
        r"""
        Compute the normalized matrix rep. for the edge
        passing between two vertices.

        INPUT:

        - ``v1`` -- 2x2 integer matrix

        - ``v2`` -- 2x2 integer matrix

        - ``normalized`` -- boolean (default: ``False``); whether the
          vertices are normalized

        OUTPUT:

        - 2x2 integer matrix, representing the edge from ``v1`` to
          ``v2``.  If ``v1`` and ``v2`` are not at distance `1`, raise
          a :exc:`ValueError`.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: v1 = T.vertex(Matrix(ZZ,2,2,[p,0,0,1])); v1
            [7 0]
            [0 1]
            sage: v2 = T.vertex(Matrix(ZZ,2,2,[p,1,0,1])); v2
            [1 0]
            [1 7]
            sage: T.edge_between_vertices(v1,v2)
            Traceback (most recent call last):
            ...
            ValueError: Vertices are not adjacent.

            sage: v3 = T.vertex(Matrix(ZZ,2,2,[1,0,0,1])); v3
            [1 0]
            [0 1]
            sage: T.edge_between_vertices(v1,v3)
            [0 1]
            [1 0]
        """
        if normalized:
            v22 = v2
        else:
            v22 = self.vertex(v2)
        for e in self.leaving_edges(v1):
            if self.target(e) == v22:
                return e
        raise ValueError('Vertices are not adjacent.')

    def leaving_edges(self, M):
        r"""
        Return edges leaving a vertex.

        INPUT:

        - ``M`` -- 2x2 integer matrix

        OUTPUT: list of size `p+1` of 2x2 integer matrices

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: T.leaving_edges(Matrix(ZZ,2,2,[1,0,0,1]))
            [
            [0 1]  [7 0]  [0 1]  [0 1]  [0 1]  [0 1]  [0 1]  [0 1]
            [7 0], [0 1], [7 1], [7 4], [7 5], [7 2], [7 3], [7 6]
            ]
        """
        return [self.edge(M * A) for A in self.edges_leaving_origin()]

    def opposite(self, e):
        r"""
        This function returns the edge oriented oppositely to a
        given edge.

        INPUT:

        - ``e`` -- 2x2 integer matrix

        OUTPUT: 2x2 integer matrix

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: e = Matrix(ZZ,2,2,[1,0,0,1])
            sage: T.opposite(e)
            [0 1]
            [7 0]
            sage: T.opposite(T.opposite(e)) == e
            True
        """
        x = copy(e)
        x.swap_columns(0, 1)
        x.rescale_col(0, self._p)
        return self.edge(x)

    def entering_edges(self, v):
        r"""
        This function returns the edges entering a given vertex.

        INPUT:

        - ``v`` -- 2x2 integer matrix

        OUTPUT: list of size `p+1` of 2x2 integer matrices

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: T.entering_edges(Matrix(ZZ,2,2,[1,0,0,1]))
            [
            [1 0]  [0 1]  [1 0]  [1 0]  [1 0]  [1 0]  [1 0]  [1 0]
            [0 1], [1 0], [1 1], [4 1], [5 1], [2 1], [3 1], [6 1]
            ]
        """
        return [self.opposite(e) for e in self.leaving_edges(v)]

    def subdivide(self, edgelist, level):
        r"""
        (Ordered) edges of ``self`` may be regarded as open balls in
        `P^1(\QQ_p)`.  Given a list of edges, this function return a list
        of edges corresponding to the level-th subdivision of the
        corresponding opens.  That is, each open ball of the input is
        broken up into `p^{\mbox{level}}` subballs of equal radius.

        INPUT:

        - ``edgelist`` -- list of edges

        - ``level`` -- integer

        OUTPUT: list of 2x2 integer matrices

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 3
            sage: T = BruhatTitsTree(p)
            sage: T.subdivide([Matrix(ZZ,2,2,[p,0,0,1])],2)
            [
            [27  0]  [0 9]  [0 9]  [0 3]  [0 3]  [0 3]  [0 3]  [0 3]  [0 3]
            [ 0  1], [3 1], [3 2], [9 1], [9 4], [9 7], [9 2], [9 5], [9 8]
            ]
        """
        if level < 0:
            return []
        if level == 0:
            return [self._Mat_22(edge) for edge in edgelist]
        else:
            newEgood = []
            for edge in edgelist:
                edge = self._Mat_22(edge)
                origin = self.origin(edge)
                newE = self.leaving_edges(self.target(edge))
                newEgood.extend([e for e in newE if self.target(e) != origin])
            return self.subdivide(newEgood, level - 1)

    def get_balls(self, center=1, level=1):
        r"""
        Return a decomposition of `P^1(\QQ_p)` into compact
        open balls.

        Each vertex in the Bruhat-Tits tree gives a decomposition of
        `P^1(\QQ_p)` into `p+1` open balls. Each of these balls may
        be further subdivided, to get a finer decomposition.

        This function returns the decomposition of `P^1(\QQ_p)`
        corresponding to ``center`` into `(p+1)p^{\mbox{level}}` balls.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 2
            sage: T = BruhatTitsTree(p)
            sage: T.get_balls(Matrix(ZZ,2,2,[p,0,0,1]),1)
            [
            [0 1]  [0 1]  [8 0]  [0 4]  [0 2]  [0 2]
            [2 0], [2 1], [0 1], [2 1], [4 1], [4 3]
            ]
        """
        return self.subdivide(self.leaving_edges(center), level)

    def find_path(self, v, boundary=None):
        r"""
        Compute a path from a vertex to a given set of so-called
        boundary vertices, whose interior must contain the origin
        vertex.  In the case that the boundary is not specified, it
        computes the geodesic between the given vertex and the origin.
        In the case that the boundary contains more than one vertex,
        it computes the geodesic to some point of the boundary.

        INPUT:

        - ``v`` -- a 2x2 matrix representing a vertex ``boundary``

        - a list of matrices (default: ``None``); if omitted, finds the
          geodesic from ``v`` to the central vertex

        OUTPUT:

        An ordered list of vertices describing the geodesic from
        ``v`` to ``boundary``, followed by the vertex in the boundary
        that is closest to ``v``.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 3
            sage: T = BruhatTitsTree(p)
            sage: T.find_path( Matrix(ZZ,2,2,[p^4,0,0,1]) )
            (
            [[81  0]
            [ 0  1], [27  0]
            [ 0  1], [9 0]
            [0 1], [3 0]      [1 0]
            [0 1]]          , [0 1]
            )
            sage: T.find_path( Matrix(ZZ,2,2,[p^3,0,134,p^2]) )
            (
            [[27  0]
            [ 8  9], [27  0]
            [ 2  3], [27  0]
            [ 0  1], [9 0]
            [0 1], [3 0]      [1 0]
            [0 1]]          , [0 1]
            )
        """
        if boundary is None:
            m = self._Mat_22(1)
            m.set_immutable()
            boundary = {m: m}
        m = self._mat_p001
        new_v = self.vertex(v)
        chain = []
        while new_v[1, 0] != 0 or new_v[0, 0].valuation(self._p) < new_v[1, 1].valuation(self._p):
            if new_v in boundary:
                return chain, boundary[new_v]
            chain.append(new_v)
            new_v = self.vertex(new_v * m)

        if new_v in boundary:
            return chain, boundary[new_v]

        while True:
            if new_v in boundary:
                return chain, boundary[new_v]
            chain.append(new_v)
            new_v = self._Mat_22([new_v[0, 0] / self._p, 0, 0, 1])
            new_v.set_immutable()
        raise RuntimeError

    def find_containing_affinoid(self, z):
        r"""
        Return the vertex corresponding to the affinoid in the
        `p`-adic upper half plane that a given (unramified!) point
        reduces to.

        INPUT:

        - ``z`` -- an element of an unramified extension of `\QQ_p`
          that is not contained in `\QQ_p`

        OUTPUT: a 2x2 integer matrix representing a vertex of ``self``

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: T = BruhatTitsTree(5)
            sage: K.<a> = Qq(5^2,20)
            sage: T.find_containing_affinoid(a)
            [1 0]
            [0 1]
            sage: z = 5*a+3
            sage: v = T.find_containing_affinoid(z).inverse(); v
            [   1    0]
            [-2/5  1/5]

        Note that the translate of ``z`` belongs to the standard
        affinoid. That is, it is a `p`-adic unit and its reduction
        modulo `p` is not in `\GF{p}`::

            sage: gz = (v[0,0]*z+v[0,1])/(v[1,0]*z+v[1,1]); gz                          # needs sage.rings.padics
            (a + 1) + O(5^19)
            sage: gz.valuation() == 0                                                   # needs sage.rings.padics
            True
        """
        # Assume z belongs to some extension of QQp.
        p = self._p
        if z.valuation() < 0:
            return self.vertex(self._Mat_22([0, 1, p, 0]) * self.find_containing_affinoid(1 / (p * z)))
        a = 0
        pn = 1
        val = z.valuation()
        L = [0 for _ in range(val)]

        L.extend(z.expansion())
        for n in range(len(L)):
            if L[n] != 0:
                if len(L[n]) > 1:
                    break
                if len(L[n]) > 0:
                    a += pn * L[n][0]
            pn *= p
        return self.vertex(self._Mat_22([pn, a, 0, 1]))

    def find_geodesic(self, v1, v2, normalized=True):
        r"""
        This function computes the geodesic between two vertices.

        INPUT:

        - ``v1`` -- 2x2 integer matrix representing a vertex

        - ``v2`` -- 2x2 integer matrix representing a vertex

        - ``normalized`` -- boolean (default: ``True``)

        OUTPUT:

        An ordered list of 2x2 integer matrices representing the vertices
        of the paths joining ``v1`` and ``v2``.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 3
            sage: T = BruhatTitsTree(p)
            sage: v1 = T.vertex( Matrix(ZZ,2,2,[p^3, 0, 1, p^1]) ); v1
            [27  0]
            [ 1  3]
            sage: v2 = T.vertex( Matrix(ZZ,2,2,[p,2,0,p]) ); v2
            [1 0]
            [6 9]
            sage: T.find_geodesic(v1,v2)
            [
            [27  0]  [27  0]  [9 0]  [3 0]  [1 0]  [1 0]  [1 0]
            [ 1  3], [ 0  1], [0 1], [0 1], [0 1], [0 3], [6 9]
            ]
        """
        if not normalized:
            v1, v2 = self.vertex(v1), self.vertex(v2)
        gamma = v2
        vv = self.vertex(gamma.adjugate() * v1)
        chain, v0 = self.find_path(vv)
        return [self.vertex(gamma * x) for x in chain + [v0]]

    def find_covering(self, z1, z2, level=0):
        r"""
        Compute a covering of `P^1(\QQ_p)` adapted to a certain
        geodesic in ``self``.

        More precisely, the `p`-adic upper half plane points ``z1``
        and ``z2`` reduce to vertices `v_1`, `v_2`.
        The returned covering consists of all the edges leaving the
        geodesic from `v_1` to `v_2`.

        INPUT:

        - ``z1``, ``z2`` -- unramified algebraic points of h_p

        OUTPUT: list of 2x2 integer matrices representing edges of self

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: from sage.modular.btquotients.btquotient import BruhatTitsTree
            sage: p = 3
            sage: K.<a> = Qq(p^2)
            sage: T = BruhatTitsTree(p)
            sage: z1 = a + a*p
            sage: z2 = 1 + a*p + a*p^2 - p^6
            sage: T.find_covering(z1,z2)
            [
            [0 1]  [3 0]  [0 1]  [0 1]  [0 1]  [0 1]
            [3 0], [0 1], [3 2], [9 1], [9 4], [9 7]
            ]

        .. NOTE::

            This function is used to compute certain Coleman integrals
            on `P^1`. That's why the input consists of two points of
            the `p`-adic upper half plane, but decomposes
            `P^1(\QQ_p)`. This decomposition is what allows us to
            represent the relevant integrand as a locally analytic
            function. The ``z1`` and ``z2`` appear in the integrand.
        """
        v1 = self.find_containing_affinoid(z1)
        v2 = self.find_containing_affinoid(z2)
        vertex_set = [self._Mat_22(0)]
        vertex_set += self.find_geodesic(v1, v2)
        vertex_set += [self._Mat_22(0)]
        E = []
        for ii in range(1, len(vertex_set) - 1):
            vv = vertex_set[ii]
            # m = vv.determinant().valuation(self._p)
            newE = self.leaving_edges(vv)
            for e in newE:
                targ = self.target(e)
                if targ != vertex_set[ii - 1] and targ != vertex_set[ii + 1]:
                    E.extend(self.subdivide([e], level))
        return E


class Vertex(SageObject):
    r"""
    This is a structure to represent vertices of quotients of the
    Bruhat-Tits tree.  It is useful to enrich the representation of
    the vertex as a matrix with extra data.

    INPUT:

    - ``p`` -- prime integer

    - ``label`` -- integer which uniquely identifies this vertex

    - ``rep`` -- a 2x2 matrix in reduced form representing this
      vertex

    - ``leaving_edges`` -- (default: empty list) list of edges
      leaving this vertex

    - ``entering_edges`` -- (default: empty list) list of edges
      entering this vertex

    - ``determinant`` -- (default: ``None``) the determinant of ``rep``,
      if known

    - ``valuation`` -- (default: ``None``) the valuation of the
      determinant of ``rep``, if known

    EXAMPLES::

        sage: from sage.modular.btquotients.btquotient import Vertex
        sage: v1 = Vertex(5,0,Matrix(ZZ,2,2,[1,2,3,18]))
        sage: v1.rep
        [ 1  2]
        [ 3 18]
        sage: v1.entering_edges
        []

    AUTHORS:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self, p, label, rep, leaving_edges=None,
                 entering_edges=None, determinant=None, valuation=None):
        """
        This initializes a structure to represent vertices of
        quotients of the Bruhat-Tits tree. It is useful to enrich the
        representation of the vertex as a matrix with extra data.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import Vertex
            sage: Y = BruhatTitsQuotient(5,13)
            sage: v1 = Vertex(5,0,Matrix(ZZ,2,2,[1,2,3,18]))
            sage: TestSuite(v1).run()
        """
        if leaving_edges is None:
            leaving_edges = []
        if entering_edges is None:
            entering_edges = []
        if determinant is None:
            determinant = rep.determinant()
        if valuation is None:
            valuation = determinant.valuation(p)
        self.p = p
        self.label = label
        self.rep = rep
        self.rep.set_immutable()
        self.determinant = determinant
        self.valuation = valuation
        self.parity = valuation % 2
        self.leaving_edges = leaving_edges
        self.entering_edges = entering_edges

    def _repr_(self):
        r"""
        Return the representation of ``self`` as a string.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,5)
            sage: X.get_vertex_list()[0]
            Vertex of Bruhat-Tits tree for p = 3
        """
        return "Vertex of Bruhat-Tits tree for p = %s" % (self.p)

    def __eq__(self, other):
        """
        Return ``self == other``.

        TESTS::

            sage: from sage.modular.btquotients.btquotient import Vertex
            sage: v1 = Vertex(7,0,Matrix(ZZ,2,2,[1,2,3,18]))
            sage: v1 == v1
            True
        """
        if self.p != other.p:
            return False
        if self.label != other.label:
            return False
        if self.rep != other.rep:
            return False
        if self.determinant != other.determinant:
            return False
        if self.valuation != other.valuation:
            return False
        if self.parity != other.parity:
            return False
        return True

    def __ne__(self, other):
        """
        Return ``self != other``.

        TESTS::

            sage: from sage.modular.btquotients.btquotient import Vertex
            sage: v1 = Vertex(7,0,Matrix(ZZ,2,2,[1,2,3,18]))
            sage: v1 != v1
            False
        """
        return not self.__eq__(other)


class Edge(SageObject):
    r"""
    This is a structure to represent edges of quotients of the
    Bruhat-Tits tree. It is useful to enrich the representation of an
    edge as a matrix with extra data.

    INPUT:

    - ``p`` -- prime integer

    - ``label`` -- integer which uniquely identifies this edge

    - ``rep`` -- a 2x2 matrix in reduced form representing this edge

    - ``origin`` -- the origin vertex of ``self``

    - ``target`` -- the target vertex of ``self``

    - ``links`` -- (default: empty list) list of elements of
      `\Gamma` which identify different edges in the Bruhat-Tits tree
      which are equivalent to ``self``

    - ``opposite`` -- (default: ``None``) the edge opposite to ``self``

    - ``determinant`` -- (default: ``None``) the determinant of ``rep``,
      if known

    - ``valuation`` -- (default: ``None``) the valuation of the
      determinant of ``rep``, if known

    EXAMPLES::

        sage: from sage.modular.btquotients.btquotient import Edge, Vertex
        sage: v1 = Vertex(7,0,Matrix(ZZ,2,2,[1,2,3,18]))
        sage: v2 = Vertex(7,0,Matrix(ZZ,2,2,[3,2,1,18]))
        sage: e1 = Edge(7,0,Matrix(ZZ,2,2,[1,2,3,18]),v1,v2)
        sage: e1.rep
        [ 1  2]
        [ 3 18]

    AUTHORS:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self, p, label, rep, origin, target, links=None,
                 opposite=None, determinant=None, valuation=None):
        """
        Representation for edges of quotients of the Bruhat-Tits
        tree. It is useful to enrich the representation of an edge as
        a matrix with extra data.

        EXAMPLES::

            sage: from sage.modular.btquotients.btquotient import Edge
            sage: Y = BruhatTitsQuotient(5,11)
            sage: el = Y.get_edge_list()
            sage: e1 = el.pop()
            sage: e2 = Edge(5,e1.label,e1.rep,e1.origin,e1.target)
            sage: TestSuite(e2).run()
        """
        if links is None:
            links = []
        if determinant is None:
            determinant = rep.determinant()
        if valuation is None:
            valuation = determinant.valuation(p)
        self.p = p
        self.label = label
        self.rep = rep
        self.rep.set_immutable()
        self.origin = origin
        self.target = target
        self.links = links
        self.opposite = opposite
        self.determinant = determinant
        self.valuation = valuation
        self.parity = valuation % 2

    def _repr_(self):
        r"""
        Return the representation of ``self`` as a string.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,5)
            sage: X.get_edge_list()[0]
            Edge of Bruhat-Tits tree for p = 3
        """
        return "Edge of Bruhat-Tits tree for p =  %s" % (self.p)

    def __eq__(self, other):
        """
        Return ``self == other``.

        TESTS::

            sage: from sage.modular.btquotients.btquotient import Edge,Vertex
            sage: v1 = Vertex(7,0,Matrix(ZZ,2,2,[1,2,3,18]))
            sage: v2 = Vertex(7,0,Matrix(ZZ,2,2,[3,2,1,18]))
            sage: e1 = Edge(7,0,Matrix(ZZ,2,2,[1,2,3,18]),v1,v2)
            sage: e1 == e1
            True
        """
        if self.p != other.p:
            return False
        if self.label != other.label:
            return False
        if self.rep != other.rep:
            return False
        if self.origin != other.origin:
            return False
        if self.target != other.target:
            return False
        if self.links != other.links:
            return False
        if self.opposite != other.opposite:
            return False
        if self.determinant != other.determinant:
            return False
        if self.valuation != other.valuation:
            return False
        if self.parity != other.parity:
            return False
        return True

    def __ne__(self, other):
        """
        Return ``self != other``.

        TESTS::

            sage: from sage.modular.btquotients.btquotient import Edge,Vertex
            sage: v1 = Vertex(7,0,Matrix(ZZ,2,2,[1,2,3,18]))
            sage: v2 = Vertex(7,0,Matrix(ZZ,2,2,[3,2,1,18]))
            sage: e1 = Edge(7,0,Matrix(ZZ,2,2,[1,2,3,18]),v1,v2)
            sage: e1 != e1
            False
        """
        return not self.__eq__(other)


class BruhatTitsQuotient(SageObject, UniqueRepresentation):
    r"""
    This function computes the quotient of the Bruhat-Tits tree
    by an arithmetic quaternionic group. The group in question is the
    group of norm 1 elements in an Eichler `\ZZ[1/p]`-order of some (tame)
    level inside of a definite quaternion algebra that is unramified
    at the prime `p`. Note that this routine relies in Magma in the case
    `p = 2` or when `N^{+} > 1`.

    INPUT:

    - ``p`` -- a prime number

    - ``Nminus`` -- squarefree integer divisible by an odd number of
      distinct primes and relatively prime to p. This is the
      discriminant of the definite quaternion algebra that one is
      quotienting by.

    - ``Nplus`` -- integer coprime to pNminus (default: 1). This is
      the tame level. It need not be squarefree! If Nplus is not 1
      then the user currently needs magma installed due to sage's
      inability to compute well with nonmaximal Eichler orders in
      rational (definite) quaternion algebras.

    - ``character`` -- a Dirichlet character (default: ``None``) of modulus
      `pN^-N^+`

    - ``use_magma`` -- boolean (default: ``False``); if True, uses Magma
      for quaternion arithmetic

    - ``magma_session`` -- (default: ``None``) if specified, the Magma session
      to use

    EXAMPLES:

    Here is an example without a Dirichlet character::

        sage: X = BruhatTitsQuotient(13, 19)
        sage: X.genus()
        19
        sage: G = X.get_graph(); G
        Multi-graph on 4 vertices

    And an example with a Dirichlet character::

        sage: f = DirichletGroup(6)[1]
        sage: X = BruhatTitsQuotient(3,2*5*7,character = f)
        sage: X.genus()
        5

    .. NOTE::

        A sage implementation of Eichler orders in rational quaternions
        algebras would remove the dependency on magma.

    AUTHORS:

    - Marc Masdeu (2012-02-20)
    """
    @staticmethod
    def __classcall__(cls, p, Nminus, Nplus=1, character=None,
                      use_magma=False, seed=None, magma_session=None):
        """
        Ensure that a canonical BruhatTitsQuotient is created.

        EXAMPLES::

            sage: BruhatTitsQuotient(3,17) is BruhatTitsQuotient(3,17,1)
            True
        """
        return super().__classcall__(cls, p, Nminus, Nplus,
                                     character, use_magma,
                                     seed, magma_session)

    def __init__(self, p, Nminus, Nplus=1, character=None,
                 use_magma=False, seed=None, magma_session=None):
        """
        Compute the quotient of the Bruhat-Tits tree by an arithmetic
        quaternionic group.

        EXAMPLES::

            sage: Y = BruhatTitsQuotient(19,11)
            sage: TestSuite(Y).run()
        """
        Nminus = Integer(Nminus)
        Nplus = Integer(Nplus)
        p = Integer(p)
        lev = p * Nminus
        self._order_is_initialized = False
        if character is not None:
            extra_level = character.conductor()
            if not extra_level.is_squarefree():
                raise ValueError("character must be of squarefree conductor")
            self._trivial_character = False
        else:
            G = DirichletGroup(lev * Nplus)
            character = G([1] * G.ngens())
            extra_level = 1
            self._trivial_character = True

        if not p.is_prime():
            raise ValueError("p must be a prime")
        if not lev.is_squarefree():
            raise ValueError("level must be squarefree")
        if (gcd(lev, Nplus) > 1):
            raise ValueError("level and conductor must be coprime")

        # if len(Nminus.factor()) % 2 != 1:
        #     raise ValueError("Nminus should be divisible by an odd number of primes")

        self._pN = p
        self._p = p
        self._Nminus = Nminus
        self._Nplus = Nplus
        if use_magma or self._Nplus != 1 or self._p == 2:
            try:
                if magma_session is None:
                    from sage.interfaces.magma import magma
                    self._magma = magma
                else:
                    self._magma = magma_session
                self._magma(p)
            except RuntimeError:
                raise NotImplementedError('Sage does not know yet how to work with the kind of orders that you are trying to use. Try installing Magma first and set it up so that Sage can use it.')

            # This is added for debugging, in order to have reproducible results
            if seed is not None:
                self._magma.function_call('SetSeed', seed, nvals=0)
            self._use_magma = True
        else:
            self._use_magma = False

        self._BT = BruhatTitsTree(p)

        self._prec = -1

        self._cached_vertices = {}
        self._cached_edges = {}
        self._cached_paths = {}
        self._cached_decomps = {}
        self._cached_equivalent = {}
        self._CM_points = {}

        self._V = (QQ ** 4).ambient_module().change_ring(ZZ)
        self._Mat_44 = MatrixSpace(ZZ, 4, 4)
        self._Mat_22 = MatrixSpace(ZZ, 2, 2)
        self._Mat_41 = MatrixSpace(ZZ, 4, 1)
        if extra_level == 1:
            self._extra_level = []
        else:
            self._extra_level = [ff[0] for ff in extra_level.factor()]
            self.get_extra_embedding_matrices()
        self._character = character
        self._Xv = [self._Mat_22([1, 0, 0, 0]),
                    self._Mat_22([0, 1, 0, 0]),
                    self._Mat_22([0, 0, 1, 0]),
                    self._Mat_22([0, 0, 0, 1])]
        self._Xe = [self._Mat_22([1, 0, 0, 0]),
                    self._Mat_22([0, 1, 0, 0]),
                    self._Mat_22([0, 0, self._p, 0]),
                    self._Mat_22([0, 0, 0, 1])]

    def _cache_key(self):
        r"""
        Return a hash of ``self``, for using in caching.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,13)
            sage: X._cache_key() == BruhatTitsQuotient(5,13)._cache_key()
            True
            sage: X._cache_key() == BruhatTitsQuotient(5,11)._cache_key()
            False

            sage: Y = BruhatTitsQuotient(5,13,use_magma = True) # optional - magma
            sage: Y._cache_key() == X._cache_key() # optional - magma
            False
        """
        return hash((self._p, self._Nminus, self._Nplus, self._character, self._use_magma))

    __hash__ = _cache_key

    def _repr_(self):
        r"""
        Return the representation of ``self`` as a string.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,13); X
            Quotient of the Bruhat Tits tree of GL_2(QQ_5) with discriminant 13 and level 1
        """
        return "Quotient of the Bruhat Tits tree of GL_2(QQ_%s) with discriminant %s and level %s" % (self.prime(), self.Nminus().factor(), self.Nplus().factor())

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,13)
            sage: Y = BruhatTitsQuotient(p = 5, Nminus = 13, Nplus=1,seed = 1231)
            sage: X == Y
            True
        """
        if self._p != other._p:
            return False
        if self._Nminus != other._Nminus:
            return False
        if self._Nplus != other._Nplus:
            return False
        if self._character != other._character:
            return False
        return True

    def __ne__(self, other):
        r"""
        Compare ``self`` with ``other``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,13)
            sage: Y = BruhatTitsQuotient(p = 5, Nminus = 13, Nplus=1,seed = 1231)
            sage: X != Y
            False
        """
        return not self.__eq__(other)

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,13); latex(X)
            X(5 \cdot 13,1)\otimes_{\Bold{Z}} \Bold{F}_{5}
        """
        return "X(%s,%s)\\otimes_{\\Bold{Z}} \\Bold{F}_{%s}" % (latex(self.level().factor()), latex(self.Nplus().factor()), latex(self.prime()))

    def get_vertex_dict(self):
        r"""
        This function returns the vertices of the quotient viewed as
        a dict.

        OUTPUT: a Python dict with the vertices of the quotient

        EXAMPLES::

            sage: X = BruhatTitsQuotient(37,3)
            sage: X.get_vertex_dict()
            {[1 0]
            [0 1]: Vertex of Bruhat-Tits tree for p = 37, [ 1  0]
            [ 0 37]: Vertex of Bruhat-Tits tree for p = 37}
        """
        try:
            return self._boundary
        except AttributeError:
            self._compute_quotient()
            return self._boundary

    def get_vertex_list(self):
        r"""
        Return a list of the vertices of the quotient.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(37,3)
            sage: X.get_vertex_list()
            [Vertex of Bruhat-Tits tree for p = 37, Vertex of Bruhat-Tits tree for p = 37]
        """
        try:
            return self._vertex_list
        except AttributeError:
            self._compute_quotient()
            return self._vertex_list

    def get_edge_list(self):
        r"""
        Return a list of ``Edge`` which represent a fundamental
        domain inside the Bruhat-Tits tree for the quotient.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(37,3)
            sage: len(X.get_edge_list())
            8
        """
        try:
            return self._edge_list
        except AttributeError:
            self._compute_quotient()
            return self._edge_list

    def get_list(self):
        r"""
        Return a list of ``Edge`` which represent a fundamental
        domain inside the Bruhat-Tits tree for the quotient,
        together with a list of the opposite edges. This is used
        to work with automorphic forms.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(37,3)
            sage: len(X.get_list())
            16
        """
        E = self.get_edge_list()
        return E + [e.opposite for e in E]

    def get_nontorsion_generators(self):
        r"""
        Use a fundamental domain in the Bruhat-Tits tree, and
        certain gluing data for boundary vertices, in order to compute
        a collection of generators for the nontorsion part
        of the arithmetic quaternionic group that one is quotienting by.
        This is analogous to using a polygonal rep. of a compact real
        surface to present its fundamental domain.

        OUTPUT:

        - A generating list of elements of an arithmetic
          quaternionic group.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,13)
            sage: len(X.get_nontorsion_generators())
            3
        """
        try:
            return list(self._nontorsion_generators)
        except AttributeError:
            self._compute_quotient()
            return list(self._nontorsion_generators)

    @cached_method
    def get_generators(self):
        r"""
        Use a fundamental domain in the Bruhat-Tits tree, and
        certain gluing data for boundary vertices, in order to compute
        a collection of generators for the arithmetic quaternionic
        group that one is quotienting by. This is analogous to using a
        polygonal rep. of a compact real surface to present its
        fundamental domain.

        OUTPUT:

        - A generating list of elements of an arithmetic
          quaternionic group.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,2)
            sage: len(X.get_generators())
            2
        """
        ans = self.get_nontorsion_generators()
        for s in self.get_vertex_stabs():
            for o in s:
                if o[2]:
                    ans.append(o[0])
                    break
        return ans

    def _compute_invariants(self):
        """
        Compute certain invariants from the level data of the quotient
        which allow one to compute the genus of the curve.

        Details to be found in Theorem 9 of [FM2014]_.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(23,11)
            sage: X._compute_invariants()
        """
        Nplus = self._Nplus
        lev = self._Nminus
        e4 = 1
        e3 = 1
        mu = Nplus
        for f in lev.factor():
            e4 *= (1 - kronecker_symbol(-4, Integer(f[0])))
            e3 *= (1 - kronecker_symbol(-3, Integer(f[0])))
            mu *= Integer(f[0]) - 1
        for f in Nplus.factor():
            if (f[1] == 1):
                e4 *= (1 + kronecker_symbol(-4, Integer(f[0])))
                e3 *= (1 + kronecker_symbol(-3, Integer(f[0])))
            else:
                if kronecker_symbol(-4, Integer(f[0])) == 1:
                    e4 *= 2
                else:
                    e4 = 0
                if kronecker_symbol(-3, Integer(f[0])) == 1:
                    e3 *= 2
                else:
                    e3 = 0
            mu *= 1 + 1 / Integer(f[0])
        self.e3 = e3
        self.e4 = e4
        self.mu = mu

    @lazy_attribute
    def e3(self):
        r"""
        Compute the `e_3` invariant defined by the formula

        .. MATH::

           e_k =\prod_{\ell\mid pN^-}\left(1-\left(\frac{-3}{\ell}\right)\right)\prod_{\ell \| N^+}\left(1+\left(\frac{-3}{\ell}\right)\right)\prod_{\ell^2\mid N^+} \nu_\ell(3)

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(31,3)
            sage: X.e3
            1
        """
        self._compute_invariants()
        return self.e3

    @lazy_attribute
    def e4(self):
        r"""
        Compute the `e_4` invariant defined by the formula

        .. MATH::

           e_k =\prod_{\ell\mid pN^-}\left(1-\left(\frac{-k}{\ell}\right)\right)\prod_{\ell \| N^+}\left(1+\left(\frac{-k}{\ell}\right)\right)\prod_{\ell^2\mid N^+} \nu_\ell(k)

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(31,3)
            sage: X.e4
            2
        """
        self._compute_invariants()
        return self.e4

    @lazy_attribute
    def mu(self):
        """
        Compute the mu invariant of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(29,3)
            sage: X.mu
            2
        """
        self._compute_invariants()
        return self.mu

    @cached_method
    def get_num_verts(self):
        r"""
        Return the number of vertices in the quotient using the formula
        `V = 2(\mu/12 + e_3/3 + e_4/4)`.

        OUTPUT:

        - An integer (the number of vertices)

        EXAMPLES::

            sage: X = BruhatTitsQuotient(29,11)
            sage: X.get_num_verts()
            4
        """
        return 2 * Integer(self.mu / 12 + self.e3 / 3 + self.e4 / 4)

    @cached_method
    def get_num_ordered_edges(self):
        """
        Return the number of ordered edges `E` in the quotient using
        the formula relating the genus `g` with the number of vertices `V`
        and that of unordered edges `E/2`: `E = 2(g + V - 1)`.

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,2)
            sage: X.get_num_ordered_edges()
            2
        """
        return 2 * (self.genus() + self.get_num_verts() - 1)

    def genus_no_formula(self):
        """
        Compute the genus of the quotient from the data of the
        quotient graph. This should agree with self.genus().

        OUTPUT: integer

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,2*3*29)
            sage: X.genus_no_formula()
            17
            sage: X.genus_no_formula() == X.genus()
            True
        """
        return ZZ(1 - len(self.get_vertex_list()) + len(self.get_edge_list()))

    @cached_method
    def genus(self):
        r"""
        Compute the genus of the quotient graph using a formula
        This should agree with self.genus_no_formula().

        Compute the genus of the Shimura curve
        corresponding to this quotient via Cerednik-Drinfeld. It is
        computed via a formula and not in terms of the quotient graph.

        INPUT:

        - ``level`` -- integer (default: ``None``); a level. By default, use that
          of ``self``.

        - ``Nplus`` -- integer (default: ``None``); a conductor. By default, use
          that of ``self``.

        OUTPUT: integer equal to the genus

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,2*5*31)
            sage: X.genus()
            21
            sage: X.genus() == X.genus_no_formula()
            True
        """
        return self.dimension_harmonic_cocycles(2)

    @cached_method
    def dimension_harmonic_cocycles(self, k, lev=None, Nplus=None,
                                    character=None):
        r"""
        Compute the dimension of the space of harmonic cocycles
        of weight `k` on ``self``.

        OUTPUT: integer equal to the dimension

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: [X.dimension_harmonic_cocycles(k) for k in range(2,20,2)]
            [1, 4, 4, 8, 8, 12, 12, 16, 16]

            sage: X = BruhatTitsQuotient(2,5) # optional - magma
            sage: [X.dimension_harmonic_cocycles(k) for k in range(2,40,2)] # optional - magma
            [0, 1, 3, 1, 3, 5, 3, 5, 7, 5, 7, 9, 7, 9, 11, 9, 11, 13, 11]

            sage: X = BruhatTitsQuotient(7, 2 * 3 * 5)
            sage: X.dimension_harmonic_cocycles(4)
            12
            sage: X = BruhatTitsQuotient(7, 2 * 3 * 5 * 11 * 13)
            sage: X.dimension_harmonic_cocycles(2)
            481
            sage: X.dimension_harmonic_cocycles(4)
            1440
        """
        k = ZZ(k)
        if lev is None:
            lev = self._p * self._Nminus
        else:
            lev = ZZ(lev)
        if Nplus is None:
            Nplus = self._Nplus
        else:
            Nplus = ZZ(Nplus)

        if character is None:
            if not self._trivial_character:
                character = self._character
                lN = lev * Nplus
                kernel = [r for r in lN.coprime_integers(lN)
                          if character(r) == 1]
            else:
                character = None
                kernel = None

        if k == 0:
            return 0

        verbose('Computing dimension for (k,level,nplus,char) = (%s, %s, %s, %s)' % (k, lev, Nplus, character), level=2)

        if lev == 1:
            return Gamma0(Nplus).dimension_cusp_forms(k=k)

        f = lev.factor()
        if any(l[1] != 1 for l in f):
            raise NotImplementedError('The level should be squarefree for '
                                      'this function to work... Sorry!')

        def GH(N, ker):
            return Gamma0(N) if character is None else GammaH_constructor(N, ker)

        def mumu(N):
            p = 1
            for _, r in ZZ(N).factor():
                if r > 2:
                    return ZZ.zero()
                elif r == 1:
                    p *= -2
            return ZZ(p)
        return sum([mumu(lev // d) * GH(d * Nplus, kernel).dimension_cusp_forms(k) for d in lev.divisors()])

    def Nplus(self):
        r"""
        Return the tame level `N^+`.

        OUTPUT: integer equal to `N^+`

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7,1)
            sage: X.Nplus()
            1
        """
        return self._Nplus

    def Nminus(self):
        r"""
        Return the discriminant of the relevant definite
        quaternion algebra.

        OUTPUT:

        An integer equal to `N^-`.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: X.Nminus()
            7
        """
        return self._Nminus

    @cached_method
    def level(self):
        r"""
        Return `p N^-`, which is the discriminant of the
        indefinite quaternion algebra that is uniformed by
        Cerednik-Drinfeld.

        OUTPUT:

        An integer equal to `p N^-`.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: X.level()
            35
        """
        return self._Nminus * self._p

    def prime(self):
        r"""
        Return the prime one is working with.

        OUTPUT: integer equal to the fixed prime `p`

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: X.prime()
            5
        """
        return self._p

    def get_graph(self):
        r"""
        Return the quotient graph (and compute it if needed).

        OUTPUT: a graph representing the quotient of the Bruhat-Tits tree

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,5)
            sage: X.get_graph()
            Multi-graph on 2 vertices
        """
        try:
            return self._S
        except AttributeError:
            self._compute_quotient()
            return self._S

    def get_fundom_graph(self):
        r"""
        Return the fundamental domain (and computes it if needed).

        OUTPUT: a fundamental domain for the action of `\Gamma`

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,5)
            sage: X.get_fundom_graph()
            Graph on 24 vertices
        """
        try:
            return self._Sfun
        except AttributeError:
            self._compute_quotient()
            return self._Sfun

    def plot(self, *args, **kwargs):
        r"""
        Plot the quotient graph.

        OUTPUT: a plot of the quotient graph

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,23)
            sage: X.plot()                                                              # needs sage.plot
            Graphics object consisting of 17 graphics primitives
        """
        S = self.get_graph()
        vertex_colors = {}
        v0 = Matrix(ZZ, 2, 2, [1, 0, 0, 1])
        v0.set_immutable()
        rainbow_color = rainbow(len(self.get_vertex_list()))
        for v in S.vertex_iterator():
            key = rainbow_color[S.get_vertex(v).label]
            if key in vertex_colors:
                vertex_colors[key].append(v)
            else:
                vertex_colors[key] = [v]

        my_args = {}
        my_args['vertex_colors'] = vertex_colors
        my_args['color_by_label'] = True
        my_args['vertex_labels'] = False
        my_args.update(kwargs)
        return S.plot(*args, **my_args)

    def plot_fundom(self, *args, **kwargs):
        r"""
        Plot a fundamental domain.

        OUTPUT: a plot of the fundamental domain

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,23)
            sage: X.plot_fundom()                                                       # needs sage.plot
            Graphics object consisting of 88 graphics primitives
        """
        S = self.get_fundom_graph()
        vertex_colors = {}
        rainbow_color = rainbow(len(self.get_vertex_list()))
        for v in S.vertex_iterator():
            key = rainbow_color[S.get_vertex(v).label]
            if key in vertex_colors:
                vertex_colors[key].append(v)
            else:
                vertex_colors[key] = [v]

        my_args = {}
        my_args['vertex_colors'] = vertex_colors
        my_args['color_by_label'] = True
        my_args['vertex_labels'] = True
        my_args.update(kwargs)
        return S.plot(*args, **my_args)

    def is_admissible(self, D):
        r"""
        Test whether the imaginary quadratic field of
        discriminant `D` embeds in the quaternion algebra. It
        furthermore tests the Heegner hypothesis in this setting
        (e.g., is `p` inert in the field, etc).

        INPUT:

        - ``D`` -- integer whose squarefree part will define the
          quadratic field

        OUTPUT: boolean describing whether the quadratic field is admissible

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: [X.is_admissible(D) for D in range(-1,-20,-1)]
            [False, True, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, True, False]
        """
        disc = fundamental_discriminant(D)
        for f in self.level().factor():
            if kronecker_symbol(disc, f[0]) != -1:
                return False
        for f in self._Nplus.factor():
            if kronecker_symbol(disc, f[0]) != 1:
                return False
        return True

    def _local_splitting_map(self, prec):
        r"""
        Return an embedding of the definite quaternion algebra
        into the algebra of 2x2 matrices with coefficients in `\QQ_p`.

        INPUT:

        - ``prec`` -- integer; the precision of the splitting

        OUTPUT: a function giving the splitting

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,3)
            sage: phi = X._local_splitting_map(10)
            sage: B.<i,j,k> = QuaternionAlgebra(3)
            sage: phi(i)**2 == QQ(i**2)*phi(B(1))
            True
        """
        I, J, K = self._local_splitting(prec)

        def phi(q):
            R = I.parent()
            v = q.coefficient_tuple()
            return R(v[0] + I * v[1] + J * v[2] + K * v[3])
        return phi

    def _local_splitting(self, prec):
        r"""
        Find an embedding of the definite quaternion algebra
        into the algebra of 2x2 matrices with coefficients in `\QQ_p`.

        INPUT:

        - ``prec`` -- integer; the precision of the splitting

        OUTPUT: matrices `I`, `J`, `K` giving the splitting

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,3)
            sage: phi = X._local_splitting_map(10)
            sage: B.<i,j,k> = QuaternionAlgebra(3)
            sage: phi(i)**2 == QQ(i**2)*phi(B(1))
            True
        """
        assert not self._use_magma
        if prec <= self._prec:
            return self._II, self._JJ, self._KK

        A = self.get_quaternion_algebra()

        ZZp = Zp(self._p, prec)
        v = A.invariants()
        a = ZZp(v[0])
        b = ZZp(v[1])
        if (A.base_ring() != QQ):
            raise ValueError("must be rational quaternion algebra")
        if (A.discriminant() % self._p == 0):
            raise ValueError("p (=%s) must be an unramified prime" % self._p)
        M = MatrixSpace(ZZp, 2)

        if a.is_square():
            alpha = a.sqrt()
            self._II = M([alpha, 0, 2 * alpha, -alpha])
            self._JJ = M([b, -b, b - 1, -b])
        else:
            self._II = M([0, a, 1, 0])
            z = 0
            self._JJ = 0
            while self._JJ == 0:
                c = a * z * z + b
                if c.is_square():
                    x = c.sqrt()
                    self._JJ = M([x, -a * z, z, -x])
                else:
                    z += 1
        self._KK = self._II * self._JJ
        return self._II, self._JJ, self._KK

    def _compute_embedding_matrix(self, prec, force_computation=False):
        r"""
        Return a matrix representing the embedding with the
        given precision.

        INPUT:

        - ``prec`` -- integer; the precision of the embedding matrix

        EXAMPLES:

        Note that the entries of the matrix are elements of Zmod::

            sage: X = BruhatTitsQuotient(3,7)
            sage: A = X.get_embedding_matrix(10)  # indirect doctest
            sage: R = A.base_ring()
            sage: B = X.get_eichler_order_basis()
            sage: R(B[0].reduced_trace()) == A[0,0]+A[3,0]
            True
        """
        if self._use_magma:
            if not force_computation:
                try:
                    return Matrix(Zmod(self._pN), 4, 4,
                                  self._cached_Iota0_matrix)
                except AttributeError:
                    pass

            Ord = self.get_eichler_order(magma=True)  # force_computation = force_computation)
            OrdMax = self.get_maximal_order(magma=True)

            OBasis = Ord.Basis()
            verbose(f'Calling magma: pMatrixRing, args = [{OrdMax}, {self._p}]')
            M, f, rho = self._magma.function_call('pMatrixRing', args=[OrdMax, self._p], params={'Precision': 2000}, nvals=3)
            v = [f.Image(OBasis[i]) for i in [1, 2, 3, 4]]

            self._cached_Iota0_matrix = [v[kk][ii, jj].sage()
                                         for ii in range(1, 3)
                                         for jj in range(1, 3)
                                         for kk in range(4)]
            return Matrix(Zmod(self._pN), 4, 4, self._cached_Iota0_matrix)
        else:
            phi = self._local_splitting_map(prec)
            B = self.get_eichler_order_basis()
            return Matrix(Zmod(self._p ** prec), 4, 4,
                          [phi(B[kk])[ii, jj] for ii in range(2)
                           for jj in range(2) for kk in range(4)])

    @cached_method
    def get_extra_embedding_matrices(self):
        r"""
        Return a list of  matrices representing the different embeddings.

        .. NOTE::

           The precision is very low (currently set to 5 digits),
           since these embeddings are only used to apply a character.

        EXAMPLES:

        This portion of the code is only relevant when working with a
        nontrivial Dirichlet character. If there is no such character
        then the code returns an empty list. Even if the character is
        not trivial it might return an empty list::

            sage: f = DirichletGroup(6)[1]
            sage: X = BruhatTitsQuotient(3,2*5*7,character = f)
            sage: X.get_extra_embedding_matrices()
            []

        ::

            sage: f = DirichletGroup(6)[1]
            sage: X = BruhatTitsQuotient(5,2,3, character = f, use_magma=True) # optional - magma
            sage: X.get_extra_embedding_matrices() # optional - magma
            [
            [1 0 2 0]
            [0 0 2 0]
            [0 0 0 0]
            [1 0 2 2]
            ]
        """
        if not self._use_magma or len(self._extra_level) == 0:
            return []
        n_iters = 0
        Ord = self.get_eichler_order(magma=True)
        OrdMax = self.get_maximal_order(magma=True)
        OBasis = Ord.Basis()
        extra_embeddings = []
        success = False
        while not success:
            success = True
            for l in self._extra_level:
                success = False
                found = False
                while not found:
                    verbose(f'Calling magma: pMatrixRing, args = [{OrdMax}, {l}]')
                    M, f, rho = self._magma.function_call('pMatrixRing', args=[OrdMax, l], params={'Precision': 20}, nvals=3)
                    v = [f.Image(OBasis[i]) for i in [1, 2, 3, 4]]
                    if all(Qp(l, 5)(v[kk][2, 1].sage()).valuation() >= 1 for kk in range(4)) and not all(Qp(l, 5)(v[kk][2, 1].sage()).valuation() >= 2 for kk in range(4)):
                        found = True
                        success = True
                    else:
                        n_iters += 1
                        verbose('Restarting magma...')
                        self._magma.quit()
                        from sage.interfaces.magma import magma
                        self._magma = magma
                        self._magma.function_call('SetSeed', n_iters, nvals=0)
                        self._order_is_initialized = False
                        self._init_order()
                        self._compute_embedding_matrix(self._prec,
                                                       force_computation=True)
                        Ord = self.get_eichler_order(magma=True)
                        OrdMax = self.get_maximal_order(magma=True)
                        OBasis = Ord.Basis()
                        extra_embeddings = []
                        success = False
                        break
                if not success:
                    break
                mat = Matrix(GF(l), 4, 4, [v[kk][ii, jj].sage()
                                           for ii in range(1, 3)
                                           for jj in range(1, 3)
                                           for kk in range(4)])
                extra_embeddings.append(mat)
        return extra_embeddings

    def _increase_precision(self, amount=1):
        r"""
        Increase the working precision.

        INPUT:

        - ``amount`` Integer (default: 1). The amount by which to
          increase the precision.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,101)
            sage: X.get_embedding_matrix()
            [    O(3) 1 + O(3) 1 + O(3) 1 + O(3)]
            [2 + O(3)     O(3) 2 + O(3) 2 + O(3)]
            [1 + O(3) 1 + O(3)     O(3) 2 + O(3)]
            [1 + O(3) 2 + O(3) 2 + O(3) 2 + O(3)]
            sage: X._increase_precision(5)
            sage: X.get_embedding_matrix()[0,0]
            2*3^3 + 2*3^5 + O(3^6)
        """
        if amount >= 1:
            self.get_embedding_matrix(prec=self._prec + amount)

    def get_embedding_matrix(self, prec=None, exact=False):
        r"""
        Return the matrix of the embedding.

        INPUT:

        - ``exact`` -- boolean (default: ``False``); if ``True``, return an
          embedding into a matrix algebra with coefficients in a
          number field. Otherwise, embed into matrices over `p`-adic
          numbers.

        - ``prec`` -- integer (default: ``None``); if specified, return the
          matrix with precision ``prec``. Otherwise, return the
          cached matrix (with the current working precision).

        OUTPUT: a 4x4 matrix representing the embedding

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,2*3*5)
            sage: X.get_embedding_matrix(4)
            [                      1 + O(7^4)         5 + 2*7 + 3*7^3 + O(7^4) 4 + 5*7 + 6*7^2 + 6*7^3 + O(7^4)       6 + 3*7^2 + 4*7^3 + O(7^4)]
            [                          O(7^4)                           O(7^4)                   3 + 7 + O(7^4) 1 + 6*7 + 3*7^2 + 2*7^3 + O(7^4)]
            [                          O(7^4)         2 + 5*7 + 6*7^3 + O(7^4) 3 + 5*7 + 6*7^2 + 6*7^3 + O(7^4)         3 + 3*7 + 3*7^2 + O(7^4)]
            [                      1 + O(7^4) 3 + 4*7 + 6*7^2 + 3*7^3 + O(7^4)                   3 + 7 + O(7^4) 1 + 6*7 + 3*7^2 + 2*7^3 + O(7^4)]
            sage: X.get_embedding_matrix(3)
            [                      1 + O(7^4)         5 + 2*7 + 3*7^3 + O(7^4) 4 + 5*7 + 6*7^2 + 6*7^3 + O(7^4)       6 + 3*7^2 + 4*7^3 + O(7^4)]
            [                          O(7^4)                           O(7^4)                   3 + 7 + O(7^4) 1 + 6*7 + 3*7^2 + 2*7^3 + O(7^4)]
            [                          O(7^4)         2 + 5*7 + 6*7^3 + O(7^4) 3 + 5*7 + 6*7^2 + 6*7^3 + O(7^4)         3 + 3*7 + 3*7^2 + O(7^4)]
            [                      1 + O(7^4) 3 + 4*7 + 6*7^2 + 3*7^3 + O(7^4)                   3 + 7 + O(7^4) 1 + 6*7 + 3*7^2 + 2*7^3 + O(7^4)]
            sage: X.get_embedding_matrix(5)
            [                              1 + O(7^5)         5 + 2*7 + 3*7^3 + 6*7^4 + O(7^5) 4 + 5*7 + 6*7^2 + 6*7^3 + 6*7^4 + O(7^5)       6 + 3*7^2 + 4*7^3 + 5*7^4 + O(7^5)]
            [                                  O(7^5)                                   O(7^5)                           3 + 7 + O(7^5)   1 + 6*7 + 3*7^2 + 2*7^3 + 7^4 + O(7^5)]
            [                                  O(7^5)         2 + 5*7 + 6*7^3 + 5*7^4 + O(7^5) 3 + 5*7 + 6*7^2 + 6*7^3 + 6*7^4 + O(7^5)         3 + 3*7 + 3*7^2 + 5*7^4 + O(7^5)]
            [                              1 + O(7^5)         3 + 4*7 + 6*7^2 + 3*7^3 + O(7^5)                           3 + 7 + O(7^5)   1 + 6*7 + 3*7^2 + 2*7^3 + 7^4 + O(7^5)]
        """
        if exact is True:
            try:
                return self._Iota_exact
            except AttributeError:
                raise RuntimeError('Exact splitting not available.')
        else:
            if prec is None:
                prec = self._prec

            if prec < 0:
                prec = 1

            if prec == self._prec:
                try:
                    return self._Iota
                except AttributeError:
                    pass

            self._pN = self._p ** prec
            self._R = Qp(self._p, prec=prec)

            if prec > self._prec:
                verbose('self._prec = %s, prec = %s' % (self._prec, prec))
                Iotamod = self._compute_embedding_matrix(prec)
                self._Iotainv_lift = Iotamod.inverse().lift()
                self._Iota = Matrix(self._R, 4, 4, [Iotamod[ii, jj]
                                                    for ii in range(4)
                                                    for jj in range(4)])

            self._prec = prec
            self._Iotainv = self._Mat_44([self._Iotainv_lift[ii, jj] % self._pN for ii in range(4) for jj in range(4)])
            return self._Iota

    def embed_quaternion(self, g, exact=False, prec=None):
        r"""
        Embed the quaternion element ``g`` into a matrix algebra.

        INPUT:

        - ``g`` -- a row vector of size `4` whose entries represent a
          quaternion in our basis

        - ``exact`` -- boolean (default: ``False``); if True, tries to embed
          ``g`` into a matrix algebra over a number field. If ``False``,
          the target is the matrix algebra over `\QQ_p`.

        OUTPUT:

        A 2x2 matrix with coefficients in `\QQ_p` if ``exact`` is
        False, or a number field if ``exact`` is True.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,2)
            sage: l = X.get_units_of_order()
            sage: len(l)
            12
            sage: l[3] # random
            [-1]
            [ 0]
            [ 1]
            [ 1]
            sage: u = X.embed_quaternion(l[3]); u # random
            [    O(7) 3 + O(7)]
            [2 + O(7) 6 + O(7)]
            sage: X._increase_precision(5)
            sage: v = X.embed_quaternion(l[3]); v # random
            [                7 + 3*7^2 + 7^3 + 4*7^4 + O(7^6)             3 + 7 + 3*7^2 + 7^3 + 4*7^4 + O(7^6)]
            [            2 + 7 + 3*7^2 + 7^3 + 4*7^4 + O(7^6) 6 + 5*7 + 3*7^2 + 5*7^3 + 2*7^4 + 6*7^5 + O(7^6)]
            sage: u == v
            True
        """
        if exact:
            return Matrix(self.get_splitting_field(), 2, 2,
                          (self.get_embedding_matrix(exact=True) * g).list())
        else:
            A = self.get_embedding_matrix(prec=prec) * g
            return Matrix(self._R, 2, 2, A.list())

    embed = embed_quaternion

    def get_embedding(self, prec=None):
        r"""
        Return a function which embeds quaternions into a matrix
        algebra.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,3)
            sage: f = X.get_embedding(prec = 4)
            sage: b = Matrix(ZZ,4,1,[1,2,3,4])
            sage: f(b)
            [2 + 3*5 + 2*5^2 + 4*5^3 + O(5^4)       3 + 2*5^2 + 4*5^3 + O(5^4)]
            [        5 + 5^2 + 3*5^3 + O(5^4)           4 + 5 + 2*5^2 + O(5^4)]
        """
        A = self.get_embedding_matrix(prec=prec)
        return lambda g: Matrix(self._R, 2, 2, (A * g).list())

    def get_edge_stabilizers(self):
        r"""
        Compute the stabilizers in the arithmetic group of all
        edges in the Bruhat-Tits tree within a fundamental domain for
        the quotient graph. The stabilizers of an edge and its
        opposite are equal, and so we only store half the data.

        OUTPUT:

        A list of lists encoding edge stabilizers. It contains one
        entry for each edge. Each entry is a list of data
        corresponding to the group elements in the stabilizer of the
        edge. The data consists of: (0) a column matrix representing
        a quaternion, (1) the power of `p` that one needs to divide
        by in order to obtain a quaternion of norm 1, and hence an
        element of the arithmetic group `\Gamma`, (2) a boolean that
        is only used to compute spaces of modular forms.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,2)
            sage: s = X.get_edge_stabilizers()
            sage: len(s) == X.get_num_ordered_edges()/2
            True
            sage: len(s[0])
            3
        """
        try:
            return self._edge_stabs
        except AttributeError:
            self._edge_stabs = [self._stabilizer(e.rep, as_edge=True)
                                for e in self.get_edge_list()]
            return self._edge_stabs

    def get_stabilizers(self):
        r"""
        Compute the stabilizers in the arithmetic group of all
        edges in the Bruhat-Tits tree within a fundamental domain for
        the quotient graph. This is similar to get_edge_stabilizers, except
        that here we also store the stabilizers of the opposites.

        OUTPUT:

        A list of lists encoding edge stabilizers. It contains one
        entry for each edge. Each entry is a list of data
        corresponding to the group elements in the stabilizer of the
        edge. The data consists of: (0) a column matrix representing
        a quaternion, (1) the power of `p` that one needs to divide
        by in order to obtain a quaternion of norm 1, and hence an
        element of the arithmetic group `\Gamma`, (2) a boolean that
        is only used to compute spaces of modular forms.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,5)
            sage: s = X.get_stabilizers()
            sage: len(s) == X.get_num_ordered_edges()
            True
            sage: gamma = X.embed_quaternion(s[1][0][0][0],prec = 20)
            sage: v = X.get_edge_list()[0].rep
            sage: X._BT.edge(gamma*v) == v
            True
        """
        S = self.get_edge_stabilizers()
        return S + S

    def get_vertex_stabs(self):
        r"""
        This function computes the stabilizers in the arithmetic
        group of all vertices in the Bruhat-Tits tree within a
        fundamental domain for the quotient graph.

        OUTPUT:

        A list of vertex stabilizers. Each vertex stabilizer is a
        finite cyclic subgroup, so we return generators for these
        subgroups.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(13,2)
            sage: S = X.get_vertex_stabs()
            sage: gamma = X.embed_quaternion(S[0][0][0],prec = 20)
            sage: v = X.get_vertex_list()[0].rep
            sage: X._BT.vertex(gamma*v) == v
            True
        """
        try:
            return self._vertex_stabs
        except AttributeError:
            self._vertex_stabs = [self._stabilizer(v.rep, as_edge=False)
                                  for v in self.get_vertex_list()]
            return self._vertex_stabs

    def get_quaternion_algebra(self):
        r"""
        Return the underlying quaternion algebra.

        OUTPUT: the underlying definite quaternion algebra

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: X.get_quaternion_algebra()
            Quaternion Algebra (-1, -7) with base ring Rational Field
        """
        try:
            return self._A
        except AttributeError:
            pass
        self._init_order()
        return self._A

    def get_eichler_order(self, magma=False, force_computation=False):
        r"""
        Return the underlying Eichler order of level `N^+`.

        OUTPUT: an Eichler order

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: X.get_eichler_order()
            Order of Quaternion Algebra (-1, -7) with base ring Rational Field with basis (1/2 + 1/2*j, 1/2*i + 1/2*k, j, k)
        """
        if magma:
            if not force_computation:
                try:
                    return self._Omagma
                except AttributeError:
                    pass
            self._init_order()
            return self._Omagma
        else:
            try:
                return self._O
            except AttributeError:
                pass
            self._init_order()
            return self._O

    def get_maximal_order(self, magma=False, force_computation=False):
        r"""
        Return the underlying maximal order containing the
        Eichler order.

        OUTPUT: a maximal order

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: X.get_maximal_order()
            Order of Quaternion Algebra (-1, -7) with base ring Rational Field with basis (1/2 + 1/2*j, 1/2*i + 1/2*k, j, k)
        """
        if magma:
            if not force_computation:
                try:
                    return self._OMaxmagma
                except AttributeError:
                    pass
            self._init_order()
            return self._OMaxmagma
        else:
            try:
                return self._OMax
            except AttributeError:
                pass
            self._init_order()
            return self._OMax

    def get_splitting_field(self):
        r"""
        Return a quadratic field that splits the quaternion
        algebra attached to ``self``. Currently requires Magma.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,11)
            sage: X.get_splitting_field()
            Traceback (most recent call last):
            ...
            NotImplementedError: Sage does not know yet how to work with the kind of orders that you are trying to use. Try installing Magma first and set it up so that Sage can use it.

        If we do have Magma installed, then it works::

            sage: X = BruhatTitsQuotient(5,11,use_magma=True) # optional - magma
            sage: X.get_splitting_field() # optional - magma
            Number Field in a with defining polynomial X1^2 + 11
        """
        if not self._use_magma:
            raise NotImplementedError('Sage does not know yet how to work with the kind of orders that you are trying to use. Try installing Magma first and set it up so that Sage can use it.')
        try:
            return self._FF
        except AttributeError:
            pass
        self._compute_exact_splitting()
        return self._FF

    def get_eichler_order_basis(self):
        r"""
        Return a basis for the global Eichler order.

        OUTPUT: basis for the underlying Eichler order of level Nplus

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,11)
            sage: X.get_eichler_order_basis()
            [1/2 + 1/2*j, 1/2*i + 1/2*k, j, k]
        """
        try:
            return self._B
        except AttributeError:
            pass
        self._init_order()
        return self._B

    def get_eichler_order_quadform(self):
        r"""
        This function return the norm form for the underlying
        Eichler order of level ``Nplus``. Required for finding elements in
        the arithmetic subgroup Gamma.

        OUTPUT: the norm form of the underlying Eichler order

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,11)
            sage: X.get_eichler_order_quadform()
            Quadratic form in 4 variables over Integer Ring with coefficients:
            [ 3 0 11 0 ]
            [ * 3 0 11 ]
            [ * * 11 0 ]
            [ * * * 11 ]
        """
        try:
            return self._OQuadForm
        except AttributeError:
            pass
        self._init_order()
        return self._OQuadForm

    def get_eichler_order_quadmatrix(self):
        r"""
        This function returns the matrix of the quadratic form of
        the underlying Eichler order in the fixed basis.

        OUTPUT: a 4x4 integral matrix describing the norm form

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,11)
            sage: X.get_eichler_order_quadmatrix()
            [ 6  0 11  0]
            [ 0  6  0 11]
            [11  0 22  0]
            [ 0 11  0 22]
        """
        try:
            return self._OM
        except AttributeError:
            pass
        self._init_order()
        return self._OM

    @cached_method
    def get_units_of_order(self):
        r"""
        Return the units of the underlying Eichler
        `\ZZ`-order. This is a finite group since the order lives in a
        definite quaternion algebra over `\QQ`.

        OUTPUT:

        A list of elements of the global Eichler `\ZZ`-order of
        level `N^+`.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,11)
            sage: X.get_units_of_order()
            [
            [ 0]  [-2]
            [-2]  [ 0]
            [ 0]  [ 1]
            [ 1], [ 0]
            ]
        """
        OM = self.get_eichler_order_quadmatrix()
        v = pari('qfminim(%s,2,0, flag = 2)' % (OM.__pari__()))
        n_units = Integer(v[0].sage() / 2)
        v = pari('qfminim(%s,2,%s, flag = 2)' % ((OM.__pari__()), n_units))
        O_units = []
        for jj in range(n_units):
            vec = Matrix(ZZ, 4, 1, [v[2][ii, jj].sage() for ii in range(4)])
            O_units.append(vec)
        return O_units

    @cached_method
    def _get_Up_data(self):
        r"""
        Return (compute if necessary) Up data.

        The Up data is a vector of length `p`, and each entry consists
        of the corresponding data for the matrix `[p,a,0,1]` where a
        varies from 0 to `p-1`. The data is a tuple (acter,edge_images),
        with edge images being of type ``DoubleCosetReduction``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: [o[0] for o in X._get_Up_data()]
            [
            [1/3   0]  [-1/3  1/3]  [-2/3  1/3]
            [  0   1], [   1    0], [   1    0]
            ]
        """
        E = self.get_edge_list()
        vec_a = self._BT.subdivide([1], 1)
        return [[alpha.inverse(),
                 [DoubleCosetReduction(self, e.rep * alpha) for e in E]
                 + [DoubleCosetReduction(self, e.opposite.rep * alpha)
                    for e in E]]
                for alpha in vec_a]

    @cached_method
    def _get_atkin_lehner_data(self, q):
        r"""
        Return (and compute if necessary) data to compute the
        Atkin-Lehner involution.

        INPUT:

        - ``q`` -- integer dividing p*Nminus*Nplus

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,5)
            sage: X._get_atkin_lehner_data(3)[0]
            [ 2]
            [ 4]
            [-3]
            [-2]
        """
        E = self.get_edge_list()

        nninc = -2
        V = []
        p = self._p
        while not V:
            nninc += 2
            V = [g for g in self._find_elements_in_order(q * self._p ** nninc)
                 if prod([self._character(ZZ((v * Matrix(ZZ, 4, 1, g))[0, 0]))
                          / self._character(p ** (nninc // 2))
                          for v in self.get_extra_embedding_matrices()]) == 1]

        beta1 = Matrix(QQ, 4, 1, V[0])

        success = False
        while not success:
            try:
                x = self.embed_quaternion(beta1)
                nn = x.determinant().valuation()
                T = [beta1,
                     [DoubleCosetReduction(self, x.adjugate() * e.rep,
                                           extrapow=nn) for e in E]]
                success = True
            except (PrecisionError, NotImplementedError):
                self._increase_precision(10)
        return T

    @cached_method
    def _get_hecke_data(self, l):
        r"""
        Return (and compute if necessary) data to compute the
        Hecke operator at a prime.

        INPUT:

        - ``l`` -- a prime l

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: len(X._get_hecke_data(5))
            2
        """
        E = self.get_edge_list()
        if (self.level() * self.Nplus()) % l == 0:
            Sset = []
        else:
            Sset = [self._p]
        BB = self._BB
        p = self._p
        T = []
        T0 = []
        V = []
        nninc = 0
        while not V:
            V = [g for g in self._find_elements_in_order(l * p ** nninc)
                 if prod([self._character(ZZ((v * Matrix(ZZ, 4, 1, g))[0, 0]))
                          / self._character((p ** (nninc // 2)))
                          for v in self.get_extra_embedding_matrices()]) == 1]
            if not V:
                nninc += 2

        alpha1 = V[0]
        alpha0 = self._conv(alpha1)

        alpha = Matrix(QQ, 4, 1, alpha1)
        alphamat = self.embed_quaternion(alpha)
        letters = self.get_nontorsion_generators()
        letters += [g for g in self._find_elements_in_order(1)
                    if prod([self._character(ZZ((v * Matrix(ZZ, 4, 1, g))[0, 0]))
                             / self._character((p ** (nninc // 2)))
                             for v in self.get_extra_embedding_matrices()]) == 1]
        n_iters = 0

        def enumerate_words(v, n=None):
            if n is None:
                n = []
            while True:
                add_new = True
                for jj in range(len(n)):
                    n[jj] += 1
                    if n[jj] != len(v):
                        add_new = False
                        break
                    else:
                        n[jj] = 0
                if add_new:
                    n.append(0)
                yield [v[x] for x in n]

        for wd in enumerate_words([self._conv(x) for x in letters]):
            if len(T) == l + 1:
                break
            v = prod(wd)
            n_iters += 1
            v0 = v * alpha0
            vinv = self.get_quaternion_algebra()(v0 ** (-1))
            new = True
            for tt in T0:
                r = vinv * tt
                r_in_order = BB * Matrix(QQ, 4, 1, r.coefficient_tuple())
                if all(a.is_S_integral(Sset) for a in r_in_order.list()):
                    new = False
                    break
            if new:
                v1 = BB * Matrix(QQ, 4, 1, v.coefficient_tuple())
                success = False
                while not success:
                    try:
                        x = self.embed_quaternion(v1, prec=max(self._prec, 40),
                                                  exact=False) * alphamat
                        nn = x.determinant().valuation()
                        dcr = [DoubleCosetReduction(self, x.adjugate() * e.rep,
                                                    extrapow=nn) for e in E]
                        T.append([v1, dcr])
                        success = True
                    except (PrecisionError, NotImplementedError):
                        self._increase_precision(10)
                        alphamat = self.embed_quaternion(alpha, prec=max(self._prec, 40), exact=False)
                T0.append(v0)
        return T, alpha

    def _find_equivalent_vertex(self, v0, V=None, valuation=None):
        r"""
        Find a vertex in ``V`` equivalent to ``v0``.

        INPUT:

        - ``v0`` -- a 2x2 matrix in `\ZZ_p` representing a
            vertex in the Bruhat-Tits tree

        - ``V`` -- list (default: ``None``); if a list of Vertex is given,
            restrict the search to the vertices in ``V``. Otherwise
            use all the vertices in a fundamental domain.

        - ``valuation`` -- integer (default: ``None``); the valuation
            of the determinant of ``v0``, if known (otherwise it is
            calculated)

        OUTPUT:

        A pair ``g``, ``v``, where ``v`` is a Vertex in ``V``
        equivalent to ``v0``, and ``g`` is such that `g\cdot v_0= v`.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: M = Matrix(ZZ,2,2,[1,3,2,7])
            sage: M.set_immutable()
            sage: X._find_equivalent_vertex(M)[-1] in X.get_vertex_list()
            True
        """
        try:
            return self._cached_vertices[v0]
        except KeyError:
            pass
        if V is None:
            V = self.get_vertex_list()
        if valuation is None:
            valuation = v0.determinant().valuation(self._p)
        parity = valuation % 2
        for v in V:
            if v.parity != parity:
                continue
            g = self._are_equivalent(v0, v.rep, False, valuation + v.valuation)
            if g is not None:
                self._cached_vertices[v0] = (g, v)
                return g, v
        return 0, None

    def _find_equivalent_edge(self, e0, E=None, valuation=None):
        r"""
        Find an edge in ``E`` equivalent to ``e0``.

        INPUT:

        - ``e0`` -- a 2x2 matrix in `\ZZ_p` representing an
            edge in the Bruhat-Tits tree

        - ``E`` -- list (default: ``None``); if a list of Edge is given,
            restrict the search to the vertices in ``E``. Otherwise
            use all the edges in a fundamental domain.

        - ``valuation`` -- integer (default: ``None``); the valuation
            of the determinant of ``e0``, if known (otherwise it is
            calculated).

        OUTPUT:

        A pair ``g``, ``e``, where ``e`` is an Edge in ``E``
        equivalent to ``e0``, and ``g`` is such that `g\cdot e_0= e`.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: M = Matrix(ZZ,2,2,[1,3,2,7])
            sage: M.set_immutable()
            sage: X._find_equivalent_edge(M)[-1] in X.get_edge_list()
            True
        """
        try:
            return self._cached_edges[e0]
        except KeyError:
            pass
        if valuation is None:
            valuation = e0.determinant().valuation(self._p)
        parity = valuation % 2
        if E is None:
            if parity == 0:
                E = self._edge_list
            else:
                E = [e.opposite for e in self._edge_list]
        for e in E:
            if e.parity != parity:
                continue
            g = self._are_equivalent(e.rep, e0, True, valuation + e.valuation)
            if g is not None:
                self._cached_edges[e0] = (g, e)
                return g, e
        return 0, None

    def fundom_rep(self, v1):
        r"""
        Find an equivalent vertex in the fundamental domain.

        INPUT:

        - ``v1`` -- a 2x2 matrix representing a normalized vertex

        OUTPUT: a ``Vertex`` equivalent to ``v1``, in the fundamental domain

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: M = Matrix(ZZ,2,2,[1,3,2,7])
            sage: M.set_immutable()
            sage: X.fundom_rep(M)
            Vertex of Bruhat-Tits tree for p = 3
        """
        try:
            tmp = self._cached_paths[v1]
            return tmp
        except KeyError:
            pass
        chain, v = self._BT.find_path(v1, self.get_vertex_dict())
        while chain:
            v0 = chain.pop()
            V = [e.target for e in v.leaving_edges]
            g, v = self._find_equivalent_vertex(v0, V)
            if v is None:
                print('Given vertex:', v0)
                print('Not equivalent to any existing vertex in the list:')
                if V is not None:
                    print([ve.label for ve in V])
                assert 0  # what the hell is that ?
            self._cached_paths[v0] = v
        return v

    def _find_lattice(self, v1, v2, as_edges, m):
        r"""
        Find the lattice attached to the pair ``v1``,``v2``.

        INPUT:

        - ``v1``, ``v2`` -- 2x2 matrices; they represent either a pair
          of normalized vertices or a pair of normalized edges

        - ``as_edges`` -- boolean; if ``True``, the inputs will be
          considered as edges instead of vertices

        - ``m`` -- integer; the valuation of the determinant of
          ``v1``*``v2``

        OUTPUT: a 4x4 integer matrix whose columns encode a lattice and a 4x4 integer matrix encoding a quadratic form

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,17)
            sage: X._find_lattice(Matrix(ZZ,2,2,[1,2,3,4]),Matrix(ZZ,2,2,[3,2,1,5]), True,0)
            (
            [1 0 0 0]  [138 204 -35 102]
            [2 3 0 0]  [204 306 -51 153]
            [0 0 1 0]  [-35 -51  12 -34]
            [0 0 0 1], [102 153 -34 102]
            )
        """
        if as_edges:
            X = self._Xe
        else:
            X = self._Xv
        if m + 1 > self._prec:
            self.get_embedding_matrix(prec=m + 1)
        v1adj = v1.adjugate()
        R = self._Mat_44
        vecM = [v2 * X[ii] * v1adj for ii in range(4)]
        M = self._Iotainv * R([[vecM[ii][jj, kk] for ii in range(4)]
                               for jj in range(2) for kk in range(2)])
        M = M.augment(R(self._pN)).transpose()
        E = M.echelon_form().submatrix(0, 0, 4, 4)
        Et = E.transpose()
        return Et, E * self.get_eichler_order_quadmatrix() * Et

    def _stabilizer(self, e, as_edge=True):
        r"""
        Find the stabilizer of an edge or vertex.

        INPUT:

        - ``e`` -- a 2x2 matrix representing an edge or vertex

        - ``as_edge`` -- boolean (default: ``True``); determines whether
          ``e`` is treated as an edge or vertex

        OUTPUT:

        A list of data describing the (finite) stabilizing subgroup
        of ``e``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,7)
            sage: X._stabilizer(Matrix(ZZ,2,2,[3,8,2,9]))[0][2]
            False
        """
        p = self._p
        m = e.determinant().valuation(p)
        twom = 2 * m
        E, A = self._find_lattice(e, e, as_edge, twom)
        n_units = len(self.get_units_of_order())
        # Using PARI to get the shortest vector in the lattice (via LLL)
        # We used to pass qfminim flag = 2
        mat = pari('qfminim(%s,,%s,flag = 2)' % (A.__pari__(), 2 * n_units))[2].sage().transpose()
        n_vecs = mat.nrows()
        stabs = []
        for jj in range(n_vecs):
            vect = mat.row(jj).row()
            vec = vect.transpose()
            nrd = Integer((vect * A * vec)[0, 0] / 2)
            if nrd == p ** twom:
                g, ans = self._nebentype_check(vec, twom, E, A, flag=0)
                if ans:
                    x = self._conv(g.transpose())
                    g.set_immutable()
                    stabs.append([g, m, x != p ** m])
        if len(stabs) <= 1:
            return [[self.B_one(), 0, False]]
        else:
            return stabs

    def _nebentype_check(self, vec, twom, E, A, flag=2):
        r"""
        Check if a quaternion maps into a subgroup of matrices
        determined by a nontrivial Dirichlet character (associated to
        self). If `N^+ = 1` then the condition is trivially satisfied.

        INPUT:

        - ``vec`` -- 4x1 integer matrix; it encodes the quaternion to
          test in the basis defined by the columns of E

        - ``twom`` -- integer

        - ``E`` -- 4x4 integer matrix; its columns should form a
          basis for an order in the quaternion algebra

        - ``A`` -- 4x4 integer matrix; it encodes the quadratic form on the
          order defined by the columns of E

        - ``flag`` -- integer (default: 0); Passed to Pari for finding
          minimal elements in a positive definite lattice

        OUTPUT:

        A pair consisting of a quaternion (represented by a 4x1 column
        matrix) and a boolean saying whether the quaternion is in the
        subgroup of `M_2(\Qp)` determined by the Dirichlet
        character. Note that if `N^+` is trivial then this function
        always outputs true.

        EXAMPLES::

            sage: f = DirichletGroup(6)[1]
            sage: X = BruhatTitsQuotient(3,2,1,f)
            sage: e = Matrix(ZZ,2,2,[1,2,5,7])
            sage: m = e.determinant().valuation(3)
            sage: twom = 2*m
            sage: E,A = X._find_lattice(e,e,True,twom)
            sage: X._nebentype_check(E**(-1)*Matrix(ZZ,4,1,[1,0,0,0]),twom,E,A)
            (
            [1]
            [0]
            [0]
            [0], True
            )
        """
        if not self._use_magma or len(self._extra_level) == 0:
            return E * vec, True
        m = ZZ(twom / 2)
        mat = pari('qfminim(%s,,%s,flag = %s)' % (A.__pari__(), 1000, flag))[2].sage().transpose()
        n_vecs = mat.nrows()
        p = self._p
        pinv = Zmod(self._character.modulus())(p) ** -1
        for jj in range(n_vecs):
            vect = mat.row(jj).row()
            vec = vect.transpose()
            nrd = Integer((vect * A * vec)[0, 0] / 2)
            if nrd == p ** twom:
                g = E * vec
                if prod([self._character(ZZ(pinv ** m * (v * g)[0, 0]))
                         for v in self.get_extra_embedding_matrices()]) == 1:
                    return g, True
        return None, False

    def _are_equivalent(self, v1, v2, as_edges=False, twom=None,
                        check_parity=False):
        r"""
        Determine whether two vertices (or edges) of the
        Bruhat-Tits tree are equivalent under the arithmetic group in
        question. The computation boils down to an application of the
        LLL short-vector algorithm to a particular lattice; for
        details see [FM2014]_.

        INPUT:

        - ``v1``, ``v2`` -- two 2x2 integral matrices representing
          either vertices or edges

        - ``as_edges`` -- boolean (default: ``False``); whether the
          matrices should be interpreted as edges (if ``True``), or as
          vertices (if ``False``)

        - ``twom`` -- integer (default: ``None``); if specified,
          indicates the valuation of the determinant of ``v1``
          `\times` ``v2``.

        OUTPUT:

        If the objects are equivalent, returns an element of
        the arithmetic group `\Gamma` that takes ``v1`` to ``v2``.
        Otherwise returns ``False``.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,5)
            sage: M1 = Matrix(ZZ,2,2,[88,3,1,1])
            sage: M1.set_immutable()
            sage: X._are_equivalent(M1,M1) == False
            False
            sage: M2 = Matrix(ZZ,2,2,[1,2,8,1]); M2.set_immutable()
            sage: X._are_equivalent(M1,M2, as_edges=True)
            sage: X._are_equivalent(M1,M2) == False
            False
        """
        try:
            return self._cached_equivalent[(v1, v2, as_edges)]
        except KeyError:
            pass
        p = self._p
        if twom is None:
            twom = v1.determinant().valuation(p) + v2.determinant().valuation(p)
        if check_parity:
            if twom % 2:
                self._cached_equivalent[(v1, v2, as_edges)] = None
                return None
        E, A = self._find_lattice(v1, v2, as_edges, twom)
        # Using PARI to get the shortest vector in the lattice (via LLL)
        vec = pari('qfminim(%s,,1,flag = 2)' % (A.__pari__()))[2].sage()

        vect = vec.transpose()
        nrd = Integer((vect * A * vec)[0, 0] / 2)
        if nrd == p ** twom:
            g, ans = self._nebentype_check(vec, twom, E, A)
            if ans:
                m = Integer(twom / 2)
                g.set_immutable()
                self._cached_equivalent[(v1, v2, as_edges)] = (g, m)
                return (g, m)
        self._cached_equivalent[(v1, v2, as_edges)] = None
        return None

    def _compute_exact_splitting(self):
        r"""
        Use Magma to calculate a splitting of the order into
        the Matrix algebra with coefficients in an appropriate
        number field.

        TESTS::

            sage: X = BruhatTitsQuotient(3,23,use_magma=True) # optional - magma
            sage: X._compute_exact_splitting() # optional - magma
        """
        # A = self.get_quaternion_algebra()
        R = self.get_maximal_order(magma=True)
        f = R.MatrixRepresentation()
        self._FF = NumberField(f.Codomain().BaseRing().DefiningPolynomial().sage(), 'a')
        allmats = []
        verbose('Calling magma, compute exact splitting')
        for kk in range(4):
            xseq = self._magma('%s(%s)' % (f.name(), R.gen(kk + 1).name())).ElementToSequence()
            allmats.append(Matrix(self._FF, 2, 2, [self._FF([QQ(xseq[ii + 1][jj + 1]) for jj in range(2)]) for ii in range(4)]))
        self._Iota_exact = Matrix(self._FF, 4, 4, [self._FF(allmats[kk][ii, jj]) for ii in range(2) for jj in range(2) for kk in range(4)])

    def _init_order(self):
        r"""
        Initialize the order of the quaternion algebra. Here we
        possibly use Magma to split it.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(3,23)
            sage: X._init_order()
        """
        if self._order_is_initialized:
            return
        if self._use_magma:
            verbose('Calling magma, init_order')
            A = self._magma.QuaternionAlgebra(self._Nminus)
            g = A.gens()
            # We store the order because we need to split it
            OMaxmagma = A.QuaternionOrder(1)
            Omagma = OMaxmagma.Order(self._Nplus)
            OBasis = Omagma.Basis()
            self._A = QuaternionAlgebra((g[0] ** 2).sage(), (g[1] ** 2).sage())
            i, j, k = self._A.gens()
            v = [1] + list(self._A.gens())
            self._B = [self._A(sum([OBasis[tt + 1][rr + 1].sage() * v[rr]
                                    for rr in range(4)])) for tt in range(4)]
            self._O = self._A.quaternion_order(self._B)
            self._Omagma = Omagma
            self._OMaxmagma = OMaxmagma
        else:
            # Note that we can't work with non-maximal orders in sage
            assert self._Nplus == 1
            self._A = QuaternionAlgebra(self._Nminus)
            v = [1] + list(self._A.gens())
            self._O = self._A.maximal_order()
            self._OMax = self._O
            OBasis = self._O.basis()
            self._B = [self._A(OBasis[tt]) for tt in range(4)]

        self._OQuadForm = QuadraticForm(self._Mat_44([(self._B[ii] * self._B[jj].conjugate()).reduced_trace() for ii in range(4) for jj in range(4)]))
        self._OM = self._OQuadForm.matrix()
        self._BB = Matrix(QQ, 4, 4, [[self._B[ii][jj] for ii in range(4)]
                                     for jj in range(4)]).inverse()
        self._order_is_initialized = True
        return

    def B_one(self):
        r"""
        Return the coordinates of `1` in the basis for the
        quaternion order.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(7,11)
            sage: v,pow = X.B_one()
            sage: X._conv(v) == 1
            True
        """
        try:
            return self._B_one
        except AttributeError:
            O = self.get_eichler_order_basis()
            self._B_one = (Matrix(ZZ, 4, 1, Matrix(QQ, 4, 4, [list(x) for x in O]).transpose().inverse().column(0).list()), 0)
            return self._B_one

    def _conv(self, v):
        r"""
        Return a quaternion having coordinates in the fixed
        basis for the order given by ``v``.

        OUTPUT: a quaternion

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: A = X.get_quaternion_algebra()
            sage: i,j,k = A.gens()
            sage: B = X.get_eichler_order_basis()
            sage: X._conv([1,2,3,4]) == B[0]+2*B[1]+3*B[2]+4*B[3]
            True
        """
        if hasattr(v, "list"):
            v = v.list()
        B = self.get_eichler_order_basis()
        return sum([v[i] * B[i] for i in range(4)])

    @cached_method
    def _find_elements_in_order(self, norm, trace=None, primitive=False):
        r"""
        Return elements in the order of the quaternion algebra
        of specified reduced norm. One may optionally choose to
        specify the reduced trace.

        INPUT:

        - ``norm`` -- integer; the required reduced norm

        - ``trace`` -- integer (default: ``None``); if specified, returns
          elements only reduced trace ``trace``

        - ``primitive`` -- boolean (default: ``False``); if ``True``, return only
          elements that cannot be divided by `p`

        EXAMPLES::

            sage: X = BruhatTitsQuotient(5,7)
            sage: X._find_elements_in_order(23)
            [[2, 9, -1, -5], [0, 8, 0, -5], [-2, 9, 1, -5], [6, 7, -3, -4], [2, 5, -1, -4], [0, 6, -1, -4], [0, 8, -1, -4], [2, 9, -1, -4], [-2, 5, 1, -4], [0, 6, 1, -4], [0, 8, 1, -4], [-2, 9, 1, -4], [-6, 7, 3, -4], [7, 6, -4, -3], [7, 6, -3, -3], [6, 7, -3, -3], [0, 8, 0, -3], [-7, 6, 3, -3], [-6, 7, 3, -3], [-7, 6, 4, -3], [0, 1, -1, -2], [0, 6, -1, -2], [0, 1, 1, -2], [0, 6, 1, -2], [9, 2, -5, -1], [6, 0, -4, -1], [8, 0, -4, -1], [5, 2, -4, -1], [9, 2, -4, -1], [1, 0, -2, -1], [6, 0, -2, -1], [0, -1, -1, -1], [-1, 0, -1, -1], [5, 2, -1, -1], [2, 5, -1, -1], [0, -1, 1, -1], [1, 0, 1, -1], [-5, 2, 1, -1], [-2, 5, 1, -1], [-6, 0, 2, -1], [-1, 0, 2, -1], [-8, 0, 4, -1], [-6, 0, 4, -1], [-9, 2, 4, -1], [-5, 2, 4, -1], [-9, 2, 5, -1], [8, 0, -5, 0], [8, 0, -3, 0]]
            sage: list(X._find_elements_in_order(23,1))
            [[1, 0, -2, -1], [1, 0, 1, -1]]
        """
        OQuadForm = self.get_eichler_order_quadform()
        if norm > 10 ** 3:
            verbose('Warning: norm (= %s) is quite large, this may take some time!' % norm)
        V = OQuadForm.vectors_by_length(norm)[norm]
        W = V if not primitive else (v for v in V
                                     if any(vi % self._p for vi in v))
        return W if trace is None else (v for v in W
                                        if self._conv(v).reduced_trace() == trace)

    def _compute_quotient(self, check=True):
        r"""
        Compute the quotient graph.

        INPUT:

        - ``check`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,2)
            sage: X.get_graph() # indirect doctest
            Multi-graph on 2 vertices

            sage: X = BruhatTitsQuotient(17,19)
            sage: X.get_graph() # indirect doctest
            Multi-graph on 4 vertices

        The following examples require magma::

            sage: X = BruhatTitsQuotient(5,7,12) # optional - magma
            sage: X.get_graph()          # optional - magma
            Multi-graph on 24 vertices
            sage: len(X._edge_list)      # optional - magma
            72

            sage: X = BruhatTitsQuotient(2,3,5)  # optional - magma
            sage: X.get_graph() # optional - magma
            Multi-graph on 4 vertices

            sage: X = BruhatTitsQuotient(2,3,35) # optional - magma
            sage: X.get_graph() # optional - magma
            Multi-graph on 16 vertices

            sage: X = BruhatTitsQuotient(53,11,2) # optional - magma
            sage: X.get_graph() # optional - magma
            Multi-graph on 6 vertices

            sage: X = BruhatTitsQuotient(2,13,9) # optional - magma
            sage: X.get_graph() # optional - magma
            Multi-graph on 24 vertices

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu
        """
        nontorsion_generators = set()
        genus = self.genus()
        num_verts = 0
        num_edges = 0
        self.get_embedding_matrix(prec=3)
        p = self._p
        v0 = Vertex(p, num_verts, self._Mat_22([1, 0, 0, 1]),
                    determinant=1, valuation=0)
        V = deque([v0])
        S = Graph(0, multiedges=True, weighted=True)  # noqa:F821
        Sfun = Graph(0)  # noqa:F821
        edge_list = []
        vertex_list = [v0]
        num_edges = 0
        num_verts += 1
        # total_verts = self.get_num_verts()
        # total_edges = genus + total_verts -1
        while len(V):
            v = V.popleft()
            E = self._BT.leaving_edges(v.rep)

            verbose('V = %s, E = %s, G = %s (target = %s), lenV = %s' % (num_verts, num_edges, 1 + num_edges - num_verts, genus, len(V)))
            for e in E:
                edge_det = e.determinant()
                edge_valuation = edge_det.valuation(p)

                g, e1 = self._find_equivalent_edge(e, v.leaving_edges,
                                                   valuation=edge_valuation)

                if e1 is not None:  # The edge is old. We just update the links
                    e1.links.append(g)
                    target = self._BT.target(e)
                    if e1.parity == 0:
                        Sfun.add_edge(v.rep, target, label=e1.label)
                    else:
                        Sfun.add_edge(v.rep, target, label=e1.opposite.label)

                    Sfun.set_vertex(target, e1.target)
                else:  # The edge is new.
                    target = self._BT.target(e)
                    target.set_immutable()
                    new_det = target.determinant()
                    new_valuation = new_det.valuation(p)
                    # new_parity = new_valuation % 2
                    g1, v1 = self._find_equivalent_vertex(target, V, valuation=new_valuation)
                    if v1 is None:
                        # The vertex is also new
                        v1 = Vertex(p, num_verts, target, determinant=new_det,
                                    valuation=new_valuation)
                        vertex_list.append(v1)
                        num_verts += 1
                        # Add the vertex to the list of pending vertices
                        V.append(v1)
                    else:
                        nontorsion_generators.add(g1[0])

                    # Add the edge to the list
                    new_e = Edge(p, num_edges, e, v, v1, determinant=edge_det,
                                 valuation=edge_valuation)
                    new_e.links.append(self.B_one())
                    Sfun.add_edge(v.rep, target, label=num_edges)
                    Sfun.set_vertex(target, v1)

                    # Add the edge to the graph
                    S.add_edge(v.rep, v1.rep, num_edges)
                    S.set_vertex(v.rep, v)
                    S.set_vertex(v1.rep, v1)

                    # Find the opposite edge
                    opp = self._BT.opposite(e)
                    # opp_det = opp.determinant()
                    new_e_opp = Edge(p, num_edges, opp, v1, v, opposite=new_e)
                    new_e.opposite = new_e_opp

                    if new_e.parity == 0:
                        edge_list.append(new_e)
                    else:
                        edge_list.append(new_e_opp)

                    v.leaving_edges.append(new_e)
                    v.entering_edges.append(new_e_opp)
                    v1.entering_edges.append(new_e)
                    v1.leaving_edges.append(new_e_opp)
                    num_edges += 1
        computed_genus = Integer(1 - len(vertex_list) + num_edges)
        if check:
            if computed_genus != genus:
                print('You found a bug! Please report!')
                print('Computed genus =', computed_genus)
                print('Theoretical genus =', genus)
                raise RuntimeError
            if self.get_num_verts() != len(vertex_list):
                raise RuntimeError('Number of vertices different '
                                   'from expected.')

        self._nontorsion_generators = nontorsion_generators
        self._boundary = {vv.rep: vv for vv in vertex_list}
        self._edge_list = edge_list
        self._vertex_list = vertex_list
        self._num_edges = num_edges
        self._S = S
        self._Sfun = Sfun

    def harmonic_cocycle_from_elliptic_curve(self, E, prec=None):
        r"""
        Return a harmonic cocycle with the same Hecke eigenvalues as ``E``.

        Given an elliptic curve `E` having a conductor `N` of the form `pN^-N^+`,
        return the harmonic cocycle over ``self`` which is attached to ``E`` via
        modularity. The result is only well-defined up to scaling.

        INPUT:

        - ``E`` -- an elliptic curve over the rational numbers

        - ``prec`` -- (default: ``None``) if specified, the harmonic cocycle will take values
          in `\QQ_p` with precision ``prec``. Otherwise it will take values in `\ZZ`.

        OUTPUT: a harmonic cocycle attached via modularity to the given elliptic curve

        EXAMPLES::

            sage: E = EllipticCurve('21a1')
            sage: X = BruhatTitsQuotient(7,3)
            sage: f = X.harmonic_cocycle_from_elliptic_curve(E,10)
            sage: T29 = f.parent().hecke_operator(29)
            sage: T29(f) == E.ap(29) * f
            True
            sage: E = EllipticCurve('51a1')
            sage: X = BruhatTitsQuotient(3,17)
            sage: f = X.harmonic_cocycle_from_elliptic_curve(E,20)
            sage: T31 = f.parent().hecke_operator(31)
            sage: T31(f) == E.ap(31) * f
            True
        """
        from .pautomorphicform import BruhatTitsHarmonicCocycles
        M = BruhatTitsHarmonicCocycles(self, 2, prec=prec)
        q = ZZ.one()
        F = E.base_ring()
        try:
            N = ZZ(E.conductor())
        except TypeError:
            try:
                N = E.conductor().norm()
            except ValueError:
                N = E.conductor().norm(QQ)
        N1 = self.level() * self.Nplus()
        K = M.base_ring() ** M.dimension()
        while K.dimension() != 1:
            q = q.next_prime()
            if N % q == 0 or N1 % q == 0:
                continue
            if F == QQ:
                Eap = E.ap(q)
            else:
                Q = F(q).factor()[0][0]
                Eap = ZZ(Q.norm() + 1 - E.reduction(Q).count_points())
            K1 = (M.hecke_matrix(q).transpose() - Eap).right_kernel()
            K = K.intersection(K1)
        col = [ZZ(o) for o in K.matrix().list()]
        return sum([a * M.gen(i) for i, a in enumerate(col) if a != 0], M(0))

    def harmonic_cocycles(self, k, prec=None, basis_matrix=None, base_field=None):
        r"""
        Compute the space of harmonic cocycles of a given even weight ``k``.

        INPUT:

        - ``k`` -- integer; the weight. It must be even.

        - ``prec`` -- integer (default: ``None``); if specified, the
          precision for the coefficient module

        - ``basis_matrix`` -- a matrix (default: ``None``)

        - ``base_field`` -- a ring (default: ``None``)

        OUTPUT: a space of harmonic cocycles

        EXAMPLES::

            sage: X = BruhatTitsQuotient(31,7)
            sage: H = X.harmonic_cocycles(2,prec=10)
            sage: H
            Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_31) with discriminant 7 and level 1
            sage: H.basis()[0]
            Harmonic cocycle with values in Sym^0 Q_31^2
        """
        from .pautomorphicform import BruhatTitsHarmonicCocycles
        return BruhatTitsHarmonicCocycles(self, k, prec=prec, basis_matrix=basis_matrix, base_field=base_field)

    def padic_automorphic_forms(self, U, prec=None, t=None, R=None, overconvergent=False):
        r"""
        The module of (quaternionic) `p`-adic automorphic forms over ``self``.

        INPUT:

        - ``U`` -- a distributions module or an integer. If ``U`` is a
          distributions module then this creates the relevant space of
          automorphic forms. If ``U`` is an integer then the coefficients
          are the (`U-2`)nd power of the symmetric representation of
          `GL_2(\QQ_p)`.

        - ``prec`` -- a precision (default: ``None``). if not ``None`` should
          be a positive integer

        - ``t`` -- (default: ``None``) the number of additional moments to store. If ``None``, determine
          it automatically from ``prec``, ``U`` and the ``overconvergent`` flag

        - ``R`` -- (default: ``None``) if specified, coefficient field of the automorphic forms.
          If not specified it defaults to the base ring of the distributions ``U``, or to `\QQ_p`
          with the working precision ``prec``.

        - ``overconvergent`` -- boolean (default: ``False``); if ``True``, will construct overconvergent
          `p`-adic automorphic forms. Otherwise it constructs the finite dimensional space of
          `p`-adic automorphic forms which is isomorphic to the space of harmonic cocycles.

        EXAMPLES::

            sage: X = BruhatTitsQuotient(11,5)
            sage: X.padic_automorphic_forms(2,prec=10)
            Space of automorphic forms on Quotient of the Bruhat Tits tree of GL_2(QQ_11) with discriminant 5 and level 1 with values in Sym^0 Q_11^2
        """
        from .pautomorphicform import pAdicAutomorphicForms
        return pAdicAutomorphicForms(self, U, prec=prec, t=t, R=R, overconvergent=overconvergent)
