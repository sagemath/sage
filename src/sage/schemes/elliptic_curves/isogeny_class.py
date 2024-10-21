# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.number_field
r"""
Isogeny class of elliptic curves over number fields

AUTHORS:

- David Roe (2012-03-29) -- initial version.
- John Cremona (2014-08) -- extend to number fields.
"""

##############################################################################
#       Copyright (C) 2012-2014 David Roe <roed.math@gmail.com>
#                          John Cremona <john.cremona@gmail.com>
#                          William Stein <wstein@gmail.com>
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
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp_method, richcmp
import sage.databases.cremona
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.flatten import flatten
from sage.misc.cachefunc import cached_method
from sage.schemes.elliptic_curves.ell_field import EllipticCurve_field
from sage.schemes.elliptic_curves.ell_number_field import EllipticCurve_number_field


@richcmp_method
class IsogenyClass_EC(SageObject):
    r"""
    Isogeny class of an elliptic curve.

    .. NOTE::

        The current implementation chooses a curve from each isomorphism
        class in the isogeny class. Over `\QQ` this is a unique reduced
        minimal model in each isomorphism class.  Over number fields the
        model chosen may change in future.
    """

    def __init__(self, E, label=None, empty=False):
        r"""
        Over `\QQ` we use curves since minimal models exist and there
        is a canonical choice of one.

        INPUT:

        - ``label`` -- string or ``None``, a Cremona or LMFDB label, used
          in printing; ignored if base field is not `\QQ`

        EXAMPLES::

            sage: cls = EllipticCurve('1011b1').isogeny_class()
            sage: print("\n".join(repr(E) for E in cls.curves))
            Elliptic Curve defined by y^2 + x*y = x^3 - 8*x - 9 over Rational Field
            Elliptic Curve defined by y^2 + x*y = x^3 - 23*x + 30 over Rational Field
        """
        self.E = E
        self._label = label
        if not empty:
            self._compute()

    def __len__(self):
        """
        The number of curves in the class.

        EXAMPLES::

            sage: E = EllipticCurve('15a')
            sage: len(E.isogeny_class()) # indirect doctest
            8
        """
        return len(self.curves)

    def __iter__(self):
        """
        Iterator over curves in the class.

        EXAMPLES::

            sage: E = EllipticCurve('15a')
            sage: all(C.conductor() == 15 for C in E.isogeny_class()) # indirect doctest
            True
        """
        return iter(self.curves)

    def __getitem__(self, i):
        """
        Return the `i`-th curve in the class.

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: iso = E.isogeny_class(order='lmfdb') # orders lexicographically on a-invariants
            sage: iso[2] == E # indirect doctest
            True
        """
        return self.curves[i]

    def index(self, C):
        """
        Return the index of a curve in this class.

        INPUT:

        - ``C`` -- an elliptic curve in this isogeny class

        OUTPUT:

        - ``i`` -- integer so that the ``i`` th curve in the class
          is isomorphic to ``C``

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: iso = E.isogeny_class(order='lmfdb') # orders lexicographically on a-invariants
            sage: iso.index(E.short_weierstrass_model())
            2
        """
        # This will need updating once we start talking about curves
        # over more general number fields
        if not isinstance(C, EllipticCurve_number_field):
            raise ValueError("x not in isogeny class")
        for i, E in enumerate(self.curves):
            if C.is_isomorphic(E):
                return i
        raise ValueError("%s is not in isogeny class %s" % (C,self))

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        If they are different, compares the sorted underlying lists of
        curves.

        Note that two isogeny classes with different orderings will
        compare as the same.  If you want to include the ordering,
        just compare the list of curves.

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: EE = EllipticCurve('990j4')
            sage: E.isogeny_class() == EE.isogeny_class() # indirect doctest
            True
        """
        if isinstance(other, IsogenyClass_EC):
            return richcmp(sorted(e.a_invariants() for e in self.curves),
                           sorted(f.a_invariants() for f in other.curves), op)
        return NotImplemented

    def __hash__(self):
        """
        Hash is based on the a-invariants of the sorted list of
        minimal models.

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: C = E.isogeny_class()
            sage: hash(C) == hash(tuple(sorted([curve.a_invariants() for curve in C.curves]))) # indirect doctest
            True
        """
        try:
            return self._hash
        except AttributeError:
            self._hash = hash(tuple(sorted(E.a_invariants() for E in self.curves)))
            return self._hash

    def _repr_(self):
        r"""
        The string representation of this isogeny class.

        .. NOTE::

            Over `\QQ`, the string representation depends on whether an
            LMFDB or Cremona label for the curve is known when this
            isogeny class is constructed.  Over general number fields,
            instead of labels the representation uses that of the curve
            initially used to create the class.

        EXAMPLES:

        If the curve is constructed from an LMFDB label then that
        label is used::

            sage: E = EllipticCurve('462.f3')
            sage: E.isogeny_class() # indirect doctest
            Elliptic curve isogeny class 462.f

        If the curve is constructed from a Cremona label then that
        label is used::

            sage: E = EllipticCurve('990j1')
            sage: E.isogeny_class()
            Elliptic curve isogeny class 990j

        Otherwise, including curves whose base field is not `\QQ`,the
        representation is determined from the curve used to create the
        class::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.isogeny_class()
            Isogeny class of Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1]); E
            Elliptic Curve defined by y^2 = x^3 + 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
            sage: C = E.isogeny_class()
            sage: C
            Isogeny class of Elliptic Curve defined by y^2 = x^3 + 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
            sage: C.curves
            [Elliptic Curve defined by y^2 = x^3 + (-27) over Number Field in i with defining polynomial x^2 + 1 with i = 1*I,
             Elliptic Curve defined by y^2 = x^3 + 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I,
             Elliptic Curve defined by y^2 + (i+1)*x*y + (i+1)*y = x^3 + i*x^2 + (-i+3)*x + 4*i over Number Field in i with defining polynomial x^2 + 1 with i = 1*I,
             Elliptic Curve defined by y^2 + (i+1)*x*y + (i+1)*y = x^3 + i*x^2 + (-i+33)*x + (-58*i) over Number Field in i with defining polynomial x^2 + 1 with i = 1*I]
        """
        if self._label:
            return "Elliptic curve isogeny class %s" % (self._label)
        else:
            return "Isogeny class of %r" % (self.E)

    def __contains__(self, x):
        """
        INPUT:

        - ``x`` -- a Python object

        OUTPUT: boolean; ``True`` iff ``x`` is an elliptic curve in this
        isogeny class

        .. NOTE::

            If the input is isomorphic but not identical to a curve in
            the class, then ``False`` will be returned.

        EXAMPLES::

            sage: cls = EllipticCurve('15a3').isogeny_class()
            sage: E = EllipticCurve('15a7'); E in cls
            True
            sage: E.short_weierstrass_model() in cls
            True
        """
        if not isinstance(x, EllipticCurve_field):
            return False
        return any(x.is_isomorphic(y) for y in self.curves)

    @cached_method
    def matrix(self, fill=True):
        """
        Return the matrix whose entries give the minimal degrees of
        isogenies between curves in this class.

        INPUT:

        - ``fill`` -- boolean (default: ``True``); if ``False`` then the
          matrix will contain only zeros and prime entries. If ``True`` it
          will fill in the other degrees.

        EXAMPLES::

            sage: isocls = EllipticCurve('15a3').isogeny_class()
            sage: isocls.matrix()
            [ 1  2  2  2  4  4  8  8]
            [ 2  1  4  4  8  8 16 16]
            [ 2  4  1  4  8  8 16 16]
            [ 2  4  4  1  2  2  4  4]
            [ 4  8  8  2  1  4  8  8]
            [ 4  8  8  2  4  1  2  2]
            [ 8 16 16  4  8  2  1  4]
            [ 8 16 16  4  8  2  4  1]
            sage: isocls.matrix(fill=False)
            [0 2 2 2 0 0 0 0]
            [2 0 0 0 0 0 0 0]
            [2 0 0 0 0 0 0 0]
            [2 0 0 0 2 2 0 0]
            [0 0 0 2 0 0 0 0]
            [0 0 0 2 0 0 2 2]
            [0 0 0 0 0 2 0 0]
            [0 0 0 0 0 2 0 0]
        """
        if self._mat is None:
            self._compute_matrix()
        mat = self._mat
        if fill and mat[0, 0] == 0:
            from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix
            mat = fill_isogeny_matrix(mat)
        if not fill and mat[0, 0] == 1:
            from sage.schemes.elliptic_curves.ell_curve_isogeny import unfill_isogeny_matrix
            mat = unfill_isogeny_matrix(mat)
        return mat

    @cached_method
    def qf_matrix(self):
        """
        Return the array whose entries are quadratic forms
        representing the degrees of isogenies between curves in this
        class (CM case only).

        OUTPUT:

        a `2x2` array (list of lists) of list, each of the form [2] or
        [2,1,3] representing the coefficients of an integral quadratic
        form in 1 or 2 variables whose values are the possible isogeny
        degrees between the i'th and j'th curve in the class.

        EXAMPLES::

            sage: pol = PolynomialRing(QQ,'x')([1,0,3,0,1])
            sage: K.<c> = NumberField(pol)
            sage: j = 1480640 + 565760*c^2
            sage: E = EllipticCurve(j=j)
            sage: C = E.isogeny_class()
            sage: C.qf_matrix()
            [[[1], [2, 2, 3]], [[2, 2, 3], [1]]]
        """
        if self._qfmat is None:
            raise ValueError("qf_matrix only defined for isogeny classes with rational CM")
        else:
            return self._qfmat

    @cached_method
    def isogenies(self, fill=False):
        r"""
        Return a list of lists of isogenies and 0s, corresponding to
        the entries of :meth:`matrix`

        INPUT:

        - ``fill`` -- boolean (default: ``False``); whether to only return
          prime degree isogenies.  Currently only implemented for
          ``fill=False``.

        OUTPUT:

        - a list of lists, where the ``j`` th entry of the ``i`` th list
          is either zero or a prime degree isogeny from the ``i`` th curve
          in this class to the ``j`` th curve.

        .. WARNING::

            The domains and codomains of the isogenies will have the same
            Weierstrass equation as the curves in this class, but they
            may not be identical python objects in the current
            implementation.

        EXAMPLES::

            sage: isocls = EllipticCurve('15a3').isogeny_class()
            sage: f = isocls.isogenies()[0][1]; f
            Isogeny of degree 2
              from Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 5*x + 2 over Rational Field
                to Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 80*x + 242 over Rational Field
            sage: f.domain() == isocls.curves[0] and f.codomain() == isocls.curves[1]
            True
        """
        if fill:
            raise NotImplementedError
        isogenies = self._maps
        if isogenies is None:
            self._compute_isogenies()
            isogenies = self._maps
        return isogenies

    @cached_method
    def graph(self):
        r"""
        Return a graph whose vertices correspond to curves in this
        class, and whose edges correspond to prime degree isogenies.

        .. NOTE::

            There are only finitely many possible isogeny graphs for
            curves over `\QQ` [Maz1978b].  This function tries to lay out
            the graph nicely by special casing each isogeny graph.
            This could also be done over other number fields, such as
            quadratic fields.

        .. NOTE::

            The vertices are labeled 1 to n rather than 0 to n-1 to
            match LMFDB and Cremona labels for curves over `\QQ`.

        EXAMPLES::

            sage: isocls = EllipticCurve('15a3').isogeny_class()
            sage: G = isocls.graph()
            sage: sorted(G._pos.items())
            [(1, [-0.8660254, 0.5]), (2, [-0.8660254, 1.5]), (3, [-1.7320508, 0]),
             (4, [0, 0]), (5, [0, -1]), (6, [0.8660254, 0.5]),
             (7, [0.8660254, 1.5]), (8, [1.7320508, 0])]
        """
        from sage.graphs.graph import Graph

        if self.E.base_field() is not QQ:
            M = self.matrix(fill=False)
            n = len(self)
            G = Graph(M, format='weighted_adjacency_matrix')
            D = {v: self.curves[v] for v in G.vertices(sort=False)}
            G.set_vertices(D)
            if self._qfmat:  # i.e. self.E.has_rational_cm():
                for i in range(n):
                    for j in range(n):
                        if M[i, j]:
                            G.set_edge_label(i, j, str(self._qfmat[i][j]))
            G.relabel(list(range(1, n + 1)))
            return G

        M = self.matrix(fill=False)
        n = M.nrows() # = M.ncols()
        G = Graph(M, format='weighted_adjacency_matrix')
        N = self.matrix(fill=True)
        D = {v: self.curves[v] for v in G.vertices(sort=False)}
        # The maximum degree classifies the shape of the isogeny
        # graph, though the number of vertices is often enough.
        # This only holds over Q, so this code will need to change
        # once other isogeny classes are implemented.

        if n == 1:
            # one vertex
            pass
        elif n == 2:
            # one edge, two vertices.  We align horizontally and put
            # the lower number on the left vertex.
            G.set_pos(pos={0:[-0.5,0],1:[0.5,0]})
        else:
            maxdegree = max(max(N))
            if n == 3:
                # o--o--o
                centervert = [i for i in range(3) if max(N.row(i)) < maxdegree][0]
                other = [i for i in range(3) if i != centervert]
                G.set_pos(pos={centervert:[0,0],other[0]:[-1,0],other[1]:[1,0]})
            elif maxdegree == 4:
                # o--o<8
                centervert = [i for i in range(4) if max(N.row(i)) < maxdegree][0]
                other = [i for i in range(4) if i != centervert]
                G.set_pos(pos={centervert:[0,0],other[0]:[0,1],other[1]:[-0.8660254,-0.5],other[2]:[0.8660254,-0.5]})
            elif maxdegree == 27:
                # o--o--o--o
                centers = [i for i in range(4) if list(N.row(i)).count(3) == 2]
                left = [j for j in range(4) if N[centers[0],j] == 3 and j not in centers][0]
                right = [j for j in range(4) if N[centers[1],j] == 3 and j not in centers][0]
                G.set_pos(pos={left:[-1.5,0],centers[0]:[-0.5,0],centers[1]:[0.5,0],right:[1.5,0]})
            elif n == 4:
                # square
                opp = [i for i in range(1,4) if not N[0,i].is_prime()][0]
                other = [i for i in range(1,4) if i != opp]
                G.set_pos(pos={0:[1,1],other[0]:[-1,1],opp:[-1,-1],other[1]:[1,-1]})
            elif maxdegree == 8:
                # 8>o--o<8
                centers = [i for i in range(6) if list(N.row(i)).count(2) == 3]
                left = [j for j in range(6) if N[centers[0],j] == 2 and j not in centers]
                right = [j for j in range(6) if N[centers[1],j] == 2 and j not in centers]
                G.set_pos(pos={centers[0]:[-0.5,0],left[0]:[-1,0.8660254],left[1]:[-1,-0.8660254],centers[1]:[0.5,0],right[0]:[1,0.8660254],right[1]:[1,-0.8660254]})
            elif maxdegree == 18:
                # two squares joined on an edge
                centers = [i for i in range(6) if list(N.row(i)).count(3) == 2]
                top = [j for j in range(6) if N[centers[0],j] == 3]
                bl = [j for j in range(6) if N[top[0],j] == 2][0]
                br = [j for j in range(6) if N[top[1],j] == 2][0]
                G.set_pos(pos={centers[0]:[0,0.5],centers[1]:[0,-0.5],top[0]:[-1,0.5],top[1]:[1,0.5],bl:[-1,-0.5],br:[1,-0.5]})
            elif maxdegree == 16:
                # tree from bottom, 3 regular except for the leaves.
                centers = [i for i in range(8) if list(N.row(i)).count(2) == 3]
                center = [i for i in centers if len([j for j in centers if N[i,j] == 2]) == 2][0]
                centers.remove(center)
                bottom = [j for j in range(8) if N[center,j] == 2 and j not in centers][0]
                left = [j for j in range(8) if N[centers[0],j] == 2 and j != center]
                right = [j for j in range(8) if N[centers[1],j] == 2 and j != center]
                G.set_pos(pos={center:[0,0],bottom:[0,-1],centers[0]:[-0.8660254,0.5],centers[1]:[0.8660254,0.5],left[0]:[-0.8660254,1.5],right[0]:[0.8660254,1.5],left[1]:[-1.7320508,0],right[1]:[1.7320508,0]})
            elif maxdegree == 12:
                # tent
                centers = [i for i in range(8) if list(N.row(i)).count(2) == 3]
                left = [j for j in range(8) if N[centers[0],j] == 2]
                right = []
                for i in range(3):
                    right.append([j for j in range(8) if N[centers[1],j] == 2 and N[left[i],j] == 3][0])
                G.set_pos(pos={centers[0]:[-0.75,0],centers[1]:[0.75,0],left[0]:[-0.75,1],right[0]:[0.75,1],left[1]:[-1.25,-0.75],right[1]:[0.25,-0.75],left[2]:[-0.25,-0.25],right[2]:[1.25,-0.25]})
        G.set_vertices(D)
        G.relabel(list(range(1, n + 1)))
        return G

    @cached_method
    def reorder(self, order):
        r"""
        Return a new isogeny class with the curves reordered.

        INPUT:

        - ``order`` -- ``None``, a string or an iterable over all curves in
          this class.  See
          :meth:`sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.isogeny_class`
          for more details.

        OUTPUT:

        Another :class:`IsogenyClass_EC` with the curves reordered (and
        matrices and maps changed as appropriate).

        EXAMPLES::

            sage: isocls = EllipticCurve('15a1').isogeny_class()
            sage: print("\n".join(repr(C) for C in isocls.curves))
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 5*x + 2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + 35*x - 28 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 135*x - 660 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 80*x + 242 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 110*x - 880 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 2160*x - 39540 over Rational Field
            sage: isocls2 = isocls.reorder('lmfdb')
            sage: print("\n".join(repr(C) for C in isocls2.curves))
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 2160*x - 39540 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 135*x - 660 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 110*x - 880 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 80*x + 242 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 5*x + 2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + 35*x - 28 over Rational Field
        """
        if order is None or isinstance(order, str) and order == self._algorithm:
            return self
        if isinstance(order, str):
            if order == "lmfdb":
                reordered_curves = sorted(self.curves,
                                          key=lambda E: E.a_invariants())
            else:
                reordered_curves = list(self.E.isogeny_class(algorithm=order))
        elif isinstance(order, (list, tuple, IsogenyClass_EC)):
            reordered_curves = list(order)
            if len(reordered_curves) != len(self.curves):
                raise ValueError("Incorrect length")
        else:
            raise TypeError("order parameter should be a string, list of curves or isogeny class")
        need_perm = self._mat is not None
        cpy = self.copy()
        curves = []
        perm = []
        for E in reordered_curves:
            try:
                j = self.curves.index(E)
            except ValueError:
                try:
                    j = self.curves.index(E.minimal_model())
                except ValueError:
                    raise ValueError("order does not yield a permutation of curves")
            curves.append(self.curves[j])
            if need_perm:
                perm.append(j+1)
        cpy.curves = tuple(curves)
        if need_perm:
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            perm = SymmetricGroup(len(self.curves))(perm)
            cpy._mat = perm.matrix() * self._mat * (~perm).matrix()
            if self._maps is not None:
                n = len(self._maps)
                cpy._maps = [self._maps[perm(i+1)-1] for i in range(n)]
                for i in range(n):
                    cpy._maps[i] = [cpy._maps[i][perm(jj + 1)-1]
                                    for jj in range(n)]
        else:
            cpy._mat = None
            cpy._maps = None
        return cpy


class IsogenyClass_EC_NumberField(IsogenyClass_EC):
    """
    Isogeny classes for elliptic curves over number fields.
    """
    def __init__(self, E, reducible_primes=None, algorithm='Billerey', minimal_models=True):
        r"""
        INPUT:

        - ``E`` -- an elliptic curve over a number field

        - ``reducible_primes`` -- list of integers, or ``None`` (default); if
          not ``None`` then this should be a list of primes; in computing the
          isogeny class, only composites isogenies of these degrees will be used.

        - ``algorithm`` -- string (default: ``'Billerey'``); the algorithm
          to use to compute the reducible primes.  Ignored for CM
          curves or if ``reducible_primes`` is provided.  Values are
          ``'Billerey'`` (default), ``'Larson'``, and ``'heuristic'``.

        - ``minimal_models`` -- boolean (default: ``True``); if ``True``,
          all curves in the class will be minimal or semi-minimal
          models.  Over fields of larger degree it can be expensive to
          compute these so set to ``False``.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1])
            sage: C = E.isogeny_class(); C
            Isogeny class of Elliptic Curve defined by y^2 = x^3 + 1
            over Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        The curves in the class (sorted)::

            sage: [E1.ainvs() for E1 in C]
            [(0, 0, 0, 0, -27),
             (0, 0, 0, 0, 1),
             (i + 1, i, i + 1, -i + 3, 4*i),
             (i + 1, i, i + 1, -i + 33, -58*i)]

        The matrix of degrees of cyclic isogenies between curves::

            sage: C.matrix()
            [1 3 6 2]
            [3 1 2 6]
            [6 2 1 3]
            [2 6 3 1]

        The array of isogenies themselves is not filled out but only
        contains those used to construct the class, the other entries
        containing the integer 0.  This will be changed when the
        class :class:`EllipticCurveIsogeny` allowed composition.  In
        this case we used `2`-isogenies to go from 0 to 2 and from 1
        to 3, and `3`-isogenies to go from 0 to 1 and from 2 to 3::

            sage: isogs = C.isogenies()
            sage: [((i,j), isogs[i][j].degree())
            ....:  for i in range(4) for j in range(4) if isogs[i][j] != 0]
            [((0, 1), 3),
             ((0, 3), 2),
             ((1, 0), 3),
             ((1, 2), 2),
             ((2, 1), 2),
             ((2, 3), 3),
             ((3, 0), 2),
             ((3, 2), 3)]
            sage: [((i,j), isogs[i][j].x_rational_map())
            ....:  for i in range(4) for j in range(4) if isogs[i][j] != 0]
            [((0, 1), (1/9*x^3 - 12)/x^2),
             ((0, 3), (-1/2*i*x^2 + i*x - 12*i)/(x - 3)),
             ((1, 0), (x^3 + 4)/x^2),
             ((1, 2), (-1/2*i*x^2 - i*x - 2*i)/(x + 1)),
             ((2, 1), (1/2*i*x^2 - x)/(x + 3/2*i)),
             ((2, 3), (x^3 + 4*i*x^2 - 10*x - 10*i)/(x^2 + 4*i*x - 4)),
             ((3, 0), (1/2*i*x^2 + x + 4*i)/(x - 5/2*i)),
             ((3, 2), (1/9*x^3 - 4/3*i*x^2 - 34/3*x + 226/9*i)/(x^2 - 8*i*x - 16))]

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([1+i, -i, i, 1, 0])
            sage: C = E.isogeny_class(); C
            Isogeny class of Elliptic Curve defined
             by y^2 + (i+1)*x*y + i*y = x^3 + (-i)*x^2 + x
             over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
            sage: len(C)
            6
            sage: C.matrix()
            [ 1  3  9 18  6  2]
            [ 3  1  3  6  2  6]
            [ 9  3  1  2  6 18]
            [18  6  2  1  3  9]
            [ 6  2  6  3  1  3]
            [ 2  6 18  9  3  1]
            sage: [E1.ainvs() for E1 in C]
            [(i + 1, i - 1, i, -i - 1, -i + 1),
            (i + 1, i - 1, i, 14*i + 4, 7*i + 14),
            (i + 1, i - 1, i, 59*i + 99, 372*i - 410),
            (i + 1, -i, i, -240*i - 399, 2869*i + 2627),
            (i + 1, -i, i, -5*i - 4, 2*i + 5),
            (i + 1, -i, i, 1, 0)]

        An example with CM by `\sqrt{-5}`::

            sage: pol = PolynomialRing(QQ,'x')([1,0,3,0,1])
            sage: K.<c> = NumberField(pol)
            sage: j = 1480640 + 565760*c^2
            sage: E = EllipticCurve(j=j)
            sage: E.has_cm()
            True
            sage: E.has_rational_cm()
            True
            sage: E.cm_discriminant()
            -20
            sage: C = E.isogeny_class()
            sage: len(C)
            2
            sage: C.matrix()
            [1 2]
            [2 1]
            sage: [E.ainvs() for E in C]
            [(0, 0, 0, 83490*c^2 - 147015, -64739840*c^2 - 84465260),
             (0, 0, 0, -161535*c^2 + 70785, -62264180*c^3 + 6229080*c)]
            sage: C.isogenies()[0][1]
            Isogeny of degree 2
              from Elliptic Curve defined by
                   y^2 = x^3 + (83490*c^2-147015)*x + (-64739840*c^2-84465260)
                   over Number Field in c with defining polynomial x^4 + 3*x^2 + 1
                to Elliptic Curve defined by
                   y^2 = x^3 + (-161535*c^2+70785)*x + (-62264180*c^3+6229080*c)
                   over Number Field in c with defining polynomial x^4 + 3*x^2 + 1

        TESTS::

            sage: TestSuite(C).run()
        """
        self._algorithm = "sage"
        self._reducible_primes = reducible_primes
        self._algorithm = algorithm
        self._minimal_models = minimal_models
        IsogenyClass_EC.__init__(self, E, label=None, empty=False)

    def copy(self):
        """
        Return a copy (mostly used in reordering).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1])
            sage: C = E.isogeny_class()
            sage: C2 = C.copy()
            sage: C is C2
            False
            sage: C == C2
            True
        """
        ans = IsogenyClass_EC_NumberField(self.E, reducible_primes=self._reducible_primes, algorithm=self._algorithm, minimal_models=self._minimal_models)
        # The following isn't needed internally, but it will keep
        # things from breaking if this is used for something other
        # than reordering.
        ans.curves = self.curves
        ans._mat = None
        ans._maps = None
        return ans

    def _compute(self, verbose=False):
        """
        Compute the list of curves, the matrix and prime-degree
        isogenies.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1])
            sage: C = E.isogeny_class()
            sage: C2 = C.copy()
            sage: C2._mat
            sage: C2._compute()
            sage: C2._mat
            [1 3 6 2]
            [3 1 2 6]
            [6 2 1 3]
            [2 6 3 1]

            sage: C2._compute(verbose=True)
            possible isogeny degrees: [2, 3] -actual isogeny degrees: {2, 3} -added curve #1 (degree 2)... -added tuple [0, 1, 2]... -added tuple [1, 0, 2]... -added curve #2 (degree 3)... -added tuple [0, 2, 3]... -added tuple [2, 0, 3]...... relevant degrees: [2, 3]... -now completing the isogeny class... -processing curve #1... -added tuple [1, 0, 2]... -added tuple [0, 1, 2]... -added curve #3... -added tuple [1, 3, 3]... -added tuple [3, 1, 3]... -processing curve #2... -added tuple [2, 3, 2]... -added tuple [3, 2, 2]... -added tuple [2, 0, 3]... -added tuple [0, 2, 3]... -processing curve #3... -added tuple [3, 2, 2]... -added tuple [2, 3, 2]... -added tuple [3, 1, 3]... -added tuple [1, 3, 3]...... isogeny class has size 4
            Sorting permutation = {0: 1, 1: 2, 2: 0, 3: 3}
            Matrix = [1 3 6 2]
            [3 1 2 6]
            [6 2 1 3]
            [2 6 3 1]

        TESTS:

        Check that :issue:`19030` is fixed (codomains of reverse isogenies were wrong)::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([1, i + 1, 1, -72*i + 8, 95*i + 146])
            sage: C = E.isogeny_class()
            sage: curves = C.curves
            sage: isos = C.isogenies()
            sage: isos[0][3].codomain() == curves[3]
            True
        """
        from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix
        from sage.matrix.matrix_space import MatrixSpace
        from sage.sets.set import Set
        self._maps = None

        if self._minimal_models:
            E = self.E.global_minimal_model(semi_global=True)
        else:
            E = self.E

        degs = self._reducible_primes
        if degs is None:
            self._reducible_primes = possible_isogeny_degrees(E, algorithm=self._algorithm)
            degs = self._reducible_primes
        if verbose:
            import sys
            sys.stdout.write(" possible isogeny degrees: %s" % degs)
            sys.stdout.flush()
        isogenies = E.isogenies_prime_degree(degs, minimal_models=self._minimal_models)
        if verbose:
            sys.stdout.write(" -actual isogeny degrees: %s" % Set(phi.degree() for phi in isogenies))
            sys.stdout.flush()
        # Add all new codomains to the list and collect degrees:
        curves = [E]
        ncurves = 1
        degs = []
        # tuples (i,j,l,phi) where curve i is l-isogenous to curve j via phi
        tuples = []

        def add_tup(t):
            for T in [t, [t[1], t[0], t[2], 0]]:
                if T not in tuples:
                    tuples.append(T)
                    if verbose:
                        sys.stdout.write(" -added tuple %s..." % T[:3])
                        sys.stdout.flush()

        for phi in isogenies:
            E2 = phi.codomain()
            d = ZZ(phi.degree())
            if not any(E2.is_isomorphic(E3) for E3 in curves):
                curves.append(E2)
                if verbose:
                    sys.stdout.write(" -added curve #%s (degree %s)..." % (ncurves,d))
                    sys.stdout.flush()
                add_tup([0,ncurves,d,phi])
                ncurves += 1
                if d not in degs:
                    degs.append(d)
        if verbose:
            sys.stdout.write("... relevant degrees: %s..." % degs)
            sys.stdout.write(" -now completing the isogeny class...")
            sys.stdout.flush()

        i = 1
        while i < ncurves:
            E1 = curves[i]
            if verbose:
                sys.stdout.write(" -processing curve #%s..." % i)
                sys.stdout.flush()

            isogenies = E1.isogenies_prime_degree(degs, minimal_models=self._minimal_models)

            for phi in isogenies:
                E2 = phi.codomain()
                d = phi.degree()
                js = [j for j,E3 in enumerate(curves) if E2.is_isomorphic(E3)]
                if js: # seen codomain already -- up to isomorphism
                    j = js[0]
                    if phi.codomain() != curves[j]:
                        phi = E2.isomorphism_to(curves[j]) * phi
                    assert phi.domain() == curves[i] and phi.codomain() == curves[j]
                    add_tup([i,j,d,phi])
                else:
                    curves.append(E2)
                    if verbose:
                        sys.stdout.write(" -added curve #%s..." % ncurves)
                        sys.stdout.flush()
                    add_tup([i,ncurves,d,phi])
                    ncurves += 1
            i += 1

        if verbose:
            print("... isogeny class has size %s" % ncurves)

        # key function for sorting
        if E.has_rational_cm():
            key_function = lambda E: (-E.cm_discriminant(),
                                      flatten([list(ai) for ai in E.ainvs()]))
        else:
            key_function = lambda E: flatten([list(ai) for ai in E.ainvs()])

        self.curves = sorted(curves, key=key_function)
        perm = {ind: self.curves.index(Ei)
                for ind, Ei in enumerate(curves)}
        if verbose:
            print("Sorting permutation = %s" % perm)

        mat = MatrixSpace(ZZ, ncurves)(0)
        self._maps = [[0] * ncurves for _ in range(ncurves)]
        for i, j, l, phi in tuples:
            if phi != 0:
                mat[perm[i], perm[j]] = l
                self._maps[perm[i]][perm[j]] = phi
        self._mat = fill_isogeny_matrix(mat)
        if verbose:
            print("Matrix = %s" % self._mat)

        if not E.has_rational_cm():
            self._qfmat = None
            return

        # In the CM case, we will have found some "horizontal"
        # isogenies of composite degree and would like to replace them
        # by isogenies of prime degree, mainly to make the isogeny
        # graph look better.  We also construct a matrix whose entries
        # are not degrees of cyclic isogenies, but rather quadratic
        # forms (in 1 or 2 variables) representing the isogeny
        # degrees.  For this we take a short cut: properly speaking,
        # when `\text{End}(E_1)=\text{End}(E_2)=O`, the set
        # `\text{Hom}(E_1,E_2)` is a rank `1` projective `O`-module,
        # hence has a well-defined ideal class associated to it, and
        # hence (using an identification between the ideal class group
        # and the group of classes of primitive quadratic forms of the
        # same discriminant) an equivalence class of quadratic forms.
        # But we currently only care about the numbers represented by
        # the form, i.e. which genus it is in rather than the exact
        # class.  So it suffices to find one form of the correct
        # discriminant which represents one isogeny degree from `E_1`
        # to `E_2` in order to obtain a form which represents all such
        # degrees.

        if verbose:
            print("Creating degree matrix (CM case)")

        allQs = {}  # keys: discriminants d
                    # values: lists of equivalence classes of
                    # primitive forms of discriminant d

        def find_quadratic_form(d, n):
            if d not in allQs:
                from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives

                allQs[d] = BinaryQF_reduced_representatives(d, primitive_only=True)
            # now test which of the Qs represents n
            for Q in allQs[d]:
                if Q.solve_integer(n):
                    return Q
            raise ValueError("No form of discriminant %d represents %s" % (d,n))

        mat = self._mat
        qfmat = [[0 for i in range(ncurves)] for j in range(ncurves)]
        for i, E1 in enumerate(self.curves):
            for j, E2 in enumerate(self.curves):
                if j < i:
                    qfmat[i][j] = qfmat[j][i]
                    mat[i,j] = mat[j,i]
                elif i == j:
                    qfmat[i][j] = [1]
                    # mat[i,j] already 1
                else:
                    d = E1.cm_discriminant()
                    if d != E2.cm_discriminant():
                        qfmat[i][j] = [mat[i,j]]
                        # mat[i,j] already unique
                    else: # horizontal isogeny
                        q = find_quadratic_form(d,mat[i,j])
                        qfmat[i][j] = list(q)
                        mat[i,j] = q.small_prime_value()

        self._mat = mat
        self._qfmat = qfmat
        if verbose:
            print("new matrix = %s" % mat)
            print("matrix of forms = %s" % qfmat)

    def _compute_matrix(self):
        """
        Compute the matrix, assuming that the list of curves is computed.

        EXAMPLES::

            sage: isocls = EllipticCurve('1225h1').isogeny_class('database')
            sage: isocls._mat
            sage: isocls._compute_matrix(); isocls._mat
            [ 0 37]
            [37  0]
        """
        self._mat = self.E.isogeny_class(order=self.curves)._mat

    def _compute_isogenies(self):
        """
        EXAMPLES::

            sage: E = EllipticCurve('15a1')
            sage: isocls = E.isogeny_class()
            sage: maps = isocls.isogenies() # indirect doctest
            sage: f = maps[0][1]
            sage: f.domain() == isocls[0] and f.codomain() == isocls[1]
            True
        """
        recomputed = self.E.isogeny_class(order=self.curves)
        self._mat = recomputed._mat
        # The domains and codomains here will be equal, but not the same Python object.
        self._maps = recomputed._maps


class IsogenyClass_EC_Rational(IsogenyClass_EC_NumberField):
    r"""
    Isogeny classes for elliptic curves over `\QQ`.
    """
    def __init__(self, E, algorithm='sage', label=None, empty=False):
        r"""
        INPUT:

        - ``E`` -- an elliptic curve over `\QQ`

        - ``algorithm`` -- string (default: ``'sage'``); one of the
          following:

          - ``'sage'`` -- use sage's implementation to compute the curves,
            matrix and isogenies

          - ``'database'`` -- use the Cremona database (only works if the
            curve is in the database)

        - ``label`` -- string; the label of this isogeny class
          (e.g. '15a' or '37.b'), used in printing

        - ``empty`` -- don't compute the curves right now (used when reordering)

        EXAMPLES::

            sage: isocls = EllipticCurve('389a1').isogeny_class(); isocls
            Elliptic curve isogeny class 389a
            sage: E = EllipticCurve([0, 0, 0, 0, 1001]) # conductor 108216108
            sage: E.isogeny_class(order='database')
            Traceback (most recent call last):
            ...
            LookupError: Cremona database does not contain entry for
            Elliptic Curve defined by y^2 = x^3 + 1001 over Rational Field
            sage: TestSuite(isocls).run()
        """
        self._algorithm = algorithm
        IsogenyClass_EC.__init__(self, E, label=label, empty=empty)

    def copy(self):
        """
        Return a copy (mostly used in reordering).

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: C = E.isogeny_class()
            sage: C2 = C.copy()
            sage: C is C2
            False
            sage: C == C2
            True
        """
        ans = IsogenyClass_EC_Rational(self.E, self._algorithm, self._label, empty=True)
        # The following isn't needed internally, but it will keep
        # things from breaking if this is used for something other
        # than reordering.
        ans.curves = self.curves
        ans._mat = None
        ans._maps = None
        return ans

    def _compute(self):
        """
        Compute the list of curves, and possibly the matrix and
        prime-degree isogenies (depending on the algorithm selected).

        EXAMPLES::

            sage: isocls = EllipticCurve('48a1').isogeny_class('sage').copy()
            sage: isocls._mat
            sage: isocls._compute(); isocls._mat
            [0 2 2 2 0 0]
            [2 0 0 0 2 2]
            [2 0 0 0 0 0]
            [2 0 0 0 0 0]
            [0 2 0 0 0 0]
            [0 2 0 0 0 0]
        """
        algorithm = self._algorithm
        from sage.matrix.matrix_space import MatrixSpace
        self._maps = None
        if algorithm == "database":
            try:
                label = self.E.cremona_label(space=False)
            except RuntimeError:
                raise RuntimeError("unable to find %s in the database" % self.E)
            db = sage.databases.cremona.CremonaDatabase()
            curves = db.isogeny_class(label)
            if not curves:
                raise RuntimeError("unable to find %s in the database" % self.E)
            # All curves will have the same conductor and isogeny class,
            # and there are most 8 of them, so lexicographic sorting is okay.
            self.curves = tuple(sorted(curves,
                                       key=lambda E: E.cremona_label()))
            self._mat = None
        elif algorithm == "sage":
            curves = [self.E.minimal_model()]
            ijl_triples = []
            l_list = None
            i = 0
            while i < len(curves):
                E = curves[i]
                isogs = E.isogenies_prime_degree(l_list)
                for phi in isogs:
                    Edash = phi.codomain()
                    l = phi.degree()
                    # look to see if Edash is new.  Note that the
                    # curves returned by isogenies_prime_degree() are
                    # standard minimal models, so it suffices to check
                    # equality rather than isomorphism here.
                    try:
                        j = curves.index(Edash)
                    except ValueError:
                        j = len(curves)
                        curves.append(Edash)
                    ijl_triples.append((i,j,l,phi))
                if l_list is None:
                    l_list = list({ZZ(f.degree()) for f in isogs})
                i += 1
            self.curves = tuple(curves)
            ncurves = len(curves)
            self._mat = MatrixSpace(ZZ,ncurves)(0)
            self._maps = [[0]*ncurves for _ in range(ncurves)]
            for i,j,l,phi in ijl_triples:
                self._mat[i,j] = l
                self._maps[i][j] = phi
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)


def isogeny_degrees_cm(E, verbose=False):
    r"""
    Return a list of primes `\ell` sufficient to generate the
    isogeny class of `E`, where `E` has CM.

    INPUT:

    - ``E`` -- an elliptic curve defined over a number field

    OUTPUT:

    A finite list of primes `\ell` such that every curve isogenous to
    this curve can be obtained by a finite sequence of isogenies of
    degree one of the primes in the list.  This list is not
    necessarily minimal.

    ALGORITHM:

    For curves with CM by the order `O` of discriminant `d`, the
    Galois representation is always non-surjective and the curve will
    admit `\ell`-isogenies for infinitely many primes `\ell`, but
    there are only finitely many codomains `E'`.  The primes can be
    divided according to the discriminant `d'` of the CM order `O'`
    associated to `E`: either `O=O'`, or one contains the other with
    index `\ell`, since `\ell O\subset O'` and vice versa.

    Case (1): `O=O'`.  The degrees of all isogenies between `E` and
    `E'` are precisely the integers represented by one of the classes
    of binary quadratic forms `Q` of discriminant `d`.  Hence to
    obtain all possible isomorphism classes of codomain `E'`, we need
    only use one prime `\ell` represented by each such class `Q`.  It
    would in fact suffice to use primes represented by forms which
    generate the class group.  Here we simply omit the principal class
    and one from each pair of inverse classes, and include a prime
    represented by each of the remaining forms.

    Case (2): `[O':O]=\ell`: so `d=\ell^2d;`.  We include all prime
    divisors of `d`.

    Case (3): `[O:O']=\ell`: we may assume that `\ell` does not divide
    `d` as we have already included these, so `\ell` either splits or
    is inert in `O`; the class numbers satisfy `h(O')=(\ell\pm1)h(O)`
    accordingly.  We include all primes `\ell` such that `\ell\pm1`
    divides the degree `[K:\QQ]`.

    For curves with only potential CM we proceed as in the CM case,
    using `2[K:\QQ]` instead of `[K:\QQ]`.

    EXAMPLES:

    For curves with CM by a quadratic order of class number greater
    than `1`, we use the structure of the class group to only give one
    prime in each ideal class::

        sage: pol = PolynomialRing(QQ,'x')([1,-3,5,-5,5,-3,1])
        sage: L.<a> = NumberField(pol)
        sage: j = hilbert_class_polynomial(-23).roots(L, multiplicities=False)[0]
        sage: E = EllipticCurve(j=j)
        sage: from sage.schemes.elliptic_curves.isogeny_class import isogeny_degrees_cm
        sage: isogeny_degrees_cm(E, verbose=True)
        CM case, discriminant = -23
        initial primes: {2}
        upward primes: {}
        downward ramified primes: {}
        downward split primes: {2, 3}
        downward inert primes: {5}
        primes generating the class group: [2]
        Set of primes before filtering: {2, 3, 5}
        List of primes after filtering: [2, 3]
        [2, 3]

    TESTS:

    Check that :issue:`36780` is fixed::

        sage: L5.<r5> = NumberField(x^2-5)
        sage: E = EllipticCurve(L5,[0,-4325477943600 *r5-4195572876000])
        sage: from sage.schemes.elliptic_curves.isogeny_class import isogeny_degrees_cm
        sage: isogeny_degrees_cm(E)
        [3, 5]
    """
    if not E.has_cm():
        raise ValueError("possible_isogeny_degrees_cm(E) requires E to be an elliptic curve with CM")
    d = E.cm_discriminant()

    if verbose:
        print("CM case, discriminant = %s" % d)

    from sage.libs.pari.all import pari
    from sage.sets.set import Set
    from sage.arith.misc import kronecker as kronecker_symbol

    n = E.base_field().absolute_degree()
    if not E.has_rational_cm():
        n *= 2
    # For discriminants with extra units there's an extra factor in the class number formula:
    if d == -4:
        n *= 2
    if d == -3:
        n *= 3
    divs = n.divisors()

    data = pari(d).quadclassunit()
    # This has 4 components: the class number, class group
    # structure (ignored), class group generators (as quadratic
    # forms) and regulator (=1 since d<0, ignored).

    h = data[0].sage()

    # We must have 2*h dividing n, and will need the quotient so
    # see if the j-invariants of any proper sub-orders could lie
    # in the same field

    n_over_2h = n//(2*h)

    # Collect possible primes.  First put in 2, and also 3 for
    # discriminant -3 (special case because of units):

    L = Set([ZZ(2), ZZ(3)]) if d == -3 else Set([ZZ(2)])
    if verbose:
        print("initial primes: %s" % L)

    # Step 1: "vertical" primes l such that the isogenous curve
    # has CM by an order whose index is l or 1/l times the index
    # of the order O of discriminant d.  The latter case can only
    # happen when l^2 divides d.

    # (a) ramified primes

    ram_l = d.odd_part().prime_factors()

    # if the CM is not rational we include all ramified primes,
    # which is simpler than using the class group later:

    if not E.has_rational_cm():
        L1 = Set(ram_l)
        L += L1
        if verbose:
            print("ramified primes: %s" % L1)

    else:

        # "Upward" primes (index divided by l):

        L1 = Set([l for l in ram_l if d.valuation(l) > 1])
        L += L1
        if verbose:
            print("upward primes: %s" % L1)

        # "Downward" ramified primes; index multiplied by l, class
        # number multiplied by l, so l must divide n/2h:

        L1 = Set([l for l in ram_l if l.divides(n_over_2h)])
        L += L1
        if verbose:
            print("downward ramified primes: %s" % L1)

    # (b) Downward split primes; the suborder has class number (l-1)*h, so
    # l-1 must divide n/2h:

    L1 = Set([lm1+1 for lm1 in divs
              if (lm1+1).is_prime() and kronecker_symbol(d,lm1+1) == +1])
    L += L1
    if verbose:
        print("downward split primes: %s" % L1)

    # (c) Downward inert primes; the suborder has class number (l+1)*h, so
    # l+1 must divide n/2h:

    L1 = Set([lp1-1 for lp1 in divs
              if (lp1-1).is_prime() and kronecker_symbol(d,lp1-1) == -1])
    L += L1
    if verbose:
        print("downward inert primes: %s" % L1)

    # Horizontal primes (rational CM only): same order, degrees are
    # all integers represented by some binary quadratic form of
    # discriminant d, so we find a prime represented by each form.

    if E.has_rational_cm():
        from sage.quadratic_forms.binary_qf import BinaryQF
        Qs = [BinaryQF(list(q)) for q in data[2]]

        L1 = [Q.small_prime_value() for Q in Qs]
        if verbose:
            print("primes generating the class group: %s" % L1)
        L += Set(L1)

    # Return sorted list

    if verbose:
        print("Set of primes before filtering: %s" % L)

    # This filter will quickly eliminate most false entries in the set
    from .gal_reps_number_field import Frobenius_filter
    L = Frobenius_filter(E, sorted(L))
    if verbose:
        print("List of primes after filtering: %s" % L)
    return L


def possible_isogeny_degrees(E, algorithm='Billerey', max_l=None,
                             num_l=None, exact=True, verbose=False):
    r"""
    Return a list of primes `\ell` sufficient to generate the
    isogeny class of `E`.

    INPUT:

    - ``E`` -- an elliptic curve defined over a number field

    - ``algorithm`` -- string (default: ``'Billerey'``); algorithm to be
      used for non-CM curves: either ``'Billerey'``, ``'Larson'``, or
      ``'heuristic'``.  Only relevant for non-CM curves and base fields
      other than `\QQ`.

    - ``max_l`` -- integer or ``None``; only relevant for non-CM curves
      and algorithms ``'Billerey'`` and ``'heuristic'``.  Controls the maximum
      prime used in either algorithm.  If ``None``, use the default
      for that algorithm.

    - ``num_l`` -- integer or ``None``; only relevant for non-CM curves
      and algorithm ``'Billerey'``.  Controls the maximum number of primes
      used in the algorithm.  If ``None``, use the default for that
      algorithm.

    - ``exact`` -- boolean (default: ``True``); if ``True``, perform an
      additional check that the primes returned are all reducible.  If
      ``False``, skip this step, in which case some of the primes
      returned may be irreducible.

    OUTPUT:

    A finite list of primes `\ell` such that every curve isogenous to
    this curve can be obtained by a finite sequence of isogenies of
    degree one of the primes in the list.

    ALGORITHM:

    For curves without CM, the set may be taken to be the finite set
    of primes at which the Galois representation is not surjective,
    since the existence of an `\ell`-isogeny is equivalent to the
    image of the mod-`\ell` Galois representation being contained in a
    Borel subgroup.  Two rigorous algorithms have been implemented to
    determine this set, due to Larson and Billeray respectively.  We
    also provide a non-rigorous 'heuristic' algorithm which only tests
    reducible primes up to a bound depending on the degree of the
    base field.

    For curves with CM see the documentation for :meth:`isogeny_degrees_cm()`.

    EXAMPLES:

    For curves without CM we determine the primes at which the mod `p`
    Galois representation is reducible, i.e. contained in a Borel
    subgroup::

        sage: from sage.schemes.elliptic_curves.isogeny_class import possible_isogeny_degrees
        sage: E = EllipticCurve('11a1')
        sage: possible_isogeny_degrees(E)
        [5]
        sage: possible_isogeny_degrees(E, algorithm='Larson')
        [5]
        sage: possible_isogeny_degrees(E, algorithm='Billerey')
        [5]
        sage: possible_isogeny_degrees(E, algorithm='heuristic')
        [5]

    We check that in this case `E` really does have rational
    `5`-isogenies::

        sage: [phi.degree() for phi in E.isogenies_prime_degree()]
        [5, 5]

    Over an extension field::

        sage: E3 = E.change_ring(CyclotomicField(3))
        sage: possible_isogeny_degrees(E3)                                              # long time (5s)
        [5]
        sage: [phi.degree() for phi in E3.isogenies_prime_degree()]
        [5, 5]

    A higher degree example (LMFDB curve 5.5.170701.1-4.1-b1)::

        sage: K.<a> = NumberField(x^5 - x^4 - 6*x^3 + 4*x + 1)
        sage: E = EllipticCurve(K, [a^3 - a^2 - 5*a + 1, a^4 - a^3 - 5*a^2 - a + 1,
        ....:                       -a^4 + 2*a^3 + 5*a^2 - 5*a - 3, a^4 - a^3 - 5*a^2 - a,
        ....:                       -3*a^4 + 4*a^3 + 17*a^2 - 6*a - 12])
        sage: possible_isogeny_degrees(E, algorithm='heuristic')
        [2]
        sage: possible_isogeny_degrees(E, algorithm='Billerey')
        [2]
        sage: possible_isogeny_degrees(E, algorithm='Larson')
        [2]

    LMFDB curve 4.4.8112.1-108.1-a5::

        sage: K.<a> = NumberField(x^4 - 5*x^2 + 3)
        sage: E = EllipticCurve(K, [a^2 - 2, -a^2 + 3, a^2 - 2, -50*a^2 + 35, 95*a^2 - 67])
        sage: possible_isogeny_degrees(E, exact=False, algorithm='Billerey')            # long time (6.5s)
        [2, 5]
        sage: possible_isogeny_degrees(E, exact=False, algorithm='Larson')
        [2, 5]
        sage: possible_isogeny_degrees(E, exact=False, algorithm='heuristic')
        [2, 5]
        sage: possible_isogeny_degrees(E)
        [2, 5]

    This function only returns the primes which are isogeny degrees::

        sage: Set(E.isogeny_class().matrix().list())                                    # long time (7s)
        {1, 2, 4, 5, 20, 10}

    For curves with CM by a quadratic order of class number greater
    than `1`, we use the structure of the class group to only give one
    prime in each ideal class::

        sage: pol = PolynomialRing(QQ,'x')([1,-3,5,-5,5,-3,1])
        sage: L.<a> = NumberField(pol)
        sage: j = hilbert_class_polynomial(-23).roots(L, multiplicities=False)[0]
        sage: E = EllipticCurve(j=j)
        sage: from sage.schemes.elliptic_curves.isogeny_class import possible_isogeny_degrees
        sage: possible_isogeny_degrees(E, verbose=True)
        CM case, discriminant = -23
        initial primes: {2}
        upward primes: {}
        downward ramified primes: {}
        downward split primes: {2, 3}
        downward inert primes: {5}
        primes generating the class group: [2]
        Set of primes before filtering: {2, 3, 5}
        List of primes after filtering: [2, 3]
        [2, 3]
    """
    if E.has_cm():
        return isogeny_degrees_cm(E, verbose)

    if E.base_field() == QQ:
        from sage.schemes.elliptic_curves.gal_reps_number_field import reducible_primes_naive
        return reducible_primes_naive(E, max_l=37, verbose=verbose)

    #  Non-CM case

    # NB The following functions first computes a finite set
    # containing the reducible primes, then checks that each is
    # reducible by computing l-isogenies.  This appears circular
    # but the computated l-isogenies for a fixed prime l is
    # cached.

    if verbose:
        print("Non-CM case, using {} algorithm".format(algorithm))

    # First we obtain a finite set of primes containing the reducible
    # ones Each of these algorithms includes application of the
    # "Frobenius filter" eliminating any ell for which there exists a
    # prime P of good reduction such that the Frobenius polynomial at
    # P does not factor modulo ell.

    if algorithm == 'Larson':
        L = E.galois_representation().isogeny_bound()

    elif algorithm == 'Billerey':
        from sage.schemes.elliptic_curves.gal_reps_number_field import reducible_primes_Billerey
        L = reducible_primes_Billerey(E, num_l=num_l, max_l=max_l, verbose=verbose)

    elif algorithm == 'heuristic':
        from sage.schemes.elliptic_curves.gal_reps_number_field import reducible_primes_naive
        L = reducible_primes_naive(E, max_l=max_l, num_P=num_l, verbose=verbose)

    else:
        raise ValueError("algorithm for possible_isogeny_degrees must be one of 'Larson', 'Billerey', 'heuristic'")

    # The set L may contain irreducible primes.  We optionally test
    # each one to see if it is actually reducible, by computing ell-isogenies:

    if exact:
        L = [l for l in L if E.isogenies_prime_degree(l, minimal_models=False)]

    return L
