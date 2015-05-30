"""
Path Algebras
"""

#*****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement

class PathAlgebra(CombinatorialFreeModule):
    r"""
    Create the path algebra of a :class:`quiver <DiGraph>` over a given field.

    Given a quiver `Q` and a field `k`, the path algebra `kQ` is defined as
    follows.  As a vector space it has basis the set of all paths in `Q`.
    Multiplication is defined on this basis and extended bilinearly.  If `p`
    is a path with terminal vertex `t` and `q` is a path with initial vertex
    `i` then the product `p*q` is defined to be the composition of the
    paths `p` and `q` if `t = i` and `0` otherwise.

    INPUT:

    - ``k`` -- field (or commutative ring), the base field of the path algebra

    - ``P`` -- the path semigroup of a quiver `Q`

    OUTPUT:

    - the path algebra `kP`

    EXAMPLES::

        sage: P = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
        sage: A = P.algebra(GF(7))
        sage: A
        Path algebra of Multi-digraph on 3 vertices over Finite Field of size 7
        sage: A.variable_names()
        ('e_1', 'e_2', 'e_3', 'a', 'b')

    Note that path algebras are uniquely defined by their quiver and field::

        sage: A is P.algebra(GF(7))
        True
        sage: A is P.algebra(RR)
        False
        sage: A is DiGraph({1:{2:['a']}}).path_semigroup().algebra(GF(7))
        False

    The path algebra of an acyclic quiver has a finite basis::

        sage: A.dimension()
        6
        sage: list(A.basis())
        [e_1, e_2, e_3, a, b, a*b]

    The path algebra can create elements from paths or from elements of the
    base ring::

        sage: A(5)
        5*e_1 + 5*e_2 + 5*e_3
        sage: S = A.semigroup()
        sage: S
        Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices
        sage: p = S((1, 2, 'a'))
        sage: r = S((2, 3, 'b'))
        sage: e2 = S((2, 2))
        sage: x = A(p) + A(e2)
        sage: x
        e_2 + a
        sage: y = A(p) + A(r)
        sage: y
        a + b

    Path algebras are graded algebras.  The grading is given by assigning
    to each basis element the length of the path corresponding to that
    basis element::

        sage: x.is_homogeneous()
        False
        sage: x.degree()
        Traceback (most recent call last):
        ...
        ValueError: Element is not homogeneous.
        sage: y.is_homogeneous()
        True
        sage: y.degree()
        1
        sage: A[1]
        Free module spanned by [a, b] over Finite Field of size 7
        sage: A[2]
        Free module spanned by [a*b] over Finite Field of size 7

    TESTS::

        sage: TestSuite(A).run()
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, k, P):
        """
        Creates a :class:`PathAlgebra` object.  Type ``PathAlgebra?`` for
        more information.

        INPUT:

        - ``k`` -- a commutative ring
        - ``P`` -- the partial semigroup formed by the paths of a quiver

        TESTS::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}}).path_semigroup()
            sage: P.algebra(GF(5))
            Path algebra of Multi-digraph on 4 vertices over Finite Field of size 5
        """
        # The following hidden methods are relevant:
        #
        # - _base
        #       The base ring of the path algebra.
        # - _basis_keys
        #       Finite enumerated set containing the QuiverPaths that form the
        #       basis.
        # - _quiver
        #       The quiver of the path algebra
        # - _semigroup
        #       Shortcut for _quiver.semigroup()

        from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
        self._quiver = P.quiver()
        self._semigroup = P
        super(PathAlgebra, self).__init__(k, self._semigroup,
                                             prefix='',
                                             element_class=self.Element,
                                             category=GradedAlgebrasWithBasis(k),
                                             bracket=False)
        self._assign_names(self._semigroup.variable_names())

    @cached_method
    def gens(self):
        """
        Return the generators of this algebra (idempotents and arrows).

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}}).path_semigroup()
            sage: A = P.algebra(GF(5))
            sage: A.variable_names()
            ('e_1', 'e_2', 'e_3', 'e_4', 'a', 'b', 'c')
            sage: A.gens()
            (e_1, e_2, e_3, e_4, a, b, c)
        """
        return tuple(self._from_dict( {index: self.base_ring().one()},
                                      remove_zeros=False )
                     for index in self._semigroup.gens())

    @cached_method
    def arrows(self):
        """
        Return the arrows of this algebra (corresponding to edges of the
        underlying quiver).

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}}).path_semigroup()
            sage: A = P.algebra(GF(5))
            sage: A.arrows()
            (a, b, c)
        """
        return tuple(self._from_dict( {index: self.base_ring().one()},
                                      remove_zeros=False )
                     for index in self._semigroup.arrows())

    @cached_method
    def idempotents(self):
        """
        Return the idempotents of this algebra (corresponding to vertices
        of the underlying quiver).

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}}).path_semigroup()
            sage: A = P.algebra(GF(5))
            sage: A.idempotents()
            (e_1, e_2, e_3, e_4)
        """
        return tuple(self._from_dict( {index: self.base_ring().one()},
                                      remove_zeros=False )
                     for index in self._semigroup.idempotents())

    def gen(self, i):
        """
        Return the `i`-th generator of this algebra.

        This is an idempotent (corresponding to a trivial path at a
        vertex) if `i < n` (where `n` is the number of vertices of the
        quiver), and a single-edge path otherwise.

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}}).path_semigroup()
            sage: A = P.algebra(GF(5))
            sage: A.gens()
            (e_1, e_2, e_3, e_4, a, b, c)
            sage: A.gen(2)
            e_3
            sage: A.gen(5)
            b
        """
        return self._from_dict( {self._semigroup.gen(i): self.base_ring().one()},
                                remove_zeros = False )

    def ngens(self):
        """
        Number of generators of this algebra.

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}}).path_semigroup()
            sage: A = P.algebra(GF(5))
            sage: A.ngens()
            7

        """
        return self._semigroup.ngens()

    def _element_constructor_(self, x):
        """
        Attempt to construct an element of ``self`` from ``x``.

        TESTS::

            sage: A = DiGraph({1:{2:['a']}}).path_semigroup().algebra(QQ)
            sage: B = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup().algebra(QQ)
            sage: x = A((1, 2, 'a')) + 1 # indirect doctest
            sage: x
            e_1 + e_2 + a
            sage: B(x) # indirect doctest
            e_1 + e_2 + a
            sage: A(1) # indirect doctest
            e_1 + e_2
        """
        from sage.quivers.paths import QuiverPath
        # If it's an element of another path algebra, do a linear combination
        # of the basis
        if isinstance(x, CombinatorialFreeModuleElement) and isinstance(x.parent(), PathAlgebra):
            coeffs = x.monomial_coefficients()
            result = self.zero()
            for key in coeffs:
                result += coeffs[key]*self.monomial(key)
            return result

        # If it's a QuiverPath return the associated basis element
        if isinstance(x, QuiverPath):
            return self.monomial(x)

        # If it's a tuple or a list try and create a QuiverPath from it and
        # then return the associated basis element
        if isinstance(x, (tuple, list)):
            return self.monomial(self._semigroup(x))

        # Otherwise let CombinatorialFreeModule try
        return super(PathAlgebra, self)._element_constructor_(x)

    def _coerce_map_from_(self, other):
        """
        Return ``True`` if there is a coercion from ``other`` to ``self``.

        The algebras that coerce into a path algebra are rings `k` or path
        algebras `kQ` such that `k` has a coercion into the base ring of
        ``self`` and `Q` is a subquiver of the quiver of ``self``.

        In particular, the path semigroup of a subquiver coerces into the
        algebra.

        TESTS::

            sage: P1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P2 = DiGraph({1:{2:['a','b']}}).path_semigroup()
            sage: A1 = P1.algebra(GF(3))
            sage: A2 = P2.algebra(GF(3))
            sage: A1.coerce_map_from(A2) # indirect doctest
            sage: A2.coerce_map_from(A1) # indirect doctest
            Conversion map:
              From: Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3
              To:   Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3
            sage: A1.coerce_map_from(ZZ) # indirect doctest
            Composite map:
                  From: Integer Ring
                  To:   Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3
                  Defn:   Natural morphism:
                          From: Integer Ring
                          To:   Finite Field of size 3
                        then
                          Generic morphism:
                          From: Finite Field of size 3
                          To:   Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3
            sage: A1.coerce_map_from(QQ) # indirect doctest
            sage: A1.coerce_map_from(ZZ)
            Composite map:
              From: Integer Ring
              To:   Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Finite Field of size 3
                    then
                      Generic morphism:
                      From: Finite Field of size 3
                      To:   Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3

        ::

            sage: A2.coerce_map_from(P1)
            Conversion map:
              From: Partial semigroup formed by the directed paths of Multi-digraph on 2 vertices
              To:   Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3
            sage: a = P1(P1.arrows()[0]); a
            a
            sage: A2.one() * a == a     # indirect doctest
            True

        """

        if isinstance(other, PathAlgebra) and self._base.has_coerce_map_from(other._base):
            OQ = other._quiver
            SQ = self._quiver
            SQE = SQ.edges()
            if all(v in SQ for v in OQ.vertices()) and all(e in SQE for e in OQ.edges()):
                return True
        if self._semigroup.has_coerce_map_from(other):
            return True
        return self._base.has_coerce_map_from(other)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: P.algebra(RR) # indirect doctest
            Path algebra of Multi-digraph on 3 vertices over Real Field with 53 bits of precision
        """
        return "Path algebra of {0} over {1}".format(self._quiver, self._base)

    ###########################################################################
    #                                                                         #
    # CATEGORY METHODS                                                        #
    #    These functions are used by the category to implement the algebra    #
    #    structure.                                                           #
    #                                                                         #
    ###########################################################################

    def _monomial(self, index):
        """
        This method makes sure that the invalid path evaluates as zero.

        TESTS::

            sage: P = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: P.algebra(ZZ)(P((1,2,'a'))*P((2,3,'b')))
            a*b
            sage: P.algebra(ZZ)(P((2,3,'b')))*P.algebra(ZZ)(P((1,2,'a')))  # indirect doctest
            0

        The following was an issue during work at :trac:`12630`::

            sage: P1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P2 = DiGraph({1:{2:['a','b']}}).path_semigroup()
            sage: A1 = P1.algebra(GF(3))
            sage: A2 = P2.algebra(GF(3))
            sage: b = P2.arrows()[1]; b
            b
            sage: A1(b)
            Traceback (most recent call last):
            ...
            ValueError: Cannot interpret b as element of Partial semigroup
            formed by the directed paths of Multi-digraph on 2 vertices
        """
        if index is not None:
            return self._from_dict( {self._semigroup(index): self.base_ring().one()},
                                    remove_zeros=False )
        return self.zero()

    def product_on_basis(self, p1, p2):
        """
        Return the product ``p1*p2`` in the path algebra.

        INPUT:

        - ``p1``, ``p2`` -- QuiverPaths

        OUTPUT:

        - :class:`~sage.combinat.free_module.CombinatorialFreeModuleElement`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p1 = Q((1, 2, 'a'))
            sage: p2 = Q([(2, 3, 'b'), (3, 4, 'c')])
            sage: A = Q.algebra(QQ)
            sage: A.product_on_basis(p1, p2)
            a*b*c
            sage: A.product_on_basis(p2, p1)
            0
        """
        PSG = self._semigroup
        p = PSG(p1)*PSG(p2)
        if p is not None:
            return self.basis()[p]
        else:
            return self.zero()

    def degree_on_basis(self, p):
        """
        Return the degree of the monomial specified by the path ``p``.

        EXAMPLES::

            sage: A = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup().algebra(QQ)
            sage: A.degree_on_basis((1, 1))
            0
            sage: A.degree_on_basis((1, 2, 'a'))
            1
            sage: A.degree_on_basis([(1, 2, 'a'), (2, 3, 'b')])
            2
        """
        return len(self._semigroup(p))

    def one(self):
        """
        Return the multiplicative identity element.

        The multiplicative identity of a path algebra is the sum of the basis
        elements corresponding to the trivial paths at each vertex.

        EXAMPLES::

            sage: A = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup().algebra(QQ)
            sage: A.one()
            e_1 + e_2 + e_3
        """
        return self.sum_of_monomials([self._semigroup([(v, v)], check=False)
                                      for v in self._quiver])

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data and subspaces of the path algebra.       #
    #                                                                         #
    ###########################################################################

    def quiver(self):
        """
        Return the quiver from which the algebra ``self`` was formed.

        OUTPUT:

        - :class:`DiGraph`, the quiver of the algebra

        EXAMPLES:

            sage: P = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: A = P.algebra(GF(3))
            sage: A.quiver() is P.quiver()
            True
        """
        return self._quiver

    def semigroup(self):
        """
        Return the (partial) semigroup from which the algebra ``self`` was
        constructed.

        .. NOTE::

            The partial semigroup is formed by the paths of a quiver,
            multiplied by concatenation. If the quiver has more than a single
            vertex, then multiplication in the path semigroup is not always
            defined.

        OUTPUT:

        - the path semigroup from which ``self`` was formed (a partial
          semigroup)

        EXAMPLES:

            sage: P = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: A = P.algebra(GF(3))
            sage: A.semigroup() is P
            True
        """
        return self._semigroup

    def homogeneous_component(self, n):
        """
        Return the `n`-th homogeneous piece of the path algebra.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        - :class:`CombinatorialFreeModule`, module spanned by the paths
          of length `n` in the quiver

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}}).path_semigroup()
            sage: A = P.algebra(GF(7))
            sage: A.homogeneous_component(2)
            Free module spanned by [a*c, b*d] over Finite Field of size 7

            sage: D = DiGraph({1: {2: 'a'}, 2: {3: 'b'}, 3: {1: 'c'}})
            sage: P = D.path_semigroup()
            sage: A = P.algebra(ZZ)
            sage: A.homogeneous_component(3)
            Free module spanned by [a*b*c, b*c*a, c*a*b] over Integer Ring
        """
        basis = []
        for v in self._semigroup._quiver:
            basis.extend(self._semigroup.iter_paths_by_length_and_startpoint(n, v))
        M = CombinatorialFreeModule(self._base, basis, prefix='', bracket=False)
        M._name = "Free module spanned by {0}".format(basis)
        return M

    __getitem__ = homogeneous_component

    def homogeneous_components(self):
        r"""
        Return the non-zero homogeneous components of ``self``.

        EXAMPLES::

            sage: Q = DiGraph([[1,2,'a'],[2,3,'b'],[3,4,'c']])
            sage: PQ = Q.path_semigroup()
            sage: A = PQ.algebra(GF(7))
            sage: A.homogeneous_components()
            [Free module spanned by [e_1, e_2, e_3, e_4] over Finite Field of size 7,
             Free module spanned by [a, b, c] over Finite Field of size 7,
             Free module spanned by [a*b, b*c] over Finite Field of size 7,
             Free module spanned by [a*b*c] over Finite Field of size 7]

        .. WARNING::

             Backward incompatible change: since :trac:`12630` and
             until :trac:`8678`, this feature was implemented under
             the syntax ``list(A)`` by means of ``A.__iter__``. This
             was incorrect since ``A.__iter__``, when defined for a
             parent, should iterate through the elements of `A`.
        """
        result = []
        i = 0
        while True:
            c = self.homogeneous_component(i)
            if not c.dimension():
                break
            result.append(c)
            i += 1
        return result

    ###########################################################################
    #                                                                         #
    # ELEMENT CLASS                                                           #
    #    The class of elements of the path algebra.                           #
    #                                                                         #
    ###########################################################################

    class Element(CombinatorialFreeModuleElement):
        def is_homogeneous(self):
            """
            Return ``True`` if and only if this element is homogeneous.

            EXAMPLES::

                sage: A = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup().algebra(QQ)
                sage: (A((1, 2, 'a')) + A((1, 2, 'b'))).is_homogeneous()
                True
                sage: (A((1, 1)) + A((1, 2, 'a'))).is_homogeneous()
                False
            """
            # Get the support, the zero element is homogeneous
            paths = self.support()
            if not paths:
                return True

            # Compare the rest of the paths, they must be the same length
            for p in paths[1:]:
                if len(p) != len(paths[0]):
                    return False

            return True

        def degree(self):
            """
            The degree of ``self``, if ``self`` is homogeneous.

            EXAMPLES::

                sage: A = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup().algebra(QQ)
                sage: A((1, 1)).degree()
                0
                sage: (A((1, 2, 'a')) + A((1, 2, 'b'))).degree()
                1

            An error is raised if the element is not homogeneous::

                sage: (A((1, 1)) + A((1, 2, 'a'))).degree()
                Traceback (most recent call last):
                ...
                ValueError: Element is not homogeneous.
            """
            # Deal with zero
            paths = self.support()
            if not paths:
                raise ValueError("The zero element does not have a well-defined degree.")

            # Check that the element is homogeneous
            for p in paths[1:]:
                if len(p) != len(paths[0]):
                    raise ValueError("Element is not homogeneous.")

            return len(paths[0])

