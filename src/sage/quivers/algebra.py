"""
Path Algebras
"""

# ****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
#                2014 Simon King <simon.king@uni-jena.de>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule
from .algebra_elements import PathAlgebraElement


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

    - ``order`` -- string; one of ``'negdegrevlex'`` (default),
      ``'degrevlex'``, ``'negdeglex'`` or ``'deglex'``, defining the monomial
      order to be used

    OUTPUT: the path algebra `kP` with the given monomial order

    .. NOTE::

        Monomial orders that are not degree orders are not supported.

    EXAMPLES::

        sage: P = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
        sage: A = P.algebra(GF(7))
        sage: A
        Path algebra of Multi-digraph on 3 vertices over Finite Field of size 7
        sage: A.variable_names()
        ('e_1', 'e_2', 'e_3', 'a', 'b')

    Note that path algebras are uniquely defined by their quiver, field and
    monomial order::

        sage: A is P.algebra(GF(7))
        True
        sage: A is P.algebra(GF(7), order='degrevlex')
        False
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
        sage: p = S([(1, 2, 'a')])
        sage: r = S([(2, 3, 'b')])
        sage: e2 = S([(2, 2)])
        sage: x = A(p) + A(e2)
        sage: x
        a + e_2
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
        ValueError: element is not homogeneous
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

    Element = PathAlgebraElement

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, k, P, order='negdegrevlex'):
        """
        Create a :class:`PathAlgebra` object.

        Type ``PathAlgebra?`` for more information.

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
        self._ordstr = order
        super().__init__(k, self._semigroup,
                         prefix='',
                         # element_class=self.Element,
                         category=GradedAlgebrasWithBasis(k),
                         bracket=False)
        self._assign_names(self._semigroup.variable_names())

    def order_string(self) -> str:
        """
        Return the string that defines the monomial order of this algebra.

        EXAMPLES::

            sage: P1 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'))
            sage: P2 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='degrevlex')
            sage: P3 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='negdeglex')
            sage: P4 = DiGraph({1:{1:['x','y','z']}}).path_semigroup().algebra(GF(25,'t'), order='deglex')
            sage: P1.order_string()
            'negdegrevlex'
            sage: P2.order_string()
            'degrevlex'
            sage: P3.order_string()
            'negdeglex'
            sage: P4.order_string()
            'deglex'
        """
        return self._ordstr

    @cached_method
    def gens(self) -> tuple:
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
        return tuple(self.gen(i) for i in range(self.ngens()))

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
        return tuple(self._from_dict({index: self.base_ring().one()},
                                     remove_zeros=False)
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
        return tuple(self._from_dict({index: self.base_ring().one()},
                                     remove_zeros=False)
                     for index in self._semigroup.idempotents())

    @cached_method
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
        return self._from_dict({self._semigroup.gen(i): self.base_ring().one()},
                               remove_zeros=False)

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

            sage: A = DiGraph({2:{3:['b']}}).path_semigroup().algebra(ZZ)
            sage: B = DiGraph({0:{1:['a']}, 1:{2:['c']}, 2:{3:['b']}}).path_semigroup().algebra(GF(5))
            sage: x = A('b') + 1 # indirect doctest
            sage: x
            e_2 + b + e_3
            sage: B(x) # indirect doctest
            e_2 + b + e_3
            sage: A(1) # indirect doctest
            e_2 + e_3
            sage: B(2) # indirect doctest
            2*e_0 + 2*e_1 + 2*e_2 + 2*e_3
            sage: B([(0,1,'a'),(1,2,'c')])  # indirect doctest
            a*c
        """
        from sage.quivers.paths import QuiverPath
        # If it's an element of another path algebra, do a linear combination
        # of the basis
        if isinstance(x, PathAlgebraElement) and isinstance(x.parent(), PathAlgebra):
            result = {}
            coeffs = x.monomial_coefficients()
            for key in coeffs:
                result[self._semigroup(key)] = coeffs[key]
            return self.element_class(self, result)

        # If it's a QuiverPath return the associated basis element
        if isinstance(x, QuiverPath):
            return self.element_class(self, {x: self.base_ring().one()})

        # If it's a scalar, return a multiple of one:
        if x in self.base_ring():
            return self.one() * x

        # If it's a tuple or a list, try and create a QuiverPath from it and
        # then return the associated basis element
        if isinstance(x, (tuple, list, str)):
            return self.element_class(self, {self._semigroup(x): self.base_ring().one()})

        if isinstance(x, dict):
            return self.element_class(self, x)

        # Otherwise let CombinatorialFreeModule try
        return super()._element_constructor_(x)

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
            Coercion map:
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
            Coercion map:
              From: Partial semigroup formed by the directed paths of Multi-digraph on 2 vertices
              To:   Path algebra of Multi-digraph on 2 vertices over Finite Field of size 3
            sage: a = P1(P1.arrows()[0]); a
            a
            sage: A2.one() * a == a     # indirect doctest
            True

        ::

            sage: A = DiGraph({2:{3:['b']}}).path_semigroup().algebra(ZZ)
            sage: B = DiGraph({0:{1:['a']}, 1:{2:['c']}, 2:{3:['b']}}).path_semigroup().algebra(GF(5))
            sage: x = A('b') + 1 # indirect doctest
            sage: x
            e_2 + b + e_3
            sage: B(2)
            2*e_0 + 2*e_1 + 2*e_2 + 2*e_3
            sage: B(2)*x*B(3)  # indirect doctest
            e_2 + b + e_3
        """
        if isinstance(other, PathAlgebra) and self._base.has_coerce_map_from(other._base):
            OQ = other._quiver
            SQ = self._quiver
            SQE = self._semigroup._sorted_edges
            if all(v in SQ for v in OQ.vertex_iterator()) and all(e in SQE for e in other._semigroup._sorted_edges):
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

    # String representation of a monomial
    def _repr_monomial(self, data):
        """
        String representation of a monomial.

        INPUT:

        - ``data`` -- list providing the indices of the path algebra
          generators occurring in the monomial

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: X         # indirect doctest
            5*e_0 + a + 2*b + 3*c + 3*e_2
        """
        # m is [list, pos, mid], where the list gives the nb of arrows, pos
        # gives the component in the module, and mid gives the length of the
        # left factor in a two-sided module.
        arrows = self.variable_names()
        return '*'.join(arrows[n] for n in data)

    def _latex_monomial(self, data):
        """
        Latex string representation of a monomial.

        INPUT:

        - ``data`` -- list providing the indices of the path algebra
          generators occurring in the monomial

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: latex(X)  # indirect doctest
            5 e_0 + a + 2 b + 3 c + 3 e_2
        """
        arrows = self.variable_names()
        return '\\cdot '.join(arrows[n] for n in data)

    @cached_method
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
        one = self.base_ring().one()
        D = {index: one for index in self._semigroup.idempotents()}
        return self._from_dict(D)

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data and subspaces of the path algebra.       #
    #                                                                         #
    ###########################################################################

    def quiver(self):
        """
        Return the quiver from which the algebra ``self`` was formed.

        OUTPUT: :class:`DiGraph`; the quiver of the algebra

        EXAMPLES::

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

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: A = P.algebra(GF(3))
            sage: A.semigroup() is P
            True
        """
        return self._semigroup

    def degree_on_basis(self, x):
        """
        Return ``x.degree()``.

        This function is here to make some methods work that are inherited
        from :class:`~sage.combinat.free_module.CombinatorialFreeModule`.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: X = a+2*b+3*c*e-a*d+5*e_0+3*e_2
            sage: X
            5*e_0 + a - a*d + 2*b + 3*e_2
            sage: X.homogeneous_component(0)   # indirect doctest
            5*e_0 + 3*e_2
            sage: X.homogeneous_component(1)
            a + 2*b
            sage: X.homogeneous_component(2)
            -a*d
            sage: X.homogeneous_component(3)
            0
        """
        return x.degree()

    def sum(self, iter_of_elements):
        """
        Return the sum of all elements in ``iter_of_elements``.

        INPUT:

        - ``iter_of_elements`` -- iterator of elements of ``self``

        .. NOTE::

            It overrides a method inherited from
            :class:`~sage.combinat.free_module.CombinatorialFreeModule`,
            which relies on a private attribute of elements---an
            implementation detail that is simply not available for
            :class:`~sage.quivers.algebra_elements.PathAlgebraElement`.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: A.sum((a, 2*b, 3*c*e, -a*d, 5*e_0, 3*e_2))
            5*e_0 + a - a*d + 2*b + 3*e_2
        """
        return sum(iter_of_elements, self.zero())

    def linear_combination(self, iter_of_elements_coeff, factor_on_left=True):
        r"""
        Return the linear combination `\lambda_1 v_1 + \cdots +
        \lambda_k v_k` (resp.  the linear combination `v_1 \lambda_1 +
        \cdots + v_k \lambda_k`) where ``iter_of_elements_coeff`` iterates
        through the sequence `((v_1, \lambda_1), ..., (v_k, \lambda_k))`.

        INPUT:

        - ``iter_of_elements_coeff`` -- iterator of pairs ``(element, coeff)``
          with ``element`` in ``self`` and ``coeff`` in ``self.base_ring()``

        - ``factor_on_left`` -- (optional) if ``True``, the coefficients are
          multiplied from the left if ``False``, the coefficients are
          multiplied from the right

        .. NOTE::

            It overrides a method inherited from
            :class:`~sage.combinat.free_module.CombinatorialFreeModule`,
            which relies on a private attribute of elements---an
            implementation detail that is simply not available for
            :class:`~sage.quivers.algebra_elements.PathAlgebraElement`.

        EXAMPLES::

            sage: A = DiGraph({0: {1: ['a'], 2: ['b']},
            ....:              1: {0: ['c'], 1: ['d']},
            ....:              2: {0: ['e'], 2: ['f']}}).path_semigroup().algebra(ZZ)
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: A.linear_combination([(a, 1), (b, 2), (c*e, 3),
            ....:                       (a*d, -1), (e_0, 5), (e_2, 3)])
            5*e_0 + a - a*d + 2*b + 3*e_2
        """
        if factor_on_left:
            return self.sum(coeff * element
                            for element, coeff in iter_of_elements_coeff)
        else:
            return self.sum(element * coeff
                            for element, coeff in iter_of_elements_coeff)

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
        Return the nonzero homogeneous components of ``self``.

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

             Backward incompatible change: since :issue:`12630` and
             until :issue:`8678`, this feature was implemented under
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
