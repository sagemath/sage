"""
Quiver Homspace
"""

# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.homset import Homset
from sage.quivers.morphism import QuiverRepHom
from sage.misc.cachefunc import cached_method


class QuiverHomSpace(Homset):
    r"""
    A homomorphism of quiver representations (of one and the same quiver)
    is given by specifying, for each vertex of the quiver, a homomorphism
    of the spaces assigned to this vertex such that these homomorphisms
    commute with the edge maps.  This class handles the set of all
    such maps, `Hom_Q(M, N)`.

    INPUT:

    - ``domain`` -- the domain of the homomorphism space

    - ``codomain`` -- the codomain of the homomorphism space

    OUTPUT:

    - :class:`QuiverHomSpace`, the homomorphism space
      ``Hom_Q(domain, codomain)``

    .. NOTE::

        The quivers of the domain and codomain must be equal or a
        ``ValueError`` is raised.

    EXAMPLES::

        sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
        sage: H = Q.S(QQ, 2).Hom(Q.P(QQ, 1))
        sage: H.dimension()
        2
        sage: H.gens()
        [Homomorphism of representations of Multi-digraph on 2 vertices,
         Homomorphism of representations of Multi-digraph on 2 vertices]
    """
    Element = QuiverRepHom

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, domain, codomain, category=None):
        """
        Initialize ``self``. Type ``QuiverHomSpace?`` for more information.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: H = Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            sage: TestSuite(H).run()
        """
        # The data in the class is stored in the following private variables:
        #
        # * _base
        #      The base ring of the representations M and N.
        # * _codomain
        #      The QuiverRep object of the codomain N.
        # * _domain
        #      The QuiverRep object of the domain M.
        # * _quiver
        #      The quiver of the representations M and N.
        # * _space
        #      A free module with ambient space.
        #
        # The free module _space is the homomorphism space.  The ambient space
        # is k^n where k is the base ring and n is the sum of the dimensions of
        # the spaces of homomorphisms between the free modules attached in M
        # and N to the vertices of the quiver.  Each coordinate represents a
        # single entry in one of those matrices.

        # Get the quiver and base ring and check that they are the same for
        # both modules
        if domain._semigroup != codomain._semigroup:
            raise ValueError("representations are not over the same quiver")
        self._quiver = domain._quiver
        self._semigroup = domain._semigroup

        # Check that the bases are compatible, and then initialise the homset:
        if codomain.base_ring() != domain.base_ring():
            raise ValueError("representations are not over the same base ring")
        Homset.__init__(self, domain, codomain, category=category,
                        base=domain.base_ring())

        # To compute the Hom Space we set up a 'generic' homomorphism where the
        # maps at each vertex are described by matrices whose entries are
        # variables.  Then the commutativity of edge diagrams gives us a
        # system of equations whose solution space is the Hom Space we're
        # looking for.  The variables will be numbered consecutively starting
        # at 0, ordered first by the vertex the matrix occurs at, then by row
        # then by column.  We'll have to keep track of which variables
        # correspond to which matrices.

        # eqs will count the number of equations in our system of equations,
        # varstart will be a list whose ith entry is the number of the
        # variable located at (0, 0) in the matrix assigned to the
        # ith vertex. (So varstart[0] will be 0.)
        eqs = 0
        verts = domain._quiver.vertices()
        varstart = [0] * (len(verts) + 1)

        # First assign to varstart the dimension of the matrix assigned to the
        # previous vertex.
        for v in verts:
            varstart[verts.index(v) + 1] = domain._spaces[v].dimension() * codomain._spaces[v].dimension()
        for e in domain._semigroup._sorted_edges:
            eqs += domain._spaces[e[0]].dimension() * codomain._spaces[e[1]].dimension()

        # After this cascading sum varstart[v] will be the sum of the
        # dimensions of the matrices assigned to vertices ordered before v.
        # This is equal to the number of the first variable assigned to v.
        for i in range(2, len(varstart)):
            varstart[i] += varstart[i - 1]

        # This will be the coefficient matrix for the system of equations.  We
        # start with all zeros and will fill in as we go.  We think of this
        # matrix as acting on the right so the columns correspond to equations,
        # the rows correspond to variables, and .kernel() will give a right
        # kernel as is needed.
        from sage.matrix.constructor import Matrix
        coef_mat = Matrix(codomain.base_ring(), varstart[-1], eqs)

        # eqn keeps track of what equation we are on.  If the maps X and Y are
        # assigned to an edge e and A and B are the matrices of variables that
        # describe the generic maps at the initial and final vertices of e
        # then commutativity of the edge diagram is described by the equation
        # AY = XB, or
        #
        #          Sum_k A_ik*Y_kj - Sum_k X_ik*B_kj == 0 for all i and j.
        #
        # Below we loop through these values of i,j,k and write the
        # coefficients of the equation above into the coefficient matrix.
        eqn = 0
        for e in domain._semigroup._sorted_edges:
            X = domain._maps[e].matrix()
            Y = codomain._maps[e].matrix()
            for i in range(X.nrows()):
                for j in range(Y.ncols()):
                    for k in range(Y.nrows()):
                        coef_mat[varstart[verts.index(e[0])] + i * Y.nrows() + k, eqn] = Y[k, j]
                    for k in range(X.ncols()):
                        coef_mat[varstart[verts.index(e[1])] + k * Y.ncols() + j, eqn] = -X[i, k]
                    eqn += 1

        # Now we can create the hom space
        self._space = coef_mat.kernel()

        # Bind identity if domain = codomain
        if domain is codomain:
            self.identity = self._identity

    @cached_method
    def zero(self):
        """
        Return the zero morphism.

        .. NOTE::

            It is needed to override the method inherited from
            the category of modules, because it would create
            a morphism that is of the wrong type and does not
            comply with :class:`~sage.quivers.morphism.QuiverRepHom`.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: H = Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            sage: H.zero() + H.an_element() == H.an_element()
            True
            sage: isinstance(H.zero(), H.element_class)
            True
        """
        return self()

    def _coerce_map_from_(self, other):
        r"""
        A coercion exists if and only if ``other`` is also a
        :class:`QuiverHomSpace` and there is a coercion from the
        domain of ``self`` to the domain of ``other`` and from the
        codomain of ``other`` to the codomain of ``self``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: S = Q.S(QQ, 1)
            sage: H1 = P.Hom(S)
            sage: H2 = (P/P.radical()).Hom(S)
            sage: H1.coerce_map_from(H2) # indirect doctest
            Coercion map:
              From: Dimension 1 QuiverHomSpace
              To:   Dimension 1 QuiverHomSpace
        """

        if not isinstance(other, QuiverHomSpace):
            return False
        if not other._domain.has_coerce_map_from(self._domain):
            return False
        if not self._codomain.has_coerce_map_from(other._codomain):
            return False
        return True

    def __call__(self, *data, **kwds):
        r"""
        A homomorphism of quiver representations (of one and the same
        quiver) is given by specifying, for each vertex of the quiver, a
        homomorphism of the spaces assigned to this vertex such that these
        homomorphisms commute with the edge maps. The domain and codomain
        of the homomorphism are required to be representations over the
        same quiver with the same base ring.

        INPUT:

        Usually, one would provide a single dict, list,
        :class:`QuiverRepElement` or :class:`QuiverRepHom` as arguments.
        The semantics is as follows:

          - list: ``data`` can be a list of images for the generators of
            the domain.  "Generators" means the output of the ``gens()``
            method.  An error will be generated if the map so defined
            is not equivariant with respect to the action of the quiver.
          - dictionary: ``data`` can be a dictionary associating to each
            vertex of the quiver either a homomorphism with domain and
            codomain the spaces associated to this vertex in the domain
            and codomain modules respectively, or a matrix defining such
            a homomorphism, or an object that sage can construct such a
            matrix from.  Not all vertices must be specified, unspecified
            vertices are assigned the zero map, and keys not corresponding
            to vertices of the quiver are ignored.  An error will be
            generated if these maps do not commute with the edge maps of
            the domain and codomain.
          - :class:`QuiverRepElement`: if the domain is a
            :class:`QuiverRep_with_path_basis` then ``data`` can be a single
            :class:`QuiverRepElement` belonging to the codomain.  The map
            is then defined by sending each path, ``p``, in the basis
            to ``data*p``.  If ``data`` is not an element of the codomain or
            the domain is not a :class:`QuiverRep_with_path_basis` then
            an error will be generated.
          - :class:`QuiverRepHom`: the input can also be a map `f : D \to C`
            such that there is a coercion from the domain of ``self`` to ``D``
            and from ``C`` to the codomain of ``self``.  The composition
            of these maps is the result.

        If there additionally are keyword arguments or if a
        :class:`QuiverRepHom` can not be created from the data, then the
        default call method of :class:`~sage.categories.homset.Homset`
        is called instead.

        OUTPUT:

        - :class:`QuiverRepHom`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: H = S.Hom(M)

        With no additional data this creates the zero map::

            sage: f = H() # indirect doctest
            sage: f.is_zero()
            True

        We must specify maps at the vertices to get a nonzero
        homomorphism.  Note that if the dimensions of the spaces assigned
        to the domain and codomain of a vertex are equal then Sage will
        construct the identity matrix from ``1``::

            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = H(maps2) # indirect doctest

        Here we create the same map by specifying images for the generators::

            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: h = H([x, y]) # indirect doctest
            sage: g == h
            True

        Here is an example of the same with a bigger identity matrix::

            sage: spaces3 = {2: QQ^2, 3: QQ^2}
            sage: maps3 = {(2, 3, 'c'): [[1, 0], [1, 0]]}
            sage: S3 = Q.representation(QQ, spaces3, maps3)
            sage: h3 = S3.Hom(M)({2: 1, 3: [[1], [0]]})
            sage: h3.get_map(2)
            Vector space morphism represented by the matrix:
            [1 0]
            [0 1]
            Domain: Vector space of dimension 2 over Rational Field
            Codomain: Vector space of dimension 2 over Rational Field

        If the domain is a module of type :class:`QuiverRep_with_path_basis`
        (for example, the indecomposable projectives) we can create maps by
        specifying a single image::

            sage: Proj = Q.P(GF(7), 3)
            sage: Simp = Q.S(GF(7), 3)
            sage: im = Simp({3: (1,)})
            sage: H2 = Proj.Hom(Simp)
            sage: H2(im).is_surjective() # indirect doctest
            True
        """
        if kwds or len(data) > 1:
            return super(Homset, self).__call__(*data, **kwds)

        if not data:
            return self.natural_map()

        data0 = data[0]
        if data0 is None or data0 == 0:
            data0 = {}
        try:
            return self.element_class(self._domain, self._codomain, data0)
        except (TypeError, ValueError):
            return super(QuiverHomSpace, self).__call__(*data, **kwds)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: Q.P(GF(3), 2).Hom(Q.S(GF(3), 2)) # indirect doctest
            Dimension 1 QuiverHomSpace
        """
        return "Dimension {} QuiverHomSpace".format(self._space.dimension())

    def natural_map(self):
        """
        The natural map from domain to codomain.

        This is the zero map.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: S.hom(M)      # indirect doctest
            Homomorphism of representations of Multi-digraph on 3 vertices
            sage: S.hom(M) == S.Hom(M).natural_map()
            True
        """
        return self.element_class(self._domain, self._codomain, {})

    def _identity(self):
        """
        Return the identity map.

        OUTPUT:

        - :class:`QuiverRepHom`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: H = P.Hom(P)
            sage: f = H.identity() # indirect doctest
            sage: f.is_isomorphism()
            True
        """
        from sage.matrix.constructor import Matrix
        maps = {v: Matrix(self._domain._spaces[v].dimension(),
                          self._domain._spaces[v].dimension(),
                          self._base.one())
                for v in self._quiver}
        return self.element_class(self._domain, self._codomain, maps)

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view and modify the representation data. #
    #                                                                         #
    ###########################################################################

    def base_ring(self):
        """
        Return the base ring of the representations.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: H = Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            sage: H.base_ring()
            Rational Field
        """
        return self._base

    def quiver(self):
        """
        Return the quiver of the representations.

        OUTPUT:

        - :class:`DiGraph`, the quiver of the representations

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: H = P.S(QQ, 2).Hom(P.P(QQ, 1))
            sage: H.quiver() is P.quiver()
            True
        """
        return self._quiver

    def domain(self):
        """
        Return the domain of the hom space.

        OUTPUT:

        - :class:`QuiverRep`, the domain of the Hom space

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: S = Q.S(QQ, 2)
            sage: H = S.Hom(Q.P(QQ, 1))
            sage: H.domain() is S
            True
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of the hom space.

        OUTPUT:

        - :class:`QuiverRep`, the codomain of the Hom space

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: H = Q.S(QQ, 2).Hom(P)
            sage: H.codomain() is P
            True
        """
        return self._codomain

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the representation.       #
    #                                                                         #
    ###########################################################################

    def dimension(self):
        """
        Return the dimension of the hom space.

        OUTPUT:

        - integer, the dimension

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: H = Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            sage: H.dimension()
            2
        """
        return self._space.dimension()

    def gens(self):
        """
        Return a list of generators of the hom space (as a `k`-vector
        space).

        OUTPUT:

        - list of :class:`QuiverRepHom` objects, the generators

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: H = Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            sage: H.gens()
            [Homomorphism of representations of Multi-digraph on 2 vertices,
             Homomorphism of representations of Multi-digraph on 2 vertices]
        """
        return [self.element_class(self._domain, self._codomain, f)
                for f in self._space.gens()]

    def coordinates(self, hom):
        """
        Return the coordinates of the map when expressed in terms of the
        generators (i. e., the output of the ``gens`` method) of the
        hom space.

        INPUT:

        - ``hom`` -- :class:`QuiverRepHom`

        OUTPUT:

        - list, the coordinates of the given map when written in terms of the
          generators of the :class:`QuiverHomSpace`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: S = Q.S(QQ, 2)
            sage: P = Q.P(QQ, 1)
            sage: H = S.Hom(P)
            sage: f = S.hom({2: [[1,-1]]}, P)
            sage: H.coordinates(f)
            [1, -1]
        """
        # Use the coordinates function on space
        return self._space.coordinates(hom._vector)

        ###########################################################################
        #                                                                         #
        # CONSTRUCTION FUNCTIONS                                                  #
        #    These functions create and return modules and homomorphisms.         #
        #                                                                         #
        ###########################################################################

    def _an_element_(self):
        """
        Return a homomorphism in the Hom space.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: S = Q.S(QQ, 2)
            sage: P = Q.P(QQ, 1)
            sage: H = S.Hom(P)
            sage: H.an_element() in H   # indirect doctest
            True
        """
        return self.element_class(self._domain, self._codomain, self._space.an_element())

    def left_module(self, basis=False):
        """
        Create the QuiverRep of ``self`` as a module over the opposite
        quiver.

        INPUT:

        - ``basis`` - bool. If ``False``, then only the module is
          returned.  If ``True``, then a tuple is returned.  The first
          element is the QuiverRep and the second element is a
          dictionary which associates to each vertex a list.  The
          elements of this list are the homomorphisms which correspond to
          the basis elements of that vertex in the module.

        OUTPUT:

        - :class:`QuiverRep` or tuple

        .. WARNING::

            The codomain of the Hom space must be a left module.

        .. NOTE::

            The left action of a path `e` on a map `f` is given by
            `(ef)(m) = ef(m)`.  This gives the Hom space its structure as
            a left module over the path algebra. This is then converted to
            a right module over the path algebra of the opposite quiver
            ``Q.reverse()`` and returned.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}}).path_semigroup()
            sage: P = Q.P(GF(3), 3)
            sage: A = Q.free_module(GF(3))
            sage: H = P.Hom(A)
            sage: H.dimension()
            6
            sage: M, basis_dict = H.left_module(true)
            sage: M.dimension_vector()
            (4, 1, 1)
            sage: Q.reverse().P(GF(3), 3).dimension_vector()
            (4, 1, 1)

        As lists start indexing at 0 the `i`-th vertex corresponds to the
        `(i-1)`-th entry of the dimension vector::

            sage: len(basis_dict[2]) == M.dimension_vector()[1]
            True
        """
        from sage.quivers.representation import QuiverRep
        if not self._codomain.is_left_module():
            raise ValueError("the codomain must be a left module")

        # Create the spaces
        spaces = {}
        for v in self._quiver:
            im_gens = [self([self._codomain.left_edge_action((v, v), f(x))
                             for x in self._domain.gens()])._vector
                       for f in self.gens()]
            spaces[v] = self._space.submodule(im_gens)

        # Create the maps
        maps = {}
        for e in self._semigroup._sorted_edges:
            e_op = (e[1], e[0], e[2])
            maps[e_op] = []
            for vec in spaces[e[1]].gens():
                vec_im = spaces[e_op[1]].coordinate_vector(self([self._codomain.left_edge_action(e, self(vec)(x))
                                                                 for x in self._domain.gens()])._vector)
                maps[e_op].append(vec_im)

        # Create and return the module (and the dict if desired)
        if basis:
            basis_dict = {}
            for v in self._quiver:
                basis_dict[v] = [self.element_class(self._domain, self._codomain, vec)
                                 for vec in spaces[v].gens()]
            return (QuiverRep(self._base, self._semigroup.reverse(), spaces, maps), basis_dict)
        else:
            return QuiverRep(self._base, self._semigroup.reverse(), spaces, maps)
