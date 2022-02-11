"""
Quiver Morphisms
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

from sage.categories.morphism import CallMorphism
from sage.matrix.constructor import Matrix


class QuiverRepHom(CallMorphism):
    r"""
    A homomorphism of quiver representations (of one and the same quiver)
    is given by specifying, for each vertex of the quiver, a homomorphism
    of the spaces assigned to this vertex such that these homomorphisms
    commute with the edge maps.  The domain and codomain of the
    homomorphism are required to be representations of the same quiver
    over the same base ring.

    INPUT:

    - ``domain`` -- :class:`QuiverRep`, the domain of the homomorphism

    - ``codomain`` -- :class:`QuiverRep`, the codomain of the homomorphism

    - ``data`` - dict, list, or :class:`QuiverRepElement`
      (default: empty dict),
      with the following meaning:

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
        :class:`QuiverRepElement` belonging to the codomain.  The map is
        then defined by sending each path, ``p``, in the basis to ``data*p``.
        If ``data`` is not an element of the codomain or the domain is not a
        :class:`QuiverRep_with_path_basis` then an error will be generated.
      - :class:`QuiverRepHom`: the input can also be a map `f : D \to C` such
        that there is a coercion from the domain of ``self`` to ``D``
        and from ``C`` to the codomain of ``self``.  The composition
        of these maps is the result.

    OUTPUT:

    - :class:`QuiverRepHom`

    EXAMPLES::

        sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
        sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
        sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
        sage: M = Q.representation(QQ, spaces, maps)
        sage: spaces2 = {2: QQ^1, 3: QQ^1}
        sage: S = Q.representation(QQ, spaces2)

    With no additional data this creates the zero map::

        sage: f = S.hom(M)
        sage: f.is_zero()
        True

    We must specify maps at the vertices to get a nonzero homomorphism.
    Note that if the dimensions of the spaces assigned to the domain and
    codomain of a vertex are equal then Sage will construct the identity
    matrix from ``1``::

        sage: maps2 = {2:[1, -1], 3:1}
        sage: g = S.hom(maps2, M)

    Here we create the same map by specifying images for the generators::

        sage: x = M({2: (1, -1)})
        sage: y = M({3: (1,)})
        sage: h = S.hom([x, y], M)
        sage: g == h
        True

    If the domain is a module of type QuiverRep_with_path_basis (for example,
    the indecomposable projectives) we can create maps by specifying a single
    image::

        sage: Proj = Q.P(GF(7), 3)
        sage: Simp = Q.S(GF(7), 3)
        sage: im = Simp({3: (1,)})
        sage: Proj.hom(im, Simp).is_surjective()
        True
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, domain, codomain, data={}):
        """
        Initialize ``self``. Type ``QuiverRepHom?`` for more information.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: f = S.hom(M)
            sage: f.is_zero()
            True
            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = S.hom(maps2, M)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: h = S.hom([x, y], M)
            sage: g == h
            True
            sage: Proj = Q.P(GF(7), 3)
            sage: Simp = Q.S(GF(7), 3)
            sage: im = Simp({3: (1,)})
            sage: Proj.hom(im, Simp).is_surjective()
            True

        ::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: H1 = Q.P(GF(3), 2).Hom(Q.S(GF(3), 2))
            sage: H2 = Q.P(GF(3), 2).Hom(Q.S(GF(3), 1))
            sage: H1.an_element() in H1   # indirect doctest
            True
        """
        # The data of a representation is held in the following private
        # variables:
        #
        # * _quiver
        #      The quiver of the representation.
        # * _base_ring
        #      The base ring of the representation.
        # * _domain
        #      The QuiverRep object that is the domain of the homomorphism.
        # * _codomain
        #      The QuiverRep object that is the codomain of the homomorphism.
        # * _vector
        #      A vector in some free module over the base ring of a length such
        #      that each coordinate corresponds to an entry in the matrix of a
        #      homomorphism attached to a vertex.
        #
        # The variable data can also be a vector of appropriate length.  When
        # this is the case it will be loaded directly into _vector and then
        # _assert_valid_hom is called.

        from sage.quivers.representation import QuiverRepElement, QuiverRep_with_path_basis

        self._domain = domain
        self._codomain = codomain
        self._quiver = domain._quiver
        self._base_ring = domain.base_ring()

        # Check that the quiver and base ring match
        if codomain._quiver != self._quiver:
            raise ValueError("the quivers of the domain and codomain must be equal")
        if codomain.base_ring() != self._base_ring:
            raise ValueError("the base ring of the domain and codomain must be equal")

        # Get the dimensions of the spaces
        mat_dims = {}
        domain_dims = {}
        codomain_dims = {}
        for v in self._quiver:
            domain_dims[v] = domain._spaces[v].dimension()
            codomain_dims[v] = codomain._spaces[v].dimension()
            mat_dims[v] = domain_dims[v] * codomain_dims[v]
        total_dim = sum(mat_dims.values())

        # Handle the case when data is a vector
        if data in self._base_ring**total_dim:
            self._vector = data
            self._assert_valid_hom()
            super(QuiverRepHom, self).__init__(domain.Hom(codomain))
            return

        # If data is not a dict, create one
        if isinstance(data, dict):
            maps_dict = data
        else:
            # If data is not a list create one, then create a dict from it
            if isinstance(data, list):
                im_list = data
            else:
                # If data is a QuiverRepHom, create a list from it
                if isinstance(data, QuiverRepHom):
                    f = data._domain.coerce_map_from(domain)
                    g = self._codomain.coerce_map_from(data._codomain)
                    im_list = [g(data(f(x))) for x in domain.gens()]

                # The only case left is that data is a QuiverRepElement
                else:
                    if not isinstance(data, QuiverRepElement):
                        raise TypeError("input data must be dictionary, list, "
                                        "QuiverRepElement or vector")
                    if not isinstance(domain, QuiverRep_with_path_basis):
                        raise TypeError("if data is a QuiverRepElement then domain "
                                        "must be a QuiverRep_with_path_basis.")
                    if data not in codomain:
                        raise ValueError("if data is a QuiverRepElement then it must "
                                         "be an element of codomain")
                    im_list = [codomain.right_edge_action(data, p) for v in domain._quiver for p in domain._bases[v]]

            # WARNING: This code assumes that the function QuiverRep.gens() returns
            # the generators ordered first by vertex and then by the order of the
            # gens() method of the space associated to that vertex.  In particular
            # this is the order that corresponds to how maps are represented via
            # matrices

            # Get the gens of the domain and check that im_list is the right length
            dom_gens = domain.gens()
            if len(im_list) != len(dom_gens):
                raise ValueError(("domain is dimension {} but only {} images"
                                  " were supplied").format(len(dom_gens), len(im_list)))

            # Get the matrices of the maps
            start_index = 0
            maps_dict = {}
            for v in self._quiver:
                maps_dict[v] = []
                dim = domain._spaces[v].dimension()
                for i in range(start_index, start_index + dim):
                    if im_list[i].support() and im_list[i].support() != [v]:
                        # If the element does not have the correct
                        # support raise an error here, otherwise we
                        # might create a valid hom that does not map
                        # the generators to the supplied images
                        txt = "generator supported at vertex {} cannot"
                        txt += " map to element with support {}"
                        raise ValueError(txt.format(v, im_list[i].support()))
                    else:
                        # If the support works out add the images coordinates
                        # as a row of the matrix
                        maps_dict[v].append(codomain._spaces[v].coordinates(im_list[i]._elems[v]))

                start_index += dim

        # Get the coordinates of the vector
        from sage.categories.map import is_Map
        vector = []
        for v in self._quiver:
            if v in maps_dict:
                if is_Map(maps_dict[v]):
                    try:
                        m = maps_dict[v].matrix()
                    except (AttributeError, ValueError):
                        gens_images = [codomain._spaces[v].coordinate_vector(maps_dict[v](x))
                                       for x in domain._spaces[v].gens()]
                        m = Matrix(self._base_ring, domain_dims[v], codomain_dims[v], gens_images)
                else:
                    m = Matrix(self._base_ring, domain_dims[v], codomain_dims[v], maps_dict[v])
            else:
                m = Matrix(self._base_ring, domain_dims[v], codomain_dims[v])
            for i in range(domain_dims[v]):
                vector += list(m[i])

        # Wrap as a vector, check it, and return
        self._vector = (self._base_ring**total_dim)(vector)
        self._assert_valid_hom()
        super(QuiverRepHom, self).__init__(domain.Hom(codomain))

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: S.hom(M) # indirect doctest
            Homomorphism of representations of Multi-digraph on 3 vertices
        """
        return "Homomorphism of representations of " + repr(self._quiver)

    def _call_(self, x):
        """
        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: h = S.hom([x, y], M)
            sage: h(S.gens()[0]) == x
            True
            sage: h(S.gens()[1]) == y
            True

        The following was an issue during work on :trac:`12630`::

            sage: Q = DiGraph({1: {}}).path_semigroup()
            sage: M = Q.I(GF(3), 1)
            sage: m = M.an_element()
            sage: R = M.quotient(M)
            sage: R(m)
            Element of quiver representation
        """
        from sage.quivers.representation import QuiverRepElement
        # Check the input
        if not isinstance(x, QuiverRepElement):
            raise ValueError("QuiverRepHom can only be called on QuiverRepElement")

        elements = dict((v, self.get_map(v)(x._elems[v])) for v in self._quiver)
        return self._codomain(elements)

    def __add__(left, right):
        """
        This function overloads the ``+`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: z = M.zero()
            sage: h = S.hom([x, z], M)
            sage: g = S.hom([z, z], M)
            sage: f = g + h
            sage: f(S.gens()[0]) == x
            True
            sage: f(S.gens()[1]) == z
            True
        """
        new_vector = left._vector + right._vector
        return left._domain.hom(new_vector, left._codomain)

    def __iadd__(self, other):
        """
        This function overloads the ``+=`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: z = M.zero()
            sage: h = S.hom([x, z], M)
            sage: g = S.hom([z, z], M)
            sage: g += h
            sage: g(S.gens()[0]) == x
            True
            sage: g(S.gens()[1]) == z
            True
        """
        self._vector += other._vector
        return self

    def __sub__(left, right):
        """
        This function overloads the ``-`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: z = M.zero()
            sage: h = S.hom([x, z], M)
            sage: g = S.hom([z, y], M)
            sage: f = h - g
            sage: f(S.gens()[0]) == x
            True
            sage: f(S.gens()[1]) == -y
            True
        """
        new_vector = left._vector - right._vector
        return left._domain.hom(new_vector, left._codomain)

    def __isub__(self, other):
        """
        This function overloads the ``-=`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: z = M.zero()
            sage: h = S.hom([x, z], M)
            sage: g = S.hom([z, y], M)
            sage: h -= g
            sage: h(S.gens()[0]) == x
            True
            sage: h(S.gens()[1]) == -y
            True
        """
        self._vector -= other._vector
        return self

    def __neg__(self):
        """
        This function overrides the unary ``-`` operator

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: h = S.hom([x, y], M)
            sage: g = -h
            sage: g(S.gens()[0]) == -x
            True
            sage: g(S.gens()[1]) == -y
            True
        """
        return self._domain.hom(-self._vector, self._codomain)

    def __pos__(self):
        """
        This function overrides the unary ``+`` operator

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: h = S.hom([x, y], M)
            sage: g = +h
            sage: g == h
            True
        """
        return self

    def __eq__(self, other):
        """
        This function overrides the ``==`` operator

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: g = S.hom([x, y], M)
            sage: h = S.hom([x, y], M)
            sage: g == h
            True
        """
        # A homomorphism can only be equal to another homomorphism between the
        # same domain and codomain
        if not isinstance(other, QuiverRepHom) or self._domain != other._domain or self._codomain != other._codomain:
            return False

        # If all that holds just check the vectors
        return self._vector == other._vector

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: g = S.hom([x, y], M)
            sage: H = hash(g)
        """
        return hash(tuple(self._vector))

    def __ne__(self, other):
        """
        This function overrides the ``!=`` operator

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: z = M.zero()
            sage: g = S.hom([x, y], M)
            sage: h = S.hom([x, z], M)
            sage: g != h
            True
        """
        # A homomorphism can only be equal to another homomorphism between the
        # same domain and codomain
        if not isinstance(other, QuiverRepHom) or self._domain != other._domain or self._codomain != other._codomain:
            return True

        # If all that holds just check the vectors
        return self._vector != other._vector

    def __bool__(self):
        """
        Return whether ``self`` is the zero morphism.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = M({2: (1, -1)})
            sage: y = M({3: (1,)})
            sage: z = M.zero()
            sage: g = S.hom([x, y], M)
            sage: h = S.hom([z, z], M)
            sage: bool(g)
            True
            sage: bool(h)
            False
        """
        return any(self._vector)

    __nonzero__ = __bool__

    def __mul__(self, other):
        """
        This function overrides the ``*`` operator

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: x = S.gens()[0]
            sage: y = S.gens()[1]
            sage: g = S.hom([x, y], S)
            sage: h = S.hom(S)
            sage: (g*h).is_zero()
            True
        """
        maps = {v: other.get_matrix(v) * self.get_matrix(v)
                for v in self._quiver}
        return other._domain.hom(maps, self._codomain)

    ###########################################################################
    #                                                                         #
    # WELL DEFINEDNESS FUNCTIONS                                              #
    #    These functions test and assert well definedness of the              #
    #    homomorphism.                                                        #
    #                                                                         #
    ###########################################################################

    def _assert_valid_hom(self):
        """
        Raise a ``ValueError`` if the homomorphism is not well defined.

        Specifically it checks that the domain and codomains of the maps are
        correct and that the edge diagrams commute.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = Q.representation(QQ, spaces2)
            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = S.hom(maps2, M) # indirect doctest
            sage: f = S.hom(maps2, S) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unable to convert {2: [1, -1], 3: 1} to an element of Dimension 2 QuiverHomSpace
        """
        # Check that the domain and codomains dimensions add correctly
        totaldim = 0
        for v in self._quiver:
            totaldim += self._domain._spaces[v].dimension() * self._codomain._spaces[v].dimension()
        if totaldim != len(self._vector):
            raise ValueError("dimensions do not match domain and codomain")

        # Check that the edge diagrams commute
        for e in self._domain._semigroup._sorted_edges:
            if self.get_matrix(e[0]) * self._codomain._maps[e].matrix() != self._domain._maps[e].matrix() * self.get_matrix(e[1]):
                raise ValueError("the diagram of edge {} does not commute".format(e))

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view the homomorphism data.              #
    #                                                                         #
    ###########################################################################

    def domain(self):
        """
        Return the domain of the homomorphism.

        OUTPUT:

        - :class:`QuiverRep`, the domain

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: S = Q.representation(QQ)
            sage: g = M.hom(S)
            sage: g.domain() is M
            True
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of the homomorphism.

        OUTPUT:

        - :class:`QuiverRep`, the codomain

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: S = Q.representation(QQ)
            sage: g = S.hom(M)
            sage: g.codomain() is M
            True
        """
        return self._codomain

    def get_matrix(self, vertex):
        """
        Return the matrix of the homomorphism attached to vertex
        ``vertex``.

        INPUT:

        - ``vertex`` -- integer, a vertex of the quiver

        OUTPUT:

        - matrix, the matrix representing the homomorphism associated to
          the given vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: I = Q.I(QQ, 3)
            sage: M = I/I.radical()
            sage: f = M.coerce_map_from(I)
            sage: f.get_matrix(1)
            [1 0]
            [0 1]
        """
        # Get dimensions
        startdim = 0
        for v in self._quiver:
            if v == vertex:
                break
            startdim += self._domain._spaces[v].dimension() * self._codomain._spaces[v].dimension()

        rows = self._domain._spaces[vertex].dimension()
        cols = self._codomain._spaces[vertex].dimension()

        # Slice out the matrix and return
        return Matrix(self._base_ring, rows, cols,
                      self._vector.list()[startdim:startdim + rows * cols])

    def get_map(self, vertex):
        """
        Return the homomorphism at the given vertex ``vertex``.

        INPUT:

        - ``vertex`` -- integer, a vertex of the quiver

        OUTPUT:

        - homomorphism, the homomorphism associated to the given vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: S = P/P.radical()
            sage: f = S.coerce_map_from(P)
            sage: f.get_map(1).is_bijective()
            True
        """
        return self._domain._spaces[vertex].hom(self.get_matrix(vertex), self._codomain._spaces[vertex])

    def quiver(self):
        """
        Return the quiver of the representations in the domain/codomain.

        OUTPUT:

        - :class:`DiGraph`, the quiver of the representations in the domain
          and codomain

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: f = P.hom({1: 1, 2: 1, 3: 1}, P)
            sage: f.quiver() is Q.quiver()
            True
        """
        return self._quiver

    def base_ring(self):
        """
        Return the base ring of the representation in the codomain.

        OUTPUT:

        - ring, the base ring of the codomain

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: f = P.hom({1: 1, 2: 1, 3: 1}, P)
            sage: f.base_ring() is QQ
            True
        """
        return self._base_ring

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the homomorphism.         #
    #                                                                         #
    ###########################################################################

    def is_injective(self):
        """
        Test whether the homomorphism is injective.

        OUTPUT:

        - bool, ``True`` if the homomorphism is injective, ``False`` otherwise

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: f = P.hom({1: 1, 2: 1, 3: 1}, P)
            sage: f.is_injective()
            True
            sage: g = P.hom(P)
            sage: g.is_injective()
            False
        """
        # The homomorphism is injective if and only if it is injective at every
        # vertex
        return not any(self.get_matrix(v).nullity() for v in self._quiver)

    def is_surjective(self):
        """
        Test whether the homomorphism is surjective.

        OUTPUT:

        - bool, ``True`` if the homomorphism is surjective, ``False`` otherwise

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: f = P.hom({1: 1, 2: 1, 3: 1}, P)
            sage: f.is_surjective()
            True
            sage: g = P.hom(P)
            sage: g.is_surjective()
            False
        """
        # The homomorphism is surjective if and only if it is surjective at
        # every vertex
        for v in self._quiver:
            m = self.get_matrix(v)
            if m.rank() != m.ncols():
                return False

        return True

    def is_isomorphism(self):
        """
        Test whether the homomorphism is an isomorphism.

        OUTPUT:

        - bool, ``True`` if the homomorphism is bijective, ``False``
          otherwise

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: f = P.hom({1: 1, 2: 1, 3: 1}, P)
            sage: f.is_isomorphism()
            True
            sage: g = P.hom(P)
            sage: g.is_isomorphism()
            False
        """
        # It's an iso if and only if it's an iso at every vertex
        return all(self.get_matrix(v).is_invertible() for v in self._quiver)

    def is_zero(self):
        """
        Test whether the homomorphism is the zero homomorphism.

        OUTPUT:

        - bool, ``True`` if the homomorphism is zero, ``False`` otherwise

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: f = P.hom({1: 1, 2: 1, 3: 1}, P)
            sage: f.is_zero()
            False
            sage: g = P.hom(P)
            sage: g.is_zero()
            True
        """
        # The homomorphism is zero if and only if it is zero at every vertex
        return all(self.get_matrix(v).is_zero() for v in self._quiver)

    def is_endomorphism(self):
        """
        Test whether the homomorphism is an endomorphism.

        OUTPUT:

        - bool, ``True`` if the domain equals the codomain, ``False``
          otherwise

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: f = P.hom({1: 1, 2: 1, 3: 1}, P)
            sage: f.is_endomorphism()
            True
            sage: S = P/P.radical()
            sage: g = S.coerce_map_from(P)
            sage: g.is_endomorphism()
            False
        """
        return self._domain == self._codomain

    def rank(self):
        """
        Return the rank of the homomorphism ``self`` (as a `k`-linear
        map).

        OUTPUT:

        - integer, the rank

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: S = P/P.radical()
            sage: f = S.coerce_map_from(P)
            sage: assert(f.rank() == 1)
        """
        # The rank is the sum of the ranks at each vertex
        return sum(self.get_matrix(v).rank() for v in self._quiver)

    ###########################################################################
    #                                                                         #
    # CONSTRUCTION FUNCTIONS                                                  #
    #    These functions create new homomorphisms, representations, and       #
    #    elements from the given homomorphism.                                #
    #                                                                         #
    ###########################################################################

    def kernel(self):
        """
        Return the kernel of ``self``.

        OUTPUT:

        - :class:`QuiverRep`, the kernel

        .. NOTE::

            To get the inclusion map of the kernel, ``K``, into the
            domain, ``D``, use ``D.coerce_map_from(K)``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^2, 3: QQ^1}
            sage: N = Q.representation(QQ, spaces2, {(2, 3, 'c'): [[1], [0]]})
            sage: maps2 = {2:[[1, 0], [0, 0]], 3:1}
            sage: g = N.hom(maps2, M)
            sage: g.kernel().dimension_vector()
            (0, 1, 0)
        """
        spaces = dict((v, self.get_map(v).kernel()) for v in self._quiver)
        return self._domain._submodule(spaces)

    def image(self):
        """
        Return the image of ``self``.

        OUTPUT:

        - :class:`QuiverRep`, the image

        .. NOTE::

            To get the inclusion map of the image, ``I``, into the
            codomain, ``C``, use ``C.coerce_map_from(I)``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^2, 3: QQ^1}
            sage: N = Q.representation(QQ, spaces2, {(2, 3, 'c'): [[1], [0]]})
            sage: maps2 = {2:[[1, 0], [0, 0]], 3:1}
            sage: g = N.hom(maps2, M)
            sage: g.image().dimension_vector()
            (0, 1, 1)
        """
        spaces = dict((v, self.get_map(v).image()) for v in self._quiver)
        return self._codomain._submodule(spaces)

    def cokernel(self):
        """
        Return the cokernel of ``self``.

        OUTPUT:

        - :class:`QuiverRep`, the cokernel

        .. NOTE::

            To get the factor map of the codomain, ``D``, onto the
            cokernel, ``C``, use ``C.coerce_map_from(D)``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: spaces2 = {2: QQ^2, 3: QQ^1}
            sage: N = Q.representation(QQ, spaces2, {(2, 3, 'c'): [[1], [0]]})
            sage: maps2 = {2:[[1, 0], [0, 0]], 3:1}
            sage: g = N.hom(maps2, M)
            sage: g.cokernel().dimension_vector()
            (2, 1, 0)
        """
        return self._codomain.quotient(self.image())

    def linear_dual(self):
        r"""
        Compute the linear dual `Df : DN \to DM` of
        ``self`` = `f : M \to N` where `D(-) = Hom_k(-, k)`.

        OUTPUT:

        - :class:`QuiverRepHom`, the map `Df : DN \to DM`

        .. NOTE::

            If `e` is an edge of the quiver `Q` and `g` is an element of
            `Hom_k(N, k)` then we let `(ga)(m) = g(ma)`.  This gives
            `Hom_k(N, k)` its structure as a module over the opposite
            quiver ``Q.reverse()``.  The map `Hom_k(N, k) \to Hom_k(M, k)`
            returned sends `g` to `gf`.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: S = P/P.radical()
            sage: f = S.coerce_map_from(P)

        The dual of a surjective map is injective and vice versa::

            sage: f.is_surjective()
            True
            sage: g = f.linear_dual()
            sage: g.is_injective()
            True

        The dual of a right module is a left module for the same quiver, Sage
        represents this as a right module for the opposite quiver::

            sage: g.quiver().path_semigroup() is Q.reverse()
            True

        The double dual of a map is the original representation::

            sage: g.linear_dual() == f
            True
        """
        # The effect of the functor D is that it just transposes the matrix of
        # a hom
        maps = dict((v, self.get_matrix(v).transpose()) for v in self._quiver)
        return self._codomain.linear_dual().hom(maps, self._domain.linear_dual())

    def algebraic_dual(self):
        r"""
        Compute the algebraic dual `f^t : N^t \to M^t` of
        ``self`` = `f : M \to N` where `(-)^t = Hom_Q(-, kQ)`.

        OUTPUT:

        - :class:`QuiverRepHom`, the map `f^t : N^t \to M^t`

        .. NOTE::

            If `e` is an edge of the quiver `Q` and `g` is an element of
            `Hom_Q(N, kQ)` then we let `(ge)(m) = eg(m)`.  This gives
            `Hom_Q(N, kQ)` its structure as a module over the opposite
            quiver ``Q.reverse()``.  The map
            `Hom_Q(N, kQ) \to Hom_Q(M, kQ)` returned sends `g` to `gf`.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b','c','d']}, 2:{4:['e','f']}, 3:{4:['g']}, 5:{2:['h','i']}}).path_semigroup()
            sage: P1 = Q.P(QQ, 4)
            sage: P1.algebraic_dual()
            Representation with dimension vector (5, 2, 1, 1, 4)

        The algebraic dual of an indecomposable projective is the indecomposable
        projective of the same vertex in the opposite quiver.

            sage: Q.reverse().P(QQ, 4)
            Representation with dimension vector (5, 2, 1, 1, 4)
        """
        # Get the domain, its basis, and the codomain
        domain, domain_gens = self._codomain.algebraic_dual(True)
        codomain, co_domain_gens = self._domain.algebraic_dual(True)

        # Find the images in the domain and create the module
        # H = QuiverHomSpace(self._domain, self._quiver.free_module(self._base_ring))
        im_gens = [codomain({v: (g * self)._vector})
                   for v in self._quiver for g in domain_gens[v]]
        return domain.hom(im_gens, codomain)

    def direct_sum(self, maps, return_maps=False, pinch=None):
        r"""
        Return the direct sum of ``self`` with the maps in the list ``maps``.

        INPUT:

        - ``maps`` -- :class:`QuiverRepHom` or list of :class:`QuiverRepHom`'s

        - ``return_maps`` -- bool (default: ``False``). If ``False``, then
          the return value is a :class:`QuiverRepHom` which is the direct sum
          of ``self`` with the :class:`QuiverRepHoms` in ``maps``.
          If ``True``, then the return value is a tuple of length either 3
          or 5.  The first entry of the tuple is the QuiverRepHom giving
          the direct sum.  If ``pinch`` is either ``None`` or
          ``'codomain'`` then the next two entries in the tuple are lists
          giving respectively the inclusion and the projection maps for
          the factors of the direct sum.  Summands are ordered as given
          in maps with ``self`` as the zeroth summand.  If ``pinch`` is
          either ``None`` or ``'domain'`` then the next two entries in the
          tuple are the inclusion and projection maps for the codomain.
          Thus if ``pinch`` is ``None`` then the tuple will have length 5.
          If ``pinch`` is either ``'domain'`` or ``'codomain'`` then the
          tuple will have length 3.

        - ``pinch`` -- string or ``None`` (default: ``None``). If this is
          equal to ``'domain'``, then the domains of ``self`` and the
          given maps must be equal.  The direct sum of `f: A \to B` and
          `g: A \to C` returned is then the map `A \to B \oplus C` defined
          by sending `x` to `(f(x), g(x))`.  If ``pinch`` equals
          ``'codomain'``, then the codomains of ``self`` and the given
          maps must be equal.  The direct sum of `f: A \to C` and
          `g: B \to C` returned is then the map `A \oplus B \to C` defined
          by sending `(x, y)` to `f(x) + g(y)`.  Finally, if ``pinch`` is
          anything other than ``'domain'`` or ``'codomain'``, then the
          direct sum of `f: A \to B` and `g: C \to D` returned is the map
          `A \oplus C \to B \oplus D` defined by sending `(x, y)` to
          `(f(x), g(y))`.

        OUTPUT:

        - :class:`QuiverRepHom` or tuple

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: P1 = Q.P(GF(3), 1)
            sage: P2 = Q.P(GF(3), 2)
            sage: S1 = P1/P1.radical()
            sage: S2 = P2/P2.radical()
            sage: pi1 = S1.coerce_map_from(P1)
            sage: pi2 = S2.coerce_map_from(P2)
            sage: f = pi1.direct_sum(pi2)
            sage: f.domain().dimension_vector() == Q.free_module(GF(3)).dimension_vector()
            True
            sage: f.is_surjective()
            True
            sage: id = P1.Hom(P1).identity()
            sage: g = pi1.direct_sum(id, pinch='domain')
            sage: g.is_surjective()
            False
        """
        # Get the list of maps to be summed
        if isinstance(maps, QuiverRepHom):
            maplist = [self, maps]
        else:
            maplist = [self] + maps

        # Check that the quivers/base rings are the same.  If pinching also
        # check that the domain/codomains are correct
        for x in maplist:
            if not isinstance(x, QuiverRepHom):
                raise TypeError("maps must be a QuiverRepHom or list of QuiverRepHoms")
            if self._quiver is not x._quiver:
                raise ValueError("Cannot direct sum maps from different quivers")
            if self._base_ring is not x._base_ring:
                raise ValueError("Base rings must be identical")
            if pinch == 'domain' and self._domain is not x._domain:
                raise ValueError("Cannot pinch maps, domains do not agree")
            if pinch == 'codomain' and self._codomain is not x._codomain:
                raise ValueError("Cannot pinch maps, codomains do not agree")

        # Get the sums and their maps
        if pinch == 'domain':
            domain = self._domain
        else:
            domain, d_incl, d_proj = self._domain.direct_sum([x._domain for x in maplist[1:]], return_maps=True)
        if pinch == 'codomain':
            codomain = self._codomain
        else:
            codomain, c_incl, c_proj = self._codomain.direct_sum([x._codomain for x in maplist[1:]], return_maps=True)

        # Start with the zero map
        result = domain.hom(codomain)

        # Add each factor
        for i in range(len(maplist)):
            if pinch == 'domain':
                result += c_incl[i] * maplist[i]
            elif pinch == 'codomain':
                result += maplist[i] * d_proj[i]
            else:
                result += c_incl[i] * maplist[i] * d_proj[i]

        # Return the results
        if return_maps:
            if pinch == 'domain':
                return (result, c_incl, c_proj)
            elif pinch == 'codomain':
                return (result, d_incl, d_proj)
            else:
                return (result, d_incl, d_proj, c_incl, c_proj)
        else:
            return result

    def lift(self, x):
        """
        Given an element `x` of the image, return an element of the domain
        that maps onto it under ``self``.

        INPUT:

        - ``x`` -- :class:`QuiverRepElement`

        OUTPUT:

        - :class:`QuiverRepElement`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c','d']}}).path_semigroup()
            sage: P = Q.P(QQ, 3)
            sage: S = P/P.radical()
            sage: proj = S.coerce_map_from(P)
            sage: x = S.an_element()
            sage: y = proj.lift(x)
            sage: proj(y) == x
            True
            sage: zero = S.hom(S, {})
            sage: zero.lift(x)
            Traceback (most recent call last):
            ...
            ValueError: element is not in the image
        """
        # Lift at each vertex
        elems = dict((v, self.get_map(v).lift(x._elems[v])) for v in self._quiver)
        return self._domain(elems)

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def scalar_mult(self, scalar):
        r"""
        Return the result of the scalar multiplication ``scalar * self``,
        where ``scalar`` is an element of the base ring `k`.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: f = M.Hom(M).an_element()
            sage: x = M.an_element()
            sage: g = f.scalar_mult(6)
            sage: g(x) == 6*f(x)
            True
        """
        return self._domain.hom(scalar * self._vector, self._codomain)

    def iscalar_mult(self, scalar):
        """
        Multiply ``self`` by ``scalar`` in place.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: f = M.Hom(M).an_element()
            sage: x = M.an_element()
            sage: y = f(x)
            sage: f.iscalar_mult(6)
            sage: f(x) == 6*y
            True
        """
        self._vector *= scalar
