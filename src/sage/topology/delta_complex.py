# -*- coding: utf-8 -*-
r"""
Finite Delta-complexes

AUTHORS:

- John H. Palmieri (2009-08)

This module implements the basic structure of finite
`\Delta`-complexes.  For full mathematical details, see Hatcher [Hat2002]_,
especially Section 2.1 and the Appendix on "Simplicial CW Structures".
As Hatcher points out, `\Delta`-complexes were first introduced by Eilenberg
and Zilber [EZ1950]_, although they called them "semi-simplicial complexes".

A `\Delta`-complex is a generalization of a :mod:`simplicial complex
<sage.homology.simplicial_complex>`; a `\Delta`-complex `X` consists
of sets `X_n` for each non-negative integer `n`, the elements of which
are called *n-simplices*, along with *face maps* between these sets of
simplices: for each `n` and for all `0 \leq i \leq n`, there are
functions `d_i` from `X_n` to `X_{n-1}`, with `d_i(s)` equal to the
`i`-th face of `s` for each simplex `s \in X_n`.  These maps must
satisfy the *simplicial identity*

  .. MATH::

    d_i d_j = d_{j-1} d_i \text{ for all } i<j.

Given a `\Delta`-complex, it has a *geometric realization*: a
topological space built by taking one topological `n`-simplex for each
element of `X_n`, and gluing them together as determined by the face
maps.

`\Delta`-complexes are an alternative to simplicial complexes.  Every
simplicial complex is automatically a `\Delta`-complex; in the other
direction, though, it seems in practice that one can often construct
`\Delta`-complex representations for spaces with many fewer simplices
than in a simplicial complex representation.  For example, the minimal
triangulation of a torus as a simplicial complex contains 14
triangles, 21 edges, and 7 vertices, while there is a `\Delta`-complex
representation of a torus using only 2 triangles, 3 edges, and 1
vertex.

.. note::

   This class derives from
   :class:`~sage.homology.cell_complex.GenericCellComplex`, and so
   inherits its methods.  Some of those methods are not listed here;
   see the :mod:`Generic Cell Complex <sage.homology.cell_complex>`
   page instead.
"""

from copy import copy
from sage.topology.cell_complex import GenericCellComplex
from sage.homology.chains import Chains, Cochains
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.matrix.constructor import matrix
from .simplicial_complex import Simplex, lattice_paths, SimplicialComplex
from sage.homology.chain_complex import ChainComplex
from sage.graphs.graph import Graph
from sage.arith.all import binomial
from sage.misc.cachefunc import cached_method


class DeltaComplex(GenericCellComplex):
    r"""
    Define a `\Delta`-complex.

    :param data: see below for a description of the options
    :param check_validity: If True, check that the simplicial identities hold.
    :type check_validity: boolean; optional, default True
    :return: a `\Delta`-complex

    Use ``data`` to define a `\Delta`-complex.  It may be in any of
    three forms:

    - ``data`` may be a dictionary indexed by simplices.  The value
      associated to a d-simplex `S` can be any of:

      - a list or tuple of (d-1)-simplices, where the ith entry is the
        ith face of S, given as a simplex,

      - another d-simplex `T`, in which case the ith face of `S` is
        declared to be the same as the ith face of `T`: `S` and `T`
        are glued along their entire boundary,

      - None or True or False or anything other than the previous two
        options, in which case the faces are just the ordinary faces
        of `S`.

      For example, consider the following::

        sage: n = 5
        sage: S5 = DeltaComplex({Simplex(n):True, Simplex(range(1,n+2)): Simplex(n)})
        sage: S5
        Delta complex with 6 vertices and 65 simplices

      The first entry in dictionary forming the argument to
      ``DeltaComplex`` says that there is an `n`-dimensional simplex
      with its ordinary boundary.  The second entry says that there is
      another simplex whose boundary is glued to that of the first
      one.  The resulting `\Delta`-complex is, of course, homeomorphic
      to an `n`-sphere, or actually a 5-sphere, since we defined `n`
      to be 5.  (Note that the second simplex here can be any
      `n`-dimensional simplex, as long as it is distinct from
      ``Simplex(n)``.)

      Let's compute its homology, and also compare it to the simplicial version::

        sage: S5.homology()
        {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: Z}
        sage: S5.f_vector()  # number of simplices in each dimension
        [1, 6, 15, 20, 15, 6, 2]
        sage: simplicial_complexes.Sphere(5).f_vector()
        [1, 7, 21, 35, 35, 21, 7]

      Both contain a single (-1)-simplex, the empty simplex; other
      than that, the `\Delta`-complex version contains fewer simplices
      than the simplicial one in each dimension.

      To construct a torus, use::

        sage: torus_dict = {Simplex([0,1,2]): True,
        ....:        Simplex([3,4,5]): (Simplex([0,1]), Simplex([0,2]), Simplex([1,2])),
        ....:        Simplex([0,1]): (Simplex(0), Simplex(0)),
        ....:        Simplex([0,2]): (Simplex(0), Simplex(0)),
        ....:        Simplex([1,2]): (Simplex(0), Simplex(0)),
        ....:        Simplex(0): ()}
        sage: T = DeltaComplex(torus_dict); T
        Delta complex with 1 vertex and 7 simplices
        sage: T.cohomology(base_ring=QQ)
        {0: Vector space of dimension 0 over Rational Field,
         1: Vector space of dimension 2 over Rational Field,
         2: Vector space of dimension 1 over Rational Field}

      This `\Delta`-complex consists of two triangles (given by
      ``Simplex([0,1,2])`` and ``Simplex([3,4,5])``); the boundary of
      the first is just its usual boundary: the 0th face is obtained
      by omitting the lowest numbered vertex, etc., and so the
      boundary consists of the edges ``[1,2]``, ``[0,2]``, and
      ``[0,1]``, in that order.  The boundary of the second is, on the
      one hand, computed the same way: the nth face is obtained by
      omitting the nth vertex.  On the other hand, the boundary is
      explicitly declared to be edges ``[0,1]``, ``[0,2]``, and
      ``[1,2]``, in that order.  This glues the second triangle to the
      first in the prescribed way.  The three edges each start and end
      at the single vertex, ``Simplex(0)``.

      .. image:: ../../media/torus_labelled.png

    - ``data`` may be nested lists or tuples.  The nth entry in the
      list is a list of the n-simplices in the complex, and each
      n-simplex is encoded as a list, the ith entry of which is its
      ith face.  Each face is represented by an integer, giving its
      index in the list of (n-1)-faces.  For example, consider this::

        sage: P = DeltaComplex( [ [(), ()],  [(1,0), (1,0), (0,0)],
        ....:                     [(1,0,2), (0, 1, 2)] ])

      The 0th entry in the list is ``[(), ()]``: there are two
      0-simplices, and their boundaries are empty.

      The 1st entry in the list is ``[(1,0), (1,0), (0,0)]``: there
      are three 1-simplices.  Two of them have boundary ``(1,0)``,
      which means that their 0th face is vertex 1 (in the list of
      vertices), and their 1st face is vertex 0.  The other edge has
      boundary ``(0,0)``, so it starts and ends at vertex 0.

      The 2nd entry in the list is ``[(1,0,2), (0,1,2)]``: there are
      two 2-simplices.  The first 2-simplex has boundary ``(1,0,2)``,
      meaning that its 0th face is edge 1 (in the list above), its 1st
      face is edge 0, and its 2nd face is edge 2; similarly for the
      2nd 2-simplex.

      If one draws two triangles and identifies them according to this
      description, the result is the real projective plane.

      .. image:: ../../media/rp2.png

      ::

        sage: P.homology(1)
        C2
        sage: P.cohomology(2)
        C2

      Closely related to this form for ``data`` is ``X.cells()``
      for a `\Delta`-complex ``X``: this is a dictionary, indexed by
      dimension ``d``, whose ``d``-th entry is a list of the
      ``d``-simplices, as a list::

        sage: P.cells()
        {-1: ((),),
         0: ((), ()),
         1: ((1, 0), (1, 0), (0, 0)),
         2: ((1, 0, 2), (0, 1, 2))}

    - ``data`` may be a dictionary indexed by integers.  For each
      integer `n`, the entry with key `n` is the list of
      `n`-simplices: this is the same format as is output by the
      :meth:`cells` method. ::

        sage: P = DeltaComplex( [ [(), ()],  [(1,0), (1,0), (0,0)],
        ....:                     [(1,0,2), (0, 1, 2)] ])
        sage: cells_dict = P.cells()
        sage: cells_dict
        {-1: ((),),
         0: ((), ()),
         1: ((1, 0), (1, 0), (0, 0)),
         2: ((1, 0, 2), (0, 1, 2))}
        sage: DeltaComplex(cells_dict)
        Delta complex with 2 vertices and 8 simplices
        sage: P == DeltaComplex(cells_dict)
        True

    Since `\Delta`-complexes are generalizations of simplicial
    complexes, any simplicial complex may be viewed as a
    `\Delta`-complex::

        sage: RP2 = simplicial_complexes.RealProjectivePlane()
        sage: RP2_delta = RP2.delta_complex()
        sage: RP2.f_vector()
        [1, 6, 15, 10]
        sage: RP2_delta.f_vector()
        [1, 6, 15, 10]

    Finally, `\Delta`-complex constructions for several familiar
    spaces are available as follows::

        sage: delta_complexes.Sphere(4)  # the 4-sphere
        Delta complex with 5 vertices and 33 simplices
        sage: delta_complexes.KleinBottle()
        Delta complex with 1 vertex and 7 simplices
        sage: delta_complexes.RealProjectivePlane()
        Delta complex with 2 vertices and 8 simplices

    Type ``delta_complexes.`` and then hit the TAB key to get the
    full list.
    """
    def __init__(self, data=None, check_validity=True):
        r"""
        Define a `\Delta`-complex.  See :class:`DeltaComplex` for more
        documentation.

        EXAMPLES::

            sage: X = DeltaComplex({Simplex(3):True, Simplex(range(1,5)): Simplex(3), Simplex(range(2,6)): Simplex(3)}); X  # indirect doctest
            Delta complex with 4 vertices and 18 simplices
            sage: X.homology()
            {0: 0, 1: 0, 2: 0, 3: Z x Z}
            sage: X == loads(dumps(X))
            True
        """
        def store_bdry(simplex, faces):
            r"""
            Given a simplex of dimension d and a list of boundaries
            (as other simplices), this stores each boundary face in
            new_data[d-1] if necessary, records the index of each
            boundary face in bdry_list, represents the simplex as
            bdry_list in new_data[d], and returns bdry_list.

            If the simplex is in the dictionary old_delayed, then it
            is already stored, temporarily, in new_data[d], so replace
            its temporary version with bdry_list.
            """
            bdry_list = []
            d = simplex.dimension()
            if d > 0:
                for f in faces:
                    if f in new_data[d-1]:
                        bdry_list.append(new_data[d-1].index(f))
                    else:
                        bdry_list.append(len(new_data[d-1]))
                        new_delayed[f] = len(new_data[d-1])
                        new_data[d-1].append(f)
                bdry_list = tuple(bdry_list)
            else:
                bdry_list = ()
            if simplex in old_delayed:
                idx = old_delayed[simplex]
                new_data[d][idx] = bdry_list
            else:
                new_data[d].append(bdry_list)
            return bdry_list

        new_data = {-1: ((),)}  # add the empty cell
        if data is None:
            pass
        else:
            if isinstance(data, (list, tuple)):
                dim = 0
                for s in data:
                    new_data[dim] = s
                    dim += 1
            elif isinstance(data, dict):
                if all(isinstance(a, (int, Integer)) for a in data):
                    # a dictionary indexed by integers
                    new_data = data
                    if -1 not in new_data:
                        new_data[-1] = ((),)  # add the empty cell
                else:
                    # else a dictionary indexed by simplices
                    dimension = max([f.dimension() for f in data])
                    old_data_by_dim = {}
                    for dim in range(0, dimension+1):
                        old_data_by_dim[dim] = []
                        new_data[dim] = []
                    for x in data:
                        if not isinstance(x, Simplex):
                            raise TypeError("Each key in the data dictionary must be a simplex.")
                        old_data_by_dim[x.dimension()].append(x)
                    old_delayed = {}
                    for dim in range(dimension, -1, -1):
                        new_delayed = {}
                        current = {}
                        for x in old_data_by_dim[dim]:
                            if x in data:
                                bdry = data[x]
                            else:
                                bdry = True
                            if isinstance(bdry, Simplex):
                                # case 1
                                # value is a simplex, so x is glued to the old
                                # one along its boundary.  So the boundary of
                                # x is the boundary of the old simplex.
                                if bdry in current:
                                    # if the old simplex is there, copy its boundary
                                    if x in old_delayed:
                                        idx = old_delayed[x]
                                        new_data[dim][idx] = current[bdry]
                                    else:
                                        new_data[dim].append(current[bdry])
                                elif bdry in data:
                                    # the old simplex has not yet been added to
                                    # new_data, but is in the data dictionary.  So
                                    # add it.
                                    current[bdry] = store_bdry(bdry, bdry.faces())
                                    new_data[dim].append(current[bdry])
                                else:
                                    raise ValueError("In the data dictionary, there is a value which is a simplex not already in the dictionary.  This is not allowed.")
                            elif isinstance(bdry, (list, tuple)):
                                # case 2
                                # boundary is a list or tuple
                                current[x] = store_bdry(x, bdry)
                            else:
                                # case 3
                                # no valid boundary specified, so the default
                                # boundary of x should be used
                                if x not in current:
                                    # x hasn't already been added, in case 1
                                    current[x] = store_bdry(x, x.faces())
                        old_delayed = new_delayed
                        if dim > 0:
                            old_data_by_dim[dim-1].extend(old_delayed.keys())
            else:
                raise ValueError("data is not a list, tuple, or dictionary")
        for n in new_data:
            new_data[n] = tuple(new_data[n])
        # at this point, new_data is a dictionary indexed by
        # dimension, with new_data[d] a list of "simplices" in
        # dimension d
        if check_validity:
            dim = max(new_data)
            for d in range(dim, 1, -1):
                for s in new_data[d]:  # s is a d-simplex
                    faces = new_data[d-1]
                    for j in range(d+1):
                        if not all(faces[s[j]][i] == faces[s[i]][j-1] for i in range(j)):
                            msg = "Simplicial identity d_i d_j = d_{j-1} d_i fails"
                            msg += " for j=%s, in dimension %s"%(j,d)
                            raise ValueError(msg)
        # self._cells_dict: dictionary indexed by dimension d: for
        # each d, have list or tuple of simplices, and for each
        # simplex, have list or tuple with its boundary (as the index
        # of an element in the list of (d-1)-simplices).
        self._cells_dict = new_data
        # self._is_subcomplex_of: if self is a subcomplex of another
        # Delta complex, record that other complex here, along with
        # data relating the cells in self to the cells in the
        # containing complex: for each dimension, a list of indices
        # specifying, for each cell in self, which cell it corresponds
        # to in the containing complex.
        self._is_subcomplex_of = None
        # self._complex: dictionary indexed by dimension d, base_ring,
        # etc.: differential from dim d to dim d-1 in the associated
        # chain complex.  thus to get the differential in the cochain
        # complex from dim d-1 to dim d, take the transpose of this
        # one.
        # self._complex = {}

    def subcomplex(self, data):
        r"""
        Create a subcomplex.

        :param data: a dictionary indexed by dimension or a list (or
          tuple); in either case, data[n] should be the list (or tuple
          or set) of the indices of the simplices to be included in
          the subcomplex.

        This automatically includes all faces of the simplices in
        ``data``, so you only have to specify the simplices which are
        maximal with respect to inclusion.

        EXAMPLES::

            sage: X = delta_complexes.Torus()
            sage: A = X.subcomplex({2: [0]})  # one of the triangles of X
            sage: X.homology(subcomplex=A)
            {0: 0, 1: 0, 2: Z}

        In the following, ``line`` is a line segment and ``ends`` is
        the complex consisting of its two endpoints, so the relative
        homology of the two is isomorphic to the homology of a circle::

            sage: line = delta_complexes.Simplex(1) # an edge
            sage: line.cells()
            {-1: ((),), 0: ((), ()), 1: ((0, 1),)}
            sage: ends = line.subcomplex({0: (0, 1)})
            sage: ends.cells()
            {-1: ((),), 0: ((), ())}
            sage: line.homology(subcomplex=ends)
            {0: 0, 1: Z}
        """
        if isinstance(data, (list, tuple)):
            data = dict(zip(range(len(data)), data))

        # new_dict: dictionary for constructing the subcomplex
        new_dict = {}
        # new_data: dictionary of all cells in the subcomplex: store
        # this with the subcomplex to make it fast to list the cells
        # in self which are not in the subcomplex.
        new_data = {}
        # max_dim: maximum dimension of cells being added
        max_dim = max(data.keys())
        # cells_to_add: in each dimension, add these cells to
        # new_dict.  start with the cells given in new_data and add
        # faces of cells one dimension higher.
        cells_to_add = data[max_dim]
        cells = self.cells()
        for d in range(max_dim, -1, -1):
            # cells_to_add is the set of indices of d-cells in self to
            # add to new_dict.
            cells_to_add = sorted(cells_to_add)
            # we add only these cells, so we need to translate their
            # indices from, for example, (0, 1, 4, 5) to (0, 1, 2, 3).
            # That is, when they appear as boundaries of (d+1)-cells,
            # we need to translate their indices in each (d+1)-cell.
            # Here is the key for that translation:
            translate = dict(zip(cells_to_add, range(len(cells_to_add))))
            new_dict[d] = []
            d_cells = cells_to_add
            new_data[d] = cells_to_add
            try:
                cells_to_add = set(new_data[d-1])  # begin to populate the (d-1)-cells
            except KeyError:
                cells_to_add = set([])
            for x in d_cells:
                if d+1 in new_dict:
                    old = new_dict[d+1]
                    new_dict[d+1] = []
                    for f in old:
                        new_dict[d+1].append(tuple([translate[n] for n in f]))
                new_dict[d].append(cells[d][x])
                cells_to_add.update(cells[d][x])
        new_cells = [new_dict[n] for n in range(0, max_dim+1)]
        sub = DeltaComplex(new_cells)
        sub._is_subcomplex_of = {self: new_data}
        return sub

    def __hash__(self):
        r"""
        TESTS::

            sage: hash(delta_complexes.Sphere(2)) == hash(delta_complexes.Sphere(2))
            True
            sage: hash(delta_complexes.Sphere(4)) == hash(delta_complexes.Sphere(4))
            True
        """
        return hash(frozenset(self._cells_dict.items()))

    def __eq__(self, right):
        r"""
        Two `\Delta`-complexes are equal, according to this, if they have
        the same ``_cells_dict``.

        EXAMPLES::

            sage: S4 = delta_complexes.Sphere(4)
            sage: S2 = delta_complexes.Sphere(2)
            sage: S4 == S2
            False
            sage: newS2 = DeltaComplex({Simplex(2):True, Simplex([8,12,17]): Simplex(2)})
            sage: newS2 == S2
            True
        """
        return self._cells_dict == right._cells_dict

    def __ne__(self, other):
        r"""
        Return ``True`` if ``self`` and ``other`` are not equal.

        EXAMPLES::

            sage: S4 = delta_complexes.Sphere(4)
            sage: S2 = delta_complexes.Sphere(2)
            sage: S4 != S2
            True
            sage: newS2 = DeltaComplex({Simplex(2):True, Simplex([8,12,17]): Simplex(2)})
            sage: newS2 != S2
            False
        """
        return not self.__eq__(other)

    def cells(self, subcomplex=None):
        r"""
        The cells of this `\Delta`-complex.

        :param subcomplex: a subcomplex of this complex
        :type subcomplex: optional, default None

        The cells of this `\Delta`-complex, in the form of a dictionary:
        the keys are integers, representing dimension, and the value
        associated to an integer d is the list of d-cells.  Each
        d-cell is further represented by a list, the ith entry of
        which gives the index of its ith face in the list of
        (d-1)-cells.

        If the optional argument ``subcomplex`` is present, then
        "return only the faces which are *not* in the subcomplex".  To
        preserve the indexing, which is necessary to compute the
        relative chain complex, this actually replaces the faces in
        ``subcomplex`` with ``None``.

        EXAMPLES::

            sage: S2 = delta_complexes.Sphere(2)
            sage: S2.cells()
            {-1: ((),),
             0: ((), (), ()),
             1: ((0, 1), (0, 2), (1, 2)),
             2: ((0, 1, 2), (0, 1, 2))}
            sage: A = S2.subcomplex({1: [0,2]}) # one edge
            sage: S2.cells(subcomplex=A)
            {-1: (None,),
             0: (None, None, None),
             1: (None, (0, 2), None),
             2: ((0, 1, 2), (0, 1, 2))}
        """
        cells = self._cells_dict.copy()
        if subcomplex is None:
            return cells
        if subcomplex._is_subcomplex_of is None or self not in subcomplex._is_subcomplex_of:
            if subcomplex == self:
                for d in range(-1, max(cells.keys())+1):
                    l = len(cells[d])
                    cells[d] = [None]*l   # get rid of all cells
                return cells
            else:
                raise ValueError("This is not a subcomplex of self.")
        else:
            subcomplex_cells = subcomplex._is_subcomplex_of[self]
            for d in range(0, max(subcomplex_cells.keys())+1):
                L = list(cells[d])
                for c in subcomplex_cells[d]:
                    L[c] = None
                cells[d] = tuple(L)
            cells[-1] = (None,)
        return cells

    def chain_complex(self, subcomplex=None, augmented=False,
                      verbose=False, check=False, dimensions=None,
                      base_ring=ZZ, cochain=False):
        r"""
        The chain complex associated to this `\Delta`-complex.

        :param dimensions: if None, compute the chain complex in all
           dimensions.  If a list or tuple of integers, compute the
           chain complex in those dimensions, setting the chain groups
           in all other dimensions to zero.  NOT IMPLEMENTED YET: this
           function always returns the entire chain complex
        :param base_ring: commutative ring
        :type base_ring: optional, default ZZ
        :param subcomplex: a subcomplex of this simplicial complex.
           Compute the chain complex relative to this subcomplex.
        :type subcomplex: optional, default empty
        :param augmented: If True, return the augmented chain complex
           (that is, include a class in dimension `-1` corresponding
           to the empty cell).  This is ignored if ``dimensions`` is
           specified or if ``subcomplex`` is nonempty.
        :type augmented: boolean; optional, default False
        :param cochain: If True, return the cochain complex (that is,
           the dual of the chain complex).
        :type cochain: boolean; optional, default False
        :param verbose: If True, print some messages as the chain
           complex is computed.
        :type verbose: boolean; optional, default False
        :param check: If True, make sure that the chain complex
           is actually a chain complex: the differentials are
           composable and their product is zero.
        :type check: boolean; optional, default False

        .. note::

           If subcomplex is nonempty, then the argument ``augmented``
           has no effect: the chain complex relative to a nonempty
           subcomplex is zero in dimension `-1`.

        EXAMPLES::

            sage: circle = delta_complexes.Sphere(1)
            sage: circle.chain_complex()
            Chain complex with at most 2 nonzero terms over Integer Ring
            sage: circle.chain_complex()._latex_()
            '\\Bold{Z}^{1} \\xrightarrow{d_{1}} \\Bold{Z}^{1}'
            sage: circle.chain_complex(base_ring=QQ, augmented=True)
            Chain complex with at most 3 nonzero terms over Rational Field
            sage: circle.homology(dim=1)
            Z
            sage: circle.cohomology(dim=1)
            Z
            sage: T = delta_complexes.Torus()
            sage: T.chain_complex(subcomplex=T)
            Trivial chain complex over Integer Ring
            sage: T.homology(subcomplex=T, algorithm='no_chomp')
            {0: 0, 1: 0, 2: 0}
            sage: A = T.subcomplex({2: [1]})  # one of the two triangles forming T
            sage: T.chain_complex(subcomplex=A)
            Chain complex with at most 1 nonzero terms over Integer Ring
            sage: T.homology(subcomplex=A)
            {0: 0, 1: 0, 2: Z}
        """
        if subcomplex is not None:
            # relative chain complex, so don't augment the chain complex
            augmented = False

        differentials = {}
        if augmented:
            empty_simplex = 1  # number of (-1)-dimensional simplices
        else:
            empty_simplex = 0
        vertices = self.n_cells(0, subcomplex=subcomplex)
        old = vertices
        old_real = [x for x in old if x is not None] # get rid of faces not in subcomplex
        n = len(old_real)
        differentials[0] = matrix(base_ring, empty_simplex, n, n*empty_simplex*[1])
        # current is list of simplices in dimension dim
        # current_real is list of simplices in dimension dim, with None filtered out
        # old is list of simplices in dimension dim-1
        # old_real is list of simplices in dimension dim-1, with None filtered out
        for dim in range(1,self.dimension()+1):
            current = list(self.n_cells(dim, subcomplex=subcomplex))
            current_real = [x for x in current if x is not None]
            i = 0
            i_real = 0
            translate = {}
            for s in old:
                if s is not None:
                    translate[i] = i_real
                    i_real += 1
                i += 1
            mat_dict = {}
            col = 0
            for s in current_real:
                sign = 1
                for row in s:
                    if old[row] is not None:
                        actual_row = translate[row]
                        if (actual_row,col) in mat_dict:
                            mat_dict[(actual_row,col)] += sign
                        else:
                            mat_dict[(actual_row,col)] = sign
                    sign *= -1
                col += 1
            differentials[dim] = matrix(base_ring, len(old_real), len(current_real), mat_dict)
            old = current
            old_real = current_real
        if cochain:
            cochain_diffs = {}
            for dim in differentials:
                cochain_diffs[dim-1] = differentials[dim].transpose()
            return ChainComplex(data=cochain_diffs, degree=1,
                                base_ring=base_ring, check=check)
        else:
            return ChainComplex(data=differentials, degree=-1,
                                base_ring=base_ring, check=check)

    def alexander_whitney(self, cell, dim_left):
        r"""
        Subdivide ``cell`` in this `\Delta`-complex into a pair of
        simplices.

        For an abstract simplex with vertices `v_0`, `v_1`, ...,
        `v_n`, then subdivide it into simplices `(v_0, v_1, ...,
        v_{dim_left})` and `(v_{dim_left}, v_{dim_left + 1}, ...,
        v_n)`. In a `\Delta`-complex, instead take iterated faces:
        take top faces to get the left factor, take bottom faces to
        get the right factor.

        INPUT:

        - ``cell`` -- a simplex in this complex, given as a pair
          ``(idx, tuple)``, where ``idx`` is its index in the list of
          cells in the given dimension, and ``tuple`` is the tuple of
          its faces

        - ``dim_left`` -- integer between 0 and one more than the
          dimension of this simplex

        OUTPUT: a list containing just the triple ``(1, left,
        right)``, where ``left`` and ``right`` are the two cells
        described above, each given as pairs ``(idx, tuple)``.

        EXAMPLES::

            sage: X = delta_complexes.Torus()
            sage: X.n_cells(2)
            [(1, 2, 0), (0, 2, 1)]
            sage: X.alexander_whitney((0, (1, 2, 0)), 1)
            [(1, (0, (0, 0)), (1, (0, 0)))]
            sage: X.alexander_whitney((0, (1, 2, 0)), 0)
            [(1, (0, ()), (0, (1, 2, 0)))]
            sage: X.alexander_whitney((1, (0, 2, 1)), 2)
            [(1, (1, (0, 2, 1)), (0, ()))]
        """
        dim = len(cell[1]) - 1
        left_cell = cell[1]
        idx_l = cell[0]
        for i in range(dim, dim_left, -1):
            idx_l = left_cell[i]
            left_cell = self.n_cells(i-1)[idx_l]
        right_cell = cell[1]
        idx_r = cell[0]
        for i in range(dim, dim - dim_left, -1):
            idx_r = right_cell[0]
            right_cell = self.n_cells(i-1)[idx_r]
        return [(ZZ.one(), (idx_l, left_cell), (idx_r, right_cell))]

    def n_skeleton(self, n):
        r"""
        The n-skeleton of this `\Delta`-complex.

        :param n: dimension
        :type n: non-negative integer

        EXAMPLES::

            sage: S3 = delta_complexes.Sphere(3)
            sage: S3.n_skeleton(1) # 1-skeleton of a tetrahedron
            Delta complex with 4 vertices and 11 simplices
            sage: S3.n_skeleton(1).dimension()
            1
            sage: S3.n_skeleton(1).homology()
            {0: 0, 1: Z x Z x Z}
        """
        if n >= self.dimension():
            return self
        else:
            data = []
            for d in range(n+1):
                data.append(self._cells_dict[d])
            return DeltaComplex(data)

    def graph(self):
        r"""
        The 1-skeleton of this `\Delta`-complex as a graph.

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: T.graph()
            Looped multi-graph on 1 vertex
            sage: S = delta_complexes.Sphere(2)
            sage: S.graph()
            Graph on 3 vertices
            sage: delta_complexes.Simplex(4).graph() == graphs.CompleteGraph(5)
            True
        """
        data = {}
        for vertex in range(len(self.n_cells(0))):
            data[vertex] = []
        for edge in self.n_cells(1):
            data[edge[0]].append(edge[1])
        return Graph(data)

    def join(self, other):
        r"""
        The join of this `\Delta`-complex with another one.

        :param other: another `\Delta`-complex (the right-hand
           factor)
        :return: the join ``self * other``

        The join of two `\Delta`-complexes `S` and `T` is the
        `\Delta`-complex `S*T` with simplices of the form `[v_0, ...,
        v_k, w_0, ..., w_n]` for all simplices `[v_0, ..., v_k]` in
        `S` and `[w_0, ..., w_n]` in `T`.  The faces are computed
        accordingly: the ith face of such a simplex is either `(d_i S)
        * T` if `i \leq k`, or `S * (d_{i-k-1} T)` if `i > k`.

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: S0 = delta_complexes.Sphere(0)
            sage: T.join(S0)  # the suspension of T
            Delta complex with 3 vertices and 21 simplices

        Compare to simplicial complexes::

            sage: K = delta_complexes.KleinBottle()
            sage: T_simp = simplicial_complexes.Torus()
            sage: K_simp = simplicial_complexes.KleinBottle()
            sage: T.join(K).homology()[3] == T_simp.join(K_simp).homology()[3] # long time (3 seconds)
            True

        The notation '*' may be used, as well::

            sage: S1 = delta_complexes.Sphere(1)
            sage: X = S1 * S1    # X is a 3-sphere
            sage: X.homology()
            {0: 0, 1: 0, 2: 0, 3: Z}
        """
        data = []
        # vertices of the join: the union of the vertices.  put the
        # vertices of self first, then the vertices of right.
        data.append(self.n_cells(0) + other.n_cells(0))
        bdries = {}
        for l_idx in range(len(self.n_cells(0))):
            bdries[(0,l_idx,-1,0)] = l_idx
        for r_idx in range(len(other.n_cells(0))):
            bdries[(-1,0,0,r_idx)] = len(self.n_cells(0)) + r_idx
        # dimension of the join:
        maxdim = self.dimension() + other.dimension() + 1
        # now for the d-cells, d>0:
        for d in range(1,maxdim+1):
            d_cells = []
            positions = {}
            new_idx = 0
            for k in range(-1,d+1):
                n = d-1-k
                # d=n+k.  need a k-cell from self and an n-cell from other
                if k == -1:
                    left = [()]
                else:
                    left = self.n_cells(k)
                l_idx = 0
                if n == -1:
                    right = [()]
                else:
                    right = other.n_cells(n)
                for l in left:
                    r_idx = 0
                    for r in right:
                        # store index of the new simplex in positions
                        positions[(k, l_idx, n, r_idx)] = new_idx
                        # form boundary of l*r and store it in d_cells
                        bdry = []
                        # first faces come from left-hand factor
                        if k == 0:
                            bdry.append(bdries[(-1, 0, n, r_idx)])
                        else:
                            for i in range(k+1):
                                bdry.append(bdries[(k-1, l[i], n, r_idx)])
                        # remaining faces come from right-hand factor
                        if n == 0:
                            bdry.append(bdries[(k, l_idx, -1, 0)])
                        else:
                            for i in range(n+1):
                                bdry.append(bdries[(k, l_idx, n-1, r[i])])
                        d_cells.append(tuple(bdry))
                        r_idx += 1
                        new_idx += 1
                    l_idx += 1
            data.append(d_cells)
            bdries = positions
        return DeltaComplex(data)

    # Use * to mean 'join':
    __mul__ = join

    def cone(self):
        r"""
        The cone on this `\Delta`-complex.

        The cone is the complex formed by adding a new vertex `C` and
        simplices of the form `[C, v_0, ..., v_k]` for every simplex
        `[v_0, ..., v_k]` in the original complex.  That is, the cone
        is the join of the original complex with a one-point complex.

        EXAMPLES::

            sage: K = delta_complexes.KleinBottle()
            sage: K.cone()
            Delta complex with 2 vertices and 14 simplices
            sage: K.cone().homology()
            {0: 0, 1: 0, 2: 0, 3: 0}
        """
        return self.join(delta_complexes.Simplex(0))

    def suspension(self, n=1):
        r"""
        The suspension of this `\Delta`-complex.

        :param n: suspend this many times.
        :type n: positive integer; optional, default 1

        The suspension is the complex formed by adding two new
        vertices `S_0` and `S_1` and simplices of the form `[S_0, v_0,
        ..., v_k]` and `[S_1, v_0, ..., v_k]` for every simplex `[v_0,
        ..., v_k]` in the original complex.  That is, the suspension
        is the join of the original complex with a two-point complex
        (the 0-sphere).

        EXAMPLES::

            sage: S = delta_complexes.Sphere(0)
            sage: S3 = S.suspension(3)  # the 3-sphere
            sage: S3.homology()
            {0: 0, 1: 0, 2: 0, 3: Z}
        """
        if n<0:
            raise ValueError("n must be non-negative.")
        if n==0:
            return self
        if n==1:
            return self.join(delta_complexes.Sphere(0))
        return self.suspension().suspension(int(n-1))

    def product(self, other):
        r"""
        The product of this `\Delta`-complex with another one.

        :param other: another `\Delta`-complex (the right-hand
           factor)
        :return: the product ``self x other``

        .. warning::

           If ``X`` and ``Y`` are `\Delta`-complexes, then ``X*Y``
           returns their join, not their product.

        EXAMPLES::

            sage: K = delta_complexes.KleinBottle()
            sage: X = K.product(K)
            sage: X.homology(1)
            Z x Z x C2 x C2
            sage: X.homology(2)
            Z x C2 x C2 x C2
            sage: X.homology(3)
            C2
            sage: X.homology(4)
            0
            sage: X.homology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 4 over Finite Field of size 2,
             2: Vector space of dimension 6 over Finite Field of size 2,
             3: Vector space of dimension 4 over Finite Field of size 2,
             4: Vector space of dimension 1 over Finite Field of size 2}
            sage: S1 = delta_complexes.Sphere(1)
            sage: K.product(S1).homology() == S1.product(K).homology()
            True
            sage: S1.product(S1) == delta_complexes.Torus()
            True
        """
        data = []
        bdries = {}
        # vertices: the vertices in the product are of the form (v,w)
        # for v a vertex in self, w a vertex in other
        vertices = []
        l_idx = 0
        for v in self.n_cells(0):
            r_idx = 0
            for w in other.n_cells(0):
                # one vertex for each pair (v,w)
                # store its indices in bdries; store its boundary in vertices
                bdries[(0,l_idx,0,r_idx, ((0,0),))] = len(vertices)
                vertices.append(()) # add new vertex (simplex with empty bdry)
                r_idx += 1
            l_idx += 1
        data.append(tuple(vertices))
        # dim of the product:
        maxdim = self.dimension() + other.dimension()
        # d-cells, d>0: these are obtained by taking products of cells
        # of dimensions k and n, where n+k >= d and n <= d, k <= d.
        simplices = []
        new = {}
        for d in range(1, maxdim+1):
            for k in range(d+1):
                for n in range(d-k,d+1):
                    k_idx = 0
                    for k_cell in self.n_cells(k):
                        n_idx = 0
                        for n_cell in other.n_cells(n):
                            # find d-dimensional faces in product of
                            # k_cell and n_cell.  to avoid repetition,
                            # only look for faces which use all
                            # vertices of each factor: the 'path'
                            # corresponding to each d-cell must hit
                            # every row and every column in the
                            # lattice.  (See the 'product' method for
                            # Simplex, as well as the function
                            # 'lattice_paths', in
                            # simplicial_complex.py.)
                            for path in lattice_paths(list(range(k + 1)),
                                                      list(range(n + 1)),
                                                      length=d+1):
                                path = tuple(path)
                                new[(k, k_idx, n, n_idx, path)] = len(simplices)
                                bdry_list = []
                                for i in range(d+1):
                                    face_path = path[:i] + path[i+1:]
                                    if ((i<d and path[i][0] == path[i+1][0]) or
                                        (i>0 and path[i][0] == path[i-1][0])):
                                        # this k-simplex
                                        k_face_idx = k_idx
                                        k_face_dim = k
                                    else:
                                        # face of this k-simplex
                                        k_face_idx = k_cell[path[i][0]]
                                        k_face_dim = k-1
                                        tail = []
                                        for j in range(i,d):
                                            tail.append((face_path[j][0]-1,
                                                       face_path[j][1]))
                                        face_path = face_path[:i] + tuple(tail)
                                    if ((i<d and path[i][1] == path[i+1][1]) or
                                        (i>0 and path[i][1] == path[i-1][1])):
                                        # this n-simplex
                                        n_face_idx = n_idx
                                        n_face_dim = n
                                    else:
                                        # face of this n-simplex
                                        n_face_idx = n_cell[path[i][1]]
                                        n_face_dim = n-1
                                        tail = []
                                        for j in range(i,d):
                                            tail.append((face_path[j][0],
                                                         face_path[j][1]-1))
                                        face_path = face_path[:i] + tuple(tail)
                                    bdry_list.append(bdries[(k_face_dim, k_face_idx,
                                                             n_face_dim, n_face_idx,
                                                             face_path)])
                                simplices.append(tuple(bdry_list))
                            n_idx += 1
                        k_idx += 1
            # add d-simplices to data, store d-simplices in bdries,
            # reset simplices
            data.append(tuple(simplices))
            bdries = new
            new = {}
            simplices = []
        return DeltaComplex(data)

    def disjoint_union(self, right):
        r"""
        The disjoint union of this `\Delta`-complex with another one.

        :param right: the other `\Delta`-complex (the right-hand factor)

        EXAMPLES::

            sage: S1 = delta_complexes.Sphere(1)
            sage: S2 = delta_complexes.Sphere(2)
            sage: S1.disjoint_union(S2).homology()
            {0: Z, 1: Z, 2: Z}
        """
        dim = max(self.dimension(), right.dimension())
        data = {}
        # in dimension n, append simplices of self with simplices of
        # right, but translate each entry of each right simplex: add
        # len(self.n_cells(n-1)) to it
        for n in range(dim, 0, -1):
            data[n] = list(self.n_cells(n))
            translate = len(self.n_cells(n-1))
            for f in right.n_cells(n):
                data[n].append(tuple([a+translate for a in f]))
        data[0] = self.n_cells(0) + right.n_cells(0)
        return DeltaComplex(data)

    def wedge(self, right):
        r"""
        The wedge (one-point union) of this `\Delta`-complex with
        another one.

        :param right: the other `\Delta`-complex (the right-hand factor)

        .. note::

            This operation is not well-defined if ``self`` or
            ``other`` is not path-connected.

        EXAMPLES::

            sage: S1 = delta_complexes.Sphere(1)
            sage: S2 = delta_complexes.Sphere(2)
            sage: S1.wedge(S2).homology()
            {0: 0, 1: Z, 2: Z}
        """
        data = self.disjoint_union(right).cells()
        left_verts = len(self.n_cells(0))
        translate = {}
        for i in range(left_verts):
            translate[i] = i
        translate[left_verts] = 0
        for i in range(left_verts + 1, left_verts + len(right.n_cells(0))):
            translate[i] = i-1
        data[0] = data[0][:-1]
        edges = []
        for e in data[1]:
            edges.append([translate[a] for a in e])
        data[1] = edges
        return DeltaComplex(data)

    def connected_sum(self, other):
        r"""
        Return the connected sum of self with other.

        :param other: another `\Delta`-complex
        :return: the connected sum ``self # other``

        .. warning::

           This does not check that self and other are manifolds.  It
           doesn't even check that their facets all have the same
           dimension.  It just chooses top-dimensional simplices from
           each complex, checks that they have the same dimension,
           removes them, and glues the remaining pieces together.
           Since a (more or less) random facet is chosen from each
           complex, this method may return random results if applied
           to non-manifolds, depending on which facet is chosen.

        ALGORITHM:

        Pick a top-dimensional simplex from each complex.  Check to
        see if there are any identifications on either simplex, using
        the :meth:`_is_glued` method.  If there are no
        identifications, remove the simplices and glue the remaining
        parts of complexes along their boundary.  If there are
        identifications on a simplex, subdivide it repeatedly (using
        :meth:`elementary_subdivision`) until some piece has no
        identifications.

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: S2 = delta_complexes.Sphere(2)
            sage: T.connected_sum(S2).cohomology() == T.cohomology()
            True
            sage: RP2 = delta_complexes.RealProjectivePlane()
            sage: T.connected_sum(RP2).homology(1)
            Z x Z x C2
            sage: T.connected_sum(RP2).homology(2)
            0
            sage: RP2.connected_sum(RP2).connected_sum(RP2).homology(1)
            Z x Z x C2
        """
        if not self.dimension() == other.dimension():
            raise ValueError("Complexes are not of the same dimension.")
        dim = self.dimension()
        # Look at the last simplex in the list of top-dimensional
        # simplices for each complex.  If there are identifications on
        # either of these simplices, subdivide until there are no more
        # identifications.
        Left = self
        while Left._is_glued():
            Left = Left.elementary_subdivision()
        Right = other
        while Right._is_glued():
            Right = Right.elementary_subdivision()
        # remove last top-dimensional face from each one and glue.
        data = {}
        for n in Left.cells():
            data[n] = list(Left.cells()[n])
        right_cells = Right.cells()
        data[dim] = data[dim][:-1]
        left_simplex = Left.n_cells(dim)[-1]
        right_simplex = Right.n_cells(dim)[-1]
        # renaming: dictionary for translating all simplices of Right
        renaming = dict(zip(right_simplex, left_simplex))
        # process_now: cells to be reindexed and added to data
        process_now = right_cells[dim][:-1]
        for n in range(dim, 0, -1):
            # glued: dictionary of just the simplices being glued
            glued = copy(renaming)
            # process_later: cells one dim lower to be added to data
            process_later = []
            old_idx = 0
            new_idx = len(data[n-1])
            # build 'renaming'
            for s in right_cells[n-1]:
                if old_idx not in renaming:
                    process_later.append(s)
                    renaming[old_idx] = new_idx
                    new_idx += 1
                old_idx += 1
            # reindex all simplices to be processed and add them to data
            for s in process_now:
                data[n].append(tuple([renaming[i] for i in s]))
            # set up for next loop, one dimension down
            renaming = {}
            process_now = process_later
            for f in glued:
                renaming.update(dict(zip(right_cells[n-1][f], data[n-1][glued[f]])))
        # deal with vertices separately.  we just need to add enough
        # vertices: all the vertices from Right, minus the number
        # being glued, which should be dim+1, the number of vertices
        # in the simplex of dimension dim being glued.
        for i in range(len(right_cells[0]) - dim - 1):
            data[0].append(())
        return DeltaComplex(data)

    def elementary_subdivision(self, idx=-1):
        r"""
        Perform an "elementary subdivision" on a top-dimensional
        simplex in this `\Delta`-complex.  If the optional argument
        ``idx`` is present, it specifies the index (in the list of
        top-dimensional simplices) of the simplex to subdivide.  If
        not present, subdivide the last entry in this list.

        :param idx: index specifying which simplex to subdivide
        :type idx: integer; optional, default -1
        :return: `\Delta`-complex with one simplex subdivided.

        *Elementary subdivision* of a simplex means replacing that
        simplex with the cone on its boundary.  That is, given a
        `\Delta`-complex containing a `d`-simplex `S` with vertices
        `v_0`, ..., `v_d`, form a new `\Delta`-complex by

        - removing `S`
        - adding a vertex `w` (thought of as being in the interior of `S`)
        - adding all simplices with vertices `v_{i_0}`, ...,
          `v_{i_k}`, `w`, preserving any identifications present
          along the boundary of `S`

        The algorithm for achieving this uses
        :meth:`_epi_from_standard_simplex` to keep track of simplices
        (with multiplicity) and what their faces are: this method
        defines a surjection `\pi` from the standard `d`-simplex to
        `S`.  So first remove `S` and add a new vertex `w`, say at the
        end of the old list of vertices.  Then for each vertex `v` in
        the standard `d`-simplex, add an edge from `\pi(v)` to `w`;
        for each edge `(v_0, v_1)` in the standard `d`-simplex, add a
        triangle `(\pi(v_0), \pi(v_1), w)`, etc.

        Note that given an `n`-simplex `(v_0, v_1, ..., v_n)` in the
        standard `d`-simplex, the faces of the new `(n+1)`-simplex are
        given by removing vertices, one at a time, from `(\pi(v_0),
        ..., \pi(v_n), w)`.  These are either the image of the old
        `n`-simplex (if `w` is removed) or the various new
        `n`-simplices added in the previous dimension.  So keep track
        of what's added in dimension `n` for use in computing the
        faces in dimension `n+1`.

        In contrast with barycentric subdivision, note that only the
        interior of `S` has been changed; this allows for subdivision
        of a single top-dimensional simplex without subdividing every
        simplex in the complex.

        The term "elementary subdivision" is taken from p. 112 in John
        M. Lee's book [Lee2011]_.

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: T.n_cells(2)
            [(1, 2, 0), (0, 2, 1)]
            sage: T.elementary_subdivision(0)  # subdivide first triangle
            Delta complex with 2 vertices and 13 simplices
            sage: X = T.elementary_subdivision(); X  # subdivide last triangle
            Delta complex with 2 vertices and 13 simplices
            sage: X.elementary_subdivision()
            Delta complex with 3 vertices and 19 simplices
            sage: X.homology() == T.homology()
            True
        """
        pi = self._epi_from_standard_simplex(idx=idx)
        cells_dict = {}
        old_cells = self.cells()
        for n in old_cells:
            cells_dict[n] = list(old_cells[n])
        dim = self.dimension()
        # cells of standard simplex of dimension dim
        std_cells = SimplicialComplex([Simplex(dim)]).delta_complex(sort_simplices=True).cells()
        # adjust zero-cells so they're distinct
        std_cells[0] = tuple([[n] for n in range(dim + 1)])
        # remove the cell being subdivided
        cells_dict[dim].pop(idx)
        # add the new vertex "w"
        cells_dict[0].append(())
        # added_cells: dict indexed by (n-1)-cells, with value the
        # corresponding new n-cell.
        added_cells = {(): len(cells_dict[0])-1}
        for n in range(0, dim):
            new_cells = {}
            # for each n-cell in the standard simplex, add an
            # (n+1)-cell to the subdivided complex.
            try:
                simplices = sorted(pi[n])
            except TypeError:
                simplices = pi[n]
            for simplex in simplices:
                # compute the faces of the new (n+1)-cell.
                cell = []
                for i in simplex:
                    if n > 0:
                        bdry = tuple(std_cells[n-1][i])
                    else:
                        bdry = ()
                    cell.append(added_cells[bdry])
                # last face is the image of the old simplex)
                cell.append(pi[n][simplex])
                cell = tuple(cell)
                cells_dict[n+1].append(cell)
                new_cells[simplex] = len(cells_dict[n+1])-1
            added_cells = new_cells
        return DeltaComplex(cells_dict)

    def _epi_from_standard_simplex(self, idx=-1, dim=None):
        r"""
        Construct an epimorphism from a standard simplex to a
        top-dimensional simplex in this `\Delta`-complex.

        If the optional argument ``dim`` is not ``None``, then
        construct the map to a simplex with this dimension.  If the
        optional argument ``idx`` is present, it specifies which
        simplex to use by giving its index in the list of simplices of
        the appropriate dimension; if not present, use the last
        simplex in this list.

        This is used by :meth:`elementary_subdivision`.

        :param idx: index specifying which simplex to examine
        :type idx: integer; optional, default -1
        :return: boolean, True if the boundary of the simplex has any
          identifications
        :param dim: dimension of simplex to consider
        :type dim: integer; optional, default = dim of complex

        Suppose that the dimension is `d`. The map is given by a
        dictionary indexed by dimension: in dimension `i`, its value
        is a dictionary specifying, for each `i`-simplex in the
        domain, the corresponding `i`-simplex in the codomain.  The
        vertices are specified as their indices in the lists of
        simplices in each complex; the same goes for all of the
        simplices in the codomain.  The simplices of dimension 1 or
        higher in the domain are listed explicitly (in the form of
        entries from the output of :meth:`cells`).

        In this function, the "standard simplex" is defined to be
        ``simplicial_complexes.Simplex(d).delta_complex(sort_simplices=True)``.

        EXAMPLES:

        The `\Delta`-complex model for a torus has two triangles and
        three edges, but only one vertex.  So a surjection from the
        standard 2-simplex to either of the triangles is a bijection
        in dimension 1, but in dimension 0, sends all three vertices
        to the same place::

            sage: T = delta_complexes.Torus()
            sage: sorted(T._epi_from_standard_simplex()[1].items())
            [((1, 0), 1), ((2, 0), 2), ((2, 1), 0)]
            sage: sorted(T._epi_from_standard_simplex()[0].items())
            [((0,), 0), ((1,), 0), ((2,), 0)]
        """
        if dim is None:
            dim = self.dimension()
        # the output is easier to read if the entries are non-negative.
        if idx == -1:
            idx = len(self.n_cells(dim)) - 1
        simplex = SimplicialComplex([Simplex(dim)]).delta_complex(sort_simplices=True)
        simplex_cells = simplex.cells()
        self_cells = self.cells()
        if dim > 0:
            map = {dim: {tuple(simplex_cells[dim][0]): idx}}
        else:
            map = {dim: {(0,): idx}}
        faces_dict = map[dim]
        for n in range(dim, 0, -1):
            n_cells = faces_dict
            faces_dict = {}
            for cell in n_cells:
                if n > 1:
                    faces = [tuple(simplex_cells[n-1][cell[j]]) for j in range(0,n+1)]
                    one_cell =  dict(zip(faces, self_cells[n][n_cells[cell]]))
                else:
                    temp =  dict(zip(cell, self_cells[n][n_cells[cell]]))
                    one_cell = {}
                    for j in temp:
                        one_cell[(j,)] = temp[j]
                for j in one_cell:
                    if j not in faces_dict:
                        faces_dict[j] = one_cell[j]
            map[n-1] = faces_dict
        return map

    def _is_glued(self, idx=-1, dim=None):
        r"""
        ``True`` if there is any gluing along the boundary of a
        top-dimensional simplex in this `\Delta`-complex.

        If the optional argument ``idx`` is present, it specifies
        which simplex to consider by giving its index in the list of
        top-dimensional simplices; if not present, look at the last
        simplex in this list.  If the optional argument ``dim`` is
        present, it specifies the dimension of the simplex to
        consider; if not present, look at a top-dimensional simplex.

        This is used by :meth:`connected_sum`.

        :param idx: index specifying which simplex to examine
        :type idx: integer; optional, default -1
        :return: boolean, True if the boundary of the simplex has any
          identifications
        :param dim: dimension of simplex to consider
        :type dim: integer; optional, default = dim of complex

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: T._is_glued()
            True
            sage: S = delta_complexes.Simplex(3)
            sage: S._is_glued()
            False
        """
        if dim is None:
            dim = self.dimension()

        simplex = self.n_cells(dim)[idx]
        i = self.dimension() - 1
        i_faces = set(simplex)
        # if there are enough i_faces, then no gluing is evident so far
        not_glued = (len(i_faces) == binomial(dim+1, i+1))
        while not_glued and i > 0:
            # count the (i-1) cells and compare to (n+1) choose i.
            old_faces = i_faces
            i_faces = set([])
            all_cells = self.n_cells(i)
            for face in old_faces:
                i_faces.update(all_cells[face])
            not_glued = (len(i_faces) == binomial(dim+1, i))
            i = i-1
        return not not_glued

    def face_poset(self):
        r"""
        The face poset of this `\Delta`-complex, the poset of
        nonempty cells, ordered by inclusion.

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: T.face_poset()
            Finite poset containing 6 elements
        """
        from sage.combinat.posets.posets import Poset
        # given the structure of self.cells(), it's easier to compute
        # the dual poset, then reverse it at the end.
        dim = self.dimension()
        covers = {}
        # store each n-simplex as a pair (n, idx).
        for n in range(dim, 0, -1):
            idx = 0
            for s in self.n_cells(n):
                covers[(n, idx)] = list(set([(n-1, i) for i in s]))
                idx += 1
        # deal with vertices separately: they have no covers (in the
        # dual poset).
        idx = 0
        for s in self.n_cells(0):
            covers[(0, idx)] = []
            idx += 1
        return Poset(Poset(covers).hasse_diagram().reverse())

    # implement using the definition?  the simplices are obtained by
    # taking chains of inclusions of simplices, etc.  have to work out
    # the faces and identifications.
    def barycentric_subdivision(self):
        r"""
        Not implemented.

        EXAMPLES::

            sage: K = delta_complexes.KleinBottle()
            sage: K.barycentric_subdivision()
            Traceback (most recent call last):
            ...
            NotImplementedError: Barycentric subdivisions are not implemented for Delta complexes.
        """
        raise NotImplementedError("Barycentric subdivisions are not implemented for Delta complexes.")

    def n_chains(self, n, base_ring=None, cochains=False):
        r"""
        Return the free module of chains in degree ``n`` over ``base_ring``.

        INPUT:

        - ``n`` -- integer
        - ``base_ring`` -- ring (optional, default `\ZZ`)
        - ``cochains`` -- boolean (optional, default ``False``); if
          ``True``, return cochains instead

        Since the list of `n`-cells for a `\Delta`-complex may have
        some ambiguity -- for example, the list of edges may look like
        ``[(0, 0), (0, 0), (0, 0)]`` if each edge starts and ends at
        vertex 0 -- we record the indices of the cells along with
        their tuples. So the basis of chains in such a case would look
        like ``[(0, (0, 0)), (1, (0, 0)), (2, (0, 0))]``.

        The only difference between chains and cochains is notation:
        the dual cochain to the chain basis element ``b`` is written
        as ``\chi_b``.

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: T.n_chains(1, QQ)
            Free module generated by {(0, (0, 0)), (1, (0, 0)), (2, (0, 0))} over Rational Field
            sage: list(T.n_chains(1, QQ, cochains=False).basis())
            [(0, (0, 0)), (1, (0, 0)), (2, (0, 0))]
            sage: list(T.n_chains(1, QQ, cochains=True).basis())
            [\chi_(0, (0, 0)), \chi_(1, (0, 0)), \chi_(2, (0, 0))]
        """
        n_cells = tuple(enumerate(self.n_cells(n)))
        if cochains:
            return Cochains(self, n, n_cells, base_ring)
        else:
            return Chains(self, n, n_cells, base_ring)

    # the second barycentric subdivision is a simplicial complex.  implement this somehow?
#     def simplicial_complex(self):
#         X = self.barycentric_subdivision().barycentric_subdivision()
#         find facets of X and return SimplicialComplex(facets)

    # This is cached for speed reasons: it can be very slow to run
    # this function.
    @cached_method
    def algebraic_topological_model(self, base_ring=None):
        r"""
        Algebraic topological model for this `\Delta`-complex with
        coefficients in ``base_ring``.

        The term "algebraic topological model" is defined by Pilarczyk
        and Réal [PR2015]_.

        INPUT:

        - ``base_ring`` - coefficient ring (optional, default
          ``QQ``). Must be a field.

        Denote by `C` the chain complex associated to this
        `\Delta`-complex. The algebraic topological model is a chain complex
        `M` with zero differential, with the same homology as `C`,
        along with chain maps `\pi: C \to M` and `\iota: M \to C`
        satisfying `\iota \pi = 1_M` and `\pi \iota` chain homotopic
        to `1_C`. The chain homotopy `\phi` must satisfy

        - `\phi \phi = 0`,
        - `\pi \phi = 0`,
        - `\phi \iota = 0`.

        Such a chain homotopy is called a *chain contraction*.

        OUTPUT: a pair consisting of

        - chain contraction ``phi`` associated to `C`, `M`, `\pi`, and
          `\iota`
        - the chain complex `M`

        Note that from the chain contraction ``phi``, one can recover the
        chain maps `\pi` and `\iota` via ``phi.pi()`` and
        ``phi.iota()``. Then one can recover `C` and `M` from, for
        example, ``phi.pi().domain()`` and ``phi.pi().codomain()``,
        respectively.

        EXAMPLES::

            sage: RP2 = delta_complexes.RealProjectivePlane()
            sage: phi, M = RP2.algebraic_topological_model(GF(2))
            sage: M.homology()
            {0: Vector space of dimension 1 over Finite Field of size 2,
             1: Vector space of dimension 1 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2}
            sage: T = delta_complexes.Torus()
            sage: phi, M = T.algebraic_topological_model(QQ)
            sage: M.homology()
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 2 over Rational Field,
             2: Vector space of dimension 1 over Rational Field}
        """
        from sage.homology.algebraic_topological_model import algebraic_topological_model_delta_complex
        if base_ring is None:
            base_ring = QQ
        return algebraic_topological_model_delta_complex(self, base_ring)

    def _string_constants(self):
        r"""
        Tuple containing the name of the type of complex, and the
        singular and plural of the name of the cells from which it is
        built.  This is used in constructing the string representation.

        EXAMPLES::

            sage: T = delta_complexes.Torus()
            sage: T._string_constants()
            ('Delta', 'simplex', 'simplices')
        """
        return ('Delta', 'simplex', 'simplices')


class DeltaComplexExamples():
    r"""
    Some examples of `\Delta`-complexes.

    Here are the available examples; you can also type
    ``delta_complexes.``  and hit TAB to get a list::

        Sphere
        Torus
        RealProjectivePlane
        KleinBottle
        Simplex
        SurfaceOfGenus

    EXAMPLES::

        sage: S = delta_complexes.Sphere(6) # the 6-sphere
        sage: S.dimension()
        6
        sage: S.cohomology(6)
        Z
        sage: delta_complexes.Torus() == delta_complexes.Sphere(3)
        False
    """

    def Sphere(self,n):
        r"""
        A `\Delta`-complex representation of the `n`-dimensional sphere,
        formed by gluing two `n`-simplices along their boundary,
        except in dimension 1, in which case it is a single 1-simplex
        starting and ending at the same vertex.

        :param n: dimension of the sphere

        EXAMPLES::

            sage: delta_complexes.Sphere(4).cohomology(4, base_ring=GF(3))
            Vector space of dimension 1 over Finite Field of size 3
        """
        if n == 1:
            return DeltaComplex([ [()], [(0, 0)] ])
        else:
            return DeltaComplex({Simplex(n): True, Simplex(range(1,n+2)): Simplex(n)})

    def Torus(self):
        r"""
        A `\Delta`-complex representation of the torus, consisting of one
        vertex, three edges, and two triangles.

        .. image:: ../../media/torus.png

        EXAMPLES::

            sage: delta_complexes.Torus().homology(1)
            Z x Z
        """
        return DeltaComplex(( ((),),  ((0,0), (0,0), (0,0)), ((1,2,0), (0, 2, 1))))

    def RealProjectivePlane(self):
        r"""
        A `\Delta`-complex representation of the real projective plane,
        consisting of two vertices, three edges, and two triangles.

        .. image:: ../../media/rp2.png

        EXAMPLES::

            sage: P = delta_complexes.RealProjectivePlane()
            sage: P.cohomology(1)
            0
            sage: P.cohomology(2)
            C2
            sage: P.cohomology(dim=1, base_ring=GF(2))
            Vector space of dimension 1 over Finite Field of size 2
            sage: P.cohomology(dim=2, base_ring=GF(2))
            Vector space of dimension 1 over Finite Field of size 2
        """
        return DeltaComplex(( ((), ()),  ((1,0), (1,0), (0,0)), ((1,0,2), (0,1,2))))

    def KleinBottle(self):
        r"""
        A `\Delta`-complex representation of the Klein bottle, consisting
        of one vertex, three edges, and two triangles.

        .. image:: ../../media/klein.png

        EXAMPLES::

            sage: delta_complexes.KleinBottle()
            Delta complex with 1 vertex and 7 simplices
        """
        return DeltaComplex(( ((),),  ((0,0), (0,0), (0,0)), ((1,2,0), (0, 1, 2))))

    def Simplex(self, n):
        r"""
        A `\Delta`-complex representation of an `n`-simplex,
        consisting of a single `n`-simplex and its faces.  (This is
        the same as the simplicial complex representation available by
        using ``simplicial_complexes.Simplex(n)``.)

        EXAMPLES::

            sage: delta_complexes.Simplex(3)
            Delta complex with 4 vertices and 16 simplices
        """
        return DeltaComplex({Simplex(n): True})

    def SurfaceOfGenus(self, g, orientable=True):
        r"""
        A surface of genus g as a `\Delta`-complex.

        :param g: the genus
        :type g: non-negative integer
        :param orientable: whether the surface should be orientable
        :type orientable: bool, optional, default ``True``

        In the orientable case, return a sphere if `g` is zero, and
        otherwise return a `g`-fold connected sum of a torus with
        itself.

        In the non-orientable case, raise an error if `g` is zero.  If
        `g` is positive, return a `g`-fold connected sum of a
        real projective plane with itself.

        EXAMPLES::

            sage: delta_complexes.SurfaceOfGenus(1, orientable=False)
            Delta complex with 2 vertices and 8 simplices
            sage: delta_complexes.SurfaceOfGenus(3, orientable=False).homology(1)
            Z x Z x C2
            sage: delta_complexes.SurfaceOfGenus(3, orientable=False).homology(2)
            0

        Compare to simplicial complexes::

            sage: delta_g4 = delta_complexes.SurfaceOfGenus(4)
            sage: delta_g4.f_vector()
            [1, 3, 27, 18]
            sage: simpl_g4 = simplicial_complexes.SurfaceOfGenus(4)
            sage: simpl_g4.f_vector()
            [1, 19, 75, 50]
            sage: delta_g4.homology() == simpl_g4.homology()
            True
        """
        try:
            g = Integer(g)
        except TypeError:
            raise ValueError("genus must be a non-negative integer")
        if g < 0:
            raise ValueError("genus must be a non-negative integer")
        if g == 0:
            if not orientable:
                raise ValueError("no non-orientable surface of genus zero")
            else:
                return delta_complexes.Sphere(2)
        if orientable:
            X = delta_complexes.Torus()
        else:
            X = delta_complexes.RealProjectivePlane()
        S = X
        for i in range(g-1):
            S = S.connected_sum(X)
        return S

delta_complexes = DeltaComplexExamples()
