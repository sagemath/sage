r"""
Generic cell complexes

AUTHORS:

- John H. Palmieri (2009-08)

This module defines a class of abstract finite cell complexes.  This
is meant as a base class from which other classes (like
:class:`~sage.homology.simplicial_complex.SimplicialComplex`,
:class:`~sage.homology.cubical_complex.CubicalComplex`, and
:class:`~sage.homology.delta_complex.DeltaComplex`) should derive.  As
such, most of its properties are not implemented.  It is meant for use
by developers producing new classes, not casual users.

.. NOTE::

    Keywords for :meth:`~GenericCellComplex.chain_complex`,
    :meth:`~GenericCellComplex.homology`, etc.: any keywords given to
    the :meth:`~GenericCellComplex.homology` method get passed on to
    the :meth:`~GenericCellComplex.chain_complex` method and also to
    the constructor for chain complexes in
    :class:`sage.homology.chain_complex.ChainComplex_class <ChainComplex>`,
    as well as its associated
    :meth:`~sage.homology.chain_complex.ChainComplex_class.homology` method.
    This means that those keywords should have consistent meaning in
    all of those situations.  It also means that it is easy to
    implement new keywords: for example, if you implement a new
    keyword for the
    :meth:`sage.homology.chain_complex.ChainComplex_class.homology` method,
    then it will be automatically accessible through the
    :meth:`~GenericCellComplex.homology` method for cell complexes --
    just make sure it gets documented.
"""


from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

class GenericCellComplex(SageObject):
    r"""
    Class of abstract cell complexes.

    This is meant to be used by developers to produce new classes, not
    by casual users.  Classes which derive from this are
    :class:`~sage.homology.simplicial_complex.SimplicialComplex`,
    :class:`~sage.homology.delta_complex.DeltaComplex`, and
    :class:`~sage.homology.cubical_complex.CubicalComplex`.

    Most of the methods here are not implemented, but probably should
    be implemented in a derived class.  Most of the other methods call
    a non-implemented one; their docstrings contain examples from
    derived classes in which the various methods have been defined.
    For example, :meth:`homology` calls :meth:`chain_complex`; the
    class :class:`~sage.homology.delta_complex.DeltaComplex`
    implements
    :meth:`~sage.homology.delta_complex.DeltaComplex.chain_complex`,
    and so the :meth:`homology` method here is illustrated with
    examples involving `\Delta`-complexes.

    EXAMPLES:

    It's hard to give informative examples of the base class, since
    essentially nothing is implemented. ::

        sage: from sage.homology.cell_complex import GenericCellComplex
        sage: A = GenericCellComplex()
    """
    def __cmp__(self,right):
        """
        Comparisons of cell complexes are not implemented.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A == B # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    ############################################################
    # self.cells() and related methods
    ############################################################

    def cells(self, subcomplex=None):
        """
        The cells of this cell complex, in the form of a dictionary:
        the keys are integers, representing dimension, and the value
        associated to an integer `d` is the set of `d`-cells.  If the
        optional argument ``subcomplex`` is present, then return only
        the faces which are *not* in the subcomplex.

        :param subcomplex: a subcomplex of this cell complex.  Return
           the cells which are not in this subcomplex.
        :type subcomplex: optional, default None

        This is not implemented in general; it should be implemented
        in any derived class.  When implementing, see the warning in
        the :meth:`dimension` method.

        This method is used by various other methods, such as
        :meth:`n_cells` and :meth:`f_vector`.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.cells()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def dimension(self):
        """
        The dimension of this cell complex: the maximum
        dimension of its cells.

        .. WARNING::

          If the :meth:`cells` method calls :meth:`dimension`,
          then you'll get an infinite loop.  So either don't use
          :meth:`dimension` or override :meth:`dimension`.

        EXAMPLES::

            sage: simplicial_complexes.RandomComplex(d=5, n=8).dimension()
            5
            sage: delta_complexes.Sphere(3).dimension()
            3
            sage: T = cubical_complexes.Torus()
            sage: T.product(T).dimension()
            4
        """
        try:
            return max([x.dimension() for x in self._facets])
        except AttributeError:
            return max(self.cells())

    def n_cells(self, n, subcomplex=None):
        """
        List of cells of dimension ``n`` of this cell complex.
        If the optional argument ``subcomplex`` is present, then
        return the ``n``-dimensional faces which are *not* in the
        subcomplex.

        :param n: the dimension
        :type n: non-negative integer
        :param subcomplex: a subcomplex of this cell complex. Return
           the cells which are not in this subcomplex.
        :type subcomplex: optional, default ``None``

        EXAMPLES::

            sage: simplicial_complexes.Simplex(2).n_cells(1)
            [(1, 2), (0, 2), (0, 1)]
            sage: delta_complexes.Torus().n_cells(1)
            [(0, 0), (0, 0), (0, 0)]
            sage: cubical_complexes.Cube(1).n_cells(0)
            [[1,1], [0,0]]
        """
        if n in self.cells(subcomplex):
            return list(self.cells(subcomplex)[n])
        else:
            # don't barf if someone asks for n_cells in a dimension where there are none
            return []

    def f_vector(self):
        """
        The `f`-vector of this cell complex: a list whose `n^{th}`
        item is the number of `(n-1)`-cells.  Note that, like all
        lists in Sage, this is indexed starting at 0: the 0th element
        in this list is the number of `(-1)`-cells (which is 1: the
        empty cell is the only `(-1)`-cell).

        EXAMPLES::

            sage: simplicial_complexes.KleinBottle().f_vector()
            [1, 8, 24, 16]
            sage: delta_complexes.KleinBottle().f_vector()
            [1, 1, 3, 2]
            sage: cubical_complexes.KleinBottle().f_vector()
            [1, 42, 84, 42]
        """
        return [self._f_dict()[n] for n in range(-1, self.dimension()+1)]

    def _f_dict(self):
        """
        The `f`-vector of this cell complex as a dictionary: the
        item associated to an integer `n` is the number of the
        `n`-cells.

        EXAMPLES::

            sage: simplicial_complexes.KleinBottle()._f_dict()[1]
            24
            sage: delta_complexes.KleinBottle()._f_dict()[1]
            3
        """
        answer = {}
        answer[-1] = 1
        for n in range(self.dimension() + 1):
            answer[n] = len(self.n_cells(n))
        return answer

    def euler_characteristic(self):
        r"""
        The Euler characteristic of this cell complex: the
        alternating sum over `n \geq 0` of the number of
        `n`-cells.

        EXAMPLES::

            sage: simplicial_complexes.Simplex(5).euler_characteristic()
            1
            sage: delta_complexes.Sphere(6).euler_characteristic()
            2
            sage: cubical_complexes.KleinBottle().euler_characteristic()
            0
        """
        return sum([(-1)**n * self.f_vector()[n+1] for n in range(self.dimension() + 1)])

    ############################################################
    # end of methods using self.cells()
    ############################################################

    def product(self, right, rename_vertices=True):
        """
        The (Cartesian) product of this cell complex with another one.

        Products are not implemented for general cell complexes.  They
        may be implemented in some derived classes (like simplicial
        complexes).

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.product(B)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def disjoint_union(self, right):
        """
        The disjoint union of this simplicial complex with another one.

        :param right: the other simplicial complex (the right-hand factor)

        Disjoint unions are not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.disjoint_union(B)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def wedge(self, right):
        """
        The wedge (one-point union) of this simplicial complex with
        another one.

        :param right: the other simplicial complex (the right-hand factor)

        Wedges are not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.wedge(B)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    ############################################################
    # self.join() and related methods
    ############################################################

    def join(self, right, **kwds):
        """
        The join of this cell complex with another one.

        :param right: the other simplicial complex (the right-hand factor)

        Joins are not implemented for general cell complexes.  They
        may be implemented in some derived classes (like simplicial
        complexes).

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.join(B)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    # for some classes, you may want * to mean join:
    ###
    # __mul__ = join

    # the cone on X is the join of X with a point.  See
    # simplicial_complex.py for one implementation.
    ###
    # def cone(self):
    #     return self.join(POINT)

    # the suspension of X is the join of X with the 0-sphere (two
    # points).  See simplicial_complex.py for one implementation.
    ###
    # def suspension(self, n=1):
    #     """
    #     The suspension of this cell complex.
    #
    #     INPUT:
    #
    #     -  ``n`` - positive integer (optional, default 1): suspend this
    #        many times.
    #     """
    #     raise NotImplementedError

    ############################################################
    # end of methods using self.join()
    ############################################################

    ############################################################
    # chain complexes, homology
    ############################################################

    def chain_complex(self, **kwds):
        """
        This is not implemented for general cell complexes.

        Some keywords to possibly implement in a derived class:

        - ``subcomplex`` -- a subcomplex: compute the relative chain complex
        - ``augmented`` -- a bool: whether to return the augmented complex
        - ``verbose`` -- a bool: whether to print informational messages as
          the chain complex is being computed
        - ``check_diffs`` -- a bool: whether to check that the each
          composite of two consecutive differentials is zero
        -  ``dimensions`` -- if ``None``, compute the chain complex in all
           dimensions.  If a list or tuple of integers, compute the
           chain complex in those dimensions, setting the chain groups
           in all other dimensions to zero.

        Definitely implement the following:

        -  ``base_ring`` -- commutative ring (optional, default ZZ)
        -  ``cochain`` -- a bool: whether to return the cochain complex

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.chain_complex()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def homology(self, dim=None, **kwds):
        r"""
        The reduced homology of this cell complex.

        :param dim: If None, then return the homology in every
           dimension.  If ``dim`` is an integer or list, return the
           homology in the given dimensions.  (Actually, if ``dim`` is
           a list, return the homology in the range from ``min(dim)``
           to ``max(dim)``.)
        :type dim: integer or list of integers or None; optional,
           default None
        :param base_ring: commutative ring, must be ZZ or a field.
        :type base_ring: optional, default ZZ
        :param subcomplex: a subcomplex of this simplicial complex.
           Compute homology relative to this subcomplex.
        :type subcomplex: optional, default empty
        :param generators: If ``True``, return generators for the homology
           groups along with the groups.  NOTE: Since :trac:`6100`, the result
           may not be what you expect when not using CHomP since its return
           is in terms of the chain complex.
        :type generators: boolean; optional, default False
        :param cohomology: If True, compute cohomology rather than homology.
        :type cohomology: boolean; optional, default False
        :param algorithm: The options are 'auto', 'dhsw', 'pari' or 'no_chomp'.
           See below for a description of what they mean.
        :type algorithm: string; optional, default 'auto'
        :param verbose: If True, print some messages as the homology is
           computed.
        :type verbose: boolean; optional, default False

        .. note::

            The keyword arguments to this function get passed on to
            :meth:``chain_complex`` and its homology.

        ALGORITHM:

        If ``algorithm`` is set to 'auto' (the default), then use
        CHomP if available.  (CHomP is available at the web page
        http://chomp.rutgers.edu/.  It is also an experimental package
        for Sage.)

        CHomP computes homology, not cohomology, and only works over
        the integers or finite prime fields.  Therefore if any of
        these conditions fails, or if CHomP is not present, or if
        ``algorithm`` is set to 'no_chomp', go to plan B: if ``self``
        has a ``_homology`` method -- each simplicial complex has
        this, for example -- then call that.  Such a method implements
        specialized algorithms for the particular type of cell
        complex.

        Otherwise, move on to plan C: compute the chain complex of
        ``self`` and compute its homology groups.  To do this: over a
        field, just compute ranks and nullities, thus obtaining
        dimensions of the homology groups as vector spaces.  Over the
        integers, compute Smith normal form of the boundary matrices
        defining the chain complex according to the value of
        ``algorithm``.  If ``algorithm`` is 'auto' or 'no_chomp', then
        for each relatively small matrix, use the standard Sage
        method, which calls the Pari package.  For any large matrix,
        reduce it using the Dumas, Heckenbach, Saunders, and Welker
        elimination algorithm: see
        :func:`sage.homology.matrix_utils.dhsw_snf` for details.

        Finally, ``algorithm`` may also be 'pari' or 'dhsw', which
        forces the named algorithm to be used regardless of the size
        of the matrices and regardless of whether CHomP is available.

        As of this writing, CHomP is by far the fastest option,
        followed by the 'auto' or 'no_chomp' setting of using the
        Dumas, Heckenbach, Saunders, and Welker elimination algorithm
        for large matrices and Pari for small ones.

        EXAMPLES::

            sage: P = delta_complexes.RealProjectivePlane()
            sage: P.homology()
            {0: 0, 1: C2, 2: 0}
            sage: P.homology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 1 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2}
            sage: S7 = delta_complexes.Sphere(7)
            sage: S7.homology(7)
            Z
            sage: cubical_complexes.KleinBottle().homology(1, base_ring=GF(2))
            Vector space of dimension 2 over Finite Field of size 2

        If CHomP is installed, Sage can compute generators of homology
        groups::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: S2.homology(dim=2, generators=True, base_ring=GF(2))  # optional - CHomP
            (Vector space of dimension 1 over Finite Field of size 2, [(0, 1, 2) + (0, 1, 3) + (0, 2, 3) + (1, 2, 3)])

        When generators are computed, Sage returns a pair for each
        dimension: the group and the list of generators.  For
        simplicial complexes, each generator is represented as a
        linear combination of simplices, as above, and for cubical
        complexes, each generator is a linear combination of cubes::

            sage: S2_cub = cubical_complexes.Sphere(2)
            sage: S2_cub.homology(dim=2, generators=True)  # optional - CHomP
            (Z, [-[[0,1] x [0,1] x [0,0]] + [[0,1] x [0,1] x [1,1]] - [[0,0] x [0,1] x [0,1]] - [[0,1] x [1,1] x [0,1]] + [[0,1] x [0,0] x [0,1]] + [[1,1] x [0,1] x [0,1]]])
        """
        from sage.interfaces.chomp import have_chomp, homcubes, homsimpl
        from sage.homology.cubical_complex import CubicalComplex
        from sage.homology.simplicial_complex import SimplicialComplex
        from sage.modules.all import VectorSpace
        from sage.homology.homology_group import HomologyGroup

        base_ring = kwds.get('base_ring', ZZ)
        cohomology = kwds.get('cohomology', False)
        subcomplex = kwds.pop('subcomplex', None)
        verbose = kwds.get('verbose', False)
        algorithm = kwds.get('algorithm', 'auto')

        if dim is not None:
            if isinstance(dim, (list, tuple)):
                low = min(dim) - 1
                high = max(dim) + 2
            else:
                low = dim - 1
                high = dim + 2
            dims = range(low, high)
        else:
            dims = None

        # try to use CHomP if computing homology (not cohomology) and
        # working over Z or F_p, p a prime.
        if (algorithm == 'auto' and cohomology is False
            and (base_ring == ZZ or (base_ring.is_prime_field()
                                     and base_ring != QQ))):
            # homcubes, homsimpl seems fastest if all of homology is computed.
            H = None
            if isinstance(self, CubicalComplex):
                if have_chomp('homcubes'):
                    H = homcubes(self, subcomplex, **kwds)
            elif isinstance(self, SimplicialComplex):
                if have_chomp('homsimpl'):
                    H = homsimpl(self, subcomplex, **kwds)

            # now pick off the requested dimensions
            if H:
                answer = {}
                if not dims:
                    dims =range(self.dimension() + 1)
                for d in dims:
                    answer[d] = H.get(d, HomologyGroup(0, base_ring))
                if dim is not None:
                    if not isinstance(dim, (list, tuple)):
                        answer = answer.get(dim, HomologyGroup(0, base_ring))
                return answer

        # Derived classes can implement specialized algorithms using a
        # _homology_ method.  See SimplicialComplex for one example.
        if hasattr(self, '_homology_'):
            return self._homology_(dim, subcomplex=subcomplex, **kwds)

        C = self.chain_complex(cochain=cohomology, augmented=True,
                               dimensions=dims, subcomplex=subcomplex, **kwds)
        answer = C.homology(**kwds)
        if dim is None:
            dim = range(self.dimension()+1)
        zero = HomologyGroup(0, base_ring)
        if isinstance(dim, (list, tuple)):
            return dict([d, answer.get(d, zero)] for d in dim)
        return answer.get(dim, zero)

    def cohomology(self, dim=None, **kwds):
        r"""
        The reduced cohomology of this cell complex.

        The arguments are the same as for the :meth:`homology` method,
        except that :meth:`homology` accepts a ``cohomology`` key
        word, while this function does not: ``cohomology`` is
        automatically true here.  Indeed, this function just calls
        :meth:`homology` with ``cohomology`` set to ``True``.

        :param dim:
        :param base_ring:
        :param subcomplex:
        :param algorithm:
        :param verbose:

        EXAMPLES::

            sage: circle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: circle.cohomology(0)
            0
            sage: circle.cohomology(1)
            Z
            sage: P2 = SimplicialComplex([[0,1,2], [0,2,3], [0,1,5], [0,4,5], [0,3,4], [1,2,4], [1,3,4], [1,3,5], [2,3,5], [2,4,5]])   # projective plane
            sage: P2.cohomology(2)
            C2
            sage: P2.cohomology(2, base_ring=GF(2))
            Vector space of dimension 1 over Finite Field of size 2
            sage: P2.cohomology(2, base_ring=GF(3))
            Vector space of dimension 0 over Finite Field of size 3

            sage: cubical_complexes.KleinBottle().cohomology(2)
            C2

        Relative cohomology::

            sage: T = SimplicialComplex([[0,1]])
            sage: U = SimplicialComplex([[0], [1]])
            sage: T.cohomology(1, subcomplex=U)
            Z

        A `\Delta`-complex example::

            sage: s5 = delta_complexes.Sphere(5)
            sage: s5.cohomology(base_ring=GF(7))[5]
            Vector space of dimension 1 over Finite Field of size 7
        """
        return self.homology(dim=dim, cohomology=True, **kwds)

    def betti(self, dim=None, subcomplex=None):
        r"""
        The Betti numbers of this simplicial complex as a dictionary
        (or a single Betti number, if only one dimension is given):
        the ith Betti number is the rank of the ith homology group.

        :param dim: If ``None``, then return every Betti number, as
           a dictionary with keys the non-negative integers.  If
           ``dim`` is an integer or list, return the Betti number for
           each given dimension.  (Actually, if ``dim`` is a list,
           return the Betti numbers, as a dictionary, in the range
           from ``min(dim)`` to ``max(dim)``.  If ``dim`` is a number,
           return the Betti number in that dimension.)
        :type dim: integer or list of integers or ``None``; optional,
           default ``None``
        :param subcomplex: a subcomplex of this cell complex.  Compute
           the Betti numbers of the homology relative to this subcomplex.
        :type subcomplex: optional, default ``None``

        EXAMPLES:

        Build the two-sphere as a three-fold join of a
        two-point space with itself::

            sage: S = SimplicialComplex([[0], [1]])
            sage: (S*S*S).betti()
            {0: 1, 1: 0, 2: 1}
            sage: (S*S*S).betti([1,2])
            {1: 0, 2: 1}
            sage: (S*S*S).betti(2)
            1

        Or build the two-sphere as a `\Delta`-complex::

            sage: S2 = delta_complexes.Sphere(2)
            sage: S2.betti([1,2])
            {1: 0, 2: 1}

        Or as a cubical complex::

            sage: S2c = cubical_complexes.Sphere(2)
            sage: S2c.betti(2)
            1
        """
        dict = {}
        H = self.homology(dim, base_ring=QQ, subcomplex=subcomplex)
        try:
            for n in H.keys():
                dict[n] = H[n].dimension()
                if n == 0:
                    dict[n] += 1
            return dict
        except AttributeError:
            return H.dimension()

    ############################################################
    # end of chain complexes, homology
    ############################################################

    def face_poset(self):
        r"""
        The face poset of this cell complex, the poset of
        nonempty cells, ordered by inclusion.

        This uses the :meth:`cells` method, and also assumes that for
        each cell ``f``, all of ``f.faces()``, ``tuple(f)``, and
        ``f.dimension()`` make sense.  (If this is not the case in
        some derived class, as happens with `\Delta`-complexes, then
        override this method.)

        EXAMPLES::

            sage: P = SimplicialComplex([[0, 1], [1,2], [2,3]]).face_poset(); P
            Finite poset containing 7 elements
            sage: P.list()
            [(3,), (2,), (2, 3), (1,), (0,), (0, 1), (1, 2)]

            sage: S2 = cubical_complexes.Sphere(2)
            sage: S2.face_poset()
            Finite poset containing 26 elements
        """
        from sage.combinat.posets.posets import Poset
        from sage.misc.flatten import flatten
        covers = {}
        # The code for posets seems to work better if each cell is
        # converted to a tuple.
        all_cells = flatten([list(f) for f in self.cells().values()])

        for C in all_cells:
            if C.dimension() >= 0:
                covers[tuple(C)] = []
        for C in all_cells:
            for face in C.faces():
                if face.dimension() >= 0:
                    covers[tuple(face)].append(tuple(C))
        return Poset(covers)

    def graph(self):
        """
        The 1-skeleton of this cell complex, as a graph.

        This is not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.graph()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def n_skeleton(self, n):
        """
        The `n`-skeleton of this cell complex: the cell
        complex obtained by discarding all of the simplices in
        dimensions larger than `n`.

        :param n: non-negative integer

        This is not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.n_skeleton(3)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def category(self):
        """
        Return the category to which this chain complex belongs: the
        category of all cell complexes.

        This is not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.category()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _string_constants(self):
        """
        Tuple containing the name of the type of complex, and the
        singular and plural of the name of the cells from which it is
        built.  This is used in constructing the string representation.

        :return: tuple of strings

        This returns ``('Cell', 'cell', 'cells')``, as in "Cell
        complex", "1 cell", and "24 cells", but in other classes it
        could be overridden, as for example with ``('Cubical', 'cube',
        'cubes')`` or ``('Delta', 'simplex', 'simplices')``.  If for a
        derived class, the basic form of the print representation is
        acceptable, you can just modify these strings.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: GenericCellComplex()._string_constants()
            ('Cell', 'cell', 'cells')
            sage: delta_complexes.Sphere(0)._string_constants()
            ('Delta', 'simplex', 'simplices')
            sage: cubical_complexes.Sphere(0)._string_constants()
            ('Cubical', 'cube', 'cubes')
        """
        return ('Cell', 'cell', 'cells')

    def _repr_(self):
        """
        Print representation.

        :return: string

        EXAMPLES::

            sage: delta_complexes.Sphere(7) # indirect doctest
            Delta complex with 8 vertices and 257 simplices
            sage: delta_complexes.Torus()._repr_()
            'Delta complex with 1 vertex and 7 simplices'
        """
        vertices = len(self.n_cells(0))
        Name, cell_name, cells_name = self._string_constants()
        if vertices != 1:
            vertex_string = "with %s vertices" % vertices
        else:
            vertex_string = "with 1 vertex"
        cells = 0
        for dim in self.cells():
            cells += len(self.cells()[dim])
        if cells != 1:
            cells_string = " and %s %s" % (cells, cells_name)
        else:
            cells_string = " and 1 %s" % cell_name
        return Name + " complex " + vertex_string + cells_string
