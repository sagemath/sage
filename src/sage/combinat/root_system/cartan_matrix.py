# sage_setup: distribution = sagemath-modules
"""
Cartan matrices

AUTHORS:

- Travis Scrimshaw (2012-04-22): Nicolas M. Thiery moved matrix creation to
  :class:`CartanType` to prepare :func:`cartan_matrix()` for deprecation.
- Christian Stump, Travis Scrimshaw (2013-04-13): Created :class:`CartanMatrix`.
- Ben Salisbury (2018-08-07): Added Borcherds-Cartan matrices.
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2012,2013 Travis Scrimshaw <tscrim at ucdavis.edu>,
#       Copyright (C) 2013 Christian Stump,
#       Copyright (C) 2018 Ben Salisbury <salis1bt at cmich.edu>,
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
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix
from sage.misc.lazy_import import lazy_import
from sage.structure.element import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.classcall_metaclass import typecall
from sage.combinat.subset import powerset
from sage.rings.integer_ring import ZZ
from sage.combinat.root_system.cartan_type import CartanType, CartanType_abstract
from sage.combinat.root_system.root_system import RootSystem
from sage.sets.family import Family

lazy_import('sage.graphs.digraph', 'DiGraph')
lazy_import('sage.combinat.root_system.dynkin_diagram', 'DynkinDiagram_class')


try:
    from sage.matrix.matrix_integer_sparse import Matrix_integer_sparse as Base
except ImportError:
    from sage.matrix.matrix_generic_sparse import Matrix_generic_sparse as Base


class CartanMatrix(Base, CartanType_abstract,
                   metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A (generalized) Cartan matrix.

    A matrix `A = (a_{ij})_{i,j \in I}` for some index set `I` is a
    generalized Cartan matrix if it satisfies the following properties:

    - `a_{ii} = 2` for all `i`,
    - `a_{ij} \leq 0` for all `i \neq j`,
    - `a_{ij} = 0` if and only if `a_{ji} = 0` for all `i \neq j`.

    Additionally some reference assume that a Cartan matrix is
    *symmetrizable* (see :meth:`is_symmetrizable`). However following Kac, we
    do not make that assumption here.

    An even, integral Borcherds--Cartan matrix is an integral matrix
    `A = (a_{ij})_{i,j \in I}` for some countable index set `I` which satisfies
    the following properties:

    - `a_{ii} \in \{2\} \cup 2\ZZ_{<0}` for all `i`,
    - `a_{ij} \leq 0` for all `i \neq j`,
    - `a_{ij} = 0` if and only if `a_{ji} = 0` for all `i \neq j`.

    INPUT:

    Can be anything which is accepted by ``CartanType`` or a matrix.

    If given a matrix, one can also use the keyword ``cartan_type`` when giving
    a matrix to explicitly state the type. Otherwise this will try to check the
    input matrix against possible standard types of Cartan matrices. To disable
    this check, use the keyword ``cartan_type_check = False``.

    If one wants to initialize a Borcherds-Cartan matrix using matrix data,
    use the keyword ``borcherds=True``. To specify the diagonal entries of
    corresponding to a Cartan type (a Cartan matrix is treated as matrix data),
    use ``borcherds`` with a list of the diagonal entries.

    EXAMPLES::

        sage: # needs sage.graphs
        sage: CartanMatrix(['A', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
        sage: CartanMatrix(['B', 6])
        [ 2 -1  0  0  0  0]
        [-1  2 -1  0  0  0]
        [ 0 -1  2 -1  0  0]
        [ 0  0 -1  2 -1  0]
        [ 0  0  0 -1  2 -1]
        [ 0  0  0  0 -2  2]
        sage: CartanMatrix(['C', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -2]
        [ 0  0 -1  2]
        sage: CartanMatrix(['D', 6])
        [ 2 -1  0  0  0  0]
        [-1  2 -1  0  0  0]
        [ 0 -1  2 -1  0  0]
        [ 0  0 -1  2 -1 -1]
        [ 0  0  0 -1  2  0]
        [ 0  0  0 -1  0  2]
        sage: CartanMatrix(['E',6])
        [ 2  0 -1  0  0  0]
        [ 0  2  0 -1  0  0]
        [-1  0  2 -1  0  0]
        [ 0 -1 -1  2 -1  0]
        [ 0  0  0 -1  2 -1]
        [ 0  0  0  0 -1  2]
        sage: CartanMatrix(['E',7])
        [ 2  0 -1  0  0  0  0]
        [ 0  2  0 -1  0  0  0]
        [-1  0  2 -1  0  0  0]
        [ 0 -1 -1  2 -1  0  0]
        [ 0  0  0 -1  2 -1  0]
        [ 0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0 -1  2]
        sage: CartanMatrix(['E', 8])
        [ 2  0 -1  0  0  0  0  0]
        [ 0  2  0 -1  0  0  0  0]
        [-1  0  2 -1  0  0  0  0]
        [ 0 -1 -1  2 -1  0  0  0]
        [ 0  0  0 -1  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0 -1  2]
        sage: CartanMatrix(['F', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -2  2 -1]
        [ 0  0 -1  2]

    This is different from MuPAD-Combinat, due to different node
    convention?

    ::

        sage: # needs sage.graphs
        sage: CartanMatrix(['G', 2])
        [ 2 -3]
        [-1  2]
        sage: CartanMatrix(['A',1,1])
        [ 2 -2]
        [-2  2]
        sage: CartanMatrix(['A', 3, 1])
        [ 2 -1  0 -1]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [-1  0 -1  2]
        sage: CartanMatrix(['B', 3, 1])
        [ 2  0 -1  0]
        [ 0  2 -1  0]
        [-1 -1  2 -1]
        [ 0  0 -2  2]
        sage: CartanMatrix(['C', 3, 1])
        [ 2 -1  0  0]
        [-2  2 -1  0]
        [ 0 -1  2 -2]
        [ 0  0 -1  2]
        sage: CartanMatrix(['D', 4, 1])
        [ 2  0 -1  0  0]
        [ 0  2 -1  0  0]
        [-1 -1  2 -1 -1]
        [ 0  0 -1  2  0]
        [ 0  0 -1  0  2]
        sage: CartanMatrix(['E', 6, 1])
        [ 2  0 -1  0  0  0  0]
        [ 0  2  0 -1  0  0  0]
        [-1  0  2  0 -1  0  0]
        [ 0 -1  0  2 -1  0  0]
        [ 0  0 -1 -1  2 -1  0]
        [ 0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0 -1  2]
        sage: CartanMatrix(['E', 7, 1])
        [ 2 -1  0  0  0  0  0  0]
        [-1  2  0 -1  0  0  0  0]
        [ 0  0  2  0 -1  0  0  0]
        [ 0 -1  0  2 -1  0  0  0]
        [ 0  0 -1 -1  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0 -1  2]
        sage: CartanMatrix(['E', 8, 1])
        [ 2  0  0  0  0  0  0  0 -1]
        [ 0  2  0 -1  0  0  0  0  0]
        [ 0  0  2  0 -1  0  0  0  0]
        [ 0 -1  0  2 -1  0  0  0  0]
        [ 0  0 -1 -1  2 -1  0  0  0]
        [ 0  0  0  0 -1  2 -1  0  0]
        [ 0  0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0  0 -1  2 -1]
        [-1  0  0  0  0  0  0 -1  2]
        sage: CartanMatrix(['F', 4, 1])
        [ 2 -1  0  0  0]
        [-1  2 -1  0  0]
        [ 0 -1  2 -1  0]
        [ 0  0 -2  2 -1]
        [ 0  0  0 -1  2]
        sage: CartanMatrix(['G', 2, 1])
        [ 2  0 -1]
        [ 0  2 -3]
        [-1 -1  2]

    Examples of Borcherds-Cartan matrices::

        sage: CartanMatrix([[2,-1],[-1,-2]], borcherds=True)                            # needs sage.graphs
        [ 2 -1]
        [-1 -2]
        sage: CartanMatrix('B3', borcherds=[-4,-6,2])                                   # needs sage.graphs
        [-4 -1  0]
        [-1 -6 -1]
        [ 0 -2  2]

    .. NOTE::

        Since this is a matrix, :meth:`row()` and :meth:`column()` will return
        the standard row and column respectively. To get the row with the
        indices as in Dynkin diagrams/Cartan types, use
        :meth:`row_with_indices()` and :meth:`column_with_indices()`
        respectively.
    """
    @staticmethod
    def __classcall_private__(cls, data=None, index_set=None,
                              cartan_type=None, cartan_type_check=True,
                              borcherds=None):
        """
        Normalize input so we can inherit from sparse integer matrix.

        .. NOTE::

            To disable the Cartan type check, use the optional argument
            ``cartan_type_check = False``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: C = CartanMatrix(['A',1,1])
            sage: C2 = CartanMatrix([[2, -2], [-2, 2]])
            sage: C3 = CartanMatrix(matrix([[2, -2], [-2, 2]]), [0, 1])
            sage: C == C2 and C == C3
            True

        TESTS:

        Check that :issue:`15740` is fixed::

            sage: # needs sage.graphs
            sage: d = DynkinDiagram()
            sage: d.add_edge('a', 'b', 2)
            sage: d.index_set()
            ('a', 'b')
            sage: cm = CartanMatrix(d)
            sage: cm.index_set()
            ('a', 'b')
        """
        # Special case with 0 args and kwds has Cartan type
        if cartan_type is not None and data is None:
            data = CartanType(cartan_type)

        if data is None:
            data = []
            n = 0
            index_set = tuple()
            cartan_type = None
            subdivisions = None
        elif isinstance(data, CartanMatrix):
            if index_set is not None:
                d = {a: index_set[i] for i,a in enumerate(data.index_set())}
                return data.relabel(d)
            return data
        else:
            dynkin_diagram = None
            subdivisions = None

            if isinstance(data, DynkinDiagram_class):
                dynkin_diagram = data
                cartan_type = data._cartan_type
            else:
                try:
                    cartan_type = CartanType(data)
                    dynkin_diagram = cartan_type.dynkin_diagram()
                except (TypeError, ValueError):
                    pass

            if dynkin_diagram is not None:
                n = dynkin_diagram.rank()
                index_set = dynkin_diagram.index_set()
                oir = dynkin_diagram.odd_isotropic_roots()
                reverse = {a: i for i,a in enumerate(index_set)}
                if isinstance(borcherds, (list, tuple)):
                    if (len(borcherds) != len(index_set)
                        and not all(val in ZZ
                                    and (val == 2 or (val % 2 == 0 and val < 0))
                                    for val in borcherds)):
                        raise ValueError("the input data is not a Borcherds-Cartan matrix")
                    data = {(i, i): val if index_set[i] not in oir else 0
                            for i,val in enumerate(borcherds)}
                else:
                    data = {(i, i): 2 if index_set[i] not in oir else 0
                            for i in range(n)}
                for (i,j,l) in dynkin_diagram.edge_iterator():
                    data[(reverse[j], reverse[i])] = -l
            else:
                M = matrix(data)
                if borcherds:
                    if not is_borcherds_cartan_matrix(M):
                        raise ValueError("the input matrix is not a Borcherds-Cartan matrix")
                else:
                    if not is_generalized_cartan_matrix(M):
                        raise ValueError("the input matrix is not a generalized Cartan matrix")
                n = M.ncols()
                data = M.dict()
                subdivisions = M._subdivisions

            if index_set is None:
                index_set = tuple(range(n))
            else:
                index_set = tuple(index_set)

        if len(index_set) != n and len(set(index_set)) != n:
            raise ValueError("the given index set is not valid")

        # We can do the Cartan type initialization later as this is not
        #   a unique representation
        mat = typecall(cls, MatrixSpace(ZZ, n, sparse=True), data, False, True)
        # FIXME: We have to initialize the CartanMatrix part separately because
        #   of the __cinit__ of the matrix. We should get rid of this workaround
        mat._CM_init(cartan_type, index_set, cartan_type_check)
        mat._subdivisions = subdivisions
        return mat

    def matrix_space(self, nrows=None, ncols=None, sparse=None):
        r"""
        Return a matrix space over the integers.

        INPUT:

        - ``nrows`` -- number of rows

        - ``ncols`` -- number of columns

        - ``sparse`` -- boolean

        EXAMPLES::

            sage: # needs sage.graphs
            sage: cm = CartanMatrix(['A', 3])
            sage: cm.matrix_space()
            Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
            sage: cm.matrix_space(2, 2)
            Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
            sage: cm[:2,1:]   # indirect doctest
            [-1  0]
            [ 2 -1]
        """
        if nrows is None:
            nrows = self.nrows()
        if ncols is None:
            ncols = self.ncols()
        if sparse is None:
            sparse = True

        if nrows == self.nrows() and ncols == self.ncols() and sparse:
            return self.parent()
        else:
            from sage.matrix.matrix_space import MatrixSpace
            return MatrixSpace(ZZ, nrows, ncols, sparse is None or bool(sparse))

    def _CM_init(self, cartan_type, index_set, cartan_type_check):
        """
        Initialize ``self`` as a Cartan matrix.

        TESTS::

            sage: C = CartanMatrix(['A',1,1])  # indirect doctest                       # needs sage.graphs
            sage: TestSuite(C).run(skip=["_test_category", "_test_change_ring"])        # needs sage.graphs

        Check that :issue:`37979` is fixed::

            sage: C = CartanMatrix([[2]], index_set=(4,))
            sage: C.index_set()
            (4,)
            sage: CartanType("A3").subtype((2,)) is CartanType("A1").relabel({1:2})
            True
        """
        self._index_set = index_set
        self.set_immutable()

        if cartan_type is not None:
            cartan_type = CartanType(cartan_type)
        elif self.nrows() == 1:
            cartan_type = CartanType(['A', 1])
            if index_set != (1,):
                cartan_type = cartan_type.relabel({1: index_set[0]})
        elif cartan_type_check:
            # Placeholder so we don't have to reimplement creating a
            #   Dynkin diagram from a Cartan matrix
            self._cartan_type = None
            cartan_type = find_cartan_type_from_matrix(self)

        self._cartan_type = cartan_type

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: CM = CartanMatrix(['A',4])                                            # needs sage.graphs
            sage: x = loads(dumps(CM))                                                  # needs sage.graphs
            sage: x._index_set                                                          # needs sage.graphs
            (1, 2, 3, 4)
        """
        if self._cartan_type:
            return (CartanMatrix, (self._cartan_type,))
        return (CartanMatrix, (self.dynkin_diagram(),))

    def root_system(self):
        """
        Return the root system corresponding to ``self``.

        EXAMPLES::

            sage: C = CartanMatrix(['A',3])                                             # needs sage.graphs
            sage: C.root_system()                                                       # needs sage.graphs
            Root system of type ['A', 3]
        """
        if self._cartan_type is not None:
            return RootSystem(self._cartan_type)
        return self.dynkin_diagram().root_system()

    def root_space(self):
        """
        Return the root space corresponding to ``self``.

        EXAMPLES::

            sage: C = CartanMatrix(['A',3])                                             # needs sage.graphs
            sage: C.root_space()                                                        # needs sage.graphs
            Root space over the Rational Field of the Root system of type ['A', 3]
        """
        return self.root_system().root_space()

    def reflection_group(self, type='matrix'):
        """
        Return the reflection group corresponding to ``self``.

        EXAMPLES::

            sage: C = CartanMatrix(['A',3])                                             # needs sage.graphs
            sage: C.reflection_group()                                                  # needs sage.graphs sage.libs.gap
            Weyl Group of type ['A', 3] (as a matrix group acting on the root space)
        """
        RS = self.root_space()

        if type == "matrix":
            return RS.weyl_group()

        if type == "permutation":
            if not self.is_finite():
                raise ValueError("only works for finite types")
            Phi = RS.roots()
            gens = {}
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            S = SymmetricGroup(len(Phi))
            for i in self.index_set():
                pi = S([ Phi.index( beta.simple_reflection(i) ) + 1 for beta in Phi ])
                gens[i] = pi
            return S.subgroup( gens[i] for i in gens )

        raise ValueError("the reflection group is only available as a matrix group or as a permutation group")

    def symmetrizer(self):
        """
        Return the symmetrizer of ``self``.

        EXAMPLES::

            sage: cm = CartanMatrix([[2,-5],[-2,2]])                                    # needs sage.graphs
            sage: cm.symmetrizer()                                                      # needs sage.graphs
            Finite family {0: 2, 1: 5}

        TESTS:

        Check that the symmetrizer computed from the Cartan matrix agrees
        with the values given by the Cartan type::

            sage: ct = CartanType(['B',4,1])
            sage: ct.symmetrizer()                                                      # needs sage.graphs
            Finite family {0: 2, 1: 2, 2: 2, 3: 2, 4: 1}
            sage: ct.cartan_matrix().symmetrizer()                                      # needs sage.graphs
            Finite family {0: 2, 1: 2, 2: 2, 3: 2, 4: 1}
        """
        sym = self.is_symmetrizable(True)
        if not sym:
            raise ValueError("the Cartan matrix is not symmetrizable")
        iset = self.index_set()
        # The result from is_symmetrizable needs to be scaled
        # to integer coefficients
        from sage.arith.functions import lcm as LCM
        from sage.rings.rational_field import QQ
        scalar = LCM([QQ(x).denominator() for x in sym])
        return Family( {iset[i]: ZZ(val*scalar) for i, val in enumerate(sym)} )

    @cached_method
    def symmetrized_matrix(self):
        """
        Return the symmetrized matrix of ``self`` if symmetrizable.

        EXAMPLES::

            sage: cm = CartanMatrix(['B',4,1])                                          # needs sage.graphs
            sage: cm.symmetrized_matrix()                                               # needs sage.graphs
            [ 4  0 -2  0  0]
            [ 0  4 -2  0  0]
            [-2 -2  4 -2  0]
            [ 0  0 -2  4 -2]
            [ 0  0  0 -2  2]
        """
        M = matrix.diagonal(list(self.symmetrizer())) * self
        M.set_immutable()
        return M

    ##########################################################################
    # Cartan type methods

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: C = CartanMatrix(['A',1,1])
            sage: C.index_set()
            (0, 1)
            sage: C = CartanMatrix(['E',6])
            sage: C.index_set()
            (1, 2, 3, 4, 5, 6)
        """
        return self._index_set

    def cartan_type(self):
        """
        Return the Cartan type of ``self`` or ``self`` if unknown.

        EXAMPLES::

            sage: C = CartanMatrix(['A',4,1])                                           # needs sage.graphs
            sage: C.cartan_type()                                                       # needs sage.graphs
            ['A', 4, 1]

        If the Cartan type is unknown::

            sage: C = CartanMatrix([[2,-1,-2], [-1,2,-1], [-2,-1,2]])                   # needs sage.graphs
            sage: C.cartan_type()                                                       # needs sage.graphs
            [ 2 -1 -2]
            [-1  2 -1]
            [-2 -1  2]
        """
        if self._cartan_type is None:
            return self
        if is_borcherds_cartan_matrix(self) and not is_generalized_cartan_matrix(self):
            return self
        return self._cartan_type

    def subtype(self, index_set):
        """
        Return a subtype of ``self`` given by ``index_set``.

        A subtype can be considered the Dynkin diagram induced from
        the Dynkin diagram of ``self`` by ``index_set``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: C = CartanMatrix(['F',4])
            sage: S = C.subtype([1,2,3])
            sage: S
            [ 2 -1  0]
            [-1  2 -1]
            [ 0 -2  2]
            sage: S.index_set()
            (1, 2, 3)
        """
        ind = self.index_set()
        I = [ind.index(i) for i in index_set]
        return CartanMatrix(self.matrix_from_rows_and_columns(I, I), index_set=index_set)

    def rank(self):
        r"""
        Return the rank of ``self``.

        EXAMPLES::

            sage: CartanMatrix(['C',3]).rank()                                          # needs sage.graphs
            3
            sage: CartanMatrix(["A2","B2","F4"]).rank()                                 # needs sage.graphs
            8
        """
        return self.ncols()

    def relabel(self, relabelling):
        """
        Return the relabelled Cartan matrix.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: CM = CartanMatrix(['C',3])
            sage: R = CM.relabel({1:0, 2:4, 3:1}); R
            [ 2  0 -1]
            [ 0  2 -1]
            [-1 -2  2]
            sage: R.index_set()
            (0, 1, 4)
            sage: CM
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
        """
        return self.dynkin_diagram().relabel(relabelling, inplace=False).cartan_matrix()

    @cached_method
    def dynkin_diagram(self):
        """
        Return the Dynkin diagram corresponding to ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: C = CartanMatrix(['A',2])
            sage: C.dynkin_diagram()
            O---O
            1   2
            A2
            sage: C = CartanMatrix(['F',4,1])
            sage: C.dynkin_diagram()
            O---O---O=>=O---O
            0   1   2   3   4
            F4~
            sage: C = CartanMatrix([[2,-4],[-4,2]])
            sage: C.dynkin_diagram()
            Dynkin diagram of rank 2
        """
        from sage.combinat.root_system.dynkin_diagram import DynkinDiagram
        if self._cartan_type is not None:
            return DynkinDiagram(self._cartan_type)
        return DynkinDiagram(self)

    def cartan_matrix(self):
        r"""
        Return the Cartan matrix of ``self``.

        EXAMPLES::

            sage: CartanMatrix(['C',3]).cartan_matrix()                                 # needs sage.graphs
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
        """
        return self

    def dual(self):
        r"""
        Return the dual Cartan matrix of ``self``, which is obtained by taking
        the transpose.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: ct = CartanType(['C',3])
            sage: M = CartanMatrix(ct); M
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
            sage: M.dual()
            [ 2 -1  0]
            [-1  2 -1]
            [ 0 -2  2]
            sage: M.dual() == CartanMatrix(ct.dual())
            True
            sage: M.dual().cartan_type() == ct.dual()
            True

        An example with arbitrary Cartan matrices::

            sage: # needs sage.graphs
            sage: cm = CartanMatrix([[2,-5], [-2, 2]]); cm
            [ 2 -5]
            [-2  2]
            sage: cm.dual()
            [ 2 -2]
            [-5  2]
            sage: cm.dual() == CartanMatrix(cm.transpose())
            True
            sage: cm.dual().dual() == cm
            True
        """
        if self._cartan_type is not None:
            return CartanMatrix(self._cartan_type.dual())
        return CartanMatrix(self.transpose())

    def is_simply_laced(self):
        """
        Implement :meth:`CartanType_abstract.is_simply_laced()`.

        A Cartan matrix is simply-laced if all non diagonal entries are `0`
        or `-1`.

        EXAMPLES::

            sage: cm = CartanMatrix([[2, -1, -1, -1], [-1, 2, -1, -1],                  # needs sage.graphs
            ....:                    [-1, -1, 2, -1], [-1, -1, -1, 2]])
            sage: cm.is_simply_laced()                                                  # needs sage.graphs
            True
        """
        for i in range(self.nrows()):
            for j in range(i+1, self.ncols()):
                if self[i, j] < -1 or self[j, i] < -1:
                    return False
        return True

    def is_crystallographic(self):
        """
        Implement :meth:`CartanType_abstract.is_crystallographic`.

        A Cartan matrix is crystallographic if it is symmetrizable.

        EXAMPLES::

            sage: CartanMatrix(['F',4]).is_crystallographic()                           # needs sage.graphs
            True
        """
        return self.is_symmetrizable()

    def column_with_indices(self, j):
        """
        Return the `j`-th column `(a_{i,j})_i` of ``self`` as a container
        (or iterator) of tuples `(i, a_{i,j})`

        EXAMPLES::

            sage: M = CartanMatrix(['B',4])                                             # needs sage.graphs
            sage: [ (i,a) for (i,a) in M.column_with_indices(3) ]                       # needs sage.graphs
            [(3, 2), (2, -1), (4, -2)]
        """
        return self.dynkin_diagram().column(j)

    def row_with_indices(self, i):
        """
        Return the `i`-th row `(a_{i,j})_j` of ``self`` as a container
        (or iterator) of tuples `(j, a_{i,j})`

        EXAMPLES::

            sage: M = CartanMatrix(['C',4])                                             # needs sage.graphs
            sage: [ (i,a) for (i,a) in M.row_with_indices(3) ]                          # needs sage.graphs
            [(3, 2), (2, -1), (4, -2)]
        """
        return self.dynkin_diagram().row(i)

    @cached_method
    def is_finite(self):
        """
        Return ``True`` if ``self`` is a finite type or ``False`` otherwise.

        A generalized Cartan matrix is finite if the determinant of all its
        principal submatrices (see :meth:`principal_submatrices`) is positive.
        Such matrices have a positive definite symmetrized matrix. Note that a
        finite matrix may consist of multiple blocks of Cartan matrices each
        having finite Cartan type.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = CartanMatrix(['C',4])
            sage: M.is_finite()
            True
            sage: M = CartanMatrix(['D',4,1])
            sage: M.is_finite()
            False
            sage: M = CartanMatrix([[2, -4], [-3, 2]])
            sage: M.is_finite()
            False
        """
        if self._cartan_type is None:
            if not self.is_symmetrizable():
                return False
            return self.symmetrized_matrix().is_positive_definite()
        return self._cartan_type.is_finite()

    @cached_method
    def is_affine(self):
        """
        Return ``True`` if ``self`` is an affine type or ``False`` otherwise.

        A generalized Cartan matrix is affine if all of its indecomposable
        blocks are either finite (see :meth:`is_finite`) or have zero
        determinant with all proper principal minors positive.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = CartanMatrix(['C',4])
            sage: M.is_affine()
            False
            sage: M = CartanMatrix(['D',4,1])
            sage: M.is_affine()
            True
            sage: M = CartanMatrix([[2, -4], [-3, 2]])
            sage: M.is_affine()
            False
        """
        if self._cartan_type is None:
            if self.det() != 0:
                return False
            for b in self.indecomposable_blocks():
                if b.det() < 0 or not all(
                    a.det() > 0 for a in b.principal_submatrices(proper=True)):
                    return False
            return True
        return self._cartan_type.is_affine()

    @cached_method
    def is_hyperbolic(self, compact=False):
        """
        Return if ``True`` if ``self`` is a (compact) hyperbolic type
        or ``False`` otherwise.

        An indecomposable generalized Cartan matrix is hyperbolic if it has
        negative determinant and if any proper connected subdiagram of its
        Dynkin diagram is of finite or affine type. It is compact hyperbolic
        if any proper connected subdiagram has finite type.

        INPUT:

        - ``compact`` -- if ``True``, check if matrix is compact hyperbolic

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = CartanMatrix([[2,-2,0],[-2,2,-1],[0,-1,2]])
            sage: M.is_hyperbolic()
            True
            sage: M.is_hyperbolic(compact=True)
            False
            sage: M = CartanMatrix([[2,-3],[-3,2]])
            sage: M.is_hyperbolic()
            True
            sage: M = CartanMatrix(['C',4])
            sage: M.is_hyperbolic()
            False
        """
        if not self.is_indefinite() or not self.is_indecomposable():
            return False

        D = self.dynkin_diagram()
        verts = tuple(D.vertex_iterator())
        for v in verts:
            l = set(verts)-set((v,))
            subg = D.subgraph(vertices=l)
            if compact and not subg.is_finite():
                return False
            elif not subg.is_finite() and not subg.is_affine():
                return False
        return True

    @cached_method
    def is_lorentzian(self):
        """
        Return ``True`` if ``self`` is a Lorentzian type or ``False`` otherwise.

        A generalized Cartan matrix is Lorentzian if it has negative determinant
        and exactly one negative eigenvalue.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = CartanMatrix([[2,-3],[-3,2]])
            sage: M.is_lorentzian()
            True
            sage: M = CartanMatrix([[2,-1],[-1,2]])
            sage: M.is_lorentzian()
            False
        """
        if self.det() >= 0:
            return False
        return sum(1 for x in self.eigenvalues() if x < 0) == 1

    @cached_method
    def is_indefinite(self):
        """
        Return if ``self`` is an indefinite type or ``False`` otherwise.

        EXAMPLES::

           sage: # needs sage.graphs
           sage: M = CartanMatrix([[2,-3],[-3,2]])
           sage: M.is_indefinite()
           True
           sage: M = CartanMatrix("A2")
           sage: M.is_indefinite()
           False
        """
        return not self.is_finite() and not self.is_affine()

    @cached_method
    def is_indecomposable(self):
        """
        Return if ``self`` is an indecomposable matrix or ``False`` otherwise.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = CartanMatrix(['A',5])
            sage: M.is_indecomposable()
            True
            sage: M = CartanMatrix([[2,-1,0],[-1,2,0],[0,0,2]])
            sage: M.is_indecomposable()
            False
        """
        comp_num = self.dynkin_diagram().connected_components_number()
        # consider the empty matrix to be indecomposable
        return comp_num <= 1

    @cached_method
    def coxeter_matrix(self):
        r"""
        Return the Coxeter matrix for ``self``.

        .. SEEALSO:: :meth:`CartanType_abstract.coxeter_matrix`

        EXAMPLES::

            sage: # needs sage.graphs
            sage: cm = CartanMatrix([[2,-5,0],[-2,2,-1],[0,-1,2]])
            sage: cm.coxeter_matrix()
            [ 1 -1  2]
            [-1  1  3]
            [ 2  3  1]
            sage: ct = CartanType([['A',2,2], ['B',3]])
            sage: ct.coxeter_matrix()
            [ 1 -1  2  2  2]
            [-1  1  2  2  2]
            [ 2  2  1  3  2]
            [ 2  2  3  1  4]
            [ 2  2  2  4  1]
            sage: ct.cartan_matrix().coxeter_matrix() == ct.coxeter_matrix()
            True
        """
        scalarproducts_to_order = {0: 2, 1: 3, 2: 4, 3: 6}
        from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
        I = self.index_set()
        n = len(I)
        M = matrix.identity(ZZ, n)
        for i in range(n):
            for j in range(i+1,n):
                val = self[i,j] * self[j,i]
                val = scalarproducts_to_order.get(val, -1)
                M[i,j] = val
                M[j,i] = val
        return CoxeterMatrix(M, index_set=self.index_set(), cartan_type=self)

    @cached_method
    def coxeter_diagram(self):
        r"""
        Construct the Coxeter diagram of ``self``.

        .. SEEALSO:: :meth:`CartanType_abstract.coxeter_diagram`

        EXAMPLES::

            sage: # needs sage.graphs
            sage: cm = CartanMatrix([[2,-5,0],[-2,2,-1],[0,-1,2]])
            sage: G = cm.coxeter_diagram(); G
            Graph on 3 vertices
            sage: G.edges(sort=True)
            [(0, 1, +Infinity), (1, 2, 3)]
            sage: ct = CartanType([['A',2,2], ['B',3]])
            sage: ct.coxeter_diagram()
            Graph on 5 vertices
            sage: ct.cartan_matrix().coxeter_diagram() == ct.coxeter_diagram()
            True
        """
        return self.coxeter_matrix().coxeter_graph()

    def principal_submatrices(self, proper=False):
        """
        Return a list of all principal submatrices of ``self``.

        INPUT:

        - ``proper`` -- if ``True``, return only proper submatrices

        EXAMPLES::

            sage: M = CartanMatrix(['A',2])                                             # needs sage.graphs
            sage: M.principal_submatrices()                                             # needs sage.graphs
            [
                          [ 2 -1]
            [], [2], [2], [-1  2]
            ]
            sage: M.principal_submatrices(proper=True)                                  # needs sage.graphs
            [[], [2], [2]]
        """
        iset = list(range(self.ncols()))
        ret = []
        for l in powerset(iset):
            if not proper or (proper and l != iset):
                ret.append(self.matrix_from_rows_and_columns(l,l))
        return ret

    @cached_method
    def indecomposable_blocks(self):
        """
        Return a tuple of all indecomposable blocks of ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = CartanMatrix(['A',2])
            sage: M.indecomposable_blocks()
            (
            [ 2 -1]
            [-1  2]
            )
            sage: M = CartanMatrix([['A',2,1],['A',3,1]])
            sage: M.indecomposable_blocks()
            (
            [ 2 -1  0 -1]
            [-1  2 -1  0]  [ 2 -1 -1]
            [ 0 -1  2 -1]  [-1  2 -1]
            [-1  0 -1  2], [-1 -1  2]
            )
        """
        subgraphs = self.dynkin_diagram().connected_components_subgraphs()
        return tuple(CartanMatrix(subg._matrix_().rows()) for subg in subgraphs)


def is_borcherds_cartan_matrix(M):
    """
    Return ``True`` if ``M`` is an even, integral Borcherds-Cartan matrix.
    For a definition of such a matrix, see :class:`CartanMatrix`.

    EXAMPLES::

        sage: from sage.combinat.root_system.cartan_matrix import is_borcherds_cartan_matrix
        sage: M = Matrix([[2,-1],[-1,2]])
        sage: is_borcherds_cartan_matrix(M)
        True
        sage: N = Matrix([[2,-1],[-1,0]])
        sage: is_borcherds_cartan_matrix(N)
        False
        sage: O = Matrix([[2,-1],[-1,-2]])
        sage: is_borcherds_cartan_matrix(O)
        True
        sage: O = Matrix([[2,-1],[-1,-3]])
        sage: is_borcherds_cartan_matrix(O)
        False
    """
    if not isinstance(M, Matrix):
        return False
    if not M.is_square():
        return False
    n = M.ncols()
    for i in range(n):
        if M[i,i] == 0:
            return False
        if M[i,i] % 2 == 1:
            return False
        for j in range(i+1, n):
            if M[i,j] > 0 or M[j,i] > 0:
                return False
            elif M[i,j] == 0 and M[j,i] != 0:
                return False
            elif M[j,i] == 0 and M[i,j] != 0:
                return False
    return True


def is_generalized_cartan_matrix(M):
    """
    Return ``True`` if ``M`` is a generalized Cartan matrix. For a definition
    of a generalized Cartan matrix, see :class:`CartanMatrix`.

    EXAMPLES::

        sage: from sage.combinat.root_system.cartan_matrix import is_generalized_cartan_matrix
        sage: M = matrix([[2,-1,-2], [-1,2,-1], [-2,-1,2]])
        sage: is_generalized_cartan_matrix(M)
        True
        sage: M = matrix([[2,-1,-2], [-1,2,-1], [0,-1,2]])
        sage: is_generalized_cartan_matrix(M)
        False
        sage: M = matrix([[1,-1,-2], [-1,2,-1], [-2,-1,2]])
        sage: is_generalized_cartan_matrix(M)
        False

    A non-symmetrizable example::

        sage: M = matrix([[2,-1,-2], [-1,2,-1], [-1,-1,2]])
        sage: is_generalized_cartan_matrix(M)
        True
    """
    if not is_borcherds_cartan_matrix(M):
        return False
    n = M.ncols()
    return all(M[i,i] == 2 for i in range(n))


def find_cartan_type_from_matrix(CM):
    r"""
    Find a Cartan type by direct comparison of Dynkin diagrams given from
    the generalized Cartan matrix ``CM`` and return ``None`` if not found.

    INPUT:

    - ``CM`` -- a generalized Cartan matrix

    EXAMPLES::

        sage: # needs sage.graphs
        sage: from sage.combinat.root_system.cartan_matrix import find_cartan_type_from_matrix
        sage: CM = CartanMatrix([[2,-1,-1], [-1,2,-1], [-1,-1,2]])
        sage: find_cartan_type_from_matrix(CM)
        ['A', 2, 1]
        sage: CM = CartanMatrix([[2,-1,0], [-1,2,-2], [0,-1,2]])
        sage: find_cartan_type_from_matrix(CM)
        ['C', 3] relabelled by {1: 0, 2: 1, 3: 2}
        sage: CM = CartanMatrix([[2,-1,-2], [-1,2,-1], [-2,-1,2]])
        sage: find_cartan_type_from_matrix(CM)

    TESTS:

    Check that :issue:`35987` is fixed::

        sage: from sage.combinat.root_system.cartan_matrix import find_cartan_type_from_matrix
        sage: cm = CartanMatrix(['A',7]).subtype([2,3,5])
        sage: find_cartan_type_from_matrix(cm)
        A2xA1 relabelled by {1: 2, 2: 3, 3: 5}

        sage: cm = CartanMatrix(['B',10,1]).subtype([0,1,2,3,5,6,8,9,10])
        sage: ct = find_cartan_type_from_matrix(cm); ct
        D4xB3xA2 relabelled by {1: 0, 2: 2, 3: 1, 4: 3, 5: 8, 6: 9, 7: 10, 8: 5, 9: 6}
        sage: ct.dynkin_diagram()
            O 3
            |
            |
        O---O---O
        0   2   1
        O---O=>=O
        8   9   10
        O---O
        5   6
        D4xB3xA2 relabelled by {1: 0, 2: 2, 3: 1, 4: 3, 5: 8, 6: 9, 7: 10, 8: 5, 9: 6}
    """
    types = []
    relabel = []
    for S in CM.dynkin_diagram().connected_components_subgraphs():
        S = DiGraph(S) # We need a simple digraph here
        n = S.num_verts()
        # Build the list to test based upon rank
        if n == 1:
            relabel.append({1: S.vertices()[0]})
            types.append(CartanType(['A', 1]))
            continue

        test = [['A', n]]
        if n >= 2:
            if n == 2:
                test += [['G',2], ['A',2,2]]
            test += [['B',n], ['A',n-1,1]]
        if n >= 3:
            if n == 3:
                test.append(['G',2,1])
            test += [['C',n], ['BC',n-1,2], ['C',n-1,1]]
        if n >= 4:
            if n == 4:
                test.append(['F',4])
            test += [['D',n], ['B',n-1,1]]
        if n >= 5:
            if n == 5:
                test.append(['F',4,1])
            test.append(['D',n-1,1])
        if n == 6:
            test.append(['E',6])
        elif n == 7:
            test += [['E',7], ['E',6,1]]
        elif n == 8:
            test += [['E',8], ['E',7,1]]
        elif n == 9:
            test.append(['E',8,1])

        # Test every possible Cartan type and its dual
        found = False
        for x in test:
            ct = CartanType(x)
            T = DiGraph(ct.dynkin_diagram()) # We need a simple digraph here
            iso, match = T.is_isomorphic(S, certificate=True, edge_labels=True)
            if iso:
                types.append(ct)
                relabel.append(match)
                found = True
                break

            if ct == ct.dual():
                continue # self-dual, so nothing more to test

            ct = ct.dual()
            T = DiGraph(ct.dynkin_diagram()) # We need a simple digraph here
            iso, match = T.is_isomorphic(S, certificate=True, edge_labels=True)
            if iso:
                types.append(ct)
                relabel.append(match)
                found = True
                break
        if not found:
            return None

    if len(types) == 1:
        # Irreducible, so just relabel
        return CartanType(types[0]).relabel(relabel[0])
    ct = CartanType(types)
    # ct._index_relabelling is a dict ``(ind, j): i``, where i is an index of
    #   ``ct``, ``ind`` is the position in the list of types, and j is the
    #   corresponding index of the type number ``ind``.
    # In other words, the j-th node of ``types[ind]`` is the i-th node of ``ct``.
    mapping = {i: relabel[d[0]][d[1]] for d, i in ct._index_relabelling.items()}
    return ct.relabel(mapping)
