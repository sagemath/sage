"""
Cartan matrices

AUTHORS:

- Travis Scrimshaw (2012-04-22): Nicolas M. Thiery moved matrix creation to
  :class:`CartanType` to prepare :func:`cartan_matrix()` for deprecation.
- Christian Stump, Travis Scrimshaw (2013-04-13): Created :class:`CartanMatrix`.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2012,2013 Travis Scrimshaw <tscrim at ucdavis.edu>,
#       Copyright (C) 2013 Chrisitan Stump,
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
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.matrix.matrix_integer_sparse import Matrix_integer_sparse
from sage.rings.all import ZZ
from sage.combinat.root_system.cartan_type import CartanType, CartanType_abstract
from sage.combinat.root_system.root_system import RootSystem
from sage.sets.family import Family

class CartanMatrix(Matrix_integer_sparse, CartanType_abstract):
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

    INPUT:

    Can be anything which is accepted by ``CartanType`` or a matrix.

    If given a matrix, one can also use the keyword ``cartan_type`` when giving
    a matrix to explicitly state the type. Otherwise this will try to check the
    input matrix against possible standard types of Cartan matrices. To disable
    this check, use the keyword ``cartan_type_check = False``.

    EXAMPLES::

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

    .. NOTE::

        Since this is a matrix, :meth:`row()` and :meth:`column()` will return
        the standard row and column respectively. To get the row with the
        indices as in Dynkin diagrams/Cartan types, use
        :meth:`row_with_indices()` and :meth:`column_with_indices()`
        respectively.
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **kwds):
        """
        Normalize input so we can inherit from spare integer matrix.

        .. NOTE::

            To disable the Cartan type check, use the optional argument
            ``cartan_type_check = False``.

        EXAMPLES::

            sage: C = CartanMatrix(['A',1,1])
            sage: C2 = CartanMatrix([[2, -2], [-2, 2]])
            sage: C3 = CartanMatrix(matrix([[2, -2], [-2, 2]]), [0, 1])
            sage: C == C2 and C == C3
            True

        TESTS:

        Check that :trac:`15740` is fixed::

            sage: d = DynkinDiagram()
            sage: d.add_edge('a', 'b', 2)
            sage: d.index_set()
            ('a', 'b')
            sage: cm = CartanMatrix(d)
            sage: cm.index_set()
            ('a', 'b')
        """
        # Special case with 0 args and kwds has Cartan type
        if "cartan_type" in kwds and len(args) == 0:
            args = (CartanType(kwds["cartan_type"]),)
        if len(args) == 0:
            data = []
            n = 0
            index_set = tuple()
            cartan_type = None
            subdivisions = None
        elif len(args) == 4 and isinstance(args[0], MatrixSpace): # For pickling
            return typecall(cls, args[0], args[1], args[2], args[3])
        elif isinstance(args[0], CartanMatrix):
            return args[0]
        else:
            cartan_type = None
            dynkin_diagram = None
            subdivisions = None

            from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
            if isinstance(args[0], DynkinDiagram_class):
                dynkin_diagram = args[0]
                cartan_type = args[0]._cartan_type
            else:
                try:
                    cartan_type = CartanType(args[0])
                    dynkin_diagram = cartan_type.dynkin_diagram()
                except (TypeError, ValueError):
                    pass

            if dynkin_diagram is not None:
                n = dynkin_diagram.rank()
                index_set = dynkin_diagram.index_set()
                reverse = dict((index_set[i], i) for i in range(len(index_set)))
                data = {(i, i): 2 for i in range(n)}
                for (i,j,l) in dynkin_diagram.edge_iterator():
                    data[(reverse[j], reverse[i])] = -l
            else:
                M = matrix(args[0])
                if not is_generalized_cartan_matrix(M):
                    raise ValueError("the input matrix is not a generalized Cartan matrix")
                n = M.ncols()
                if "cartan_type" in kwds:
                    cartan_type = CartanType(kwds["cartan_type"])
                elif n == 1:
                    cartan_type = CartanType(['A', 1])
                elif kwds.get("cartan_type_check", True):
                    cartan_type = find_cartan_type_from_matrix(M)
                data = M.dict()
                subdivisions = M._subdivisions

            if len(args) == 1:
                if cartan_type is not None:
                    index_set = tuple(cartan_type.index_set())
                elif dynkin_diagram is None:
                    index_set = tuple(range(n))
            elif len(args) == 2:
                index_set = tuple(args[1])
                if len(index_set) != n and len(set(index_set)) != n:
                    raise ValueError("the given index set is not valid")
            else:
                raise ValueError("too many arguments")

        mat = typecall(cls, MatrixSpace(ZZ, n, sparse=True), data, cartan_type, index_set)
        mat._subdivisions = subdivisions
        return mat

    def __init__(self, parent, data, cartan_type, index_set):
        """
        Initialize ``self``.

        TESTS::

            sage: C = CartanMatrix(['A',1,1])
            sage: TestSuite(C).run(skip=["_test_category", "_test_change_ring"])
        """
        Matrix_integer_sparse.__init__(self, parent, data, False, False)
        self._cartan_type = cartan_type
        self._index_set = index_set
        self.set_immutable()

    def root_system(self):
        """
        Return the root system corresponding to ``self``.

        EXAMPLES::

            sage: C = CartanMatrix(['A',3])
            sage: C.root_system()
            Root system of type ['A', 3]
        """
        if self._cartan_type is not None:
            return RootSystem(self._cartan_type)
        return self.dynkin_diagram().root_system()

    def root_space(self):
        """
        Return the root space corresponding to ``self``.

        EXAMPLES::

            sage: C = CartanMatrix(['A',3])
            sage: C.root_space()
            Root space over the Rational Field of the Root system of type ['A', 3]
        """
        return self.root_system().root_space()

    def reflection_group(self, type="matrix"):
        """
        Return the reflection group corresponding to ``self``.

        EXAMPLES::

            sage: C = CartanMatrix(['A',3])
            sage: C.reflection_group()
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

        raise ValueError("The reflection group is only available as a matrix group or as a permutation group.")

    def symmetrizer(self):
        """
        Return the symmetrizer of ``self``.

        EXAMPLES::

            sage: cm = CartanMatrix([[2,-5],[-2,2]])
            sage: cm.symmetrizer()
            Finite family {0: 2, 1: 5}

        TESTS:

        Check that the symmetrizer computed from the Cartan matrix agrees
        with the values given by the Cartan type::

            sage: ct = CartanType(['B',4,1])
            sage: ct.symmetrizer()
            Finite family {0: 2, 1: 2, 2: 2, 3: 2, 4: 1}
            sage: ct.cartan_matrix().symmetrizer()
            Finite family {0: 2, 1: 2, 2: 2, 3: 2, 4: 1}
        """
        sym = self.is_symmetrizable(True)
        if not sym:
            raise ValueError("the Cartan matrix is not symmetrizable")
        iset = self.index_set()
        # The result from is_symmetrizable needs to be scaled
        # to integer coefficients
        from sage.rings.arith import LCM
        from sage.rings.all import QQ
        scalar = LCM(map(lambda x: QQ(x).denominator(), sym))
        return Family( {iset[i]: ZZ(val*scalar) for i, val in enumerate(sym)} )

    @cached_method
    def symmetrized_matrix(self):
        """
        Return the symmetrized matrix of ``self`` if symmetrizable.

        EXAMPLES::

            sage: cm = CartanMatrix(['B',4,1])
            sage: cm.symmetrized_matrix()
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

            sage: C = CartanMatrix(['A',4,1])
            sage: C.cartan_type()
            ['A', 4, 1]

        If the Cartan type is unknown::

            sage: C = CartanMatrix([[2,-1,-2], [-1,2,-1], [-2,-1,2]])
            sage: C.cartan_type()
            [ 2 -1 -2]
            [-1  2 -1]
            [-2 -1  2]
        """
        if self._cartan_type is None:
            return self
        return self._cartan_type

    def rank(self):
        r"""
        Return the rank of ``self``.

        EXAMPLES::

            sage: CartanMatrix(['C',3]).rank()
            3
            sage: CartanMatrix(["A2","B2","F4"]).rank()
            8
        """
        return self.ncols()

    def relabel(self, relabelling):
        """
        Return the relabelled Cartan matrix.

        EXAMPLES::

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

            sage: CartanMatrix(['C',3]).cartan_matrix()
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
        Implements :meth:`CartanType_abstract.is_simply_laced()`.

        A Cartan matrix is simply-laced if all non diagonal entries are `0`
        or `-1`.

        EXAMPLES::

            sage: cm = CartanMatrix([[2, -1, -1, -1], [-1, 2, -1, -1], [-1, -1, 2, -1], [-1, -1, -1, 2]])
            sage: cm.is_simply_laced()
            True
        """
        for i in range(self.nrows()):
            for j in range(i+1, self.ncols()):
                if self[i, j] < -1 or self[j, i] < -1:
                    return False
        return True

    def is_crystallographic(self):
        """
        Implements :meth:`CartanType_abstract.is_crystallographic`.

        A Cartan matrix is crystallographic if it is symmetrizable.

        EXAMPLES::

            sage: CartanMatrix(['F',4]).is_crystallographic()
            True
        """
        return self.is_symmetrizable()

    def column_with_indices(self, j):
        """
        Return the `j^{th}` column `(a_{i,j})_i` of ``self`` as a container
        (or iterator) of tuples `(i, a_{i,j})`

        EXAMPLES::

            sage: M = CartanMatrix(['B',4])
            sage: [ (i,a) for (i,a) in M.column_with_indices(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        return self.dynkin_diagram().column(j)

    def row_with_indices(self, i):
        """
        Return the `i^{th}` row `(a_{i,j})_j` of ``self`` as a container
        (or iterator) of tuples `(j, a_{i,j})`

        EXAMPLES::

            sage: M = CartanMatrix(['C',4])
            sage: [ (i,a) for (i,a) in M.row_with_indices(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        return self.dynkin_diagram().row(i)

    def is_finite(self):
        """
        Return if ``self`` is a finite type or ``False`` if unknown.

        EXAMPLES::

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
            return self.det() > 0
        return self._cartan_type.is_finite()

    def is_affine(self):
        """
        Return if ``self`` is an affine type or ``False`` if unknown.

        EXAMPLES::

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
            return self.det() == 0
        return self._cartan_type.is_affine()

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
    if not is_Matrix(M):
        return False
    if not M.is_square():
        return False
    n = M.ncols()
    for i in xrange(n):
        if M[i,i] != 2:
            return False
        for j in xrange(i+1, n):
            if M[i,j] > 0 or M[j,i] > 0:
                return False
            elif M[i,j] == 0 and M[j,i] != 0:
                return False
            elif M[j,i] == 0 and M[i,j] != 0:
                return False
    return True

def find_cartan_type_from_matrix(CM):
    """
    Find a Cartan type by direct comparison of matrices given from the
    generalized Cartan matrix ``CM`` and return ``None`` if not found.

    INPUT:

    - ``CM`` -- A generalized Cartan matrix

    EXAMPLES::

        sage: from sage.combinat.root_system.cartan_matrix import find_cartan_type_from_matrix
        sage: M = matrix([[2,-1,-1], [-1,2,-1], [-1,-1,2]])
        sage: find_cartan_type_from_matrix(M)
        ['A', 2, 1]
        sage: M = matrix([[2,-1,0], [-1,2,-2], [0,-1,2]])
        sage: find_cartan_type_from_matrix(M)
        ['C', 3]
        sage: M = matrix([[2,-1,-2], [-1,2,-1], [-2,-1,2]])
        sage: find_cartan_type_from_matrix(M)
    """
    n = CM.ncols()
    # Build the list to test based upon rank
    if n == 1:
        return CartanType(['A', 1])

    test = [['A', n]]
    if n >= 2:
        if n == 2:
            test += [['G',2], ['A',2,2]]
        test += [['B',n], ['A',n-1,1]]
    if n >= 3:
        if n == 3:
            test += [['G',2,1], ['D',4,3]]
        test += [['C',n], ['BC',n-1,2], ['C',n-1,1]]
    if n >= 4:
        if n == 4:
            test += [['F',4], ['G',2,1], ['D',4,3]]
        test += [['D',n], ['B',n-1,1]]
    if n == 5:
        test += [['F',4,1], ['D',n-1,1]]
    elif n == 6:
        test.append(['E',6])
    elif n == 7:
        test += [['E',7], ['E',6,1]]
    elif n == 8:
        test += [['E',8], ['E',7,1]]
    elif n == 9:
        test.append(['E',8,1])

    # Test every possible Cartan type and its dual
    for x in test:
        ct = CartanType(x)
        if ct.cartan_matrix() == CM:
            return ct
        if ct == ct.dual():
            continue
        ct = ct.dual()
        if ct.cartan_matrix() == CM:
            return ct
    return None

def cartan_matrix(t):
    """
    Return the Cartan matrix of type `t`.

    .. NOTE::

        This function is deprecated in favor of
        ``CartanMatrix(...)``, to avoid polluting the
        global namespace.

    EXAMPLES::

        sage: cartan_matrix(['A', 4])
        doctest:...: DeprecationWarning: cartan_matrix() is deprecated. Use CartanMatrix() instead
        See http://trac.sagemath.org/14137 for details.
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
    """
    from sage.misc.superseded import deprecation
    deprecation(14137, 'cartan_matrix() is deprecated. Use CartanMatrix() instead')
    return CartanMatrix(t)

