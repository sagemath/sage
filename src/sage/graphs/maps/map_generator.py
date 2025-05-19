from collections import deque
import random
import numpy as np
from sage.graphs.maps.custom_swap import CustomSwap
from sage.graphs.maps.map_permutation import MapPermutation
from sage.graphs.maps.rooted_map import RootedMap
from sage.graphs.maps.primitive_mutable_labelled_map import PrimitiveMutableLabelledMap
from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge


class MapGenerator:
    """
    This class represents an abstraction containing
    methods to generate Map.
    """

    def __init__(self):
        """
        Init a MapGenerator instance

        EXAMPLES::

            sage: mg = MapGenerator()

        NOTE:

            Complexity is O(1)
        """
        # Set it to true when in production
        # during debugging to False
        self._production = True

    def cube(self) -> RootedMap:
        """
        OUTPUT:

        Returns the standard cube map.

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: mg.cube().pretty_print()
                        Alpha: [(1, 3), (2, 5), (4, 7), (6, 10), (8, 13), (9, 14), (11, 17), (12, 18), (15, 16), (19, 22), (20, 23), (21, 24)]
                        Sigma (Node): [(1, 2, 4), (3, 6, 9), (5, 8, 12), (7, 11, 16), (10, 15, 20), (13, 14, 19), (17, 18, 21), (22, 23, 24)]
                        Phi (Face): [(1, 6, 15, 7), (2, 8, 14, 3), (4, 11, 18, 5), (9, 19, 23, 10), (12, 21, 22, 13), (16, 20, 24, 17)]

        NOTE:

            Complexity is O(1)
        """
        return RootedMap(
            adj=[
                (5, 4, 2),
                (1, 3, 6),
                (4, 7, 2),
                (8, 3, 1),
                (8, 1, 6),
                (5, 2, 7),
                (3, 8, 6),
                (7, 4, 5),
            ],
            trust=self._production,
        )

    def complete_map(self, n: int) -> RootedMap:
        """

        INPUT:

         - ``n`` --  int ; ``n`` >=1

        OUTPUT:

            Returns an arbitrary rooted map corresponding to the complete
            graph with n nodes. The genus is guaranteed to be zero if the
            graph is planar (i.e., n <= 4).

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: mg.complete_map(3)
            Labelled map | Sigma : [2, 1, 4, 3, 6, 5], Alpha : [6, 3, 2, 5, 4, 1]

        NOTE:

            Complexity is O(n^2)
        """
        adj = list(tuple((j + i) % n + 1 for j in range(1, n))
                   for i in range(n))
        m = RootedMap(
            adj=adj,
            trust=self._production,
        )
        if n <= 4:
            m = m.force_planar()
        return m

    def getRandomDyckPath(self, n: int, seed: int | None = None) -> list[int]:
        """
        Returns a random Dyck path of size n (uniform random generation).

        INPUT:

        - ``n`` -- int; size of the path
        - ``seed`` -- int | None; A random seed; if None is used, no random seed will be set.

        OUTPUT:

            A list of size 2*n with +1 for up and
            -1 for down steps in the Dyck path.

        EXAMPLES::

            sage: MapGenerator().getRandomDyckPath(10,seed=42)
            [1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, -1, -1]

        TESTS::

            sage: dyckPath = MapGenerator().getRandomDyckPath(50)
            sage: level = 0
            sage: for step in dyckPath:
            ....:     level += step
            ....:     assert level >= 0
            ....:
            sage: assert level == 0

        NOTE:

            Complexity is O(n)
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))
        N = 2 * n + 1
        dyck = [1] * n + [-1] * (n + 1)
        rng.shuffle(dyck)
        level = 0
        minlevel = 0
        posmin = 0
        for i in range(N):
            level += dyck[i]
            if level < minlevel:
                posmin = i + 1
                minlevel = level
        Dyckfinal = dyck[posmin:] + dyck[:posmin]
        return Dyckfinal[:-1]

    def getRandomPermutation(self, n: int, seed: int | None = None) -> MapPermutation:
        """
        INPUT:

        - ``n`` -- int ; The size of the permutation.
        - ``seed`` -- int|None ; A random seed; if None is used, no random seed will be set.

        OUTPUT:

            A random permutation of size n, as a MapPermutation

        EXAMPLES::

            sage: MapGenerator().getRandomPermutation(4,seed=42)
            [3, 2, 4, 1]

        NOTE:

            Complexity is O(n)
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))
        lst = [i + 1 for i in range(n)]
        rng.shuffle(lst)
        return MapPermutation(lst)

    def isValidDyckPath(self, dyckPathCandidate: list[int]) -> bool:
        """
        Checks whether the given Dyck path candidate is valid.

        INPUT:

        - ``dyckPathCandidate`` -- List[int] ; A list representing a potential Dyck path.

        OUTPUT:

            A boolean indicating whether or not dyckPathCandidate is a
            correct Dyck path.

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: dyckPath = mg.getRandomDyckPath(10,seed =42)
            sage: mg.isValidDyckPath(dyckPath)
            True

        NOTE:

            Complexity is O(n)
        """
        if len(dyckPathCandidate) == 0 or len(dyckPathCandidate) % 2 == 1:
            return False

        for e in dyckPathCandidate:
            if e != -1 and e != 1:
                return False

        S = 0
        for e in dyckPathCandidate:
            S += e
            if S < 0:
                return False
        return S == 0

    def getTreeFromDyckPath(self, dyckPath: list[int], trust=False) -> RootedMap:
        """
        Given a Dyck path, this function returns the associated rooted tree.

        INPUT:

        - ``dyckPath`` -- List[int] ; A list representing a Dyck path, with +1 for up and -1 for down.
        - ``trust`` -- bool ; A boolean indicating whether to trust that we have a dyckPath


        OUTPUT:

            The corresponding rooted plane tree if dyckPath is valid;
            otherwise, raises an error.

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: dyckPath = mg.getRandomDyckPath(10,seed =42)
            sage: tree = mg.getTreeFromDyckPath(dyckPath).pretty_print()
                        Alpha: [(1, 3), (2, 5), (4, 6), (7, 9), (8, 10), (11, 13), (12, 15), (14, 17), (16, 18), (19, 20)]
                        Sigma (Node): [(1, 2, 4), (3,), (5,), (6, 7, 8), (9, 11, 12, 14), (10,), (13, 16), (15,), (17, 19), (18,), (20,)]
                        Phi (Face): [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]


        NOTE:

            O(k), where k = len(dyckPath)
        """
        if not trust and not self.isValidDyckPath(dyckPath):
            raise ValueError("The given list isn't a Dyck path")

        phiCycle = []
        alphaCycle = []
        p = []

        for i in range(len(dyckPath)):
            phiCycle.append(i + 1)
            if dyckPath[i] < 0:
                otherDemiEdge = p.pop()
                alphaCycle.append((i + 1, otherDemiEdge))
            else:
                p.append(i + 1)

        phi = MapPermutation([tuple(phiCycle)])
        alpha = MapPermutation(alphaCycle)
        sigma = phi.left_action_product(alpha)

        return RootedMap(
            alpha=alpha,
            sigma=sigma,
            trust=self._production,
        )

    def getRandomLabellingTree(self, tree: RootedMap, seed: int | None = None) -> list[int]:
        """
        Generates a uniformly random correct labelling of a tree.
        A function on the nodes of the tree considered up to translation
        such that if u and v are adjacent f(u) and f(v) differs by at most one

        INPUT:

        - ``tree`` - RootedMap ; The input rooted tree.
        - ``seed`` - int | None ;  A random seed; if None is used, no random seed will be set.

        OUTPUT:

            A list of size 2*tree.m + 1 where labelling[i] (for i >= 1)
            represents the label of demi-edge i. The first value
            (labelling[0]) is set to -1 but has no meaning.

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: tree = mg.getRandomTree(10,seed=42)
            sage: mg.getRandomLabellingTree(tree,seed = 42)
            [-1, 0, 0, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -2, -1, -2, -2, -2, -2, -2, -3]


        NOTE:

            O(m), where m is the number of edges in the tree.
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))

        sigma = tree.sigma
        alpha = tree.alpha

        p = deque()
        nodes = tree.nodes()

        demiEdgeToNodeId = [-1 for _ in range(2 * tree.m + 1)]

        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                demiEdgeToNodeId[nodes[i][j]] = i

        labelling = [-1 for _ in range(2 * tree.m + 1)]
        startNodeId = demiEdgeToNodeId[1]

        p.append(startNodeId)
        labellingNodes = [-1 for _ in range(len(nodes))]
        seen = [False for _ in range(len(nodes))]

        seen[startNodeId] = True
        transition = [-1, 1, 0]
        labellingNodes[startNodeId] = 0

        while len(p) > 0:
            nodeId = p.popleft()

            for demiEdge in nodes[nodeId]:
                labelling[demiEdge] = labellingNodes[nodeId]

                alphaDemiEdge = alpha(demiEdge)
                alphaNodeId = demiEdgeToNodeId[alphaDemiEdge]
                if not seen[alphaNodeId]:
                    dLabel = rng.sample(transition, 1)[0]
                    labellingNodes[alphaNodeId] = dLabel + \
                        labellingNodes[nodeId]
                    seen[alphaNodeId] = True
                    p.append(alphaNodeId)

        return labelling

    def getRandomTree(self, numberOfEdge: int, seed: int | None = None) -> RootedMap:
        """
        Generates a uniformly random rooted tree.

        INPUT:

        - ``numberOfEdge`` -- int ; The number of edges in the tree.
        - ``seed`` -- int | None ; A random seed; if None is used, no random seed will be set.

        OUTPUT:

            A randomly selected rooted tree with numberOfEdge edges.

        EXAMPLES::

            sage: MapGenerator().getRandomTree(10,seed=42)
            Rooted map | Sigma : [2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20] Alpha : [3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19]

        NOTE:

            O(numberOfEdge)
        """
        return self.getTreeFromDyckPath(
            self.getRandomDyckPath(numberOfEdge, seed=seed), trust=self._production
        )

    def getRandomLabelledTree(self, numberOfEdge: int, seed: int | None = None) -> tuple[RootedMap, list[int]]:
        """
        Generates a uniformly random rooted tree along with a labelling.

        INPUT:

        - ``numberOfEdge`` -- int ; The number of edges in the tree.
        - ``seed`` -- int ; A random seed if None is used, no random seed will be set.

        OUTPUT:

            A tuple (tree, labelling) where:
            - tree : A randomly selected rooted tree with numberOfEdge edges.
            - labelling : A list of labels for the tree’s demi-edges.

        EXAMPLES::

            sage: MapGenerator().getRandomLabelledTree(10,seed=42)
            (Rooted map | Sigma : [2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20] Alpha : [3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19],
             [-1,
              0,
              0,
              0,
              0,
              -1,
              -1,
              -1,
              -1,
              -1,
              0,
              -1,
              -1,
              -2,
              -1,
              -2,
              -2,
              -2,
              -2,
              -2,
              -3])

        NOTE:

            O(numberOfEdge)
        """
        tree = self.getRandomTree(numberOfEdge, seed=seed)
        return tree, self.getRandomLabellingTree(tree, seed=seed)

    def getRandomPlanarQuadrangulation(self, numberOfFace: int, seed: int | None = None) -> RootedMap:
        """
        Generates a uniformly random rooted planar quadrangulation with a
        specified number of faces.

        INPUT:

        - ``numberOfFace`` -- int ; The number of faces in the quadrangulation.
        - ``seed`` -- int | None ; A random seed; if None is used, no random seed will be set.

        OUTPUT:

            A randomly selected rooted planar quadrangulation with
            numberOfFace faces.

        EXAMPLES::

            sage: MapGenerator().getRandomPlanarQuadrangulation(10,seed=42).faces()
            [(1, 6, 9, 7),
             (2, 8, 13, 3),
             (4, 10, 16, 5),
             (11, 19, 24, 12),
             (14, 18, 17, 15),
             (20, 27, 37, 21),
             (22, 29, 34, 23),
             (25, 28, 38, 26),
             (30, 36, 35, 31),
             (32, 40, 39, 33)]

        NOTE:

            O(numberOfFace)
        """
        tree, labelling = self.getRandomLabelledTree(numberOfFace, seed=seed)
        quadA, quadB = tree.inverseShaefferTree(
            returnMarkedDemiEdge=False, labelled=labelling
        )

        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))

        if rng.random() < 0.5:
            return quadB
        return quadA

    def getRandomPlanarMap(self, numberOfEdge: int, seed: int | None = None) -> RootedMap:
        """
        Generates a uniformly random rooted planar map with a specified
        number of edges.

        INPUT:

        - ``numberOfEdge`` -- int; The number of edges in the rooted map.
        - ``seed`` -- int | None; A random seed; if None is used, no random seed will be set.

        OUTPUT:

            A randomly selected rooted planar map with numberOfEdge edges.


        EXAMPLES::

            sage: MapGenerator().getRandomPlanarMap(10,seed=42)
            Rooted map | Sigma : [2, 3, 4, 1, 6, 7, 8, 10, 11, 5, 12, 14, 16, 9, 13, 18, 19, 17, 20, 15] Alpha : [2, 1, 5, 6, 3, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19]

        NOTE:

            O(numberOfEdge)
        """
        quad = self.getRandomPlanarQuadrangulation(numberOfEdge, seed=seed)

        return quad.inverseQuadrangulation()

    def generateRandomBaseTwoLeafBitString(self, n: int, seed: int | None = None) -> list[int]:
        """
        INPUT:

        - ``n`` -- int ; ``n``>=1

        OUTPUT:

            A uniformly generated base two leaf bit string i.e a bit string of size 4n-2, with n-1  1

        EXAMPLES:

            sage: len(MapGenerator().generateRandomBaseTwoLeafBitString(10,seed=42))
            38
            sage: MapGenerator().generateRandomBaseTwoLeafBitString(5, seed=42)
            [1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]

        NOTE:

            O(n)
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))

        L = 4 * n - 2
        weight = n - 1
        bits = [0] * L
        one_positions = rng.sample(range(L), weight)
        for pos in one_positions:
            bits[pos] = 1
        return bits

    def checkPrefixCondition(self, bits: list[int]) -> bool:
        """
        INPUT:

        - ``bits`` -- List[int] ; a list containing 0 and 1

        OUTPUT:

            A boolean indicating if for every prefix 3*n_1-n_0>-2 where n_1 is the number of 1 in the prefix
            and n_0 the number of 0 in the prefix

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: mg.getRandomTwoLeafBitString(4,seed=42)
            [0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]
            sage: mg.checkPrefixCondition(mg.getRandomTwoLeafBitString(4,seed=42))
            True

        NOTE:

            O(len(bits))
        """
        current_sum = 0
        for j in range(len(bits)-1):
            bit = bits[j]
            current_sum += 3 if bit == 1 else -1
            if current_sum <= -2:
                return False
        return True

    def cyclicShift(self, bits: list[int], shift: int) -> list[int]:
        """
        INPUT:

        -``bits`` -- List[int]
        - ``shift`` -- int ; a positive integer < len(shift)

        OUTPUT:

            bits shifted by shift

        EXAMPLES::
            sage: mg = MapGenerator()
            sage: lst = mg.getRandomDyckPath(4,seed=42)
            sage: lst
            [1, -1, 1, 1, 1, -1, -1, -1]
            sage: mg.cyclicShift(lst,2)
            [1, 1, 1, -1, -1, -1, 1, -1]

        NOTE:
            O(len(bits))
        """
        return bits[shift:] + bits[:shift]

    def getRandomTwoLeafBitString(self, n: int, seed: int | None = None) -> list[int]:
        """
        INPUT:

        - ``n`` -- int; ``n``>=1

        OUTPUT:

            A random two leaf bit string i.e  a sequence of size 4n-2 of 0 and 1 such
            that for every prefix 3*n_1-n_0>-2 where n_1 is the number of 1 in the prefix
            and n_0 the number of 0 in the prefix

        EXAMPLES::

            sage: MapGenerator().getRandomTwoLeafBitString(4,seed=42)
            [0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]

        NOTE:

            O(n)
        """
        L = 4 * n - 2
        b = self.generateRandomBaseTwoLeafBitString(n, seed=seed)
        elems = [4 * bit - 1 for bit in b]          # maps 0 to -1 and 1 to 3

        q = deque()
        s = 0
        # our goal is to find the valid shifts of b, i. e. such that the value of every prefix is > -2 (and = -2 for the whole string)
        # note that here, the value of a string is 3 * (number of bits equal ro 1) - (number of bits equal to 0)

        # at each iteration i, s is the value of b from index 0 to i excluded
        # q is a deque of pairs (index, value) such that, at each iteration i:
        #   - if (j, v) is in q, the value of the bitstring from i to j included is v - s
        #   - the deque is ordered both by increasing indices and values (ie if (j, v) is before (j', v') in q, it holds that j<j' and v<=v')

        # note that we aren't interested in (j, v), (j', v') such that j < j' and v > v' because
        # the sum from i to j' will always be greater than the sum from i to j

        # at each iteration i, we just have to check that the value of the minimum, ie. the first element of the deque - s, is equal to -2
        # and that this minimum is reached at the end of the deque

        def add(j: int, v: int) -> None:
            while q and q[-1][1] > v:
                q.pop()
            q.append((j, v))

        # initialize the deque
        for i in range(L):
            s += elems[i]
            add(i, s)

        s = 0

        valid_shifts = []

        for i in range(L):
            if q[0][1] - s == -2 and q[0][0] == i + L - 1:
                valid_shifts.append(i)

            if q[0][0] == i:
                q.popleft()
            s += elems[i]

            add(i + L, -2 + s)

        assert len(valid_shifts) == 2

        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))

        return self.cyclicShift(b, rng.choice(valid_shifts))

    def getRandomRootedTwoLeafTree(self, n: int, seed: int | None = None) -> RootedMap:
        """
        INPUT:

        - ``n`` -- int; ``n``>=1

        OUTPUT:

            A randomly generated rooted on one leaf two leaf tree

        EXAMPLES::

            sage: MapGenerator().getRandomRootedTwoLeafTree(4,seed=42)
            Rooted map | Sigma : [1, 3, 4, 6, 5, 2, 9, 10, 11, 13, 7, 12, 16, 18, 15, 8, 17, 20, 19, 14, 21, 22] Alpha : [2, 1, 5, 7, 3, 8, 4, 6, 12, 14, 15, 9, 17, 10, 11, 19, 13, 21, 16, 22, 18, 20]

        NOTE:

            O(n)
        """
        b = self.getRandomTwoLeafBitString(n, seed=seed)

        return self.rootedTwoLeafTreeFromBit(b)

    def rootedTwoLeafTreeFromBit(self, b: list[int]) -> RootedMap:
        """
        INPUT:

            - ``b`` -- LabelledMap; a two leaf bit string

        OUTPUT:

            The two leaf tree (a tree where each internal node has 2 leaf) associated to b rooted at a leaf

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: bits = mg.getRandomTwoLeafBitString(4,seed = 42)
            sage: mg.rootedTwoLeafTreeFromBit(bits)
            Rooted map | Sigma : [1, 3, 4, 6, 5, 2, 9, 10, 11, 13, 7, 12, 16, 18, 15, 8, 17, 20, 19, 14, 21, 22] Alpha : [2, 1, 5, 7, 3, 8, 4, 6, 12, 14, 15, 9, 17, 10, 11, 19, 13, 21, 16, 22, 18, 20]

        NOTE:

            O(len(b))
        """
        n = (len(b)+2)//4

        alphaCycle = [(i, i+1) for i in range(1, 6*n-1, 2)]
        sigmaCycle = []

        p = [1]  # pile contenant les points fixes non couplés
        h = 0  # hauteur dans l'arbre

        H = [0]  # pile des hauteurs associés aux elements de p
        D = [0]*(6*n-1)  # hauteur de toutes les demi-edges
        Q = [[] for _ in range(n)]
        # Is the node full (in term of leaf)
        isFull = [False for _ in range(n)]
        Q[0].append(2)
        D[1] = 0
        D[2] = 0
        Q.append([1])
        j = 4
        i = 0
        while i < len(b) and j <= 6 * n - 2:
            if b[i] == 1:
                D[j] = h+1
                D[j-1] = h
                Q[h].append(j-1)
                Q[h+1].append(j)
                j += 2
                h += 1
                i += 1
                continue
            # b[i] == 0
            if H and H[-1] == h:
                D[j] = h
                D[j-1] = h
                Q[h].append(j-1)
                Q.append([j])
                p.pop()
                H.pop()

                isFull[h] = True

                j += 2
                i += 1
                continue

            if isFull[h]:
                if Q[h]:
                    sigmaCycle.append(tuple(Q[h]))
                    Q[h] = []
                isFull[h] = False
                i += 1
                h -= 1
                continue

            D[j] = h
            D[j-1] = h

            Q[h].append(j-1)
            p.append(j)
            Q.append([j])
            H.append(h)

            j += 2
            i += 1

        for level in range(len(Q)):
            if Q[level]:
                sigmaCycle.append(tuple(Q[level]))

        sigma = MapPermutation(sigmaCycle)
        alpha = MapPermutation(alphaCycle)

        tree = RootedMap(sigma=sigma, alpha=alpha)

        return tree

    def randomTreeToTriangulation(self, tree: RootedMap, seed: int | None = None) -> RootedMap:
        """
        INPUT:

        - ``tree`` -- LabelledMap ; A two leaf tree rooted tree rooted at a leaf

        OUTPUT:

            A triangulation between the two associated to the tree with equal probability

        EXAMPLES::

            sage: mg = MapGenerator()
            sage: tree = mg.getRandomRootedTwoLeafTree(4,seed=42)
            sage: tri = mg.randomTreeToTriangulation(tree,seed=42)
            sage: tri.isTriangulation()
            True

        NOTE:

            O(n) where n is the size of the tree
        """
        def isOnInnerEdge(Z: TopologicalDemiEdge) -> bool:
            return Z.n != Z and (Z.c).n != Z.c

        def isOnLeafEdge(Z: TopologicalDemiEdge) -> bool:
            return not isOnInnerEdge(Z)

        def extremeOnEdgeAfter(Z: TopologicalDemiEdge) -> bool:
            return ((Z.c).n).c

        def isSpecial(Z: TopologicalDemiEdge) -> bool:
            return isOnLeafEdge(extremeOnEdgeAfter(Z)) and isOnLeafEdge(Z)

        triangulation = PrimitiveMutableLabelledMap(lmap=tree)

        root = triangulation.X(1)

        outerList = []
        isLeaf = np.zeros(2*triangulation.m+1)

        for A in root.face():
            if A.n == A:
                outerList.append(A)
                isLeaf[A.raw] = True

        cntSpecial = 0
        for A in outerList:
            cntSpecial += isSpecial(A)

        A = triangulation.X(1)
        B, C = A.f, (A.f).f

        while cntSpecial > 2:
            if isOnInnerEdge(A) and isOnInnerEdge(B) and isOnLeafEdge(C):

                cntSpecial -= isSpecial(C.c)

                N, _ = A.link(C.c)
                N.contract()
                A = (C.c).pf
                B, C = A.f, (A.f).f
            else:
                A = A.f
                B = B.f
                C = C.f

        specialList = []
        for A in outerList:
            if isSpecial(A):
                specialList.append((A, extremeOnEdgeAfter(A)))

        A, B = specialList[0]
        C, D = specialList[1]

        def closeOp(U: TopologicalDemiEdge, V: TopologicalDemiEdge) -> TopologicalDemiEdge:
            W = U.addEdgeAfter()
            W.link(V)[0].contract()
            W.contract()
            return V

        def enclose(A: TopologicalDemiEdge, D: TopologicalDemiEdge) -> None:
            V = A.c.pf
            U = A
            while True:
                nxt = V.pf
                if V.n == V:
                    U = closeOp(U, V)
                if V == D:
                    break
                V = nxt
        if A != D:
            enclose(A, D)
            enclose(C, B)

        X, Y = A.link(C)

        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))

        tau = CustomSwap([(X.raw, 1)])

        if rng.random() < 0.5:
            tau = CustomSwap([(Y.raw, 1)])

        triangulation = triangulation.relabel(tau)
        return RootedMap(triangulation)

    def getRandomTriangulation(self, n: int, seed: int | None = None) -> RootedMap:
        """
        INPUT:

        - ``n`` --  int ; ``n`` >=1

        OUTPUT:

            A random rooted triangulation of size n (i.e with 2n faces, 3n edge and  n+2 node)
            uniformly

        EXAMPLES::

            sage: MapGenerator().getRandomTriangulation(22,seed=42).isTriangulation()
            True

        NOTE:

            O(n)
        """
        return self.randomTreeToTriangulation(self.getRandomRootedTwoLeafTree(n, seed=seed), seed=seed)
