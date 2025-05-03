import random
from sage.graphs.planar_maps.MapPermutation import MapPermutation
from sage.graphs.planar_maps.RootedMap import RootedMap
import numpy as np
from sage.graphs.planar_maps.CustomSwap import CustomSwap
from queue import deque
from sage.graphs.planar_maps.PrimitiveMutableLabelledMap import PrimitiveMutableLabelledMap


class MapGenerator:
    """
    This class represents an abstraction containing
    methods to generate a Map.
    """

    def __init__(self):
        # Set it to true when in production
        # during debugging to False
        self._production = True

    def cube(self):
        """Returns the standard cube map."""
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

    def complete_map(self, n):
        """
        Returns an arbitrary rooted map corresponding to the complete
        graph with n nodes. The genus is guaranteed to be zero if the
        graph is planar (i.e., n <= 4).
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

    def getRandomDyckPath(self, n, seed=None):
        """
        Returns a random Dyck path of size n (uniform random generation).

        INPUT:
            - ``n`` -- int; size of the path
            - ``seed`` -- int; A random seed; if None is used, no random seed will be set.

        OUTPUT:
            A list of size 2*n with +1 for up and
            -1 for down steps in the Dyck path.

        EXAMPLE::
            sage: dyckPath = MapGenerator().getRandomDyckPath(10, seed=42)
            sage: dyckPath
            [1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, -1, -1]

        TESTS::
            sage: dyckPath = MapGenerator().getRandomDyckPath(50)
            sage: level = 0
            sage: for step in dyckPath:
            ....:     level += step
            ....:     assert level >= 0
            ....:
            sage: assert level == 0
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

    def getRandomPermutation(self, n, seed=None):
        """
        Returns a random permutation of size n.

        Args:
            n : The size of the permutation.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A random permutation of size n.
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))
        lst = [i + 1 for i in range(n)]
        rng.shuffle(lst)
        return MapPermutation(lst)

    def isValidDyckPath(self, dyckPathCandidate):
        """
        Checks whether the given Dyck path candidate is valid.

        Args:
            dyckPathCandidate : A list representing a potential Dyck path.

        Returns:
            A boolean indicating whether or not dyckPathCandidate is a
            correct Dyck path.
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

    def getTreeFromDyckPath(self, dyckPath, trust=False):
        """
        Given a Dyck path, this function returns the associated rooted tree.

        Args:
            dyckPath : A list representing a Dyck path, with +1 for up and -1 for down.
            trust: A boolean indicating whether to trust that we have a dyckPath

        Returns:
            The corresponding rooted plane tree if dyckPath is valid;
            otherwise, raises an error.

        Complexity:
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

    def getRandomLabellingTree(self, tree, seed=None):
        """
        Generates a uniformly random labelling of a tree.

        Args:
            tree : The input rooted tree.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A list of size 2*tree.m + 1 where labelling[i] (for i >= 1)
            represents the label of demi-edge i. The first value
            (labelling[0]) is set to -1 but has no meaning.

        Complexity:
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

    def getRandomTree(self, numberOfEdge, seed=None):
        """
        Generates a uniformly random rooted tree.

        Args:
            numberOfEdge : The number of edges in the tree.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A randomly selected rooted tree with numberOfEdge edges.

        Complexity:
            O(numberOfEdge)
        """
        return self.getTreeFromDyckPath(
            self.getRandomDyckPath(numberOfEdge, seed=seed), trust=self._production
        )

    def getRandomLabelledTree(self, numberOfEdge, seed=None):
        """
        Generates a uniformly random rooted tree along with a labelling.

        Args:
            numberOfEdge : The number of edges in the tree.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A tuple (tree, labelling) where:
            - tree : A randomly selected rooted tree with numberOfEdge edges.
            - labelling : A list of labels for the tree’s demi-edges.

        Complexity:
            O(numberOfEdge)
        """
        tree = self.getRandomTree(numberOfEdge, seed=seed)
        return tree, self.getRandomLabellingTree(tree, seed=seed)

    def getRandomPlanarQuadrangulation(self, numberOfFace, seed=None):
        """
        Generates a uniformly random rooted planar quadrangulation with a
        specified number of faces.

        Args:
            numberOfFace : The number of faces in the quadrangulation.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A randomly selected rooted planar quadrangulation with
            numberOfFace faces.

        Complexity:
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

    def getRandomPlanarMap(self, numberOfEdge, seed=None):
        """
        Generates a uniformly random rooted planar map with a specified
        number of edges.

        Args:
            numberOfEdge : The number of edges in the rooted map.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A randomly selected rooted planar map with numberOfEdge edges.

        Complexity:
            O(numberOfEdge)
        """
        quad = self.getRandomPlanarQuadrangulation(numberOfEdge, seed=seed)

        return quad.inverseQuadrangulation()

    def generateRandomBitstring(self, n, seed=None):
        """
        Args:
            n>=1
        Returns:
            A uniformly generated bit string of size 4*n-2 such that for every prefix
            3*n_1-n_0>-2 where n_1 is the number of 1 in the prefix
            and n_0 the number of 0 in the prefix 
        ------
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

    def checkPrefixCondition(self, bits):
        """
        Args:
            bits a list containing 0 and 1
        Returns:
            A boolean indicating if for every prefix 3*n_1-n_0>-2 where n_1 is the number of 1 in the prefix
            and n_0 the number of 0 in the prefix
        """
        current_sum = 0
        for j in range(len(bits)-1):
            bit = bits[j]
            current_sum += 3 if bit == 1 else -1
            if current_sum <= -2:
                return False
        return True

    def cyclicShift(self, bits, shift):
        """
        Args:
            -bits a list
            -shift a positive integer < len(shift)
        Returns:
            bits shifted by shift
        """
        return bits[shift:] + bits[:shift]

    def fastGenerateValidCodeword(self, n, seed=None):
        """
        Args:
            n>=1
        Returns:
            A random valid codeword i.e a sequence of size 4n-2 of 0 and 1 such
            that for every prefix 3*n_1-n_0>-2 where n_1 is the number of 1 in the prefix
            and n_0 the number of 0 in the prefix
        """
        L = 4 * n - 2
        b = self.generateRandomBitstring(n, seed=seed)
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

        def add(j, v):
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

        return self.cyclicShift(b, random.choice(valid_shifts))

    def getRandomRootedTwoLeafTree(self, n, seed=None):
        """
        Args:
            n>=1
        Returns:
            A randomly generated rooted on one leaf two leaf tree
        """
        b = self.fastGenerateValidCodeword(n, seed=seed)

        return self.rootedTwoLeafTreeFromBit(b)

    def rootedTwoLeafTreeFromBit(self, b):
        """
        Args:
            b a correct code word such of size of the form 4n-2, such that the sum of every prefix > -2
        Returns:
            The two leaf tree (a tree where each internal node has 2 leaf) associated to b rooted at a leaf
        -----
        O(n) where n is such that len(b)=4n-2
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

    def randomTreeToTriangulation(self, tree, seed=None):
        """
        Args:
            tree a two leaf tree rooted tree rooted at a leaf
        Returns:
            A triangulation between the two associated to the tree with equal probability
        -----
        O(n) where n is the size of the tree
        """
        def isOnInnerEdge(Z):
            return Z.n != Z and (Z.c).n != Z.c

        def isOnLeafEdge(Z):
            return not isOnInnerEdge(Z)

        def extremeOnEdgeAfter(Z):
            return ((Z.c).n).c

        def isSpecial(Z):
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

        def closeOp(U, V):
            W = U.addEdgeAfter()
            W.link(V)[0].contract()
            W.contract()
            return V

        def enclose(A, D):
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

    def getRandomTriangulation(self, n, seed=None):
        """
        Args:
            n>=1
        Returns: 
            A random rooted triangulation of size n (i.e with 2n faces, 3n edge and  n+2 node)
            uniformly
        ----
        O(n)
        """
        return self.randomTreeToTriangulation(self.getRandomRootedTwoLeafTree(n, seed=seed), seed=seed)
