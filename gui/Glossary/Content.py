# There is likely a better solution to this, but for now this is how content will be stored
def outputContent(topic):
    if topic == "GT":
        return {"Graph" : "A <b>graph</b> G is a pair of sets (V,E) where the members of V are called <b>vertices</b>"
                +" and the members of E are unordered pairs of distinct vertices called <b>edges</b>." 
                "We say two vertices are <b>adjacent</b> if they are the elements of an edge. We then say the edge is <b>incident</b> to those vertices.",
                "Degree":"For a vertex x, the set of vertices adjacent to x is denoted N(x). We say the <b>degree</b> of x is the size of N(x), denoted |N(x)|.",
                "Subgraph":"A graph S is a <b>subgraph</b> of a graph G if V(S) ⊆ V(G) and E(S) ⊆ E(G).",
                "Complete Graph":"A <b>complete graph</b> on n vertices, denoted K(n), is a graph on n vertices such that every vertex is adjacent to every other vertex.",
                "Walk":"A <b>walk</b> is an alternating sequence of vertices and edges.",
                "Path":"A <b>path</b> is a walk where all vertices are distinct.",
                "Trail":"If the edges in a walk are distinct, the walk is called a <b>trail</b>.",
                "Circuit":"A closed trail where end vertices coincide, is called a <b>circuit</b>.",
                "Hamiltonian":"A graph is <b>Hamiltonian</b> if there is a cycle that traverses ever vertex exactly once.",
                "Eularian":"A graph is <b>Eularian</b> if there is a circuit that contains all the edges of G.",
                "Planar":"A graph is <b>planar</b> if it can be drawn on the plane with no edge crossings.",
                "Subdivision":"A graph G is a <b>subdivision</b> of a graph H if edges in H are replaced with paths where two such paths can have only their end vertices in common.",
                "Connected":"A graph is <b>connected</b> if for every pair of distinct vertices, there is a path joining those vertices.",
                "Handshaking Lemma":"<i>Let G be a simple graph where</i> a,b,c,...,x,y,z <i>are vertices of</i> G, <i>then</i> |N(a)|+|N(b)|+|N(c)|+...+|N(x)|+|N(y)|+|N(z)|=2|E(G)|",
                "Euler's Formula":"<i>If G is a connected planar graph with</i> n <i>vertices,</i> m <i>edges, and</i> f <i>faces, then </i>n-m+f=2.",
                "Bipartite Graph":"A graph is <b>bipartite</b> if the vertex set can be partitioned into two disjoint sets A and B and for every edge e={x,y}, we have that x is in A and y is in B.",
                "Kuratowski's Theorem":"<i>A graph is planar if and only if it has no subgraph that is a subdivision of K(5) or K(3,3)</i>",
                "Density":"The <b>density</b> of a graph G with n vertices and m edges is given by D=2m/n(n-1)."}
    elif topic == "LA":
        return {"Vector": "A <b>vector</b> is a quantity that has both direction and magnitude.",
                 "Matrix": "Let a(i,j) be a real number where 1≤i≤n and 1≤j≤m. A n×m <b>(real) matrix</b> A a function such that A=(a(i,j)) for i=1,...,n and j=1,...,m. We say a matrix is <b>square</b> if n=m.",
                 "Linearly independent":"A set of vectors {x,y,z} are <b>linearly independent</b> if for all real numbers a,b, and c we have that ax+by+cz=0" ,
                 "Basis":"The <b>basis</b> of a vector space X is the set of all linearly independent vectors that form all linear combinations of vectors in X",
                 "Dimenseion":"The dimension of a space X, is the size of basis.",
                 "Eigenvalue":"Let A be an n×m matrix and let v be an m-dimensional column vector. The eigenvalue of a matrix A is a real number k such that Av=kv.",
                 "Eigenvector":"Let A be an n×m matrix and let k be an eigenvalue of A. The <b>eigenvector</b> is the vector v such that Av=kv.",
                 "Identity matrix":"The <b>identity matrix</b>, denoted I, is a matrix such that a(i,j)=1 for i=j and a(i,j)=0 otherwise.",
                 "Inverse matrix":"The <b>inverse matrix</b> of a matrix A is denoted A⁻¹ and is a matrix such that AA⁻¹=A⁻¹A=I.",
                 "Diagonal matrix": "A <b>diagonal matrix</b> is a matrix D such that a(i,j)=0 for all i≠j.",
                 "Diagonalizable": "A square matrix A is <b>diagonalizable</b> fi there exists an invertible matrix P and a diagonal matrix D such that P⁻¹AP=D.",
                 "Matrix transpose":"The <b>transpose</b> of a matrix A is a matrix Aᵀ such that aᵀ(i,j)=a(j,i).",
                 "Dot product":"Let u=ax+by+cz and v=dx+ey+fz, the <b>dot product</b> of u and v, denoted u•v, is given by u•v=ad+be+cf.",
                 "Orthogonal matrix":"An <b>orthogonal matrix</b> Q is a matrix such that QᵀQ=QQᵀ=I." }
    elif topic == "ST":
        return {"Set": "A <b>set</b> is a collection of distinct objects, considered as an object in its own right. Sets are typically denoted using curly brackets, e.g., {1, 2, 3}.",
                "Subset":"A set A is a <b>subset</b> of a set B if every element of A is also an element of B. This is denoted as A ⊆ B.",
                "Union":"The <b>union</b> of two sets A and B, denoted A ∪ B, is the set of elements that are in A, in B, or in both.",
                "Intersection":"The <b>intersection</b> of two sets A and B, denoted A ∩ B, is the set of elements that are in both A and B.",
                "Difference":"The <b>difference</b> of two sets A and B, denoted A - B, is the set of elements that are in A but not in B.",
                "Compliment":"Let A⊆U. The <b>complement</b> of a set A, denoted A', is the set of all elements not in A, but in U.",
                "Power Set":"The <b>power set</b> of a set A, denoted P(A), is the set of all possible subsets of A, including the empty set and A itself.",
                "Cartesian Product":"The <b>Cartesian product</b> of two sets A and B, denoted A × B, is the set of all ordered pairs (a, b) where a ∈ A and b ∈ B.",
                "Partition":"A <b>partition</b> of a set A is a collection of non-empty, disjoint subsets of A whose union is A.",
                "Cardinality":'The <b>cardinality</b> of a set A, denoted |A|, is a measure of the "number of elements" in A.',
                "Set Equality":"Two sets A and B are <b>equal</b>, denoted A=B, if they contain exactly the same elements.",
                "Disjoint Sets":"Two sets A and B are <b>disjoint</b> if their intersection is the empty set, i.e., A ∩ B = ∅.",
                "Symmetric Difference":"The <b>symmetric difference</b> of two sets A and B, denoted A Δ B, is the set of elements that are in either of the sets A or B but not in their intersection.",
                "Mutual Containment": "Two sets A and B are equal if and only if A⊆B and B⊆A."}
    else:
        return {f"key_{i}": f"value_{i}" for i in range(1000)}