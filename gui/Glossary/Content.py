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
                "Curcuit":"A closed trail where end vertices coincide, is called a <b>circuit</b>.",
                "Hamiltonian":"A graph is <b>Hamiltonian</b> if there is a cycle that traverses ever vertex exactly once.",
                "Eularian":"A graph is <b>Eularian</b> if there is a circuit that contains all the edges of G.",
                "Planar":"A graph is <b>planar</b> if it can be drawn on the plane with no edge crossings.",
                "Subdivision":"A graph G is a <b>subdivision</b> of a graph H if edges in H are replaced with paths where two such paths can have only their end vertices in common.",
                "Connected":"A graph is <b>connected</b> if for every pair of distinct vertices, there is a path joining those vertices.",
                "Handshaking Lemma":"<i>Let G be a simple graph where</i> a,b,c,...,x,y,z <i>are vertices of</i> G, <i>then</i> |N(a)|+|N(b)|+|N(c)|+...+|N(x)|+|N(y)|+|N(z)|=2|E(G)|",
                "Euler's Formula":"<i>If G is a connected planar graph with</i> n <i>vertices,</i> m <i>edges, and</i> f <i>faces, then </i>n-m+f=2.",
                "Bipartite Graph":"A graph is <b>bipartite</b> if the vertex set can be partitioned into two disjoint sets A and B and for every edge e={x,y}, we have that x is in A and y is in B.",
                "Kuratowski's Theorem":"<i>A graph is planar if and only if it has no subgraph that is a subdivision of K(5) or K(3,3)</i>"}
    elif topic == "LA":
        return {"Linear Map": "",
                 "Vector space": "",
                 "Orthogonality":""}
    else:
        return {f"key_{i}": f"value_{i}" for i in range(1000)}