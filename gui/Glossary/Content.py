# There is likely a better solution to this, but for now this is how content will be stored
def outputContent(topic):
    if topic == "GT":
        return {"Graph" : "A graph G is a pair of sets (V,E) where the members of V are called vertices"
                +"\nand the members of E are unordered pairs of distinct vertices called edges." 
                "\nWe say two vertices are adjacent if they are the elements of an edge. We \nthen say the edge is incident to those vertices.",
                "Degree":"For a vertex x, the set of vertices adjacent to x is denoted N(x). We say the \ndegree of x is the size of N(x), denoted |N(x)|.",
                "Subgraph":"A graph S is a subgraph of a graph G if V(S) ⊆ V(G) and E(S) ⊆ E(G).",
                "Complete Graph":"A complete graph on n vertices, denoted K(n), is a graph on n vertices \nsuch that every vertex is adjacent to every other vertex.",
                "Walk":"A walk is an alternating sequence of vertices and edges.",
                "Path":"A path is a walk where all vertices are distinct.",
                "Trail":"If the edges in a walk are distinct, the walk is called a trail.",
                "Curcuit":"A closed trail that is one where end vertices coincide, is called a circuit.",
                "Hamiltonian":"A graph is hamiltonian if there is a cycle that traverses ever vertex exactly \nonce.",
                "Eularian":"A graph is Eularian if there is a circuit that contains all the edges of G.",
                "Planar":"A graph is planar if it can be drawn on the plane with no edge crossings.",
                "Subdivision":"A graph G is a subdivision of a graph H if edges in H are replaced with \npaths where two such paths can have only their end vertices in common.",
                "Connected":"A graph is connected if for ever pair of distinct vertices, there is a path \njoining those vertices.",
                "Handshaking Lemma":"Let G be a simple graph where a,b,c,...,x,y,z are vertices of G, \nthen |N(a)|+|N(b)|+|N(c)|+...+|N(x)|+|N(y)|+|N(z)|=2|E(G)|",
                "Euler's Formula":"If G is a connected planar graph with n vertices, m edges, and f faces, \nthen n-m+f=2.",
                "Bipartite Graph":"A graph is bipartite if the vertex set can be partitioned into two disjoint \nsets A and B and for every edge e={x,y}, we have that x is in A and y is in B.",
                "Kuratowski's Theorem":"A graph is planar if and only if it has no subgraph that is a subdivision of \nK(5) or K(3,3)"}
    elif topic == "LA":
        return {"Linear Map": "",
                 "Vector space": "",
                 "Orthogonality":""}
    else:
        return {f"key_{i}": f"value_{i}" for i in range(1000)}