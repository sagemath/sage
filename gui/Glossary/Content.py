# There is likely a better solution to this, but for now this is how content will be stored
def outputContent(topic):
    if topic == "GT":
        return {"Graph" : "A " + '{:s}'.format('\u0332'.join('Graph'))+ " G is a pair of sets (V,E) where the members of V are called "+ '{:}'.format('\u0333'.join('vertices'))+
                "\nand the members of E are unordered pairs of distinct vertices called "+'{:s}'.format('\u0332'.join('edges'))+"." + 
                "\nWe say two vertices are adjacent if they are the elements of an edge. We then say the \nedge is incident to those vertices.",
                "Degree":"For a vertex x, the set of vertices adjacent to x is denoted N(x). We say the \ndegree of x is the size of N(x).",
                "Subgraph":"A graph S is a subgraph of a graph G if V(S) ⊆ V(G) and E(S) ⊆ E(G).",
                "Complete Graph":"A complete graph on n vertices, denoted K(n), is a graph on n vertices \nsuch that every vertex is adjacent to every other vertex.",
                "Walk":"A walk is an alternating sequence of vertices and edges.",
                "Path":"A path is a walk where all vertices are distinct.",
                "Trail":"If the edges in a walk are distinct, the walk is called a trail.",
                "Curcuit":"A closed trail that is one where end vertices coincide, is called a circuit.",
                "Hamiltonian":"A graph is hamiltonian if there is a cycle that traverses ever vertex exactly once.",
                "Eularian":"A graph is Eularian if there is a circuit that contains all the edges of G"}
    elif topic == "LA":
        return {"Matrix": "stuff", "Vector space": "other stuff", "Orthogonality":"some definition involving inner products and other stuff like that"}
    else:
        return {f"key_{i}": f"value_{i}" for i in range(1000)}