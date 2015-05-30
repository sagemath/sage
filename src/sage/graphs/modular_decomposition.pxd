cdef extern from "modular_decomposition.h":

    ctypedef struct noeud:
        pass

    ctypedef struct c_adj "adj":

        int s # number of the vertex
        c_adj * suiv # adress of next pair in the list, NULL if last

    ctypedef struct c_fils "fils":

        noeud * pointe # adress of the node in the tree
        c_fils * suiv # adress of the next pair in the list, NULL if last

    ctypedef struct c_noeud "noeud":

        int type # is FEUILLE, SERIE, PARALLELE or PREMIER
        int ps # internal use
        int bg # internal use
        int ds # internal use
        int bd # internal use
        int sommet # internal use
        int nom # if type=FEUILLE, number of the corresponding vertex of the graph
        int id # internal use (node unique ID)

        c_noeud *pere # adress of parent node, NULL if root
        c_fils *fpere # points the head of the linked list of sons (if type is not FEUILLE, else is NULL)
        c_fils *fils # points the head of the linked list of sons
        c_fils *lastfils # internal use (points the last item in the listed list of sons)

    ctypedef struct c_graphe "graphe":

        int n
        c_adj ** G

    c_noeud *decomposition_modulaire(c_graphe G)
