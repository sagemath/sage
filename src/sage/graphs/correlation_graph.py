from sage.graphs.all import Graph
from numpy import corrcoef
from sage.all import Matrix

class CorrelationGraph(Graph):
    """
    TODO: docstring
    EXAMPLES:

        sage: data=[[1,2,3],[4,5,6],[7,8,9999]]
        sage: CG = CorrelationGraph(data, 0.9)
        sage: CG
        Looped graph on 3 vertices
    """
    
    def __init__(self, seqs, alpha, include_anticorrelation=False):

        # compute correlations
        corrs = corrcoef(seqs)

        # compare against alpha to get adjacency matrix
        if include_anticorrelation:
            boolean_adjacency_matrix = (abs(corrs)>=alpha)
        else:
            boolean_adjacency_matrix = (corrs>=alpha)
        adjacency_matrix=Matrix(boolean_adjacency_matrix.astype(int))

        # call graph constructor
        super().__init__(adjacency_matrix, format="adjacency_matrix")    
    
    
    def plot(self):

        # hack, hopefully temporary, see trac for discussion
        self.__class__=Graph
        p=self.plot()
        self.__class__=CorrelationGraph

        return p