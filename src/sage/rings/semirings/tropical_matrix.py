from sage.rings.infinity import infinity
from sage.matrix.constructor import matrix
from sage.matrix.matrix_generic_dense import Matrix_generic_dense

class Matrix_tropical_dense(Matrix_generic_dense):
    def extremum_mean_weight(self):
        # Karp algorithm
        T = self.base_ring()
        n = self.ncols()
        if self.nrows() != n:
            raise TypeError("matrix must be square")
        v = matrix(1, n, n*[T.one()])  # ???
        vs = [v]
        for _ in range(n):
            v = v * self
            vs.append(v)
        w = [vs[n][0,j].lift() for j in range(n)]
        if T._use_min:
            return min(max((w[j] - vs[k][0,j].lift()) / (n-k) for k in range(n))
                       for j in range(n) if w[j] is not infinity)
        else:
            return max(min((w[j] - vs[k][0,j].lift()) / (n-k) for k in range(n))
                       for j in range(n) if w[j] is not infinity)

    def weak_transitive_closure(self):
        # Floyd-Warshall algorithm
        T = self.base_ring()
        n = self.ncols()
        if self.nrows() != n:
            raise TypeError("matrix must be square")
        G = self.__copy__()
        for p in range(n):
            for i in range(n):
                if i == p:
                    continue
                for j in range(n):
                    if j == p:
                        continue
                    G[i,j] += G[i,p] * G[p,j]
                    if i == j:
                        if T._use_min and G[i,i].lift() < 0:
                            raise ValueError("negative cycle exists")
                        if not T._use_min and G[i,i].lift() > 0:
                            raise ValueError("positive cycle exists")
        return G

    def strong_transitive_closure(self):
        return self.parent().identity_matrix() + self.weak_transitive_closure()

    kleene_star = strong_transitive_closure
