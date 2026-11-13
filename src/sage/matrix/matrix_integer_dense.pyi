from sage.matrix.matrix_dense import Matrix_dense

class Matrix_integer_dense(Matrix_dense):
    def LLL(self, delta=None, eta=None, algorithm='fpLLL:wrapper', fp=None, prec=0, early_red=False, use_givens=False, use_siegel=False, transformation=False, **kwds) -> Matrix_integer_dense: ...
