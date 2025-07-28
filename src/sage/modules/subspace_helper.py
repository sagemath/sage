from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.fields import Fields
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.matrix.matrix_polynomial_dense import Matrix_polynomial_dense
from sage.matrix.constructor import matrix


class SubspaceHelper(metaclass=ClasscallMetaclass):
    def __classcall_private__(self, mat, saturate=False):
        base = mat.base_ring()
        if base in Fields():
            cls = SubspaceHelper_field
        elif (isinstance(mat, Matrix_polynomial_dense)
          and base.base_ring() in Fields()):
            cls = SubspaceHelper_polynomial_ring
        elif base in PrincipalIdealDomains():
            cls = SubspaceHelper_PID
        else:
            raise NotImplementedError("subspaces and quotients are only implemented over PID")
        return cls.__call__(mat, saturate)

    def __hash__(self):
        return hash(self.basis)

    def __eq__(self, other):
        return self.basis == other.basis

    def __repr__(self):
        return self.basis.__repr__()


class SubspaceHelper_field(SubspaceHelper):
    def __init__(self, mat, saturate):
        base = mat.base_ring()
        n = mat.ncols()
        basis = mat.echelon_form()
        self.rank = r = basis.rank()
        pivots = basis.pivots()
        self.basis = basis.matrix_from_rows(range(r))
        self.basis.set_immutable()
        self.complement = matrix(base, n-r, n)
        self.coordinates = matrix(base, n, n)
        indices = []
        i = 0
        for j in range(n):
            if i < r and pivots[i] == j:
                self.coordinates[j, i] = base.one()
                i += 1
            else:
                indices.append(j)
                self.complement[j-i, j] = base.one()
                self.coordinates[j, j-i+r] = base.one()
        for i in range(r):
            for j in range(n-r):
                self.coordinates[pivots[i], j+r] = -basis[i, indices[j]]
        self.is_saturated = True


class SubspaceHelper_PID(SubspaceHelper):
    def __init__(self, mat, saturate):
        base = mat.base_ring()
        n = mat.ncols()
        S, U, V = mat.smith_form()
        r = 0
        for i in range(min(S.nrows(), S.ncols())):
            if S[i,i] == 0:
                break
            r += 1
        self.rank = r
        W = V.inverse().change_ring(base)
        if saturate:
            self.basis = W.matrix_from_rows(range(r))
            self.complement = W.matrix_from_rows(range(r, n))
            self.coordinates = V
            self.is_saturated = True
        else:
            S = S.matrix_from_rows(range(r))
            self.basis = matrix(base, [[S[i,i]*W[i,j] for j in range(n)]
                                        for i in range(r)])
            self.complement = W.matrix_from_rows(range(r, n))
            K = base.fraction_field()
            scalars = [~S[i,i] for i in range(r)] + (n-r)*[K.one()]
            self.coordinates = matrix(K, [[scalars[j]*V[i,j] for j in range(n)]
                                           for i in range(n)])
            self.is_saturated = all(S[i,i].is_unit() for i in range(r))
        self.basis.set_immutable()


class SubspaceHelper_polynomial_ring(SubspaceHelper):
    def __init__(self, mat, saturate):
        base = mat.base_ring()
        if saturate:
            S, _, V = mat.smith_form()
            W = V.inverse().change_ring(base)
            r = 0
            for i in range(min(S.nrows(), S.ncols())):
                if S[i,i] == 0:
                    break
                r += 1
            mat = W.matrix_from_rows(range(r))
            self.is_saturated = True
        self.basis = mat.popov_form(include_zero_vectors=False)
        self.basis.set_immutable()
        self.rank = self.basis.nrows()

    @lazy_attribute
    def _popov(self):
        return self.basis.stack(self.complement).popov_form(transformation=True)

    @lazy_attribute
    def complement(self):
        return self.basis.basis_completion().popov_form()

    @lazy_attribute
    def coordinates(self):
        P, T = self._popov
        if P.is_one():
            return T
        else:
            return self.basis.stack(self.complement).inverse()

    @lazy_attribute
    def is_saturated(self):
        P, _ = self._popov
        return P.is_one()
