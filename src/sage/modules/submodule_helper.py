r"""
Helper classes for handling submodules over various base rings.

Currently, it is only used for Ore modules.

AUTHOR:

- Xavier Caruso (2025-08)
"""

# ***************************************************************************
#    Copyright (C) 2025 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.fields import Fields
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.matrix.matrix_polynomial_dense import Matrix_polynomial_dense
from sage.matrix.constructor import matrix


class SubmoduleHelper(metaclass=ClasscallMetaclass):
    r"""
    A class for manipulating submodules at the level of matrices.
    The class provides the arguments:

    - ``basis``: a matrix in normal form whose rows form a
      basis of the submodule

    - ``complement``: a matrix in normal form whose rows form
      a basis of a complement of the submodule

    - ``coordinates``: a change-of-basis matrix from the canonical
      basis to the basis obtained by concatening ``basis`` and
      ``complement``

    - ``is_saturated``: a boolean; whether this submodule is
      saturated in the ambient space
    """
    def __classcall_private__(self, mat, saturate=False):
        r"""
        Dispatch to the appropriate class.

        TESTS::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: Fq = GF(5)
            sage: SH = SubmoduleHelper(matrix(Fq, [[1, 3]]))
            sage: type(SH)
            <class 'sage.modules.submodule_helper.SubmoduleHelper_field'>

        ::

            sage: A.<t> = Fq[]
            sage: SH = SubmoduleHelper(matrix(A, [[t, 3 + t]]))
            sage: type(SH)
            <class 'sage.modules.submodule_helper.SubmoduleHelper_polynomial_ring'>

        ::

            sage: SH = SubmoduleHelper(matrix(ZZ, [[1, 3]]))
            sage: type(SH)
            <class 'sage.modules.submodule_helper.SubmoduleHelper_PID'>

        ::

            sage: R.<x,y> = Fq[]
            sage: SH = SubmoduleHelper(matrix(R, [[x, 3 + y]]))
            Traceback (most recent call last):
            ...
            NotImplementedError: submodules and quotients are only implemented over PIDs
        """
        base = mat.base_ring()
        if base in Fields():
            cls = SubmoduleHelper_field
        elif (isinstance(mat, Matrix_polynomial_dense)
          and base.base_ring() in Fields()):
            cls = SubmoduleHelper_polynomial_ring
        elif base in PrincipalIdealDomains():
            cls = SubmoduleHelper_PID
        else:
            raise NotImplementedError("submodules and quotients are only implemented over PIDs")
        return cls.__call__(mat, saturate)

    def __hash__(self):
        r"""
        Return a hash of this submodule.

        Two differents instances corresponding to the same submodule
        have identical hashes.

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: Fq = GF(5)
            sage: SH1 = SubmoduleHelper(matrix(Fq, [[1, 3]]))
            sage: h1 = hash(SH1)
            sage: h1  # random
            -8066549207035083627

        ::

            sage: Fq = GF(5)
            sage: SH2 = SubmoduleHelper(matrix(Fq, [[2, 6]]))
            sage: h2 = hash(SH2)
            sage: h2  # random
            -8066549207035083627
            sage: h1 == h2
            True
        """
        return hash(self.basis)

    def __eq__(self, other):
        r"""
        Return ``True`` is ``self`` and ``other`` define the
        same submodule.

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: Fq = GF(5)
            sage: SH1 = SubmoduleHelper(matrix(Fq, [[1, 3]]))
            sage: SH2 = SubmoduleHelper(matrix(Fq, [[2, 6]]))
            sage: SH1 == SH2
            True

        ::

            sage: SHZ = SubmoduleHelper(matrix(ZZ, [[1, 3]]))
            sage: SHZ == SH1
            False
        """
        sb = self.basis
        ob = other.basis
        return sb.parent() is ob.parent() and sb == ob


class SubmoduleHelper_field(SubmoduleHelper):
    r"""
    Submodules over fields.
    """
    def __init__(self, mat, saturate):
        r"""
        Initialize this submodule.

        INPUT:

        - ``mat`` -- a matrix whose rows span the submodule

        - ``saturate`` -- a boolean, ignored

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: Fq = GF(5)
            sage: SH1 = SubmoduleHelper(matrix(Fq, [[1, 2, 3], [1, 3, 4]]))
            sage: SH1.basis
            [1 0 1]
            [0 1 1]
            sage: SH1.complement
            [0 0 1]
            sage: SH1.coordinates
            [1 0 4]
            [0 1 4]
            [0 0 1]

        We check that the previous outputs are coherent::

            sage: SH1.basis.stack(SH1.complement) * SH1.coordinates == 1
            True
        """
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


class SubmoduleHelper_PID(SubmoduleHelper):
    r"""
    Submodules over principal ideal domains (except
    polynomial rings to which a special class is dedicated).
    """
    def __init__(self, mat, saturate):
        r"""
        Initialize this submodule.

        INPUT:

        - ``mat`` -- a matrix whose rows span the submodule

        - ``saturate`` -- a boolean

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: SH = SubmoduleHelper(matrix(ZZ, [[1, 2, 3], [1, 3, 4]]))
            sage: SH.basis
            [1 0 1]
            [0 1 1]
            sage: SH.complement
            [1 0 0]
            sage: SH.coordinates
            [ 0  0  1]
            [-1  1  1]
            [ 1  0 -1]

        We check that the previous outputs are coherent::

            sage: SH.basis.stack(SH.complement) * SH.coordinates == 1
            True

        When ``saturate=True``, the saturation of the span of the given
        matrix is created::

            sage: SH1 = SubmoduleHelper(matrix(ZZ, [[1, 2, 3], [2, 6, 8]]))
            sage: SH1.basis
            [1 0 1]
            [0 2 2]
            sage: SH2 = SubmoduleHelper(matrix(ZZ, [[1, 2, 3], [2, 6, 8]]), True)
            sage: SH2.basis
            [1 0 1]
            [0 1 1]
        """
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
            basis = W.matrix_from_rows(range(r))
            complement = W.matrix_from_rows(range(r, n))
            self.is_saturated = True
        else:
            S = S.matrix_from_rows(range(r))
            basis = matrix(base, [[S[i,i]*W[i,j] for j in range(n)]
                                        for i in range(r)])
            complement = W.matrix_from_rows(range(r, n))
            self.is_saturated = all(S[i,i].is_unit() for i in range(r))
        self.basis = basis.echelon_form()
        self.basis.set_immutable()
        self.complement = complement.echelon_form()
        self.coordinates = self.basis.stack(self.complement).inverse()


class SubmoduleHelper_polynomial_ring(SubmoduleHelper):
    r"""
    Submodules over polynomial rings.
    """
    def __init__(self, mat, saturate):
        r"""
        Initialize this submodule.

        INPUT:

        - ``mat`` -- a matrix whose rows span the submodule

        - ``saturate`` -- a boolean

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: A.<t> = GF(5)[]
            sage: SH = SubmoduleHelper(matrix(A, [[t, t^2, t^3], [t, t+1, t+2]]))
            sage: SH.basis
            [    t^3 + 3*t^2 + 3*t t^3 + 3*t^2 + 2*t + 4                     3]
            [                    t                 t + 1                 t + 2]

        When ``saturate=True``, the saturation of the span of the given
        matrix is created::

            sage: SHsat = SubmoduleHelper(matrix(A, [[t, t^2, t^3], [t, t+1, t+2]]), True)
            sage: SHsat.basis
            [t^2 + 3*t + 4 t^2 + 3*t + 3             1]
            [            t         t + 1         t + 2]
        """
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
        r"""
        Return a Popov form and the corresponding transformation matrix
        of the matrix whose rows are the concatenation of the rows of
        ``self.basis`` and ``self.complement``.

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: A.<t> = GF(5)[]
            sage: SH = SubmoduleHelper(matrix(A, [[1, t, t^2], [t, t+1, t+2]]))
            sage: SH._popov
            (
            [1 0 0]  [2*t^3 + 4*t^2 + 4*t + 4 2*t^3 + 4*t^2 + 4*t + 3             3*t^3 + 4*t]
            [0 1 0]  [  3*t^3 + t^2 + 4*t + 1   3*t^3 + t^2 + 4*t + 2         2*t^3 + 3*t + 1]
            [0 0 1], [        2*t^2 + 2*t + 2         2*t^2 + 2*t + 2         3*t^2 + 2*t + 2]
            )
        """
        return self.basis.stack(self.complement).popov_form(transformation=True)

    @lazy_attribute
    def complement(self):
        r"""
        Return a basis of a complement submodule of this submodule.

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: A.<t> = GF(5)[]
            sage: SH = SubmoduleHelper(matrix(A, [[1, t, t^2], [t, t+1, t+2]]))
            sage: SH.complement
            [t^2 + t + 1 t^2 + t + 1           t]
        """
        return self.basis.basis_completion().popov_form()

    @lazy_attribute
    def coordinates(self):
        r"""
        Return a change-of-basis matrix from the canonical basis
        to the basis obtained by concatening ``self.basis`` and
        ``self.complement``.

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: A.<t> = GF(5)[]
            sage: SH = SubmoduleHelper(matrix(A, [[1, t, t^2], [t, t+1, t+2]]))
            sage: SH.coordinates
            [2*t^3 + 4*t^2 + 4*t + 4 2*t^3 + 4*t^2 + 4*t + 3             3*t^3 + 4*t]
            [  3*t^3 + t^2 + 4*t + 1   3*t^3 + t^2 + 4*t + 2         2*t^3 + 3*t + 1]
            [        2*t^2 + 2*t + 2         2*t^2 + 2*t + 2         3*t^2 + 2*t + 2]

        ::

            sage: SH.basis.stack(SH.complement) * SH.coordinates == 1
            True
        """
        P, T = self._popov
        if P.is_one():
            return T
        else:
            return self.basis.stack(self.complement).inverse()

    @lazy_attribute
    def is_saturated(self):
        r"""
        Return ``True`` if this submodule is saturated; ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.submodule_helper import SubmoduleHelper
            sage: A.<t> = GF(5)[]
            sage: SH = SubmoduleHelper(matrix(A, [[1, t, t^2], [t, t+1, t+2]]))
            sage: SH.is_saturated
            True

        ::

            sage: SH = SubmoduleHelper(matrix(A, [[t, t^2, t^3], [t, t+1, t+2]]))
            sage: SH.is_saturated
            False
        """
        P, _ = self._popov
        return P.is_one()
