import pytest

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import random_vector
from sage.rings.cc import CC
from sage.rings.complex_double import CDF
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.qqbar import AA, QQbar
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RR
from sage.symbolic.ring import SR

QF = QuadraticField(-1, 'I')

# A list of suitable rings to test. The first component is the ring
# itself, the second is a friendly name, and the third component is
# the maximum size of the matrices we'll generate over that ring.
rings = [(ZZ,    "ZZ",    5),
         (QQ,    "QQ",    5),
         (AA,    "AA",    2),  # Can be ungodly slow
         (QQbar, "QQbar", 2),  # Can be ungodly slow
         (QF,    "QQ(I)", 5),
         (RDF,   "RDF",   5),
         (RR,    "RR",    5),
         (CDF,   "CDF",   5),
         (CC,    "CC",    5),
         (SR,    "SR",    5)]

# Not chosen rigorously, but a Hilbert matrix of order 5
# factors quite nicely, so what could possibly go wrong?
tolerance = 1e-6

# Generate an id from a (ring, name, size) triple.
id_gen = lambda s: f"{s[1]}^{s[2]}x{s[2]}"

A_and_I_params = ( (r[0],r[1],n)
                   for r in rings
                   for n in range(r[2]+1) )
@pytest.fixture(params=A_and_I_params, ids=id_gen)
def A_and_I(request):
    r"""
    A generic random Hermitian matrix over any base ring, of any size,
    paired with the corresponding identity matrix.
    """
    R,_,n = request.param
    _A = matrix.random(R,n)
    _I = matrix.identity(R,n)
    return (_A + _A.H, _I)

@pytest.fixture
def A(A_and_I):
    return A_and_I[0]

@pytest.fixture
def I(A_and_I):
    return A_and_I[1]


# SR is inexact, but it's a special case because most arithmetic can
# be done exactly, and new users are predisposed to use it. Its
# random_element() method is relatively tame, otherwise we wouldn't be
# able to lump it in with the others.
B_params = ( (r[0],r[1],n)
             for r in rings
             for n in range(r[2]+1)
             if r[0].is_exact() or r[0] is SR)
@pytest.fixture(params=B_params, ids=id_gen)
def B(request):
    r"""
    A random positive-definite Hermitian matrix over an
    exact/symbolic base ring, of any size.
    """
    R,_,n = request.param
    _B = matrix.random(R,n)
    return (_B * _B.H) + matrix.identity(R,n)


def test_identity_factorization_is_trivial(I):
    r"""
    All three factors should be the identity when the input
    matrix is.
    """
    assert(I.block_ldlt() == (I,I,I))

def test_diagonal_factorization_is_trivial(A_and_I):
    r"""
    The factorization should be (essentially) a no-op on diagonal
    matrices.
    """
    A,I = A_and_I
    A = matrix.diagonal(A.diagonal())
    assert(A.block_ldlt() == (I,I,A))

def test_reconstruction_is_accurate(A):
    r"""
    Ensure that a random Hermitian matrix can be reconstructed
    from its factorization.
    """
    P,L,D = A.block_ldlt()
    difference = P.T*A*P - L*D*L.H
    if A.base_ring().is_exact() or A.base_ring() is SR:
        assert(difference.is_zero())
    else:
        assert(difference.norm() < tolerance)

def test_base_rings_of_factors_are_unchanged(A):
    r"""
    All three factors have the same base ring, even when they're
    trivial.
    """
    P,L,D = A.block_ldlt()
    assert(P.base_ring() == L.base_ring())
    assert(L.base_ring() == D.base_ring())

def test_diagonal_is_1x1_for_positive_definite(B):
    r"""
    The diagonal blocks should be 1x1 if the input is
    positive-definite and if the greedy algorithm is used.

    We only test exact/symbolic rings because this property isn't
    guaranteed with Bunch-Kaufman, and the greedy algorithm may not
    work well in other rings.
    """
    P,L,D = B.block_ldlt(algorithm="greedy")
    assert(matrix.diagonal(D.diagonal()) == D)

def test_classical_and_indefinite_factorizations_agree(B):
    r"""
    For a positive-definite matrix, the classical `LDL^{T}` and
    indefinite factorizations should agree.

    We skip over ``SR`` because ``indefinite_factorization()``
    requires an exact ring.
    """
    if B.base_ring() is SR:
        return
    _,L,D = B.block_ldlt(algorithm="classical")
    l,d = B.indefinite_factorization(algorithm='hermitian')
    assert(L == l)
    assert(D == matrix.diagonal(d))
