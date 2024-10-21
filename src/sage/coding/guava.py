# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.libs.gap sage.modules sage.rings.finite_rings
r"""
Constructions of generator matrices using the GUAVA package for GAP

This module only contains Guava wrappers (GUAVA is an optional GAP package).

AUTHORS:

- David Joyner (2005-11-22, 2006-12-03): initial version

- Nick Alexander (2006-12-10): factor GUAVA code to guava.py

- David Joyner (2007-05): removed Golay codes, toric and trivial codes and
  placed them in code_constructions; renamed RandomLinearCode to
  RandomLinearCodeGuava

- David Joyner (2008-03): removed QR, XQR, cyclic and ReedSolomon codes

- David Joyner (2009-05): added "optional package" comments, fixed some
  docstrings to be sphinx compatible

- Dima Pasechnik (2019-11): port to libgap
"""
#*****************************************************************************
#       Copyright (C) 2007 David Joyner <wdj@usna.edu>
#                     2006 Nick Alexander <ncalexan@math.uci.edu>
#                     2019 Dima Pasechnik <dima@pasechnik.info>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gap.libgap import libgap
from sage.misc.randstate import current_randstate
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from .linear_code import LinearCode
from sage.features.gap import GapPackage


def QuasiQuadraticResidueCode(p):
    r"""
    A (binary) quasi-quadratic residue code (or QQR code).

    Follows the definition of Proposition 2.2 in [BM2003]_. The code has a generator
    matrix in the block form `G=(Q,N)`. Here `Q` is a `p \times p` circulant
    matrix whose top row is `(0,x_1,...,x_{p-1})`, where `x_i=1` if and only if
    `i` is a quadratic residue mod `p`, and `N` is a `p \times p` circulant
    matrix whose top row is `(0,y_1,...,y_{p-1})`, where `x_i+y_i=1` for all
    `i`.

    INPUT:

    - ``p`` -- a prime `>2`

    OUTPUT: a QQR code of length `2p`

    EXAMPLES::

        sage: C = codes.QuasiQuadraticResidueCode(11); C   # optional - gap_package_guava
        [22, 11] linear code over GF(2)

    These are self-orthogonal in general and self-dual when `p \equiv 3 \pmod 4`.

    AUTHOR: David Joyner (11-2005)
    """
    GapPackage('guava', spkg='gap_packages').require()
    libgap.load_package('guava')
    C = libgap.QQRCode(p)
    G = C.GeneratorMat()
    MS = MatrixSpace(GF(2), len(G), len(G[0]))
    return LinearCode(MS(G))


def RandomLinearCodeGuava(n, k, F):
    r"""
    The method used is to first construct a `k \times n` matrix of the block
    form `(I,A)`, where `I` is a `k \times k` identity matrix and `A` is a
    `k \times (n-k)` matrix constructed using random elements of `F`. Then
    the columns are permuted using a randomly selected element of the symmetric
    group `S_n`.

    INPUT:

    - ``n``, ``k`` -- integers with `n>k>1`

    OUTPUT: a "random" linear code with length `n`, dimension `k` over field `F`

    EXAMPLES::

        sage: C = codes.RandomLinearCodeGuava(30,15,GF(2)); C      # optional - gap_package_guava
        [30, 15] linear code over GF(2)
        sage: C = codes.RandomLinearCodeGuava(10,5,GF(4,'a')); C   # optional - gap_package_guava
        [10, 5] linear code over GF(4)

    AUTHOR: David Joyner (11-2005)
    """
    current_randstate().set_seed_gap()

    GapPackage('guava', spkg='gap_packages').require()
    libgap.load_package('guava')
    C = libgap.RandomLinearCode(n,k,F)
    G = C.GeneratorMat()
    MS = MatrixSpace(F, len(G), len(G[0]))
    return LinearCode(MS(G))
