# -*- coding: utf-8 -*-
r"""
Generalized Tamari lattices

These lattices depend on three parameters `a`, `b` and `m`, where `a`
and `b` are coprime positive integers and `m` is a nonnegative
integer.

The elements are :func:`Dyck paths<sage.combinat.dyck_word.DyckWord>`
in the `(a \times b)`-rectangle. The order relation depends on `m`.

To use the provided functionality, you should import Generalized
Tamari lattices by typing::

    sage: from sage.combinat.tamari_lattices import GeneralizedTamariLattice

Then, ::

    sage: GeneralizedTamariLattice(3,2)
    Finite lattice containing 2 elements
    sage: GeneralizedTamariLattice(4,3)
    Finite lattice containing 5 elements

The classical **Tamari lattices** are special cases of this construction and
are also available directly using the catalogue of posets, as follows::

    sage: posets.TamariLattice(3)
    Finite lattice containing 5 elements

.. SEEALSO::

    For more detailed information see :meth:`TamariLattice`,
    :meth:`GeneralizedTamariLattice`.
"""
# ****************************************************************************
#    Copyright (C) 2012-2018 Frédéric Chapoton <chapoton@math.unistra.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.combinat.posets.lattices import LatticePoset, MeetSemilattice
from sage.arith.all import gcd


def paths_in_triangle(i, j, a, b):
    r"""
    Return all Dyck paths from `(0,0)` to `(i,j)` in the `(a \times
    b)`-rectangle.

    This means that at each step of the path, one has `a y \geq b x`.

    A path is represented by a sequence of `0` and `1`, where `0` is an
    horizontal step `(1,0)` and `1` is a vertical step `(0,1)`.

    INPUT:

    - `a` and `b` -- coprime integers with `a \geq b`

    - `i` and `j` -- nonnegative integers with `1 \geq \frac{j}{b} \geq
      \frac{bi}{a} \geq 0`

    OUTPUT:

    - a list of paths

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import paths_in_triangle
        sage: paths_in_triangle(2,2,2,2)
        [(1, 0, 1, 0), (1, 1, 0, 0)]
        sage: paths_in_triangle(2,3,4,4)
        [(1, 0, 1, 0, 1), (1, 1, 0, 0, 1), (1, 0, 1, 1, 0),
        (1, 1, 0, 1, 0), (1, 1, 1, 0, 0)]
        sage: paths_in_triangle(2,1,4,4)
        Traceback (most recent call last):
        ...
        ValueError: the endpoint is not valid
        sage: paths_in_triangle(3,2,5,3)
        [(1, 0, 1, 0, 0), (1, 1, 0, 0, 0)]
    """
    if not(b >= j and j * a >= i * b and i >= 0):
        raise ValueError("the endpoint is not valid")

    if i == 0:
        return [tuple([1] * j)]

    if (j - 1) * a >= (i) * b:
        result = [u + tuple([1]) for u in paths_in_triangle(i, j - 1, a, b)]
        result += [u + tuple([0]) for u in paths_in_triangle(i - 1, j, a, b)]
        return result

    return [u + tuple([0]) for u in paths_in_triangle(i - 1, j, a, b)]


def swap(p, i, m=1):
    r"""
    Perform a covering move in the `(a,b)`-Tamari lattice of parameter `m`.

    The letter at position `i` in `p` must be a `0`, followed by at
    least one `1`.

    INPUT:

    - ``p`` -- a Dyck path in the `(a \times b)`-rectangle

    - ``i`` -- an integer between `0` and `a+b-1`

    OUTPUT:

    - a Dyck path in the `(a \times b)`-rectangle

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import swap
        sage: swap((1,0,1,0,0),1)
        (1, 1, 0, 0, 0)
        sage: swap((1,1,0,0,1,1,0,0,0),3)
        (1, 1, 0, 1, 1, 0, 0, 0, 0)

    TESTS::

        sage: swap((1,0,1,0),6)
        Traceback (most recent call last):
        ...
        ValueError: the index is greater than the length of the path
        sage: swap((1,1,0,0,1,1,0,0),2)
        Traceback (most recent call last):
        ...
        ValueError: there is no such covering move
    """
    if i >= len(p):
        raise ValueError("the index is greater than the length of the path")
    if i == len(p) - 1 or p[i + 1] == 0:
        raise ValueError("there is no such covering move")

    found = False
    height = 0
    j = i
    while not found and j <= len(p) - 2:
        j += 1
        if p[j]:
            height += m
        else:
            height -= 1
        if height == 0:
            found = True
    q = [k for k in p]
    for k in range(i, j):
        q[k] = p[k + 1]
    q[j] = 0
    return tuple(q)


def GeneralizedTamariLattice(a, b, m=1, check=True):
    r"""
    Return the `(a,b)`-Tamari lattice of parameter `m`.

    INPUT:

    - `a` and `b` -- coprime integers with `a \geq b`

    - `m` -- a nonnegative integer such that `a \geq b m`

    OUTPUT:

    - a finite lattice (the lattice property is only conjectural in general)

    The elements of the lattice are
    :func:`Dyck paths<sage.combinat.dyck_word.DyckWord>` in the
    `(a \times b)`-rectangle.

    The parameter `m` (slope) is used only to define the covering relations.
    When the slope `m` is `0`, two paths are comparable if and only if
    one is always above the other.

    The usual :wikipedia:`Tamari lattice<Tamari_lattice>` of index `b`
    is the special case `a=b+1` and `m=1`.

    Other special cases give the `m`-Tamari lattices studied in [BMFPR]_.

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import GeneralizedTamariLattice
        sage: GeneralizedTamariLattice(3,2)
        Finite lattice containing 2 elements
        sage: GeneralizedTamariLattice(4,3)
        Finite lattice containing 5 elements
        sage: GeneralizedTamariLattice(4,4)
        Traceback (most recent call last):
        ...
        ValueError: the numbers a and b must be coprime with a>=b
        sage: GeneralizedTamariLattice(7,5,2)
        Traceback (most recent call last):
        ...
        ValueError: the condition a>=b*m does not hold
        sage: P = GeneralizedTamariLattice(5,3);P
        Finite lattice containing 7 elements

    TESTS::

        sage: P.coxeter_transformation()**18 == 1
        True

    REFERENCES:

    .. [BMFPR] \M. Bousquet-Melou, E. Fusy, L.-F. Preville Ratelle.
       *The number of intervals in the m-Tamari lattices*. :arxiv:`1106.1498`
    """
    if not(gcd(a, b) == 1 and a >= b):
        raise ValueError("the numbers a and b must be coprime with a>=b")
    if not(a >= b * m):
        raise ValueError("the condition a>=b*m does not hold")

    def covers(p):
        return [swap(p, i, m) for i in range(len(p) - 1)
                if not p[i] and p[i + 1]]
    return LatticePoset({p: covers(p)
                         for p in paths_in_triangle(a, b, a, b)}, check=check)


def TamariLattice(n, m=1):
    r"""
    Return the `n`-th Tamari lattice.

    Using the slope parameter `m`, one can also get the `m`-Tamari lattices.

    INPUT:

    - `n` -- a nonnegative integer (the index)

    - `m` -- an optional nonnegative integer (the slope, default to 1)

    OUTPUT:

    a finite lattice

    In the usual case, the elements of the lattice are :func:`Dyck
    paths<sage.combinat.dyck_word.DyckWord>` in the `(n+1 \times
    n)`-rectangle. For a general slope `m`, the elements are Dyck
    paths in the `(m n+1 \times n)`-rectangle.

    See :wikipedia:`Tamari lattice<Tamari_lattice>` for mathematical
    background.

    EXAMPLES::

        sage: posets.TamariLattice(3)
        Finite lattice containing 5 elements

        sage: posets.TamariLattice(3, 2)
        Finite lattice containing 12 elements

    REFERENCES:

    - [BMFPR]_
    """
    return GeneralizedTamariLattice(m * n + 1, n, m, check=False)


# a variation : the Dexter meet-semilattices


def swap_dexter(p, i):
    r"""
    Perform covering moves in the `(a,b)`-Dexter posets.

    The letter at position `i` in `p` must be a `0`, followed by at
    least one `1`.

    INPUT:

    - ``p`` -- a Dyck path in the `(a \times b)`-rectangle

    - ``i`` -- an integer between `0` and `a+b-1`

    OUTPUT:

    - a list of Dyck paths in the `(a \times b)`-rectangle

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import swap_dexter
        sage: swap_dexter((1,0,1,0,0),1)
        [(1, 1, 0, 0, 0)]
        sage: swap_dexter((1,1,0,0,1,1,0,0,0),3)
        [(1, 1, 0, 1, 1, 0, 0, 0, 0), (1, 1, 1, 1, 0, 0, 0, 0, 0)]
        sage: swap_dexter((1,1,0,1,0,0,0),2)
        []

    TESTS::

        sage: swap_dexter((1,0,1,0,0),6)
        Traceback (most recent call last):
        ...
        ValueError: the index is greater than the length of the path

        sage: swap_dexter((1,1,0,0,1,1,0,0,0),2)
        Traceback (most recent call last):
        ...
        ValueError: there is no such covering move
    """
    m = 1

    if i >= len(p):
        raise ValueError("the index is greater than the length of the path")
    if i == len(p) - 1 or p[i + 1] == 0:
        raise ValueError("there is no such covering move")

    found = False
    height = 0
    j = i
    while not found and j <= len(p) - 2:
        j += 1
        if p[j]:
            height += m
        else:
            height -= 1
        if height == 0:
            found = True
    if not (j == len(p) - 2 or p[j + 1]):  # forbidden moves
        return []

    resu = []
    tp = tuple(p)
    for deb in range(i, 0, -1):
        if not p[deb]:
            q = tp[:deb] + tp[i + 1: j + 1] + tp[deb: i + 1] + tp[j + 1:]
            resu.append(q)
        else:
            break
    return resu


def DexterSemilattice(n):
    r"""
    Return the `n`-th Dexter meet-semilattice.

    INPUT:

    - ``n`` -- a nonnegative integer (the index)

    OUTPUT:

    a finite meet-semilattice

    The elements of the semilattice are :func:`Dyck
    paths<sage.combinat.dyck_word.DyckWord>` in the `(n+1 \times
    n)`-rectangle.

    EXAMPLES::

        sage: posets.DexterSemilattice(3)
        Finite meet-semilattice containing 5 elements

        sage: P = posets.DexterSemilattice(4); P
        Finite meet-semilattice containing 14 elements
        sage: len(P.maximal_chains())
        15
        sage: len(P.maximal_elements())
        4
        sage: P.chain_polynomial()
        q^5 + 19*q^4 + 47*q^3 + 42*q^2 + 14*q + 1

    REFERENCES:

    - [Cha18]_
    """
    a = n + 1
    b = n

    def covers_dexter(p):
        data = [swap_dexter(p, i) for i in range(len(p) - 1)
                if not p[i] and p[i + 1]]
        return [cov for L in data for cov in L]
    return MeetSemilattice({p: covers_dexter(p)
                            for p in paths_in_triangle(a, b, a, b)},
                           check=False)
