# sage.doctest: needs sage.modules
r"""
Generalized Tamari lattices

These lattices depend on three parameters `a`, `b` and `m`, where `a`
is a positive integer, `b` is a nonnegative integer, and `m` is a
nonnegative rational number.

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
from __future__ import annotations
from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.posets.lattices import LatticePoset, MeetSemilattice
from sage.rings.finite_rings.integer_mod import IntegerMod_abstract
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


def _coerce_nonnegative_integer(value, name):
    r"""
    Return ``value`` as a nonnegative integer.

    INPUT:

    - ``value`` -- the value to coerce

    - ``name`` -- string; the parameter name to use in error messages

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import _coerce_nonnegative_integer
        sage: _coerce_nonnegative_integer(QQ(3), 'n')
        3
        sage: _coerce_nonnegative_integer(RDF(3), 'n')                                  # needs sage.rings.real_double
        3
        sage: _coerce_nonnegative_integer(Zmod(2)(1), 'n')                              # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: n must be a nonnegative integer
        sage: _coerce_nonnegative_integer(3/2, 'n')
        Traceback (most recent call last):
        ...
        ValueError: n must be a nonnegative integer
    """
    if isinstance(value, IntegerMod_abstract):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        value = ZZ(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a nonnegative integer") from None
    if value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _coerce_positive_integer(value, name):
    r"""
    Return ``value`` as a positive integer.

    INPUT:

    - ``value`` -- the value to coerce

    - ``name`` -- string; the parameter name to use in error messages

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import _coerce_positive_integer
        sage: _coerce_positive_integer(QQ(3), 'a')
        3
        sage: _coerce_positive_integer(RDF(3), 'a')                                     # needs sage.rings.real_double
        3
        sage: _coerce_positive_integer(Zmod(2)(1), 'a')                                 # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: a must be a positive integer
        sage: _coerce_positive_integer(0, 'a')
        Traceback (most recent call last):
        ...
        ValueError: a must be a positive integer
    """
    try:
        value = _coerce_nonnegative_integer(value, name)
    except ValueError:
        raise ValueError(f"{name} must be a positive integer") from None
    if value == 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _coerce_nonnegative_rational(value, name):
    r"""
    Return ``value`` as a nonnegative rational number.

    INPUT:

    - ``value`` -- the value to coerce

    - ``name`` -- string; the parameter name to use in error messages

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import _coerce_nonnegative_rational
        sage: _coerce_nonnegative_rational(RDF(2)/5, 'm') == QQ(RDF(2)/5)               # needs sage.rings.real_double
        True
        sage: from fractions import Fraction
        sage: _coerce_nonnegative_rational(Fraction('1/2'), 'm')
        1/2
        sage: _coerce_nonnegative_rational(Zmod(2)(1), 'm')                             # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: m must be a nonnegative rational number
    """
    if isinstance(value, IntegerMod_abstract):
        raise ValueError(f"{name} must be a nonnegative rational number")
    try:
        value = QQ(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a nonnegative rational number") from None
    if value < 0:
        raise ValueError(f"{name} must be a nonnegative rational number")
    return value


def paths_in_triangle(i, j, a, b) -> list[tuple[int, ...]]:
    r"""
    Return all Dyck paths from `(0,0)` to `(i,j)` in the `(a \times
    b)`-rectangle.

    This means that at each step of the path, one has `a y \geq b x`.

    A path is represented by a sequence of `0` and `1`, where `0` is an
    horizontal step `(1,0)` and `1` is a vertical step `(0,1)`.

    INPUT:

    - ``a``, ``b`` -- integers with `a \geq b`

    - ``i``, ``j`` -- nonnegative integers with `1 \geq \frac{j}{b} \geq
      \frac{i}{a} \geq 0`

    OUTPUT: list of paths

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
    if not (b >= j and j * a >= i * b and i >= 0):
        raise ValueError("the endpoint is not valid")

    if i == 0:
        return [tuple([1] * j)]

    if (j - 1) * a >= (i) * b:
        result = [u + (1,) for u in paths_in_triangle(i, j - 1, a, b)]
        result += [u + (0,) for u in paths_in_triangle(i - 1, j, a, b)]
        return result

    return [u + (0,) for u in paths_in_triangle(i - 1, j, a, b)]


def swap(p, i, m=1) -> tuple[int, ...]:
    r"""
    Perform a covering move in the `(a,b)`-Tamari lattice of slope parameter `m`.

    The letter at position `i` in `p` must be a `0`, followed by at
    least one `1`.

    INPUT:

    - ``p`` -- a Dyck path in the `(a \times b)`-rectangle

    - ``i`` -- integer between `0` and `a+b-1`

    OUTPUT: a Dyck path in the `(a \times b)`-rectangle

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import swap
        sage: swap((1,0,1,0,0),1)
        (1, 1, 0, 0, 0)
        sage: swap((1,1,0,0,1,1,0,0,0),3)
        (1, 1, 0, 1, 1, 0, 0, 0, 0)
        sage: swap((1,0,1,0,1,0,0,0), 1, 1)
        (1, 1, 0, 0, 1, 0, 0, 0)
        sage: swap((1,0,1,0,1,0,0,0), 1, 5/3)
        (1, 1, 0, 1, 0, 0, 0, 0)


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
        if height <= 0:
            found = True
    q = list(p)
    for k in range(i, j):
        q[k] = p[k + 1]
    q[j] = 0
    return tuple(q)


def GeneralizedTamariLattice(a, b, m=1):
    r"""
    Return the `(a,b)`-Tamari lattice of parameter `m`.

    INPUT:

    - ``a`` -- positive integer

    - ``b`` -- nonnegative integer with `a \geq b`, unless `m = 0`

    - ``m`` -- a nonnegative rational number such that `a \geq b m`

    OUTPUT:

    - a finite lattice (special case of the alt `\nu`-Tamari lattices
      in [CC2023]_)

    The elements of the lattice are
    :func:`Dyck paths<sage.combinat.dyck_word.DyckWord>` in the
    `(a \times b)`-rectangle.

    The parameter `m` (slope) is used only to define the covering relations.
    When the slope `m` is `0`, two paths are comparable if and only if
    one is always above the other.

    The usual :wikipedia:`Tamari lattice<Tamari_lattice>` of index `b`
    is the special case `a=b+1` and `m=1`.

    Other special cases give the `m`-Tamari lattices studied in [BMFPR2011]_,
    or the rational Tamari lattices when a and b are coprime and m = a/b (see [PRV2017]_).

    EXAMPLES::

        sage: from sage.combinat.tamari_lattices import GeneralizedTamariLattice
        sage: GeneralizedTamariLattice(3,2)
        Finite lattice containing 2 elements
        sage: GeneralizedTamariLattice(4,3)
        Finite lattice containing 5 elements
        sage: GeneralizedTamariLattice(7,5,2)
        Traceback (most recent call last):
        ...
        ValueError: the condition a>=b*m does not hold
        sage: P = GeneralizedTamariLattice(5,3); P
        Finite lattice containing 7 elements
        sage: P = GeneralizedTamariLattice(5, 3, m=5/3); P
        Finite lattice containing 7 elements

        sage: GeneralizedTamariLattice(3, 2, -1)
        Traceback (most recent call last):
        ...
        ValueError: m must be a nonnegative rational number
        sage: GeneralizedTamariLattice(0, 0)
        Traceback (most recent call last):
        ...
        ValueError: a must be a positive integer
        sage: GeneralizedTamariLattice(3, -1)
        Traceback (most recent call last):
        ...
        ValueError: b must be a nonnegative integer
        sage: GeneralizedTamariLattice(3/2, 1)
        Traceback (most recent call last):
        ...
        ValueError: a must be a positive integer
        sage: GeneralizedTamariLattice(2, 4, 1/4)
        Traceback (most recent call last):
        ...
        ValueError: the condition a>=b does not hold
        sage: GeneralizedTamariLattice(1, 3, 0)
        Finite lattice containing 1 elements

    TESTS::

        sage: P.coxeter_transformation()**18 == 1                                       # needs sage.libs.flint
        True

        sage: GeneralizedTamariLattice(3, 2, sqrt(2))                                   # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: m must be a nonnegative rational number
        sage: GeneralizedTamariLattice(3, QQ(2))
        Finite lattice containing 2 elements
        sage: GeneralizedTamariLattice(3, RDF(2))                                      # needs sage.rings.real_double
        Finite lattice containing 2 elements
        sage: GeneralizedTamariLattice(Zmod(2)(1), 1)                                  # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: a must be a positive integer
        sage: GeneralizedTamariLattice(7, 7, RDF(2)/5).is_lattice()                    # needs sage.rings.real_double
        True
        sage: from fractions import Fraction
        sage: GeneralizedTamariLattice(3, 2, Fraction('1/2'))
        Finite lattice containing 2 elements

    REFERENCES:

    - [BMFPR2011]_

    - [PRV2017]_

    - [CC2023]_
    """
    a = _coerce_positive_integer(a, "a")
    b = _coerce_nonnegative_integer(b, "b")
    m = _coerce_nonnegative_rational(m, "m")

    if m != 0 and a < b:
        raise ValueError("the condition a>=b does not hold")
    if a < b * m:
        raise ValueError("the condition a>=b*m does not hold")

    def covers(p):
        return [swap(p, i, m) for i in range(len(p) - 1)
                if not p[i] and p[i + 1]]

    # TODO check the exact best categories to use
    if m == 0:  # generalized Dyck lattices
        cat = FiniteLatticePosets().Distributive()
    elif m == 1 and a == b + 1:  # Tamari lattices
        cat = FiniteLatticePosets().Trim().CongruenceUniform()
    else:
        cat = FiniteLatticePosets()

    return LatticePoset({p: covers(p)
                         for p in paths_in_triangle(a, b, a, b)},
                        check=False, category=cat)


def TamariLattice(n, m=1):
    r"""
    Return the `n`-th Tamari lattice.

    Using the slope parameter `m`, one can also get the `m`-Tamari lattices.

    INPUT:

    - ``n`` -- nonnegative integer (the index)

    - ``m`` -- nonnegative integer (the slope, default: 1)

    OUTPUT: a finite lattice

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
        sage: posets.TamariLattice(-1)
        Traceback (most recent call last):
        ...
        ValueError: n must be a nonnegative integer
        sage: posets.TamariLattice(3, -1)
        Traceback (most recent call last):
        ...
        ValueError: m must be a nonnegative integer
        sage: posets.TamariLattice(3, 1/2)
        Traceback (most recent call last):
        ...
        ValueError: m must be a nonnegative integer
        sage: posets.TamariLattice(QQ(3))
        Finite lattice containing 5 elements
        sage: posets.TamariLattice(3, RDF(2))                                          # needs sage.rings.real_double
        Finite lattice containing 12 elements
        sage: posets.TamariLattice(Zmod(5)(3))                                         # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: n must be a nonnegative integer
        sage: posets.TamariLattice(3, 0)
        Finite lattice containing 1 elements

    REFERENCES:

    - [BMFPR2011]_
    """
    n = _coerce_nonnegative_integer(n, "n")
    m = _coerce_nonnegative_integer(m, "m")

    return GeneralizedTamariLattice(m * n + 1, n, m)


# a variation : the Dexter meet-semilattices


def swap_dexter(p, i) -> list[tuple[int, ...]]:
    r"""
    Perform covering moves in the `(a,b)`-Dexter posets.

    The letter at position `i` in `p` must be a `0`, followed by at
    least one `1`.

    INPUT:

    - ``p`` -- a Dyck path in the `(a \times b)`-rectangle

    - ``i`` -- integer between `0` and `a+b-1`

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

    - ``n`` -- nonnegative integer (the index)

    OUTPUT: a finite meet-semilattice

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
        sage: posets.DexterSemilattice(-1)
        Traceback (most recent call last):
        ...
        ValueError: n must be a nonnegative integer
        sage: posets.DexterSemilattice(1/2)
        Traceback (most recent call last):
        ...
        ValueError: n must be a nonnegative integer
        sage: posets.DexterSemilattice(RDF(3))                                         # needs sage.rings.real_double
        Finite meet-semilattice containing 5 elements
        sage: posets.DexterSemilattice(Zmod(5)(3))                                     # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: n must be a nonnegative integer

    REFERENCES:

    - [Cha18]_
    """
    n = _coerce_nonnegative_integer(n, "n")

    a = n + 1
    b = n

    def covers_dexter(p):
        data = [swap_dexter(p, i) for i in range(len(p) - 1)
                if not p[i] and p[i + 1]]
        return [cov for L in data for cov in L]
    return MeetSemilattice({p: covers_dexter(p)
                            for p in paths_in_triangle(a, b, a, b)},
                           check=False)
