# sage_setup: distribution = sagemath-graphs
# sage.doctest: needs sage.combinat
r"""
`\nu`-Tamari lattice

A class of the `\nu`-Tamari lattice, and `(\delta,\nu)`-Tamari lattice (or alt
`\nu`-Tamari lattices), see [PRV2017]_ and [CC2023]_ for details.

These lattices depend on a parameter `\nu` where `\nu` is a path of North
and East steps. The alt `\nu`-Tamari lattice depends on an additional parameter
`\delta`, which is an increment vector with respect to `\nu`.

The elements are :func:`\nu-Dyck paths<sage.combinat.nu_dyck_word.NuDyckWord>`
which are weakly above `\nu`.

To use the provided functionality, you should import `\nu`-Tamari lattices by
typing::

    sage: from sage.combinat.nu_tamari_lattice import NuTamariLattice, AltNuTamariLattice

Then, ::

    sage: NuTamariLattice([1,1,1,0,0,1,1,0])
    Finite lattice containing 6 elements
    sage: NuTamariLattice([0,0,0,1,1,0,0,1])
    Finite lattice containing 40 elements
    sage: AltNuTamariLattice([0,0,0,1,1,0,0,1])
    Finite lattice containing 40 elements
    sage: AltNuTamariLattice([0,0,0,1,1,0,0,1], [0,1,0])
    Finite lattice containing 40 elements


The classical **Tamari lattices** and the **Generalized Tamari lattices** are
special cases of this construction and are also available with this poset::

    sage: NuTamariLattice([1,0,1,0,1,0])
    Finite lattice containing 5 elements

    sage: NuTamariLattice([1,0,0,1,0,0,1,0,0])
    Finite lattice containing 12 elements

.. SEEALSO::

    For more detailed information see :meth:`NuTamariLattice` and
    :meth:`AltNuTamariLattice`. For more
    information on the standard Tamari lattice see
    :meth:`sage.combinat.tamari_lattices.TamariLattice`,
    :meth:`sage.combinat.tamari_lattices.GeneralizedTamariLattice`

AUTHORS:

- Aram Dermenjian (2020-09-26): initial version

- Clément Chenevière (2024-02-01): added the alt `\nu`-Tamari lattices
"""
# ****************************************************************************
#    Copyright (C) 2020-2020 Aram Dermenjian <aram.dermenjian@gmail.com>
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
from sage.combinat.nu_dyck_word import NuDyckWords, NuDyckWord
from sage.combinat.posets.lattices import LatticePoset


def NuTamariLattice(nu):
    r"""
    Return the `\nu`-Tamari lattice.

    INPUT:

    - `\nu` -- list of 0s and 1s or a string of 0s and 1s

    OUTPUT: a finite lattice

    The elements of the lattice are
    :func:`\nu-Dyck paths<sage.combinat.nu_dyck_word.NuDyckWord>` weakly above
    `\nu`.

    The usual :wikipedia:`Tamari lattice<Tamari_lattice>` is the special case
    where `\nu = (NE)^h` where `h` is the height.

    Other special cases give the `m`-Tamari lattices studied in [BMFPR]_.

    EXAMPLES::

        sage: from sage.combinat.nu_tamari_lattice import NuTamariLattice
        sage: NuTamariLattice([1,0,1,0,0,1,0])
        Finite lattice containing 7 elements
        sage: NuTamariLattice([1,0,1,0,1,0])
        Finite lattice containing 5 elements
        sage: NuTamariLattice([1,0,1,0,1,0,1,0])
        Finite lattice containing 14 elements
        sage: NuTamariLattice([1,0,1,0,1,0,0,0,1])
        Finite lattice containing 24 elements
    """
    NDW = NuDyckWords(nu)
    covers = []
    elements = []
    height = NDW[0].height()
    for ndw in NDW:
        elements.append(ndw)
        for i in range(1, height + 1):
            new_ndw = ndw.mutate(i)
            if new_ndw is not None:
                covers.append([ndw, new_ndw])
    return LatticePoset([elements, covers])


def delta_swap(p, k, delta):
    r"""
    Perform a covering move in the `(\delta,\nu)`-Tamari lattice (or alt
    `\nu`-Tamari lattice, see [CC2023]_).

    The letter at position `k` is a North step of the `\nu`-Dyck word `p`, and
    must be preceded by an East step.

    The vector `\delta = (\delta_1, \dots, \delta_n)` is an increment vector
    with respect to the path `\nu`, that is to say `\delta_i \leq \nu_i`, where
    `\nu_i` is the number of East steps following the `i`-th North step of
    `\nu`.

    INPUT:

    - ``p`` -- a `\nu`-Dyck word

    - ``k`` -- integer between `0` and ``p.length()-1``

    - ``delta`` -- list of nonnegative integers of length ``p.height()``

    OUTPUT: a `\nu`-Dyck word

    EXAMPLES::

        sage: from sage.combinat.nu_tamari_lattice import delta_swap
        sage: delta_swap(NuDyckWord('0101', '0101'), 3, delta = [1, 0])
        [0, 1, 1, 0]
        sage: delta_swap(NuDyckWord('1001110100', '0100010111'), 3, [3, 1, 0, 0, 0])
        [1, 0, 1, 1, 1, 0, 0, 1, 0, 0]
        sage: delta_swap(NuDyckWord('10100101000', '01001000110'), 2, [2, 3, 0, 1])
        [1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0]
        sage: delta_swap(NuDyckWord('10100101000', '01001000110'), 2, [1, 1, 0, 0])
        [1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0]


    TESTS::

        sage: delta_swap(NuDyckWord('10011101000', '01000101110'), 0, [3, 1, 0, 0, 1])
        Traceback (most recent call last):
        ...
        ValueError: there is no such covering move
        sage: delta_swap(NuDyckWord('10011101000', '01000101110'), 1, [3, 1, 0, 0, 1])
        Traceback (most recent call last):
        ...
        ValueError: there is no such covering move
        sage: delta_swap(NuDyckWord('10011101000', '01000101110'), 11, [3, 1, 0, 0, 1])
        Traceback (most recent call last):
        ...
        ValueError: the index is greater than the length of the path
    """
    if k >= p.length():
        raise ValueError("the index is greater than the length of the path")
    # if delta is None:
    #     delta = [len(_) for _ in str(p._nu).split(sep='1')[1:]]
    if k == 0 or p[k - 1] == 1:
        raise ValueError("there is no such covering move")
    found = False
    i = p[:k].count(1)
    j = k
    alt = 0
    while not found and j <= p.length() - 1:
        if p[j]:
            alt += delta[i]
            i += 1
        else:
            alt -= 1
        if alt == 0:
            found = True
        j += 1
    q = p[:k - 1] + p[k:j] + [p[k - 1]] + p[j:]
    return NuDyckWord(q, p._nu)


def AltNuTamariLattice(nu, delta=None):
    r"""
    Return the `(\delta,\nu)`-Tamari lattice (or alt `\nu`-Tamari lattice).

    For more information, see [CC2023]_.

    The path `\nu` is a path of North steps (represented as `1` s) and East
    steps (represented as `0` s).

    The vector `\delta = (\delta_1, \dots, \delta_n)` is an increment vector
    with respect to the path `\nu`, that is to say `\delta_i \leq \nu_i`, where
    `\nu_i` is the number of `0` s following the `i`-th `1` of `\nu`. If not
    provided, `\delta` is set by default to produce the classical `\nu`-Tamari
    lattice.

    INPUT:

    - `\nu` -- list of 0s and 1s or a string of 0s and 1s

    - `\delta` -- list of nonnegative integers

    OUTPUT: a finite lattice

    EXAMPLES::

        sage: from sage.combinat.nu_tamari_lattice import AltNuTamariLattice, NuTamariLattice
        sage: AltNuTamariLattice('01001', [0, 0])
        Finite lattice containing 7 elements
        sage: AltNuTamariLattice('01001', [1, 0])
        Finite lattice containing 7 elements
        sage: AltNuTamariLattice('01001') == AltNuTamariLattice('01001', [2, 0])
        True
        sage: nu = '00100100101'; P = AltNuTamariLattice(nu); Q = NuTamariLattice(nu); P == Q
        True

    TESTS::

        sage: AltNuTamariLattice('012', [0,0])
        Traceback (most recent call last):
        ...
        ValueError: nu must be a list or a string of 0s and 1s
        sage: AltNuTamariLattice([0,10,0,11], [2,0,0])
        Traceback (most recent call last):
        ...
        ValueError: nu must be a list or a string of 0s and 1s
        sage: AltNuTamariLattice('01001', [0, 1, 0])
        Traceback (most recent call last):
        ...
        ValueError: delta is not a valid increment vector
        sage: AltNuTamariLattice('0100101', [3, 0, 0])
        Traceback (most recent call last):
        ...
        ValueError: delta is not a valid increment vector

    REFERENCES:

    - [PRV2017]_

    - [CC2023]_
    """
    if not ((isinstance(nu, (list, tuple)) and all(x in [0, 1] for x in nu)) or
            (isinstance(nu, str) and all(x in ['0', '1'] for x in nu))):
        raise ValueError("nu must be a list or a string of 0s and 1s")
    nu = [int(a) for a in nu]
    # transforms nu in a sequence of 0s and 1s if it is a list
    nu = ''.join(str(a) for a in nu)
    # produces delta if delta is None, and check that delta is valid otherwise
    deltamax = [len(a) for a in nu.split(sep='1')[1:]]
    if delta is None:
        delta = deltamax
    elif len(delta) != len(deltamax) or any(delta[i] > deltamax[i] for i in range(len(delta))):
        raise ValueError("delta is not a valid increment vector")

    def covers(p):
        return [delta_swap(p, k, delta=delta) for k in range(1, p.length())
                if not p[k - 1] and p[k]]
    return LatticePoset({p: covers(p) for p in NuDyckWords(nu)},
                        check=False)
