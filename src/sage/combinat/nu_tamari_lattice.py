# -*- coding: utf-8 -*-
r"""
`\nu`-Tamari lattice

A class of the `\nu`-Tamari lattice, see [PRV2017]_ for details.

These lattices depend on one parameter `\nu` where `\nu` is a path of North
and East steps.

The elements are :func:`\nu-Dyck paths<sage.combinat.nu_dyck_word.NuDyckWord>`
which are weakly above `\nu`.

To use the provided functionality, you should import `\nu`-Tamari lattices by
typing::

    sage: from sage.combinat.nu_tamari_lattice import NuTamariLattice

Then, ::

    sage: NuTamariLattice([1,1,1,0,0,1,1,0])
    Finite lattice containing 6 elements
    sage: NuTamariLattice([0,0,0,1,1,0,0,1])
    Finite lattice containing 40 elements

The classical **Tamari lattices** and the **Generalized Tamari lattices** are
special cases of this construction and are also available with this poset::

    sage: NuTamariLattice([1,0,1,0,1,0])
    Finite lattice containing 5 elements

    sage: NuTamariLattice([1,0,0,1,0,0,1,0,0])
    Finite lattice containing 12 elements

.. SEEALSO::

    For more detailed information see :meth:`NuTamariLattice`. For more
    information on the standard Tamari lattice see
    :meth:`sage.combinat.tamari_lattices.TamariLattice`,
    :meth:`sage.combinat.tamari_lattices.GeneralizedTamariLattice`
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
from sage.combinat.nu_dyck_word import NuDyckWords
from sage.combinat.posets.lattices import LatticePoset


def NuTamariLattice(nu):
    r"""
    Return the `\nu`-Tamari lattice.

    INPUT:

    - `\nu` -- a list of 0s and 1s or a string of 0s and 1s.

    OUTPUT:

    a finite lattice

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
