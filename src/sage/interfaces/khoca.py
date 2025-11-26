r"""
Interface to Khoca

Khoca is computer program written by Lukas Lewark to calculate sl(N)-homology
of knots and links. It calculates the following:

* Khovanov sl(2)-homology of arbitrary links, given as a braid or in PD code.
* Khovanov-Rozansky `sl(N)`-homology with `N > 2` of bipartite knots, given by a
  certain encoding of a matched diagram of the knot.
* Homology over the integers, the rationals or a prime field.
* Either equivariant homology, or homology with an arbitrary fixed potential.
* All pages of the spectral sequence of filtered homology over a field.
* Reduced and unreduced homology.
* Homology of sums and mirror images of knots.

For more details please have a look at the `Khoca repository <https://github.com/LLewark/khoca>`__.
If you are using khoca for a project or publication, please cite the web page or the literature
given there. Furthermore, note that not all functionality listed above is available through
this interface. Especially this is true for `sl(N)` homology for `N > 2` since we don't have
a method to obtain a *matched diagram* for bipartite knots.

The Khoca interface will only work if the optional Sage package Khoca is installed.

AUTHORS:

- Sebastian Oehms (2025):
"""

##############################################################################
#       Copyright (C) 2025 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.misc.cachefunc import cached_function
from enum import Enum


class KnownKeywords(Enum):
    r"""
    Enum class to specify if a keyword belongs to the interface.

    EXAMPLES::

        sage: from sage.interfaces.khoca import KnownKeywords
        sage: [kwd for kwd in KnownKeywords]
        [<KnownKeywords.frobenius_algebra: 'frobenius_algebra'>,
        <KnownKeywords.root: 'root'>,
        <KnownKeywords.equivariant: 'equivariant'>,
        <KnownKeywords.reduced: 'reduced'>,
        <KnownKeywords.code: 'code'>]
    """
    frobenius_algebra = 'frobenius_algebra'
    root = 'root'
    equivariant = 'equivariant'
    reduced = 'reduced'
    code = 'code'


@cached_function
def check_kwds(**kwds):
    r"""
    Return the keys of ``kwds`` as a list of elements of ``KnownKeywords``.
    If ``kwds`` contains a key that is not in ``KnownKeywords`` a ``KeyError``
    is raised.

    EXAMPLES::

        sage: from sage.interfaces.khoca import check_kwds
        sage: check_kwds(frobenius_algebra=(0,0), root=0)
        [<KnownKeywords.frobenius_algebra: 'frobenius_algebra'>,
        <KnownKeywords.root: 'root'>]
        sage: check_kwds(frobenius=(0,0), root=0)
        Traceback (most recent call last):
        ...
        ValueError: 'frobenius' is not a valid KnownKeywords
    """
    keylist = []
    for k in kwds:
        keylist.append(KnownKeywords(k))
    return keylist


@cached_function
def khoca_interface(ring, **kwds):
    r"""
    Return an instance of ``InteractiveCalculator` of the ``Khoca``.
    This is a calculator for Khovanov homology written by Lukas
    Lewark. For more information see :class:`~sage.features.khoca.Khoca`.

    EXAMPLES::

        sage: # optional khoca
        sage: from sage.interfaces.khoca import khoca_interface
        sage: khoca_interface(ZZ)
        Khovanov homology calculator for Frobenius algebra: Z[X] / (1*X^2).
        sage: khoca_interface(QQ)
        Khovanov homology calculator for Frobenius algebra: Q[X] / (1*X^2).
        sage: khoca_interface(GF(3))
        Khovanov homology calculator for Frobenius algebra: F_3[X] / (1*X^2).
        sage: khoca_interface(ZZ, frobenius_algebra=(1,-2), root=1)
        Khovanov homology calculator for Frobenius algebra: Z[X] / (1*X^2 + -2*X + 1).
        sage: khoca_interface(QQ, equivariant=3)
        Traceback (most recent call last):
        ...
        NotImplementedError: keyword equivariant is not implemented yet
    """
    from sage.features.khoca import Khoca
    Khoca().require()
    keys = check_kwds(**kwds)
    ch = ring.characteristic()
    rg = ch
    if rg == 0 and ring.is_field():
        rg = 1
    from khoca import InteractiveCalculator
    frobenius_algebra = (0, 0)
    if KnownKeywords.frobenius_algebra in keys:
        frobenius_algebra = kwds[KnownKeywords.frobenius_algebra.value]
    root = 0
    if KnownKeywords.root in keys:
        root = kwds[KnownKeywords.root.value]
    equivariant = None
    if KnownKeywords.equivariant in keys:
        raise NotImplementedError('keyword %s is not implemented yet' % KnownKeywords.equivariant.value)
    return InteractiveCalculator(coefficient_ring=rg,
                                 frobenius_algebra=frobenius_algebra,
                                 root=root,
                                 equivariant=equivariant)


@cached_function
def khoca_raw_data(link, ring, red_typ=True, **kwds):
    r"""
    Return the raw data for the Khovanov homology from the ``Khoca``
    calculator. This needs the optional feature ``khoca`` be present.

    INPUT:

    - ``link`` -- :class:`~sage.knots.link.Link`
    - ``ring`` -- the coefficient ring
    - ``kwds`` -- dictionary of options to be passes to ``Khoca``

    OUTPUT:

    A list of quadruples ``[degree, height, torsion, rank]`` each of which
    represents a summand of a homology group for the given degree and height.

    EXAMPLES::

        sage: # optional - khoca
        sage: from sage.interfaces.khoca import khoca_raw_data
        sage: B2 = BraidGroup(2)
        sage: b2 = B2((1,1))
        sage: L2 = Link(b2)
        sage: khoca_raw_data(L2, ZZ)
        {(0, 0, 0): 1, (2, 0, 0): 1, (4, 2, 0): 1, (6, 2, 0): 1}
        sage: khoca_raw_data(L2, ZZ, reduced=True)
        {(1, 0, 0): 1, (5, 2, 0): 1}
        sage: b3 = B2((1,1,1))
        sage: K3 = Link(b3)
        sage: khoca_raw_data(K3, ZZ, code='pd') == khoca_raw_data(K3, ZZ, code='braid')
        True
        sage: khoca_raw_data(L2, ZZ, equivariant=3)
        Traceback (most recent call last):
        ...
        NotImplementedError: keyword equivariant is not implemented yet
        sage: khoca_raw_data(K3, QQ, code='gauss')
        Traceback (most recent call last):
        ...
        ValueError: unknown code gauss, must be one of (pd, braid)
    """
    def prepare_data(data):
        r"""
        compress and adapt data to Sage
        """
        data_as_dict = {}

        for i in data:
            # i[0]: degree, i[1]: height, i[2]: torsion i[3]: rank
            d, h, t, r = i
            # make d compatible with Sage
            d = int(t / 2) - d
            if (h, d, t) in data_as_dict:
                data_as_dict[(h, d, t)] += r
            else:
                data_as_dict[(h, d, t)] = r
        return data_as_dict

    if not red_typ:
        keys = check_kwds(**kwds)
        arg = link.braid().Tietze()
        if KnownKeywords.code in keys:
            code = kwds[KnownKeywords.code.value]
            if code == 'pd':
                arg = link.pd_code()
            elif code != 'braid':
                raise ValueError('unknown code %s, must be one of (pd, braid)' % code)

        KH = khoca_interface(ring, **kwds)
        khres = KH(arg)
        return {'red': prepare_data(khres[0]), 'unred': prepare_data(khres[1])}

    raw_data = khoca_raw_data(link, ring, red_typ=False, **kwds)
    if 'reduced' in kwds:
        red = kwds['reduced']
        if isinstance(red, bool):
            if red:
                return raw_data['red']
        else:
            raise TypeError('reduced must be a boolean')
    return raw_data['unred']
