r"""
Strata of differentials on Riemann surfaces

The space of Abelian (or quadratic) differentials is stratified by the
degrees of the zeroes (and simple poles for quadratic
differentials). Each stratum has one, two or three connected
components and each is associated to an (extended) Rauzy class. The
:meth:`~sage.dynamics.flat_surfaces.strata.AbelianStratum.connected_components`
method (only available for Abelian stratum) give the decomposition of
a stratum (which corresponds to the SAGE object
:class:`~sage.dynamics.flat_surfaces.strata.AbelianStratum`).

The work for Abelian differentials was done by Maxim Kontsevich and Anton
Zorich in [KonZor03]_ and for quadratic differentials by Erwan Lanneau in
[Lan08]_. Zorich gave an algorithm to pass from a connected component of a
stratum to the associated Rauzy class (for both interval exchange
transformations and linear involutions) in [Zor08]_ and is implemented for
Abelian stratum at different level (approximately one for each component):

- for connected stratum :meth:`~ConnectedComponentOfAbelianStratum.representative`
- for hyperellitic component :meth:`~HypConnectedComponentOfAbelianStratum.representative`
- for non hyperelliptic component, the algorithm is the same as for connected
  component
- for odd component :meth:`~OddConnectedComponentOfAbelianStratum.representative`
- for even component :meth:`~EvenConnectedComponentOfAbelianStratum.representative`

The inverse operation (pass from an interval exchange transformation to
the connected component) is partially written in [KonZor03]_ and
simply named here
:meth:`~sage.dynamics.interval_exchanges.template.PermutationIET.connected_component`.

All the code here was first available on Mathematica [ZS]_.

REFERENCES:

.. [KonZor03] M. Kontsevich, A. Zorich "Connected components of the moduli space
   of Abelian differentials with prescripebd singularities" Invent. math. 153,
   631-678 (2003)

.. [Lan08] E. Lanneau "Connected components of the strata of the moduli spaces
   of quadratic differentials", Annales sci. de l'ENS, serie 4, fascicule 1,
   41, 1-56 (2008)

.. [Zor08] A. Zorich "Explicit Jenkins-Strebel representatives of all strata of
   Abelian and quadratic differentials", Journal of Modern Dynamics, vol. 2,
   no 1, 139-185 (2008) (http://www.math.psu.edu/jmd)

.. [ZS] Anton Zorich, "Generalized Permutation software"
   (http://perso.univ-rennes1.fr/anton.zorich/Software/software_en.html)

.. NOTE::

    The quadratic strata are not yet implemented.

AUTHORS:

- Vincent Delecroix (2009-09-29): initial version


EXAMPLES:

Construction of a stratum from a list of singularity degrees::

    sage: a = AbelianStratum(1,1)
    sage: print a
    H(1, 1)
    sage: print a.genus()
    2
    sage: print a.nintervals()
    5

::

    sage: a = AbelianStratum(4,3,2,1)
    sage: print a
    H(4, 3, 2, 1)
    sage: print a.genus()
    6
    sage: print a.nintervals()
    15

By convention, the degrees are always written in decreasing order::

    sage: a1 = AbelianStratum(4,3,2,1)
    sage: a1
    H(4, 3, 2, 1)
    sage: a2 = AbelianStratum(2,3,1,4)
    sage: a2
    H(4, 3, 2, 1)
    sage: a1 == a2
    True

It is also possible to consider stratum with an incoming or an
outgoing separatrix marked (the aim of this consideration is to
attach a specified degree at the left or the right of the associated
interval exchange transformation)::

    sage: a_out = AbelianStratum(1, 1, marked_separatrix='out')
    sage: a_out
    H^out(1, 1)
    sage: a_in = AbelianStratum(1, 1, marked_separatrix='in')
    sage: a_in
    H^in(1, 1)
    sage: a_out == a_in
    False

Get a list of strata with constraints on genus or on the number of intervals
of a representative::

    sage: for a in AbelianStrata(genus=3):
    ....:     print a
    H(4)
    H(3, 1)
    H(2, 2)
    H(2, 1, 1)
    H(1, 1, 1, 1)

::

    sage: for a in AbelianStrata(nintervals=5):
    ....:     print a
    H^out(0, 2)
    H^out(2, 0)
    H^out(1, 1)
    H^out(0, 0, 0, 0)

::

    sage: for a in AbelianStrata(genus=2, nintervals=5):
    ....:     print a
    H^out(0, 2)
    H^out(2, 0)
    H^out(1, 1)

Obtains the connected components of a stratum::

    sage: a = AbelianStratum(0)
    sage: print a.connected_components()
    [H_hyp(0)]

::

    sage: a = AbelianStratum(6)
    sage: cc = a.connected_components()
    sage: print cc
    [H_hyp(6), H_odd(6), H_even(6)]
    sage: for c in cc:
    ....:     print c, "\n", c.representative(alphabet=range(1,9))
    H_hyp(6)
    1 2 3 4 5 6 7 8
    8 7 6 5 4 3 2 1
    H_odd(6)
    1 2 3 4 5 6 7 8
    4 3 6 5 8 7 2 1
    H_even(6)
    1 2 3 4 5 6 7 8
    6 5 4 3 8 7 2 1

::

    sage: a = AbelianStratum(1, 1, 1, 1)
    sage: print a.connected_components()
    [H_c(1, 1, 1, 1)]
    sage: c = a.connected_components()[0]
    sage: print c.representative(alphabet="abcdefghi")
    a b c d e f g h i
    e d c f i h g b a

The zero attached on the left of the associated Abelian permutation
corresponds to the first singularity degree::

    sage: a = AbelianStratum(4, 2, marked_separatrix='out')
    sage: b = AbelianStratum(2, 4, marked_separatrix='out')
    sage: print a == b
    False
    sage: print a, ":", a.connected_components()
    H^out(4, 2) : [H_odd^out(4, 2), H_even^out(4, 2)]
    sage: print b, ":", b.connected_components()
    H^out(2, 4) : [H_odd^out(2, 4), H_even^out(2, 4)]
    sage: a_odd, a_even = a.connected_components()
    sage: b_odd, b_even = b.connected_components()

The representatives are hence different::

    sage: print a_odd.representative(alphabet=range(1,10))
    1 2 3 4 5 6 7 8 9
    4 3 6 5 7 9 8 2 1
    sage: print b_odd.representative(alphabet=range(1,10))
    1 2 3 4 5 6 7 8 9
    4 3 5 7 6 9 8 2 1

::

    sage: print a_even.representative(alphabet=range(1,10))
    1 2 3 4 5 6 7 8 9
    6 5 4 3 7 9 8 2 1
    sage: print b_even.representative(alphabet=range(1,10))
    1 2 3 4 5 6 7 8 9
    7 6 5 4 3 9 8 2 1

You can retrieve the decomposition of the irreducible Abelian permutations into
Rauzy diagrams from the classification of strata::

    sage: a = AbelianStrata(nintervals=4)
    sage: l = sum([stratum.connected_components() for stratum in a], [])
    sage: n = map(lambda x: x.rauzy_diagram().cardinality(), l)
    sage: for c,i in zip(l,n):
    ....:     print c, ":", i
    H_hyp^out(2) : 7
    H_hyp^out(0, 0, 0) : 6
    sage: print sum(n)
    13

::

    sage: a = AbelianStrata(nintervals=5)
    sage: l = sum([stratum.connected_components() for stratum in a], [])
    sage: n = map(lambda x: x.rauzy_diagram().cardinality(), l)
    sage: for c,i in zip(l,n):
    ....:     print c, ":", i
    H_hyp^out(0, 2) : 11
    H_hyp^out(2, 0) : 35
    H_hyp^out(1, 1) : 15
    H_hyp^out(0, 0, 0, 0) : 10
    sage: print sum(n)
    71

::

    sage: a = AbelianStrata(nintervals=6)
    sage: l = sum([stratum.connected_components() for stratum in a], [])
    sage: n = map(lambda x: x.rauzy_diagram().cardinality(), l)
    sage: for c,i in zip(l,n):
    ....:     print c, ":", i
    H_hyp^out(4) : 31
    H_odd^out(4) : 134
    H_hyp^out(0, 2, 0) : 66
    H_hyp^out(2, 0, 0) : 105
    H_hyp^out(0, 1, 1) : 20
    H_hyp^out(1, 1, 0) : 90
    H_hyp^out(0, 0, 0, 0, 0) : 15
    sage: print sum(n)
    461
"""
#*****************************************************************************
#       Copyright (C) 2009 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject

from sage.combinat.combinat import CombinatorialClass
from sage.combinat.combinat import InfiniteAbstractCombinatorialClass
from sage.combinat.partition import Partitions

from sage.rings.integer import Integer
from sage.rings.rational import Rational


def AbelianStrata(genus=None, nintervals=None, marked_separatrix=None):
    r"""
    Abelian strata.

    INPUT:

    - ``genus`` - a non negative integer or ``None``

    - ``nintervals`` - a non negative integer or ``None``

    - ``marked_separatrix`` - 'no' (for no marking), 'in' (for marking an
      incoming separatrix) or 'out' (for marking an outgoing separatrix)

    EXAMPLES:

    Abelian strata with a given genus::

        sage: for s in AbelianStrata(genus=1): print s
        H(0)

    ::

        sage: for s in AbelianStrata(genus=2): print s
        H(2)
        H(1, 1)

    ::

        sage: for s in AbelianStrata(genus=3): print s
        H(4)
        H(3, 1)
        H(2, 2)
        H(2, 1, 1)
        H(1, 1, 1, 1)

    ::

        sage: for s in AbelianStrata(genus=4): print s
        H(6)
        H(5, 1)
        H(4, 2)
        H(4, 1, 1)
        H(3, 3)
        H(3, 2, 1)
        H(3, 1, 1, 1)
        H(2, 2, 2)
        H(2, 2, 1, 1)
        H(2, 1, 1, 1, 1)
        H(1, 1, 1, 1, 1, 1)

    Abelian strata with a given number of intervals::

        sage: for s in AbelianStrata(nintervals=2): print s
        H^out(0)

    ::

        sage: for s in AbelianStrata(nintervals=3): print s
        H^out(0, 0)

    ::

        sage: for s in AbelianStrata(nintervals=4): print s
        H^out(2)
        H^out(0, 0, 0)

    ::

        sage: for s in AbelianStrata(nintervals=5): print s
        H^out(0, 2)
        H^out(2, 0)
        H^out(1, 1)
        H^out(0, 0, 0, 0)

    Abelian strata with both constraints::

        sage: for s in AbelianStrata(genus=2, nintervals=4): print s
        H^out(2)

    ::

        sage: for s in AbelianStrata(genus=5, nintervals=12): print s
        H^out(8, 0, 0)
        H^out(0, 8, 0)
        H^out(0, 7, 1)
        H^out(1, 7, 0)
        H^out(7, 1, 0)
        H^out(0, 6, 2)
        H^out(2, 6, 0)
        H^out(6, 2, 0)
        H^out(1, 6, 1)
        H^out(6, 1, 1)
        H^out(0, 5, 3)
        H^out(3, 5, 0)
        H^out(5, 3, 0)
        H^out(1, 5, 2)
        H^out(2, 5, 1)
        H^out(5, 2, 1)
        H^out(0, 4, 4)
        H^out(4, 4, 0)
        H^out(1, 4, 3)
        H^out(3, 4, 1)
        H^out(4, 3, 1)
        H^out(2, 4, 2)
        H^out(4, 2, 2)
        H^out(2, 3, 3)
        H^out(3, 3, 2)
    """
    if genus is None:
        if nintervals is None:
            return AbelianStrata_all()
        else:
            return AbelianStrata_d(
                nintervals=nintervals,
                marked_separatrix=marked_separatrix)
    else:
        if nintervals is None:
            return AbelianStrata_g(
                genus=genus,
                marked_separatrix=marked_separatrix)
        else:
            return AbelianStrata_gd(
                genus=genus,
                nintervals=nintervals,
                marked_separatrix=marked_separatrix)


class AbelianStrata_g(CombinatorialClass):
    r"""
    Stratas of genus g surfaces.

    INPUT:

    - ``genus`` - a non negative integer

    - ``marked_separatrix`` - 'no', 'out' or 'in'
    """
    def __init__(self, genus=None, marked_separatrix=None):
        r"""
        TESTS::

            sage: s = AbelianStrata(genus=3)
            sage: s == loads(dumps(s))
            True

            sage: AbelianStrata(genus=-3)
            Traceback (most recent call last):
            ...
            ValueError: genus must be positive

            sage: AbelianStrata(genus=3, marked_separatrix='yes')
            Traceback (most recent call last):
            ...
            ValueError: marked_separatrix must be no, out or in
        """
        genus = Integer(genus)

        if not(genus >= 0):
            raise ValueError('genus must be positive')

        if marked_separatrix is None:
            marked_separatrix = 'no'
        if not marked_separatrix in ['no', 'out', 'in']:
            raise ValueError("marked_separatrix must be no, out or in")
        self._marked_separatrix = marked_separatrix

        self._genus = genus

    def _repr_(self):
        r"""
        TESTS::

            sage: repr(AbelianStrata(genus=3))   #indirect doctest
            'Abelian strata of genus 3 surfaces'
        """
        if self._marked_separatrix == 'no':
            return "Abelian strata of genus %d surfaces" % (self._genus)
        elif self._marked_separatrix == 'in':
            return "Abelian strata of genus %d surfaces and a marked incoming separatrix" % (self._genus)
        else:
            return "Abelian strata of genus %d surfaces and a marked outgoing separatrix" % (self._genus)

    def __iter__(self):
        r"""
        TESTS::

            sage: list(AbelianStrata(genus=1))
            [H(0)]
        """
        if self._genus == 0:
            pass
        elif self._genus == 1:
            yield AbelianStratum(0, marked_separatrix=self._marked_separatrix)
        else:
            if self._marked_separatrix == 'no':
                for p in Partitions(2*self._genus-2):
                    yield AbelianStratum(p)
            else:
                for p in Partitions(2*self._genus-2):
                    l = list(p)
                    for t in set(l):
                        i = l.index(t)
                        yield AbelianStratum([t] + l[:i] + l[i+1:],
                                             marked_separatrix=self._marked_separatrix)


class AbelianStrata_d(CombinatorialClass):
    r"""
    Strata with constraint number of intervals.

    INPUT:

    - ``nintervals`` - an integer greater than 1

    - ``marked_separatrix`` - 'no', 'out' or 'in'
    """
    def __init__(self, nintervals=None, marked_separatrix=None):
        r"""
        TESTS::

            sage: s = AbelianStrata(nintervals=10)
            sage: s == loads(dumps(s))
            True

            sage: AbelianStrata(nintervals=1)
            Traceback (most recent call last):
            ...
            ValueError: number of intervals must be at least 2

            sage: AbelianStrata(nintervals=4, marked_separatrix='maybe')
            Traceback (most recent call last):
            ...
            ValueError: marked_separatrix must be no, out or in
        """
        nintervals = Integer(nintervals)

        if not(nintervals > 1):
            raise ValueError("number of intervals must be at least 2")

        self._nintervals = nintervals

        if marked_separatrix is None:
            marked_separatrix = 'out'
        if not marked_separatrix in ['no', 'out', 'in']:
            raise ValueError("marked_separatrix must be no, out or in")
        self._marked_separatrix = marked_separatrix

    def _repr_(self):
        r"""
        TESTS::

            sage: repr(AbelianStrata(nintervals=2,marked_separatrix='no')) #indirect doctest
            'Abelian strata with 2 intervals IET'
        """
        if self._marked_separatrix == 'no':
            return "Abelian strata with %d intervals IET" % (self._nintervals)
        elif self._marked_separatrix == 'in':
            return "Abelian strata with %d intervals IET and a marked incoming separatrix" % (self._nintervals)
        else:
            return "Abelian strata with %d intervals IET and a marked outgoing separatrix" % (self._nintervals)

    def __iter__(self):
        r"""
        TESTS::

            sage: for a in AbelianStrata(nintervals=4): print a
            H^out(2)
            H^out(0, 0, 0)
        """
        n = self._nintervals
        for s in range(1+n % 2, n, 2):
            for p in Partitions(n-1, length=s):
                l = [k-1 for k in p]
                if self._marked_separatrix == 'no':
                    yield AbelianStratum(l, marked_separatrix='no')
                else:
                    for t in set(l):
                        i = l.index(t)
                        yield AbelianStratum([t] + l[:i] + l[i+1:],
                                             marked_separatrix=self._marked_separatrix)


class AbelianStrata_gd(CombinatorialClass):
    r"""
    Abelian strata of prescribed genus and number of intervals.

    INPUT:

    - ``genus`` - integer: the genus of the surfaces

    - ``nintervals`` - integer: the number of intervals

    - ``marked_separatrix`` - 'no', 'in' or 'out'
    """
    def __init__(self, genus=None, nintervals=None, marked_separatrix=None):
        r"""
        TESTS::

            sage: s = AbelianStrata(genus=4,nintervals=10)
            sage: s == loads(dumps(s))
            True

            sage: AbelianStrata(genus=-1)
            Traceback (most recent call last):
            ...
            ValueError: genus must be positive

            sage: AbelianStrata(genus=1, nintervals=1)
            Traceback (most recent call last):
            ...
            ValueError: number of intervals must be at least 2

            sage: AbelianStrata(genus=1, marked_separatrix='so')
            Traceback (most recent call last):
            ...
            ValueError: marked_separatrix must be no, out or in
        """
        genus = Integer(genus)

        if not(genus >= 0):
            raise ValueError("genus must be positive")
        self._genus = genus

        nintervals = Integer(nintervals)
        if not(nintervals > 1):
            raise ValueError("number of intervals must be at least 2")
        self._nintervals = nintervals

        if marked_separatrix is None:
            marked_separatrix = 'out'
        if not marked_separatrix in ['no', 'out', 'in']:
            raise ValueError("marked_separatrix must be no, out or in")

        self._marked_separatrix = marked_separatrix

    def __repr__(self):
        r"""
        TESTS::

            sage: a = AbelianStrata(genus=2,nintervals=4,marked_separatrix='no')
            sage: repr(a)   #indirect doctest
            'Abelian strata of genus 2 surfaces and 4 intervals'
        """
        if self._marked_separatrix == 'no':
            return "Abelian strata of genus %d surfaces and %d intervals" % (self._genus, self._nintervals)
        elif self._marked_separatrix == 'in':
            return "Abelian strata of genus %d surfaces and %d intervals and a marked incoming ihorizontal separatrix" % (self._genus, self._nintervals)
        else:
            return "Abelian strata of genus %d surfaces and %d intervals and a marked outgoing horizontal separatrix" % (self._genus, self._nintervals)

    def __iter__(self):
        r"""
        TESTS::

            sage: list(AbelianStrata(genus=2,nintervals=4))
            [H^out(2)]
        """
        if self._genus == 0:
            pass
        elif self._genus == 1:
            if self._nintervals >= 2:
                yield AbelianStratum([0]*(self._nintervals-1),
                                     marked_separatrix='out')
        else:
            s = self._nintervals - 2*self._genus + 1
            for p in Partitions(2*self._genus - 2 + s, length=s):
                l = [k-1 for k in p]
                for t in set(l):
                    i = l.index(t)
                    yield AbelianStratum([t] + l[:i] +
                                         l[i+1:], marked_separatrix='out')


class AbelianStrata_all(InfiniteAbstractCombinatorialClass):
    r"""
    Abelian strata.
    """
    def __repr__(self):
        r"""
        TESTS::

            sage: repr(AbelianStrata())   #indirect doctest
            'Abelian strata'
        """
        return "Abelian strata"

    def _infinite_cclass_slice(self, g):
        r"""
        TESTS::

            sage: AbelianStrata()[0]
            H(0)
            sage: AbelianStrata()[1]
            H(2)
            sage: AbelianStrata()[2]
            H(1, 1)

        ::

            sage: a = AbelianStrata()
            sage: a._infinite_cclass_slice(0) == AbelianStrata(genus=0)
            True
            sage: a._infinite_cclass_slice(10) == AbelianStrata(genus=10)
            True
        """
        return AbelianStrata_g(g)


class AbelianStratum(SageObject):
    """
    Stratum of Abelian differentials.

    A stratum with a marked outgoing separatrix corresponds to Rauzy diagram
    with left induction, a stratum with marked incoming separatrix correspond
    to Rauzy diagram with right induction.
    If there is no marked separatrix, the associated Rauzy diagram is the
    extended Rauzy diagram (consideration of the
    :meth:`sage.dynamics.interval_exchanges.template.Permutation.symmetric`
    operation of Boissy-Lanneau).

    When you want to specify a marked separatrix, the degree on which it is is
    the first term of your degrees list.

    INPUT:

    - ``marked_separatrix`` - ``None`` (default) or 'in' (for incoming
      separatrix) or 'out' (for outgoing separatrix).

    EXAMPLES:

    Creation of an Abelian stratum and get its connected components::

        sage: a = AbelianStratum(2, 2)
        sage: print a
        H(2, 2)
        sage: a.connected_components()
        [H_hyp(2, 2), H_odd(2, 2)]

    Specification of marked separatrix:

    ::

        sage: a = AbelianStratum(4,2,marked_separatrix='in')
        sage: print a
        H^in(4, 2)
        sage: b = AbelianStratum(2,4,marked_separatrix='in')
        sage: print b
        H^in(2, 4)
        sage: a == b
        False

    ::

        sage: a = AbelianStratum(4,2,marked_separatrix='out')
        sage: print a
        H^out(4, 2)
        sage: b = AbelianStratum(2,4,marked_separatrix='out')
        sage: print b
        H^out(2, 4)
        sage: a == b
        False

    Get a representative of a connected component::

        sage: a = AbelianStratum(2,2)
        sage: a_hyp, a_odd = a.connected_components()
        sage: print a_hyp.representative()
        1 2 3 4 5 6 7
        7 6 5 4 3 2 1
        sage: print a_odd.representative()
        0 1 2 3 4 5 6
        3 2 4 6 5 1 0

    You can choose the alphabet::

        sage: print a_odd.representative(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        A B C D E F G
        D C E G F B A

    By default, you get a reduced permutation, but you can specify
    that you want a labelled one::

        sage: p_reduced = a_odd.representative()
        sage: p_labelled = a_odd.representative(reduced=False)
    """
    def __init__(self, *l, **d):
        """
        TESTS::

            sage: s = AbelianStratum(0)
            sage: s == loads(dumps(s))
            True
            sage: s = AbelianStratum(1,1,1,1)
            sage: s == loads(dumps(s))
            True

            sage: AbelianStratum('no','way')
            Traceback (most recent call last):
            ...
            ValueError: input must be a list of integers

            sage: AbelianStratum([1,1,1,1], marked_separatrix='full')
            Traceback (most recent call last):
            ...
            ValueError: marked_separatrix must be one of 'no', 'in', 'out'
        """
        if l == ():
            pass

        elif hasattr(l[0], "__iter__") and len(l) == 1:
            l = l[0]

        if not all(isinstance(i, (Integer, int)) for i in l):
            raise ValueError("input must be a list of integers")

        if 'marked_separatrix' in d:
            m = d['marked_separatrix']

            if m is None:
                m = 'no'

            if (m != 'no' and m != 'in' and m != 'out'):
                raise ValueError("marked_separatrix must be one of 'no', "
                                 "'in', 'out'")
            self._marked_separatrix = m

        else:  # default value
            self._marked_separatrix = 'no'

        self._zeroes = list(l)

        if not self._marked_separatrix == 'no':
            self._zeroes[1:] = sorted(self._zeroes[1:], reverse=True)
        else:
            self._zeroes.sort(reverse=True)

        self._genus = sum(l)/2 + 1

        self._genus = Integer(self._genus)

        zeroes = sorted(x for x in self._zeroes if x > 0)

        if self._genus == 1:
            self._cc = (HypCCA,)

        elif self._genus == 2:
            self._cc = (HypCCA,)

        elif self._genus == 3:
            if zeroes == [2, 2] or zeroes == [4]:
                self._cc = (HypCCA, OddCCA)
            else:
                self._cc = (CCA,)

        elif len(zeroes) == 1:
            # just one zeros [2g-2]
            self._cc = (HypCCA, OddCCA, EvenCCA)

        elif zeroes == [self._genus-1, self._genus-1]:
            # two similar zeros [g-1, g-1]
            if self._genus % 2 == 0:
                self._cc = (HypCCA, NonHypCCA)

            else:
                self._cc = (HypCCA, OddCCA, EvenCCA)

        elif len([x for x in zeroes if x % 2]) == 0:
            # even zeroes [2 l_1, 2 l_2, ..., 2 l_n]
            self._cc = (OddCCA, EvenCCA)

        else:
            self._cc = (CCA, )

    def _repr_(self):
        """
        TESTS::

            sage: repr(AbelianStratum(1,1))   #indirect doctest
            'H(1, 1)'
        """
        if self._marked_separatrix == 'no':
            return "H(" + str(self._zeroes)[1:-1] + ")"
        else:
            return ("H" +
                    '^' + self._marked_separatrix +
                    "(" + str(self._zeroes)[1:-1] + ")")

    def __str__(self):
        r"""
        TESTS::

            sage: str(AbelianStratum(1,1))
            'H(1, 1)'
        """
        if self._marked_separatrix == 'no':
            return "H(" + str(self._zeroes)[1:-1] + ")"
        else:
            return ("H" +
                    '^' + self._marked_separatrix +
                    "(" + str(self._zeroes)[1:-1] + ")")

    def __eq__(self, other):
        r"""
        TESTS:

            sage: a = AbelianStratum(1,3)
            sage: b = AbelianStratum(3,1)
            sage: c = AbelianStratum(1,3,marked_separatrix='out')
            sage: d = AbelianStratum(3,1,marked_separatrix='out')
            sage: e = AbelianStratum(1,3,marked_separatrix='in')
            sage: f = AbelianStratum(3,1,marked_separatrix='in')
            sage: a == b  # no difference for unmarked
            True
            sage: c == d  # difference for out mark
            False
            sage: e == f  # difference for in mark
            False
            sage: a == c  # difference between no mark and out mark
            False
            sage: a == e  # difference between no mark and in mark
            False
            sage: c == e  # difference between out mark adn in mark
            False

            sage: a == False
            Traceback (most recent call last):
            ...
            TypeError: the right member must be a stratum
        """
        if not isinstance(self, type(other)):
            raise TypeError("the right member must be a stratum")

        return (self._marked_separatrix == other._marked_separatrix and
                self._zeroes == other._zeroes)

    def __ne__(self, other):
        r"""
        TESTS::

            sage: a = AbelianStratum(1,3)
            sage: b = AbelianStratum(3,1)
            sage: c = AbelianStratum(1,3,marked_separatrix='out')
            sage: d = AbelianStratum(3,1,marked_separatrix='out')
            sage: e = AbelianStratum(1,3,marked_separatrix='in')
            sage: f = AbelianStratum(3,1,marked_separatrix='in')
            sage: a != b  # no difference for unmarked
            False
            sage: c != d  # difference for out mark
            True
            sage: e != f  # difference for in mark
            True
            sage: a != c  # difference between no mark and out mark
            True
            sage: a != e  # difference between no mark and in mark
            True
            sage: c != e  # difference between out mark adn in mark
            True
            sage: a != False
            Traceback (most recent call last):
            ...
            TypeError: the right member must be a stratum
        """
        if not isinstance(self, type(other)):
            raise TypeError("the right member must be a stratum")

        return (self._marked_separatrix != other._marked_separatrix or
                self._zeroes != other._zeroes)

    def __cmp__(self, other):
        r"""
        The order is given by the natural:

        self < other iff adherance(self) c adherance(other)

        TESTS::

            sage: a3 = AbelianStratum(3,2,1)
            sage: a3_out = AbelianStratum(3,2,1,marked_separatrix='out')
            sage: a3_in = AbelianStratum(3,2,1,marked_separatrix='in')
            sage: a3 == a3_out
            False
            sage: a3 == a3_in
            False
            sage: a3_out == a3_in
            False
        """
        if (not isinstance(self, type(other)) or
            self._marked_separatrix != other._marked_separatrix):
            raise TypeError("the other must be a stratum with same marking")

        if self._zeroes < other._zeroes:
            return 1
        elif self._zeroes > other._zeroes:
            return -1
        return 0

    def connected_components(self):
        """
        Lists the connected components of the Stratum.

        OUTPUT:

        list -- a list of connected components of stratum

        EXAMPLES:

        ::

            sage: AbelianStratum(0).connected_components()
            [H_hyp(0)]

        ::

            sage: AbelianStratum(2).connected_components()
            [H_hyp(2)]
            sage: AbelianStratum(1,1).connected_components()
            [H_hyp(1, 1)]

        ::

            sage: AbelianStratum(4).connected_components()
            [H_hyp(4), H_odd(4)]
            sage: AbelianStratum(3,1).connected_components()
            [H_c(3, 1)]
            sage: AbelianStratum(2,2).connected_components()
            [H_hyp(2, 2), H_odd(2, 2)]
            sage: AbelianStratum(2,1,1).connected_components()
            [H_c(2, 1, 1)]
            sage: AbelianStratum(1,1,1,1).connected_components()
            [H_c(1, 1, 1, 1)]
        """
        return map(lambda x: x(self), self._cc)

    def is_connected(self):
        r"""
        Tests if the strata is connected.

        OUTPUT:

        boolean -- ``True`` if it is connected else ``False``

        EXAMPLES:

        ::

            sage: AbelianStratum(2).is_connected()
            True
            sage: AbelianStratum(2).connected_components()
            [H_hyp(2)]

        ::

            sage: AbelianStratum(2,2).is_connected()
            False
            sage: AbelianStratum(2,2).connected_components()
            [H_hyp(2, 2), H_odd(2, 2)]
        """
        return len(self._cc) == 1

    def genus(self):
        r"""
        Returns the genus of the stratum.

        OUTPUT:

        integer -- the genus

        EXAMPLES:

        ::

            sage: AbelianStratum(0).genus()
            1
            sage: AbelianStratum(1,1).genus()
            2
            sage: AbelianStratum(3,2,1).genus()
            4
        """
        return self._genus

    def nintervals(self):
        r"""
        Returns the number of intervals of any iet of the strata.

        OUTPUT:

        integer -- the number of intervals for any associated iet

        EXAMPLES:

        ::

            sage: AbelianStratum(0).nintervals()
            2
            sage: AbelianStratum(0,0).nintervals()
            3
            sage: AbelianStratum(2).nintervals()
            4
            sage: AbelianStratum(1,1).nintervals()
            5
        """
        return 2 * self.genus() + len(self._zeroes) - 1


class ConnectedComponentOfAbelianStratum(SageObject):
    r"""
    Connected component of Abelian stratum.

    .. warning::

        Internal class! Do not use directly!

    TESTS:

    Tests for outgoing marked separatrices::

        sage: a = AbelianStratum(4,2,0,marked_separatrix='out')
        sage: a_odd, a_even = a.connected_components()
        sage: a_odd.representative().attached_out_degree()
        4
        sage: a_even.representative().attached_out_degree()
        4

    ::

        sage: a = AbelianStratum(2,4,0,marked_separatrix='out')
        sage: a_odd, a_even = a.connected_components()
        sage: a_odd.representative().attached_out_degree()
        2
        sage: a_even.representative().attached_out_degree()
        2

    ::

        sage: a = AbelianStratum(0,4,2,marked_separatrix='out')
        sage: a_odd, a_even = a.connected_components()
        sage: a_odd.representative().attached_out_degree()
        0
        sage: a_even.representative().attached_out_degree()
        0

    ::

        sage: a = AbelianStratum(3,2,1,marked_separatrix='out')
        sage: a_c = a.connected_components()[0]
        sage: a_c.representative().attached_out_degree()
        3

    ::

        sage: a = AbelianStratum(2,3,1,marked_separatrix='out')
        sage: a_c = a.connected_components()[0]
        sage: a_c.representative().attached_out_degree()
        2

    ::

        sage: a = AbelianStratum(1,3,2,marked_separatrix='out')
        sage: a_c = a.connected_components()[0]
        sage: a_c.representative().attached_out_degree()
        1

    Tests for incoming separatrices::

        sage: a = AbelianStratum(4,2,0,marked_separatrix='in')
        sage: a_odd, a_even = a.connected_components()
        sage: a_odd.representative().attached_in_degree()
        4
        sage: a_even.representative().attached_in_degree()
        4

    ::

        sage: a = AbelianStratum(2,4,0,marked_separatrix='in')
        sage: a_odd, a_even = a.connected_components()
        sage: a_odd.representative().attached_in_degree()
        2
        sage: a_even.representative().attached_in_degree()
        2

    ::

        sage: a = AbelianStratum(0,4,2,marked_separatrix='in')
        sage: a_odd, a_even = a.connected_components()
        sage: a_odd.representative().attached_in_degree()
        0
        sage: a_even.representative().attached_in_degree()
        0

    ::

        sage: a = AbelianStratum(3,2,1,marked_separatrix='in')
        sage: a_c = a.connected_components()[0]
        sage: a_c.representative().attached_in_degree()
        3

    ::

        sage: a = AbelianStratum(2,3,1,marked_separatrix='in')
        sage: a_c = a.connected_components()[0]
        sage: a_c.representative().attached_in_degree()
        2

    ::

        sage: a = AbelianStratum(1,3,2,marked_separatrix='in')
        sage: a_c = a.connected_components()[0]
        sage: a_c.representative().attached_in_degree()
        1
    """
    _name = 'c'

    def __init__(self, parent):
        r"""
        TESTS::

            sage: a = AbelianStratum([1]*10).connected_components()[0]
            sage: a == loads(dumps(a))
            True
        """
        self._parent = parent

    def _repr_(self):
        r"""
        TESTS::

            sage: a = AbelianStratum([1]*8).connected_components()[0]
            sage: repr(a)   #indirect doctest
            'H_c(1, 1, 1, 1, 1, 1, 1, 1)'
        """
        if self._parent._marked_separatrix == 'no':
            return ("H" +
                    "_" + self._name +
                    "(" + str(self._parent._zeroes)[1:-1] + ")")

        else:
            return ("H" +
                    "_" + self._name +
                    "^" + self._parent._marked_separatrix +
                    "(" + str(self._parent._zeroes)[1:-1] + ")")

    def __str__(self):
        r"""
        TESTS::

            sage: str(AbelianStratum([1]*8))
            'H(1, 1, 1, 1, 1, 1, 1, 1)'
        """
        if self._parent._marked_separatrix == 'no':
            return ("H" +
                    "_" + self._name +
                    "(" + str(self._parent._zeroes)[1:-1] + ")")

        else:
            return ("H" +
                    "_" + self._name +
                    "^" + self._parent._marked_separatrix +
                    "(" + str(self._parent._zeroes)[1:-1] + ")")

    def parent(self):
        r"""
        The stratum of this component

        OUTPUT:

        stratum - the stratum where this component leaves

        EXAMPLES::

            sage: p = iet.Permutation('a b','b a')
            sage: c = p.connected_component()
            sage: c.parent()
            H(0)
        """
        return self._parent

    def representative(self, reduced=True, alphabet=None):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (default: ``True``): whether you
          obtain a reduced or labelled permutation

        - ``alphabet`` - an alphabet or ``None``: whether you want to
          specify an alphabet for your permutation

        OUTPUT:

        permutation -- a permutation which lives in this component

        EXAMPLES:

        ::

            sage: c = AbelianStratum(1,1,1,1).connected_components()[0]
            sage: print c
            H_c(1, 1, 1, 1)
            sage: p = c.representative(alphabet=range(9))
            sage: print p
            0 1 2 3 4 5 6 7 8
            4 3 2 5 8 7 6 1 0
            sage: p.connected_component()
            H_c(1, 1, 1, 1)
        """
        g = self._parent._genus
        zeroes = [x for x in self._parent._zeroes if x > 0]
        n = self._parent._zeroes.count(0)

        l0 = range(0, 4*g-3)
        l1 = [4, 3, 2]
        for k in range(5, 4*g-6, 4):
            l1 += [k, k+3, k+2, k+1]
        l1 += [1, 0]
        k = 3
        for d in zeroes:
            for i in range(d-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 2
            k += 2

        if n != 0:
            interval = range(4*g-3, 4*g-3+n)

            if self._parent._zeroes[0] == 0:
                k = l0.index(4)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)


    def genus(self):
        r"""
        Returns the genus of the surfaces in this connected component.

        OUTPUT:

        integer -- the genus of the surface

        EXAMPLES:

        ::

            sage: a = AbelianStratum(6,4,2,0,0)
            sage: c_odd, c_even = a.connected_components()
            sage: c_odd.genus()
            7
            sage: c_even.genus()
            7

        ::

            sage: a = AbelianStratum([1]*8)
            sage: c = a.connected_components()[0]
            sage: c.genus()
            5
        """
        return self._parent.genus()

    def nintervals(self):
        r"""
        Returns the number of intervals of the representative.

        OUTPUT:

        integer -- the number of intervals in any representative

        EXAMPLES:

        ::

            sage: a = AbelianStratum(6,4,2,0,0)
            sage: c_odd, c_even = a.connected_components()
            sage: c_odd.nintervals()
            18
            sage: c_even.nintervals()
            18

        ::

            sage: a = AbelianStratum([1]*8)
            sage: c = a.connected_components()[0]
            sage: c.nintervals()
            17
        """
        return self.parent().nintervals()

    def rauzy_diagram(self, reduced=True):
        r"""
        Returns the Rauzy diagram associated to this connected component.

        OUTPUT:

        rauzy diagram -- the Rauzy diagram associated to this stratum

        EXAMPLES:

        ::

            sage: c = AbelianStratum(0).connected_components()[0]
            sage: r = c.rauzy_diagram()
        """
        return self.representative(reduced=reduced).rauzy_diagram()

    def __cmp__(self, other):
        r"""
        TESTS::

            sage: a1 = AbelianStratum(1,1,1,1)
            sage: c1 = a1.connected_components()[0]
            sage: a2 = AbelianStratum(3,1)
            sage: c2 = a2.connected_components()[0]
            sage: c1 == c1
            True
            sage: c1 == c2
            False
            sage: a1 = AbelianStratum(1,1,1,1)
            sage: c1 = a1.connected_components()[0]
            sage: a2 = AbelianStratum(2, 2)
            sage: c2_hyp, c2_odd = a2.connected_components()
            sage: c1 != c1
            False
            sage: c1 != c2_hyp
            True
            sage: c2_hyp != c2_odd
            True
            sage: c1 == True
            Traceback (most recent call last):
            ...
            TypeError: other must be a connected component
        """
        if not isinstance(other, CCA):
            raise TypeError("other must be a connected component")

        if isinstance(self, type(other)):
            if self._parent._zeroes > other._parent._zeroes:
                return 1
            elif self._parent._zeroes < other._parent._zeroes:
                return -1
            return 0

        return cmp(type(self), type(other))

CCA = ConnectedComponentOfAbelianStratum


class HypConnectedComponentOfAbelianStratum(CCA):
    """
    Hyperelliptic component of Abelian stratum.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'hyp'

    def representative(self, reduced=True, alphabet=None):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (defaut: ``True``): whether you obtain
          a reduced or labelled permutation

        - ``alphabet`` - alphabet or ``None`` (defaut: ``None``):
          whether you want to specify an alphabet for your
          representative

        EXAMPLES:

        ::

            sage: c = AbelianStratum(0).connected_components()[0]
            sage: c
            H_hyp(0)
            sage: p = c.representative(alphabet="01")
            sage: p
            0 1
            1 0
            sage: p.connected_component()
            H_hyp(0)

        ::

            sage: c = AbelianStratum(0,0).connected_components()[0]
            sage: c
            H_hyp(0, 0)
            sage: p = c.representative(alphabet="abc")
            sage: p
            a b c
            c b a
            sage: p.connected_component()
            H_hyp(0, 0)

        ::

            sage: c = AbelianStratum(2).connected_components()[0]
            sage: c
            H_hyp(2)
            sage: p = c.representative(alphabet="ABCD")
            sage: p
            A B C D
            D C B A
            sage: p.connected_component()
            H_hyp(2)

        ::

            sage: c = AbelianStratum(1,1).connected_components()[0]
            sage: c
            H_hyp(1, 1)
            sage: p = c.representative(alphabet="01234")
            sage: p
            0 1 2 3 4
            4 3 2 1 0
            sage: p.connected_component()
            H_hyp(1, 1)
        """
        g = self._parent._genus
        n = self._parent._zeroes.count(0)
        m = len(self._parent._zeroes) - n

        if m == 0:  # on the torus
            if n == 1:
                l0 = [0, 1]
                l1 = [1, 0]
            elif n == 2:
                l0 = [0, 1, 2]
                l1 = [2, 1, 0]
            else:
                l0 = range(1, n+2)
                l1 = [n+1] + range(1, n+1)

        elif m == 1:  # H(2g-2,0^n) or H(0,2g-2,0^(n-1))
            l0 = range(1, 2*g+1)
            l1 = range(2*g, 0, -1)
            interval = range(2*g+1, 2*g+n+1)

            if self._parent._zeroes[0] == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        else:  # H(g-1,g-1,0^n) or H(0,g-1,g-1,0^(n-1))
            l0 = range(1, 2*g+2)
            l1 = range(2*g+1, 0, -1)
            interval = range(2*g+2, 2*g+n+2)

            if self._parent._zeroes[0] == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)

HypCCA = HypConnectedComponentOfAbelianStratum


class NonHypConnectedComponentOfAbelianStratum(CCA):
    """
    Non hyperelliptic component of Abelian stratum.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'nonhyp'

NonHypCCA = NonHypConnectedComponentOfAbelianStratum


class EvenConnectedComponentOfAbelianStratum(CCA):
    """
    Connected component of Abelian stratum with even spin structure.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'even'

    def representative(self, reduced=True, alphabet=None):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        EXAMPLES:

        ::

            sage: c = AbelianStratum(6).connected_components()[2]
            sage: c
            H_even(6)
            sage: p = c.representative(alphabet=range(8))
            sage: p
            0 1 2 3 4 5 6 7
            5 4 3 2 7 6 1 0
            sage: p.connected_component()
            H_even(6)

        ::

            sage: c = AbelianStratum(4,4).connected_components()[2]
            sage: c
            H_even(4, 4)
            sage: p = c.representative(alphabet=range(11))
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            5 4 3 2 6 8 7 10 9 1 0
            sage: p.connected_component()
            H_even(4, 4)
        """
        zeroes = [x for x in self._parent._zeroes if x > 0]
        n = self._parent._zeroes.count(0)
        g = self._parent._genus

        l0 = range(3*g-2)
        l1 = [6, 5, 4, 3, 2, 7, 9, 8]
        for k in range(10, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in zeroes:
            for i in range(d/2-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # if there are marked points we transform 0 in [3g-2, 3g-3, ...]
        if n != 0:
            interval = range(3*g-2, 3*g - 2 + n)

            if self._parent._zeroes[0] == 0:
                k = l0.index(6)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)

EvenCCA = EvenConnectedComponentOfAbelianStratum


class OddConnectedComponentOfAbelianStratum(CCA):
    r"""
    Connected component of an Abelian stratum with odd spin parity.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'odd'

    def representative(self, reduced=True, alphabet=None):
        """
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        EXAMPLES:

        ::

            sage: a = AbelianStratum(6).connected_components()[1]
            sage: print a.representative(alphabet=range(8))
            0 1 2 3 4 5 6 7
            3 2 5 4 7 6 1 0

        ::

            sage: a = AbelianStratum(4,4).connected_components()[1]
            sage: print a.representative(alphabet=range(11))
            0 1 2 3 4 5 6 7 8 9 10
            3 2 5 4 6 8 7 10 9 1 0
        """
        zeroes = [x//2 for x in self._parent._zeroes if x > 0]

        n = self._parent._zeroes.count(0)
        g = self._parent._genus

        l0 = range(3*g-2)
        l1 = [3, 2]
        for k in range(4, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in zeroes:
            for i in range(d-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # marked points
        if n != 0:
            interval = range(3*g-2, 3*g-2+n)

            if self._parent._zeroes[0] == 0:
                k = l0.index(3)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)

OddCCA = OddConnectedComponentOfAbelianStratum
