# sage_setup: distribution = sagemath-graphs
# sage.doctest: needs sage.graphs sage.groups
r"""
Access to the KnotInfo database

This module contains the class :class:`KnotInfoBase` which is derived from
:class:`Enum` and provides knots and links listed in the databases at the
web-pages `KnotInfo <https://knotinfo.math.indiana.edu/>`__
and `LinkInfo <https://linkinfo.sitehost.iu.edu/>`__ as its items.

This interface contains a set of about twenty knots and links statically as
demonstration cases. The complete database can be installed as an optional Sage
package using

- ``sage -i database_knotinfo`` (does not install if the current version is present)
- ``sage -f database_knotinfo`` (installs even if the current version is present)

This will install a `Python wrapper <https://github.com/soehms/database_knotinfo#readme>`__
for the original databases in Sage. This wrapper perfoms an automatic progress
of version numbers. For more details and further install instructions please see
the correspondig web-page.

To perform all the doctests concerning the usage of the database on the installation
add the option ``-c``. In this case (for instance ``sage -f -c database_knotinfo``)
the installation breaks on failing tests.

The installation of the complete database  will be necessary in order to have
access to all the properties recorded in the databases, as well.

If the entire database is installed as explained above, the import instructions
for :class:`KnotInfo` and :class:`KnotInfoSeries`, which can be seen in the opening
lines of the examples, are unnecessary.

Be aware that there are a couple of conventions used differently on KnotInfo as
in Sage.

For different conventions regarding normalization of the polynomial invariants see
the according documentation of :meth:`KnotInfoBase.homfly_polynomial`,
:meth:`KnotInfoBase.jones_polynomial` and :meth:`KnotInfoBase.alexander_polynomial`.

Also, note that the braid notation is used according to Sage, even thought in
the source where it is taken from, the braid generators are assumed to have a
negative crossing which would be opposite to the convention in Sage (see definition
3 of
:arxiv:`Gittings, T., "Minimum Braids: A Complete Invariant of Knots and Links" <math/0401051>`).

Furthermore, note that not all columns available in the database are visible on the web
pages. It is planned to remove non-visible columns from the database in the future (see
the `Python Wrapper <https://github.com/soehms/database_knotinfo#readme>`__ for
updated information).

EXAMPLES::

    sage: L = KnotInfo.L4a1_0
    sage: L.pd_notation()
    [[6, 1, 7, 2], [8, 3, 5, 4], [2, 5, 3, 6], [4, 7, 1, 8]]
    sage: L.pd_notation(original=True)
    '{{6, 1, 7, 2}, {8, 3, 5, 4}, {2, 5, 3, 6}, {4, 7, 1, 8}}'
    sage: L.is_knot()
    False
    sage: L.num_components()
    2

Items for knots need a leading ``K`` for technical reason::

    sage: K = KnotInfo.K4_1
    sage: K.is_knot()
    True

Injecting the variable name into the namespace::

    sage: KnotInfo.K5_1.inject()
    Defining K5_1
    sage: K5_1.dt_notation()
    [6, 8, 10, 2, 4]

Defining a link from the original name string::

    sage: KnotInfo('L6a1{1}').inject()
    Defining L6a1_1
    sage: L6a1_1.is_alternating()
    True

Obtaining an instance of :class:`~sage.groups.braid.Braid`::

    sage: L.braid()
    s1^-2*s0^-1*s1*s0^-1
    sage: type(_)
    <class 'sage.groups.braid.BraidGroup_class_with_category.element_class'>

Obtaining an instance of :class:`Link`::

    sage: l = L.link(); l
    Link with 2 components represented by 4 crossings
    sage: type(l)
    <class 'sage.knots.link.Link'>

If you have `SnapPy <https://snappy.math.uic.edu/index.html>`__ installed inside
Sage, you can obtain an instance of :class:`~spherogram.links.links_base.Link`,
too::

    sage: # optional - snappy
    sage: L6 = KnotInfo.L6a1_0
    sage: l6s = L6.link(snappy=True); l6s
    ...
    <Link: 2 comp; 6 cross>
    sage: type(l6s)
    <class 'spherogram.links.invariants.Link'>
    sage: l6  = L6.link()
    sage: l6 == l6s.sage_link()
    True
    sage: L6.link(L6.items.name, snappy=True)
    <Link L6a1: 2 comp; 6 cross>
    sage: l6sn = _
    sage: l6s == l6sn
    False
    sage: l6m = l6.mirror_image()
    sage: l6sn.sage_link().is_isotopic(l6m)
    True

But observe that the name conversion to SnapPy does not distinguish orientation
types::

    sage: L6b = KnotInfo.L6a1_1
    sage: L6b.link(L6b.items.name, snappy=True)  # optional - snappy
    <Link L6a1: 2 comp; 6 cross>
    sage: _.PD_code() == l6sn.PD_code()          # optional - snappy
    True

Obtaining the HOMFLY-PT polynomial::

    sage: L.homfly_polynomial()
    -v^-1*z - v^-3*z - v^-3*z^-1 + v^-5*z^-1
    sage: _ == l.homfly_polynomial(normalization='vz')
    True


Obtaining the original string from the database for an arbitrary property::

    sage: K[K.items.classical_conway_name]  # optional - database_knotinfo
    '4_1'

Further methods::

    sage: K.crossing_number()
    4
    sage: K.gauss_notation()
    [-1, 2, -3, 1, -4, 3, -2, 4]
    sage: K.dt_notation()
    [4, 6, 8, 2]
    sage: K.determinant()
    5
    sage: K.symmetry_type()
    'fully amphicheiral'
    sage: _ == K[K.items.symmetry_type]
    True
    sage: K.is_reversible()
    True
    sage: K.is_amphicheiral()
    True
    sage: K.jones_polynomial()                                                          # needs sage.symbolic
    t^2 - t - 1/t + 1/t^2 + 1
    sage: K.kauffman_polynomial()
    a^2*z^2 + a*z^3 - a^2 - a*z + 2*z^2 + a^-1*z^3 - 1 - a^-1*z + a^-2*z^2 - a^-2
    sage: K.alexander_polynomial()
    t^2 - 3*t + 1

Using the ``column_type`` of a property::

    sage: def select_column(i):
    ....:     return i.column_type() != i.types.OnlyLinks and K[i] == 'Y'
    sage: [i.column_name() for i in K.items if select_column(i)]  # optional - database_knotinfo
    ['Alternating', 'Fibered', 'Quasialternating', 'Adequate']

You can launch web-pages attached to the links::

    sage: # not tested
    sage: K.diagram()
    True
    sage: L.diagram(single=True)
    True
    sage: L.knot_atlas_webpage()
    True
    sage: K.knotilus_webpage()
    True

and the description web-pages of the properties::

    sage: K.items.positive.description_webpage()  # not tested
    True

To see all the properties available in this interface you can use "tab-completion".
For example type ``K.items.`` and than hit the :kbd:`Tab` key. You can select the item
you want from the list. If you know some first letters type them first to obtain a
reduced selection list.

In a similar way you may select the knots and links. Here you have to type ``KnotInfo.``
or ``KnotInfo.L7`` before stroking the :kbd:`Tab` key. In the latter case  the selection list
will be reduced to proper links with 7 crossings.

Finally there is a method :meth:`Link.get_knotinfo` of class :class:`Link` to find an instance
in the KnotInfo database::

    sage: L = Link([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
    ....:           [17,19,8,18], [9,10,11,14], [10,12,13,11],
    ....:           [12,19,15,13], [20,16,14,15], [16,20,17,2]])
    sage: L.get_knotinfo()
    KnotInfo['K0_1']


REFERENCES:

- `KnotInfo <https://knotinfo.math.indiana.edu/>`__
- `LinkInfo <https://linkinfo.sitehost.iu.edu/>`__


AUTHORS:

- Sebastian Oehms August 2020: initial version
- Sebastian Oehms June   2022: add :meth:`conway_polynomial` and :meth:`khovanov_polynomial` (:issue:`33969`)

Thanks to Chuck Livingston and Allison Moore for their support. For further acknowledgments see the correspondig hompages.
"""


##############################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################


from enum import Enum
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.sage_eval import sage_eval
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.knots.knot import Knots
from sage.databases.knotinfo_db import KnotInfoColumns, db

lazy_import('sage.groups.braid', 'BraidGroup')


def eval_knotinfo(string, locals={}, to_tuple=True):
    r"""
    Preparse a string from the KnotInfo database and evaluate it by ``sage_eval``.

    INPUT:

    - ``string`` -- string that gives a value of some database entry
    - ``locals`` -- dictionary of locals passed to ``sage_eval``

    EXAMPLES::

        sage: from sage.knots.knotinfo import eval_knotinfo
        sage: L = KnotInfo.L4a1_0
        sage: L.braid_notation(original=True)
        '{3, {-2, -2, -1, 2, -1}}'
        sage: eval_knotinfo(_)
        (3, (-2, -2, -1, 2, -1))
        sage: KnotInfo.K13a_1.kauffman_polynomial()  # optional - database_knotinfo # indirect doctest
        Traceback (most recent call last):
        ...
        NotImplementedError: this value is not provided by the database
    """
    if not string:
        # An empty string in the database Excel spreadsheet indicates that
        # the property is not provided for that particular knot or link.
        raise NotImplementedError('this value is not provided by the database')
    if to_tuple:
        new_string = string.replace('{', '(')
        new_string = new_string.replace('}', ')')
    else:
        new_string = string.replace('{', '[')
        new_string = new_string.replace('}', ']')
    new_string = new_string.replace(';', ',')
    return sage_eval(new_string, locals=locals)


def knotinfo_int(string):
    r"""
    Preparse a string from the KnotInfo database representing an integer.

    INPUT:

    - ``string`` -- string that gives a value of some database entry

    EXAMPLES::

        sage: from sage.knots.knotinfo import knotinfo_int
        sage: knotinfo_int('7')
        7
        sage: KnotInfo.K13a_1.braid_index() # optional - database_knotinfo # indirect doctest
        Traceback (most recent call last):
        ...
        NotImplementedError: this integer is not provided by the database
    """
    if not string:
        # an empty string in the Excel sheet of the database indicates that
        # the property is not provided for this special knot or link.
        raise NotImplementedError('this integer is not provided by the database')
    else:
        return int(string)


def knotinfo_bool(string):
    r"""
    Preparse a string from the KnotInfo database representing a boolean.

    INPUT:

    - ``string`` -- string that gives a value of some database entry

    EXAMPLES::

        sage: from sage.knots.knotinfo import knotinfo_bool
        sage: knotinfo_bool('Y')
        True
        sage: KnotInfo.K13a_1.is_almost_alternating() # optional - database_knotinfo # indirect doctest
        Traceback (most recent call last):
        ...
        NotImplementedError: this boolean is not provided by the database
    """
    if not string:
        # an empty string in the Excel sheet of the database indicates that
        # the property is not provided for this special knot or link.
        raise NotImplementedError('this boolean is not provided by the database')
    if string == 'Y':
        return True
    elif string == 'N':
        return False
    raise ValueError('%s is not a KnotInfo boolean')


class SymmetryMutant(Enum):
    r"""
    Enum to specify the symmetry mutant link of the prime link listed in the
    KnotInfo and LinkInfo databases. From the KnotInfo description page:

        If a knot is viewed as the oriented diffeomorphism
        class of an oriented pair, `K = (S_3, S_1)`, with `S_i`
        diffeomorphic to `S^i`, there are four oriented knots
        associated to any particular knot `K`. In addition to
        `K` itself, there is the reverse, `K^r = (S_3, -S_1)`,
        the concordance inverse, `-K = (-S_3, -S_1)`, and the
        mirror image, `K^m = (-S_3, S_1)`.
    """
    itself = 's'
    reverse = 'r'
    concordance_inverse = 'c'
    mirror_image = 'm'
    mixed = 'x'  # to be used in connection with KnotInfoSeries
    unknown = '?'

    def __gt__(self, other):
        r"""
        Implement comparison of different items in order to have ``sorted`` work.

        EXAMPLES::

            sage: from sage.knots.knotinfo import SymmetryMutant
            sage: sorted(SymmetryMutant)        # indirect doctest
            [<SymmetryMutant.mixed: 'x'>,
            <SymmetryMutant.itself: 's'>,
            <SymmetryMutant.reverse: 'r'>,
            <SymmetryMutant.mirror_image: 'm'>,
            <SymmetryMutant.concordance_inverse: 'c'>,
            <SymmetryMutant.unknown: '?'>]
        """
        # We use the reversal of the alphabetical order of the values so that
        # `itself` occurs before the mirrored cases
        return self.value < other.value

    def rev(self):
        r"""
        Return the reverse of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import SymmetryMutant
            sage: all( sym.rev().rev() == sym for sym in SymmetryMutant)
            True
        """
        if self is SymmetryMutant.itself:
            return SymmetryMutant.reverse
        elif self is SymmetryMutant.reverse:
            return SymmetryMutant.itself
        elif self is SymmetryMutant.mirror_image:
            return SymmetryMutant.concordance_inverse
        elif self is SymmetryMutant.concordance_inverse:
            return SymmetryMutant.mirror_image
        return self

    def mir(self):
        r"""
        Return the mirror image of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import SymmetryMutant
            sage: all( sym.mir().mir() == sym for sym in SymmetryMutant)
            True
        """
        if self is SymmetryMutant.itself:
            return SymmetryMutant.mirror_image
        elif self is SymmetryMutant.reverse:
            return SymmetryMutant.concordance_inverse
        elif self is SymmetryMutant.mirror_image:
            return SymmetryMutant.itself
        elif self is SymmetryMutant.concordance_inverse:
            return SymmetryMutant.reverse
        return self

    def matches(self, link):
        r"""
        Return the list of other symmetry mutants that give isotopic links
        with respect to ``link`` and ``self``. For ``self`` is
        ``SymmetryMutant.unknown`` a boolean is returned which is ``True``
        if the chirality of ``link`` is unknown.

        EXAMPLES::

            sage: from sage.knots.knotinfo import SymmetryMutant
            sage: SymmetryMutant.itself.matches(KnotInfo.K6_1)
            [<SymmetryMutant.reverse: 'r'>]
            sage: SymmetryMutant.mirror_image.matches(KnotInfo.K6_1)
            [<SymmetryMutant.concordance_inverse: 'c'>]
        """
        rev = link.is_reversible()
        achp = link.is_amphicheiral(positive=True)
        ach = link.is_amphicheiral()
        if self is SymmetryMutant.unknown:
            if rev is None or ach is None or achp is None:
                return True
            else:
                return False
        res = []
        if rev:
            res.append(self.rev())
        if achp:
            res.append(self.mir())
        if ach:
            res.append(self.rev().mir())
        return res

    def is_minimal(self, link):
        r"""
        Return whether ``self`` is minimal among its matching mutants.

        EXAMPLES::

            sage: from sage.knots.knotinfo import SymmetryMutant
            sage: SymmetryMutant.itself.is_minimal(KnotInfo.K6_1)
            True
            sage: SymmetryMutant.concordance_inverse.is_minimal(KnotInfo.K6_1)
            False
        """
        if self in [SymmetryMutant.unknown, SymmetryMutant.mixed]:
            return False
        matches = self.matches(link)
        return all(self < other for other in matches)


# ---------------------------------------------------------------------------------
# KnotInfoBase
# ---------------------------------------------------------------------------------
class KnotInfoBase(Enum):
    r"""
    Enum class to select the knots and links listed in the databases at the web-pages
    `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and `LinkInfo <https://linkinfo.sitehost.iu.edu/>`__.

    EXAMPLES::

        sage: [knot.name for knot in KnotInfo if knot.crossing_number() < 5]
        ['K0_1', 'K3_1', 'K4_1', 'L2a1_0', 'L2a1_1', 'L4a1_0', 'L4a1_1']

    More examples and information can be seen in the module header
    :mod:`~sage.knots.knotinfo` (by typing)::

        sage: import sage.knots.knotinfo   # not tested
        sage: sage.knots.knotinfo?         # not tested

    TESTS:

        sage: KnotInfo.K7_1.inject()
        Defining K7_1
        sage: TestSuite(K7_1).run()
    """

    def __gt__(self, other):
        r"""
        Implement comparison of different items in order to have ``sorted`` work.

        EXAMPLES::

            sage: KnotInfo.L4a1_0 < KnotInfo.L4a1_1 # indirect doctest
            True
            sage: KnotInfo.L2a1_0 < KnotInfo.K3_1   # indirect doctest
            False
            sage: KnotInfo.K10_3 > KnotInfo.K3_1    # optional - database_knotinfo
            True
        """
        if self.__class__ is other.__class__:
            tups = (not self.is_knot(), self.crossing_number(), self.value)
            tupo = (not other.is_knot(), other.crossing_number(), other.value)
            return tups > tupo
        return NotImplemented

    @property
    def items(self):
        r"""
        Return an Enum class to select a column item of the KnotInfo database.

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: it = L.items
            sage: [i.name for i in it if i.name.startswith('braid')]
            ['braid_index', 'braid_length', 'braid_notation', 'braid_notation_old']
            sage: L.items.dt_notation.column_name()
            'DT Notation'

        To check if the item is available for proper links or only knots type::

            sage: it.gauss_notation.column_type()
            <KnotInfoColumnTypes.KnotsAndLinks: 'B'>
            sage: it.dt_notation.column_type()
            <KnotInfoColumnTypes.OnlyKnots: 'K'>

        To see the description of the item in your web browser type::

            sage: it.gauss_notation.description_webpage()    # not tested
            True
        """
        return db.columns()

    @cached_method
    def __getitem__(self, item):
        r"""
        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L[L.items.alternating]
            'Y'
            sage: L[L.items.arc_notation]
            '{{6, 4}, {3, 5}, {4, 2}, {1, 3}, {2, 6}, {5, 1}}'
            sage: L[L.items.braid_notation]
            '{3, {-2, -2, -1, 2, -1}}'
            sage: L[0]
            Traceback (most recent call last):
            ...
            KeyError: "item must be an instance of <enum 'KnotInfoColumns'>"
        """
        if not isinstance(item, KnotInfoColumns):
            raise KeyError('item must be an instance of %s' % (KnotInfoColumns))
        if item.column_type() == item.types.OnlyLinks and self.is_knot():
            raise KeyError('item not available for knots' % (KnotInfoColumns))
        if item.column_type() == item.types.OnlyKnots and not self.is_knot():
            raise KeyError('item not available for links' % (KnotInfoColumns))

        l = db.read(item)
        ind = db.read_row_dict()[self.name][0]
        offset = 0
        if item.column_type() == item.types.OnlyLinks:
            offset = self._offset_knots()

        return l[ind - offset]

    def _offset_knots(self):
        r"""
        Return the list index of the first proper link in a combined
        list containing knots and proper links together which is the
        case for columns used for KnotInfo and LinkInfo in common.
        This index is exactly the total number of knots recorded
        in KnotInfo.

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L._offset_knots()          # optional - database_knotinfo
            12966
        """
        return db.read_num_knots()

    @cached_method
    def _braid_group(self):
        r"""
        Return the braid group corresponding to the braid index
        of ``self``.

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L._braid_group()
            Braid group on 3 strands
        """
        try:
            n = self.braid_index()
        except NotImplementedError:
            bn = self.braid_notation()
            n = max(abs(i) for i in bn) + 1

        if n == 1:
            return BraidGroup(2)
        else:
            return BraidGroup(n)

    @cached_method
    def _homfly_pol_ring(self, var1, var2):
        r"""
        Return the parent Laurent polynomial ring for the HOMFLY-PT
        polynomial according to Sage's internal one.

        EXAMPLES::

            sage: L = KnotInfo.L4a1_1
            sage: L._homfly_pol_ring('u', 'v')
            Multivariate Laurent Polynomial Ring in u, v over Integer Ring
        """
        K3_1 = Knots().from_table(3, 1)
        return K3_1.homfly_polynomial(var1=var1, var2=var2).parent()

    @cached_method
    def pd_notation(self, original=False):
        r"""
        Return the value of column ``pd_notation`` for this
        link as a Python list of Python lists. For more information
        type ``KnotInfo.K0_1.items.pd_notation.description_webpage()``.

        INPUT:

        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        Python list of python lists each entry of the outer list
        representing a crossing.

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L.pd_notation()
            [[6, 1, 7, 2], [8, 3, 5, 4], [2, 5, 3, 6], [4, 7, 1, 8]]
            sage: L.pd_notation(original=True)
            '{{6, 1, 7, 2}, {8, 3, 5, 4}, {2, 5, 3, 6}, {4, 7, 1, 8}}'
            sage: K = KnotInfo.K4_1
            sage: K.pd_notation()
            [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]
        """
        if self.is_knot():
            pd_notation = self[self.items.pd_notation]
        else:
            pd_notation = self[self.items.pd_notation_vector]

        if original:
            return pd_notation

        if not pd_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(pd_notation, to_tuple=False)

    @cached_method
    def dt_notation(self, original=False):
        r"""
        Return the value of column ``dt_notation`` for this
        link as a Python list of Python lists. For more information
        type ``KnotInfo.K0_1.items.dt_notation.description_webpage()``.

        INPUT:

        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        Python list of python lists each entry of the outer list
        representing a crossing.

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L.dt_notation()
            [[6, 8], [2, 4]]
            sage: L.dt_notation(original=True)
            '[{6, 8}, {2, 4}]'
            sage: L = KnotInfo.L4a1_0
            sage: K = KnotInfo.K4_1
            sage: K.dt_notation()
            [4, 6, 8, 2]
        """
        if self.is_knot():
            dt_notation = self[self.items.dt_notation]
        else:
            dt_notation = self[self.items.dt_code]

        if original:
            return dt_notation

        if not dt_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(dt_notation, to_tuple=False)

    @cached_method
    def gauss_notation(self, original=False):
        r"""
        Return the value of column ``gauss_notation`` for this
        link as a Python list of Python lists. For more information
        type ``KnotInfo.K0_1.items.gauss_notation.description_webpage()``.

        INPUT:

        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L.gauss_notation()
            [[1, -3, 2, -4], [3, -1, 4, -2]]
            sage: L.gauss_notation(original=True)
            '{{1, -3, 2, -4}, {3, -1, 4, -2}}'
        """
        gauss_notation = self[self.items.gauss_notation]
        if original:
            return gauss_notation

        if not gauss_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(gauss_notation, to_tuple=False)

    @cached_method
    def braid_notation(self, original=False):
        r"""
        Return the value of column ``braid_notation`` for this
        link as a Python tuple (Tietze form). For more information
        type ``KnotInfo.K0_1.items.braid_notation.description_webpage()``.

        INPUT:

        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        Python tuple representing the braid whose closure is ``self``
        in Tietze form.

        ..NOTE::

            There has been a major change to braid representatives for
            proper links since version 2021.10.1. The former braid
            reresentatives can be obtained by the column
            ``braid_notation_old`` (see the final example below).

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L.braid_notation()
            (-2, -2, -1, 2, -1)
            sage: L.braid_notation(original=True)
            '{3, {-2, -2, -1, 2, -1}}'
            sage: L[L.items.braid_notation_old]
            '{4, {1, -2, 3, -2, -1, -2, -3, -2}}'

        TESTS:

        Check that :issue:`33966` is fixed::

            sage: KnotInfo.K0_1.braid_notation()
            (1,)
        """
        braid_notation = self[self.items.braid_notation]
        if original:
            return braid_notation

        if not braid_notation:
            # don't forget the unknot
            return (1, )

        braid_notation = eval_knotinfo(braid_notation)
        if type(braid_notation) is list:
            # in some cases there are a pair of braid representations
            # in the database. If this is the case we select the
            # corresponding to the braid index.
            if type(braid_notation[0]) is tuple:
                i = self.braid_index()
                for b in braid_notation:
                    if -i < min(b) and max(b) < i:
                        braid_notation = b
                        break

        if not self.is_knot():
            # in linkinfo the braid_notation includes the braid_index as
            # first item of a pair
            braid_notation = braid_notation[1]
        return braid_notation

    @cached_method
    def braid_index(self):
        r"""
        Return the value of column ``braid_index`` for this
        link as a Python int.

        OUTPUT:

        Python int giving the minimum of strands needed to
        represent ``self`` as closure of a braid.

        EXAMPLES::

            sage: L = KnotInfo.L4a1_0
            sage: L.braid_index()
            3
            sage: KnotInfo.K13a_1.inject()    # optional - database_knotinfo
            Defining K13a_1
            sage: K13a_1.braid_index()        # optional - database_knotinfo
            Traceback (most recent call last):
            ...
            NotImplementedError: this integer is not provided by the database
        """
        if self.is_knot():
            return knotinfo_int(self[self.items.braid_index])
        else:
            braid_notation = self[self.items.braid_notation]
            braid_notation = eval_knotinfo(braid_notation)
            return knotinfo_int(braid_notation[0])

    @cached_method
    def braid_length(self):
        r"""
        Return the value of column ``braid_length`` for this
        link as a Python int.

        OUTPUT:

        Python int giving the minimum length of a braid word
        needed to represent ``self`` as closure of a braid.

        EXAMPLES::

            sage: K = KnotInfo.K3_1
            sage: K.braid_length()
            3
        """
        return knotinfo_int(self[self.items.braid_length])

    @cached_method
    def braid(self):
        r"""
        Return the braid notation of ``self`` as an instance of :class:`~sage.groups.braid.Braid`.

        EXAMPLES::

            sage: K = KnotInfo.K3_1
            sage: K.braid()
            s^3
            sage: K.braid_notation()
            (1, 1, 1)
            sage: KnotInfo.K13n_1448.braid()    # optional - database_knotinfo
            s0^-1*s1*s2*s3*s4*s3^2*s2^-1*s1^-1*s0*s2^-1*s1*(s3*s2)^2*s4^-1*s3*s2*s1^-1*s3*s2^-1*s3
        """
        return self._braid_group()(self.braid_notation())

    @cached_method
    def num_components(self):
        r"""
        Return the number of components of ``self``.

        EXAMPLES::

            sage: KnotInfo.L6a1_0.num_components()
            2
        """
        return db.read_row_dict()[self.name][1]

    @cached_method
    def crossing_number(self):
        r"""
        Return the minimal number of crossings of ``self``.

        .. NOTE::

           In contrast to the number of crossings displayed for instances
           of :class:`Link` this number is the minimum over all possible
           diagrams of the link. The number of crossings displayed in
           the representation string of :class:`Link` refers to the
           special diagram which could be larger.

        EXAMPLES::

            sage: KnotInfo.L4a1_0.crossing_number()
            4
            sage: KnotInfo.K3_1.crossing_number()
            3
            sage: Link(KnotInfo.L4a1_0.braid())
            Link with 2 components represented by 5 crossings
        """
        return knotinfo_int(self[self.items.crossing_number])

    @cached_method
    def determinant(self):
        r"""
        Return the determinant of ``self``.

        From the KnotInfo description page:

            The determinant of a knot is `\det(V + V^t)`, where `V` is a Seifert
            matrix for the knot.

        To read the complete description type
        ``KnotInfo.K0_1.items.determinant.description_webpage()``.

        .. NOTE::

           KnotInfo's value for the unknot ``0_1`` is zero. This is not
           compatible whith Sage's result (the value of the Alexander
           polynomial at -1). Since this method is needed to identify
           Sage links we take the according value in that case.

        EXAMPLES::

            sage: KnotInfo.L4a1_0.determinant()
            4
            sage: KnotInfo.K3_1.determinant()
            3
            sage: KnotInfo.K0_1.determinant()
            1
        """
        if self.crossing_number() == 0:
            # see note above
            return 1
        return knotinfo_int(self[self.items.determinant])

    @cached_method
    def three_genus(self):
        r"""
        Return the three genus of ``self``.

        From the KnotInfo description page:

            The three-genus of a knot is defined to be the minimal genus of
            a Seifert surface for a knot.

        To read the complete description type
        ``KnotInfo.K0_1.items.three_genus.description_webpage()``.

        EXAMPLES::

            sage: KnotInfo.K5_2.three_genus()     # optional - database_knotinfo
            1

        Note that this differs from the corresponding result in Sage
        since the latter is obtained for a Seifert surface that does not
        have the minimal genus::

            sage: KnotInfo.K5_2.link().genus()
            3
        """
        return knotinfo_int(self[self.items.three_genus])

    @cached_method
    def signature(self):
        r"""
        Return the signature of ``self``.

        From the KnotInfo description page:

            The signature of a knot, `\sigma (K)`, is equal to `\sigma (V + V^t)`,
            the signature of `V + V^t` where `V` is a Seifert matrix for the knot
            and `V^t` is its transpose.

        To read the complete description type
        ``KnotInfo.K0_1.items.signatur.description_webpage()``.

        EXAMPLES::

            sage: KnotInfo.K5_2.signature()       # optional - database_knotinfo
            -2
        """
        return knotinfo_int(self[self.items.signature])

    @cached_method
    def is_knot(self):
        r"""
        Return whether ``self`` is a knot or a proper link.

        EXAMPLES::

            sage: KnotInfo.L7a1_0.is_knot()      # optional - database_knotinfo
            False
            sage: KnotInfo.K6_3.is_knot()
            True
        """
        return self.num_components() == 1

    @cached_method
    def name_unoriented(self):
        r"""
        Return the part of the name of ``self`` which is independent on the
        orientation.

        EXAMPLES::

            sage: KnotInfo.L10a122_1_0.name_unoriented()  # optional - database_knotinfo
            'L10a122'
        """
        return self[self.items.name_unoriented]

    @cached_method
    def symmetry_type(self):
        r"""
        Return the symmetry type of ``self``.

        From the KnotInfo description page:

            If a knot is viewed as the oriented diffeomorphism
            class of an oriented pair, `K = (S_3, S_1)`, with `S_i`
            diffeomorphic to `S^i`, there are four oriented knots
            associated to any particular knot `K`. In addition to
            `K` itself, there is the reverse, `K^r = (S_3, -S_1)`,
            the concordance inverse, `-K = (-S_3, -S_1)`, and the
            mirror image, `K^m = (-S_3, S_1)`. A knot is called
            reversible if `K = K^r`, negative amphicheiral if
            `K = -K`, and positive amphicheiral if `K = K^m`.

            A knot possessing any two of these types of symmetry
            has all three. Thus, in the table, a knot is called
            reversible if that is the only type of symmetry it has,
            and likewise for negative amphicheiral. If it has none
            of these types of symmetry it is called chiral, and if
            it has all three it is called fully amphicheiral.

            For prime knots with fewer than 12 crossings, all
            amphicheiral knots are negative amphicheiral.

        EXAMPLES::

            sage: KnotInfo.K6_1.series().inject()
            Defining K6
            sage: [(K.name, K.symmetry_type()) for K in K6]
            [('K6_1', 'reversible'),
             ('K6_2', 'reversible'),
             ('K6_3', 'fully amphicheiral')]
        """
        if not self.is_knot():
            raise NotImplementedError('this is only available for knots')

        symmetry_type = self[self.items.symmetry_type].strip()  # for example K10_88 is a case with trailing whitespaces
        if not symmetry_type and self.crossing_number() == 0:
            return 'fully amphicheiral'
        return symmetry_type

    @cached_method
    def is_reversible(self):
        r"""
        Return whether ``self`` is reversible.

        EXAMPLES::

            sage: KnotInfo.K6_3.is_reversible()
            True

        TESTS::

            sage: KnotInfo.K10_67.is_reversible() # optional - database_knotinfo
            False
            sage: KnotInfo.L7a4_0.is_reversible() # optional - database_knotinfo
        """
        if self.is_knot():
            symmetry_type = self.symmetry_type()
            if symmetry_type == 'reversible':
                return True
            if symmetry_type == 'fully amphicheiral':
                return True
            return False

        # revert orientation
        b = self.braid()
        bt = list(b.Tietze())
        bt.reverse()
        br = b.parent()(tuple(bt))
        if b.is_conjugated(br):
            return True
        return None

    @cached_method
    def is_amphicheiral(self, positive=False):
        r"""
        Return whether ``self`` is amphicheiral.

        INPUT:

        - ``positive`` -- boolean (default: ``False``); whether to check
          if ``self`` is positive or negative amphicheiral (see documentation
          of :meth:`symmetry_type`)

        OUTPUT: boolean or ``None`` if this cannot be determined

        ``True`` if ``self`` is fully or negative amphicheiral per default. If
        ``positive`` is set to ``True`` than fully and positive amphicheiral
        links give ``True``.

        .. NOTE::

            For proper links this property is not provided in the database.
            Anyway, we support it here in this case, as well, except for a few
            items where it cannot be determined easily and where ``None``
            is returned as answer.

        EXAMPLES::

            sage: # optional - database_knotinfo
            sage: Kp = KnotInfo.K12a_427
            sage: Kp.is_amphicheiral()
            False
            sage: Kp.is_amphicheiral(positive=True)
            True
            sage: Kn = KnotInfo.K10_88
            sage: Kn.is_amphicheiral()
            True
            sage: Kn.is_amphicheiral(positive=True)
            False
            sage: KnotInfo.L4a1_0.is_amphicheiral()
            False
            sage: KnotInfo.L10n59_1.is_amphicheiral()
            True
            sage: KnotInfo.L10n36_0.inject()
            Defining L10n36_0
            sage: L10n36_0.is_amphicheiral() is None
            True
        """
        if self.is_knot():
            symmetry_type = self.symmetry_type()
            if positive:
                if symmetry_type == 'positive amphicheiral':
                    return True
            else:
                if symmetry_type == 'negative amphicheiral':
                    return True

            if symmetry_type == 'fully amphicheiral':
                return True
            return False

        h = self.homfly_polynomial()
        v, z = h.parent().gens()
        hm = h.subs(v=~v, z=-z)
        if h != hm:
            return False

        k = self.kauffman_polynomial()
        a, z = k.parent().gens()
        km = k.subs(a=~a)
        if k != km:
            return False

        b = self.braid()
        bi = ~b
        if b.is_conjugated(bi):
            # at least negative amphicheiral
            if not positive:
                return True

        # revert orientation (back)
        bit = list(bi.Tietze())
        bit.reverse()
        bm = b.parent()(tuple(bit))
        if b.is_conjugated(bm):
            if positive:
                return True

        return None

    @cached_method
    def is_alternating(self):
        r"""
        Return whether ``self`` is alternating.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_alternating()
            True
        """
        return knotinfo_bool(self[self.items.alternating])

    @cached_method
    def is_almost_alternating(self):
        r"""
        Return whether ``self`` is almost alternating.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_almost_alternating() # optional - database_knotinfo
            False
        """
        db._feature.require()    # column not available in demo-version
        return knotinfo_bool(self[self.items.almost_alternating])

    @cached_method
    def is_quasi_alternating(self):
        r"""
        Return whether ``self`` is quasi alternating.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_quasi_alternating() # optional - database_knotinfo
            True
        """
        db._feature.require()    # column not available in demo-version
        return knotinfo_bool(self[self.items.quasi_alternating])

    @cached_method
    def is_adequate(self):
        r"""
        Return whether ``self`` is adequate.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_adequate()         # optional - database_knotinfo
            True
        """
        db._feature.require()    # column not available in demo-version
        return knotinfo_bool(self[self.items.adequate])

    @cached_method
    def is_positive(self):
        r"""
        Return whether ``self`` is positive.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_positive()
            True
        """
        return knotinfo_bool(self[self.items.positive])

    @cached_method
    def is_quasipositive(self):
        r"""
        Return whether ``self`` is quasi-positive.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_quasipositive()     # optional - database_knotinfo
            True
        """
        db._feature.require()    # column not available in demo-version
        return knotinfo_bool(self[self.items.quasipositive])

    @cached_method
    def is_strongly_quasipositive(self):
        r"""
        Return whether ``self`` is strongly quasi-positive.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_strongly_quasipositive() # optional - database_knotinfo
            True
        """
        db._feature.require()    # column not available in demo-version
        return knotinfo_bool(self[self.items.strongly_quasipositive])

    @cached_method
    def is_positive_braid(self):
        r"""
        Return whether ``self`` is a positive braid.

        EXAMPLES::

            sage: KnotInfo.K5_2.is_positive_braid()         # optional - database_knotinfo
            False
        """
        db._feature.require()    # column not available in demo-version
        return knotinfo_bool(self[self.items.positive_braid])

    @cached_method
    def is_fibered(self):
        r"""
        Return whether ``self`` is fibered.

        EXAMPLES::

            sage: KnotInfo.K6_3.is_fibered()
            True
        """
        return knotinfo_bool(self[self.items.fibered])

    @cached_method
    def is_oriented(self):
        r"""
        Return whether ``self`` is oriented.

        EXAMPLES::

            sage: KnotInfo.L6a2_1.is_oriented()
            True
        """
        return not knotinfo_bool(self[self.items.unoriented])

    @cached_method
    def homfly_polynomial(self, var1='v', var2='z', original=False):
        r"""
        Return the HOMFLY-PT polynomial according to the value of column
        ``homfly_polynomial`` for this knot or link (in the latter case the
        column ``homflypt_polynomial`` is used) as an instance of the
        element class according to the output of :meth:`Link.homfly_polynomial`
        of :class:`Link`.

        The HOMFLY-PT polynomial `P(L)` of a link `L` satisfies the following skein
        relation (see the corresponding `KnotInfo description page
        <https://knotinfo.math.indiana.edu/descriptions/jones_homfly_kauffman_description/polynomial_defn.html)>`__):

        .. MATH::

            P(O) = 1,\,\,\,   v^{-1} P(L_+) -  v P(L_-) = z P(L_0)

        INPUT:

        - ``var1`` -- string (default: ``'v'``); for the name of the first variable
        - ``var2`` -- string (default: ``'z'``); for the name of the second variable
        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        A Laurent polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.
        If ``original`` is set to ``True`` then a string is returned.

        .. NOTE::

            The skein-relation for the HOMFLY-PT polynomial given on KnotInfo
            does not match the default used in Sage. For comparison you have
            to use the keyword argument ``normalization='vz'`` on the side
            of Sage.

        EXAMPLES::

            sage: K3_1 = KnotInfo.K3_1
            sage: PK3_1 = K3_1.homfly_polynomial(); PK3_1
            -v^4 + v^2*z^2 + 2*v^2
            sage: K3_1.homfly_polynomial(original=True)
            '(2*v^2-v^4)+v^2*z^2'
            sage: PK3_1 == K3_1.link().homfly_polynomial(normalization='vz')
            True

        for proper links::

            sage: L4a1_1 = KnotInfo.L4a1_1
            sage: PL4a1_1 = L4a1_1.homfly_polynomial(var1='x', var2='y'); PL4a1_1
            -x^5*y + x^3*y^3 - x^5*y^-1 + 3*x^3*y + x^3*y^-1
            sage: _ == L4a1_1.link().homfly_polynomial('x', 'y', 'vz')
            True

        check the skein-relation from the KnotInfo description page (applied to one
        of the positive crossings of the right-handed trefoil)::

            sage: R = PK3_1.parent()
            sage: PO = R.one()
            sage: L2a1_1 = KnotInfo.L2a1_1
            sage: PL2a1_1 = L2a1_1.homfly_polynomial()
            sage: v, z = R.gens()
            sage: ~v*PK3_1 -v*PO == z*PL2a1_1
            True

        TESTS::

            sage: H = KnotInfo.L11n459_1_1_1.homfly_polynomial()   # optional - database_knotinfo
            sage: all(L.homfly_polynomial() == L.link().homfly_polynomial(normalization='vz')\
                      for L in KnotInfo if L.crossing_number() < 7)
            True

        REFERENCES:

        - :wikipedia:`HOMFLY_polynomial`
        """
        if self.is_knot():
            homfly_polynomial = self[self.items.homfly_polynomial]
        else:
            homfly_polynomial = self[self.items.homflypt_polynomial]

        if original:
            return homfly_polynomial

        R = self._homfly_pol_ring(var1, var2)
        if not homfly_polynomial and self.crossing_number() == 0:
            return R.one()

        # As of February 2021 there is a wrong character for the link in the
        # last row of the database. This is removed here (see SPKG.rst and
        # the test above). Once this is fixed upstream, the following three
        # lines of code can be removed again:
        if self.value == 'L11n459{1,1,1}':
            if homfly_polynomial.endswith('}'):
                homfly_polynomial = homfly_polynomial.strip('}')

        L, M = R.gens()
        lc = {'v': L, 'z': M}
        return eval_knotinfo(homfly_polynomial, locals=lc)

    @cached_method
    def kauffman_polynomial(self, var1='a', var2='z', original=False):
        r"""
        Return the Kauffman polynomial according to the value of column
        ``kauffman_polynomial`` for this knot or link as an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.

        The Kauffman polynomial `F(L)` respectivlely its corresponding invariant
        under regular isotopy `\Delta (L) = a^{w(L)} F(L)` where `w(L)` is the
        writhe of the link `L` satisfies the following skein relation
        (see the corresponding `KnotInfo description page
        <https://knotinfo.math.indiana.edu/descriptions/jones_homfly_kauffman_description/polynomial_defn.html)>`__):

        .. MATH::

            \Delta(O) = 1,\,\,\,   \Delta(L_+) -  \Delta(L_-) = z (\Delta(L_0 + \Delta(L_{\infty}))

        Furthermore, removing a curl of sign `\epsilon` leads to a multiplication
        of `\Delta(L)` with `a^{\epsilon}`.

        INPUT:

        - ``var1`` -- (default: ``'a'``) the first variable
        - ``var2`` -- (default: ``'z'``) the second variable
        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        A Laurent polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.
        If ``original`` is set to ``False`` then a string is returned.

        EXAMPLES::

            sage: L = KnotInfo.L2a1_1
            sage: K = KnotInfo.K4_1

            sage: L.kauffman_polynomial()
            a^-1*z - a^-1*z^-1 + a^-2 + a^-3*z - a^-3*z^-1
            sage: K.kauffman_polynomial()
            a^2*z^2 + a*z^3 - a^2 - a*z + 2*z^2 + a^-1*z^3 - 1 - a^-1*z + a^-2*z^2 - a^-2

        Comparison with Jones polynomial::

            sage: # needs sage.symbolic
            sage: k    = _
            sage: a, z = k.variables()
            sage: j    = K.jones_polynomial(skein_normalization=True)
            sage: t,   = j.variables()
            sage: k.subs(a=-t^3, z=~t+t) == j.subs(t=t^4)
            True

        Check the skein relation::

            sage: K3_1    = KnotInfo.K3_1
            sage: FK3_1   = K3_1.kauffman_polynomial()
            sage: FL2a1_1 = L.kauffman_polynomial()
            sage: z, a    = FK3_1.variables()
            sage: ΔK3_1   = FK3_1   * a**K3_1.link().writhe()
            sage: ΔL2a1_1 = FL2a1_1 * a**L.link().writhe()
            sage: ΔO1p    = a          # unknot with one positive curl
            sage: ΔO2n    = a**-2      # unknot with two negative curls
            sage: ΔK3_1 + ΔO1p == z*(ΔL2a1_1 + ΔO2n)
            True

        REFERENCES:

        - :wikipedia:`Kauffman_polynomial`
        """
        kauffman_polynomial = self[self.items.kauffman_polynomial]

        if original:
            return kauffman_polynomial

        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        R = LaurentPolynomialRing(ZZ, (var1, var2))
        if not kauffman_polynomial and self.crossing_number() == 0:
            return R.one()

        a, z = R.gens()
        lc = {'a': a, 'z': z}
        return R(eval_knotinfo(kauffman_polynomial, locals=lc))

    @cached_method
    def jones_polynomial(self, variab=None, skein_normalization=False, puiseux=False, original=False, use_sqrt=False):
        r"""
        Return the Jones polynomial according to the value of column
        ``jones_polynomial`` for this knot or link as an element of the symbolic
        ring :class:`~sage.symbolic.ring.SR` or an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial`
        depending on the keyword ``skein_normalization``. Using the keyword
        ``puiseux`` instead of an element of the symbolic ring an instance of
        :class:`~sage.rings.puiseux_series_ring_element.PuiseuxSeries` can be
        returned.

        The Jones polynomial `V(L)` of a link `L` satisfies the following skein
        relation (see the corresponding `KnotInfo description page
        <https://knotinfo.math.indiana.edu/descriptions/jones_homfly_kauffman_description/polynomial_defn.html)>`__):

        .. MATH::

            V(O) = 1,\,\,\,   t^{-1} V(L_+) -  t V(L_-) = (t^{\frac{1}{2}} - t^{-\frac{1}{2}}) V(L_0)

        INPUT:

        - ``variab`` -- variable (default: ``None``); used according to
          :meth:`Link.jones_polynomial`
        - ``skein_normalization`` -- boolean (default: ``False``); used
          according to :meth:`Link.jones_polynomial`
        - ``puiseux`` -- boolean (default: ``True``); only used in case
          ``skein_normalization=False``. If set to ``True`` instead of an element
          of the symbolic ring an instance of :class:`~sage.rings.puiseux_series_ring_element.PuiseuxSeries`
          is returned
        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string
        - ``use_sqrt`` -- boolean (default: ``False``); see the note below

        OUTPUT:

        Depends on the keywords (in excluding order):

        - ``original=True`` a string according to the original value from the
          database
        - ``skein_normalization=True`` a Laurent polynomial over the integers,
          more precisely an instance of :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial`
        - ``puiseux=True`` a puiseux series over the integers, more precisely an
          instance of :class:`~sage.rings.puiseux_series_ring_element.PuiseuxSeries`

        In all other cases an element of the symbolic ring :class:`~sage.symbolic.ring.SR`.

        .. NOTE::

            There is a difference to Sage's conventions concerning the Jones
            polynomial in the case of proper links. KnotInfo does not display
            these polynomials in the indeterminate `t` used in the skein relation.
            Instead a variable `x` is used defined by `x^2 = t`. Sage uses `t` in
            both cases, knots and proper links. Thus, to obtain the Jones polynomial
            for a proper link in `t` you have to set the keyword ``use_sqrt``
            to ``True``.

        EXAMPLES::

            sage: K = KnotInfo.K4_1
            sage: Kj = K.jones_polynomial(); Kj                                         # needs sage.symbolic
            t^2 - t - 1/t + 1/t^2 + 1
            sage: Kjs = K.jones_polynomial(skein_normalization=True); Kjs
            A^-8 - A^-4 + 1 - A^4 + A^8
            sage: Kjp = K.jones_polynomial(puiseux=True); Kjp
            t^-2 - t^-1 + 1 - t + t^2

        for proper links::

            sage: L = KnotInfo.L2a1_1
            sage: Lj = L.jones_polynomial(); Lj                                         # needs sage.symbolic
            -x^5 - x
            sage: Ljt = L.jones_polynomial(use_sqrt=True); Ljt                          # needs sage.symbolic
            -t^(5/2) - sqrt(t)
            sage: Ljp = L.jones_polynomial(puiseux=True); Ljp
            -t^(1/2) - t^(5/2)
            sage: Ljs = L.jones_polynomial(skein_normalization=True); Ljs
            -A^2 - A^10
            sage: Lj.parent()                                                           # needs sage.symbolic
            Symbolic Ring
            sage: Ljt.parent()                                                          # needs sage.symbolic
            Symbolic Ring
            sage: Ljp.parent()
            Puiseux Series Ring in t over Integer Ring
            sage: Ljs.parent()
            Univariate Laurent Polynomial Ring in A over Integer Ring

        Comparison with Sage's results::

            sage: k = K.link()
            sage: kj = k.jones_polynomial()                                             # needs sage.symbolic
            sage: bool(Kj == kj)                                                        # needs sage.symbolic
            True
            sage: kjs = k.jones_polynomial(skein_normalization=True)
            sage: Kjs == kjs
            True
            sage: l = L.link()
            sage: lj = l.jones_polynomial()                                             # needs sage.symbolic
            sage: bool(Lj == lj)                                                        # needs sage.symbolic
            False
            sage: bool(Ljt == lj)   # see note above                                    # needs sage.symbolic
            True
            sage: ljs = l.jones_polynomial(skein_normalization=True)
            sage: Ljs == ljs
            True

        Check the skein-relation from the KnotInfo description page (applied to one
        of the positive crossings of the right-handed trefoil)::

            sage: # needs sage.symbolic
            sage: K3_1  = KnotInfo.K3_1

            sage: # needs sage.symbolic
            sage: K3_1j = K3_1.jones_polynomial()
            sage: L2a1_1j = Ljt     # see note above
            sage: R = L2a1_1j.parent()
            sage: Oj = R(1)
            sage: t = R('t')
            sage: lhs = expand(~t*K3_1j - t*Oj)
            sage: rhs = expand((sqrt(t) - ~sqrt(t))*L2a1_1j)
            sage: bool(lhs == rhs)
            True

        The same with the Puiseux series version::

            sage: K3_1jp = K3_1.jones_polynomial(puiseux=True)
            sage: L2a1_1jp = Ljp
            sage: R = L2a1_1jp.parent()
            sage: Ojp = R(1)
            sage: t = R('t')
            sage: ~t*K3_1jp - t*Ojp == (t^(1/2)-~t^(1/2))*L2a1_1jp
            True

        The same in the case of skein normalization (using `t = A^4`)::

            sage: K3_1js = K3_1.jones_polynomial(skein_normalization=True)
            sage: L2a1_1js = L.jones_polynomial(skein_normalization=True)
            sage: Rs = K3_1js.parent()
            sage: Ojs = Rs.one()
            sage: A, = Rs.gens()
            sage: ~A^4*K3_1js - A^4*Ojs == (A^2-~A^2)*L2a1_1js
            True

        REFERENCES:

        - :wikipedia:`Jones_polynomial`
        """
        jones_polynomial = self[self.items.jones_polynomial]

        if original:
            return jones_polynomial

        if skein_normalization:
            if not variab:
                variab = 'A'
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
            R = LaurentPolynomialRing(ZZ, variab)
        else:
            if not variab:
                if use_sqrt or self.is_knot() or puiseux:
                    variab = 't'
                else:
                    variab = 'x'
            if puiseux:
                from sage.rings.puiseux_series_ring import PuiseuxSeriesRing  # since PuiseuxPolynomial is not available, so far
                R = PuiseuxSeriesRing(ZZ, variab)
            else:
                from sage.symbolic.ring import SR
                R = SR

        if not jones_polynomial and self.crossing_number() == 0:
            return R.one()

        t = R(variab)
        if skein_normalization:
            if self.is_knot():
                lc = {'t': t**4}
            else:
                lc = {'x': t**2}
        else:
            if self.is_knot():
                lc = {'t': t}
            elif puiseux:
                lc = {'x': t**(1/2)}
            elif use_sqrt:
                from sage.misc.functional import sqrt
                lc = {'x': sqrt(t)}
            else:
                lc = {'x': t}

        return R(eval_knotinfo(jones_polynomial, locals=lc))

    @cached_method
    def alexander_polynomial(self, var='t', original=False, laurent_poly=False):
        r"""
        Return the Alexander polynomial according to the value of column
        ``alexander_polynomial`` for this knot as an instance of
        :class:`~sage.rings.polynomial.polynomial_element.Polynomial`.

        It is obtained from the Seifert matrix `V` of ``self`` by the following
        formula (see the KnotInfo description web-page; to launch it see the
        example below):

        .. MATH::

            A(L) = \det(V -t V^t)

        Here `V^t` stands for the transpose of `V`.


        INPUT:

        - ``var`` -- (default: ``'t'``) the variable
        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string
        - ``laurent_poly`` -- boolean (default: ``False``); see the note below

        OUTPUT:

        A polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.polynomial_element.Polynomial`.
        If ``laurent_poly`` is set to ``True`` a Laurent polynomial
        over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial`
        is returned. If ``original`` is set to ``True`` then a string
        is returned.

        .. NOTE::

            As an invariant the Alexander polynomial is only unique up to
            a unit factor in the Laurent polynomial ring over the integers
            in the indeterminate `t`. While the normalization of the exponents
            in KnotInfo guarantees it to be a proper polynomial, this is
            not the case for the implementation in Sage. Use the keyword
            ``laurent_poly`` to achiev a normalization according to Sage.
            But, still there may be a difference in sign (see the example below).

        EXAMPLES::

            sage: K = KnotInfo.K4_1
            sage: Ka = K.alexander_polynomial(); Ka
            t^2 - 3*t + 1

        Comparison with Sage's results::

            sage: k = K.link()
            sage: ka = k.alexander_polynomial(); ka
            -t^-1 + 3 - t
            sage: K.alexander_polynomial(laurent_poly=True)
            t^-1 - 3 + t
            sage: _ == -ka
            True

        Launch the KnotInfo description web-page::

            sage: K.items.alexander_polynomial.description_webpage() # not tested
            True
        """
        alexander_polynomial = self[self.items.alexander_polynomial]

        if original:
            return alexander_polynomial

        if laurent_poly:
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
            R = LaurentPolynomialRing(ZZ, var)
        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(ZZ, var)

        if not alexander_polynomial and self.crossing_number() == 0:
            return R.one()

        t, = R.gens()
        lc = {'t': t}
        ap = R(eval_knotinfo(alexander_polynomial, locals=lc))
        if not laurent_poly or ap.is_constant():
            return ap

        exp = ap.exponents()
        return t ** ((-max(exp) - min(exp)) // 2) * ap

    @cached_method
    def conway_polynomial(self, var='t', original=False):
        r"""
        Return the Conway polynomial according to the value of column
        ``conway_polynomial`` for this knot or link as an instance of
        :class:`~sage.rings.polynomial.polynomial_element.Polynomial`.

        It is obtained from the Seifert matrix `V` of ``self`` by the following
        formula (see the KnotInfo description web-page; to launch it see the
        example below):

        .. MATH::

            \nabla(L) = \det(t^{\frac{1}{2}} V -t^{\frac{-1}{2}} V^t)

        Here `V^t` stands for the transpose of `V`.


        INPUT:

        - ``var`` -- (default: ``'t'``) the variable
        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        A polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.polynomial_element.Polynomial`.
        If ``original`` is set to ``True`` then a string is returned.

        EXAMPLES::

            sage: K = KnotInfo.K4_1
            sage: Kc = K.conway_polynomial(); Kc
            -t^2 + 1
            sage: L = KnotInfo.L5a1_0
            sage: Lc = L.conway_polynomial(); Lc
            t^3

        Comparison to Sage's results::

            sage: Kc == K.link().conway_polynomial()
            True
            sage: Lc == L.link().conway_polynomial()
            True

        Launch the KnotInfo description web-page::

            sage: K.items.conway_polynomial.description_webpage()  # not tested
            True
        """
        conway_polynomial = self[self.items.conway_polynomial]

        if original:
            return conway_polynomial

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(ZZ, var)

        if not conway_polynomial and self.crossing_number() == 0:
            return R.one()

        t, = R.gens()
        lc = {'z': t}
        return R(eval_knotinfo(conway_polynomial, locals=lc))

    @cached_method
    def khovanov_polynomial(self, var1='q', var2='t', base_ring=ZZ, original=False, reduced=False, odd=False, KhoHo=False):
        r"""
        Return the Khovanov polynomial according to the value of column
        ``khovanov_polynomial`` for this knot or link as an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.

        INPUT:

        - ``var1`` -- (default: ``'q'``) the first variable
        - ``var2`` -- (default: ``'t'``) the second variable
        - ``base_ring`` -- (default: ``ZZ``) the ring of the polynomial's
          coefficients
        - ``original`` -- boolean (default: ``False``); if set to
          ``True`` the original table entry is returned as a string
        - ``reduced`` -- boolean (default: ``False``); if set to ``True``
          the reduced version of the homology is used
        - ``odd`` -- boolean (default: ``False``); if set to ``True``
          the odd version of the homology is used
        - ``KhoHo`` -- boolean (deprecated). The corresponding values have
          disappeared from the database since January 2024

        OUTPUT:

        A Laurent polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.
        If ``original`` is set to ``True`` then a string is returned.

        .. NOTE ::

            The data used for multi-component links were calculated with the program
            `KhoHo <https://github.com/AShumakovitch/KhoHo>`__.which uses the ``DT``
            notation. For knots data calculated with
            `KnotJob <https://www.maths.dur.ac.uk/users/dirk.schuetz/knotjob.html>`__
            are used. The latter program is more accurate in terms of orientation
            and reflection as it is based on ``PD`` code.

        EXAMPLES::

            sage: K = KnotInfo.K6_3
            sage: Kk = K.khovanov_polynomial(); Kk
            q^7*t^3 + q^5*t^2 + q^3*t^2 + q^3*t + q*t + 2*q + 2*q^-1 + q^-1*t^-1
            + q^-3*t^-1 + q^-3*t^-2 + q^-5*t^-2 + q^-7*t^-3
            sage: Kk2 = K.khovanov_polynomial(var1='p', base_ring=GF(2)); Kk2
            p^7*t^3 + p^5*t^3 + p^5*t^2 + p^3*t + p^-1 + p^-1*t^-1 + p^-3*t^-2 + p^-7*t^-3

            sage: L = KnotInfo.L5a1_0
            sage: Lk = L.khovanov_polynomial(); Lk
            q^4*t^2 + t + 2 + 2*q^-2 + q^-2*t^-1 + q^-4*t^-2 + q^-6*t^-2 + q^-8*t^-3
            sage: L.khovanov_polynomial(original=True)
             '2 + 2/q^2 + 1/(q^8*t^3) + 1/(q^6*t^2) + 1/(q^4*t^2) + 1/(q^2*t) + t + q^4*t^2'

        Obtaining the reduced homology (for knots only)::

            sage: Kkr = K.khovanov_polynomial(reduced=True); Kkr
            q^6*t^3 + 2*q^4*t^2 + 2*q^2*t + 3 + 2*q^-2*t^-1 + 2*q^-4*t^-2 + q^-6*t^-3
            sage: K.khovanov_polynomial(base_ring=QQ, reduced=True) == Kkr
            True
            sage: Kkr2 = K.khovanov_polynomial(base_ring=GF(2), reduced=True); Kkr2
            q^6*t^3 + 1 + q^-6*t^-3
            sage: KnotInfo.K8_19.inject()                               # optional database_knotinfo
            Defining K8_19
            sage: K8kr = K8_19.khovanov_polynomial(reduced=True); K8kr  # optional database_knotinfo
            q^16*t^5 + q^12*t^4 + q^12*t^3 + q^10*t^2 + q^6

        Obtaining the odd Khovanov homology (for knots only)::

            sage: K.khovanov_polynomial(odd=True) == Kkr
            True
            sage: K.khovanov_polynomial(base_ring=QQ, odd=True) == Kkr
            True
            sage: K.khovanov_polynomial(base_ring=GF(2), odd=True) == Kkr2
            True
            sage: K8ko = K8_19.khovanov_polynomial(odd=True); K8ko     # optional database_knotinfo
            q^16*t^5 + q^10*t^2 + q^6
            sage: K8kr == K8ko                                         # optional database_knotinfo
            False


        Comparison to Sage's results::

            sage: Kk == K.link().khovanov_polynomial()
            True
            sage: Kk2 == K.link().khovanov_polynomial(var1='p', base_ring=GF(2))
            True
            sage: Lk == L.link().khovanov_polynomial()
            True

        TESTS::

            sage: KnotInfo.K0_1.inject()
            Defining K0_1
            sage: K0_1.khovanov_polynomial()
            q + q^-1
            sage: K0_1.khovanov_polynomial(reduced=True)
            1
            sage: K0_1.khovanov_polynomial(odd=True)
            1
            sage: K0_1.khovanov_polynomial(base_ring=GF(3), reduced=True)
            Traceback (most recent call last):
            ...
            ValueError: characteristic 3 of base ring is not valid
            sage: K0_1.khovanov_polynomial(base_ring=GF(3), odd=True)
            Traceback (most recent call last):
            ...
            ValueError: characteristic 3 of base ring is not valid
            sage: L.khovanov_polynomial(base_ring=GF(2))
            Traceback (most recent call last):
            ...
            NotImplementedError: Khovanov polynomial available only for knots in characteristic 2

        REFERENCES:

        - :wikipedia:`Khovanov_homology`
        - :wikipedia:`Reduced_homology`
        - [ORS2013]_
        """
        ch = base_ring.characteristic()
        integral = ch == 0 and base_ring.is_field()
        if not self.is_knot():
            # KnotJob calculated results only available for knots
            khovanov_polynomial = self[self.items.khovanov_polynomial]
            KhoHo = True
        else:
            if KhoHo:
                KhoHo = False
                from sage.misc.superseded import deprecation
                deprecation(37014, "the KhoHo option is deprecated and ignored.")
            if reduced:
                if integral:
                    khovanov_polynomial = self[self.items.khovanov_reduced_integral_polynomial]
                elif ch == 0:
                    khovanov_polynomial = self[self.items.khovanov_reduced_rational_polynomial]
                elif ch == 2:
                    khovanov_polynomial = self[self.items.khovanov_reduced_mod2_polynomial]
                else:
                    raise ValueError('characteristic %s of base ring is not valid' % ch)
            elif odd:
                if integral:
                    khovanov_polynomial = self[self.items.khovanov_odd_integral_polynomial]
                elif ch == 0:
                    khovanov_polynomial = self[self.items.khovanov_odd_rational_polynomial]
                elif ch == 2:
                    khovanov_polynomial = self[self.items.khovanov_odd_mod2_polynomial]
                else:
                    raise ValueError('characteristic %s of base ring is not valid' % ch)
            else:
                khovanov_polynomial = self[self.items.khovanov_unreduced_integral_polynomial]

        if original:
            return khovanov_polynomial

        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        var_names = [var1, var2]
        R = LaurentPolynomialRing(base_ring, var_names)

        if not khovanov_polynomial and self.crossing_number() == 0:
            if reduced or odd:
                return R.one()
            else:
                return R({(1, 0): 1, (-1, 0): 1})

        if ch == 2:
            if not self.is_knot():
                raise NotImplementedError('Khovanov polynomial available only for knots in characteristic 2')
            if KhoHo:
                khovanov_torsion_polynomial = self[self.items.khovanov_torsion_polynomial]
                khovanov_torsion_polynomial = khovanov_torsion_polynomial.replace('Q', 'q')
                khovanov_polynomial = '%s + %s' % (khovanov_polynomial, khovanov_torsion_polynomial)

        if not khovanov_polynomial:
            # given just for links with less than 12 crossings
            raise NotImplementedError('Khovanov polynomial not available for this link')

        from sage.repl.preparse import implicit_mul
        # since implicit_mul does not know about the choice of variable names
        # we have to insert * between them separately
        for i in ['q', 't', 'T', ')']:
            for j in ['q', 't', 'T', '(']:
                khovanov_polynomial = khovanov_polynomial.replace('%s%s' % (i, j), '%s*%s' % (i, j))
        khovanov_polynomial = implicit_mul(khovanov_polynomial)
        gens = R.gens_dict()
        lc = {}
        lc['q'] = gens[var1]
        lc['t'] = gens[var2]
        if ch == 2:
            lc['T'] = 1
        else:
            lc['T'] = 0

        return R(eval_knotinfo(khovanov_polynomial, locals=lc))

    @cached_method
    def link(self, use_item=db.columns().pd_notation, snappy=False):
        r"""
        Return ``self`` as an instance of :class:`Link` or optional
        ``spherogram.links.invariants.Link``  (SnapPy).

        INPUT:

        - ``use_item`` -- (default: ``self.items.pd_notation``)
          instance of :class:`KnotInfoColumns` to choose the column
          that should be used to construct the link. Allowed values
          are:

          - ``self.items.pd_notation``
          - ``self.items.braid_notation``
          - ``self.items.name``           (only for ``snappy=True``)
          - ``self.items.dt_notation``    (only for knots and ``snappy=False``)
          - ``self.items.gauss_notation`` (only for knots and ``snappy=False``)

        - ``snappy`` -- boolean (default: ``False``); if set to ``True``
          the target of the conversion is the ``pip`` installable
          package `SnapPy <https://snappy.math.uic.edu/index.html>`__
          (explicitely, ``spherogram.links.invariants.Link``).
          If SnapPy is not installed an :exc:`ImportError` is raised. To
          install SnapPy use ``sage -pip install snappy``.

        .. NOTE::

            We use the PD-notation to construct ``self`` as
            default. This ensures that the number of crossings
            displayed in the representation string of the link
            coincides with the crossing number as a topological
            invariant.

            Furthermore, note that the mirror version may depend
            on the used KnotInfo-notation. For instance, regarding to
            the knot ``5_1`` the Gauss- and the DT-notation refer to
            the mirror image (see example below).

        EXAMPLES::

            sage: K = KnotInfo.K3_1
            sage: K.link()
            Knot represented by 3 crossings
            sage: _.braid()
            s^3
            sage: _ == K.braid()
            True

        using ``dt_notation``::

            sage: K.link(use_item=K.items.dt_notation)
            Knot represented by 3 crossings
            sage: _.braid()
            s^3

            sage: L = KnotInfo.L4a1_0
            sage: L.link()
            Link with 2 components represented by 4 crossings

            sage: L.link(use_item=L.items.dt_notation)
            Traceback (most recent call last):
            ...
            ValueError: link construction using Columns.dt_notation not possible

        using ``snappy``::

            sage: L2  = KnotInfo.L2a1_1
            sage: l2  = L2.link()
            sage: l2s = L2.link(snappy=True).sage_link()  # optional -  snappy
            sage: l2 == l2s                               # optional -  snappy
            True

        but observe::

            sage: K7   = KnotInfo.K7_2
            sage: k7s  = K7.link(snappy=True); k7s        # optional - snappy
            <Link: 1 comp; 7 cross>
            sage: k7sn = K7.link(K7.items.name, snappy=True); k7sn     # optional - snappy
            <Link 7_2: 1 comp; 7 cross>
            sage: k7s.sage_link().is_isotopic(k7sn)       # optional - snappy
            False
            sage: k7snm = k7sn.sage_link().mirror_image() # optional - snappy
            sage: k7s.sage_link().is_isotopic(k7snm)      # optional - snappy
            True

        using ``braid_notation``::

            sage: L2.link(use_item=L.items.braid_notation) == l2
            True

        observe::

            sage: L.link(use_item=L.items.braid_notation)
            Link with 2 components represented by 5 crossings

            sage: K6_1 = KnotInfo.K6_1
            sage: K6_1.link().braid() == K6_1.braid()
            False

        also observe::

            sage: K4_1 = KnotInfo.K4_1
            sage: K4_1.link().pd_code()
            [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]
            sage: K4_1.pd_notation()
            [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]

            sage: K5_1 = KnotInfo.K5_1
            sage: K5_1.link().braid()
            s^5
            sage: K5_1.link(K5_1.items.dt_notation).braid()
            s^-5
            sage: K5_1.link(K5_1.items.gauss_notation).braid()
            s^-5
        """
        if not isinstance(use_item, KnotInfoColumns):
            raise TypeError('%s must be an instance of %s' % (use_item, KnotInfoColumns))

        if snappy:
            try:
                from snappy import Link
            except ImportError:
                raise ImportError('this option demands snappy to be installed')
        elif self.is_knot():
            from sage.knots.knot import Knot as Link
        else:
            from sage.knots.link import Link

        if use_item == self.items.pd_notation:
            return Link(self.pd_notation())
        elif use_item == self.items.braid_notation:
            return Link(self.braid())
        elif use_item == self.items.name and snappy:
            if not self.is_knot():
                use_item = self.items.name_unoriented
            return Link(self[use_item])
        elif self.is_knot() and not snappy:
            # Construction via Gauss and DT-Code only possible for knots
            from sage.knots.knot import Knots
            if use_item == self.items.dt_notation:
                return Knots().from_dowker_code(self.dt_notation())
            elif use_item == self.items.gauss_notation:
                return Knots().from_gauss_code(self.gauss_notation())

        raise ValueError('link construction using %s not possible' % use_item)

    @cached_method
    def is_unique(self):
        r"""
        Return whether there is no other isotopic link in the database or not.

        OUTPUT: boolean or ``None`` if this cannot be determined

        EXAMPLES::

            sage: KnotInfo.L4a1_0.is_unique()
            True
            sage: KnotInfo.L5a1_0.is_unique()
            False
            sage: L = KnotInfo.L9a43_0_1             # optional - database_knotinfo
            sage: L.series(oriented=True).inject()   # optional - database_knotinfo
            Defining L9a43
            sage: [(L,L.is_unique()) for L in L9a43] # optional - database_knotinfo
            [(<KnotInfo.L9a43_0_0: 'L9a43{0,0}'>, True),
             (<KnotInfo.L9a43_1_0: 'L9a43{1,0}'>, False),
             (<KnotInfo.L9a43_0_1: 'L9a43{0,1}'>, None),
             (<KnotInfo.L9a43_1_1: 'L9a43{1,1}'>, False)]
        """
        # an isotopic pair must have the same unoriented name. So, we can focus
        # on such series
        if self.is_knot():
            return True
        S = self.series(oriented=True)
        hp = self.homfly_polynomial()
        Sl = S.list(homfly=hp)
        if len(Sl) == 1:
            return True
        kp = self.kauffman_polynomial()
        Sl = [L for L in Sl if L != self and L.kauffman_polynomial() == kp]
        if not Sl:
            return True

        b = self.braid()
        for L in Sl:
            Lb = L.braid()
            if L.braid() == b:
                return False
            if Lb.is_conjugated(b):
                return False

        return None

    @cached_method
    def is_recoverable(self, unique=True):
        r"""
        Return if ``self`` can be recovered from its conversion to Sage links
        using the ``pd_notation`` and the ``braid_notation`` and their
        mirror images.

        The method is indirectly used by the ``TestSuite`` of the series of ``self``.

        INPUT:

        - ``unique`` -- boolean (default: ``True``); if set to ``False``
          it is only checked if ``self`` is among the recovered items

        EXAMPLES::

            sage: KnotInfo.L4a1_0.inject()
            Defining L4a1_0
            sage: L4a1_0.is_recoverable()
            True
            sage: L4a1_0.is_recoverable(unique=False)
            True
            sage: KnotInfo.L5a1_0.inject()
            Defining L5a1_0
            sage: L5a1_0.is_recoverable()
            False
            sage: L5a1_0.is_recoverable(unique=False)
            True

        TESTS:

            sage: KnotInfo.K12a_165.is_recoverable(unique=False)  # optional - database_knotinfo, long time
            True
        """
        def recover(sym_mut, braid):
            r"""
            Check if ``self`` can be recovered form its associated
            Sage link.
            """
            if braid:
                l = self.link(self.items.braid_notation)
            else:
                l = self.link()
            if sym_mut is SymmetryMutant.mirror_image:
                l = l.mirror_image()
            elif sym_mut is SymmetryMutant.reverse:
                l = l.reverse()
            elif sym_mut is SymmetryMutant.concordance_inverse:
                l = l.mirror_image().reservse()

            def check_result(res):
                r"""
                Check a single result from ``get_knotinfo``.
                """
                if type(res) is tuple:
                    L, s = res
                else:
                    L, s = res.to_knotinfo()[0]
                if not isinstance(L, KnotInfoBase):
                    return False
                if L != self:
                    return False
                return s == sym_mut

            try:
                res = l.get_knotinfo(unique=unique)
            except NotImplementedError:
                return False
            if unique:
                return check_result(res)
            else:
                return any(check_result(r) for r in res)

        from sage.misc.misc import some_tuples
        if SymmetryMutant.unknown.matches(self):
            sym_muts = [SymmetryMutant.unknown]
        else:
            sym_muts = [s for s in SymmetryMutant if s.is_minimal(self)]
        return all(recover(sym, braid) for sym, braid in some_tuples(sym_muts, 2, 8))

    def inject(self, verbose=True):
        """
        Inject ``self`` with its name into the namespace of the
        Python code from which this function is called.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``); whether to suppress
          the message printed on the invocation

        EXAMPLES::

            sage: KnotInfo.K5_2.inject()
            Defining K5_2
            sage: K5_2.is_alternating()
            True
        """
        name = self.name
        if verbose:
            print("Defining %s" % (name))
        from sage.repl.user_globals import set_global
        set_global(name, self)

    @cached_method
    def series(self, oriented=False):
        r"""
        Return the series of links ``self`` belongs to.

        INPUT:

        - ``oriented`` -- boolean (default: ``False``); it only affects proper
          links. By default the items of the series will be again series of
          links collecting all orientation mutants of an unoriented name. To
          obtain the series of the individual links this keyword has to be set
          to ``True``.

        EXAMPLES::

            sage: K5 = KnotInfo.K5_2.series()
            sage: K5(1)
            <KnotInfo.K5_1: '5_1'>
            sage: KnotInfo.L4a1_1.series().inject()
            Defining L4a
            sage: L4a(1)
            Series of links L4a1
            sage: KnotInfo.L4a1_1.series(oriented=True).inject()
            Defining L4a1
            sage: L4a(1) == L4a1
            True
            sage: L4a1(1)
            <KnotInfo.L4a1_1: 'L4a1{1}'>
        """
        if oriented:
            S = KnotInfoSeries(self.crossing_number(), self.is_knot(), self.is_alternating(), self.name_unoriented())
        else:
            S = KnotInfoSeries(self.crossing_number(), self.is_knot(), self.is_alternating())
        return S

    def diagram(self, single=False, new=0, autoraise=True):
        r"""
        Launch the diagram of ``self`` given on the KnotInfo web-page.

        INPUT:

        - ``single`` -- boolean (default: ``False``); if set to ``True`` only one
          diagram is shown
        - ``new`` -- integer according to :func:`open` of :mod:`webbrowser`
          (``0`` default, ``1`` new window, ``2`` new tab)
        - ``autoraise`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: K = KnotInfo.K3_1
            sage: K.diagram()            # not tested
            True
            sage: K.diagram(single=True) # not tested
            True
        """
        import webbrowser
        if self.is_knot():
            filename = db.filename.knots
        else:
            filename = db.filename.links

        if single:
            return webbrowser.open(filename.diagram_url(self[self.items.diagram], single=single), new=new, autoraise=autoraise)
        else:
            return webbrowser.open(filename.diagram_url(self[self.items.name]), new=new, autoraise=autoraise)

    def knot_atlas_webpage(self, new=0, autoraise=True):
        r"""
        Launch the Knot Atlas web-page for ``self``.

        INPUT:

        - ``new`` -- integer according to :func:`open` of :mod:`webbrowser`
          (``0`` default, ``1`` new window, ``2`` new tab)
        - ``autoraise`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: K = KnotInfo.K3_1
            sage: K.knot_atlas_webpage()        # not tested
            True
        """
        import webbrowser
        return webbrowser.open(self[self.items.knot_atlas_anon], new=new, autoraise=autoraise)

    def knotilus_webpage(self, new=0, autoraise=True):
        r"""
        Launch the Knotilus web-page for ``self``.

        INPUT:

        - ``new`` -- integer according to :func:`open` of :mod:`webbrowser`
          (``0`` default, ``1`` new window, ``2`` new tab)
        - ``autoraise`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: K = KnotInfo.K3_1
            sage: K.knotilus_webpage(new=1)   # not tested
            True
        """
        import webbrowser
        return webbrowser.open(self[self.items.knotilus_page_anon], new=new, autoraise=autoraise)


# --------------------------------------------------------------------------------------------
# KnotInfoSeries
# --------------------------------------------------------------------------------------------
class KnotInfoSeries(UniqueRepresentation, SageObject):
    r"""
    This class can be used to access knots and links via their index
    according to the series they belong to.

    INPUT:

    - ``crossing_number`` -- integer giving the crossing number of this series
      of links
    - ``is_knot`` -- boolean; whether this series is a series of knots
      or proper links
    - ``is_alternating`` -- boolean; whether this series is restricted to
      alternating links or not
      This is not relevant for knots with less than 11 crossings
    - ``name_unoriented`` -- string restricting the series to all links with
      that ``name_unoriented``

    EXAMPLES::

        sage: from sage.knots.knotinfo import KnotInfoSeries
        sage: K6 = KnotInfoSeries(6, True, True); K6
        Series of knots K6
        sage: K6(3)
        <KnotInfo.K6_3: '6_3'>
        sage: list(K6)
        [<KnotInfo.K6_1: '6_1'>, <KnotInfo.K6_2: '6_2'>, <KnotInfo.K6_3: '6_3'>]
        sage: L6a = KnotInfoSeries(6, False, True); L6a
        Series of links L6a
        sage: L6a(2)
        Series of links L6a2
        sage: _.inject()
        Defining L6a2
        sage: list(L6a2)
        [<KnotInfo.L6a2_0: 'L6a2{0}'>, <KnotInfo.L6a2_1: 'L6a2{1}'>]
        sage: L6a2(0).series() == L6a
        True
        sage: L6a2(0) == L6a2('0')
        True
    """

    def __init__(self, crossing_number, is_knot, is_alternating, name_unoriented=None):
        r"""
        Python constructor.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: L6a = KnotInfoSeries(6, False, True); L6a
            Series of links L6a
        """
        self._crossing_number = crossing_number
        self._is_knot = is_knot
        self._is_alternating = is_alternating
        self._name_unoriented = name_unoriented

    @cached_method
    def list(self, oriented=False, comp=None, det=None, homfly=None):
        r"""
        Return this series as a Python list.

        INPUT:

        - ``oriented`` -- boolean (default: ``False``); it only affects
          series of proper links. By default the list items of a series of proper
          links are again series of links collecting all orientation types of an
          unoriented name. To obtain the list of the individual links this
          keyword has to be set to ``True``.

        - ``comp`` -- (default: ``None``) if given an integer for this
          keyword the list is restriced to links having the according number
          of components. This keyword implies ``oriented=True``.

        - ``det`` -- (default: ``None``) if given an integer for this
          keyword the list is restriced to links having the according value
          for its determinant. This keyword implies ``oriented=True``.

        - ``homfly`` -- (default: ``None``) if given a HOMFLY-PT polynomial
          having ``normalization='vz'`` for this keyword the list is restriced to
          links having the according value for its HOMFLY-PT polynomial. This
          keyword implies ``oriented=True``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: K6 = KnotInfoSeries(6, True, True); K6
            Series of knots K6
            sage: K6.list()
            [<KnotInfo.K6_1: '6_1'>, <KnotInfo.K6_2: '6_2'>, <KnotInfo.K6_3: '6_3'>]
            sage: KnotInfoSeries(2, False, True).inject()
            Defining L2a
            sage: L2a.list()
            [Series of links L2a1]
            sage: L2a.list(oriented=True)
            [<KnotInfo.L2a1_0: 'L2a1{0}'>, <KnotInfo.L2a1_1: 'L2a1{1}'>]
        """
        if homfly is not None:
            # additional restriction to number of components, determinant and
            # HOMFLY-PT polynomial
            l = self.list(oriented=True, comp=comp, det=det)
            return [L for L in l if L.homfly_polynomial() == homfly]

        if det is not None:
            # additional restriction to number of components and determinant
            l = self.list(oriented=True, comp=comp)
            return [L for L in l if L.determinant() == det]

        if comp is not None:
            # additional restriction to number of components
            l = self.list(oriented=True)
            return [L for L in l if L.num_components() == comp]

        # default case
        is_knot = self._is_knot
        cross_nr = self._crossing_number
        is_alt = self._is_alternating
        n_unori = self._name_unoriented

        res = []
        curr_n_unori = None
        for K in KnotInfo:
            if K.is_knot() != is_knot:
                continue
            if K.crossing_number() != cross_nr:
                continue
            if not is_knot or cross_nr > 10:
                if K.is_alternating() != is_alt:
                    continue
            if is_knot or oriented:
                res.append(K)
            else:
                this_n_unori = K.name_unoriented()
                if n_unori:
                    if this_n_unori != n_unori:
                        continue
                    res.append(K)
                elif this_n_unori != curr_n_unori:
                    if curr_n_unori:
                        res.append(KnotInfoSeries(cross_nr, is_knot, is_alt, curr_n_unori))
                    curr_n_unori = this_n_unori
                else:
                    continue

        if curr_n_unori:
            res.append(KnotInfoSeries(cross_nr, is_knot, is_alt, curr_n_unori))
        return res

    @cached_method
    def lower_list(self, oriented=False, comp=None, det=None, homfly=None):
        r"""
        Return this series together with all series with smaller crossing number
        as a Python list.

        INPUT:

        - ``oriented`` -- boolean (default: ``False``); see the
          description for :meth:`list`

        - ``comp`` -- (default: ``None``) see the description for
          :meth:`list`

        - ``det`` -- (default: ``None``) see the description for
          :meth:`list`

        - ``homfly`` -- (default: ``None``) see the description for
          :meth:`list`

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(5, True, True).lower_list()
            [<KnotInfo.K0_1: '0_1'>,
             <KnotInfo.K3_1: '3_1'>,
             <KnotInfo.K4_1: '4_1'>,
             <KnotInfo.K5_1: '5_1'>,
             <KnotInfo.K5_2: '5_2'>]
            sage: KnotInfoSeries(4, False, True).lower_list()
            [Series of links L2a1, Series of links L4a1]
            sage: KnotInfoSeries(4, False, True).lower_list(oriented=True)
            [<KnotInfo.L2a1_0: 'L2a1{0}'>,
             <KnotInfo.L2a1_1: 'L2a1{1}'>,
             <KnotInfo.L4a1_0: 'L4a1{0}'>,
             <KnotInfo.L4a1_1: 'L4a1{1}'>]
        """
        l = []
        cr = self._crossing_number
        if cr > 0:
            LS = type(self)(cr - 1, self._is_knot, self._is_alternating, self._name_unoriented)
            l = LS.lower_list(oriented=oriented, comp=comp, det=det, homfly=homfly)
        return l + self.list(oriented=oriented, comp=comp, det=det, homfly=homfly)

    def __repr__(self):
        r"""
        Return the representation string of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True)
            Series of knots K6
            sage: _.__repr__()
            'Series of knots K6'
        """
        if self._is_knot:
            return 'Series of knots %s' % (self._name())
        else:
            return 'Series of links %s' % (self._name())

    def __getitem__(self, item):
        r"""
        Return the given ``item`` from the list of ``self``
        (making the Python build-in ``list`` work).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True).inject()
            Defining K6
            sage: list(K6)                      # indirect doctest
            [<KnotInfo.K6_1: '6_1'>, <KnotInfo.K6_2: '6_2'>, <KnotInfo.K6_3: '6_3'>]
        """
        from sage.rings.integer import Integer
        if type(item) not in (int, Integer):
            raise ValueError('item must be an integer')
        l = self.list()
        max_item = len(l)
        if item < 0 or item > max_item:
            raise ValueError('item must be nonnegative and smaller than %s' % (max_item))

        return l[item]

    def __call__(self, item):
        r"""
        Return the given ``item`` from the list of ``self``
        (making the function call for ``self`` work).
        In contrast to ``__getitem__`` the first ``item``
        has to be ``1`` (not ``0``).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True).inject()
            Defining K6
            sage: K6(2)                         # indirect doctest
            <KnotInfo.K6_2: '6_2'>

            sage: # optional - database_knotinfo
            sage: KnotInfo.L8a21_0_1_0.inject()
            Defining L8a21_0_1_0
            sage: L8a21_0_1_0.series().inject()
            Defining L8a
            sage: L8a(1)
            Series of links L8a1
            sage: L8a(21)(2)     == L8a21_0_1_0
            True
            sage: L8a(21)('010') == L8a21_0_1_0
            True
        """
        if self._name_unoriented:
            if isinstance(item, str):
                # allow input as dual number according to naming
                item = int(item, 2)
            return self[item]

        from sage.rings.integer import Integer
        if type(item) not in (int, Integer):
            raise ValueError('item must be an integer')
        l = self.list()
        max_item = len(l) + 1
        if item < 1 or item > max_item:
            raise ValueError('item must be positive and smaller than %s' % (max_item))

        return l[item-1]

    def _name(self):
        r"""
        Return the name of the series.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True)._name()
            'K6'
        """
        is_knot = self._is_knot
        cross_nr = self._crossing_number
        is_alt = self._is_alternating
        n_unori = self._name_unoriented

        alt = 'a'
        if not is_alt:
            alt = 'n'

        if is_knot:
            if cross_nr > 10:
                res = 'K%s%s' % (cross_nr, alt)
            else:
                res = 'K%s' % (cross_nr)
        elif n_unori:
            res = '%s' % (n_unori)
        else:
            res = 'L%s%s' % (cross_nr, alt)
        return res

    def is_recoverable(self, unique=True, max_samples=8):
        r"""
        Return if all items of ``self`` can be recovered from its conversion to
        Sage links using the ``pd_notation`` and the ``braid_notation`` and their
        mirror images.

        The method is indirectly used by the ``TestSuite``.

        INPUT:

        - ``unique`` -- boolean (default: ``True``); see
          :meth:`KnotInfoBase.is_recoverable`
        - ``max_samples`` -- nonnegative integer or ``infinity``
          (default: `8`); limits the number of items to check (random sample).
          If set to ``infinity`` then no limit is set.

        EXAMPLES::

            sage: KnotInfo.L4a1_0.series().inject()
            Defining L4a
            sage: L4a.is_recoverable()
            True
            sage: L4a.is_recoverable(unique=False)
            True
            sage: KnotInfo.L5a1_0.series().inject()
            Defining L5a
            sage: L5a.is_recoverable()
            False
            sage: L5a.is_recoverable(unique=False)
            True
        """
        from sage.misc.misc import some_tuples
        l = self.list(oriented=True)
        bound = len(l)
        return all(L.is_recoverable(unique=unique) for L, in some_tuples(l, 1, bound, max_samples=max_samples))

    def _test_recover(self, **options):
        r"""
        Method used by ``TestSuite``. Tests if all links of the series can be
        recovered from their conversion to Sage links. It uses :meth:`is_recoverable`.
        Thus, per default maximal `8` items (random sample) are tested. Use the
        option ``max_samples`` to choose another limit or test all
        (``max_samples=infinity``)

        EXAMPLES::

            sage: TestSuite(KnotInfo.L5a1_0.series()).run(verbose=True)  # indirect doctest
            running ._test_category() . . . pass
            running ._test_new() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass
            running ._test_recover() . . . pass
            sage: TestSuite(KnotInfo.K6_1.series()).run(max_samples=infinity)  # indirect doctest
        """
        tester = options['tester']
        max_samples = tester._max_samples
        try:
            if max_samples:
                tester.assertTrue(self.is_recoverable(unique=False, max_samples=max_samples))
            else:
                tester.assertTrue(self.is_recoverable(unique=False))
        except ImportError:
            pass

    def inject(self, verbose=True):
        r"""
        Inject ``self`` with its name into the namespace of the
        Python code from which this function is called.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``); to suppress
          the message printed on the invocation

        EXAMPLES::

            sage: KnotInfoSeries(6, True, True).inject()
            Defining K6
            sage: K6(2)
            <KnotInfo.K6_2: '6_2'>
        """
        name = self._name()
        if verbose:
            print("Defining %s" % (name))
        from sage.repl.user_globals import set_global
        set_global(name, self)


KnotInfo = KnotInfoBase('KnotInfo', db.row_names())
