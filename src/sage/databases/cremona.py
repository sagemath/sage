"""
Cremona's tables of elliptic curves

Sage includes John Cremona's tables of elliptic curves in an
easy-to-use format. An instance of the class CremonaDatabase()
gives access to the database.

If the optional full CremonaDatabase is not installed, a mini-version
is included by default with Sage.  It contains Weierstrass equations,
rank, and torsion for curves up to conductor 10000.

The large database includes all curves in John Cremona's tables. It also
includes data related to the BSD conjecture and modular degrees for all of
these curves, and generators for the Mordell-Weil groups. To install it via the
optional :ref:`database_cremona_ellcurve <spkg_database_cremona_ellcurve>`
package, run the following command in the shell ::

    sage -i database_cremona_ellcurve

This causes the latest version of the database to be downloaded from
the internet.

Both the mini and full versions of John Cremona's tables are stored in
SAGE_SHARE/cremona as SQLite databases. The mini version has the layout::

    CREATE TABLE t_class(conductor INTEGER, class TEXT PRIMARY KEY, rank INTEGER);
    CREATE TABLE t_curve(class TEXT, curve TEXT PRIMARY KEY, eqn TEXT UNIQUE, tors INTEGER);
    CREATE INDEX i_t_class_conductor ON t_class(conductor);
    CREATE INDEX i_t_curve_class ON t_curve(class);

while the full version has the layout::

    CREATE TABLE t_class(conductor INTEGER, class TEXT PRIMARY KEY, rank INTEGER, L REAL, deg INTEGER);
    CREATE TABLE t_curve(class TEXT, curve TEXT PRIMARY KEY, eqn TEXT UNIQUE, gens TEXT, tors INTEGER, cp INTEGER, om REAL, reg REAL, sha);
    CREATE INDEX i_t_class_conductor ON t_class(conductor);
    CREATE INDEX i_t_curve_class ON t_curve(class);
"""
# ****************************************************************************
#       Copyright (C) 2014 John Cremona <john.cremona@gmail.com>
#       Copyright (C) 2011 R. Andrew Ohana <andrew.ohana@gmail.com>
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
from sage.misc.prandom import randint

import sage.schemes.elliptic_curves.constructor as elliptic
from .sql_db import SQLDatabase, verify_column
from sage.features.databases import DatabaseCremona
from sage.misc.timing import walltime

import re
import string

_cremonaSkeleton = {
    't_class': {
        'conductor': {'sql':'INTEGER', 'index':True},
        'class':     {'sql':'TEXT',    'primary_key':True},
        'rank':      {'sql':'INTEGER'},
        'L':         {'sql':'REAL'},
        'deg':       {'sql':'INTEGER'}
    },
    't_curve': {
        'class':    {'sql':'TEXT', 'index':True},
        'curve':    {'sql':'TEXT', 'primary_key':True},
        'eqn':      {'sql':'TEXT', 'unique':True},
        'gens':     {'sql':'TEXT'},
        'tors':     {'sql':'INTEGER'},
        'cp':       {'sql':'INTEGER'},
        'om':       {'sql':'REAL'},
        'reg':      {'sql':'REAL'},
        'sha':      {'sql':'NOTYPE'}
    }
}
_miniCremonaSkeleton = {
    't_class': {
        'conductor': {'sql':'INTEGER', 'index':True},
        'class':     {'sql':'TEXT',    'primary_key':True},
        'rank':      {'sql':'INTEGER'}
    },
    't_curve': {
        'class':    {'sql':'TEXT', 'index':True},
        'curve':    {'sql':'TEXT', 'primary_key':True},
        'eqn':      {'sql':'TEXT', 'unique':True},
        'tors':     {'sql':'INTEGER'}
    }
}

for t in _cremonaSkeleton:
    for c in _cremonaSkeleton[t]:
        _cremonaSkeleton[t][c] = verify_column(_cremonaSkeleton[t][c])
    for c in _miniCremonaSkeleton[t]:
        _miniCremonaSkeleton[t][c] = verify_column(_miniCremonaSkeleton[t][c])


def build(name, data_tgz, largest_conductor=0, mini=False, decompress=True):
    """
    Build the CremonaDatabase with given name from scratch
    using the data_tgz tarball.

    .. NOTE::

           For data up to level 350000, this function takes about
           3m40s.  The resulting database occupies 426MB disk space.

    To create the large Cremona database from Cremona's data_tgz
    tarball, obtainable from
    http://homepages.warwick.ac.uk/staff/J.E.Cremona/ftp/data/, run
    the following command::

        sage: d = sage.databases.cremona.build('cremona','ecdata.tgz')   # not tested
    """
    from sage.env import SAGE_SHARE
    db_path = os.path.join(SAGE_SHARE,'cremona',name.replace(' ','_')+'.db')
    if os.path.exists(db_path):
        raise RuntimeError('Please (re)move %s before building ' % db_path
                + 'database')
    if not os.path.exists(data_tgz):
        raise OSError("The data file is not at %s" % data_tgz)
    t = walltime()

    if decompress:
        cmd = "tar zxf {} ecdata/allcurves ecdata/allbsd ecdata/degphi ecdata/allgens".format(data_tgz)
        print("Extracting files from {}...".format(data_tgz))
        n = os.system(cmd)
        if n:
            raise RuntimeError("Error extracting tarball.")
        else:
            print("...finished file extraction")
    if mini:
        c = MiniCremonaDatabase(name,False,True)
    else:
        c = LargeCremonaDatabase(name,False,True)
    # The following line assumes that the tarball extracts to a
    # directory called 'ecdata'
    c._init_from_ftpdata('ecdata', largest_conductor)
    print("Total time: ", walltime(t))


def is_optimal_id(id):
    """
    Return ``True`` if the Cremona id refers to an optimal curve, and
    ``False`` otherwise.

    The curve is optimal if the id, which is of the
    form [letter code][number] has number 1.

    .. NOTE::

       990h3 is the optimal curve in that class, so doesn't obey
       this rule.

    INPUT:

    - ``id`` -- string of form letter code followed by an
      integer, e.g., a3, bb5, etc.

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.databases.cremona import is_optimal_id
        sage: is_optimal_id('b1')
        True
        sage: is_optimal_id('bb1')
        True
        sage: is_optimal_id('c1')
        True
        sage: is_optimal_id('c2')
        False
    """
    return id[-1] == '1' and not id[-2].isdigit()


def cremona_letter_code(n):
    """
    Return the Cremona letter code corresponding to an integer.

    For example, 0 - a 25 - z 26 - ba 51 - bz 52 - ca 53 - cb etc.

    .. NOTE::

       This is just the base 26 representation of n, where a=0, b=1,
       ..., z=25. This extends the old Cremona notation (counting from
       0) for the first 26 classes, and is different for classes above
       26.

    INPUT:

    - ``n`` -- nonnegative integer

    OUTPUT: string

    EXAMPLES::

        sage: from sage.databases.cremona import cremona_letter_code
        sage: cremona_letter_code(0)
        'a'
        sage: cremona_letter_code(26)
        'ba'
        sage: cremona_letter_code(27)
        'bb'
        sage: cremona_letter_code(521)
        'ub'
        sage: cremona_letter_code(53)
        'cb'
        sage: cremona_letter_code(2005)
        'czd'

    TESTS::

        sage: cremona_letter_code(QQ)
        Traceback (most recent call last):
        ...
        ValueError: Cremona letter codes are only defined for nonnegative integers
        sage: cremona_letter_code(x)                                                    # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: Cremona letter codes are only defined for nonnegative integers
        sage: cremona_letter_code(-1)
        Traceback (most recent call last):
        ...
        ValueError: Cremona letter codes are only defined for nonnegative integers
        sage: cremona_letter_code(3.14159)
        Traceback (most recent call last):
        ...
        ValueError: Cremona letter codes are only defined for nonnegative integers
    """
    try:
        m = int(n)
        if n == m:
            n = m
        else:
            n = -1
    except (ValueError, TypeError):
        n = -1

    if n < 0:
        raise ValueError("Cremona letter codes are only defined for nonnegative integers")

    if n == 0:
        return "a"
    s = ""
    while n != 0:
        s = chr(n % 26+97) + s
        n //= 26
    return s


def old_cremona_letter_code(n):
    r"""
    Return the *old* Cremona letter code corresponding to an integer.

    For example::

        1  --> A
        26 --> Z
        27 --> AA
        52 --> ZZ
        53 --> AAA
        etc.

    INPUT:

    - ``n`` -- integer

    OUTPUT: string

    EXAMPLES::

        sage: from sage.databases.cremona import old_cremona_letter_code
        sage: old_cremona_letter_code(1)
        'A'
        sage: old_cremona_letter_code(26)
        'Z'
        sage: old_cremona_letter_code(27)
        'AA'
        sage: old_cremona_letter_code(521)
        'AAAAAAAAAAAAAAAAAAAAA'
        sage: old_cremona_letter_code(53)
        'AAA'
        sage: old_cremona_letter_code(2005)
        'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    """
    n -= 1
    k = n % 26 + 65
    label = chr(k)*int(n//26 + 1)
    return label


old_cremona_label_regex = re.compile(r'(\d+)([A-Z]*)(\d*)$')
cremona_label_regex = re.compile(r'(\d+)([a-z]*)(\d*)$')
lmfdb_label_regex = re.compile(r'(\d+)\.([a-z]+)(\d*)$')


def parse_cremona_label(label, numerical_class_code=False):
    """
    Given a Cremona label that defines an elliptic
    curve, e.g., 11a1 or 37b3, parse the label and return the
    conductor, isogeny class label, and number.

    For this function, the curve number may be omitted, in which case
    it defaults to 1.  If the curve number and isogeny class are both
    omitted (label is just a string representing a conductor), then
    the isogeny class defaults to 'a' and the number to 1.  Valid
    labels consist of one or more digits, followed by zero or more
    letters (either all in upper case for an old Cremona label, or all
    in lower case), followed by zero or more digits.

    INPUT:

    - ``label`` -- string; a valid Cremona elliptic curve label

    - ``numerical_class_code`` -- boolean (default: ``False``); if ``True``,
      convert the isogeny class label from a letter code in base 26
      to an integer.  This is useful for sorting.

    OUTPUT:

    -  integer; the conductor
    -  string or integer; the isogeny class label
    -  integer; the number

    EXAMPLES::

        sage: from sage.databases.cremona import parse_cremona_label
        sage: parse_cremona_label('37a2')
        (37, 'a', 2)
        sage: parse_cremona_label('37b1')
        (37, 'b', 1)
        sage: parse_cremona_label('10bb2')
        (10, 'bb', 2)
        sage: parse_cremona_label('11a')
        (11, 'a', 1)
        sage: parse_cremona_label('11')
        (11, 'a', 1)

    Valid old Cremona labels are allowed::

        sage: parse_cremona_label('17CCCC')
        (17, 'dc', 1)
        sage: parse_cremona_label('5AB2')
        Traceback (most recent call last):
        ...
        ValueError: 5AB2 is not a valid Cremona label

    When ``numerical_class_code`` is ``True``, the output is a triple of integers::

        sage: from sage.databases.cremona import parse_cremona_label
        sage: parse_cremona_label('100800hj2')
        (100800, 'hj', 2)
        sage: parse_cremona_label('100800hj2', numerical_class_code=True)
        (100800, 191, 2)

    TESTS::

        sage: from sage.databases.cremona import parse_cremona_label
        sage: parse_cremona_label('x11')
        Traceback (most recent call last):
        ...
        ValueError: x11 is not a valid Cremona label
    """
    m = cremona_label_regex.match(str(label))
    if m is None:
        m = old_cremona_label_regex.match(str(label))
        if m is None:
            raise ValueError(label + " is not a valid Cremona label")

    conductor, iso, num = m.groups()
    if len(iso) == 0:
        iso = "a"
    if len(num) == 0:
        num = "1"

    # convert old cremona labels to new ones
    if iso.upper() == iso and iso[0]*len(iso) == iso:
        iso = cremona_letter_code((len(iso)-1)*26+ord(iso[0])-ord('A'))

    # verify cremona label is valid
    if iso.lower() != iso:
        raise ValueError('%s is not a valid Cremona label' % label)

    # convert class label to an int if requested
    if numerical_class_code:
        iso = class_to_int(iso)

    return int(conductor), iso, int(num)


def parse_lmfdb_label(label, numerical_class_code=False):
    """
    Given an LMFDB label that defines an elliptic curve, e.g., 11.a1
    or 37.b3, parse the label and return the conductor, isogeny class
    label, and number.

    The LMFDB label (named after the `L`-functions and modular forms
    database), is determined by the following two orders:

    - Isogeny classes with the same conductor are ordered
      lexicographically by the coefficients in the `q`-expansion of the
      associated modular form.

    - Curves within the same isogeny class are ordered
      lexicographically by the a-invariants of the minimal model.

    The format is <conductor>.<iso><curve>, where the isogeny class is
    encoded using the same base-26 encoding into letters used in
    Cremona's labels.  For example, 990.h3 is the same as Cremona's 990j1

    For this function, the curve number may be omitted, in which case
    it defaults to 1.  If the curve number and isogeny class are both
    omitted (label is just a string representing a conductor), then
    the isogeny class defaults to 'a' and the number to 1.

    INPUT:

    - ``label`` -- str

    - ``numerical_class_code`` -- boolean (default: ``False``); if ``True``,
      convert the isogeny class label from a letter code in base 26
      to an integer.  This is useful for sorting.

    OUTPUT:

    - ``int`` -- the conductor
    - ``str`` or ``int`` -- the isogeny class label
    - ``int`` -- the number

    EXAMPLES::

        sage: from sage.databases.cremona import parse_lmfdb_label
        sage: parse_lmfdb_label('37.a2')
        (37, 'a', 2)
        sage: parse_lmfdb_label('37.b')
        (37, 'b', 1)
        sage: parse_lmfdb_label('10.bb2')
        (10, 'bb', 2)

    When ``numerical_class_code`` is ``True``, the output is a triple of integers::

        sage: from sage.databases.cremona import parse_lmfdb_label
        sage: parse_lmfdb_label('100800.bg4')
        (100800, 'bg', 4)
        sage: parse_lmfdb_label('100800.bg4', numerical_class_code=True)
        (100800, 32, 4)
    """
    m = lmfdb_label_regex.match(str(label).lower())
    if m is None:
        raise ValueError(label + " is not a valid LMFDB label")
    conductor, iso, num = m.groups()
    if len(iso) == 0:
        iso = "a"
    if len(num) == 0:
        num = "1"
    # convert class label to an int if requested
    if numerical_class_code:
        iso = class_to_int(iso)

    return int(conductor), iso, int(num)


_class_curve_re = re.compile(r'(?P<class>[a-z]+)(?P<curve>\d+)')


def split_code(key):
    """
    Split class + curve id string into its two parts.

    EXAMPLES::

        sage: import sage.databases.cremona as cremona
        sage: cremona.split_code('ba2')
        ('ba', '2')
        sage: cremona.split_code('42')
        Traceback (most recent call last):
        ...
        ValueError: invalid curve ID: '42'
    """
    m = _class_curve_re.match(key)

    if not m:
        raise ValueError("invalid curve ID: '{0}'".format(key))

    return (m.group('class'), m.group('curve'))


def class_to_int(k):
    """
    Convert class id string into an integer.

    Note that this is the inverse of :func:`cremona_letter_code`.

    EXAMPLES::

        sage: import sage.databases.cremona as cremona
        sage: cremona.class_to_int('ba')
        26
        sage: cremona.class_to_int('cremona')
        821863562
        sage: cremona.cremona_letter_code(821863562)
        'cremona'
    """
    kk = [string.ascii_lowercase.index(ch) for ch in list(k)]
    kk.reverse()
    return sum(kk[i] * 26 ** i for i in range(len(kk)))


def sort_key(key1):
    """
    Comparison key for curve id strings.

    .. NOTE::

       Not the same as standard lexicographic order!

    EXAMPLES::

        sage: from sage.databases.cremona import sort_key
        sage: l = ['ba1', 'z1']
        sage: sorted(l, key=sort_key)
        ['z1', 'ba1']
    """
    cl1, cu1 = split_code(key1)
    return (class_to_int(cl1), cu1)


def cremona_to_lmfdb(cremona_label, CDB=None):
    """
    Convert a Cremona label into an LMFDB label.

    See :func:`parse_lmfdb_label` for an explanation of LMFDB labels.

    INPUT:

    - ``cremona_label`` -- string, the Cremona label of a curve; this can be
      the label of a curve (e.g. '990j1') or of an isogeny class (e.g. '990j')
    - ``CDB`` -- the Cremona database in which to look up the isogeny
      classes of the same conductor

    OUTPUT: ``lmfdb_label``; string, the corresponding LMFDB label

    EXAMPLES::

        sage: from sage.databases.cremona import cremona_to_lmfdb, lmfdb_to_cremona
        sage: cremona_to_lmfdb('990j1')
        '990.h3'
        sage: lmfdb_to_cremona('990.h3')
        '990j1'

    TESTS::

        sage: for label in ['5077a1','66a3','102b','420c2']:
        ....:     assert(lmfdb_to_cremona(cremona_to_lmfdb(label)) == label)
        sage: for label in ['438.c2','306.b','462.f3']:
        ....:     assert(cremona_to_lmfdb(lmfdb_to_cremona(label)) == label)
    """
    from sage.libs.pari.all import pari
    m = cremona_label_regex.match(cremona_label)
    if m is None:
        raise ValueError("Invalid Cremona label")
    N, cremona_iso, cremona_number = m.groups()
    if CDB is None:
        CDB = CremonaDatabase()
    classes = CDB.isogeny_classes(N)
    ft = 53
    tff = 255 # This should be enough to distinguish between curves (using heuristics from Sato-Tate for example)
    isos = []
    for i, iso in enumerate(classes):
        alist = iso[0][0]
        E = pari(alist).ellinit(precision=ft)
        isos.append((E.ellan(tff, python_ints=True), cremona_letter_code(i)))
    isos.sort()
    sorted_letters = [iso[1] for iso in isos]
    lmfdb_iso = cremona_letter_code(sorted_letters.index(cremona_iso))
    if len(cremona_number) > 0:
        iso_class = sorted([(curve[0],str(i+1)) for i,curve in enumerate(classes[class_to_int(cremona_iso)])])
        sorted_numbers = [curve[1] for curve in iso_class]
        lmfdb_number = str(sorted_numbers.index(cremona_number)+1)
        return N + '.' + lmfdb_iso + lmfdb_number
    else:
        return N + '.' + lmfdb_iso


def lmfdb_to_cremona(lmfdb_label, CDB=None):
    """
    Convert an LMFDB label into a Cremona label.

    See :func:`parse_lmfdb_label` for an explanation of LMFDB labels.

    INPUT:

    - ``lmfdb_label`` -- string, the LMFDB label of a curve; this can be the
      label of a curve (e.g. '990.j1') or of an isogeny class (e.g. '990.j')
    - ``CDB`` -- the Cremona database in which to look up the isogeny
      classes of the same conductor

    OUTPUT: ``cremona_label``; a string, the corresponding Cremona label

    EXAMPLES::

        sage: from sage.databases.cremona import cremona_to_lmfdb, lmfdb_to_cremona
        sage: lmfdb_to_cremona('990.h3')
        '990j1'
        sage: cremona_to_lmfdb('990j1')
        '990.h3'
    """
    from sage.libs.pari.all import pari
    m = lmfdb_label_regex.match(lmfdb_label)
    if m is None:
        raise ValueError("Invalid LMFDB label")
    N, lmfdb_iso, lmfdb_number = m.groups()
    if CDB is None:
        CDB = CremonaDatabase()
    classes = CDB.isogeny_classes(N)
    ft = 53
    tff = 255 # This should be enough to distinguish between curves (using heuristics from Sato-Tate for example)
    isos = []
    for i, iso in enumerate(classes):
        alist = iso[0][0]
        E = pari(alist).ellinit(precision=ft)
        isos.append((E.ellan(tff, python_ints=True), cremona_letter_code(i)))
    isos.sort()
    cremona_iso = isos[class_to_int(lmfdb_iso)][1]
    if len(lmfdb_number) > 0:
        iso_class = sorted([(curve[0],i+1) for i,curve in enumerate(classes[class_to_int(cremona_iso)])])
        cremona_number = str(iso_class[int(lmfdb_number)-1][1])
        return N + cremona_iso + cremona_number
    else:
        return N + cremona_iso


class MiniCremonaDatabase(SQLDatabase):
    """
    The Cremona database of elliptic curves.

    EXAMPLES::

        sage: c = CremonaDatabase()
        sage: c.allcurves(11)
        {'a1': [[0, -1, 1, -10, -20], 0, 5],
         'a2': [[0, -1, 1, -7820, -263580], 0, 1],
         'a3': [[0, -1, 1, 0, 0], 0, 5]}
    """
    _expected_skeleton = _miniCremonaSkeleton

    def __init__(self, name, read_only=True, build=False):
        """
        Initialize the database.

        TESTS::

            sage: c = CremonaDatabase('cremona mini')
            sage: c.name
            'cremona mini'
            sage: c = CremonaDatabase('cremona')    # optional - database_cremona_ellcurve
            sage: c.name                            # optional - database_cremona_ellcurve
            'cremona'
        """
        self.name = name
        name = name.replace(' ', '_')
        db_path = DatabaseCremona(name=name).absolute_filename()
        if build:
            if read_only:
                raise RuntimeError('The database must not be read_only.')
            SQLDatabase.__init__(self, db_path, read_only=read_only,
                    skeleton=self._expected_skeleton)
            return
        SQLDatabase.__init__(self, db_path, read_only=read_only)
        if self.get_skeleton() != self._expected_skeleton:
            raise RuntimeError('Database at %s does ' % (self.__dblocation__)
              + 'not appear to be a valid SQL Cremona database.')

    def __iter__(self):
        """
        Return an iterator through all EllipticCurve objects in the
        Cremona database.

        TESTS::

            sage: it = CremonaDatabase().__iter__()
            sage: next(it).label()
            '11a1'
            sage: next(it).label()
            '11a2'
            sage: next(it).label()
            '11a3'
            sage: next(it).label()
            '14a1'
            sage: skip = [next(it) for _ in range(100)]
            sage: next(it).label()
            '45a3'
        """
        query = "SELECT curve FROM t_curve,t_class USING(class) ORDER BY conductor"
        for c in self.__connection__.cursor().execute(query):
            yield self.elliptic_curve(c[0])

    def __getitem__(self, N):
        """
        If `N` is an integer, return all data about level `N` in the database.
        If `N` is a string it must be a Cremona label, in which case return
        the corresponding elliptic curve, if it is in the database.

        INPUT:

        - ``N`` -- integer or string

        OUTPUT: dictionary (if `N` is an integer) or EllipticCurve (if `N` is
        a string)

        TESTS::

            sage: c = CremonaDatabase()
            sage: c[11]['allcurves']['a2']
            [[0, -1, 1, -7820, -263580], 0, 1]
            sage: c['11a2']
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        """
        if isinstance(N, str):
            return self.elliptic_curve(N)

        try:
            N = int(N)
        except ValueError:
            raise KeyError("N (=%s) must be a string or positive integer." % N)

        if N <= 0:
            raise KeyError("N (=%s) must be a string or positive integer." % N)

        ret = {'allcurves': self.allcurves(N)}
        if hasattr(self, 'allbsd'):
            ret['allbsd'] = self.allbsd(N)
        if hasattr(self, 'degphi'):
            ret['degphi'] = self.degphi(N)
        if hasattr(self, 'allgens'):
            ret['allgens'] = self.allgens(N)
        return ret

    def __repr__(self):
        """
        String representation of this database.

        TESTS::

            sage: c = CremonaDatabase('cremona mini')
            sage: c.__repr__()
            "Cremona's database of elliptic curves with conductor at most 9999"
        """
        return "Cremona's database of elliptic curves with conductor at most "\
            + str(self.largest_conductor())

    def allcurves(self, N):
        """
        Return the allcurves table of curves of conductor N.

        INPUT:

        - ``N`` -- integer; the conductor

        OUTPUT: dictionary; id:[ainvs, rank, tor], ...

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.allcurves(11)['a3']
            [[0, -1, 1, 0, 0], 0, 5]
            sage: c.allcurves(12)
            {}
            sage: c.allcurves(12001)['a1']   # optional - database_cremona_ellcurve
            [[1, 0, 0, -101, 382], 1, 1]
        """
        ret = {}
        for c in self.__connection__.cursor().execute('SELECT curve,eqn,'
            + 'rank,tors FROM t_curve,t_class USING(class) WHERE '
            + 'conductor=?', (int(N),)):
            N,iso,num = parse_cremona_label(c[0])
            ret[iso+str(num)] = [eval(c[1]), c[2], c[3]]
        return ret

    def curves(self, N):
        """
        Return the curves table of all *optimal* curves of conductor N.

        INPUT:

        - ``N`` -- integer; the conductor

        OUTPUT: dictionary; id:[ainvs, rank, tor], ...

        EXAMPLES:

        Optimal curves of conductor 37::

            sage: CremonaDatabase().curves(37)
            {'a1': [[0, 0, 1, -1, 0], 1, 1], 'b1': [[0, 1, 1, -23, -50], 0, 3]}

        Note the 'h3', which is the unique case in the tables where
        the optimal curve doesn't have label ending in 1::

            sage: sorted(CremonaDatabase().curves(990))
            ['a1', 'b1', 'c1', 'd1', 'e1', 'f1', 'g1', 'h3', 'i1', 'j1', 'k1', 'l1']

        TESTS::

            sage: c = CremonaDatabase()
            sage: c.curves(12001)['a1']   # optional - database_cremona_ellcurve
            [[1, 0, 0, -101, 382], 1, 1]
        """
        ret = {}
        for c in self.__connection__.cursor().execute('SELECT curve,eqn,'
            + 'rank,tors FROM t_curve,t_class USING(class) WHERE '
            + 'curve=class||1 AND conductor=?',(int(N),)):
            N,iso,num = parse_cremona_label(c[0])
            ret[iso+str(num)] = [eval(c[1]),c[2],c[3]]
        if N == 990:
            del ret['h1']
            ret['h3'] = [[1,-1,1,-1568,-4669],1,6]
        return ret

    def coefficients_and_data(self, label):
        """
        Return the Weierstrass coefficients and other data for the
        curve with given label.

        EXAMPLES::

            sage: c, d = CremonaDatabase().coefficients_and_data('144b1')
            sage: c
            [0, 0, 0, 6, 7]
            sage: d['conductor']
            144
            sage: d['cremona_label']
            '144b1'
            sage: d['rank']
            0
            sage: d['torsion_order']
            2

        Check that :issue:`17904` is fixed::

            sage: 'gens' in CremonaDatabase().coefficients_and_data('100467a2')[1] # optional - database_cremona_ellcurve
            True
        """
        # There are two possible strings: the Cremona label and the LMFDB label.
        # They are distinguished by the presence of a period.
        if label.find('.') == -1:
            cremona_label = label
            lmfdb_label = None
        else:
            cremona_label = lmfdb_to_cremona(label)
            lmfdb_label = label

        N, iso, num = parse_cremona_label(cremona_label)
        label = str(N)+iso+str(num)
        if self.get_skeleton() == _miniCremonaSkeleton:
            q = self.__connection__.cursor().execute("SELECT eqn,rank,tors "
                + 'FROM t_curve,t_class USING(class) WHERE curve=?', (label,))
        else:
            q = self.__connection__.cursor().execute("SELECT eqn,rank,tors,"
                + "deg,gens,cp,om,L,reg,sha FROM t_curve,t_class "
                + "USING(class) WHERE curve=?",(label,))
        try:
            c = next(q)
        except StopIteration:
            if N < self.largest_conductor():
                message = "There is no elliptic curve with label " + label \
                    + " in the database"
            elif DatabaseCremona().is_present():
                message = "There is no elliptic curve with label " + label \
                    + " in the currently available databases"
            else:
                message = "There is no elliptic curve with label " \
                    + label + " in the default database; try installing " \
                    + "the optional package database_cremona_ellcurve which " \
                    + "contains the complete Cremona database"
            raise ValueError(message)
        ainvs = eval(c[0])
        data = {'cremona_label': label,
                'rank': c[1],
                'torsion_order': c[2],
                'conductor': N}
        if lmfdb_label:
            data['lmfdb_label'] = lmfdb_label
        if len(c) > 3:
            data['modular_degree'] = (c[3])
            data['gens'] = eval(c[4])
            data['db_extra'] = list(c[5:])
        return ainvs, data

    def data_from_coefficients(self, ainvs):
        """
        Return elliptic curve data for the curve with given
        Weierstrass coefficients.

        EXAMPLES::

            sage: d = CremonaDatabase().data_from_coefficients([1, -1, 1, 31, 128])
            sage: d['conductor']
            1953
            sage: d['cremona_label']
            '1953c1'
            sage: d['rank']
            1
            sage: d['torsion_order']
            2

        Check that :issue:`17904` is fixed::

            sage: ai = EllipticCurve('100467a2').ainvs() # optional - database_cremona_ellcurve
            sage: 'gens' in CremonaDatabase().data_from_coefficients(ai) # optional - database_cremona_ellcurve
            True
        """
        ainvs = str(list(ainvs))
        if self.get_skeleton() == _miniCremonaSkeleton:
            q = self.__connection__.cursor().execute("SELECT curve,rank,tors "
                + 'FROM t_curve,t_class USING(class) WHERE eqn=?',
                (ainvs.replace(' ', ''),))
        else:
            q = self.__connection__.cursor().execute("SELECT curve,rank,tors,"
                + "deg,gens,cp,om,L,reg,sha FROM t_curve,t_class "
                + "USING(class) WHERE eqn=?",
                (ainvs.replace(' ', ''),))
        try:
            c = next(q)
        except StopIteration:
            raise RuntimeError("There is no elliptic curve with coefficients "
                               + ainvs + " in the database")
        label = str(c[0])
        N, iso, num = parse_cremona_label(label)
        data = {'cremona_label': label,
                'rank': c[1],
                'torsion_order': c[2],
                'conductor': N}
        if len(c) > 3:
            data['modular_degree'] = (c[3])
            data['gens'] = eval(c[4])
            data['db_extra'] = list(c[5:])
        return data

    def elliptic_curve_from_ainvs(self, ainvs):
        """
        Return the elliptic curve in the database of with minimal ``ainvs``
        if it exists.

        This raises a :exc:`RuntimeError` exception otherwise.

        INPUT:

        - ``ainvs`` -- list (5-tuple of int's); the minimal
          Weierstrass model for an elliptic curve

        OUTPUT: EllipticCurve

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.elliptic_curve_from_ainvs([0, -1, 1, -10, -20])
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: c.elliptic_curve_from_ainvs([1, 0, 0, -101, 382])  # optional - database_cremona_ellcurve
            Elliptic Curve defined by y^2 + x*y = x^3 - 101*x + 382 over Rational Field

        Old (pre-2006) Cremona labels are also allowed::

            sage: c.elliptic_curve('9450KKKK1')
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 5*x + 7 over Rational Field

        Make sure :issue:`12565` is fixed::

            sage: c.elliptic_curve('10a1')
            Traceback (most recent call last):
            ...
            ValueError: There is no elliptic curve with label 10a1 in the database
        """
        data = self.data_from_coefficients(ainvs)
        return elliptic.EllipticCurve(ainvs, **data)

    def elliptic_curve(self, label):
        """
        Return an elliptic curve with given label with some data about it
        from the database pre-filled in.

        INPUT:

        - ``label`` -- string (Cremona or LMFDB label)

        OUTPUT: an :class:`sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field`

        .. NOTE::

            For more details on LMFDB labels see :func:`parse_lmfdb_label`.

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.elliptic_curve('11a1')
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: c.elliptic_curve('12001a1')    # optional - database_cremona_ellcurve
            Elliptic Curve defined by y^2 + x*y = x^3 - 101*x + 382 over Rational Field
            sage: c.elliptic_curve('48c1')
            Traceback (most recent call last):
            ...
            ValueError: There is no elliptic curve with label 48c1 in the database

        You can also use LMFDB labels::

            sage: c.elliptic_curve('462.f3')
            Elliptic Curve defined by y^2 + x*y = x^3 - 363*x + 1305 over Rational Field
        """
        ainvs, data = self.coefficients_and_data(label)
        return elliptic.EllipticCurve(ainvs, **data)

    def iter(self, conductors):
        """
        Return an iterator through all curves in the database with given
        conductors.

        INPUT:

        - ``conductors`` -- list or generator of ints

        OUTPUT: generator that iterates over EllipticCurve objects

        EXAMPLES::

            sage: [e.cremona_label() for e in CremonaDatabase().iter([11..15])]
            ['11a1', '11a2', '11a3', '14a1', '14a2', '14a3', '14a4', '14a5',
             '14a6', '15a1', '15a2', '15a3', '15a4', '15a5', '15a6', '15a7', '15a8']
        """
        for N in conductors:
            for c in self.__connection__.cursor().execute('SELECT curve '
                + 'FROM t_curve,t_class USING(class) WHERE conductor=?',
                (int(N),)):
                yield self.elliptic_curve(c[0])

    def isogeny_classes(self, conductor):
        """
        Return the allcurves data (ainvariants, rank and torsion) for the
        elliptic curves in the database of given conductor as a list of
        lists, one for each isogeny class. The curve with number 1 is
        always listed first.

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.isogeny_classes(11)
            [[[[0, -1, 1, -10, -20], 0, 5],
             [[0, -1, 1, -7820, -263580], 0, 1],
             [[0, -1, 1, 0, 0], 0, 5]]]
            sage: c.isogeny_classes(12001)   # optional - database_cremona_ellcurve
            [[[[1, 0, 0, -101, 382], 1, 1]],
             [[[0, 0, 1, -247, 1494], 1, 1]],
             [[[0, 0, 1, -4, -18], 1, 1]],
             [[[0, 1, 1, -10, 18], 1, 1]]]
        """
        conductor = int(conductor)
        classes = []
        A = self.allcurves(conductor)
        K = sorted(A, key=sort_key)
        for k in K:
            v = A[k]
            # test if not first curve in class
            if not (k[-1] == '1' and k[-2].isalpha()):
                classes[len(classes)-1].append(v)
            else:
                classes.append([v])
        return classes

    def isogeny_class(self, label):
        """
        Return the isogeny class of elliptic curves that are
        isogenous to the curve with given Cremona label.

        INPUT:

        - ``label`` -- string

        OUTPUT: list of EllipticCurve objects

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.isogeny_class('11a1')
            [Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 - x^2 over Rational Field]
            sage: c.isogeny_class('12001a1')   # optional - database_cremona_ellcurve
            [Elliptic Curve defined by y^2 + x*y = x^3 - 101*x + 382 over Rational Field]
        """
        conductor, iso, num = parse_cremona_label(label)
        q = self.__connection__.cursor().execute("SELECT curve FROM t_curve "
            + "WHERE class=?",(str(conductor)+iso,))
        return [self.elliptic_curve(c[0]) for c in q]

    def iter_optimal(self, conductors):
        """
        Return an iterator through all optimal curves in the database
        with given conductors.

        INPUT:

        - ``conductors`` -- list or generator of ints

        OUTPUT: generator that iterates over EllipticCurve objects

        EXAMPLES:

        We list optimal curves with conductor up to 20::

            sage: [e.cremona_label() for e in CremonaDatabase().iter_optimal([11..20])]
            ['11a1', '14a1', '15a1', '17a1', '19a1', '20a1']

        Note the unfortunate 990h3 special case::

            sage: [e.cremona_label() for e in CremonaDatabase().iter_optimal([990])]
            ['990a1', '990b1', '990c1', '990d1', '990e1', '990f1', '990g1', '990h3', '990i1', '990j1', '990k1', '990l1']
        """
        for N in conductors:
            if N == 990:
                for c in self.__connection__.cursor().execute('SELECT class '
                    + 'FROM t_class WHERE conductor=990'):
                    if c[0][-1] == 'h':
                        yield self.elliptic_curve(c[0]+'3')
                    else:
                        yield self.elliptic_curve(c[0]+'1')
                continue
            for c in self.__connection__.cursor().execute('SELECT curve '
                + 'FROM t_curve,t_class USING(class) WHERE curve=class||1 '
                + 'AND conductor=?',(int(N),)):
                yield self.elliptic_curve(c[0])

    def list(self, conductors):
        """
        Return a list of all curves with given conductors.

        INPUT:

        - ``conductors`` -- list or generator of ints

        OUTPUT: list of EllipticCurve objects

        EXAMPLES::

            sage: CremonaDatabase().list([37])
            [Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 3*x + 1 over Rational Field]
        """
        return list(self.iter(conductors))

    def list_optimal(self, conductors):
        """
        Return a list of all optimal curves with given conductors.

        INPUT:

        - ``conductors`` -- list or generator of ints list of EllipticCurve
          objects

        OUTPUT: list of EllipticCurve objects

        EXAMPLES::

            sage: CremonaDatabase().list_optimal([37])
            [Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field]
        """
        return list(self.iter_optimal(conductors))

    def largest_conductor(self):
        """
        The largest conductor for which the database is complete.

        OUTPUT: integer; largest conductor

        EXAMPLES::

            sage: c = CremonaDatabase('cremona mini')
            sage: c.largest_conductor()
            9999
        """
        if hasattr(self, '__largest_conductor__'):
            return self.__largest_conductor__
        q = self.__connection__.cursor().execute('SELECT conductor FROM '
            + 't_class ORDER BY conductor DESC LIMIT 1')
        self.__largest_conductor__ = next(q)[0]
        return self.__largest_conductor__

    def smallest_conductor(self):
        """
        The smallest conductor for which the database is complete: always 1.

        OUTPUT: integer; smallest conductor

        .. NOTE::

            This always returns the integer 1, since that is the
            smallest conductor for which the database is complete,
            although there are no elliptic curves of conductor 1.  The
            smallest conductor of a curve in the database is 11.

        EXAMPLES::

            sage: CremonaDatabase().smallest_conductor()
            1
        """
        return 1

    def conductor_range(self):
        """
        Return the range of conductors that are covered by the database.

        OUTPUT: tuple of ints (N1,N2+1) where N1 is the smallest and
        N2 the largest conductor for which the database is complete.

        EXAMPLES::

            sage: c = CremonaDatabase('cremona mini')
            sage: c.conductor_range()
            (1, 10000)
        """
        return 1, self.largest_conductor()+1

    def number_of_curves(self,  N=0, i=0):
        """
        Return the number of curves stored in the database with conductor
        `N`. If `N = 0`, returns the total number of curves in the database.

        If `i` is nonzero, returns the number of curves in the `i`-th isogeny
        class. If `i` is a Cremona letter code, e.g., 'a' or 'bc', it is
        converted to the corresponding number.

        INPUT:

        - ``N`` -- integer
        - ``i`` -- integer or string

        OUTPUT: integer

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.number_of_curves(11)
            3
            sage: c.number_of_curves(37)
            4
            sage: c.number_of_curves(990)
            42
            sage: num = c.number_of_curves()
        """
        if N == 0:
            if hasattr(self, '__number_of_curves__'):
                return self.__number_of_curves__
            q = self.__connection__.cursor().execute('SELECT COUNT(curve) '
                + 'FROM t_curve')
            self.__number_of_curves__ = next(q)[0]
            return self.__number_of_curves__
        if i == 0:
            q = self.__connection__.cursor().execute('SELECT COUNT(curve) '
                + 'FROM t_curve,t_class USING(class) WHERE conductor=?',
                (int(N),))
            return next(q)[0]
        if not isinstance(i, str):
            i = cremona_letter_code(i)
        q = self.__connection__.cursor().execute('SELECT COUNT(curve) FROM '
            + 't_curve WHERE class=?',(str(N)+i,))
        return next(q)[0]

    def number_of_isogeny_classes(self, N=0):
        """
        Return the number of isogeny classes of curves in the database of
        conductor N. If N is 0, return the total number of isogeny classes
        of curves in the database.

        INPUT:

        - ``N`` -- integer

        OUTPUT: integer

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.number_of_isogeny_classes(11)
            1
            sage: c.number_of_isogeny_classes(37)
            2
            sage: num = c.number_of_isogeny_classes()
        """
        if N == 0:
            if hasattr(self, '__number_of_isogeny_classes__'):
                return self.__number_of_isogeny_classes__
            q = self.__connection__.cursor().execute('SELECT COUNT(class) '
                + 'FROM t_class')
            self.__number_of_isogeny_classes__ = next(q)[0]
            return self.__number_of_isogeny_classes__
        q = self.__connection__.cursor().execute('SELECT COUNT(class) FROM '
            + 't_class WHERE conductor=?',(int(N),))
        return next(q)[0]

    def random(self):
        """
        Return a random curve from the database.

        EXAMPLES::

            sage: CremonaDatabase().random() # random -- depends on database installed
            Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 - 224*x + 3072 over Rational Field
        """
        N = randint(11, self.largest_conductor())
        q = self.__connection__.cursor().execute('SELECT conductor FROM '
            + 't_class WHERE conductor>=? ORDER BY conductor',(int(N),))
        try:
            N = next(q)[0]
        except StopIteration:
            N = 11
        iso = randint(0, self.number_of_isogeny_classes(N)-1)
        iso = cremona_letter_code(iso)
        num = randint(1, self.number_of_curves(N,iso))
        return self.elliptic_curve(str(N)+iso+str(num))

    ###############################################################################
    # Functions for loading data from Cremona's ftpdata directory.
    ###############################################################################
    def _init_from_ftpdata(self, ftpdata, largest_conductor=0):
        """Create the SQL Cremona Database from the Cremona elliptic curve
        data repository ecdata, which is available from
        https://github.com/JohnCremona/ecdata.

        To create the large database from Cremona's text files, see
        sage.databases.cremona.build.  Alternatively:

        If the Cremona database has already been installed, remove
        `SAGE_DATA/cremona/cremona.db`. Then run::

            sage: C = sage.databases.cremona.LargeCremonaDatabase('cremona',False, True)  # not tested
            sage: C._init_from_ftpdata('path/to/ecdata/',0)                               # not tested

        EXAMPLES::

            sage: d = sage.databases.cremona.MiniCremonaDatabase(name='cremona', read_only=False, rebuild=True)   # not tested
            sage: d._init_from_ftpdata('/home/jec/ecdata')                                                        # not tested
        """
        if self.__read_only__:
            raise RuntimeError("The database must not be read_only.")

        if not os.path.exists(ftpdata):
            raise RuntimeError("The cremona ftpdata directory '" + ftpdata
                + "' does not exist.")

        if largest_conductor:
            print("largest conductor =", largest_conductor)
            self.__largest_conductor__ = largest_conductor

        # Since July 2014 the data files have been arranged in
        # subdirectories (see trac #16903).
        allcurves_dir = os.path.join(ftpdata,'allcurves')
        allbsd_dir = os.path.join(ftpdata,'allbsd')
        allgens_dir = os.path.join(ftpdata,'allgens')
        degphi_dir = os.path.join(ftpdata,'degphi')
        num_curves, num_iso_classes = self._init_allcurves(allcurves_dir, largest_conductor)
        self.__number_of_curves__ = num_curves
        self.__number_of_isogeny_classes__ = num_iso_classes
        if hasattr(self, 'degphi'):
            self._init_degphi(degphi_dir, largest_conductor)
        if hasattr(self, 'allbsd'):
            self._init_allbsd(allbsd_dir, largest_conductor)
        if hasattr(self, 'allgens'):
            self._init_allgens(allgens_dir, largest_conductor)
        self.vacuum()

    def _init_allcurves(self, ftpdata, largest_conductor=0):
        """
        Initialize the allcurves table by reading the corresponding ftpdata
        files and importing them into the database.

        To create the large database from Cremona's text files, see
        sage.databases.cremona.build, do NOT run this method directly.

        INPUT:

        - ``ftpdata`` -- string; the name of the directory in which the data is

        - ``largest_conductor`` -- integer (default: 0); if 0,
          then only include data up to that conductor

        OUTPUT:

        -  integer; number_of_curves
        -  integer; number_of_isogeny_classes

        EXAMPLES::

            sage: d = sage.databases.cremona.MiniCremonaDatabase(name='cremona', read_only=False, rebuild=True)   # not tested
            sage: d._init_allcurves('.', 11)    # not tested
            (3, 1)
        """
        if self.__read_only__:
            raise RuntimeError("The database must not be read_only.")
        files = sorted(os.listdir(ftpdata))
        name = 'allcurves'
        num_curves = 0
        num_iso_classes = 0
        con = self.get_connection()
        for F in files:
            if not F[:len(name)] == name:
                continue
            print("Inserting", F)
            class_data = []
            curve_data = []
            for L in open(ftpdata + "/" + F).readlines():
                N, iso, num, ainvs, r, tor = L.split()
                if largest_conductor and int(N) > largest_conductor:
                    break
                cls = N+iso
                cur = cls+num
                if num == "1":
                    class_data.append((N,cls,r))
                    num_iso_classes += 1
                curve_data.append((cur,cls,ainvs,tor))
                num_curves += 1
            con.executemany('INSERT INTO t_class (conductor,class,rank) '
                + 'VALUES (?,?,?)', class_data)
            con.executemany('INSERT INTO t_curve (curve,class,eqn,tors) '
                + 'VALUES (?,?,?,?)', curve_data)
            print("Committing...")
            print("num_iso_classes =", num_iso_classes)
            self.commit()
            if largest_conductor and int(N) > largest_conductor:
                break
        return num_curves, num_iso_classes


class LargeCremonaDatabase(MiniCremonaDatabase):
    """
    The Cremona database of elliptic curves.

    EXAMPLES::

        sage: c = CremonaDatabase('cremona')  # optional - database_cremona_ellcurve
        sage: c.allcurves(11)                 # optional - database_cremona_ellcurve
        {'a1': [[0, -1, 1, -10, -20], 0, 5],
        'a2': [[0, -1, 1, -7820, -263580], 0, 1],
        'a3': [[0, -1, 1, 0, 0], 0, 5]}
    """
    _expected_skeleton = _cremonaSkeleton

    def allbsd(self, N):
        r"""
        Return the allbsd table for conductor N. The entries are::

            [id, tamagawa_product, Omega_E, L, Reg_E, Sha_an(E)]

        where id is the isogeny class (letter) followed by a number, e.g.,
        b3, and L is `L^r(E,1)/r!`, where E has rank r.

        INPUT:

        - ``N`` -- integer; the conductor

        OUTPUT: dictionary containing the allbsd table for each isogeny class
        in conductor N

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.allbsd(12)            # optional - database_cremona_ellcurve
            {}
            sage: c.allbsd(19)['a3']      # optional - database_cremona_ellcurve
            [1, 4.07927920046493, 0.453253244496104, 1.0, 1]
            sage: c.allbsd(12001)['a1']   # optional - database_cremona_ellcurve
            [2, 3.27608135248722, 1.54910143090506, 0.236425971187952, 1.0]
        """
        ret = {}
        for c in self.__connection__.cursor().execute('SELECT curve,cp,om,L,'
            + 'reg,sha FROM t_curve,t_class USING(class) WHERE conductor=?',
            (int(N),)):
            N,iso,num = parse_cremona_label(c[0])
            ret[iso+str(num)] = list(c[1:])
        return ret

    def allgens(self, N):
        """
        Return the allgens table for conductor N.

        INPUT:

        - ``N`` -- integer; the conductor

        OUTPUT: dictionary; id:[points, ...], ...

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.allgens(12)            # optional - database_cremona_ellcurve
            {}
            sage: c.allgens(1001)['a1']    # optional - database_cremona_ellcurve
            [[61, 181, 1]]
            sage: c.allgens(12001)['a1']   # optional - database_cremona_ellcurve
            [[7, 2, 1]]
        """
        ret = {}
        for c in self.__connection__.cursor().execute('SELECT curve,gens '
            + 'FROM t_curve,t_class USING(class) WHERE conductor=?',(int(N),)):
            N,iso,num = parse_cremona_label(c[0])
            ret[iso+str(num)] = eval(c[1])
        return ret

    def degphi(self, N):
        """
        Return the degphi table for conductor N.

        INPUT:

        - ``N`` -- integer; the conductor

        OUTPUT: dictionary; id:degphi, ...

        EXAMPLES::

            sage: c = CremonaDatabase()
            sage: c.degphi(11)            # optional - database_cremona_ellcurve
            {'a1': 1}
            sage: c.degphi(12001)['c1']   # optional - database_cremona_ellcurve
            1640
        """
        ret = {}
        for c in self.__connection__.cursor().execute('SELECT curve,deg FROM'
            + ' t_curve,t_class USING(class) WHERE curve=class||1 AND '
            + 'conductor=?', (int(N),)):
            N,iso,num = parse_cremona_label(c[0])
            ret[iso+str(num)] = c[1]
        return ret

    def _init_degphi(self, ftpdata, largest_conductor=0):
        """
        Initialize the degphi table by reading the corresponding ftpdata
        files and importing them into the database.

        To create the large database from Cremona's text files, see
        sage.databases.cremona.build, do NOT run this method directly.

        EXAMPLES::

            sage: d = sage.databases.cremona.LargeCremonaDatabase(name='cremona', read_only=False, rebuild=True)   # not tested
            sage: d._init_degphi('.')           # not tested
        """
        if self.__read_only__:
            raise RuntimeError("The database must not be read_only.")
        files = sorted(os.listdir(ftpdata))
        name = "degphi"
        con = self.get_connection()
        for F in files:
            if not F[:len(name)] == name:
                continue
            print("Inserting", F)
            class_data = []
            for L in open(ftpdata + "/" + F).readlines():
                N, iso, num, degree, primes, curve = L.split()
                if largest_conductor and int(N) > largest_conductor:
                    break
                class_data.append((degree,N+iso))
            con.executemany('UPDATE t_class SET deg=? WHERE class=?',
                class_data)
            print("Committing...")
            self.commit()
            if largest_conductor and int(N) > largest_conductor:
                break

    def _init_allbsd(self, ftpdata, largest_conductor=0):
        """
        Initialize the allbsd table by reading the corresponding ftpdata
        files and importing them into the database.

        To create the large database from Cremona's text files, see
        sage.databases.cremona.build, do NOT run this method directly.

        EXAMPLES::

            sage: d = sage.databases.cremona.LargeCremonaDatabase(name='cremona', read_only=False, rebuild=True)   # not tested
            sage: d._init_allbsd('.')           # not tested
        """
        if self.__read_only__:
            raise RuntimeError("The database must not be read_only.")
        files = sorted(os.listdir(ftpdata))
        name = "allbsd"
        con = self.get_connection()
        for F in files:
            if not F[:len(name)] == name:
                continue
            print("Inserting", F)
            curve_data = []
            class_data = []
            for L in open(ftpdata + "/" + F).readlines():
                N, iso, num, eqn, rank, tor, cp, om, L, reg, sha = L.split()
                if largest_conductor and int(N) > largest_conductor:
                    break
                cls = N+iso
                if num == "1":
                    class_data.append((L,cls))
                curve_data.append((cp,om,reg,eval(sha),cls+num))
            con.executemany("UPDATE t_class SET L=? WHERE class=?", class_data)
            con.executemany("UPDATE t_curve SET cp=?,om=?,reg=?,sha=? WHERE "
                    + "curve=?", curve_data)
            print("Committing...")
            self.commit()
            if largest_conductor and int(N) > largest_conductor:
                break

    def _init_allgens(self, ftpdata, largest_conductor=0):
        """
        Initialize the allgens table by reading the corresponding ftpdata
        files and importing them into the database.

        To create the large database from Cremona's text files, see
        sage.databases.cremona.build, do NOT run this method directly.

        EXAMPLES::

            sage: d = sage.databases.cremona.LargeCremonaDatabase(name='cremona', read_only=False, rebuild=True)   # not tested
            sage: d._init_allgens('.')          # not tested
        """
        if self.__read_only__:
            raise RuntimeError("The database must not be read_only.")
        files = sorted(os.listdir(ftpdata))
        name = "allgens"
        con = self.get_connection()
        for F in files:
            if not F[:len(name)] == name:
                continue
            print("Inserting", F)
            curve_data = []
            for L in open(ftpdata + "/" + F).readlines():
                v = L.split()
                if largest_conductor and int(v[0]) > largest_conductor:
                    break
                gens = '['+','.join(v[6:6+int(v[4])]).replace(':',',')+']'
                curve_data.append((gens,''.join(v[:3])))
            con.executemany("UPDATE t_curve SET gens=? WHERE curve=?",
                curve_data)
            print("Committing...")
            self.commit()
            if largest_conductor and int(v[0]) > largest_conductor:
                break


_db = None


def CremonaDatabase(name=None, mini=None, set_global=None):
    """
    Initialize the Cremona database with name ``name``. If ``name`` is
    ``None`` it instead initializes large Cremona database (named 'cremona'),
    if available or default mini Cremona database (named 'cremona mini').

    If the Cremona database in question is in the format of the mini database,
    you must set ``mini=True``, otherwise it must be set to ``False``.

    TESTS::

        sage: c = CremonaDatabase()
        sage: isinstance(c, sage.databases.cremona.MiniCremonaDatabase)
        True
        sage: isinstance(c, sage.databases.cremona.LargeCremonaDatabase)  # optional - database_cremona_ellcurve
        True

    Verify that :issue:`12341` has been resolved::

        sage: c = CremonaDatabase('should not exist', mini=True)
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: database_should_not_exist_ellcurve is not available.
        '...db' not found in any of [...]
        ...Further installation instructions might be available at https://github.com/JohnCremona/ecdata.
        sage: c = CremonaDatabase('should not exist',mini=False)
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: database_should_not_exist_ellcurve is not available.
        '...db' not found in any of [...]
        ...Further installation instructions might be available at https://github.com/JohnCremona/ecdata.

    Verify that :issue:`39072` has been resolved::

        sage: C = CremonaDatabase(mini=False)  # optional - !database_cremona_ellcurve
        Traceback (most recent call last):
        ...
        ValueError: full Cremona database is not available; consider using mini Cremona database by mini=True
    """
    if set_global is not None:
        from sage.misc.superseded import deprecation
        deprecation(25825, "the set_global argument for CremonaDatabase is deprecated and ignored")

    if name is None:
        if mini is None:
            if DatabaseCremona().is_present():
                name = 'cremona'
                mini = False
            else:
                name = 'cremona mini'
                mini = True
        elif mini:
            name = 'cremona mini'
        else:
            if not DatabaseCremona().is_present():
                raise ValueError('full Cremona database is not available; '
                                 'consider using mini Cremona database by mini=True')
            name = 'cremona'
    elif name == 'cremona mini':
        mini = True
    elif name == 'cremona':
        mini = False
    else:
        if mini is None:
            raise ValueError('mini must be set as either True or False')

    if mini:
        return MiniCremonaDatabase(name)
    else:
        return LargeCremonaDatabase(name)
