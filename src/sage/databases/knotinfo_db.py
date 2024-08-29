r"""
KnotInfo database

This module contains the class :class:`KnotInfoDataBase` and auxiliary classes
for it, which serves as an interface to the lists of named knots and links provided
at the web-pages `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and
`LinkInfo <https://linkinfo.sitehost.iu.edu>`__.

To use the database, you need to install the optional :ref:`database_knotinfo
<spkg_database_knotinfo>` package by the Sage command ::

    sage -i database_knotinfo

EXAMPLES::

    sage: # optional - database_knotinfo
    sage: from sage.databases.knotinfo_db import KnotInfoDataBase
    sage: ki_db = KnotInfoDataBase()
    sage: ki_db
    <sage.databases.knotinfo_db.KnotInfoDataBase object at ...>

AUTHORS:

- Sebastian Oehms (2020-08): initial version
"""
##############################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
##############################################################################
import os
from enum import Enum

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.persist import save, load
from sage.misc.verbose import verbose
from sage.misc.cachefunc import cached_method

columns_white_list = ['knot_atlas_anon', 'knotilus_page_anon']
columns_black_list = ['homfly_polynomial_old']

class KnotInfoColumnTypes(Enum):
    r"""
    Enum class to specify if a column from the table of knots and links provided
    at the web-pages `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and
    `LinkInfo <https://linkinfo.sitehost.iu.edu>`__.  is used for knots only,
    links only or both.

    EXAMPLES::

        sage: from sage.databases.knotinfo_db import KnotInfoColumnTypes
        sage: [col_type for col_type in KnotInfoColumnTypes]
        [<KnotInfoColumnTypes.OnlyKnots: 'K'>,
        <KnotInfoColumnTypes.OnlyLinks: 'L'>,
        <KnotInfoColumnTypes.KnotsAndLinks: 'B'>]
    """
    OnlyKnots = 'K'       # column that is only used in the KnotInfo table
    OnlyLinks = 'L'       # column that is only used in the LinkInfo table
    KnotsAndLinks = 'B'       # column that is only used in both tables


class KnotInfoColumns(Enum):
    r"""
    Enum class to select a column from the table of knots and links provided
    at the web-pages `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and
    `LinkInfo <https://linkinfo.sitehost.iu.edu>`__.

    EXAMPLES::

        sage: from sage.databases.knotinfo_db import KnotInfoDataBase
        sage: ki_db = KnotInfoDataBase()
        sage: cols = ki_db.columns(); cols
        <enum 'Columns'>
        sage: from sage.databases.knotinfo_db import KnotInfoColumns
        sage: isinstance(cols.name, KnotInfoColumns)
        True

        sage: def only_links(c):
        ....:     return c.column_type() == c.types.OnlyLinks
        sage: [c.column_name() for c in cols if only_links(c)]  # optional - database_knotinfo
        ['Name - Unoriented',
         'Orientation',
         'Unoriented Rank',
         'PD Notation (vector)',
         'PD Notation (KnotTheory)',
         'Braid Notation',
         'Quasipositive Braid',
         'Multivariable Alexander Polynomial',
         'HOMFLYPT Polynomial',
         'Khovanov Polynomial',
         'Unoriented',
         'Arc Notation',
         'Linking Matrix',
         'Rolfsen Name',
         'Components',
         'DT code',
         'Splitting Number',
         'Nullity',
         'Unlinking Number',
         'Weak Splitting Number']
    """
    @property
    def types(self):
        r"""
        Return :class:`KnotInfoColumnTypes` to be used for checks.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: cols = ki_db.columns()
            sage: cols.dt_code.column_type() == cols.dt_code.types.OnlyLinks
            True
        """
        return KnotInfoColumnTypes

    def column_name(self):
        r"""
        Return the name of ``self`` displayed on the KnotInfo web-page.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: cols = ki_db.columns()
            sage: cols.dt_code.column_name()
            'DT code'
        """
        return self.value[0]

    def column_type(self):
        r"""
        Return the type of ``self``. That is an instance of ``Enum``
        :class:`KnotInfoColumnTypes`.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: cols = ki_db.columns()
            sage: cols.homfly_polynomial.column_type()
            <KnotInfoColumnTypes.OnlyKnots: 'K'>
            sage: cols.homflypt_polynomial.column_type()
            <KnotInfoColumnTypes.OnlyLinks: 'L'>
            sage: cols.name.column_type()
            <KnotInfoColumnTypes.KnotsAndLinks: 'B'>
        """
        return self.value[1]

    def description_webpage(self, new=0, autoraise=True):
        r"""
        Launch the description page of ``self`` in the standard web browser.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: cols = ki_db.columns()
            sage: cols.pd_notation.description_webpage()            # not tested
            True
            sage: cols.homflypt_polynomial.description_webpage()    # not tested
            True
        """
        import webbrowser
        if self.column_type() == self.types.OnlyLinks:
            url = KnotInfoFilename.links.description_url(self)
        else:
            url = KnotInfoFilename.knots.description_url(self)
        return webbrowser.open(url, new=new, autoraise=autoraise)


class KnotInfoFilename(Enum):
    r"""
    Enum for the different data files. The following choices are possible:

    - ``knots`` -- contains the data from KnotInfo
    - ``links`` -- contains the data for proper links from LinkInfo

    Examples::

        sage: from sage.databases.knotinfo_db import KnotInfoDataBase
        sage: ki_db = KnotInfoDataBase()
        sage: ki_db.filename
        <enum 'KnotInfoFilename'>
    """

    def url(self):
        r"""
        Return the URL to download the data from the web-page.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.url()
            'https://knotinfo.math.indiana.edu/'
        """
        if self == KnotInfoFilename.knots:
            return self.value[0]
        else:
            return self.value[0]

    def excel(self):
        r"""
        Return the Excel-file name to download the data from the web-page.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.excel()
            'knotinfo_data_complete.xls'
        """
        if self == KnotInfoFilename.knots:
            return '%s.xls' % (self.value[1])
        else:
            return '%s.xlsx' % (self.value[1])

    def csv(self):
        r"""
        Return the file name under which the data from the web-page
        are stored as csv file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.csv()
            'knotinfo_data_complete.csv'
        """
        return '%s.csv' % (self.value[1])

    def num_knots(self, version):
        r"""
        Return the file name under which the number of knots is stored
        in an sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.num_knots('21.7')
            'num_knots_21.7.sobj'
        """
        return 'num_knots_%s.sobj' % version

    def sobj_row(self):
        r"""
        Return the file name under which the row-data of the csv-File
        is stored as python dictionary in a sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.sobj_row()
            'row_dict.sobj'
        """
        return 'row_dict.sobj'

    def sobj_column(self):
        r"""
        Return the file name under which the column-data of the csv-File
        is stored as python dictionary in a sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.sobj_column()
            'column_dict.sobj'
        """
        return 'column_dict.sobj'

    def sobj_data(self, column):
        r"""
        Return the file name under which the data of the given
        column is stored as python list in a sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.sobj_data(ki_db.columns().braid_notation)
            'knotinfo_braid_notation'
        """
        if column.column_type() == column.types.OnlyLinks:
            return 'linkinfo_%s' % (column.name)
        else:
            return 'knotinfo_%s' % (column.name)

    def description_url(self, column):
        r"""
        Return the url of the description page of the given column.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.description_url(ki_db.columns().braid_notation)
            'https://knotinfo.math.indiana.edu/descriptions/braid_notation.html'
        """
        return '%sdescriptions/%s.html' % (self.url(), column.name)

    def diagram_url(self, fname, single=False):
        r"""
        Return the url of the diagram page of the given link.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.diagram_url('3_1-50.png')
            'https://knotinfo.math.indiana.edu/diagram_display.php?3_1-50.png'
            sage: ki_db.filename.knots.diagram_url('3_1', single=True)
            'https://knotinfo.math.indiana.edu/diagrams/3_1'
        """
        if single:
            return '%sdiagrams/%s' % (self.url(), fname)
        else:
            return '%sdiagram_display.php?%s' % (self.url(), fname)

    knots = ['https://knotinfo.math.indiana.edu/', 'knotinfo_data_complete']
    links = ['https://linkinfo.sitehost.iu.edu/',  'linkinfo_data_complete']


#----------------------------------------------------------------------------------------------------------------------------
# Class to provide data for knots and links from the KnotInfo web-page
#----------------------------------------------------------------------------------------------------------------------------
class KnotInfoDataBase(SageObject, UniqueRepresentation):
    r"""
    Database interface to KnotInfo.

    The original data are obtained from KnotInfo web-page (URL see the example
    below). In order to have these data installed during the build process as
    a sage-package they are converted as csv files into a tarball. This tarball
    has been created using the method :meth:`create_spkg_tarball`.

    EXAMPLES::

        sage: from sage.databases.knotinfo_db import KnotInfoDataBase
        sage: ki_db = KnotInfoDataBase()
        sage: ki_db.filename.knots
        <KnotInfoFilename.knots: ['https://knotinfo.math.indiana.edu/',
                                  'knotinfo_data_complete']>
    """

    filename = KnotInfoFilename

    def __init__(self, install=False):
        r"""
        Python constructor.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.links
            <KnotInfoFilename.links: ['https://linkinfo.sitehost.iu.edu/',
                                      'linkinfo_data_complete']>
        """
        # some constants
        self._names_column = 'name'
        self._components_column = 'components'
        self._knot_prefix = 'K'

        self._knot_list = None
        self._link_list = None
        self._demo = None
        self._num_knots = None

        from sage.features.databases import DatabaseKnotInfo
        from sage.env import DOT_SAGE
        self._feature = DatabaseKnotInfo()
        self._sobj_path = os.path.join(DOT_SAGE, 'knotinfo')

    def create_filecache(self, force=False):
        r"""
        Create the internal files containing the database.

        INPUT:

        - ``force`` -- boolean (default: ``False``); if set to ``True`` the
          existing file-cache is overwritten

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.create_filecache()    # optional - database_knotinfo
        """
        if not self._feature.is_present():
            return

        if os.path.isdir(self._sobj_path):
            # if it exists then remove it if it belongs to an older version of
            # the database or should be reset by the user because it is damaged.
            test_version = os.path.join(self._sobj_path, self.filename.knots.num_knots(self.version()))
            if force or not os.path.isfile(test_version):
                import shutil
                shutil.rmtree(self._sobj_path)

        from sage.misc.temporary_file import atomic_dir
        with atomic_dir(self._sobj_path) as d:
            sobj_path = d.name
            num_knots_file = os.path.join(sobj_path, self.filename.knots.num_knots(self.version()))
            knot_list = self.knot_list()
            num_knots = len(knot_list) - 1
            save(num_knots, num_knots_file)
            self._num_knots = num_knots
            self._create_col_dict_sobj(sobj_path=sobj_path)
            self._create_data_sobj(sobj_path=sobj_path)
        return

    def version(self):
        r"""
        Return the version of the database currently installed on the device.

        .. NOTE::

            The development of the original databases on the KnotInfo and
            LinkInfo web-pages is in a continuous flow. The installed version
            can be behind the current available state of these databases. Every
            month a cronjob on the
            `GitHub repository <https://github.com/soehms/database_knotinfo/>`__
            searches for differences and creates a new release on
            `PyPI <https://pypi.org/project/database-knotinfo/>`__ in case of
            success.

            If you note that your version is behind the version on PyPI
            and would like to have Sage working with that release you should
            first try to upgrade using ``sage -i database_knotinfo``. If this
            is not successful even though you are on the latest Sage release
            please create an issue for that in the GitHub repository.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.version()   >= '2021.10.1'   # optional database_knotinfo
            True
        """
        self._feature.require()
        from database_knotinfo import version
        return version()

    def demo_version(self):
        r"""
        Return whether the KnotInfo databases are installed completely or
        just the demo version is used.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.demo_version()       # optional - database_knotinfo
            False
        """
        if self._demo is None:
            if self._feature.is_present():
                num_knots_file = os.path.join(self._sobj_path, self.filename.knots.num_knots(self.version()))
                from builtins import FileNotFoundError
                try:
                    self._num_knots = load(num_knots_file)
                except FileNotFoundError:
                    self.create_filecache()
                self._demo = False
            else:
                self._demo = True
                self._num_knots = len([v for v in row_demo_sample.values() if v[1] == 1])
        return self._demo

    def knot_list(self):
        r"""
        Return the KnotInfo table rows as Python list.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: len(ki_db.knot_list())  # not tested (just used on installation)
        """
        if self._knot_list:
            return self._knot_list

        from database_knotinfo import link_list
        self._knot_list = link_list()
        return self._knot_list

    def link_list(self):
        r"""
        Return the LinkInfo table rows as Python list.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: len(ki_db.link_list())  # not tested (just used on installation)
        """
        if self._link_list:
            return self._link_list

        from database_knotinfo import link_list
        self._link_list = link_list(proper_links=True)
        return self._link_list

    def _create_col_dict_sobj(self, sobj_path=None):
        r"""
        Create ``sobj`` files containing the number of knots and a dictionary
        for the columns of the table.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db._create_col_dict_sobj() # not tested (used on installation)
        """
        knot_list = self.knot_list()
        knot_column_names = knot_list[0]

        link_list = self.link_list()
        link_column_names = link_list[0]

        column_dict = {}

        if not sobj_path:
            sobj_path = self._sobj_path

        # ----------------------------------------------------------------
        # Columns that exist for knots and links
        # ----------------------------------------------------------------
        for col in knot_column_names:

            name = knot_column_names[col]
            if not name and col not in columns_white_list:
                # not of interest
                continue

            if col in columns_black_list:
                # not of interest
                continue

            col_type = KnotInfoColumnTypes.OnlyKnots
            if col in link_column_names:
                col_type = KnotInfoColumnTypes.KnotsAndLinks
            column_dict[col] = [name, col_type]

        # ----------------------------------------------------------------
        # Columns that exist for links only
        # ----------------------------------------------------------------
        for col in link_column_names:

            name = link_column_names[col]
            if not name and col not in columns_white_list:
                # not of interest
                continue

            if col in columns_black_list:
                # not of interest
                continue

            if col in knot_column_names:
                # already used
                continue

            col_type = KnotInfoColumnTypes.OnlyLinks
            column_dict[col] = [name, col_type]

        save(column_dict, '%s/%s' % (sobj_path, self.filename.knots.sobj_column()))

    def _create_data_sobj(self, sobj_path=None):
        r"""
        Create ``sobj`` files containing the contents of the whole table.
        To each column there is created one file containing a list of
        strings giving the entries of the database table.

        The length of these lists depends on the type of the corresponding
        column. If a column is used in both tables
        (``KnotInfoColumnTypes.KnotsAndLinks``) the list of proper links
        is appended to the list of knots.  In both other cases the lenght
        of the list corresponds to the number of listed knots and proper
        links respectively.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db._create_data_sobj() # not tested (just used on installation)
        """
        knot_list = self.knot_list()
        link_list = self.link_list()
        len_knots = len(knot_list)
        len_links = len(link_list)

        row_dict = {}

        if not sobj_path:
            sobj_path = self._sobj_path

        # ----------------------------------------------------------------
        # Columns that exist for knots and links
        # ----------------------------------------------------------------
        column_dict = load('%s/%s' % (sobj_path, self.filename.knots.sobj_column()))
        cols = KnotInfoColumns('ColsTemp', column_dict)
        for col in cols:
            val_list = []

            if col.column_type() != col.types.OnlyLinks:
                for i in range(1, len_knots):
                    if col.name == self._names_column:
                        row_dict[self._knot_prefix + knot_list[i][col.name]] = [i - 1, 1]

                    val_list.append(knot_list[i][col.name])

            if col.column_type() != col.types.OnlyKnots:
                for i in range(1, len_links):
                    if col.name == self._names_column:
                        link_name = link_list[i][col.name]
                        link_name = link_name.replace('{', '_')
                        link_name = link_name.replace(',', '_')
                        link_name = link_name.replace('}', '')

                        num_comp = int(link_list[i][self._components_column])
                        row_dict[link_name] = [i + len_knots - 2, num_comp]

                    val_list.append(link_list[i][col.name])

            if val_list:
                save(val_list, '%s/%s' % (sobj_path, self.filename.knots.sobj_data(col)))

        save(row_dict,    '%s/%s' % (sobj_path, self.filename.knots.sobj_row()))

    @cached_method
    def columns(self):
        r"""
        Return the columns ot the database table as instances of enum class
        :class:`KnotInfoColumns`.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: cols = ki_db.columns()
            sage: [col.column_name() for col in cols if col.name.startswith('pd')]   # optional - database_knotinfo
            ['PD Notation', 'PD Notation (vector)', 'PD Notation (KnotTheory)']
        """
        column_dict = self.read_column_dict()
        return KnotInfoColumns('Columns', column_dict)

    # -------------------------------------------------------------------------------------------------------------
    # read the dictionary for the column names from sobj-file
    # -------------------------------------------------------------------------------------------------------------
    @cached_method
    def read_column_dict(self):
        r"""
        Read the dictionary for the column names from the according sobj-file.

        OUTPUT: a Python dictionary containing the column names and types

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: len(ki_db.read_column_dict()) > 120       # optional - database_knotinfo
            True
        """
        if self.demo_version():
            return column_demo_sample
        sobj_path = self._sobj_path
        filename = self.filename.knots.sobj_column()
        return load('%s/%s' % (sobj_path, filename))

    # -------------------------------------------------------------------------------------------------------------
    # read the dictionary for the row names that is the knot and link names from sobj-file
    # -------------------------------------------------------------------------------------------------------------
    @cached_method
    def read_row_dict(self):
        r"""
        Read the dictionary for the row names that is the knot and link names
        from the according sobj-file

        OUTPUT:

        A python dictionary containing the names of the knots and links
        together with their table index and the corresponding number of
        components

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.read_row_dict()['K7_1']
            [8, 1]
        """
        if self.demo_version():
            return row_demo_sample
        sobj_path = self._sobj_path
        filename = self.filename.knots.sobj_row()
        return load('%s/%s' % (sobj_path, filename))

    # -------------------------------------------------------------------------------------------------------------
    # return a dictionary to obtain the original name to a row_dict key
    # -------------------------------------------------------------------------------------------------------------
    @cached_method
    def row_names(self):
        r"""
        Return a dictionary to obtain the original name to a row_dict key.

        OUTPUT: a Python dictionary containing the names of the knots and links
        together with their original names from the database

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.row_names()['K7_1']      # optional - database_knotinfo
            '7_1'
        """
        row_dict = self.read_row_dict()
        names = self.read(self.columns().name)
        return {k:names[v[0]] for k, v in row_dict.items()}

    # -------------------------------------------------------------------------------------------------------------
    # read the number of knots contained in the database (without proper links) from the according sobj-file.
    # -------------------------------------------------------------------------------------------------------------
    def read_num_knots(self):
        r"""
        Read the number of knots contained in the database (without
        proper links) from the according sobj-file.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.read_num_knots()              # optional - database_knotinfo
            12966
        """
        if not self._num_knots:
            self.demo_version()
        return self._num_knots

    # -------------------------------------------------------------------------------------------------------------
    # read an sobj-file obtained from KnotInfo database
    # -------------------------------------------------------------------------------------------------------------
    @cached_method
    def read(self, column):
        r"""
        Access a column of KnotInfo / LinkInfo.

        INPUT:

        - ``column`` -- instance of enum :class:`KnotInfoColumns`
          to select the data to be read in

        OUTPUT: a Python list containing the data corresponding to the column

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
        """
        if not isinstance(column, KnotInfoColumns):
            raise TypeError('column must be an instance of enum %s' % (KnotInfoColumns))

        if self.demo_version():
            return data_demo_sample[column]

        sobj_path = self._sobj_path
        if column.column_type() == column.types.OnlyLinks:
            filename = self.filename.links.sobj_data(column)
        else:
            filename = self.filename.knots.sobj_data(column)

        verbose('loading data library %s ...' % (filename))
        res = load('%s/%s' % (sobj_path, filename))
        verbose('... finished!')

        return res

    def _test_database(self, **options):
        r"""
        Method used by TestSuite. Performs :meth:`KnotInfoBase.is_recoverable`.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: TestSuite(ki_db).run()    # long time indirect doctest
        """
        from sage.knots.knotinfo import KnotInfo
        from sage.misc.misc import some_tuples
        tester = options['tester']
        max_samples = tester._max_samples
        if not max_samples:
            max_samples = 20
        l = list(KnotInfo)
        sample = some_tuples(l, 1, len(l), max_samples=max_samples)
        tester.assertTrue(all(L.is_recoverable(unique=False) for L, in sample))


column_demo_sample = {
    'name':                 ['Name',                 KnotInfoColumnTypes.KnotsAndLinks],
    'name_unoriented':      ['Name - Unoriented',    KnotInfoColumnTypes.OnlyLinks],
    'dt_notation':          ['DT Notation',          KnotInfoColumnTypes.OnlyKnots],
    'gauss_notation':       ['Gauss Notation',       KnotInfoColumnTypes.KnotsAndLinks],
    'pd_notation':          ['PD Notation',          KnotInfoColumnTypes.OnlyKnots],
    'pd_notation_vector':   ['PD Notation (vector)', KnotInfoColumnTypes.OnlyLinks],
    'crossing_number':      ['Crossing Number',      KnotInfoColumnTypes.KnotsAndLinks],
    'braid_index':          ['Braid Index',          KnotInfoColumnTypes.OnlyKnots],
    'braid_length':         ['Braid Length',         KnotInfoColumnTypes.OnlyKnots],
    'braid_notation':       ['Braid Notation',       KnotInfoColumnTypes.KnotsAndLinks],
    'braid_notation_old':   ['Braid Notation',       KnotInfoColumnTypes.OnlyLinks],
    'alternating':          ['Alternating',          KnotInfoColumnTypes.KnotsAndLinks],
    'alexander_polynomial': ['Alexander',            KnotInfoColumnTypes.OnlyKnots],
    'jones_polynomial':     ['Jones',                KnotInfoColumnTypes.KnotsAndLinks],
    'conway_polynomial':    ['Conway',               KnotInfoColumnTypes.KnotsAndLinks],
    'homfly_polynomial':    ['HOMFLY',               KnotInfoColumnTypes.OnlyKnots],
    'homflypt_polynomial':  ['HOMFLYPT Polynomial',  KnotInfoColumnTypes.OnlyLinks],
    'kauffman_polynomial':  ['Kauffman',             KnotInfoColumnTypes.KnotsAndLinks],
    'khovanov_polynomial':  ['Khovanov',             KnotInfoColumnTypes.OnlyLinks],
    'khovanov_unreduced_integral_polynomial': ['KH Unred Z Poly', KnotInfoColumnTypes.OnlyKnots],
    'khovanov_reduced_integral_polynomial': ['KH Red Z Poly', KnotInfoColumnTypes.OnlyKnots],
    'khovanov_reduced_rational_polynomial': ['KH Red Q Poly', KnotInfoColumnTypes.OnlyKnots],
    'khovanov_reduced_mod2_polynomial': ['KH Red Mod2 Poly', KnotInfoColumnTypes.OnlyKnots],
    'khovanov_odd_integral_polynomial': ['KH Odd Red Z Poly', KnotInfoColumnTypes.OnlyKnots],
    'khovanov_odd_rational_polynomial': ['KH Odd Red Q Poly', KnotInfoColumnTypes.OnlyKnots],
    'khovanov_odd_mod2_polynomial': ['KH Red Odd Mod2 Poly', KnotInfoColumnTypes.OnlyKnots],
    'determinant':          ['Determinant',          KnotInfoColumnTypes.KnotsAndLinks],
    'positive':             ['Positive',             KnotInfoColumnTypes.OnlyKnots],
    'fibered':              ['Fibered',              KnotInfoColumnTypes.OnlyKnots],
    'unoriented':           ['Unoriented',           KnotInfoColumnTypes.OnlyLinks],
    'symmetry_type':        ['Symmetry Type',        KnotInfoColumnTypes.OnlyKnots],
    'width':                ['Width',                KnotInfoColumnTypes.OnlyKnots],
    'arc_notation':         ['Arc Notation',         KnotInfoColumnTypes.OnlyLinks],
    'dt_code':              ['DT code',              KnotInfoColumnTypes.OnlyLinks]
}


row_demo_sample = {
    'K0_1':   [0, 1],
    'K3_1':   [1, 1],
    'K4_1':   [2, 1],
    'K5_1':   [3, 1],
    'K5_2':   [4, 1],
    'K6_1':   [5, 1],
    'K6_2':   [6, 1],
    'K6_3':   [7, 1],
    'K7_1':   [8, 1],
    'K7_2':   [9, 1],
    'L2a1_0': [10, 2],
    'L2a1_1': [11, 2],
    'L4a1_0': [12, 2],
    'L4a1_1': [13, 2],
    'L5a1_0': [14, 2],
    'L5a1_1': [15, 2],
    'L6a1_0': [16, 2],
    'L6a1_1': [17, 2],
    'L6a2_0': [18, 2],
    'L6a2_1': [19, 2]
}

db = KnotInfoDataBase()
dc = db.columns()


data_demo_sample = {
    dc.name: ['0_1', '3_1', '4_1', '5_1', '5_2', '6_1', '6_2', '6_3', '7_1', '7_2',
              'L2a1{0}', 'L2a1{1}', 'L4a1{0}', 'L4a1{1}', 'L5a1{0}', 'L5a1{1}',
              'L6a1{0}', 'L6a1{1}', 'L6a2{0}', 'L6a2{1}', 'L6a3{0}'
             ],
    dc.name_unoriented: ['L2a1', 'L2a1', 'L4a1', 'L4a1', 'L5a1', 'L5a1', 'L6a1', 'L6a1', 'L6a2', 'L6a2', 'L6a3'],
    dc.crossing_number: ['0', '3', '4', '5', '5', '6', '6', '6', '7', '7', '2', '2', '4', '4', '5', '5', '6', '6', '6', '6', '6'],
    dc.braid_notation: [
        '',
        '{1,1,1}',
        '{1,-2,1,-2}',
        '{1,1,1,1,1}',
        '{1,1,1,2,-1,2}',
        '{1,1,2,-1,-3,2,-3}',
        '{1,1,1,-2,1,-2}',
        '{1,1,-2,1,-2,-2}',
        '{1,1,1,1,1,1,1}',
        '{1,1,1,2,-1,2,3,-2,3}',
        '{2, {-1, -1}}',
        '{2, {1, 1}}',
        '{3, {-2, -2, -1, 2, -1}}',
        '{2, {1, 1, 1, 1}}',
        '{3, {-1, 2, -1, 2, -1}}',
        '{3, {-1, 2, -1, 2, -1}}',
        '{4, {1, -2, 3, -2, 1, -2, -3, -2}}',
        '{3, {2, 2, 2, 1, 1, -2, 1}}',
        '{3, {-1, 2, -1, -2, -2, -1, -1}}',
        '{3, {1, -2, 1, 2, 2, 1, 1}}',
        '{2, {-1, -1, -1, -1, -1, -1}}'
        ],
    dc.braid_notation_old: [
        '{2, {-1, -1}}',
        '{2, {1, 1}}',
        '{4, {1, -2, 3, -2, -1, -2, -3, -2}}',
        '{2, {1, 1, 1, 1}}',
        '{3, {-1, 2, -1, 2, -1}}',
        '{3, {-1, 2, -1, 2, -1}}',
        '{4, {1, -2, 3, -2, 1, -2, -3, -2}}',
        '{4, {1, 2, 3, 2, 2, -1, 2, 2, -3, 2}}',
        '{4, {1, -2, -2, -2, 3, -2, -1, -2, -3, -2}}',
        '{4, {1, 2, -3, 2, -1, 2, 3, 2, 2, 2}}',
        '{2, {-1, -1, -1, -1, -1, -1}}'
        ],
    dc.braid_index: ['1', '2', '3', '2', '3', '4', '3', '3', '2', '4'],
    dc.braid_length: ['', '3', '4', '5', '6', '7', '6', '6', '7', '9'],
    dc.determinant: ['0', '3', '5', '5', '7', '9', '11', '13', '7', '11', '2', '2', '4', '4', '8', '8', '12', '12', '10', '10', '6'],
    dc.positive: ['', 'Y', 'N', 'Y', 'Y', 'N', 'N', 'N', 'Y', 'Y'],
    dc.fibered: ['', 'Y', 'Y', 'Y', 'N', 'N', 'Y', 'Y', 'Y', 'N'],
    dc.unoriented: ['Y', 'N', 'Y', 'N', 'Y', 'N', 'Y', 'N', 'Y', 'N', 'Y'],
    dc.pd_notation: [
        '',
        '[[1,5,2,4],[3,1,4,6],[5,3,6,2]]',
        '[[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]',
        '[[2,8,3,7],[4,10,5,9],[6,2,7,1],[8,4,9,3],[10,6,1,5]]',
        '[[1,5,2,4],[3,9,4,8],[5,1,6,10],[7,3,8,2],[9,7,10,6]]',
        '[[1,7,2,6],[3,10,4,11],[5,3,6,2],[7,1,8,12],[9,4,10,5],[11,9,12,8]]',
        '[[1,8,2,9],[3,11,4,10],[5,1,6,12],[7,2,8,3],[9,7,10,6],[11,5,12,4]]',
        '[[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]]',
        '[[1,9,2,8],[3,11,4,10],[5,13,6,12],[7,1,8,14],[9,3,10,2],[11,5,12,4],[13,7,14,6]]',
        '[[2,10,3,9],[4,14,5,13],[6,12,7,11],[8,2,9,1],[10,8,11,7],[12,6,13,5],[14,4,1,3]]',
        ],
    dc.pd_notation_vector: [
        '{{4, 1, 3, 2}, {2, 3, 1, 4}}',
        '{{4, 2, 3, 1}, {2, 4, 1, 3}}',
        '{{6, 1, 7, 2}, {8, 3, 5, 4}, {2, 5, 3, 6}, {4, 7, 1, 8}}',
        '{{6, 2, 7, 1}, {8, 4, 5, 3}, {2, 8, 3, 7}, {4, 6, 1, 5}}',
        '{{6, 1, 7, 2}, {10, 7, 5, 8}, {4, 5, 1, 6}, {2, 10, 3, 9}, {8, 4, 9, 3}}',
        '{{8, 2, 9, 1}, {10, 7, 5, 8}, {4, 10, 1, 9}, {2, 5, 3, 6}, {6, 3, 7, 4}}',
        '{{6, 1, 7, 2}, {10, 3, 11, 4}, {12, 8, 5, 7}, {8, 12, 9, 11}, {2, 5, 3, 6}, {4, 9, 1, 10}}',
        '{{10, 2, 11, 1}, {6, 4, 7, 3}, {12, 10, 5, 9}, {8, 6, 9, 5}, {2, 12, 3, 11}, {4, 8, 1, 7}}',
        '{{8, 1, 9, 2}, {12, 5, 7, 6}, {10, 3, 11, 4}, {4, 11, 5, 12}, {2, 7, 3, 8}, {6, 9, 1, 10}}',
        '{{10, 2, 11, 1}, {12, 6, 7, 5}, {8, 4, 9, 3}, {4, 8, 5, 7}, {2, 12, 3, 11}, {6, 10, 1, 9}}',
        '{{8, 1, 9, 2}, {2, 9, 3, 10}, {10, 3, 11, 4}, {12, 5, 7, 6}, {6, 7, 1, 8}, {4, 11, 5, 12}}'
        ],
    dc.dt_notation: [
        '',
        '[4, 6, 2]',
        '[4, 6, 8, 2]',
        '[6, 8, 10, 2, 4]',
        '[4, 8, 10, 2, 6]',
        '[4, 8, 12, 10, 2, 6]',
        '[4, 8, 10, 12, 2, 6]',
        '[4, 8, 10, 2, 12, 6]',
        '[8, 10, 12, 14, 2, 4, 6]',
        '[4, 10, 14, 12, 2, 8, 6]'
        ],
    dc.dt_code: [
        '[{4}, {2}]',
        '[{4}, {2}]',
        '[{6, 8}, {2, 4}]',
        '[{6, 8}, {4, 2}]',
        '[{6, 8}, {4, 10, 2}]',
        '[{8, 6}, {2, 10, 4}]',
        '[{6, 10}, {2, 12, 4, 8}]',
        '[{10, 6}, {8, 4, 12, 2}]',
        '[{8, 10, 12}, {2, 6, 4}]',
        '[{10, 8, 12}, {4, 6, 2}]',
        '[{8, 10, 12}, {6, 2, 4}]'
        ],
    dc.gauss_notation: [
        '',
        '{1, -2, 3, -1, 2, -3}',
        '{-1, 2, -3, 1, -4, 3, -2, 4}',
        '{-1, 2, -3, 4, -5, 1, -2, 3, -4, 5}',
        '{1, -2, 3, -1, 4, -5, 2, -3, 5, -4}',
        '{1, -2, 3, -4, 2, -1, 5, -6, 4, -3, 6, -5}',
        '{1, -2, 3, -4, 5, -6, 2, -1, 6, -3, 4, -5}',
        '{-1, 2, -3, 1, -4, 5, -2, 3, -6, 4, -5, 6}',
        '{1, -2, 3, -4, 5, -6, 7, -1, 2, -3, 4, -5, 6, -7}',
        '{-1, 2, -3, 4, -5, 6, -7, 1, -2, 7, -6, 5, -4, 3}',
        '{{1, -2}, {2, -1}}',
        '{{1, -2}, {2, -1}}',
        '{{1, -3, 2, -4}, {3, -1, 4, -2}}',
        '{{1, -3, 2, -4}, {4, -1, 3, -2}}',
        '{{1, -4, 5, -3}, {3, -1, 2, -5, 4, -2}}',
        '{{1, -4, 5, -3}, {4, -5, 2, -1, 3, -2}}',
        '{{1, -5, 2, -6}, {5, -1, 3, -4, 6, -2, 4, -3}}',
        '{{1, -5, 2, -6}, {4, -2, 6, -4, 3, -1, 5, -3}}',
        '{{1, -5, 3, -4, 2, -6}, {5, -1, 6, -3, 4, -2}}',
        '{{1, -5, 3, -4, 2, -6}, {4, -3, 6, -1, 5, -2}}',
        '{{1, -2, 3, -6, 4, -5}, {5, -1, 2, -3, 6, -4}}'
        ],
    dc.arc_notation: [
        '{{4, 2}, {3, 1}, {4, 2}, {1, 3}}',
        '{{2, 4}, {3, 1}, {2, 4}, {3, 1}}',
        '{{6, 4}, {3, 5}, {4, 2}, {1, 3}, {2, 6}, {5, 1}}',
        '{{3, 6}, {2, 5}, {6, 4}, {1, 3}, {5, 2}, {4, 1}}',
        '{{6, 2}, {1, 4}, {3, 5}, {4, 7}, {2, 6}, {7, 3}, {5, 1}}',
        '{{3, 5}, {6, 4}, {5, 2}, {7, 3}, {1, 6}, {2, 7}, {4, 1}}',
        '{{8, 4}, {3, 5}, {4, 2}, {6, 3}, {5, 7}, {1, 6}, {2, 8}, {7, 1}}',
        '{{2, 8}, {1, 7}, {8, 4}, {5, 3}, {4, 2}, {3, 6}, {7, 5}, {6, 1}}',
        '{{8, 3}, {2, 7}, {3, 1}, {4, 8}, {5, 2}, {6, 4}, {7, 5}, {1, 6}}',
        '{{3, 8}, {2, 7}, {8, 4}, {1, 3}, {5, 2}, {4, 6}, {7, 5}, {6, 1}}',
        '{{8, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 6}, {5, 7}, {6, 8}, {7, 1}}'
        ],
    dc.alternating: ['Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y'],
    dc.symmetry_type: [
        '',
        'reversible',
        'fully amphicheiral',
        'reversible',
        'reversible',
        'reversible',
        'reversible',
        'fully amphicheiral',
        'reversible',
        'reversible'
        ],
    dc.homfly_polynomial: [
        '',
        '(2*v^2-v^4)+v^2*z^2',
        '(v^(-2)-1+v^2)-z^2',
        '(3*v^4-2*v^6)+(4*v^4-v^6)*z^2+v^4*z^4',
        '(v^2+v^4-v^6)+(v^2+v^4)*z^2',
        '(v^(-2)-v^2+v^4)+(-1-v^2)*z^2',
        '(2-2*v^2+v^4)+(1-3*v^2+v^4)*z^2-v^2*z^4',
        '(-v^(-2)+3-v^2)+(-v^(-2)+3-v^2)*z^2+z^4',
        '(4*v^6-3*v^8)+(10*v^6-4*v^8)*z^2+(6*v^6-v^8)*z^4+v^6*z^6',
        '(v^2+v^6-v^8)+(v^2+v^4+v^6)*z^2'
        ],
    dc.homflypt_polynomial: [
        '1/(v^3*z)-1/(v*z)-z/v',
        'v/z-v^3/z + v*z',
        '1/(v^5*z)-1/(v^3*z)-z/v^3-z/v',
        'v^3/z-v^5/z + 3*v^3*z-v^5*z + v^3*z^3',
        '1/(v*z)-v/z-z/v^3 + (2*z)/v-v*z + z^3/v',
        '1/(v*z)-v/z-z/v^3 + (2*z)/v-v*z + z^3/v',
        '1/(v^5*z)-1/(v^3*z)-(2*z)/v^3 + z/v-v*z + z^3/v',
        'v^3/z-v^5/z + 2*v^3*z + v^5*z-v^7*z + v^3*z^3 + v^5*z^3',
        '1/(v^7*z)-1/(v^5*z) + z/v^7-(2*z)/v^5-(2*z)/v^3-z^3/v^5-z^3/v^3',
        'v^5/z-v^7/z + 2*v^3*z + 2*v^5*z-v^7*z + v^3*z^3 + v^5*z^3',
        '1/(v^7*z)-1/(v^5*z) + (3*z)/v^7-(6*z)/v^5 + z^3/v^7-(5*z^3)/v^5-z^5/v^5'
        ],
    dc.kauffman_polynomial: [
        '',
        '(-a^(-4)-2*a^(-2))*z^(0)+ (a^(-5)+ a^(-3))*z^(1)+ (a^(-4)+ a^(-2))*z^(2)',
        '(-a^(-2)-1-a^2)*z^(0)+ (-a^(-1)-a)*z^(1)+ (a^(-2)+ 2+ a^2)*z^(2)+ (a^(-1)+ a)*z^(3)',
        '(2*a^(-6)+ 3*a^(-4))*z^(0)+ (a^(-9)-a^(-7)-2*a^(-5))*z^(1)+ (a^(-8)-3*a^(-6)-4*a^(-4))*z^(2)+ (a^(-7)+ a^(-5))*z^(3)+ (a^(-6)+ a^(-4))*z^(4)',
        '(a^(-6)+ a^(-4)-a^(-2))*z^(0)+ (-2*a^(-7)-2*a^(-5))*z^(1)+ (-2*a^(-6)-a^(-4)+ a^(-2))*z^(2)+ (a^(-7)+ 2*a^(-5)+ a^(-3))*z^(3)+ (a^(-6)+ a^(-4))*z^(4)',
        '(a^(-4)+ a^(-2)-a^2)*z^(0)+ (2*a^(-3)+ 2*a^(-1))*z^(1)+ (-3*a^(-4)-4*a^(-2)+ a^2)*z^(2)+ (-3*a^(-3)-2*a^(-1)+ a)*z^(3)+ (a^(-4)+ 2*a^(-2)+ 1)*z^(4)+ (a^(-3)+ a^(-1))*z^(5)',
        '(a^(-4)+ 2*a^(-2)+ 2)*z^(0)+ (-a^(-5)-a^(-3))*z^(1)+ (a^(-6)-2*a^(-4)-6*a^(-2)-3)*z^(2)+ (2*a^(-5)-2*a^(-1))*z^(3)+ (2*a^(-4)+ 3*a^(-2)+ 1)*z^(4)+ (a^(-3)+ a^(-1))*z^(5)',
        '(a^(-2)+ 3+ a^2)*z^(0)+ (-a^(-3)-2*a^(-1)-2*a-a^3)*z^(1)+ (-3*a^(-2)-6-3*a^2)*z^(2)+ (a^(-3)+ a^(-1)+ a+ a^3)*z^(3)+ (2*a^(-2)+ 4+ 2*a^2)*z^(4)+ (a^(-1)+ a)*z^(5)',
        '(-3*a^(-8)-4*a^(-6))*z^(0)+ (a^(-13)-a^(-11)+ a^(-9)+ 3*a^(-7))*z^(1)+ (a^(-12)-2*a^(-10)+ 7*a^(-8)+ 10*a^(-6))*z^(2)+ (a^(-11)-3*a^(-9)-4*a^(-7))*z^(3)+ (a^(-10)-5*a^(-8)-6*a^(-6))*z^(4)+ (a^(-9)+ a^(-7))*z^(5)+ (a^(-8)+ a^(-6))*z^(6)',
        '(-a^(-8)-a^(-6)-a^(-2))*z^(0)+ (3*a^(-9)+ 3*a^(-7))*z^(1)+ (4*a^(-8)+ 3*a^(-6)+ a^(-2))*z^(2)+ (-4*a^(-9)-6*a^(-7)-a^(-5)+ a^(-3))*z^(3)+ (-4*a^(-8)-3*a^(-6)+ a^(-4))*z^(4)+ (a^(-9)+ 2*a^(-7)+ a^(-5))*z^(5)+ (a^(-8)+ a^(-6))*z^(6)',
        'a^2-a/z-a^3/z + a*z + a^3*z',
        'a^(-2)-1/(a^3*z)-1/(a*z) + z/a^3 + z/a',
        '-a^4 + a^3/z + a^5/z + a*z-2*a^3*z-3*a^5*z + a^2*z^2 + a^4*z^2 + a^3*z^3 + a^5*z^3',
        '-a^(-4) + 1/(a^5*z) + 1/(a^3*z) + z/a^7-(2*z)/a^5-(3*z)/a^3 + z^2/a^6 + z^2/a^4 + z^3/a^5 + z^3/a^3',
        '-1 + 1/(a*z) + a/z-(2*z)/a-4*a*z-2*a^3*z-z^2 + a^4*z^2 + z^3/a + 3*a*z^3 + 2*a^3*z^3 + z^4 + a^2*z^4',
        '-1 + 1/(a*z) + a/z-(2*z)/a-4*a*z-2*a^3*z-z^2 + a^4*z^2 + z^3/a + 3*a*z^3 + 2*a^3*z^3 + z^4 + a^2*z^4',
        '-a^4 + a^3/z + a^5/z-z/a-a^3*z-2*a^5*z-3*z^2-3*a^2*z^2 + z^3/a + a^5*z^3 + 2*z^4 + 3*a^2*z^4 + a^4*z^4 + a*z^5 + a^3*z^5',
        '-a^(-4) + 1/(a^5*z) + 1/(a^3*z)-z/a^9-z/a^5-(2*z)/a^3-(3*z^2)/a^8-(3*z^2)/a^6 + z^3/a^9 + z^3/a^3 + (2*z^4)/a^8 + (3*z^4)/a^6 + z^4/a^4 + z^5/a^7 + z^5/a^5',
        'a^6-a^5/z-a^7/z-2*a^3*z + 3*a^5*z + 3*a^7*z-2*a^9*z-a^4*z^2-2*a^6*z^2-a^8*z^2 + a^3*z^3-2*a^5*z^3-2*a^7*z^3 + a^9*z^3 + a^4*z^4 + 2*a^6*z^4 + a^8*z^4 + a^5*z^5 + a^7*z^5',
        'a^(-6)-1/(a^7*z)-1/(a^5*z)-(2*z)/a^9 + (3*z)/a^7 + (3*z)/a^5-(2*z)/a^3-z^2/a^8-(2*z^2)/a^6-z^2/a^4 + z^3/a^9-(2*z^3)/a^7-(2*z^3)/a^5 + z^3/a^3 + z^4/a^8 + (2*z^4)/a^6 + z^4/a^4 + z^5/a^7 + z^5/a^5',
        'a^6-a^5/z-a^7/z + 6*a^5*z + 4*a^7*z-a^9*z + a^11*z-3*a^6*z^2-2*a^8*z^2 + a^10*z^2-5*a^5*z^3-4*a^7*z^3 + a^9*z^3 + a^6*z^4 + a^8*z^4 + a^5*z^5 + a^7*z^5'
        ],
    dc.jones_polynomial: [
        '1',
        't+ t^3-t^4',
        't^(-2)-t^(-1)+ 1-t+ t^2',
        't^2+ t^4-t^5+ t^6-t^7',
        't-t^2+ 2*t^3-t^4+ t^5-t^6',
        't^(-2)-t^(-1)+ 2-2*t+ t^2-t^3+ t^4',
        't^(-1)-1+ 2*t-2*t^2+ 2*t^3-2*t^4+ t^5',
        '-t^(-3)+ 2*t^(-2)-2*t^(-1)+ 3-2*t+ 2*t^2-t^3',
        't^3+ t^5-t^6+ t^7-t^8+ t^9-t^10',
        't-t^2+ 2*t^3-2*t^4+ 2*t^5-t^6+ t^7-t^8',
        '-x^(-5)-x^(-1)',
        '-x-x^5',
        '-x^(-9)-x^(-5) + x^(-3)-x^(-1)',
        '-x^3-x^7 + x^9-x^11',
        'x^(-7)-2/x^5 + x^(-3)-2/x + x-x^3',
        'x^(-7)-2/x^5 + x^(-3)-2/x + x-x^3',
        '-x^(-9) + x^(-7)-3/x^5 + 2/x^3-2/x + 2*x-x^3',
        '-x^3 + x^5-3*x^7 + 2*x^9-2*x^11 + 2*x^13-x^15',
        '-x^(-15) + x^(-13)-2/x^11 + 2/x^9-2/x^7 + x^(-5)-x^(-3)',
        '-x^3 + x^5-2*x^7 + 2*x^9-2*x^11 + x^13-x^15',
        '-x^(-17) + x^(-15)-x^(-13) + x^(-11)-x^(-9)-x^(-5)'
        ],
    dc.alexander_polynomial: [
        '1',
        '1-t+ t^2',
        '1-3*t+ t^2',
        '1-t+ t^2-t^3+ t^4',
        '2-3*t+ 2*t^2',
        '2-5*t+ 2*t^2',
        '1-3*t+ 3*t^2-3*t^3+ t^4',
        '1-3*t+ 5*t^2-3*t^3+ t^4',
        '1-t+ t^2-t^3+ t^4-t^5+ t^6',
        '3-5*t+ 3*t^2'],
    dc.conway_polynomial: [
        '1',
        '1+z^2',
        '1-z^2',
        '1+3*z^2+z^4',
        '1+2*z^2',
        '1-2*z^2',
        '1-z^2-z^4',
        '1+z^2+z^4',
        '1+6*z^2+5*z^4+z^6',
        '1+3*z^2',
        '-z',
        'z',
        '-2*z',
        '2*z + z^3',
        'z^3',
        'z^3',
        '-2*z + z^3',
        '2*z + 2*z^3',
        '-3*z-2*z^3',
        '3*z + 2*z^3',
        '-3*z-4*z^3-z^5'],
    dc.khovanov_polynomial: [
        '1 + q^(-2) + 1/(q^6*t^2) + 1/(q^4*t^2)',
        '1 + q^2 + q^4*t^2 + q^6*t^2',
        '1 + q^(-2) + 1/(q^10*t^4) + 1/(q^8*t^4) + 1/(q^6*t^2) + 1/(q^2*t)',
        'q^2 + q^4 + q^6*t^2 + q^10*t^3 + q^10*t^4 + q^12*t^4',
        '2 + 2/q^2 + 1/(q^8*t^3) + 1/(q^6*t^2) + 1/(q^4*t^2) + 1/(q^2*t) + t + q^4*t^2',
        '2 + 2/q^2 + 1/(q^8*t^3) + 1/(q^6*t^2) + 1/(q^4*t^2) + 1/(q^2*t) + t + q^4*t^2',
        '1 + 2/q^2 + 1/(q^10*t^4) + 1/(q^8*t^4) + 1/(q^8*t^3) + 2/(q^6*t^2) + 1/(q^4*t^2) + 2/(q^2*t) + t + q^2*t + q^4*t^2',
        'q^2 + q^4 + q^4*t + 2*q^6*t^2 + q^8*t^2 + 2*q^10*t^3 + 2*q^10*t^4 + q^12*t^4 + q^12*t^5 + q^14*t^5 + q^16*t^6',
        'q^(-4) + q^(-2) + 1/(q^16*t^6) + 1/(q^14*t^6) + 1/(q^14*t^5) + 1/(q^12*t^4) + 1/(q^10*t^4) + 1/(q^10*t^3) + 1/(q^8*t^3) + 1/(q^8*t^2) + 1/(q^6*t^2) + 1/(q^4*t)',
        'q^2 + q^4 + q^4*t + q^6*t^2 + q^8*t^2 + q^8*t^3 + q^10*t^3 + q^10*t^4 + q^12*t^4 + q^14*t^5 + q^14*t^6 + q^16*t^6',
        'q^(-6) + q^(-4) + 1/(q^18*t^6) + 1/(q^16*t^6) + 1/(q^16*t^5) + 1/(q^12*t^4) + 1/(q^12*t^3) + 1/(q^8*t^2)'],
    dc.khovanov_unreduced_integral_polynomial: [
        '',
        'q + q^(3) + t^(2) q^(5) + t^(3) q^(9) + t^(3) q^(7) T^(2)',
        't^(-2) q^(-5) + t^(-1) q^(-1) + q^(-1) + q + t q + t^(2) q^(5) + t^(-1) q^(-3) T^(2) + t^(2) q^(3) T^(2)',
        'q^(3) + q^(5) + t^(2) q^(7) + t^(3) q^(11) + t^(4) q^(11) + t^(5) q^(15) + t^(3) q^(9) T^(2) + t^(5) q^(13) T^(2)',
        'q + q^(3) + t q^(3) + t^(2) q^(5) + t^(2) q^(7) + t^(3) q^(9) + t^(4) q^(9) + t^(5) q^(13) + t^(2) q^(5) T^(2) + t^(3) q^(7) T^(2) + t^(5) q^(11) T^(2)',
        't^(-2) q^(-5) + t^(-1) q^(-1) + 2 q^(-1) + q + t q + t q^(3) + t^(2) q^(5) + t^(3) q^(5) + t^(4) q^(9) + t^(-1) q^(-3) T^(2) + t q T^(2) + t^(2) q^(3) T^(2) + t^(4) q^(7) T^(2)',
        't^(-2) q^(-3) + t^(-1) q + 2 q + q^(3) + t q^(3) + t q^(5) + t^(2) q^(5) + t^(2) q^(7) + t^(3) q^(7) + t^(3) q^(9) + t^(4) q^(11) + t^(-1) q^(-1) T^(2) + t q^(3) T^(2) + t^(2) q^(5) T^(2) + t^(3) q^(7) T^(2) + t^(4) q^(9) T^(2)',
        't^(-3) q^(-7) + t^(-2) q^(-5) + t^(-2) q^(-3) + t^(-1) q^(-3) + t^(-1) q^(-1) + 2 q^(-1) + 2 q + t q + t q^(3) + t^(2) q^(3) + t^(2) q^(5) + t^(3) q^(7) + t^(-2) q^(-5) T^(2) + t^(-1) q^(-3) T^(2) + q^(-1) T^(2) + t q T^(2) + t^(2) q^(3) T^(2) + t^(3) q^(5) T^(2)',
        'q^(5) + q^(7) + t^(2) q^(9) + t^(3) q^(13) + t^(4) q^(13) + t^(5) q^(17) + t^(6) q^(17) + t^(7) q^(21) + t^(3) q^(11) T^(2) + t^(5) q^(15) T^(2) + t^(7) q^(19) T^(2)',
        'q + q^(3) + t q^(3) + t^(2) q^(5) + t^(2) q^(7) + t^(3) q^(7) + t^(3) q^(9) + t^(4) q^(9) + t^(4) q^(11) + t^(5) q^(13) + t^(6) q^(13) + t^(7) q^(17) + t^(2) q^(5) T^(2) + t^(3) q^(7) T^(2) + t^(4) q^(9) T^(2) + t^(5) q^(11) T^(2) + t^(7) q^(15) T^(2)'],
    dc.khovanov_reduced_integral_polynomial: [
        '',
        'q^(2) + t^(2) q^(6) + t^(3) q^(8)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 1 + t q^(2) + t^(2) q^(4)',
        'q^(4) + t^(2) q^(8) + t^(3) q^(10) + t^(4) q^(12) + t^(5) q^(14)',
        'q^(2) + t q^(4) + 2 t^(2) q^(6) + t^(3) q^(8) + t^(4) q^(10) + t^(5) q^(12)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 2 + 2 t q^(2) + t^(2) q^(4) + t^(3) q^(6) + t^(4) q^(8)',
        't^(-2) q^(-2) + t^(-1) + 2 q^(2) + 2 t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + t^(4) q^(10)',
        't^(-3) q^(-6) + 2 t^(-2) q^(-4) + 2 t^(-1) q^(-2) + 3 + 2 t q^(2) + 2 t^(2) q^(4) + t^(3) q^(6)',
        'q^(6) + t^(2) q^(10) + t^(3) q^(12) + t^(4) q^(14) + t^(5) q^(16) + t^(6) q^(18) + t^(7) q^(20)',
        'q^(2) + t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + 2 t^(4) q^(10) + t^(5) q^(12) + t^(6) q^(14) + t^(7) q^(16)'],
    dc.khovanov_reduced_rational_polynomial: [
        '',
        ' q^(2) + t^(2) q^(6) + t^(3) q^(8)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 1 + t q^(2) + t^(2) q^(4)',
        ' q^(4) + t^(2) q^(8) + t^(3) q^(10) + t^(4) q^(12) + t^(5) q^(14)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + t^(3) q^(8) + t^(4) q^(10) + t^(5) q^(12)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 2 + 2 t q^(2) + t^(2) q^(4) + t^(3) q^(6) + t^(4) q^(8)',
        't^(-2) q^(-2) + t^(-1) + 2 q^(2) + 2 t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + t^(4) q^(10)',
        't^(-3) q^(-6) + 2 t^(-2) q^(-4) + 2 t^(-1) q^(-2) + 3 + 2 t q^(2) + 2 t^(2) q^(4) + t^(3) q^(6)',
        ' q^(6) + t^(2) q^(10) + t^(3) q^(12) + t^(4) q^(14) + t^(5) q^(16) + t^(6) q^(18) + t^(7) q^(20)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + 2 t^(4) q^(10) + t^(5) q^(12) + t^(6) q^(14) + t^(7) q^(16)'],
    dc.khovanov_reduced_mod2_polynomial: [
        '',
        ' q^(2) + t^(2) q^(6) + t^(3) q^(8)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 1 + t q^(2) + t^(2) q^(4)',
        ' q^(4) + t^(2) q^(8) + t^(3) q^(10) + t^(4) q^(12) + t^(5) q^(14)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + t^(3) q^(8) + t^(4) q^(10) + t^(5) q^(12)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 2 + 2 t q^(2) + t^(2) q^(4) + t^(3) q^(6) + t^(4) q^(8)',
        't^(-2) q^(-2) + t^(-1) + 2 q^(2) + 2 t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + t^(4) q^(10)',
        't^(-3) q^(-6) + 2 t^(-2) q^(-4) + 2 t^(-1) q^(-2) + 3 + 2 t q^(2) + 2 t^(2) q^(4) + t^(3) q^(6)',
        ' q^(6) + t^(2) q^(10) + t^(3) q^(12) + t^(4) q^(14) + t^(5) q^(16) + t^(6) q^(18) + t^(7) q^(20)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + 2 t^(4) q^(10) + t^(5) q^(12) + t^(6) q^(14) + t^(7) q^(16)'],
    dc.khovanov_odd_integral_polynomial: [
        '',
        'q^(2) + t^(2) q^(6) + t^(3) q^(8)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 1 + t q^(2) + t^(2) q^(4)',
        'q^(4) + t^(2) q^(8) + t^(3) q^(10) + t^(4) q^(12) + t^(5) q^(14)',
        'q^(2) + t q^(4) + 2 t^(2) q^(6) + t^(3) q^(8) + t^(4) q^(10) + t^(5) q^(12)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 2 + 2 t q^(2) + t^(2) q^(4) + t^(3) q^(6) + t^(4) q^(8)',
        't^(-2) q^(-2) + t^(-1) + 2 q^(2) + 2 t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + t^(4) q^(10)',
        't^(-3) q^(-6) + 2 t^(-2) q^(-4) + 2 t^(-1) q^(-2) + 3 + 2 t q^(2) + 2 t^(2) q^(4) + t^(3) q^(6)',
        'q^(6) + t^(2) q^(10) + t^(3) q^(12) + t^(4) q^(14) + t^(5) q^(16) + t^(6) q^(18) + t^(7) q^(20)',
        'q^(2) + t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + 2 t^(4) q^(10) + t^(5) q^(12) + t^(6) q^(14) + t^(7) q^(16)'],
    dc.khovanov_odd_rational_polynomial: [
        '',
        ' q^(2) + t^(2) q^(6) + t^(3) q^(8)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 1 + t q^(2) + t^(2) q^(4)',
        ' q^(4) + t^(2) q^(8) + t^(3) q^(10) + t^(4) q^(12) + t^(5) q^(14)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + t^(3) q^(8) + t^(4) q^(10) + t^(5) q^(12)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 2 + 2 t q^(2) + t^(2) q^(4) + t^(3) q^(6) + t^(4) q^(8)',
        't^(-2) q^(-2) + t^(-1) + 2 q^(2) + 2 t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + t^(4) q^(10)',
        't^(-3) q^(-6) + 2 t^(-2) q^(-4) + 2 t^(-1) q^(-2) + 3 + 2 t q^(2) + 2 t^(2) q^(4) + t^(3) q^(6)',
        ' q^(6) + t^(2) q^(10) + t^(3) q^(12) + t^(4) q^(14) + t^(5) q^(16) + t^(6) q^(18) + t^(7) q^(20)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + 2 t^(4) q^(10) + t^(5) q^(12) + t^(6) q^(14) + t^(7) q^(16)'],
    dc.khovanov_odd_mod2_polynomial: [
        '',
        ' q^(2) + t^(2) q^(6) + t^(3) q^(8)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 1 + t q^(2) + t^(2) q^(4)',
        ' q^(4) + t^(2) q^(8) + t^(3) q^(10) + t^(4) q^(12) + t^(5) q^(14)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + t^(3) q^(8) + t^(4) q^(10) + t^(5) q^(12)',
        't^(-2) q^(-4) + t^(-1) q^(-2) + 2 + 2 t q^(2) + t^(2) q^(4) + t^(3) q^(6) + t^(4) q^(8)',
        't^(-2) q^(-2) + t^(-1) + 2 q^(2) + 2 t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + t^(4) q^(10)',
        't^(-3) q^(-6) + 2 t^(-2) q^(-4) + 2 t^(-1) q^(-2) + 3 + 2 t q^(2) + 2 t^(2) q^(4) + t^(3) q^(6)',
        ' q^(6) + t^(2) q^(10) + t^(3) q^(12) + t^(4) q^(14) + t^(5) q^(16) + t^(6) q^(18) + t^(7) q^(20)',
        ' q^(2) + t q^(4) + 2 t^(2) q^(6) + 2 t^(3) q^(8) + 2 t^(4) q^(10) + t^(5) q^(12) + t^(6) q^(14) + t^(7) q^(16)']
}
