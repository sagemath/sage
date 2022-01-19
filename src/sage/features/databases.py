# -*- coding: utf-8 -*-
r"""
Features for testing the presence of various databases
"""


from . import StaticFile, PythonModule
from sage.env import (
    CONWAY_POLYNOMIALS_DATA_DIR,
    CREMONA_MINI_DATA_DIR, CREMONA_LARGE_DATA_DIR)


class DatabaseConwayPolynomials(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of Frank Luebeck's
    database of Conway polynomials.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseConwayPolynomials
        sage: DatabaseConwayPolynomials().is_present()
        FeatureTestResult('conway_polynomials', True)
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseConwayPolynomials
            sage: isinstance(DatabaseConwayPolynomials(), DatabaseConwayPolynomials)
            True
        """
        if CONWAY_POLYNOMIALS_DATA_DIR:
            search_path = [CONWAY_POLYNOMIALS_DATA_DIR]
        else:
            search_path = []
        StaticFile.__init__(self, "conway_polynomials",
                            filename='conway_polynomials.p',
                            search_path=search_path,
                            spkg='conway_polynomials',
                            description="Frank Luebeck's database of Conway polynomials")


CREMONA_DATA_DIRS = set([CREMONA_MINI_DATA_DIR, CREMONA_LARGE_DATA_DIR])


class DatabaseCremona(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of John Cremona's
    database of elliptic curves.

    INPUT:

    - ``name`` -- either ``'cremona'`` (the default) for the full large
      database or ``'cremona_mini'`` for the small database.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseCremona
        sage: DatabaseCremona('cremona_mini').is_present()
        FeatureTestResult('database_cremona_mini_ellcurve', True)
        sage: DatabaseCremona().is_present()  # optional - database_cremona_ellcurve
        FeatureTestResult('database_cremona_ellcurve', True)
    """
    def __init__(self, name="cremona", spkg="database_cremona_ellcurve"):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseCremona
            sage: isinstance(DatabaseCremona(), DatabaseCremona)
            True
        """
        StaticFile.__init__(self, f"database_{name}_ellcurve",
                            filename='{}.db'.format(name.replace(' ', '_')),
                            search_path=CREMONA_DATA_DIRS,
                            spkg=spkg,
                            url="https://github.com/JohnCremona/ecdata",
                            description="Cremona's database of elliptic curves")


class DatabaseJones(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of John Jones's tables of number fields.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseJones
        sage: bool(DatabaseJones().is_present())  # optional - database_jones_numfield
        True
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseJones
            sage: isinstance(DatabaseJones(), DatabaseJones)
            True
        """
        StaticFile.__init__(self, "database_jones_numfield",
                            filename='jones/jones.sobj',
                            spkg="database_jones_numfield",
                            description="John Jones's tables of number fields")


class DatabaseKnotInfo(PythonModule):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the databases at the
    web-pages `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and
    `LinkInfo <https://linkinfo.sitehost.iu.edu>`__.



    EXAMPLES::

        sage: from sage.features.databases import DatabaseKnotInfo
        sage: DatabaseKnotInfo().is_present()  # optional - database_knotinfo
        FeatureTestResult('database_knotinfo', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseKnotInfo
            sage: isinstance(DatabaseKnotInfo(), DatabaseKnotInfo)
            True
        """
        PythonModule.__init__(self, 'database_knotinfo', spkg='database_knotinfo')


def all_features():
    return [DatabaseConwayPolynomials(),
            DatabaseCremona(), DatabaseCremona('cremona_mini'),
            DatabaseJones(),
            DatabaseKnotInfo()]
