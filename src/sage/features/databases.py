r"""
Features for testing the presence of various databases
"""

# *****************************************************************************
#       Copyright (C) 2016      Julian Rüth
#                     2018-2019 Jeroen Demeyer
#                     2018      Timo Kaufmann
#                     2020-2022 Matthias Koeppe
#                     2020-2022 Sebastian Oehms
#                     2021      Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************


from . import StaticFile, PythonModule
from sage.env import (
    CONWAY_POLYNOMIALS_DATA_DIR,
    CREMONA_MINI_DATA_DIR, CREMONA_LARGE_DATA_DIR,
    POLYTOPE_DATA_DIR)


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

class DatabaseCubicHecke(PythonModule):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the databases at the
    web-page `Cubic Hecke algebra on 4 strands <http://www.lamfa.u-picardie.fr/marin/representationH4-en.html>`__
    of Ivan Marin.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseCubicHecke
        sage: DatabaseCubicHecke().is_present()  # optional - database_cubic_hecke
        FeatureTestResult('database_cubic_hecke', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseCubicHecke
            sage: isinstance(DatabaseCubicHecke(), DatabaseCubicHecke)
            True
        """
        PythonModule.__init__(self, 'database_cubic_hecke', spkg='database_cubic_hecke')

class DatabaseReflexivePolytopes(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the PALP database
    of reflexive lattice polytopes.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseReflexivePolytopes
        sage: bool(DatabaseReflexivePolytopes().is_present())                              # optional - polytopes_db
        True
        sage: bool(DatabaseReflexivePolytopes('polytopes_db_4d', 'Hodge4d').is_present())  # optional - polytopes_db_4d
        True
    """
    def __init__(self, name='polytopes_db', dirname='Full3D'):
        """
        TESTS::

            sage: from sage.features.databases import DatabaseReflexivePolytopes
            sage: isinstance(DatabaseReflexivePolytopes(), DatabaseReflexivePolytopes)
            True
        """
        StaticFile.__init__(self, name, dirname,
                            search_path=[POLYTOPE_DATA_DIR])


def all_features():
    return [DatabaseConwayPolynomials(),
            DatabaseCremona(), DatabaseCremona('cremona_mini'),
            DatabaseJones(),
            DatabaseKnotInfo(),
            DatabaseCubicHecke(),
            DatabaseReflexivePolytopes(),
            DatabaseReflexivePolytopes('polytopes_db_4d', 'Hodge4d')]
