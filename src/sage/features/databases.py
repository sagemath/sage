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

import atexit
from contextlib import ExitStack
from functools import cache
from importlib.resources import as_file, files
from os.path import join

from sage.env import sage_data_paths
from sage.features import PythonModule, StaticFile
from sage.features.join_feature import JoinFeature

# Keep extracted resources alive for callers that retain filenames or open
# database connections. With installed wheels, as_file uses the original path.
_data_resource_contexts = ExitStack()
atexit.register(_data_resource_contexts.close)


@cache
def _packaged_data_path(package, subdirectory):
    """
    Return a packaged data directory, retaining any temporary extraction.
    """
    resource = files(package).joinpath(subdirectory)
    return str(_data_resource_contexts.enter_context(as_file(resource)))


class DatabaseCremona(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of :ref:`John Cremona's
    database of elliptic curves <spkg_database_cremona_ellcurve>`.

    INPUT:

    - ``name`` -- either ``'cremona'`` (the default) for the full large
      database or ``'cremona_mini'`` for the small database

    EXAMPLES::

        sage: from sage.features.databases import DatabaseCremona
        sage: DatabaseCremona('cremona_mini', type='standard').is_present()
        FeatureTestResult('database_cremona_mini_ellcurve', True)
        sage: DatabaseCremona().is_present()                                    # optional - database_cremona_ellcurve
        FeatureTestResult('database_cremona_ellcurve', True)
    """

    def __init__(
        self, name="cremona", spkg="database_cremona_ellcurve", type="optional"
    ):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseCremona
            sage: isinstance(DatabaseCremona(), DatabaseCremona)
            True
        """
        self.filename = f"{name}.db"
        if name == "cremona_mini":
            feature = PythonModule(
                "elliptic_curve_databases", spkg="elliptic_curves", type="standard"
            )
        else:
            from sage.env import CREMONA_LARGE_DATA_DIR

            feature = StaticFile(
                f"database_{name}_ellcurve",
                filename=self.filename,
                search_path=CREMONA_LARGE_DATA_DIR or tuple(sage_data_paths("cremona")),
                spkg="database_cremona_ellcurve",
                type="optional",
                url="https://github.com/JohnCremona/ecdata",
            )

        JoinFeature.__init__(
            self,
            f"database_{name}_ellcurve",
            [feature],
            url="https://github.com/JohnCremona/ecdata",
            description="Cremona's database of elliptic curves",
        )

    def absolute_filename(self):
        """
        Return the database filename, retaining any extracted resource.

        EXAMPLES::

            sage: from pathlib import Path
            sage: from sage.features.databases import DatabaseCremona
            sage: Path(DatabaseCremona('cremona_mini').absolute_filename()).is_file()
            True
            sage: Path(DatabaseCremona().absolute_filename()).is_file()  # optional - database_cremona_ellcurve
            True
        """
        self.require()
        if self.filename == "cremona_mini.db":
            return join(_packaged_data_path("elliptic_curve_databases.data", "cremona"), self.filename)
        return self._features[0].absolute_filename()


class DatabaseEllcurves(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of
    William Stein's database of interesting curves.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseEllcurves
        sage: bool(DatabaseEllcurves().is_present())  # optional - database_ellcurves
        True
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseEllcurves
            sage: isinstance(DatabaseEllcurves(), DatabaseEllcurves)
            True
        """
        self.filename = "rank0"
        JoinFeature.__init__(
            self,
            "database_ellcurves",
            [PythonModule("elliptic_curve_databases", spkg="elliptic_curves", type="standard")],
            description="William Stein's database of interesting curves",
        )

    def absolute_filename(self):
        """
        Return the rank-zero data filename, retaining its extracted directory.

        EXAMPLES::

            sage: from pathlib import Path
            sage: from sage.features.databases import DatabaseEllcurves
            sage: Path(DatabaseEllcurves().absolute_filename()).is_file()  # optional - database_ellcurves
            True
        """
        self.require()
        return join(_packaged_data_path("elliptic_curve_databases.data", "ellcurves"), self.filename)


class DatabaseGraphs(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of
    the graphs database.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseGraphs
        sage: bool(DatabaseGraphs().is_present())  # optional - database_graphs
        True
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseGraphs
            sage: isinstance(DatabaseGraphs(), DatabaseGraphs)
            True
        """
        self.filename = "graphs.db"
        JoinFeature.__init__(
            self,
            "database_graphs",
            [PythonModule("graph_databases", spkg="graphs", type="standard")],
            description="A database of graphs",
        )

    def absolute_filename(self):
        """
        Return the database filename, retaining its extracted directory.

        EXAMPLES::

            sage: from pathlib import Path
            sage: from sage.features.databases import DatabaseGraphs
            sage: Path(DatabaseGraphs().absolute_filename()).is_file()  # optional - database_graphs
            True
        """
        self.require()
        return join(_packaged_data_path("graph_databases.data", ""), self.filename)


class DatabaseJones(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of
    :ref:`John Jones's tables of number fields <spkg_database_jones_numfield>`.

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
        StaticFile.__init__(
            self,
            "database_jones_numfield",
            filename="jones.sobj",
            search_path=sage_data_paths("jones"),
            spkg="database_jones_numfield",
            description="John Jones's tables of number fields",
        )


class DatabaseKnotInfo(PythonModule):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the
    :ref:`package providing the KnotInfo and LinkInfo databases <spkg_database_knotinfo>`.

    The homes of these databases are the
    web-pages `KnotInfo <https://knotinfo.org/>`__ and
    `LinkInfo <https://link-info-repo.onrender.com/>`__.

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
        PythonModule.__init__(self, "database_knotinfo", spkg="database_knotinfo")


class DatabaseMatroids(PythonModule):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of
    :ref:`Yoshitake Matsumoto's Database of Matroids <spkg_matroid_database>`.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseMatroids
        sage: DatabaseMatroids().is_present()                                           # optional - matroid_database
        FeatureTestResult('matroid_database', True)

    REFERENCES:

    [Mat2012]_
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseMatroids
            sage: isinstance(DatabaseMatroids(), DatabaseMatroids)
            True
        """
        PythonModule.__init__(self, "matroid_database", spkg="matroid_database")


class DatabaseCubicHecke(PythonModule):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the
    :ref:`Cubic Hecke algebra database package <spkg_database_cubic_hecke>`.

    The home of this database is the
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
        PythonModule.__init__(self, "database_cubic_hecke", spkg="database_cubic_hecke")


class DatabaseReflexivePolytopes(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the
    :ref:`PALP databases of reflexive three-dimensional <spkg_polytopes_db>`
    and :ref:`four-dimensional lattice polytopes <spkg_polytopes_db_4d>`.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseReflexivePolytopes
        sage: bool(DatabaseReflexivePolytopes().is_present())                   # optional - polytopes_db
        True
        sage: bool(DatabaseReflexivePolytopes('polytopes_db_4d').is_present())  # optional - polytopes_db_4d
        True
    """

    def __init__(self, name="polytopes_db"):
        """
        TESTS::

            sage: from sage.features.databases import DatabaseReflexivePolytopes
            sage: isinstance(DatabaseReflexivePolytopes(), DatabaseReflexivePolytopes)
            True
            sage: DatabaseReflexivePolytopes().filename
            'Full3d'
            sage: DatabaseReflexivePolytopes('polytopes_db_4d').filename
            'Hodge4d'
        """
        if name == "polytopes_db_4d":
            from sage.env import POLYTOPE_DATA_DIR

            self.filename = "Hodge4d"
            feature = StaticFile(
                name,
                filename=self.filename,
                search_path=POLYTOPE_DATA_DIR or tuple(sage_data_paths("reflexive_polytopes")),
                spkg="polytopes_db_4d",
            )
        else:
            self.filename = "Full3d"
            feature = PythonModule(
                "polytope_databases", spkg="polytopes_db", type="standard"
            )

        JoinFeature.__init__(self, name, [feature])

    def absolute_filename(self):
        """
        Return the PALP data directory, retaining any extracted resource.

        EXAMPLES::

            sage: from pathlib import Path
            sage: from sage.features.databases import DatabaseReflexivePolytopes
            sage: Path(DatabaseReflexivePolytopes().absolute_filename()).is_dir()  # optional - polytopes_db
            True
            sage: Path(DatabaseReflexivePolytopes('polytopes_db_4d').absolute_filename()).is_dir()  # optional - polytopes_db_4d
            True
        """
        self.require()
        if self.filename == "Full3d":
            return join(_packaged_data_path("polytope_databases", "data"), self.filename)
        return self._features[0].absolute_filename()


def all_features():
    return [
        PythonModule("conway_polynomials", spkg="conway_polynomials", type="standard"),
        DatabaseCremona(),
        DatabaseCremona("cremona_mini", type="standard"),
        DatabaseEllcurves(),
        DatabaseGraphs(),
        DatabaseJones(),
        DatabaseKnotInfo(),
        DatabaseMatroids(),
        DatabaseCubicHecke(),
        DatabaseReflexivePolytopes(),
        DatabaseReflexivePolytopes("polytopes_db_4d"),
    ]
