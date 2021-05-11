# -*- coding: utf-8 -*-
r"""
Testing for databases at runtime
"""

import os

from . import StaticFile
from sage.env import CREMONA_MINI_DATA_DIR, CREMONA_LARGE_DATA_DIR

CREMONA_DATA_DIRS = set([CREMONA_MINI_DATA_DIR, CREMONA_LARGE_DATA_DIR])


class DatabaseCremona(StaticFile):
    r"""
    A :class:`Feature` which describes the presence of John Cremona's
    database of elliptic curves.

    INPUT:

    - ``name`` -- either ``'cremona'`` (the default) for the full large
      database or ``'cremona_mini'`` for the small database.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseCremona
        sage: DatabaseCremona('cremona_mini').is_present()
        FeatureTestResult("Cremona's database of elliptic curves", True)
        sage: DatabaseCremona().is_present()  # optional: database_cremona_ellcurve
        FeatureTestResult("Cremona's database of elliptic curves", True)
    """
    def __init__(self, name="cremona", spkg="database_cremona_ellcurve"):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseCremona
            sage: isinstance(DatabaseCremona(), DatabaseCremona)
            True
        """
        StaticFile.__init__(self, "Cremona's database of elliptic curves",
                            filename='{}.db'.format(name.replace(' ', '_')),
                            search_path=CREMONA_DATA_DIRS,
                            spkg=spkg,
                            url="https://github.com/JohnCremona/ecdata")

class DatabaseJones(StaticFile):
    r"""
    A :class:`Feature` which describes the presence of John Jones's tables of number fields.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseJones
        sage: bool(DatabaseJones().is_present())  # optional: database_jones_numfield
        True
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseJones
            sage: isinstance(DatabaseJones(), DatabaseJones)
            True
        """
        StaticFile.__init__(self, "John Jones's tables of number fields",
                            filename='jones/jones.sobj',
                            spkg="database_jones_numfield")
