# -*- coding: utf-8 -*-
r"""
Testing for databases at runtime
"""

import os

from . import StaticFile


class DatabaseCremona(StaticFile):
    r"""
    A :class:`Feature` which describes the presence of John Cremona's database
    of elliptic curves.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseCremona
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
        filename = name.replace(' ', '_') + ".db"
        StaticFile.__init__(self, "Cremona's database of elliptic curves",
                            filename=os.path.join("cremona", filename),
                            spkg=spkg,
                            url="https://github.com/JohnCremona/ecdata")
