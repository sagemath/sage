# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of Sloane Online Encyclopedia of Integer Sequences
"""

# ****************************************************************************
#       Copyright (C) 2024 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Feature


class SloaneOEIS(Feature):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of
    the Sloane Online Encyclopedia of Integer Sequences.

    EXAMPLES::

        sage: from sage.features.sloane_database import SloaneOEIS
        sage: bool(SloaneOEIS().is_present())  # optional - sloane_database
        True
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sloane_database import SloaneOEIS
            sage: isinstance(SloaneOEIS(), SloaneOEIS)
            True
        """
        Feature.__init__(self, name='sloane_database',
                         description='Sloane Online Encyclopedia of Integer Sequences')

    def _is_present(self):
        r"""
        Return whether the database is available.

        EXAMPLES::

            sage: from sage.features.sloane_database import SloaneOEIS
            sage: bool(SloaneOEIS().is_present())  # optional - !sloane_database
            False
        """
        try:
            from sage.databases.sloane import SloaneEncyclopedia
        except ImportError:
            return False
        return SloaneEncyclopedia.is_installed()


def all_features():
    return [SloaneOEIS()]
