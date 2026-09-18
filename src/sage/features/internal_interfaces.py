r"""
Features for testing whether pip installable interfaces to ``mathics``, ``regina``, etc. are functional

EXAMPLES::

    sage: from sage.features.internal_interfaces import Mathics
    sage: F = Mathics
    sage: F.hide()
    sage: mathics(~7)
    Traceback (most recent call last):
    ...
    sage.features.FeatureNotPresentError: mathics is not available.
    Feature `mathics` is hidden.
    Use method `unhide` to make it available again.
    sage: F.unhide()

    sage: from sage.features.internal_interfaces import Regina
    sage: F = Regina
    sage: F.hide()
    sage: regina(~7)
    Traceback (most recent call last):
    ...
    sage.features.FeatureNotPresentError: regina is not available.
    Feature `regina` is hidden.
    Use method `unhide` to make it available again.
    sage: F.unhide()

    sage: from sage.features.internal_interfaces import SnapPy
    sage: F = SnapPy
    sage: F.hide()
    sage: snappy(~7)
    Traceback (most recent call last):
    ...
    sage.features.FeatureNotPresentError: snappy is not available.
    Feature `snappy` is hidden.
    Use method `unhide` to make it available again.
    sage: F.unhide()
"""

# ****************************************************************************
#       Copyright (C) 2026 Sebastian Oehms
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import PythonModule

Mathics = PythonModule('mathics', spkg='mathics')
Regina = PythonModule('regina', spkg='regina')
SnapPy = PythonModule('snappy', spkg='snappy')


def all_features():
    r"""
    Return features corresponding to internal interfaces.

    EXAMPLES::

        sage: from sage.features.internal_interfaces import all_features
        sage: list(all_features())
        [Feature('mathics'), Feature('regina'), Feature('snappy')]
    """
    return [Mathics, Regina, SnapPy]
