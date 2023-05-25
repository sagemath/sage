r"""
Features for testing the presence of GAP packages
"""
# *****************************************************************************
#       Copyright (C) 2016 Julian Rüth
#                     2018 Jeroen Demeyer
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import Feature, FeatureTestResult, PythonModule
from .join_feature import JoinFeature


class GapPackage(Feature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of a GAP package.

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: GapPackage("grape", spkg="gap_packages")
        Feature('gap_package_grape')
    """
    def __init__(self, package, **kwds):
        r"""
        TESTS::

            sage: from sage.features.gap import GapPackage
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), GapPackage)
            True
        """
        Feature.__init__(self, f"gap_package_{package}", **kwds)
        self.package = package

    def _is_present(self):
        r"""
        Return whether the package is available in GAP.

        This does not check whether this package is functional.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape", spkg="gap_packages")._is_present()  # optional - gap_packages
            FeatureTestResult('gap_package_grape', True)
        """
        from sage.libs.gap.libgap import libgap
        command = 'TestPackageAvailability("{package}")'.format(package=self.package)
        presence = libgap.eval(command)
        if presence:
            return FeatureTestResult(self, True,
                    reason="`{command}` evaluated to `{presence}` in GAP.".format(command=command, presence=presence))
        else:
            return FeatureTestResult(self, False,
                    reason="`{command}` evaluated to `{presence}` in GAP.".format(command=command, presence=presence))


class sage__libs__gap(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.gap`
    (the library interface to GAP) and :mod:`sage.interfaces.gap` (the pexpect
    interface to GAP). By design, we do not distinguish between these two, in order
    to facilitate the conversion of code from the pexpect interface to the library
    interface.

    EXAMPLES::

        sage: from sage.features.gap import sage__libs__gap
        sage: sage__libs__gap().is_present()                       # optional - sage.libs.gap
        FeatureTestResult('sage.libs.gap', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.gap import sage__libs__gap
            sage: isinstance(sage__libs__gap(), sage__libs__gap)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.gap',
                             [PythonModule('sage.libs.gap.libgap'),
                              PythonModule('sage.interfaces.gap')])


def all_features():
    return [sage__libs__gap()]
