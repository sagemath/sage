# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``lrslib``
"""

# *****************************************************************************
#       Copyright (C) 2016      Julian Rüth
#                     2018      Jeroen Demeyer
#                     2021-2022 Matthias Koeppe
#                     2021      Kwankyu Lee
#                     2022      Sébastien Labbé
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import subprocess

from . import Executable, FeatureTestResult
from .join_feature import JoinFeature


class Lrs(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the ``lrs``
    binary which comes as a part of ``lrslib``.

    EXAMPLES::

        sage: from sage.features.lrs import Lrs
        sage: Lrs().is_present()  # optional - lrslib
        FeatureTestResult('lrs', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.lrs import Lrs
            sage: isinstance(Lrs(), Lrs)
            True
        """
        Executable.__init__(self, "lrs", executable='lrs', spkg='lrslib',
                            url='http://cgm.cs.mcgill.ca/~avis/C/lrs.html')

    def is_functional(self):
        r"""
        Test whether ``lrs`` works on a trivial input.

        EXAMPLES::

            sage: from sage.features.lrs import Lrs
            sage: Lrs().is_functional()  # optional - lrslib
            FeatureTestResult('lrs', True)
        """
        from sage.misc.temporary_file import tmp_filename

        tf_name = tmp_filename()
        with open(tf_name, 'w') as tf:
            tf.write("V-representation\nbegin\n 1 1 rational\n 1 \nend\nvolume")
        command = [self.absolute_filename(), tf_name]
        try:
            result = subprocess.run(command, capture_output=True, text=True)
        except OSError as e:
            return FeatureTestResult(self, False, reason='Running command "{}" '
                        'raised an OSError "{}" '.format(' '.join(command), e))

        if result.returncode:
            return FeatureTestResult(self, False,
                reason="Call to `{command}` failed with exit code {result.returncode}.".format(command=" ".join(command), result=result))

        expected_list = ["Volume= 1", "Volume=1"]
        if all(result.stdout.find(expected) == -1 for expected in expected_list):
            return FeatureTestResult(self, False,
                reason="Output of `{command}` did not contain the expected result {expected}; output: {result.stdout}".format(
                    command=" ".join(command),
                    expected=" or ".join(expected_list),
                    result=result))

        return FeatureTestResult(self, True)


class LrsNash(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the ``lrsnash``
    binary which comes as a part of ``lrslib``.

    EXAMPLES::

        sage: from sage.features.lrs import LrsNash
        sage: LrsNash().is_present()  # optional - lrslib
        FeatureTestResult('lrsnash', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.lrs import LrsNash
            sage: isinstance(LrsNash(), LrsNash)
            True
        """
        Executable.__init__(self, "lrsnash", executable='lrsnash', spkg='lrslib',
                            url='http://cgm.cs.mcgill.ca/~avis/C/lrs.html')

    def is_functional(self):
        r"""
        Test whether ``lrsnash`` works on a trivial input.

        EXAMPLES::

            sage: from sage.features.lrs import LrsNash
            sage: LrsNash().is_functional()  # optional - lrslib
            FeatureTestResult('lrsnash', True)
        """
        from sage.misc.temporary_file import tmp_filename

        # Checking whether `lrsnash` can handle the new input format
        # This test is currently done in build/pkgs/lrslib/spkg-configure.m4
        tf_name = tmp_filename()
        with open(tf_name, 'w') as tf:
            tf.write("1 1\n \n 0\n \n 0\n")
        command = [self.absolute_filename(), tf_name]
        try:
            result = subprocess.run(command, capture_output=True, text=True)
        except OSError as e:
            return FeatureTestResult(self, False, reason='Running command "{}" '
                        'raised an OSError "{}" '.format(' '.join(command), e))
        if result.returncode:
            return FeatureTestResult(self, False, reason='Running command "{}" '
                        'returned nonzero exit status "{}" with stderr '
                        '"{}" and stdout "{}".'.format(' '.join(result.args),
                                                        result.returncode,
                                                        result.stderr.strip(),
                                                        result.stdout.strip()))

        return FeatureTestResult(self, True)


class Lrslib(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the executables
    :class:`lrs <Lrs>` and :class:`lrsnash <LrsNash>` provided by the :ref:`lrslib <spkg_lrslib>` package.

    EXAMPLES::

        sage: from sage.features.lrs import Lrslib
        sage: Lrslib().is_present()  # optional - lrslib
        FeatureTestResult('lrslib', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.lrs import Lrslib
            sage: isinstance(Lrslib(), Lrslib)
            True
        """
        JoinFeature.__init__(self, "lrslib",
                             (Lrs(), LrsNash()))


def all_features():
    return [Lrslib()]
