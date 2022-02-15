# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of ``lrslib``
"""

import subprocess

from . import Executable, FeatureTestResult


class Lrs(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the ``lrs``
    binary which comes as a part of ``lrslib``.

    EXAMPLES::

        sage: from sage.features.lrs import Lrs
        sage: Lrs().is_present()  # optional - lrslib
        FeatureTestResult('lrslib', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.lrs import Lrs
            sage: isinstance(Lrs(), Lrs)
            True
        """
        Executable.__init__(self, "lrslib", executable="lrs", spkg="lrslib",
                            url="http://cgm.cs.mcgill.ca/~avis/C/lrs.html")

    def is_functional(self):
        r"""
        Test whether ``lrs`` works on a trivial input.

        EXAMPLES::

            sage: from sage.features.lrs import Lrs
            sage: Lrs().is_functional()  # optional - lrslib
            FeatureTestResult('lrslib', True)
        """
        from sage.misc.temporary_file import tmp_filename

        # Check #1
        tf_name = tmp_filename()
        with open(tf_name, 'w') as tf:
            tf.write("V-representation\nbegin\n 1 1 rational\n 1 \nend\nvolume")
        command = ['lrs', tf_name]
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

        # Check #2
        # Checking whether `lrsnash` can handle the new input format
        # This test is currently done in build/pkgs/lrslib/spkg-configure.m4
        tf_name = tmp_filename()
        with open(tf_name, 'w') as tf:
            tf.write("1 1\n \n 0\n \n 0\n")
        command = ['lrsnash', tf_name]
        try:
            result = subprocess.run(command, capture_output=True, text=True)
        except OSError as e:
            return FeatureTestResult(self, False, reason='Running command "{}" '
                        'raised an OSError "{}" '.format(' '.join(command), e))
        if result.returncode:
            return FeatureTestResult(self, False, reason='Running command "{}" '
                        'returned non-zero exit status "{}" with stderr '
                        '"{}" and stdout "{}".'.format(' '.join(result.args),
                                                        result.returncode,
                                                        result.stderr.strip(),
                                                        result.stdout.strip()))

        return FeatureTestResult(self, True)


def all_features():
    return [Lrs()]
