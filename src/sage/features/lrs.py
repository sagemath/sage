# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of ``lrslib``
"""

import os
import subprocess

from . import Executable, FeatureTestResult
from sage.cpython.string import str_to_bytes, bytes_to_str


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
        tf_name = tmp_filename()
        with open(tf_name, 'wb') as tf:
            tf.write(str_to_bytes("V-representation\nbegin\n 1 1 rational\n 1 \nend\nvolume"))
        devnull = open(os.devnull, 'wb')
        command = ['lrs', tf_name]
        try:
            lines = bytes_to_str(subprocess.check_output(command, stderr=devnull))
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                reason="Call to `{command}` failed with exit code {e.returncode}.".format(command=" ".join(command), e=e))

        expected_list = ["Volume= 1", "Volume=1"]
        if all(lines.find(expected) == -1 for expected in expected_list):
            print(lines)
            return FeatureTestResult(self, False,
                reason="Output of `{command}` did not contain the expected result {expected}.".format(
                    command=" ".join(command),
                    expected=" or ".join(expected_list)))

        return FeatureTestResult(self, True)


def all_features():
    return [Lrs()]
