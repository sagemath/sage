# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of ``csdp``
"""

import os
import re
import subprocess
from sage.cpython.string import bytes_to_str

from . import Executable, FeatureTestResult


class CSDP(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for the ``theta`` binary
    of CSDP.

    EXAMPLES::

        sage: from sage.features.csdp import CSDP
        sage: CSDP().is_present()  # optional - csdp
        FeatureTestResult('csdp', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.csdp import CSDP
            sage: isinstance(CSDP(), CSDP)
            True
        """
        Executable.__init__(self, name="csdp", spkg="csdp", executable="theta",
                                url="https://github.com/dimpase/csdp")

    def is_functional(self):
        r"""
        Check whether ``theta`` works on a trivial example.

        EXAMPLES::

            sage: from sage.features.csdp import CSDP
            sage: CSDP().is_functional()  # optional - csdp
            FeatureTestResult('csdp', True)
        """
        from sage.misc.temporary_file import tmp_filename
        tf_name = tmp_filename()
        with open(tf_name, 'wb') as tf:
            tf.write("2\n1\n1 1".encode())
        with open(os.devnull, 'wb') as devnull:
            command = ['theta', tf_name]
            try:
                lines = subprocess.check_output(command, stderr=devnull)
            except subprocess.CalledProcessError as e:
                return FeatureTestResult(self, False,
                    reason="Call to `{command}` failed with exit code {e.returncode}."
                                             .format(command=" ".join(command), e=e))

        result = bytes_to_str(lines).strip().split('\n')[-1]
        match = re.match("^The Lovasz Theta Number is (.*)$", result)
        if match is None:
            return FeatureTestResult(self, False,
                reason="Last line of the output of `{command}` did not have the expected format."
                                         .format(command=" ".join(command)))

        return FeatureTestResult(self, True)


def all_features():
    return [CSDP()]
