# -*- coding: utf-8 -*-
r"""
Check for Kenzo
"""

from sage.libs.ecl import ecl_eval
from . import Feature, FeatureTestResult

class Kenzo(Feature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``Kenzo``.

    EXAMPLES::

        sage: from sage.features.kenzo import Kenzo
        sage: Kenzo().is_present()  # optional - kenzo
        FeatureTestResult('Kenzo', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.kenzo import Kenzo
            sage: isinstance(Kenzo(), Kenzo)
            True
        """
        Feature.__init__(self, name="Kenzo", spkg="kenzo",
                         url="https://github.com/miguelmarco/kenzo/")

    def _is_present(self):
        r"""
        Check whether Kenzo is installed and works.

        EXAMPLES::

            sage: from sage.features.kenzo import Kenzo
            sage: Kenzo()._is_present()  # optional - kenzo
            FeatureTestResult('Kenzo', True)
        """
        # Redirection of ECL and Maxima stdout to /dev/null
        # This is also done in the Maxima library, but we
        # also do it here for redundancy.
        ecl_eval(r"""(defparameter *dev-null* (make-two-way-stream
                      (make-concatenated-stream) (make-broadcast-stream)))""")
        ecl_eval("(setf original-standard-output *standard-output*)")
        ecl_eval("(setf *standard-output* *dev-null*)")

        try:
            from sage.env import KENZO_FAS
            if KENZO_FAS:
                ecl_eval("(require :kenzo \"{}\")".format(KENZO_FAS))
            else:
                ecl_eval("(require :kenzo)")

        except RuntimeError:
            return FeatureTestResult(self, False, reason="Unable to make ECL require kenzo")
        return FeatureTestResult(self, True)

