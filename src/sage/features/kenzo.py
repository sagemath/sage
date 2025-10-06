# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``kenzo``
"""

# ****************************************************************************
#       Copyright (C) 2020 Travis Scrimshaw
#                     2021 Matthias Koeppe
#                     2021 Michael Orlitzky
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from . import Feature, FeatureTestResult


class Kenzo(Feature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`Kenzo <spkg_kenzo>`.

    EXAMPLES::

        sage: from sage.features.kenzo import Kenzo
        sage: Kenzo().is_present()  # optional - kenzo
        FeatureTestResult('kenzo', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.kenzo import Kenzo
            sage: isinstance(Kenzo(), Kenzo)
            True
        """
        Feature.__init__(self, name='kenzo', spkg='kenzo',
                         url='https://github.com/miguelmarco/kenzo/')

    def _is_present(self):
        r"""
        Check whether Kenzo is installed and works.

        EXAMPLES::

            sage: from sage.features.kenzo import Kenzo
            sage: Kenzo()._is_present()  # optional - kenzo
            FeatureTestResult('kenzo', True)
        """
        try:
            from sage.libs.ecl import ecl_eval
        except ImportError:
            return FeatureTestResult(self, False, reason="sage.libs.ecl is not available")
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


def all_features():
    return [Kenzo()]
