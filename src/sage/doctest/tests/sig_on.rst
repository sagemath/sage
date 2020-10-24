Test a bad value of ``sig_on_count``::

    sage: from cysignals.tests import _sig_on
    sage: _sig_on()

The following test should succeed as usual::

    sage: 3**12 - 2**19
    7153
