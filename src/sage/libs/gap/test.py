"""
Short tests for GAP
"""

from sage.libs.gap.libgap import libgap
from sage.misc.temporary_file import tmp_filename


def check_write_to_file():
    """
    Test that libgap can write to files.

    See :issue:`16502`, :issue:`15833`.

    EXAMPLES::

        sage: from sage.libs.gap.test import check_write_to_file
        sage: check_write_to_file()
    """
    fname = tmp_filename()
    message = "Ceci n'est pas une groupe"
    libgap.PrintTo(fname, message)
    with open(fname) as f:
        assert f.read() == message
    SystemFile = libgap.function_factory('StringFile')
    assert SystemFile(fname).sage() == message
