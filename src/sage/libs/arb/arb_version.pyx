# -*- coding: utf-8
from sage.libs.arb.arb cimport arb_version 
from sage.cpython.string cimport char_to_str

def version():
    """
    Get arb version

    TESTS::

        sage: from sage.libs.arb.arb_version import version
        sage: version().split('.')[0]
        '2'
    """
    try:
        py_string = char_to_str(arb_version)
    finally:
        pass
    return py_string
