"""
Tests Deprecation

EXAMPLES::

    sage: import sage.tests.deprecation_test
    sage: sage.tests.deprecation_test.function_old()
    doctest:...: DeprecationWarning: function_old is deprecated. Please
    use sage.tests.deprecation_test.function_new instead.
    See http://trac.sagemath.org/12345 for details.
"""
from sage.misc.superseded import deprecated_function_alias


def function_new():
    """
    New function, deprecating ``old_function``.

    EXAMPLES::

        sage: from sage.tests.deprecation_test import function_old
        sage: function_old()
        doctest:...: DeprecationWarning: function_old is deprecated. Please
        use sage.tests.deprecation_test.function_new instead.
        See http://trac.sagemath.org/12345 for details.
    """
    pass


function_old = deprecated_function_alias(12345, function_new)
