# sage_setup: distribution = sagemath-repl
"""
Tests Deprecation

EXAMPLES::

    sage: import sage.tests.test_deprecation
    sage: sage.tests.test_deprecation.function_old()
    doctest:...: DeprecationWarning: function_old is deprecated. Please
    use sage.tests.test_deprecation.function_new instead.
    See https://github.com/sagemath/sage/issues/12345 for details.
"""
from sage.misc.superseded import deprecated_function_alias


def function_new():
    """
    New function, deprecating ``old_function``.

    EXAMPLES::

        sage: from sage.tests.test_deprecation import function_old
        sage: function_old()
        doctest:...: DeprecationWarning: function_old is deprecated. Please
        use sage.tests.test_deprecation.function_new instead.
        See https://github.com/sagemath/sage/issues/12345 for details.
    """
    pass


function_old = deprecated_function_alias(12345, function_new)
