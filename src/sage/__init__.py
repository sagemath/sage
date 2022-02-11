# Do not add anything to this file.
# It will be removed soon in order to turn 'sage' into a native namespace package.
# See https://trac.sagemath.org/ticket/29705



# Deprecated leftover of monkey-patching inspect.isfunction() to support Cython functions.
# We cannot use lazy_import for the deprecation here.
def isfunction(obj):
    """
    Check whether something is a function.

    This is a variant of ``inspect.isfunction``:
    We assume that anything which has a genuine ``__code__``
    attribute (not using ``__getattr__`` overrides) is a function.
    This is meant to support Cython functions.

    This function is deprecated.  Most uses of ``isfunction``
    can be replaced by ``callable``.

    EXAMPLES::

        sage: from sage import isfunction
        sage: def f(): pass
        sage: isfunction(f)
        doctest:warning...
        DeprecationWarning: sage.isfunction is deprecated; use callable or sage.misc.sageinspect.is_function_or_cython_function instead
        See https://trac.sagemath.org/32479 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(32479, "sage.isfunction is deprecated; use callable or sage.misc.sageinspect.is_function_or_cython_function instead")
    from sage.misc.sageinspect import is_function_or_cython_function
    return is_function_or_cython_function(obj)
