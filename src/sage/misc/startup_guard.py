from contextlib import contextmanager

IS_STARTUP: bool = False

@contextmanager
def startup():
    """
    Simple context manager to indicate whether Sage is currently starting,
    e.g. importing `sage.all` etc.
    
    EXAMPLES::

        >>> import sage.misc.startup_guard as startup_guard
        >>> print(startup_guard.IS_STARTUP)
        False
        >>> with startup_guard.startup():
        >>>     print(startup_guard.IS_STARTUP)
        True
        >>> print(startup_guard.IS_STARTUP)
        False
    """
    global IS_STARTUP
    IS_STARTUP = True
    yield
    IS_STARTUP = False