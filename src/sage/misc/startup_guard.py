from contextlib import contextmanager

IS_STARTUP: bool = False

@contextmanager
def startup():
    """
    Simple context manager to indicate whether Sage is currently starting,
    i.e. importing `sage.all`.
    
    EXAMPLES::

        sage: import sage.misc.startup_guard as startup_guard
        sage: print(startup_guard.IS_STARTUP)
        False
        sage: with startup_guard.startup():
        ....:     print(startup_guard.IS_STARTUP)
        True
        sage: print(startup_guard.IS_STARTUP)
        False
    """
    global IS_STARTUP
    IS_STARTUP = True
    yield
    IS_STARTUP = False
