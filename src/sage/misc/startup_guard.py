from contextlib import contextmanager
from enum import Enum


class StartupState(Enum):
    """
    Possible states of the startup guard.
    """

    UNINITIALIZED = 0
    """
    The startup process is not yet started.
    This either that we are at the very beginning of the startup process,
    or that another global environment other than :mod:`~sage.all` is used.
    """

    RUNNING = 1
    """
    The startup is in process.
    """

    FINISHED = 2
    """
    The startup process is finished.
    """


startup_state: StartupState = StartupState.UNINITIALIZED


@contextmanager
def startup():
    """
    Simple context manager to indicate whether Sage is currently starting,
    i.e. importing `sage.all`.

    EXAMPLES::

        sage: import sage.misc.startup_guard as startup_guard
        sage: print(startup_guard.startup_state) # This returns FINISHED, as sage.all is loaded for doctests.
        StartupState.FINISHED
        sage: with startup_guard.startup():
        ....:     print(startup_guard.startup_state)
        StartupState.RUNNING
        sage: print(startup_guard.startup_state)
        StartupState.FINISHED
    """
    global startup_state
    startup_state = StartupState.RUNNING
    yield
    startup_state = StartupState.FINISHED
