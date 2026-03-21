# sage.doctest: optional - giac
r"""
Context manager and local wrapper for giacsettings.
"""

from sage.libs.giac.giac import giacsettings, libgiac  # noqa: F401


class GiacSettingsDefaultContext:
    r"""
    Context preserve libgiac settings.
    """

    def __enter__(self):
        """
        EXAMPLES::

           sage: from sage.libs.giac.context import GiacSettingsDefaultContext
           sage: from sage.libs.giac.giac import giacsettings
           sage: giacsettings.proba_epsilon = 1e-16
           sage: with GiacSettingsDefaultContext(): giacsettings.proba_epsilon = 1e-12
           sage: giacsettings.proba_epsilon < 1e-14
           True

        """
        self.proba_epsilon = giacsettings.proba_epsilon
        self.threads = giacsettings.threads
        # Change the debug level at the end to not have messages at each modification
        self.debuginfolevel = libgiac('debug_infolevel()')

    def __exit__(self, typ, value, tb):
        """
        EXAMPLES::

           sage: from sage.libs.giac.context import GiacSettingsDefaultContext
           sage: from sage.libs.giac.giac import giacsettings
           sage: giacsettings.proba_epsilon = 1e-16
           sage: with GiacSettingsDefaultContext(): giacsettings.proba_epsilon = 1e-30
           sage: giacsettings.proba_epsilon < 1e-20
           False
        """
        # Restore the debug level first to not have messages at each modification
        libgiac('debug_infolevel')(self.debuginfolevel)
        # NB: giacsettings.epsilon has a different meaning that giacsettings.proba_epsilon.
        giacsettings.proba_epsilon = self.proba_epsilon
        giacsettings.threads = self.threads


def local_giacsettings(func):
    """
    Decorator to preserve Giac's proba_epsilon and threads settings.

    EXAMPLES::

        sage: def testf(a, b):
        ....:     giacsettings.proba_epsilon = a/100
        ....:     giacsettings.threads = b+2
        ....:     return (giacsettings.proba_epsilon, giacsettings.threads)

        sage: from sage.libs.giac.giac import giacsettings
        sage: from sage.libs.giac.context import local_giacsettings
        sage: gporig, gtorig = (giacsettings.proba_epsilon,giacsettings.threads)
        sage: gp, gt = local_giacsettings(testf)(giacsettings.proba_epsilon,giacsettings.threads)
        sage: gporig == giacsettings.proba_epsilon
        True
        sage: gtorig == giacsettings.threads
        True
        sage: gp<gporig, gt-gtorig
        (True, 2)
    """
    from sage.misc.decorators import sage_wraps

    @sage_wraps(func)
    def wrapper(*args, **kwds):
        """
        Execute function in ``GiacSettingsDefaultContext``.
        """
        with GiacSettingsDefaultContext():
            return func(*args, **kwds)

    return wrapper
