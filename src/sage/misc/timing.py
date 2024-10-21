# sage_setup: distribution = sagemath-objects
r"""
Timing functions
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 Gonzalo Tornaria
#                     2008 Martin Albrecht
#                     2009 Mike Hansen
#                     2018 Frédéric Chapoton
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import time


def cputime(t=0, subprocesses=False):
    """
    Return the time in CPU seconds since Sage started, or with
    optional argument ``t``, return the time since ``t``. This is how
    much time Sage has spent using the CPU.  If ``subprocesses=False``
    this does not count time spent in subprocesses spawned by Sage
    (e.g., Gap, Singular, etc.). If ``subprocesses=True`` this
    function tries to take all subprocesses with a working
    ``cputime()`` implementation into account.

    The measurement for the main Sage process is done via a call to
    :func:`resource.getrusage()`.

    INPUT:

    - ``t`` -- (optional) time in CPU seconds, if ``t`` is a result
      from an earlier call with ``subprocesses=True``, then
      ``subprocesses=True`` is assumed.

    - ``subprocesses`` -- boolean (default: ``False``); include subprocesses

    OUTPUT:

    - ``float`` -- time in CPU seconds if ``subprocesses=False``

    - :class:`GlobalCputime` -- object which holds CPU times of
      subprocesses otherwise

    EXAMPLES::

        sage: t = cputime()
        sage: F = gp.factor(2^199-1)                                                    # needs sage.libs.pari
        sage: cputime(t)            # somewhat random
        0.010999000000000092

        sage: t = cputime(subprocesses=True)
        sage: F = gp.factor(2^199-1)                                                    # needs sage.libs.pari
        sage: cputime(t)            # somewhat random
        0.091999

        sage: w = walltime()
        sage: F = gp.factor(2^199-1)                                                    # needs sage.libs.pari
        sage: walltime(w)           # somewhat random
        0.58425593376159668

    .. NOTE::

        Even with ``subprocesses=True`` there is no guarantee that the
        CPU time is reported correctly because subprocesses can be
        started and terminated at any given time.
    """
    try:
        import resource
    except ImportError:
        # The module 'resource' is removed in Pyodide to browser limitations.
        return walltime(t)

    if isinstance(t, GlobalCputime):
        subprocesses = True

    if not subprocesses:
        try:
            t = float(t)
        except TypeError:
            t = 0.0
        u, s = resource.getrusage(resource.RUSAGE_SELF)[:2]
        return u + s - t
    else:
        try:
            from sage.interfaces.quit import expect_objects
        except ImportError:
            expect_objects = ()
        if t == 0:
            ret = GlobalCputime(cputime())
            for s in expect_objects:
                S = s()
                if S and S.is_running():
                    try:
                        ct = S.cputime()
                        ret.total += ct
                        ret.interfaces[s] = ct
                    except NotImplementedError:
                        pass
            return ret
        else:
            if not isinstance(t, GlobalCputime):
                t = GlobalCputime(t)
            ret = GlobalCputime(cputime() - t.local)
            for s in expect_objects:
                S = s()
                if S and S.is_running():
                    try:
                        ct = S.cputime() - t.interfaces.get(s, 0.0)
                        ret.total += ct
                        ret.interfaces[s] = ct
                    except NotImplementedError:
                        pass
            return ret


class GlobalCputime:
    """
    Container for CPU times of subprocesses.

    AUTHOR:

    - Martin Albrecht - (2008-12): initial version

    EXAMPLES:

    Objects of this type are returned if ``subprocesses=True`` is
    passed to :func:`cputime`::

        sage: cputime(subprocesses=True)        # indirect doctest, output random
        0.2347431

    We can use it to keep track of the CPU time spent in Singular for
    example::

        sage: t = cputime(subprocesses=True)
        sage: P = PolynomialRing(QQ,7,'x')
        sage: I = sage.rings.ideal.Katsura(P)                                           # needs sage.libs.singular
        sage: gb = I.groebner_basis()  # calls Singular                                 # needs sage.libs.singular
        sage: cputime(subprocesses=True) - t    # output random
        0.462987

    For further processing we can then convert this container to a
    float::

        sage: t = cputime(subprocesses=True)
        sage: float(t)                          # output somewhat random
        2.1088339999999999

    .. SEEALSO::

      :func:`cputime`
    """
    def __init__(self, t):
        """
        Create a new CPU time object which also keeps track of
        subprocesses.

        EXAMPLES::

            sage: from sage.misc.timing import GlobalCputime
            sage: ct = GlobalCputime(0.0); ct
            0.0...
        """
        self.total = t
        self.local = t
        self.interfaces = {}

    def __repr__(self):
        """
        EXAMPLES::

            sage: cputime(subprocesses=True)    # indirect doctest, output random
            0.2347431
        """
        return str(self.total)

    def __add__(self, other):
        """
        EXAMPLES::

            sage: t = cputime(subprocesses=True)
            sage: P = PolynomialRing(QQ,7,'x')
            sage: I = sage.rings.ideal.Katsura(P)                                       # needs sage.libs.singular
            sage: gb = I.groebner_basis()  # calls Singular                             # needs sage.libs.singular
            sage: cputime(subprocesses=True) + t # output random
            2.798708
        """
        if not isinstance(other, GlobalCputime):
            other = GlobalCputime(other)
        ret = GlobalCputime(self.total + other.total)
        return ret

    def __sub__(self, other):
        """
        EXAMPLES::

            sage: t = cputime(subprocesses=True)
            sage: P = PolynomialRing(QQ,7,'x')
            sage: I = sage.rings.ideal.Katsura(P)                                       # needs sage.libs.singular
            sage: gb = I.groebner_basis()  # calls Singular                             # needs sage.libs.singular
            sage: cputime(subprocesses=True) - t # output random
            0.462987
        """
        if not isinstance(other, GlobalCputime):
            other = GlobalCputime(other)
        ret = GlobalCputime(self.total - other.total)
        return ret

    def __float__(self):
        """
        EXAMPLES::

            sage: t = cputime(subprocesses=True)
            sage: float(t)                          # output somewhat random
            2.1088339999999999
        """
        return float(self.total)


def walltime(t=0):
    """
    Return the wall time in second, or with optional argument ``t``, return
    the wall time since time ``t``. "Wall time" means the time on a wall
    clock, i.e., the actual time.

    INPUT:

    - ``t`` -- (optional) float, time in CPU seconds

    OUTPUT: ``float`` -- time in seconds

    EXAMPLES::

        sage: w = walltime()
        sage: F = factor(2^199-1)                                                       # needs sage.libs.pari
        sage: walltime(w)   # somewhat random
        0.8823847770690918
    """
    return time.time() - t
