Test raising slow doctest warning (cputime instead of walltime is checked so we need busy loop)::

    sage: t = walltime()
    sage: while walltime(t) < 2: pass
