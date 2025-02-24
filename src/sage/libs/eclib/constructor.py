"Cremona modular symbols"

def CremonaModularSymbols(level, sign=0, cuspidal=False, verbose=0):
    """
    Return the space of Cremona modular symbols with given level, sign, etc.

    INPUT:

    - ``level`` -- integer; at least 2, not just positive!
    - ``sign`` -- integer (default: 0); either 0 or 1 or -1
    - ``cuspidal`` -- boolean (default: ``False``); if ``True``, compute only
      the cuspidal subspace
    - ``verbose`` -- boolean (default: ``False``); if ``True``, print verbose
      information while creating space

    EXAMPLES::

        sage: M = CremonaModularSymbols(43); M
        Cremona Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0
        sage: M = CremonaModularSymbols(43, sign=1); M
        Cremona Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1
        sage: M = CremonaModularSymbols(43, cuspidal=True); M
        Cremona Cuspidal Modular Symbols space of dimension 6 for Gamma_0(43) of weight 2 with sign 0
        sage: M = CremonaModularSymbols(43, cuspidal=True, sign=1); M
        Cremona Cuspidal Modular Symbols space of dimension 3 for Gamma_0(43) of weight 2 with sign 1

    When run interactively, the following command will display verbose output::

        sage: M = CremonaModularSymbols(43, verbose=1)
        After 2-term relations, ngens = 22
        ngens     = 22
        maxnumrel = 32
        relation matrix has = 704 entries...
        Finished 3-term relations: numrel = 16 ( maxnumrel = 32)
        relmat has 42 nonzero entries (density = 0.0596591)
        Computing kernel...
        time to compute kernel =  (... seconds)
        rk = 7
        Number of cusps is 2
        ncusps = 2
        About to compute matrix of delta
        delta matrix done: size 2x7.
        About to compute kernel of delta
        done
        Finished constructing homspace.
        sage: M
        Cremona Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0

    The input must be valid or a :exc:`ValueError` is raised::

        sage: M = CremonaModularSymbols(-1)
        Traceback (most recent call last):
        ...
        ValueError: the level (= -1) must be at least 2
        sage: M = CremonaModularSymbols(0)
        Traceback (most recent call last):
        ...
        ValueError: the level (= 0) must be at least 2

    The sign can only be 0 or 1 or -1::

        sage: M = CremonaModularSymbols(10, sign = -2)
        Traceback (most recent call last):
        ...
        ValueError: sign (= -2) is not supported; use 0, +1 or -1

    We do allow -1 as a sign (see :issue:`9476`)::

        sage: CremonaModularSymbols(10, sign = -1)
        Cremona Modular Symbols space of dimension 0 for Gamma_0(10) of weight 2 with sign -1
    """
    from .homspace import ModularSymbols
    return ModularSymbols(level=level, sign=sign, cuspidal=cuspidal, verbose=verbose)
