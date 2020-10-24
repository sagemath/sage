def create_RealField(prec=53, type="MPFR", rnd="RNDN", sci_not=0):
    """
    Create a real field with given precision, type, rounding mode and
    scientific notation.

    Some options are ignored for certain types (RDF for example).

    INPUT:

    - ``prec`` -- a positive integer

    - ``type`` -- type of real field:

      - ``'RDF'`` -- the Sage real field corresponding to native doubles
      - ``'Interval'`` -- real fields implementing interval arithmetic
      - ``'RLF'`` -- the real lazy field
      - ``'MPFR'`` -- floating point real numbers implemented using the MPFR
        library

    - ``rnd`` -- rounding mode:

      - ``'RNDN'`` -- round to nearest
      - ``'RNDZ'`` -- round toward zero
      - ``'RNDD'`` -- round down
      - ``'RNDU'`` -- round up

    - ``sci_not`` -- boolean, whether to use scientific notation for printing

    OUTPUT:

    the appropriate real field

    EXAMPLES::

        sage: from sage.rings.real_field import create_RealField
        sage: create_RealField(30)
        Real Field with 30 bits of precision
        sage: create_RealField(20, 'RDF') # ignores precision
        Real Double Field
        sage: create_RealField(60, 'Interval')
        Real Interval Field with 60 bits of precision
        sage: create_RealField(40, 'RLF') # ignores precision
        Real Lazy Field
    """
    if type == "RDF":
        from .real_double import RDF
        return RDF
    elif type == "Interval":
        from .real_mpfi import RealIntervalField
        return RealIntervalField(prec, sci_not)
    elif type == "Ball":
        from .real_arb import RealBallField
        return RealBallField(prec)
    elif type == "RLF":
        from .real_lazy import RLF
        return RLF
    else:
        from .real_mpfr import RealField
        return RealField(prec, sci_not, rnd)


