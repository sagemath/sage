# cython: old_style_globals=True
"""
Dependency usage tracking for citations
"""

from sage.misc.temporary_file import tmp_filename
from sage.env import SAGE_LOCAL, SAGE_VENV

systems = {}
systems['PARI'] = ['cypari2', 'sage.interfaces.gp']
systems['Singular'] = ['sage.interfaces.singular', '_libsingular',
                       'sage.libs.singular']
systems['Maxima'] = ['sage.interfaces.maxima']
systems['GAP'] = ['sage.interfaces.gap']
systems['Magma'] = ['sage.interfaces.magma', 'sage.interfaces.magma_free']
systems['Axiom'] = ['sage.interfaces.axiom']
systems['ECM'] = ['sage.interfaces.ecm']
systems['scipy'] = ['scipy']
systems['numpy'] = ['numpy']
systems['ginac'] = ['sage.symbolic']
systems['Maple'] = ['sage.interfaces.maple']
systems['Mathematica'] = ['sage.interfaces.mathematica']
systems['MuPAD'] = ['sage.interfaces.mupad']
systems['Octave'] = ['sage.interfaces.octave']
systems['povray'] = ['sage.interfaces.povray']
systems['Macaulay2'] = ['sage.interfaces.macaulay2']
systems['mwrank'] = ['sage.interfaces.mwrank', 'sage.libs.eclib']
systems['matlab'] = ['sage.interfaces.matlab']
systems['LiE'] = ['sage.interfaces.lie']
systems['Tachyon'] = ['sage.interfaces.tachyon']
systems['Frobby'] = ['sage.interfaces.frobby']
systems['gfan'] = ['sage.interfaces.gfan']
systems['R'] = ['sage.interfaces.r']
systems['KASH'] = ['sage.interfaces.kash']
systems['Linbox'] = ['sage.libs.linbox']
systems['Symmetrica'] = ['sage.libs.symmetrica']
systems['NTL'] = ['sage.libs.ntl',
                  'sage.rings.finite_rings.element_ntl_gf2e']
systems['FLINT'] = ['_flint']
systems['GMP'] = ['sage.rings.integer.Integer']
systems['MPFR'] = ['sage.rings.real_mpfr',
                   'sage.rings.complex_mpfr']
systems['MPFI'] = ['sage.rings.real_mpfi',
                   'sage.rings.complex_interval']
systems['M4RI'] = ['sage.matrix.matrix_mod2_dense']
systems['Givaro'] = ['sage.rings.finite_rings.element_givaro']
systems['PolyBoRi'] = ['sage.rings.polynomial.pbori']


def get_systems(cmd):
    """
    Return a list of the systems used in running the command ``cmd``.

    Note that the results can sometimes include systems that did not
    actually contribute to the computation. Due to caching, it
    could miss some dependencies as well.

    INPUT:

    - ``cmd`` -- string to run

    .. WARNING::

        In order to properly support Cython code, this requires that
        Sage was compiled with the environment variable
        ``SAGE_PROFILE=yes``. If this was not the case, a warning will
        be given when calling this function.

    EXAMPLES::

        sage: from sage.misc.citation import get_systems
        sage: get_systems('print("hello")')  # random (may print warning)
        []
        sage: integrate(x^2, x)  # Priming coercion model                               # needs sage.symbolic
        1/3*x^3
        sage: get_systems('integrate(x^2, x)')                                          # needs sage.symbolic
        ['Maxima', 'ginac']
        sage: R.<x,y,z> = QQ[]
        sage: I = R.ideal(x^2+y^2, z^2+y)
        sage: get_systems('I.primary_decomposition()')
        ['Singular']
    """
    import cProfile
    import pstats
    import re

    if not cython_profile_enabled():
        from warnings import warn
        warn("get_systems() requires Cython profiling to be enabled, "
             "otherwise the results will be very unreliable. "
             "Rebuild Sage with the environment variable 'SAGE_PROFILE=yes' "
             "to enable profiling.")

    if not isinstance(cmd, str):
        raise TypeError("command must be a string")

    from sage.repl.preparse import preparse
    cmd = preparse(cmd)

    # Run the command and get the stats
    filename = tmp_filename()
    cProfile.runctx(cmd, globals(), {}, filename)
    stats = pstats.Stats(filename)

    # Strings is a list of method names and modules which get run
    def string_from_stat(a):
        s = a[0]
        if SAGE_LOCAL:
            s = s.replace(SAGE_LOCAL, "")
        if SAGE_VENV:
            s = s.replace(SAGE_VENV, "")
        return s + " " + a[2]

    strings = [string_from_stat(a)
               for a in stats.stats]

    # Remove trivial functions
    bad_res = [re.compile(r'is_.*Element'), re.compile("is_[a-z_]*_type")]
    for bad_re in bad_res:
        i = 0
        while i < len(strings):
            if bad_re.findall(strings[i]):
                strings.pop(i)
            else:
                i += 1

    # Check to see which systems appear in the profiled run
    systems_used = []
    for system in systems:
        if any((r in s) or (r.replace('.', '/') in s)
               for r in systems[system] for s in strings):
            systems_used.append(system)
    return sorted(systems_used)


cdef extern from *:
    int CYTHON_PROFILE """
        #ifdef CYTHON_PROFILE
        CYTHON_PROFILE
        #else
        0
        #endif
        """


cpdef inline bint cython_profile_enabled() noexcept:
    """
    Return whether Cython profiling is enabled.

    EXAMPLES::

        sage: from sage.misc.citation import cython_profile_enabled
        sage: cython_profile_enabled()  # random
        False
        sage: type(cython_profile_enabled()) is bool
        True
    """
    return CYTHON_PROFILE
