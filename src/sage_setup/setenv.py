# Set some environment variables in the running process

import os
import sage.env
from pathlib import Path

def _environ_prepend(var, value, separator=':'):
    if value:
        if var in os.environ:
            os.environ[var] = value + separator + os.environ[var]
        else:
            os.environ[var] = value

def setenv():
    from sage.env import UNAME, SAGE_LOCAL, SAGE_VENV, SAGE_ARCHFLAGS, SAGE_PKG_CONFIG_PATH

    ##
    ## from sage-env:
    ##

    # not done: CC, CXX, FC, OBJC, OBJCXX, F77, F90, F95
    if 'ARCHFLAGS' not in os.environ and SAGE_ARCHFLAGS != "unset":
        os.environ['ARCHFLAGS'] = SAGE_ARCHFLAGS
    _environ_prepend('PKG_CONFIG_PATH', SAGE_PKG_CONFIG_PATH)
    # Issue #32057: As sage.env gives SAGE_LOCAL a fallback value from SAGE_VENV,
    # SAGE_LOCAL is never unset.  So we only set it if it differs from SAGE_VENV.
    # We assume that compiling/linking against libraries installed in SAGE_VENV
    # works -- that's the responsibility of the venv activation, not ours.
    if SAGE_LOCAL and Path(SAGE_VENV).resolve() != Path(SAGE_LOCAL).resolve():
        _environ_prepend('PATH',         f'{SAGE_LOCAL}/bin')
        _environ_prepend('LIBRARY_PATH', f'{SAGE_LOCAL}/lib')
        _environ_prepend('CPATH',        f'{SAGE_LOCAL}/include')
        _environ_prepend('LDFLAGS',      f'-L{SAGE_LOCAL}/lib -Wl,-rpath,{SAGE_LOCAL}/lib',
                         separator=' ')
        if UNAME == 'Linux':
            _environ_prepend('LDFLAGS',      f'-Wl,-rpath-link,{SAGE_LOCAL}/lib',
                             separator=' ')
        if Path(SAGE_VENV).resolve() != Path(SAGE_LOCAL).resolve():
            # This condition is always true, but we are keeping it for clarity.
            _environ_prepend('PATH',         f'{SAGE_VENV}/bin')
            # the following two are not done by sage-env
            #_environ_prepend('LIBRARY_PATH', f'{SAGE_VENV}/lib')
            #_environ_prepend('CPATH',        f'{SAGE_VENV}/include')

    # not done: PATH prepend of SAGE_SRC/bin, SAGE_ROOT/build/bin
    # not done: MACOSX_DEPLOYMENT_TARGET
    # not done: PATH prepend for ccache & CCACHE_BASEDIR
    # not done: OPENBLAS_NUM_THREADS
