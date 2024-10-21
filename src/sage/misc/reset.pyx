# sage_setup: distribution = sagemath-objects
# cython: old_style_globals=True
"""
Interpreter reset
"""
import sys

# Exclude these from the reset command.
# DATA, base64 -- needed by the notebook
# Add exit and quit to EXCLUDE to resolve Issue #22529 and Issue #16704
EXCLUDE = set(['sage_mode', '__DIR__', 'DIR', 'DATA', 'base64', 'exit', 'quit'])


def reset(vars=None, attached=False):
    """
    Delete all user-defined variables, reset all global variables
    back to their default states, and reset all interfaces to other
    computer algebra systems.

    If vars is specified, just restore the value of vars and leave
    all other variables alone (i.e., call restore).

    Note that the variables in the set :obj:`sage.misc.reset.EXCLUDE` are
    excluded from being reset.

    INPUT:

    - ``vars`` -- list or space or comma separated string (default:
      ``None``); variables to restore

    - ``attached`` -- boolean (default: ``False``); if ``vars`` is not ``None``,
      whether to detach all attached files

    EXAMPLES::

        sage: x = 5
        sage: reset()
        sage: x                                                                         # needs sage.symbolic
        x

        sage: fn = tmp_filename(ext='.py')
        sage: sage.misc.reset.EXCLUDE.add('fn')
        sage: with open(fn, 'w') as f:
        ....:     _ = f.write('a = 111')
        sage: attach(fn)
        sage: af = attached_files(); len(af)
        1
        sage: af == [fn]
        True
        sage: reset()
        sage: af = attached_files(); len(af)
        1
        sage: reset(attached=True)
        sage: af = attached_files(); len(af)
        0
        sage: sage.misc.reset.EXCLUDE.remove('fn')

    TESTS:

    Confirm that assumptions do not survive a reset (:issue:`10855`)::

        sage: # needs sage.symbolic
        sage: assume(x > 3)
        sage: assumptions()
        [x > 3]
        sage: bool(x > 3)
        True
        sage: reset()
        sage: assumptions()
        []
        sage: bool(x > 3)
        False
    """
    if vars is not None:
        restore(vars)
        return
    G = globals()  # this is the reason the code must be in Cython.
    T = type(sys)
    for k in list(G):
        if k[0] != '_' and not isinstance(k, T) and k not in EXCLUDE:
            try:
                del G[k]
            except KeyError:
                pass
    restore()
    try:
        from sage.symbolic.assumptions import forget
    except ImportError:
        pass
    else:
        forget()
    reset_interfaces()
    if attached:
        import sage.repl.attach
        sage.repl.attach.reset()


def restore(vars=None):
    """
    Restore predefined global variables to their default values.

    INPUT:

    - ``vars`` -- string or list (default: ``None``); if not ``None``, restores
      just the given variables to the default value

    EXAMPLES::

        sage: x = 10; y = 15/3; QQ='red'
        sage: QQ
        'red'
        sage: restore('QQ')
        sage: QQ
        Rational Field
        sage: x
        10
        sage: y = var('y')                                                              # needs sage.symbolic
        sage: restore('x y')
        sage: x                                                                         # needs sage.symbolic
        x
        sage: y
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined
        sage: x = 10; y = 15/3; QQ='red'
        sage: ww = 15
        sage: restore()
        sage: x, QQ, ww                                                                 # needs sage.symbolic
        (x, Rational Field, 15)
        sage: restore('ww')
        sage: ww
        Traceback (most recent call last):
        ...
        NameError: name 'ww' is not defined
    """
    G = globals()  # this is the reason the code must be in Cython.
    if 'sage_mode' not in G:
        import sage.all
        D = sage.all.__dict__
    else:
        mode = G['sage_mode']
        if mode == 'cmdline':
            import sage.all_cmdline
            D = sage.all_cmdline.__dict__
        else:
            import sage.all
            D = sage.all.__dict__
    _restore(G, D, vars)
    try:
        import sage.calculus.calculus
    except ImportError:
        pass
    else:
        _restore(sage.calculus.calculus.syms_cur, sage.calculus.calculus.syms_default, vars)


def _restore(G, D, vars):
    if vars is None:
        for k, v in D.iteritems():
            G[k] = v
    else:
        if isinstance(vars, str):
            if ',' in vars:
                vars = vars.split(',')
            else:
                vars = vars.split()
        for k in vars:
            if k in D:
                G[k] = D[k]
            else:
                try:
                    del G[k]      # the default value was "unset"
                except KeyError:
                    pass


def reset_interfaces():
    try:
        from sage.interfaces.quit import expect_quitall
    except ImportError:
        pass
    else:
        expect_quitall()
