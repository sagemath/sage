# sage_setup: distribution = sagemath-repl
r"""
Evaluating a String in Sage
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy
import sage.repl.preparse as preparser


def sage_eval(source, locals=None, cmds='', preparse=True):
    r"""
    Obtain a Sage object from the input string by evaluating it using
    Sage. This means calling eval after preparsing and with globals
    equal to everything included in the scope of ``from sage.all
    import *``.).

    INPUT:

    - ``source`` -- string or object with a ``_sage_``
      method

    - ``locals`` -- evaluate in namespace of :mod:`sage.all` plus
      the locals dictionary

    - ``cmds`` -- string; sequence of commands to be run
      before source is evaluated

    - ``preparse`` -- boolean (default: ``True``); if ``True``, preparse the
      string expression

    EXAMPLES: This example illustrates that preparsing is applied::

        sage: eval('2^3')
        1
        sage: sage_eval('2^3')
        8

    However, preparsing can be turned off::

        sage: sage_eval('2^3', preparse=False)
        1

    Note that you can explicitly define variables and pass them as the
    second option::

        sage: x = PolynomialRing(RationalField(),"x").gen()
        sage: sage_eval('x^2+1', locals={'x':x})
        x^2 + 1

    This example illustrates that evaluation occurs in the context of
    ``from sage.all import *``. Even though ``bernoulli`` has
    been redefined in the local scope, when calling
    :func:`sage_eval` the default value meaning of :func:`bernoulli`
    is used. Likewise for ``QQ`` below::

        sage: bernoulli = lambda x : x^2
        sage: bernoulli(6)
        36
        sage: eval('bernoulli(6)')
        36
        sage: sage_eval('bernoulli(6)')                                                 # needs sage.libs.flint
        1/42

    ::

        sage: QQ = lambda x : x^2
        sage: QQ(2)
        4
        sage: sage_eval('QQ(2)')
        2
        sage: parent(sage_eval('QQ(2)'))
        Rational Field

    This example illustrates setting a variable for use in evaluation::

        sage: x = 5
        sage: eval('4//3 + x', {'x': 25})
        26
        sage: sage_eval('4/3 + x',  locals={'x': 25})
        79/3

    You can also specify a sequence of commands to be run before the
    expression is evaluated::

        sage: sage_eval('p', cmds='K.<x> = QQ[]\np = x^2 + 1')
        x^2 + 1

    If you give commands to execute and a dictionary of variables, then
    the dictionary will be modified by assignments in the commands::

        sage: vars = {}
        sage: sage_eval('None', cmds='y = 3', locals=vars)
        sage: vars['y'], parent(vars['y'])
        (3, Integer Ring)

    You can also specify the object to evaluate as a tuple. A 2-tuple
    is assumed to be a pair of a command sequence and an expression; a
    3-tuple is assumed to be a triple of a command sequence, an
    expression, and a dictionary holding local variables. (In this
    case, the given dictionary will not be modified by assignments in
    the commands.)

    ::

        sage: sage_eval(('f(x) = x^2', 'f(3)'))                                         # needs sage.symbolic
        9
        sage: vars = {'rt2': sqrt(2.0)}
        sage: sage_eval(('rt2 += 1', 'rt2', vars))
        2.41421356237309
        sage: vars['rt2']
        1.41421356237310

    This example illustrates how :mod:`sage_eval` can be
    useful when evaluating the output of other computer algebra
    systems::

        sage: # needs sage.libs.gap
        sage: R.<x> = PolynomialRing(RationalField())
        sage: gap.eval('R:=PolynomialRing(Rationals,["x"]);')
        'Rationals[x]'
        sage: ff = gap.eval('x:=IndeterminatesOfPolynomialRing(R);; f:=x^2+1;'); ff
        'x^2+1'
        sage: sage_eval(ff, locals={'x':x})
        x^2 + 1
        sage: eval(ff)
        Traceback (most recent call last):
        ...
        RuntimeError: Use ** for exponentiation, not '^', which means xor
        in Python, and has the wrong precedence.

    Here you can see that :func:`eval` simply will not work but
    :func:`sage_eval` will.

    TESTS:

    We get a nice minimal error message for syntax errors, that still
    points to the location of the error (in the input string)::

        sage: sage_eval('RR(22/7]')
        Traceback (most recent call last):
        ...
         File "<string>", line 1
            RR(Integer(22)/Integer(7)]
                                     ^
        SyntaxError: ...

    ::

        sage: sage_eval('None', cmds='$x = $y[3] # Does Perl syntax work?')
        Traceback (most recent call last):
        ...
         File "<string>", line 1
            $x = $y[Integer(3)] # Does Perl syntax work?
            ^
        SyntaxError: invalid ...
    """
    if isinstance(source, (list, tuple)):
        cmds = source[0]
        if len(source) > 2:
            locals = copy(source[2])
        source = source[1]

    if not isinstance(source, str):
        raise TypeError("source must be a string.")

    if locals is None:
        locals = {}

    try:
        import sage.all as toplevel
    except ImportError:
        try:
            import sage.all__sagemath_modules as toplevel
        except ImportError:
            try:
                import sage.all__sagemath_categories as toplevel
            except ImportError:
                import sage.all__sagemath_objects as toplevel
    if cmds:
        cmd_seq = cmds + '\n_sage_eval_returnval_ = ' + source
        if preparse:
            cmd_seq = preparser.preparse_file(cmd_seq)
    else:
        if preparse:
            source = preparser.preparse(source)

    if cmds:
        exec(cmd_seq, toplevel.__dict__, locals)
        return locals['_sage_eval_returnval_']
    else:
        return eval(source, toplevel.__dict__, locals)


def sageobj(x, vars=None):
    """
    Return a native Sage object associated to ``x``, if possible and
    implemented.

    If the object has a ``_sage_`` method it is called and the value is
    returned. Otherwise, :func:`str` is called on the object, and all preparsing
    is applied and the resulting expression is evaluated in the context
    of ``from sage.all import *``. To evaluate the
    expression with certain variables set, use the ``vars`` argument, which
    should be a dictionary.

    EXAMPLES::

        sage: type(sageobj(gp('34/56')))                                                # needs sage.libs.pari
        <class 'sage.rings.rational.Rational'>

        sage: n = 5/2
        sage: sageobj(n) is n
        True
        sage: k = sageobj('Z(8^3/1)', {'Z':ZZ}); k
        512
        sage: type(k)
        <class 'sage.rings.integer.Integer'>

    This illustrates interfaces::

        sage: # needs sage.libs.pari
        sage: f = gp('2/3')
        sage: type(f)
        <class 'sage.interfaces.gp.GpElement'>
        sage: f._sage_()
        2/3
        sage: type(f._sage_())
        <class 'sage.rings.rational.Rational'>

        sage: # needs sage.libs.gap
        sage: a = gap(939393/2433)
        sage: a._sage_()
        313131/811
        sage: type(a._sage_())
        <class 'sage.rings.rational.Rational'>
    """
    try:
        return x._sage_()
    except (TypeError, NotImplementedError, AttributeError):
        return sage_eval(str(x), vars)
